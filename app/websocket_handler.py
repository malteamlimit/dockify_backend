import asyncio

from fastapi import WebSocket, WebSocketDisconnect
from sqlmodel import Session, select

from .db import db as _db
from .models import DockingJob, JobStatus, DockingJobWComp


# Every connected stream client gets its own queue of changed job_ids.
_client_queues: set[asyncio.Queue] = set()

# The main event loop, captured at startup so worker threads (where docking
# runs) can schedule broadcasts onto it.
_main_loop: asyncio.AbstractEventLoop | None = None


def set_event_loop(loop: asyncio.AbstractEventLoop) -> None:
    """Store the main event loop so background threads can broadcast updates."""
    global _main_loop
    _main_loop = loop


def notify_job_update(job_id: str) -> None:
    """Broadcast a job change to every connected stream client.

    Safe to call from any thread: the queue put is scheduled onto the main
    event loop. A no-op while no loop has been registered yet.
    """
    if _main_loop is None:
        return
    for queue in list(_client_queues):
        asyncio.run_coroutine_threadsafe(queue.put(job_id), _main_loop)


def recover_orphaned_jobs() -> None:
    """Fail jobs left mid-flight by a previous backend run.

    A QUEUED/RUNNING status in the DB at startup means the process died before
    the job finished -- no task is executing it anymore, so it would otherwise
    hang forever. Mark such jobs as failed so the frontend reflects reality.
    """
    with Session(_db.engine) as session:
        statement = select(DockingJob).where(
            DockingJob.job_status.in_([JobStatus.QUEUED, JobStatus.RUNNING])
        )
        for job in session.exec(statement).all():
            job.job_status = JobStatus.FAILED
            job.error = "Job interrupted by a backend restart."
            job.runs = len(job.complexes)
            session.add(job)
        session.commit()


async def job_status_stream(websocket: WebSocket) -> None:
    """Single global stream of job updates across all jobs.

    On connect the client receives a snapshot of every job (one message per
    job); afterwards it receives one message per job whenever that job changes.
    A deleted job is signalled with ``{"deleted": <job_id>}``.
    """
    await websocket.accept()

    queue: asyncio.Queue = asyncio.Queue()
    _client_queues.add(queue)

    # The client never sends data; this task exists only to notice a
    # disconnect, even while we are idle waiting on the queue.
    async def watch_disconnect() -> None:
        try:
            while True:
                message = await websocket.receive()
                if message["type"] == "websocket.disconnect":
                    return
        except WebSocketDisconnect:
            return

    disconnect_task = asyncio.create_task(watch_disconnect())

    try:
        # Initial snapshot so a freshly (re)connected client is consistent
        # even if updates happened between its REST fetch and this connect.
        with Session(_db.engine) as session:
            for job in session.exec(select(DockingJob)).all():
                payload = DockingJobWComp.model_validate(job).model_dump_json()
                await websocket.send_text(payload)

        # Stream deltas: one job per message, until the client disconnects.
        while not disconnect_task.done():
            next_update = asyncio.ensure_future(queue.get())
            done, _ = await asyncio.wait(
                {next_update, disconnect_task},
                return_when=asyncio.FIRST_COMPLETED,
            )
            if next_update not in done:
                next_update.cancel()
                break

            job_id = next_update.result()
            with Session(_db.engine) as session:
                job = session.get(DockingJob, job_id)
                if job is None:
                    await websocket.send_json({"deleted": job_id})
                else:
                    payload = DockingJobWComp.model_validate(job).model_dump_json()
                    await websocket.send_text(payload)

    except WebSocketDisconnect:
        pass
    finally:
        disconnect_task.cancel()
        _client_queues.discard(queue)
