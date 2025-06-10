import asyncio
from fastapi import WebSocket, WebSocketDisconnect

from sqlmodel import Session

from .db.db import engine
from .models import DockingJob, JobStatus, DockingJobWComp


job_update_queues: dict[str, asyncio.Queue] = {}

async def get_job_status(websocket: WebSocket, job_id: str):
    await websocket.accept()

    job_queue = job_update_queues.setdefault(job_id, asyncio.Queue())

    with Session(engine) as session:
        job = session.get(DockingJob, job_id)
        if job:
            job_dict = DockingJobWComp.model_validate(job)
            await websocket.send_text(job_dict.model_dump_json())

    try:
        while True:
            await job_queue.get()

            with Session(engine) as session:
                job = session.get(DockingJob, job_id)
                if not job:
                    await websocket.send_json({"error": "Job not found"})
                    await websocket.close()
                    break

                job_dict = DockingJobWComp.model_validate(job)
                await websocket.send_text(job_dict.model_dump_json())

                if job.job_status in (JobStatus.COMPLETED, JobStatus.FAILED):
                    print(f"Job break for job_id={job_id}")
                    break

    except WebSocketDisconnect:
        print(f"Client disconnected: {job_id}")
    finally:
        await websocket.close()
