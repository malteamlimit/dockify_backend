import os

from fastapi import APIRouter, BackgroundTasks, Depends, WebSocket, status, HTTPException
from sqlmodel import Session, select

from .util import generate_sdf_from_smiles
from ..db.db import get_session
from ..dependencies import get_docking_wrapper
from .. import docking
from ..models import *
from ..util import draw2D
from ..websocket_handler import job_status_stream, notify_job_update

router = APIRouter()


@router.post("/jobs/create", tags=['jobs'], response_model=DockingJobWComp)
def create_job(request: DockingJob, session: Session = Depends(get_session)):
    """
    Create a new docking job.
    """
    session.add(request)
    session.flush()
    if not request.sdf:
        try:
            request.sdf = generate_sdf_from_smiles(request.smiles)
        except ValueError as e:
            raise HTTPException(status_code=400, detail=str(e))

    session.commit()
    session.refresh(request)
    draw2D(request.job_id, request.smiles)

    return request


@router.patch("/jobs/{job_id}/name", tags=['jobs'])
def update_job_name(job_id: str, new_name: str, session: Session = Depends(get_session)):
    """
    Update the name of a docking job.
    """
    job = session.get(DockingJob, job_id)
    if not job:
        return {"error": "Job not found"}
    job.name = new_name
    session.add(job)
    session.commit()
    return {"job_id": job_id, "new_name": new_name}


@router.patch("/jobs/{job_id}/thresholds", tags=['jobs'], response_model=DockingJobWComp)
def update_job_thresholds(
        job_id: str,
        delta_g_threshold: float,
        atom_pair_cst_threshold: float,
        session: Session = Depends(get_session)
):
    """
    Update a job's violation thresholds and re-analyze its existing results.

    No re-docking happens: the best valid complex and the pose RMSD are
    recomputed from the already docked poses.
    """
    job = session.get(DockingJob, job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    job.delta_g_threshold = delta_g_threshold
    job.atom_pair_cst_threshold = atom_pair_cst_threshold

    if job.complexes:
        docking.compute_best_and_rmsd(job, recompute_all=True)

    session.add(job)
    session.commit()
    session.refresh(job)
    return job


@router.delete("/jobs/{job_id}", tags=['jobs'])
def delete_job_by_id(job_id: str, session: Session = Depends(get_session)):
    statement = select(DockingJob).where(DockingJob.job_id == job_id)
    result = session.exec(statement).first()
    if not result:
        return {"error": "Job not found"}

    session.delete(result)
    session.commit()
    notify_job_update(job_id)
    # delete files in poses and previews with the name containing job_id
    poses_path = 'app/static/poses/'
    previews_path = 'app/static/previews/'
    for filename in os.listdir(poses_path):
        if job_id in filename:
            os.remove(os.path.join(poses_path, filename))
    for filename in os.listdir(previews_path):
        if job_id in filename:
            os.remove(os.path.join(previews_path, filename))
    return {"message": "Job deleted successfully"}


@router.post("/jobs/{job_id}/run", tags=['jobs'], status_code=status.HTTP_202_ACCEPTED)
async def run_job(
        job_id: str,
        runs: int,
        background_tasks: BackgroundTasks,
        dw: docking.DockingWrapper = Depends(get_docking_wrapper),
        session: Session = Depends(get_session)
):
    job = session.get(DockingJob, job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # QUEUED until a worker thread picks the job up; run_docking changes it to
    # RUNNING
    job.job_status = JobStatus.QUEUED
    job.error = None
    job.progress_info = "Queued..."
    job.progress = 0
    session.add(job)
    session.commit()
    notify_job_update(job_id)

    background_tasks.add_task(dw.run_docking, job_id, runs)

    return {"queued": True}


@router.post("/jobs/{job_id}/cancel", tags=['jobs'])
def cancel_job(job_id: str, session: Session = Depends(get_session)):
    """
    Request cancellation of a queued or running job.

    The flag is picked up between docking rounds, so a running
    job stops before its next round rather than instantly. A queued job is
    cancelled immediately since it has not started yet.
    """
    job = session.get(DockingJob, job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    if job.job_status not in (JobStatus.QUEUED, JobStatus.RUNNING):
        raise HTTPException(status_code=409, detail="Job is not running")

    docking.request_cancel(job_id)

    if job.job_status == JobStatus.QUEUED:
        # Not started yet -- reflect the cancellation right away. run_docking
        # will still see the flag and skip the job when the worker reaches it.
        job.job_status = JobStatus.CANCELLED
        job.progress_info = "Cancelled."
        session.add(job)
        session.commit()
        notify_job_update(job_id)

    return {"cancelling": True}


@router.websocket("/jobs/status")
async def job_status(websocket: WebSocket):
    await job_status_stream(websocket)





@router.get("/jobs/", tags=['jobs'], response_model=list[DockingJobWComp])
def get_jobs(session: Session = Depends(get_session)):
    """
    Get a list of all docking jobs.
    """
    statement = select(DockingJob)
    results = session.exec(statement).all()
    return results