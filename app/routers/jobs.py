import asyncio
import os

from fastapi import APIRouter, BackgroundTasks, Depends, WebSocket, status
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from sqlmodel import Session, select

from ..db.db import get_session
from ..dependencies import get_docking_wrapper
from .. import docking
from ..models import *
from ..websocket_handler import get_job_status

router = APIRouter()


@router.post("/jobs/create", tags=['jobs'])
def create_job(request: DockingJob, session: Session = Depends(get_session)):
    """
    Create a new docking job.
    """
    session.add(request)
    session.commit()
    drawer = rdMolDraw2D.MolDraw2DSVG(1000, 1000)
    drawer.DrawMolecule(rdMolDraw2D.PrepareMolForDrawing(Chem.MolFromSmiles(request.smiles)))
    drawer.FinishDrawing()
    with open('app/static/previews/' + request.job_id + '.svg', 'w') as f:
        f.write(drawer.GetDrawingText())
    return {request.model_dump_json()}


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


@router.delete("/jobs/{job_id}", tags=['jobs'])
def delete_job_by_id(job_id: str, session: Session = Depends(get_session)):
    statement = select(DockingJob).where(DockingJob.job_id == job_id)
    result = session.exec(statement).first()
    if not result:
        return {"error": "Job not found"}

    session.delete(result)
    session.commit()
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
        return False

    job.job_status = JobStatus.RUNNING
    job.error = None
    job.progress_info = "Preparing Ligand..."
    job.progress = 0
    session.add(job)
    session.commit()

    background_tasks.add_task(dw.run_docking, job_id, runs)
    # asyncio.create_task(dw.run_docking(job_id, runs))

    return True


@router.websocket("/jobs/{job_id}/status")
async def check_job(websocket: WebSocket, job_id: str):
    await get_job_status(websocket, job_id)





@router.get("/jobs/", tags=['jobs'], response_model=list[DockingJobWComp])
def get_jobs(session: Session = Depends(get_session)):
    """
    Get a list of all docking jobs.
    """
    statement = select(DockingJob)
    results = session.exec(statement).all()
    return results