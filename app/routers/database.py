import datetime
import os
import shutil

from fastapi import APIRouter, HTTPException, UploadFile, File
from sqlmodel import Session, select
from starlette.responses import FileResponse

from app.db.db import init_db, recreate_engine, engine
from app.models import DockingJob
from app.util import draw2D

router = APIRouter()

@router.get("/database/export", tags=['database'])
def export_db():
    db_path = "data/dockify.db"

    if not os.path.exists(db_path):
        raise HTTPException(status_code=404, detail="Database file not found")

    return FileResponse(
        path=db_path,
        filename="database_export.db",
        media_type="application/octet-stream"
    )


@router.post("/database/import", tags=['database'])
async def import_db(file: UploadFile = File(...)):
    db_path = "data/dockify.db"
    timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")

    if not file.filename.endswith('.db'):
        raise HTTPException(status_code=400, detail="Only .db files are allowed")

    os.makedirs(os.path.dirname(db_path), exist_ok=True)

    backup_path = None
    if os.path.exists(db_path):
        backup_path = f"data/dockify_backup_{timestamp}.db"
        shutil.copy2(db_path, backup_path)

    try:
        # Stream the uploaded file to disk in chunks to avoid large memory usage
        with open(db_path, "wb") as buffer:
            while True:
                chunk = await file.read(1024 * 1024)
                if not chunk:
                    break
                buffer.write(chunk)

        # recreate all 2d preview images
        recreate_engine()
        with Session(engine) as session:
            jobs = session.exec(select(DockingJob)).all()
            for job in jobs:
                try:
                    draw2D(job.job_id, job.smiles)
                except Exception as e:
                    print(f"Failed to generate preview for job {job.job_id}: {e}")

        return {
            "message": "Database imported successfully",
            "filename": file.filename,
            "backup": backup_path,
        }
    except (OSError, IOError) as e:
        if backup_path and os.path.exists(backup_path):
            shutil.move(backup_path, db_path)
        raise HTTPException(status_code=500, detail=f"Failed to import database: {str(e)}")
    finally:
        await file.close()


@router.post("/database/reset", tags=['database'])
async def reset_db():
    db_path = "data/dockify.db"
    timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")

    os.makedirs(os.path.dirname(db_path), exist_ok=True)

    backup_path = None
    old_db_exists = os.path.exists(db_path)

    try:
        # backup existing db
        if old_db_exists:
            backup_path = f"data/dockify_backup_{timestamp}.db"
            shutil.copy2(db_path, backup_path)
            os.rename(db_path, f"{db_path}.old")

        # recreate db
        recreate_engine()
        init_db()

        # remove old db
        if old_db_exists:
            os.remove(f"{db_path}.old")

        return {
            "message": "Database reset successfully",
            "backup": backup_path,
        }

    except Exception as e:
        # rollback old db
        if old_db_exists and os.path.exists(f"{db_path}.old"):
            if os.path.exists(db_path):
                os.remove(db_path)
            os.rename(f"{db_path}.old", db_path)

        raise HTTPException(
            status_code=500,
            detail=f"Database reset failed: {str(e)}"
        )