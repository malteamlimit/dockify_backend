import datetime
import io
import json
import os
import shutil
import zipfile

from fastapi import APIRouter, HTTPException, UploadFile, File
from sqlmodel import Session, select
from starlette.responses import StreamingResponse

from app.db.db import init_db, recreate_engine, engine
from app.models import DockingJob
from app.util import draw2D

router = APIRouter()

DB_PATH = "data/dockify.db"
POSES_DIR = "app/static/poses"
POSE_EXTS = (".pdb",)
ARCHIVE_VERSION = 1


@router.get("/database/export", tags=['database'])
def export_db():
    if not os.path.exists(DB_PATH):
        raise HTTPException(status_code=404, detail="Database file not found")

    timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, mode='w', compression=zipfile.ZIP_DEFLATED) as zf:
        zf.write(DB_PATH, arcname="database.db")

        pose_count = 0
        if os.path.isdir(POSES_DIR):
            for name in sorted(os.listdir(POSES_DIR)):
                full = os.path.join(POSES_DIR, name)
                if os.path.isfile(full) and name.endswith(POSE_EXTS):
                    zf.write(full, arcname=f"poses/{name}")
                    pose_count += 1

        manifest = {
            "version": ARCHIVE_VERSION,
            "created_at": datetime.datetime.now().isoformat(),
            "pose_files": pose_count,
        }
        zf.writestr("manifest.json", json.dumps(manifest, indent=2))

    size = buf.tell()
    buf.seek(0)
    return StreamingResponse(
        buf,
        media_type="application/zip",
        headers={
            "Content-Disposition": f'attachment; filename="dockify_backup_{timestamp}.zip"',
            "Content-Length": str(size),
        },
    )


@router.post("/database/import", tags=['database'])
async def import_db(file: UploadFile = File(...)):
    if not file.filename.endswith('.zip'):
        raise HTTPException(status_code=400, detail="Only .zip backup files are allowed")

    timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")

    data = await file.read()
    await file.close()

    try:
        zf = zipfile.ZipFile(io.BytesIO(data))
    except zipfile.BadZipFile:
        raise HTTPException(status_code=400, detail="Uploaded file is not a valid zip archive")

    names = set(zf.namelist())
    if "database.db" not in names:
        zf.close()
        raise HTTPException(status_code=400, detail="Backup archive is missing database.db")

    os.makedirs(os.path.dirname(DB_PATH), exist_ok=True)
    os.makedirs(POSES_DIR, exist_ok=True)

    db_backup = None
    poses_backup = None
    if os.path.exists(DB_PATH):
        db_backup = f"data/dockify_backup_{timestamp}.db"
        shutil.copy2(DB_PATH, db_backup)
    if os.listdir(POSES_DIR):
        poses_backup = f"app/static/poses_backup_{timestamp}"
        shutil.copytree(POSES_DIR, poses_backup)

    try:
        with zf.open("database.db") as src, open(DB_PATH, "wb") as dst:
            shutil.copyfileobj(src, dst)

        for name in os.listdir(POSES_DIR):
            full = os.path.join(POSES_DIR, name)
            if os.path.isfile(full):
                os.remove(full)

        for member in zf.namelist():
            if not member.startswith("poses/") or member.endswith("/"):
                continue
            base = os.path.basename(member)
            if not base or not base.endswith(POSE_EXTS):
                continue
            with zf.open(member) as src, open(os.path.join(POSES_DIR, base), "wb") as dst:
                shutil.copyfileobj(src, dst)

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
            "backup": db_backup,
            "poses_backup": poses_backup,
        }
    except (OSError, IOError, zipfile.BadZipFile) as e:
        if db_backup and os.path.exists(db_backup):
            shutil.move(db_backup, DB_PATH)
        if poses_backup and os.path.isdir(poses_backup):
            shutil.rmtree(POSES_DIR, ignore_errors=True)
            shutil.move(poses_backup, POSES_DIR)
        raise HTTPException(status_code=500, detail=f"Failed to import database: {str(e)}")
    finally:
        zf.close()


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