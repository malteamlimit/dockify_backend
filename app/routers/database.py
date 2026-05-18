import datetime
import io
import json
import os
import shutil
import sqlite3
import tempfile
import zipfile

from fastapi import APIRouter, HTTPException, UploadFile, File
from sqlmodel import Session, select
from starlette.responses import StreamingResponse

from app.db.db import init_db, recreate_engine, engine
from app.docking import clear_all_cancels
from app.models import DockingJob, TargetConfig
from app.targets import AVAILABLE_TARGETS
from app.util import draw2D

router = APIRouter()

DB_PATH = "data/dockify.db"


def _extract_zip_db_to_tempfile(zf: zipfile.ZipFile) -> str:
    """Write the ZIP's database.db to a temp file and return its path."""
    with zf.open("database.db") as src:
        db_bytes = src.read()
    with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as tmp:
        tmp.write(db_bytes)
        return tmp.name


def _check_db_integrity(tmp_path: str) -> None:
    """Raise HTTPException 400 if the SQLite file is malformed."""
    conn = sqlite3.connect(tmp_path)
    try:
        result = conn.execute("PRAGMA integrity_check").fetchone()
        if result is None or result[0] != "ok":
            raise HTTPException(
                status_code=400,
                detail=f"The backup database is corrupted and cannot be imported (integrity_check: {result}).",
            )
    finally:
        conn.close()


def _read_target_from_zip(zf: zipfile.ZipFile) -> str | None:
    """Return the target_id stored in the ZIP's database.db, or None if absent."""
    try:
        tmp_path = _extract_zip_db_to_tempfile(zf)
        try:
            conn = sqlite3.connect(tmp_path)
            try:
                row = conn.execute(
                    "SELECT target_id FROM targetconfig WHERE id = 1"
                ).fetchone()
                return str(row[0]) if row else None
            except sqlite3.OperationalError:
                return None  # Old backup without targetconfig table
            finally:
                conn.close()
        finally:
            os.unlink(tmp_path)
    except HTTPException:
        raise
    except Exception:
        return None
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

    # Validate the archive DB before touching the live DB.
    try:
        tmp_path = _extract_zip_db_to_tempfile(zf)
        try:
            _check_db_integrity(tmp_path)
        finally:
            os.unlink(tmp_path)
    except HTTPException:
        zf.close()
        raise

    # Check whether the backup's target matches the currently active target.
    incoming_target_id = _read_target_from_zip(zf)
    # Normalise to str and reject unknown/legacy values (e.g. old integer "0").
    incoming_target_id = str(incoming_target_id) if incoming_target_id is not None else None
    if incoming_target_id and incoming_target_id not in AVAILABLE_TARGETS:
        incoming_target_id = None  # treat unrecognised legacy ID as "no target"

    with Session(engine) as check_session:
        current_config = check_session.get(TargetConfig, 1)
    current_target_id = current_config.target_id if current_config else None

    if incoming_target_id is not None and current_target_id is not None and incoming_target_id != current_target_id:
        zf.close()
        current_name = AVAILABLE_TARGETS.get(current_target_id, {}).get("name", current_target_id)
        incoming_name = AVAILABLE_TARGETS.get(incoming_target_id, {}).get("name", incoming_target_id)
        raise HTTPException(
            status_code=409,
            detail=(
                f"This backup was created for a different target ({incoming_name} · {incoming_target_id}), "
                f"but the current target is {current_name} · {current_target_id}. "
                f"Use 'Reset Database' to switch targets before importing this backup."
            ),
        )

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
            # Remove any targetconfig row with an unrecognised target_id
            # (e.g. legacy integer "0") so the frontend shows the selection dialog.
            tc = session.get(TargetConfig, 1)
            if tc and tc.target_id not in AVAILABLE_TARGETS:
                session.delete(tc)
                session.commit()

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

    old_db_exists = os.path.exists(db_path)

    try:
        if old_db_exists:
            os.remove(db_path)

        # recreate db and clear any in-memory state tied to the old db
        clear_all_cancels()
        recreate_engine()
        init_db()

        # clear poses (app/static/poses) and previews (app/static/previews)
        for directory in ("app/static/poses", "app/static/previews"):
            if os.path.isdir(directory):
                for name in os.listdir(directory):
                    full = os.path.join(directory, name)
                    if os.path.isfile(full):
                        os.remove(full)

        return {"message": "Database reset successfully"}

    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Database reset failed: {str(e)}"
        )