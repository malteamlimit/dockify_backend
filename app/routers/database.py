import datetime
import os
import shutil

from fastapi import APIRouter, HTTPException, UploadFile, File
from starlette.responses import FileResponse

from app.db.db import init_db, recreate_engine

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
        backup_path = f"data/dockify_{timestamp}.backup"
        shutil.copy2(db_path, backup_path)

    try:
        # Stream the uploaded file to disk in chunks to avoid large memory usage
        with open(db_path, "wb") as buffer:
            while True:
                chunk = await file.read(1024 * 1024)
                if not chunk:
                    break
                buffer.write(chunk)

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
    if os.path.exists(db_path):
        backup_path = f"data/dockify_{timestamp}.backup"
        shutil.copy2(db_path, backup_path)
        os.remove(db_path)

    recreate_engine()
    init_db()

    return {
        "message": "Database reset successfully",
        "backup": backup_path,
    }

