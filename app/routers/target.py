import os

from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import FileResponse
from sqlmodel import Session

from app.db.db import get_session
from app.models import TargetConfig
from app.targets import AVAILABLE_TARGETS, get_target_public_info

router = APIRouter()


@router.get("/targets", tags=["target"])
def list_targets():
    """Return summary info for all available targets (for the selection dialog)."""
    return [
        {
            "id": t["id"],
            "name": t["name"],
            "full_name": t["full_name"],
            "pdb_code": t["pdb_code"],
            "description": t["description"],
        }
        for t in AVAILABLE_TARGETS.values()
    ]


@router.get("/target", tags=["target"])
def get_active_target(session: Session = Depends(get_session)):
    """Return the active target with core_smiles and constraints, or 404 if none selected."""
    config = session.get(TargetConfig, 1)
    if not config:
        raise HTTPException(status_code=404, detail="No target selected")
    info = get_target_public_info(config.target_id)
    if not info:
        raise HTTPException(status_code=404, detail="Target not found")
    return info


@router.post("/target/{target_id}", tags=["target"])
def select_target(target_id: str, session: Session = Depends(get_session)):
    """Select a target by ID. Upserts the singleton TargetConfig row."""
    if target_id not in AVAILABLE_TARGETS:
        raise HTTPException(status_code=400, detail=f"Unknown target: {target_id}")

    config = session.get(TargetConfig, 1)
    if config:
        config.target_id = target_id
    else:
        config = TargetConfig(id=1, target_id=target_id)

    session.add(config)
    session.commit()
    return get_target_public_info(target_id)


@router.get("/target/{target_id}/preview", tags=["target"])
def get_target_preview(target_id: str):
    """Serve the target's relaxed PDB file for the 3D card preview."""
    target = AVAILABLE_TARGETS.get(target_id)
    if not target:
        raise HTTPException(status_code=404, detail="Target not found")
    pose_path = target["pose_path"]
    if not os.path.exists(pose_path):
        raise HTTPException(status_code=404, detail="Preview file not found")
    return FileResponse(pose_path, media_type="text/plain", filename=f"{target_id}.pdb")