from fastapi import APIRouter, Depends, HTTPException
from openbabel import pybel
from pydantic import BaseModel
from sqlmodel import Session
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, QED

from app.db.db import get_session
from app.models import DockingJob
from app.util import draw2D

router = APIRouter()

class ConfBase(BaseModel):
    smiles: str
    job_id: str


def generate_sdf_from_smiles(smiles: str) -> str:
    """Utility: SMILES → SDF"""
    if not smiles:
        return ""

    mol_rdkit = Chem.MolFromSmiles(smiles)
    if not mol_rdkit:
        raise ValueError(f"Invalid SMILES: {smiles}")
    mol_rdkit = Chem.AddHs(mol_rdkit)
    result = AllChem.EmbedMolecule(mol_rdkit)

    if result != 0:
        raise ValueError("Failed to generate 3D conformer")

    return Chem.MolToMolBlock(mol_rdkit)

@router.post("/util/genConf", tags=['util'])
def generate_conformer(conf: ConfBase, session: Session = Depends(get_session)):
    """Generate conformer and update job"""
    try:
        sdf = generate_sdf_from_smiles(conf.smiles)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

    job = session.get(DockingJob, conf.job_id)
    if job:
        job.smiles = conf.smiles
        job.sdf = sdf

        session.commit()
        draw2D(job.job_id, job.smiles)

    return {"sdf": sdf}

@router.post("/util/props", tags=['util'])
def generate_props(conf: ConfBase, session: Session = Depends(get_session)):
    mol_rdkit = Chem.MolFromSmiles(conf.smiles)

    pybel_mol = next(pybel.readfile("sdf", 'input/ref_ligand_core.sdf'))
    pybel_sdf = pybel_mol.write('sdf')
    ref_mol = Chem.MolFromMolBlock(pybel_sdf)

    is_sub = mol_rdkit.HasSubstructMatch(ref_mol)

    weight = Descriptors.MolWt(mol_rdkit)
    hbond_acc = Descriptors.NOCount(mol_rdkit)
    hbond_don = Descriptors.NHOHCount(mol_rdkit)
    logp = Descriptors.MolLogP(mol_rdkit)
    qed = QED.qed(mol_rdkit)

    job = session.get(DockingJob, conf.job_id)
    if job:
        job.weight = weight
        job.hbond_acc = hbond_acc
        job.hbond_don = hbond_don
        job.logp = logp
        job.qed = qed
        job.is_sub = is_sub
        session.add(job)
        session.commit()

    return {"is_sub": is_sub, "weight": weight, "hbond_acc": hbond_acc, "hbond_don": hbond_don, "logp": logp, "qed": qed}