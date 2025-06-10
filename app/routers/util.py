from fastapi import APIRouter, Depends
from openbabel import pybel
from pydantic import BaseModel
from sqlmodel import Session
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, QED

from app.db.db import get_session
from app.models import DockingJob

router = APIRouter()

class ConfBase(BaseModel):
    smiles: str
    job_id: str

@router.post("/util/genConf", tags=['util'])
def generate_conformer(conf: ConfBase, session: Session = Depends(get_session)):
    if conf.smiles == "":
        return {"sdf": ""}
    mol_rdkit = Chem.MolFromSmiles(conf.smiles)
    mol_rdkit = Chem.AddHs(mol_rdkit)
    AllChem.EmbedMolecule(mol_rdkit)
    sdf = Chem.MolToMolBlock(mol_rdkit)

    job = session.get(DockingJob, conf.job_id)
    if job:
        job.sdf = sdf
        session.add(job)
        session.commit()

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