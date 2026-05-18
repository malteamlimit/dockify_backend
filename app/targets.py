"""Definitions of all available docking targets."""
from typing import Any


def _smiles_from_sdf(path: str) -> str:
    """Derive canonical SMILES from an SDF file using the same pybel→RDKit pipeline as docking."""
    try:
        from openbabel import pybel
        from rdkit import Chem
        pybel_mol = next(pybel.readfile("sdf", path))
        pybel_sdf = pybel_mol.write("sdf")
        mol = Chem.MolFromMolBlock(pybel_sdf, removeHs=True)
        if mol is None:
            return ""
        return Chem.MolToSmiles(mol)
    except Exception as exc:
        print(f"[targets] Could not extract SMILES from {path}: {exc}")
        return ""


_7F83_CONSTRAINTS = [
    [364, "HG",  [-6.752,  -0.1555, 13.0855], 1.8,  0.125],
    [65,  "OD2", [-7.1638,  5.8368, 16.5862], 3.23, 0.25],
    [65,  "OD2", [-7.5181,  3.1143, 15.5623], 3.25, 0.25],
    [89,  "CB",  [-6.0966,  5.3594, 15.7673], 3.7,  0.25],
    [86,  "CD",  [-7.1638,  5.8368, 16.5862], 5.11, 0.5],
]

AVAILABLE_TARGETS: dict[str, dict[str, Any]] = {
    "7f83": {
        "id": "7f83",
        "name": "GHSR",
        "full_name": "Ghrelin Receptor (GHSR-1a)",
        "pdb_code": "7F83",
        "description": "G protein-coupled receptor regulating appetite and growth hormone secretion. Target for obesity and metabolic disorder therapies.",
        "pose_path": "input/new-targets/Ghrelin_Receptor/Final/7F83_A_relax.pdb",
        "core_ligand_path": "input/new-targets/Ghrelin_Receptor/Final/7f83_X_LIG_clean_core.sdf",
        "constraints": _7F83_CONSTRAINTS,
    },
    "cox1": {
        "id": "cox1",
        "name": "COX-1",
        "full_name": "Cyclooxygenase-1",
        "pdb_code": "3N8Z",
        "description": "Enzyme responsible for prostaglandin synthesis. Primary target for anti-inflammatory drugs (NSAIDs).",
        "pose_path": "input/new-targets/COX_1/Final/3N8Z_A_relax.pdb",
        "core_ligand_path": "input/new-targets/COX_1/Final/3n8z_L_FLP_clean_core.sdf",
        "constraints": None,
    },
    "egfr": {
        "id": "egfr",
        "name": "EGFR",
        "full_name": "Epidermal Growth Factor Receptor",
        "pdb_code": "1M17",
        "description": "Receptor tyrosine kinase mutated in ~15% of non-small cell lung cancers. Major oncology target.",
        "pose_path": "input/new-targets/EGFR/Final/1M17_A_relax.pdb",
        "core_ligand_path": "input/new-targets/EGFR/Final/1m17_B_AQ4_clean_core.sdf",
        "constraints": None,
    },
    "estrogen": {
        "id": "estrogen",
        "name": "ERα",
        "full_name": "Estrogen Receptor α",
        "pdb_code": "7KBS",
        "description": "Nuclear receptor mediating estrogen signaling. Central target for hormone-sensitive breast cancer therapy.",
        "pose_path": "input/new-targets/Estrogen_Receptor/Final/7KBS_A_relax.pdb",
        "core_ligand_path": "input/new-targets/Estrogen_Receptor/Final/7kbs_C_RAL_clean_core.sdf",
        "constraints": None,
    },
    "sert": {
        "id": "sert",
        "name": "SERT",
        "full_name": "Serotonin Transporter",
        "pdb_code": "6AWP",
        "description": "Membrane transporter regulating serotonin reuptake. The primary target for SSRI antidepressants.",
        "pose_path": "input/new-targets/Serotonin_Transporter/Final/6AWP_A_relax.pdb",
        "core_ligand_path": "input/new-targets/Serotonin_Transporter/Final/6awp_D_FVX_clean_core.sdf",
        "constraints": None,
    },
}


def get_target_public_info(target_id: str) -> dict | None:
    """Return public target info including dynamically derived core_smiles."""
    target = AVAILABLE_TARGETS.get(target_id)
    if not target:
        return None
    return {
        "id": target["id"],
        "name": target["name"],
        "full_name": target["full_name"],
        "pdb_code": target["pdb_code"],
        "description": target["description"],
        "core_smiles": _smiles_from_sdf(target["core_ligand_path"]),
        "constraints": target["constraints"],
    }