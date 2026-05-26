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

_COX1_CONSTRAINTS = [
    [354, "CZ",  [-24.005, 54.479, 10.638], 3.70, 0.25],
    [318, "CG2", [-22.027, 53.465,  9.688], 4.30, 0.25],
    [500, "CG",  [-20.234, 49.478,  8.793], 3.90, 0.25],
    [487, "CZ",  [-21.029, 50.916, 12.567], 4.20, 0.25],
]

_EGFR_CONSTRAINTS = [
    [97, "N",   [21.746, -1.973, 54.834], 3.10, 0.25],
    [31, "CG1", [21.545,  0.701, 52.319], 4.40, 0.25],
    [23, "CD1", [19.960, -1.225, 53.519], 3.60, 0.25],
]

_ESTROGEN_CONSTRAINTS = [
    [81, "O",   [-1.263, -24.794, 28.161], 3.40, 0.25],
    [98, "CD1", [-2.981, -27.861, 27.237], 3.40, 0.25],
    [98, "CE1", [-0.673, -30.794, 27.049], 4.30, 0.25],
]

_SERT_CONSTRAINTS = [
    [365, "CB", [32.410, 186.874, 144.154], 3.90, 0.25],
    [99,  "CB", [30.547, 184.819, 143.960], 4.30, 0.25],
]

AVAILABLE_TARGETS: dict[str, dict[str, Any]] = {
    "7f83": {
        "id": "7f83",
        "name": "GHSR",
        "full_name": "Ghrelin Receptor (GHSR-1a)",
        "pdb_code": "7F83",
        "description": "G protein-coupled receptor regulating appetite and growth hormone secretion. Target for obesity and metabolic disorder therapies.",
        "pose_path": "input/targets/Ghrelin_Receptor/Final/7F83_A_relax_noligand.pdb",
        "preview_path": "input/targets/Ghrelin_Receptor/Final/7F83_A_relax_noligand.pdb",
        "core_ligand_path": "input/targets/Ghrelin_Receptor/Final/7f83_X_LIG_clean_core.sdf",
        "constraints": _7F83_CONSTRAINTS,
    },
    "cox1": {
        "id": "cox1",
        "name": "COX-1",
        "full_name": "Cyclooxygenase-1",
        "pdb_code": "3N8Z",
        "description": "Enzyme responsible for prostaglandin synthesis. Primary target for anti-inflammatory drugs (NSAIDs).",
        "pose_path": "input/targets/COX_1/Final/3N8Z_A_relax_noligand.pdb",
        "preview_path": "input/targets/COX_1/Final/3N8Z_A_relax_noligand.pdb",
        "core_ligand_path": "input/targets/COX_1/Final/3n8z_L_FLP_clean_core.sdf",
        "constraints": _COX1_CONSTRAINTS,
    },
    "egfr": {
        "id": "egfr",
        "name": "EGFR",
        "full_name": "Epidermal Growth Factor Receptor",
        "pdb_code": "1M17",
        "description": "Receptor tyrosine kinase mutated in ~15% of non-small cell lung cancers. Major oncology target.",
        "pose_path": "input/targets/EGFR/Final/1M17_A_relax_noligand.pdb",
        "preview_path": "input/targets/EGFR/Final/1M17_A_relax_noligand.pdb",
        "core_ligand_path": "input/targets/EGFR/Final/1m17_B_AQ4_clean_core.sdf",
        "constraints": _EGFR_CONSTRAINTS,
    },
    "estrogen": {
        "id": "estrogen",
        "name": "ERα",
        "full_name": "Estrogen Receptor α",
        "pdb_code": "7KBS",
        "description": "Nuclear receptor mediating estrogen signaling. Central target for hormone-sensitive breast cancer therapy.",
        "pose_path": "input/targets/Estrogen_Receptor/Final/7KBS_A_relax_noligand.pdb",
        "preview_path": "input/targets/Estrogen_Receptor/Final/7KBS_A_relax_noligand.pdb",
        "core_ligand_path": "input/targets/Estrogen_Receptor/Final/7kbs_C_RAL_clean_core.sdf",
        "constraints": _ESTROGEN_CONSTRAINTS,
    },
    "sert": {
        "id": "sert",
        "name": "SERT",
        "full_name": "Serotonin Transporter",
        "pdb_code": "6AWP",
        "description": "Membrane transporter regulating serotonin reuptake. The primary target for SSRI antidepressants.",
        "pose_path": "input/targets/Serotonin_Transporter/Final/6AWP_A_relax_noligand.pdb",
        "preview_path": "input/targets/Serotonin_Transporter/Final/6AWP_A_relax_noligand.pdb",
        "core_ligand_path": "input/targets/Serotonin_Transporter/Final/6awp_D_FVX_clean_core.sdf",
        "constraints": _SERT_CONSTRAINTS,
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