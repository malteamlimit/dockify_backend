from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D


def parse_ligand_atoms(pdb_path: str, resname: str = "UNK", heavy_only: bool = True) -> dict[str, tuple[float, float, float]]:
    """Parse the ligand atoms of a docked pose PDB file, keyed by atom name."""
    coords: dict[str, tuple[float, float, float]] = {}
    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if line[17:20].strip() != resname:
                continue
            atom_name = line[12:16].strip()
            element = line[76:78].strip() if len(line) >= 78 else ""
            if heavy_only and (element == "H" or (not element and atom_name.startswith("H"))):
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords[atom_name] = (x, y, z)
    return coords


def rmsd_from_pdb(pdb_a: str, pdb_b: str, resname: str = "UNK", heavy_only: bool = True) -> float:
    """Compute ligand-atom RMSD between two docked poses without superposition.
    Receptor frame is fixed across poses, so atoms can be compared directly by name.
    """
    a = parse_ligand_atoms(pdb_a, resname, heavy_only)
    b = parse_ligand_atoms(pdb_b, resname, heavy_only)
    common = a.keys() & b.keys()
    if not common:
        raise ValueError(f"No common ligand atoms between {pdb_a} and {pdb_b}")
    sq = 0.0
    for name in common:
        ax, ay, az = a[name]
        bx, by, bz = b[name]
        sq += (ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2
    return (sq / len(common)) ** 0.5


def draw2D(job_id, smiles):
    """Generate a 2D depiction of the ligand from its SMILES string, matching the reference core structure if possible."""
    pybel_mol = next(pybel.readfile("sdf", 'input/ref_ligand_core.sdf'))
    pybel_sdf = pybel_mol.write('sdf')
    template = Chem.MolFromMolBlock(pybel_sdf)
    query = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(template)
    if query.HasSubstructMatch(template):
        AllChem.GenerateDepictionMatching2DStructure(query, template)

    drawer = rdMolDraw2D.MolDraw2DSVG(1000, 1000)
    drawer.SetLineWidth(5)
    drawer.drawOptions().clearBackground = False
    drawer.DrawMolecule(rdMolDraw2D.PrepareMolForDrawing(query))
    drawer.FinishDrawing()
    with open('app/static/previews/' + job_id + '.svg', 'w') as f:
        f.write(drawer.GetDrawingText())
