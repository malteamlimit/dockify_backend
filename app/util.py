from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D


def draw2D(job_id, smiles):
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
