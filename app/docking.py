import asyncio
import pickle

from pyrosetta.rosetta.core.scoring import residue_rmsd_nosuper
from rdkit import Chem
from rdkit.Chem import Lipinski, AllChem, Descriptors, QED
from rdkit.Chem.Draw import rdMolDraw2D
from sqlmodel import Session

from openbabel import pybel

from .db.db import engine
from .models import DockingJob, JobStatus, ComplexResult
from .websocket_handler import job_update_queues

class DockingWrapper:
    def __init__(self, loop=None):
        import pyrosetta
        self.pyrosetta = pyrosetta
        self.pyrosetta.init()

        self.loop = loop



    def rdkit_to_mutable_res(self, mol):
        chem_manager = self.pyrosetta.rosetta.core.chemical.ChemicalManager.get_instance()

        tag = "fa_standard"
        lig_restype = self.pyrosetta.rosetta.core.chemical.MutableResidueType(
            chem_manager.atom_type_set(tag),
            chem_manager.element_set("default"),
            chem_manager.mm_atom_type_set(tag),
            chem_manager.orbital_type_set(tag)
        )
        lig_restype.name("UNKNOWN")
        lig_restype.name3("UNK")
        lig_restype.name1("X")
        lig_restype.interchangeability_group("UNK")

        index_to_vd = {}

        conf = mol.GetConformer(0)

        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            element_name = atom.GetSymbol()
            charge = atom.GetFormalCharge()

            vd_atom = lig_restype.add_atom("")
            restype_atom = lig_restype.atom(vd_atom)
            restype_atom.element_type(lig_restype.element_set().element(element_name))
            restype_atom.formal_charge(charge)
            restype_atom.mm_name("VIRT")

            atom_pos = conf.GetAtomPosition(i)
            xyz = self.pyrosetta.rosetta.numeric.xyzVector_double_t(atom_pos.x, atom_pos.y, atom_pos.z)
            restype_atom.ideal_xyz(xyz)

            index_to_vd[i] = vd_atom

        for bond in mol.GetBonds():
            bond_name = bond.GetBondType()
            if bond_name == Chem.rdchem.BondType.SINGLE:
                bond_name = self.pyrosetta.rosetta.core.chemical.BondName.SingleBond
            elif bond_name == Chem.rdchem.BondType.DOUBLE:
                bond_name = self.pyrosetta.rosetta.core.chemical.BondName.DoubleBond
            elif bond_name == Chem.rdchem.BondType.TRIPLE:
                bond_name = self.pyrosetta.rosetta.core.chemical.BondName.TripleBond
            elif bond_name == Chem.rdchem.BondType.AROMATIC:
                bond_name = self.pyrosetta.rosetta.core.chemical.BondName.AromaticBond
            else:
                print("ERROR: encountered unknown bond type", bond_name)
                bond_name = self.pyrosetta.rosetta.core.chemical.BondName.UnknownBond

            lig_restype.add_bond(
                index_to_vd[bond.GetBeginAtom().GetIdx()],
                index_to_vd[bond.GetEndAtom().GetIdx()],
                bond_name
            )

        self.pyrosetta.rosetta.core.chemical.rename_atoms(lig_restype, True)
        self.pyrosetta.rosetta.core.chemical.rosetta_retype_fullatom(lig_restype, True)
        self.pyrosetta.rosetta.core.chemical.rosetta_recharge_fullatom(lig_restype)

        self.pyrosetta.rosetta.core.chemical.find_bonds_in_rings(lig_restype)

        nbr_vd = 0
        shortest_nbr_dist = 999999.99
        for vd in index_to_vd.values():
            if lig_restype.atom(vd).element_type().get_chemical_symbol() == "H":
                continue
            tmp_dist = self.pyrosetta.rosetta.core.chemical.find_nbr_dist(lig_restype, vd)
            if tmp_dist < shortest_nbr_dist:
                shortest_nbr_dist = tmp_dist
                nbr_vd = vd

        lig_restype.nbr_radius(shortest_nbr_dist)
        lig_restype.nbr_atom(nbr_vd)
        lig_restype.assign_internal_coordinates()
        lig_restype.autodetermine_chi_bonds()

        index_to_name = {}
        pos_to_name = {}

        for idx in range(mol.GetNumAtoms()):
            vd = index_to_vd[idx]
            atm_name = lig_restype.atom_name(vd)
            atom_pos = conf.GetAtomPosition(idx)
            binned_pos = (int(atom_pos.x * 10), int(atom_pos.y * 10), int(atom_pos.z * 10))
            index_to_name[idx] = atm_name
            pos_to_name[binned_pos] = atm_name

        return lig_restype, index_to_vd, index_to_name, pos_to_name


    def add_confs_to_res(self, mol, mutable_res, index_to_vd):
        rotamers_spec = self.pyrosetta.rosetta.core.chemical.rotamers.StoredRotamerLibrarySpecification()

        for i in range(mol.GetNumConformers()):
            conf = mol.GetConformer(i)
            single_conf_spec = self.pyrosetta.rosetta.std.map_std_string_numeric_xyzVector_double_t_std_allocator_std_pair_const_std_string_numeric_xyzVector_double_t()
            for idx, atm_vd in index_to_vd.items():
                rdkit_atm_pos = conf.GetAtomPosition(idx)
                single_conf_spec[mutable_res.atom_name(atm_vd)] = self.pyrosetta.rosetta.numeric.xyzVector_double_t(
                    rdkit_atm_pos.x,
                    rdkit_atm_pos.y,
                    rdkit_atm_pos.z)

            rotamers_spec.add_rotamer(single_conf_spec)

        mutable_res.rotamer_library_specification(rotamers_spec)
        return mutable_res


    def mutable_res_to_res(self, mutable_res):
        lig_restype_non_mutable = self.pyrosetta.rosetta.core.chemical.ResidueType.make(mutable_res)
        return self.pyrosetta.rosetta.core.conformation.Residue(lig_restype_non_mutable, True)


    def prepare_pose(self, smiles: str, pose, ref_mol):
        # Generate RDKit Mol from SMILES
        new_mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if new_mol is None:
            raise ValueError("Invalid SMILES string provided. (RDKit failed to generate mol object)")
        else:
            try:
                Chem.SanitizeMol(new_mol)
            except:
                raise ValueError("Invalid SMILES string provided. (Sanitization failed)")

        # Determine number of different conformers based on molecule flexibility
        nrotbonds = Lipinski.NumRotatableBonds(new_mol)
        nconf = 50
        if nrotbonds > 12:
            nconf = 300
        elif nrotbonds >= 8:
            nconf = 200
        new_mol = Chem.AddHs(new_mol)

        # Align new Mol with reference Mol and generate conformers
        if not new_mol.HasSubstructMatch(ref_mol):
            print('ERROR: New mol has no substructure match with reference mol')
            return
        for i in range(nconf):
            AllChem.ConstrainedEmbed(new_mol, ref_mol, clearConfs=False, coreConfId=0, randomseed=42 + i)
        AllChem.MMFFOptimizeMoleculeConfs(new_mol, numThreads=0)
        match = new_mol.GetSubstructMatch(ref_mol)
        atom_map = list(zip(match, range(ref_mol.GetNumAtoms())))

        for cid in range(new_mol.GetNumConformers()):
            AllChem.AlignMol(new_mol, ref_mol, atomMap=atom_map, prbCid=cid)

        # with Chem.SDWriter('constraints/rdkit_conf.sdf') as w:
        #    for cid in range(new_mol.GetNumConformers()):
        #        w.write(new_mol, confId=cid)

        # # Align new Mol with reference Mol and generate conformers
        # AllChem.ConstrainedEmbed(new_mol, ref_mol)
        # AllChem.EmbedMultipleConfs(new_mol, numConfs=nconf, clearConfs=False, maxAttempts=30)
        # AllChem.AlignMolConformers(new_mol)

        # Turn RDKit Mol into Rosetta Residue
        mutres, index_to_vd, index_to_name, pos_to_name = self.rdkit_to_mutable_res(new_mol)
        mutres = self.add_confs_to_res(new_mol, mutres, index_to_vd)
        res = self.mutable_res_to_res(mutres)

        # # Examplified creation of 2D image with Rosetta atom names
        # d = rdMolDraw2D.MolDraw2DCairo(750, 750)  # or MolDraw2DSVG to get SVGs
        # for i in range(new_mol.GetNumAtoms()):
        #     new_mol.GetAtomWithIdx(i).SetProp('atomNote', index_to_name[i])
        # AllChem.Compute2DCoords(new_mol)
        # d.DrawMolecule(new_mol)
        # d.FinishDrawing()
        # d.WriteDrawingText('output/atom_map.png')

        # copy receptor pose (should not contain paul ligand!)
        new_pose = self.pyrosetta.rosetta.core.pose.Pose()
        new_pose.detached_copy(pose)

        # Add ligand residue and update information
        new_pose.append_residue_by_jump(res, 1, "", "", True)
        new_pose.pdb_info().chain(new_pose.total_residue(), 'X')
        new_pose.update_pose_chains_from_pdb_chains()

        return new_pose, new_mol, pos_to_name


    def find_closest(self, pos, pos_to_name):
        min_dist = float("inf")
        min_pos = None
        for ref_pos in pos_to_name:
            d_sq = 0.0
            for i in range(len(pos)):
                d_sq += (pos[i] - ref_pos[i]) ** 2
            if d_sq < min_dist:
                min_dist = d_sq
                min_pos = ref_pos

        print('Mapped', pos, 'to', min_pos, 'with minimal distance', min_dist, '- returning', pos_to_name[min_pos])
        return pos_to_name[min_pos]


    def add_constraints(self, pose, constraints, pos_to_name):
        for constraint in constraints:
            atom1 = self.pyrosetta.rosetta.core.id.AtomID(pose.residue(constraint[0]).type().atom_index(constraint[1]),
                                                          constraint[0])
            pos = (int(constraint[2][0] * 10), int(constraint[2][1] * 10), int(constraint[2][2] * 10))
            name = self.find_closest(pos, pos_to_name)
            atom2 = self.pyrosetta.rosetta.core.id.AtomID(
                pose.residue(pose.total_residue()).type().atom_index(name),
                pose.total_residue())
            harmonic = self.pyrosetta.rosetta.core.scoring.func.HarmonicFunc(constraint[3], constraint[4])
            cst = self.pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint(atom1, atom2, harmonic)
            pose.add_constraint(cst)


    def dock(self, job: DockingJob, dbsession, pose, runs, mover, scfx):
        work_poses = []
        for current_repeat in range(runs):

            job.progress_info = "Round " + str(current_repeat + 1) + "/" + str(runs) + "..."
            job.progress = round(10 + 80 * (current_repeat + 1) / runs)
            update(self.loop, dbsession, job)

            work_pose = self.pyrosetta.rosetta.core.pose.Pose()
            work_pose.detached_copy(pose)
            mover.apply(work_pose)

            temp_total_score = work_pose.energies().total_energy()
            temp_atom_pair_cst = work_pose.energies().total_energies()[
                self.pyrosetta.rosetta.core.scoring.ScoreType.atom_pair_constraint]

            work_pose.remove_constraints()

            interface_scores = self.pyrosetta.rosetta.protocols.ligand_docking.get_interface_deltas('X', work_pose, scfx)
            idelta_scores = {}
            score_types = ['if_X_fa_atr', 'if_X_fa_elec', 'if_X_fa_rep', 'if_X_fa_sol', 'if_X_hbond_bb_sc',
                           'if_X_hbond_sc',
                           'interface_delta_X', 'if_X_fa_pair']
            for score_type in score_types:
                idelta_scores[score_type] = interface_scores[score_type]

            result: ComplexResult = ComplexResult(
                job_id=job.job_id,
                id=current_repeat + job.runs,
                total_score=temp_total_score,
                atom_pair_cst=temp_atom_pair_cst,
                atom_attraction=idelta_scores['if_X_fa_atr'],
                electrostatic=idelta_scores['if_X_fa_elec'],
                atom_repulsion=idelta_scores['if_X_fa_rep'],
                solvation=idelta_scores['if_X_fa_sol'],
                hbond=idelta_scores['if_X_hbond_bb_sc'] + idelta_scores['if_X_hbond_sc'],
                delta_g=idelta_scores['interface_delta_X'],
                pairwise_energy=idelta_scores['if_X_fa_pair'],
            )
            work_poses.append(work_pose)

            work_poses[result.id - job.runs].dump_pdb("app/static/poses/" + job.job_id + "_" + str(result.id) + ".pdb")


            job.complexes.append(result)
            update(self.loop, dbsession, job)

        return work_poses

    def analyze_results(self, job: DockingJob, dbsession, work_poses, rdkit_mol):

        for result in filter(lambda r: r.id >= job.runs, job.complexes):
            v = []
            if result.atom_pair_cst > 15:
                v.append("ATOM_PAIR_CST")
            if result.delta_g > 0.0:
                v.append("DELTA_G")
            result.violation = ''.join(v)


        if job.runs == 0:
            any_valid = any(len(result.violation) == 0 for result in job.complexes)
            if not any_valid:
                job.complexes = []
                raise Exception("No valid docking results found. Please check your constraints and try again.")

        best_index, best_pose, best_score = 0, None, float('inf')
        for result in job.complexes:
            print("#############resultid, violation", result.id, result.violation)
            if result.delta_g < best_score and len(result.violation) == 0:
                best_score = result.delta_g
                best_pose = work_poses[result.id - job.runs] if result.id - job.runs >= 0 else None
                best_index = result.id

        # if best_pose is None:
        #     with open("app/static/poses/" + job.job_id + "#best.pkl", "rb") as f:
        #         best_pose = pickle.load(f)
        #     calc_for = filter(lambda r: r.id >= job.runs, job.complexes)
        # else:
        #     with open("app/static/poses/" + job.job_id + "#best.pkl", "wb") as f:
        #         pickle.dump(best_pose, f)
        #     calc_for = job.complexes


        all_poses = []
        if best_pose is None:
            with open("app/static/poses/" + job.job_id + "#" + str(best_index) + ".pkl", "rb") as f:
                best_pose = pickle.load(f)
            for i in range(job.runs, len(job.complexes)):
                with open("app/static/poses/" + job.job_id + "#" + str(i) + ".pkl", "wb") as f:
                    pickle.dump(work_poses[i - job.runs], f)
            calc_for_all = False
        else:
            if job.runs > 0:
                for i in range(job.runs):
                    with open("app/static/poses/" + job.job_id + "#" + str(i) + ".pkl", "rb") as f:
                        all_poses.append(pickle.load(f))
            for i in range(job.runs, len(job.complexes)):
                with open("app/static/poses/" + job.job_id + "#" + str(i) + ".pkl", "wb") as f:
                    pickle.dump(work_poses[i - job.runs], f)
            all_poses.extend(work_poses)
            calc_for_all = True


        res_selector = self.pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(
            best_pose.total_residue())
        rmsd_calc = self.pyrosetta.rosetta.core.simple_metrics.metrics.RMSDMetric()
        rmsd_calc.set_residue_selector(res_selector)
        rmsd_calc.set_comparison_pose(best_pose)

        if calc_for_all:
            for result in job.complexes:
                result.rmsd = rmsd_calc.calculate(all_poses[result.id])
        else:
            for result in filter(lambda r: r.id >= job.runs, job.complexes):
                result.rmsd = rmsd_calc.calculate(work_poses[result.id - job.runs])


        header = ['name']
        for key in job.complexes[0].model_dump().keys():
            if key not in ('pose', 'violation', 'pose_path'):
                header.append(key)
        padding = 16

        # TODO: simplify / remove debug output
        weight = Descriptors.MolWt(rdkit_mol)
        hbond_acc = Descriptors.NOCount(rdkit_mol)
        hbond_don = Descriptors.NHOHCount(rdkit_mol)
        logp = Descriptors.MolLogP(rdkit_mol)
        qed = QED.qed(rdkit_mol)

        job.weight = weight
        job.hbond_acc = hbond_acc
        job.hbond_don = hbond_don
        job.logp = logp
        job.qed = qed
        job.best_complex_nr = best_index
        update(self.loop, dbsession, job)


        def safe_format(value):
            if value is None:
                return 'N/A'.rjust(padding)
            try:
                return f'{float(value):.4f}'.rjust(padding)
            except (TypeError, ValueError):
                return str(value).rjust(padding)

        print(''.rjust(padding), ''.join(h.rjust(padding) for h in ['weight', 'hbond_acc', 'hbond_don', 'logp', 'qed']),
              sep='')
        print('Ligand Prop'.rjust(padding),
              ''.join(f'{v:.4f}'.rjust(padding) for v in [weight, hbond_acc, hbond_don, logp, qed]), sep='')
        print()

        print("".join([h.rjust(padding) for h in header]) + " violation".rjust(padding))
        print('Best'.rjust(padding),
              ''.join(safe_format(job.complexes[best_index].model_dump()[h]) for h in header[1:]),
              str(', '.join(job.complexes[best_index].violation)).rjust(padding), sep='')
        for i in range(len(job.complexes)):
            print(str(i + 1).rjust(padding),
                  ''.join(safe_format(job.complexes[i].model_dump()[h]) for h in header[1:]),
                  str(', '.join(job.complexes[i].violation)).rjust(padding), sep='')


    def generate_assets(self, job: DockingJob, work_poses, rdkit_mol):
        withoutHs = Chem.RemoveHs(rdkit_mol)
        AllChem.Compute2DCoords(withoutHs)
        drawer = rdMolDraw2D.MolDraw2DSVG(500, 500)
        drawer.DrawMolecule(rdMolDraw2D.PrepareMolForDrawing(withoutHs))
        drawer.FinishDrawing()
        with open('app/static/previews/' + job.job_id + '.svg', 'w') as f:
            f.write(drawer.GetDrawingText())


    def run_docking(self, job_id: str, runs: int):
        with Session(engine) as dbsession:
            job = dbsession.get(DockingJob, job_id)

            try:

                pybel_mol = next(pybel.readfile("sdf", 'input/ref_ligand_core.sdf'))
                pybel_sdf = pybel_mol.write('sdf')

                # RDKits sanitization step automatically detects aromaticity
                # I am assuming that it is best to provide molecules in kekulized form with all hydrogens present
                ref_mol = Chem.MolFromMolBlock(pybel_sdf)

                pose = self.pyrosetta.pose_from_pdb("input/7f83_relax.pdb")

                # DEFAULT: "COc5cc(OCc1ccncc1)cc6nc(c4cc(CCc2ccnnc2)c(OCCO)c(Cc3cnccn3)c4)[nH]c(=O)c56"
                # DEFAULT: "CN(C(=O)CN3CC2(CCN(C(=O)c1cccnc1)CC2)C3)c5ccc4COCc4c5"
                smiles = job.smiles
                lig_pose, new_mol, pos_to_name = self.prepare_pose(smiles, pose, ref_mol)

                # DEFAULT "[
                #         [92, 'HE2', 'O5', 2.5, 0.25],
                #         [45, 'HH', 'N1', 2.5, 0.25],
                #         [87, 'OH', 'O4', 2.5, 0.25],
                #     ]"
                if job.constraints:
                    self.add_constraints(lig_pose, job.constraints, pos_to_name)

                job.progress_info = "Preparing Docking Protocol..."
                job.progress = 5
                update(self.loop, dbsession, job)

                protocol_xml = self.pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_file(
                    "input/transform_std.xml")
                mover = protocol_xml.get_mover("ParsedProtocol")
                scfx = protocol_xml.get_score_function("hard_rep")

                work_poses = self.dock(job, dbsession, lig_pose, runs, mover, scfx)

                job.progress_info = "Finished Docking Process..."
                job.progress = 90
                update(self.loop, dbsession, job)

                self.analyze_results(job, dbsession, work_poses, new_mol)

                self.generate_assets(job, work_poses, new_mol)

                job.progress_info = "Finalized Calculation..."
                job.progress = 100
                job.job_status = JobStatus.COMPLETED
                job.runs = job.runs + runs
                update(self.loop, dbsession, job)

            except Exception as e:
                job.job_status = JobStatus.FAILED
                job.error = str(e)
                update(self.loop, dbsession, job)
                raise e


def update(loop: asyncio.AbstractEventLoop, dbsession: Session, job: DockingJob):
    dbsession.add(job)
    dbsession.commit()
    if job.job_id in job_update_queues:
        asyncio.run_coroutine_threadsafe(job_update_queues[job.job_id].put(True), loop)
    print("UPDATE", job.job_id, job.job_status, job.progress_info, job.progress)
