import traceback

from pyrosetta.rosetta.cppdb import session
from rdkit import Chem
from rdkit.Chem import AllChem, Lipinski, Descriptors, QED
from openbabel import pybel
import matplotlib
from rdkit.Chem.Draw import rdMolDraw2D
from sqlmodel import Session

from .config import settings
from .db.db import engine
from .models import DockingJob, ComplexResult

matplotlib.use("agg")


class DockingWrapper:

    def __init__(self):
        import pyrosetta
        self.pyrosetta = pyrosetta
        self.pyrosetta.init()

        self.current_job: DockingJob | None = None
        self.current_results: list[ComplexResult] = []

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
            # Receptor atom selected via Residue index and atom name
            atom1 = self.pyrosetta.rosetta.core.id.AtomID(pose.residue(constraint[0]).type().atom_index(constraint[1]),
                                                          constraint[0])
            # ligand atom selected from last residue (aka ligand residue) and atom name
            pos = (int(constraint[2][0] * 10), int(constraint[2][1] * 10), int(constraint[2][2] * 10))
            name = self.find_closest(pos, pos_to_name)
            atom2 = self.pyrosetta.rosetta.core.id.AtomID(
                pose.residue(pose.total_residue()).type().atom_index(name),
                pose.total_residue())
            harmonic = self.pyrosetta.rosetta.core.scoring.func.HarmonicFunc(constraint[3], constraint[4])
            cst = self.pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint(atom1, atom2, harmonic)
            pose.add_constraint(cst)

    def dock(self, pose, runs, all_results: list[ComplexResult], mover, scfx, repeats):
        for current_repeat in range(repeats):

            self.current_job.update(
                progress_info="Round " + str(current_repeat + 1) + "/" + str(runs) + "...",
                progress=round(10 + 80 * (current_repeat + 1) / runs))

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

            results: ComplexResult = ComplexResult(
                job_id=self.current_job.job_id,
                id=current_repeat + self.current_job.runs,
                pose=work_pose,
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

            all_results.append(results)
            self.current_results.append(results)

            with Session(engine) as session:
                session.add(results)
                session.commit()

    def analyze_results(self, results: list[ComplexResult], rdkit_mol):

        # TODO: HANDLE FALSE RESULTS FRONTEND
        # filtered_results = []
        # cst_penalty_violations = 0
        # delta_g_violations = 0
        # for result in results:
        #     keep = True
        #     if result['atom_pair_cst'] > max_cst_penalty:
        #         keep = False
        #         cst_penalty_violations += 1
        #     if result['delta_g'] > 0.0:
        #         keep = False
        #         delta_g_violations += 1
        #     if keep:
        #         filtered_results.append(result)
        # print(len(filtered_results), 'of', len(results), 'results were successful.')
        # print(cst_penalty_violations, 'results were discarde due to constraint penalties above', max_cst_penalty)
        # print(delta_g_violations, 'results were discarde due to positive delta gs')
        # if len(filtered_results) == 0:
        #     print('ERROR: No successful results. Skipping analysis.')
        #     return

        for index, result in enumerate(results):
            if result.atom_pair_cst > 15:
                result.violation.append("ATOM_PAIR_CST")
            if result.delta_g > 0.0:
                result.violation.append("DELTA_G")

        any_valid = any(len(result.violation) == 0 for result in results)
        if not any_valid:
            self.current_job.update(
                error="No valid docking results found. Please check your constraints and try again.")
            print("No valid docking results found. Please check your constraints and try again.")
            return

        best_index, best_pose, best_score = 0, None, float('inf')
        for idx, result in enumerate(results):
            if result.delta_g < best_score and len(result.violation) == 0:
                best_score = result.delta_g
                best_pose = result.pose
                best_index = idx

        res_selector = self.pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(
            best_pose.total_residue())
        rmsd_calc = self.pyrosetta.rosetta.core.simple_metrics.metrics.RMSDMetric()
        rmsd_calc.set_residue_selector(res_selector)
        rmsd_calc.set_comparison_pose(best_pose)

        for result in results:
            result.rmsd = rmsd_calc.calculate(result.pose)

        # TODO: HANDLE THIS IN FRONTEND
        # x = []
        # y = []
        # for result in results:
        #     rmsd = result['rmsd']
        #     score = result['delta_g']
        #     x.append(rmsd)
        #     y.append(score)
        #
        # plt.scatter(x, y)
        # plt.ylabel('Delta G (REU)')
        # plt.xlabel('RMSD')
        # plt.savefig('output/funnel.png')
        #
        # best_pose.dump_pdb("output/best_pose.pdb")

        header = ['name']
        for key in results[0].model_dump().keys():
            if key != 'pose' and key != 'violation' and key != 'pose_path':
                header.append(key)
        padding = 16

        # TODO: simplify / remove debug output
        weight = Descriptors.MolWt(rdkit_mol)
        hbond_acc = Descriptors.NOCount(rdkit_mol)
        hbond_don = Descriptors.NHOHCount(rdkit_mol)
        logp = Descriptors.MolLogP(rdkit_mol)
        qed = QED.qed(rdkit_mol)

        with Session(engine) as session:
            job = session.get(DockingJob, self.current_job.job_id)

            for i, db_complex in enumerate(job.complexes):
                if i < len(results):
                    db_complex.rmsd = results[i].rmsd
            session.commit()

        self.current_job.update(
            weight=weight,
            hbond_acc=hbond_acc,
            hbond_don=hbond_don,
            logp=logp,
            qed=qed,
            best_complex_id=best_index
        )

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
              ''.join(safe_format(results[best_index].model_dump()[h]) for h in header[1:]),
              str(', '.join(results[best_index].violation)).rjust(padding), sep='')
        for i in range(len(results)):
            print(str(i + 1).rjust(padding),
                  ''.join(safe_format(results[i].model_dump()[h]) for h in header[1:]),
                  str(', '.join(results[i].violation)).rjust(padding), sep='')

    def generate_assets(self, job_id, results, rdkit_mol):
        withoutHs = Chem.RemoveHs(rdkit_mol)
        AllChem.Compute2DCoords(withoutHs)
        drawer = rdMolDraw2D.MolDraw2DSVG(500, 500)
        drawer.DrawMolecule(rdMolDraw2D.PrepareMolForDrawing(withoutHs))
        drawer.FinishDrawing()
        with open('app/static/previews/' + job_id + '.svg', 'w') as f:
            f.write(drawer.GetDrawingText())

        for result in results:
            result.pose.dump_pdb("app/static/poses/" + result.job_id + "#" + result.id + ".pdb")


    def run_docking(self, job: DockingJob, runs: int):

        try:
            self.current_job = job

            self.current_job.update(progress_info="Preparing Ligand...", progress=0)

            pybel_mol = next(pybel.readfile("sdf", 'input/ref_ligand_core.sdf'))
            pybel_sdf = pybel_mol.write('sdf')

            # RDKits sanitization step automatically detects aromaticity
            # I am assuming that it is best to provide molecules in kekulized form with all hydrogens present
            ref_mol = Chem.MolFromMolBlock(pybel_sdf)

            pose = self.pyrosetta.pose_from_pdb("input/7f83_relax.pdb")

            # DEFAULT: "COc5cc(OCc1ccncc1)cc6nc(c4cc(CCc2ccnnc2)c(OCCO)c(Cc3cnccn3)c4)[nH]c(=O)c56"
            # DEFAULT: "CN(C(=O)CN3CC2(CCN(C(=O)c1cccnc1)CC2)C3)c5ccc4COCc4c5"
            smiles = self.current_job.smiles
            lig_pose, new_mol, pos_to_name = self.prepare_pose(smiles, pose, ref_mol)

            # DEFAULT "[
            #         [92, 'HE2', 'O5', 2.5, 0.25],
            #         [45, 'HH', 'N1', 2.5, 0.25],
            #         [87, 'OH', 'O4', 2.5, 0.25],
            #     ]"
            if self.current_job.constraints:
                self.add_constraints(lig_pose, self.current_job.constraints, pos_to_name)

            self.current_job.update(progress_info="Preparing Docking Protocol...", progress=5)

            protocol_xml = self.pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_file(
                "input/transform_std.xml")
            mover = protocol_xml.get_mover("ParsedProtocol")
            scfx = protocol_xml.get_score_function("hard_rep")

            results = []
            self.dock(lig_pose, runs, results, mover, scfx, repeats=runs)

            self.current_job.update(progress_info="Finished Docking Process...", progress=90)

            self.analyze_results(results, new_mol)

            with Session(engine) as dbsession:
                job = dbsession.get(DockingJob, job.job_id)
                self.current_job.update(
                    job_status="COMPLETED",
                    runs=self.current_job.runs + runs,
                )
                job.job_status = self.current_job.job_status
                job.runs = self.current_job.runs

                job.best_complex_nr = self.current_job.best_complex_nr
                job.weight = self.current_job.weight
                job.hbond_acc = self.current_job.hbond_acc
                job.hbond_don = self.current_job.hbond_don
                job.logp = self.current_job.logp
                job.qed = self.current_job.qed
                job.error = self.current_job.error

                dbsession.add(job)
                dbsession.commit()

            self.generate_assets(job.job_id, results, new_mol)

            self.current_job.update(status="Finalized Calculation...", progress=100)

        except Exception as e:
            print("Docking failed with error:", e)
            traceback.print_exc()
            self.current_job.update(error=str(e))

