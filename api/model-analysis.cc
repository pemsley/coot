#include "molecules-container.hh"

std::vector <positioned_atom_spec_t>
molecules_container_t::get_atom_differences(int imol1, int imol2) {

    std::vector <positioned_atom_spec_t> v;
    if (is_valid_model_molecule(imol1)) {
        if (is_valid_model_molecule(imol2)) {
            mmdb::Manager *mol1 = get_mol(imol1);
            mmdb::Manager *mol2 = get_mol(imol2);
            if (mol1 && mol2) {
                int nmod1 = mol1->GetNumberOfModels();
                for (int imod = 1; imod <= nmod1; imod++) {
                    mmdb::Model *model1 = mol1->GetModel(imod);
                    mmdb::Model *model2 = mol2->GetModel(imod);
                    if (model1) {
                        if (model2) {
                            int nchains1 = model1->GetNumberOfChains();
                            for(int ichain=0; ichain<nchains1; ichain++) {
                                mmdb::Chain *chain1 = model1->GetChain(ichain);
                                mmdb::Chain *chain2 = model2->GetChain(ichain);
                                if (chain1) {
                                    if (chain2) {
                                        int nres1 = chain1->GetNumberOfResidues();
                                        for(int ires=0; ires<nres1; ires++) {
                                            mmdb::Residue *res1 = chain1->GetResidue(ires);
                                            if (res1) {
                                                int res_no = res1->GetSeqNum();
                                                const char *ins_code = res1->GetInsCode();
                                                mmdb::Residue *res2 = chain2->GetResidue(res_no, ins_code);
                                                if (res2) {
                                                    int n_atoms1 = res1->GetNumberOfAtoms();
                                                    for(int iatom=0; iatom<n_atoms1; iatom++) {
                                                        mmdb::Atom *atom1 = res1->GetAtom(iatom);
                                                        if (atom1) {
                                                            const char *atom_name = atom1->GetAtomName();
                                                            mmdb::Atom *atom2 = res2->GetAtom(atom_name);
                                                            if (atom2) {
                                                                coot::atom_spec_t as1 = coot::atom_spec_t(atom1);
                                                                coot::Cartesian pos1 = coot::Cartesian(atom1->x, atom1->y, atom1->z);
                                                                coot::Cartesian pos2 = coot::Cartesian(atom2->x, atom2->y, atom2->z);
                                                                positioned_atom_spec_t pas;
                                                                pas.atom_spec = as1;
                                                                pas.pos1 = pos1;
                                                                pas.pos2 = pos2;
                                                                v.push_back(pas);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return v;
}
