
#include "rigid-body-fit.hh"
#include "ligand/ligand.hh"
#include "coot-utils/atom-selection-container.hh"
#include "coot-utils/coot-coord-utils.hh"

coot::minimol::molecule
coot::api::rigid_body_fit_inner(const coot::minimol::molecule &mol_without_moving_atoms,
                                const coot::minimol::molecule &mol_for_moving_atoms,
                                const clipper::Xmap<float> &xmap) {

   bool debug = false;
   // debugging
   if (debug) {
      mol_without_moving_atoms.write_file("rigid-body-moving-atoms-mol.pdb", 44);
      mol_without_moving_atoms.write_file("rigid-body-without-moving-atoms.pdb", 44);
   }

   coot::ligand lig;
   bool mask_water_flag = false;
   float sigma = 0.4; // FIXME
   lig.import_map_from(xmap, sigma);

   lig.install_ligand(mol_for_moving_atoms);
   lig.find_centre_by_ligand(0); // don't test ligand size
   // lig.set_map_atom_mask_radius(0.5);  // why did I have this?
   lig.mask_map(mol_without_moving_atoms, mask_water_flag);
   if (debug)
      lig.output_map("rigid-body.map");
   lig.set_dont_write_solutions();
   lig.set_dont_test_rotations();
   lig.set_acceptable_fit_fraction(0.5);
   lig.fit_ligands_to_clusters(1);
   unsigned int iclust = 0;
   unsigned int isol   = 0;
   coot::minimol::molecule moved_mol = lig.get_solution(iclust, isol);
   return moved_mol;
}

#include "geometry/residue-and-atom-specs.hh"

void
coot::api::rigid_body_fit(mmdb::Manager *mol, int udd_atom_selection_fitting_atoms,
                          const clipper::Xmap<float> &xmap) {

   // make a minimol that is the atoms of the atom selection
   // and
   // make a minimol that is the atoms that are not in the atom selection

   if (true) { // debugging
      mmdb::PAtom *atoms = NULL;
      int n_atoms;
      mol->GetSelIndex(udd_atom_selection_fitting_atoms, atoms, n_atoms);
      std::cout << "----------- debug:: in rigid_body_fit() we selected " << n_atoms << " atoms " << std::endl;
   }

   bool fill_masking_molecule_flag = true;
   coot::ligand lig;
   std::pair<coot::minimol::molecule, coot::minimol::molecule> p = coot::make_mols_from_atom_selection(mol,
                                                                                                       udd_atom_selection_fitting_atoms,
                                                                                                       fill_masking_molecule_flag);
   const minimol::molecule &mol_without_moving_atoms = p.first;
   const minimol::molecule &mol_for_moving_atoms     = p.second;

   minimol::molecule moved_atoms_mol = rigid_body_fit_inner(mol_without_moving_atoms, mol_for_moving_atoms, xmap);

   mol_without_moving_atoms.write_file("mol_without_moving_atoms.pdb", true);
   mol_for_moving_atoms.write_file("mol_for_moving_atoms.pdb", true);

   // now transfer the new positions of the atoms in moved_atoms_mol into mol

   unsigned int n_moved = 0;
   for(unsigned int ifrag=0; ifrag<moved_atoms_mol.fragments.size(); ifrag++) {
      for(int ires=moved_atoms_mol[ifrag].min_res_no(); ires<=moved_atoms_mol[ifrag].max_residue_number(); ires++) {
         std::string ins_code = "";
         residue_spec_t this_residue_spec(moved_atoms_mol[ifrag].fragment_id, ires, ins_code);
         mmdb::Residue *this_residue = util::get_residue(this_residue_spec, mol);
         if (this_residue) {
            for (unsigned int iat=0; iat<moved_atoms_mol[ifrag][ires].atoms.size(); iat++) {

               if (false)
                  std::cout << "replacing " << moved_atoms_mol[ifrag].fragment_id << " " << moved_atoms_mol[ifrag][ires]
                            << " " << moved_atoms_mol[ifrag][ires][iat].name
                            << " " << moved_atoms_mol[ifrag][ires][iat].pos.format() << std::endl;

               mmdb::Atom **residue_atoms = 0;
               int n_residue_atoms = 0;
               this_residue->GetAtomTable(residue_atoms, n_residue_atoms);
               for (int jat=0; jat<n_residue_atoms; jat++) {
                  mmdb::Atom *at = residue_atoms[jat];
                  if (! at->isTer()) {
                     std::string this_atom_name(at->GetAtomName());
                     // std::cout << "   comparing names " << moved_atoms_mol[ifrag][ires][iat].name << " "
                     // << this_atom_name <<std::endl;
                     if (moved_atoms_mol[ifrag][ires][iat].name == this_atom_name) {
                        std::string this_atom_alt_conf(at->altLoc);
                        if (this_atom_alt_conf == moved_atoms_mol[ifrag][ires][iat].altLoc) {
                           at->x = moved_atoms_mol[ifrag][ires][iat].pos.x();
                           at->y = moved_atoms_mol[ifrag][ires][iat].pos.y();
                           at->z = moved_atoms_mol[ifrag][ires][iat].pos.z();
                           n_moved++;
                           // std::cout << "      moving " << coot::atom_spec_t(at) << std::endl;
                        }
                     }
                  }
               }
            }
	 }
      }
   }

   if (true)
      std::cout << "DEBUG:: in rigid_body_fit() moved " << n_moved << " atoms " << std::endl;
}
