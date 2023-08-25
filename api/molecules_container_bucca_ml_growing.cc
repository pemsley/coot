
#include "molecules_container.hh"
#include "buccaneer_ml_growing/ml-grow.h"
#include "mini-mol/mini-mol-utils.hh"

//! buccaneer building, called by the above
int
molecules_container_t::add_terminal_residue_directly_using_bucca_ml_growing(int imol, const coot::residue_spec_t &spec) {

/*
   auto delete_ter = [] (mmdb::Residue *residue_p) {
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (at->isTer()) {
            delete at;
            break;
         }
      }
   };

   auto add_or_insert_residue = [delete_ter](mmdb::Manager *mol, mmdb::Chain *chain_p, mmdb::Residue *residue_p,
                                   int serial_number) {

      if (serial_number == -1) {
         mmdb::Residue *last_residue = coot::util::get_last_residue_in_chain(chain_p);
         if (last_residue) {
            chain_p->AddResidue(residue_p);
            delete_ter(last_residue);
         }
      } else {
         chain_p->InsResidue(residue_p, serial_number);
      }
   };

   int status = 0;
   int  imol_map = imol_refinement_map;
   if (is_valid_map_molecule(imol_map)) {
      if (is_valid_model_molecule(imol)) {
         mmdb::Manager *mol = get_mol(imol);
         mmdb::Residue *residue_p = molecules[imol].get_residue(spec);
         if (residue_p) {
            mmdb::Chain *chain_p = residue_p->GetChain();
            int res_no_this = spec.res_no;
            int res_no_next = spec.res_no + 1;
            int res_no_prev = spec.res_no - 1;
            const std::string &chain_id = spec.chain_id;
            // Ca_group next_ca_group(const Ca_group &ca_group, const clipper::Xmap<float> &xmap);
            // Ca_group prev_ca_group(const Ca_group &ca_group, const clipper::Xmap<float> &xmap);

            coot::residue_spec_t spec_next = spec.next();
            coot::residue_spec_t spec_prev = spec.previous();

            mmdb::Residue *residue_next_p = molecules[imol].get_residue(spec_next);
            mmdb::Residue *residue_prev_p = molecules[imol].get_residue(spec_prev);

            clipper::Xmap<float> &xmap = molecules[imol_map].xmap;

            // bad cell? test here
            // ::api::cell_t c = get_cell(imol_map);
            // std::cout << "debug in add_terminal_residue_directly_using_bucca_ml_growing(); cell: "
            // << c.a << " " << c.b << " " << c.c << " " << c.alpha << " " << c.beta << " " << c.gamma << std::endl;

            bool found_n_pos  = false;
            bool found_ca_pos = false;
            bool found_c_pos  = false;
            clipper::Coord_orth c_pos;
            clipper::Coord_orth n_pos;
            clipper::Coord_orth ca_pos;
            mmdb::Atom **residue_atoms = 0;
            int n_residue_atoms = 0;
            mmdb::Atom *o_at_this_residue = nullptr; // O atom need to be added after addition of next residue
            residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
            for (int iat=0; iat<n_residue_atoms; iat++) {
               mmdb::Atom *at = residue_atoms[iat];
               std::string name(at->name);
               if (! at->isTer()) {
                  if (name == " CA ") { found_ca_pos = true;  ca_pos = coot::co(at); }
                  if (name == " C  ") { found_c_pos  = true;  c_pos  = coot::co(at); }
                  if (name == " N  ") { found_n_pos  = true;  n_pos  = coot::co(at); }
                  if (name == " O  ") { o_at_this_residue = at; }
               }
            }
            if (found_n_pos && found_ca_pos && found_c_pos) {

               if (false) {
                  std::cout << "r_this_group:  n_pos " <<  n_pos.format() << std::endl;
                  std::cout << "r_this_group: ca_pos " << ca_pos.format() << std::endl;
                  std::cout << "r_this_group:  c_pos " <<  c_pos.format() << std::endl;
               }
               Ca_group r_this_group(n_pos, ca_pos, c_pos);
               if (residue_prev_p == nullptr) {
                  Ca_group r_next = ml_grow::prev_ca_group(r_this_group, xmap);
               }
               if (residue_next_p == nullptr) {
                  Ca_group r_next = ml_grow::next_ca_group(r_this_group, xmap);
                  clipper::Coord_orth ca_next_pos = r_next.coord_ca();
                  clipper::Coord_orth  c_next_pos = r_next.coord_c();
                  clipper::Coord_orth  n_next_pos = r_next.coord_n();
                  clipper::Coord_orth cb_next_pos = r_next.coord_cb();

                  coot::minimol::residue r_this_mm(residue_p);
                  coot::minimol::residue r_next_mm(res_no_next);

                  // the CB may not be there - so don't try to fit it to the density as much as we do
                  // for the other atoms.
                  r_next_mm.addatom(" CA ", " C", ca_next_pos.x(), ca_next_pos.y(), ca_next_pos.z(), "", 1.0, 30.0);
                  r_next_mm.addatom(" C  ", " C",  c_next_pos.x(),  c_next_pos.y(),  c_next_pos.z(), "", 1.0, 30.0);
                  r_next_mm.addatom(" N  ", " N",  n_next_pos.x(),  n_next_pos.y(),  n_next_pos.z(), "", 1.0, 30.0);
                  r_next_mm.addatom(" CB ", " C", cb_next_pos.x(), cb_next_pos.y(), cb_next_pos.z(), "", 0.13, 30.0);

                  std::pair<bool, clipper::Coord_orth> o_position_pair = coot::o_position(r_this_mm, r_next_mm);

                  if (o_position_pair.first) {

                     const auto &o_pos = o_position_pair.second;
                     if (o_at_this_residue == nullptr) {
                        // add the O to this residue
                        mmdb::Atom *new_atom = new mmdb::Atom;
                        new_atom->SetAtomName(" O  ");
                        new_atom->SetElementName(" O");
                        new_atom->SetCoordinates(o_pos.x(), o_pos.y(), o_pos.z(), 1.0, 30.0);
                        residue_p->AddAtom(new_atom);
                     } else {
                        // or move just move it, if it already exists
                        o_at_this_residue->x = o_pos.x();
                        o_at_this_residue->y = o_pos.y();
                        o_at_this_residue->z = o_pos.z();
                     }
                     r_next_mm.name = "ALA";
                     mmdb::Residue *r_next_p = r_next_mm.make_residue();

                     std::string ins_code;
                     int serial_number = find_serial_number_for_insert(res_no_next, ins_code, chain_p);
                     add_or_insert_residue(mol, chain_p, r_next_p, serial_number);
                     coot::copy_segid(residue_p, r_next_p);
                     mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
                     mol->FinishStructEdit();
                     coot::util::pdbcleanup_serial_residue_numbers(mol);
                     status = 1;
                  }
               }
            }
         }
      }
   }
   return status;
   */
   return 0;
}
