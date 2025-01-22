
#include "coot-utils/atom-tree.hh"
#include "coot-utils/coot-coord-extras.hh"
#include "coot-molecule.hh"

int
coot::molecule_t::replace_residue(const std::string &residue_cid, const std::string &new_residue_type,
                                  int imol_enc, const coot::protein_geometry &geom) {

   int status = 0;

   mmdb::Residue *residue_p = cid_to_residue(residue_cid);
   if (residue_p) {

      std::pair<bool, dictionary_residue_restraints_t> rp = geom.get_monomer_restraints(new_residue_type, imol_enc);
      if (rp.first) {
      const auto &restraints_new_type = rp.second;

         std::string current_residue_type = residue_p->GetResName();
         std::pair<bool, dictionary_residue_restraints_t> rp_current = geom.get_monomer_restraints(current_residue_type, imol_enc);
         if (rp_current.first) {
            const auto &restraints_current_type = rp_current.second;
            status = util::mutate_by_overlap(residue_p, atom_sel.mol, restraints_current_type, restraints_new_type);
         }
      }
   }
   return status;
}


int
coot::molecule_t::rotate_around_bond(const std::string &residue_cid,
                                     const std::string &alt_conf,
                                     coot::atom_name_quad quad,
                                     double torsion_angle,
                                     coot::protein_geometry &geom) {

   int status = 0;
   double r = -999.9;

   mmdb::Residue *residue_p = cid_to_residue(residue_cid);
   if (residue_p) {
      std::string res_name(residue_p->GetResName());
      std::pair<bool, coot::dictionary_residue_restraints_t> restraints_info =
         geom.get_monomer_restraints(res_name, imol_no);
      if (! restraints_info.first) {
         std::cout << "WARNING:: set_torsion: No restraints for " << res_name << std::endl;
      } else {
         coot::atom_tree_t tree(restraints_info.second, residue_p, alt_conf);

         try {
            r = tree.set_dihedral(quad.atom_name(0),
                                  quad.atom_name(1),
                                  quad.atom_name(2),
                                  quad.atom_name(3),
                                  torsion_angle);
            atom_sel.mol->FinishStructEdit();
         }
         catch(const std::runtime_error &rte) {
            std::cout << "in set_torsion:: set_dihedral() error: " << rte.what() << std::endl;
         }

      }
   } else {
      std::cout << "failed to find residue " << residue_cid << std::endl;
   }
   return status;
}


//! copy chain using NCS matrix
bool
coot::molecule_t::copy_ncs_chain(const std::string &from_chain_id, const std::string &to_chain_id) {

   bool status = false;

   return status;

}
