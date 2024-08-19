
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
