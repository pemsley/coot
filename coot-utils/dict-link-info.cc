
#include "dict-link-info.hh"
#include "coot-coord-utils.hh"

#ifdef USE_BACKWARD
#include <utils/backward.hpp>
#endif

// this can throw a std::runtime_error
coot::dict_link_info_t::dict_link_info_t(mmdb::Residue *residue_ref,
                                         mmdb::Residue *residue_new,
                                         const std::string &link_type,
                                         const coot::protein_geometry &geom) {

   if (false) {
      std::cout << "dict_link_info_t() constructor start with link_type: "
                << link_type << std::endl;
      std::cout << "dict_link_info_t() constructor start with residue_ref: "
                << coot::residue_spec_t(residue_ref) << std::endl;
      std::cout << "dict_link_info_t() constructor start with residue_new: "
                << coot::residue_spec_t(residue_new) << std::endl;
   }

   bool ifound = false;
   if (! residue_ref) {
      throw (std::runtime_error("Null residue_ref"));
   } else {
      if (! residue_ref) {
         throw (std::runtime_error("Null residue_new"));
      } else {
         coot::dictionary_residue_link_restraints_t rr = geom.link(link_type);
         if (rr.link_id == "") {
            std::string rn1 = residue_ref->GetResName();
            std::string rn2 = residue_new->GetResName();
            std::string mess = "Link not found in dictionary " + link_type + " between " + rn1 + " " + rn2;

            throw (std::runtime_error(mess));
         } else {

            bool order_switch_flag = check_for_order_switch(residue_ref,
                                                            residue_new,
                                                            link_type, geom);

            mmdb::Residue *res_1 = residue_ref;
            mmdb::Residue *res_2 = residue_new;

            if (order_switch_flag) {
               std::swap(res_1, res_2);
            }

            // we found it (i.e. not null object)
            coot::residue_spec_t res_spec_ref(res_1);
            coot::residue_spec_t res_spec_new(res_2);
            for (unsigned int ibond=0; ibond<rr.link_bond_restraint.size(); ibond++) {
               mmdb::PPAtom residue_atoms_1 = 0;
               int n_residue_atoms_1;
               res_1->GetAtomTable(residue_atoms_1, n_residue_atoms_1);
               for (int iat1=0; iat1<n_residue_atoms_1; iat1++) {
                  std::string atom_name_1(residue_atoms_1[iat1]->name);
                  if (atom_name_1 == rr.link_bond_restraint[ibond].atom_id_1_4c()) {
                     // OK so the first atom matched
                     mmdb::PPAtom residue_atoms_2 = 0;
                     int n_residue_atoms_2;
                     res_2->GetAtomTable(residue_atoms_2, n_residue_atoms_2);
                     for (int iat2=0; iat2<n_residue_atoms_2; iat2++) {
                        std::string atom_name_2(residue_atoms_2[iat2]->name);
                        if (atom_name_2 == rr.link_bond_restraint[ibond].atom_id_2_4c()) {
                           ifound = 1;
                           spec_ref = coot::atom_spec_t(res_spec_ref.chain_id,
                                                        res_spec_ref.res_no,
                                                        res_spec_ref.ins_code,
                                                        atom_name_1, "");
                           spec_new = coot::atom_spec_t(res_spec_new.chain_id,
                                                        res_spec_new.res_no,
                                                        res_spec_new.ins_code,
                                                        atom_name_2, "");

                           dist = coot::distance(residue_atoms_1[iat1],
                                                 residue_atoms_2[iat2]);
                           break;
                        }
                        if (ifound)
                           break;
                     }
                  }
                  if (ifound)
                     break;
               }
               if (ifound)
                  break;
            }

            if (! ifound)
               throw std::runtime_error("Dictionary links atom not found in link residues");
         }
      }
   }
}

bool
coot::dict_link_info_t::check_for_order_switch(mmdb::Residue *residue_ref,
                                               mmdb::Residue *residue_new,
                                               const std::string &link_type,
                                               const coot::protein_geometry &geom) const {

   bool order_switch_flag = false;
   std::string comp_id_ref = residue_ref->GetResName();
   std::string comp_id_new = residue_new->GetResName();

   try {
      std::string group_ref = geom.get_group(residue_ref);
      // std::cout << "got group_ref: " << group_ref << std::endl;
      std::string group_new = geom.get_group(residue_new);
      // std::cout << "got group_new:o " << group_new << std::endl;
      std::vector<std::pair<coot::chem_link, bool> > link_infos;
      std::vector<coot::chem_link> link_infos_f = geom.matching_chem_links(comp_id_ref, group_ref, comp_id_new, group_new);
      std::vector<coot::chem_link> link_infos_b = geom.matching_chem_links(comp_id_new, group_new, comp_id_ref, group_ref);
      for (const auto &link : link_infos_f) link_infos.push_back(std::make_pair(link, false));
      for (const auto &link : link_infos_b) link_infos.push_back(std::make_pair(link,  true));

      if (false)
         std::cout << "DEBUG:: in check_for_order_switch() found " << link_infos.size()
                   << " link infos " << std::endl;

      for (unsigned int ilink=0; ilink<link_infos.size(); ilink++) {
         // std::cout << "   chem_link: " << ilink << " " << link_infos[ilink].first
         // << " " << link_infos[ilink].second << std::endl;
         if (link_infos[ilink].first.Id() == link_type) {
            order_switch_flag = link_infos[ilink].second;
            break;
         }
      }
   }
   catch (const std::runtime_error &rte) {
      std::cout << "WARNING:: check_for_order_switch() exception: " << rte.what() << std::endl;
   }
   return order_switch_flag;
}
