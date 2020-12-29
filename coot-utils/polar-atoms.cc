
#include <iostream>
#include <map>

#include "geometry/protein-donor-acceptors.hh"
#include "polar-atoms.hh"

#include "contacts-by-bricks.hh"
#include "geometry/residue-and-atom-specs.hh"

#include "analysis/daca.hh"

// This molecule should have hydrogen atoms.
void
coot::buried_unsatisfied_polar_atoms(mmdb::Manager *mol) {

   std::set<unsigned int> fixed_atom_indices; // empty

   int SelectionHandle = mol->NewSelection(); // d
   mol->SelectAtoms (SelectionHandle, 0, "*",
                     mmdb::ANY_RES, // starting resno, an int
                     "*", // any insertion code
                     mmdb::ANY_RES, // ending resno
                     "*", // ending insertion code
                     "*", // any residue name
                     "*", // atom name
                     "*", // elements
                     "*"  // alt loc.
                     );

   // -: ASP GLU
   // +: LYS ARG HIS
   // H: ASN GLN SER THR TYR

   mmdb::Atom **atom_selection = 0;
   int n_selected_atoms = 0;
   mol->GetSelIndex(SelectionHandle, atom_selection, n_selected_atoms);
   std::map<std::string, bool> is_polar_map;

   coot::quick_protein_donor_acceptors pda;

   std::cout << "selected " << n_selected_atoms << " atoms " << std::endl;
   if (n_selected_atoms > 0) {
      contacts_by_bricks contacts(atom_selection, n_selected_atoms, fixed_atom_indices);
      // max distance between hydrogen atom and a hydrogen bond acceptor
      contacts.set_dist_max(3.0); // Hmmm
      std::vector<std::set<unsigned int> > vec;
      contacts.find_the_contacts(&vec);

#if 0
      daca d;
      std::vector<std::pair<mmdb::Residue *, float> > se = d.solvent_exposure(mol);
      std::map<mmdb::Residue *, float> se_map;
      for (auto p : se)
         se_map[p.first] = p.second;
#endif

      std::vector<bool> is_polar(n_selected_atoms, false);
      for (int i=0; i<n_selected_atoms; i++) {
         mmdb::Atom *at = atom_selection[i];
         std::string res_name(at->residue->GetResName());
         std::string at_name(at->GetAtomName());
         quick_protein_donor_acceptors::key k(res_name, at_name);
         hb_t hb_type = pda.get_type(k);
         if (hb_type == HB_BOTH || hb_type == HB_DONOR || hb_type == HB_ACCEPTOR || hb_type == HB_HYDROGEN)
            is_polar[i] = true;
      }

      if (!vec.empty()) {
         int vs = vec.size();
         if (vs != n_selected_atoms) {
            std::cout << "size problem " << vs << " " << n_selected_atoms << std::endl;
         } else {
            for (int i=0; i<n_selected_atoms; i++) {
               if (is_polar[i]) {
                  if (! vec[i].empty()) {
                     mmdb::Atom *at = atom_selection[i];
                     if (at) {
                        if (! at->isTer()) {
                           bool found_something = false;
                           std::string res_name(at->residue->GetResName());
                           std::string at_name(at->GetAtomName());
                           quick_protein_donor_acceptors::key key_1(res_name, at_name);
                           std::set<unsigned int>::const_iterator it;
                           for (it=vec[i].begin(); it!=vec[i].end(); ++it) {
                              mmdb::Atom *at_neighb = atom_selection[*it];
                              if (at_neighb) {
                                 if (! at_neighb->isTer()) {
                                    std::string res_name_n(at->residue->GetResName());
                                    std::string at_name_n(at->GetAtomName());
                                    quick_protein_donor_acceptors::key key_2(res_name_n, at_name_n);
                                    std::pair<bool, bool> is_valid = pda.is_hydrogen_bond_by_types(key_1, key_2);
                                    if (is_valid.first)
                                       if (is_valid.second)
                                          found_something = true;
                                 }
                              }
                           }
                           if (! found_something) {
                              // std::cout << "Nothing was found for atom " << atom_spec_t(at) << std::endl;
                           }

                           // mmdb::Residue *residue_p = at->residue;
                           // std::cout << atom_spec_t(at) << " " << se_map[residue_p] << " " << found_something
                           // << std::endl;
                        }
                     }
                  }
               }
            }
         }
      } else {
         std::cout << "empty vec - sad face" << std::endl;
      }
   }
   
}
