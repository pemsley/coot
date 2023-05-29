
#include "protein-geometry.hh"

void
coot::dictionary_residue_restraints_t::move_3GP_atoms() {

   std::map<std::string, clipper::Coord_orth> position_map;
   auto add_ref_pos = [&position_map] (const std::string &atom_name, const clipper::Coord_orth &co) {
      position_map[atom_name] = co;
   };

   auto position = [] (double x, double y, double z) {
     return clipper::Coord_orth(x,y,z);
   };

   add_ref_pos(" O6 ", position(-0.191 ,  0.192 ,  0.457 ));
   add_ref_pos(" C6 ", position( 0.721 ,  1.013 ,  0.364 ));
   add_ref_pos(" C5 ", position( 2.046 ,  0.587 ,  0.335 ));
   add_ref_pos(" C4 ", position( 3.070 ,  1.523 ,  0.230 ));
   add_ref_pos(" N3 ", position( 2.781 ,  2.849 ,  0.156 ));
   add_ref_pos(" C2 ", position( 1.465 ,  3.282 ,  0.185 ));
   add_ref_pos(" N2 ", position( 1.189 ,  4.581 ,  0.113 ));
   add_ref_pos("HN21", position( 0.295 ,  4.888 , -0.258 ));
   add_ref_pos("HN22", position( 1.832 ,  5.258 ,  0.510 ));
   add_ref_pos(" N1 ", position( 0.436 ,  2.364 ,  0.289 ));
   add_ref_pos(" HN1", position(-0.552 ,  2.689 ,  0.311 ));
   add_ref_pos(" N7 ", position( 2.607 , -0.640 ,  0.390 ));
   add_ref_pos(" C8 ", position( 3.951 , -0.507 ,  0.323 ));
   add_ref_pos(" H8 ", position( 4.670 , -1.325 ,  0.345 ));
   add_ref_pos(" N9 ", position( 4.255 ,  0.827 ,  0.224 ));
   add_ref_pos(" C1'", position( 5.607 ,  1.416 ,  0.126 ));
   add_ref_pos(" H1'", position( 5.574 ,  2.306 , -0.502 ));
   add_ref_pos(" O4'", position( 6.528 ,  0.487 , -0.417 ));
   add_ref_pos(" C4'", position( 7.768 ,  0.519 ,  0.268 ));
   add_ref_pos(" H4'", position( 8.520 ,  0.983 , -0.370 ));
   add_ref_pos(" C5'", position( 8.234 , -0.897 ,  0.604 ));
   add_ref_pos("H5'1", position( 9.028 , -0.857 ,  1.349 ));
   add_ref_pos("H5'2", position( 8.639 , -1.368 , -0.291 ));
   add_ref_pos(" O5'", position( 7.162 , -1.670 ,  1.097 ));
   add_ref_pos("HO5'", position( 6.941 , -1.384 ,  2.007 ));
   add_ref_pos(" C2'", position( 6.131 ,  1.812 ,  1.500 ));
   add_ref_pos(" H2'", position( 5.590 ,  1.255 ,  2.265 ));
   add_ref_pos(" O2'", position( 5.989 ,  3.198 ,  1.714 ));
   add_ref_pos("HO2'", position( 5.441 ,  3.351 ,  2.511 ));
   add_ref_pos(" C3'", position( 7.587 ,  1.380 ,  1.513 ));
   add_ref_pos(" H3'", position( 7.771 ,  0.783 ,  2.406 ));
   add_ref_pos(" O3'", position( 8.455 ,  2.493 ,  1.479 ));
   add_ref_pos(" P  ", position( 9.987 ,  2.376 ,  1.962 ));
   add_ref_pos(" O1P", position(10.003 ,  1.895 ,  3.393 ));
   add_ref_pos(" O3P", position(10.730 ,  1.394 ,  1.088 ));
   add_ref_pos(" O2P", position(10.649 ,  3.730 ,  1.876 ));

   for (unsigned int i=0; i<atom_info.size(); i++) {
      auto &atom = atom_info[i];
      const std::string &atom_name = atom.atom_id_4c;
      std::map<std::string, clipper::Coord_orth>::const_iterator it;
      it = position_map.find(atom_name);
      if (it != position_map.end()) {
         auto prev_1 = atom.pdbx_model_Cartn_ideal.second;
         auto prev_2 = atom.model_Cartn.second;
         atom.pdbx_model_Cartn_ideal.second = it->second;
         atom.model_Cartn.second = it->second;
         std::cout << "replacing " << atom_name << " " << atom.pdbx_model_Cartn_ideal.second.format()
                   << " was " << prev_1.format() << " " << prev_2.format() << std::endl;
      } else {
         std::cout << "3GP move fail! " << atom_name << std::endl;
      }
   }

   for (unsigned int i=0; i<atom_info.size(); i++) {
      auto &atom = atom_info[i];
      const std::string &atom_name = atom.atom_id_4c;
      std::cout << "After move: " << atom_name << " " << atom.model_Cartn.first << " " << atom.model_Cartn.second.format()
                << std::endl;
   }

}
