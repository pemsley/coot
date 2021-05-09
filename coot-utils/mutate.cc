
#include <string.h> // bleugh
#include "compat/coot-sysdep.h"
#include "geometry/main-chain.hh"
#include "coot-coord-utils.hh"

void
coot::util::mutate_internal(mmdb::Residue *residue,
                            mmdb::Residue *std_residue,
                            const std::string &alt_conf,
                            short int is_from_shelx_ins_flag,
                            float b_factor) {

   mmdb::PPAtom residue_atoms;
   int nResidueAtoms;
   residue->GetAtomTable(residue_atoms, nResidueAtoms);
   std::string res_name(residue->GetResName());

   mmdb::PPAtom std_residue_atoms;
   int n_std_ResidueAtoms;
   std_residue->GetAtomTable(std_residue_atoms, n_std_ResidueAtoms);

   //
   std::string old_seg_id_for_residue_atoms;
   bool use_old_seg_id = 0;
   try {
      old_seg_id_for_residue_atoms = coot::residue_atoms_segid(residue);
      use_old_seg_id = 1;
   }
   catch (const std::runtime_error &mess) {
   } 

   bool verb = false;
   if (verb) { 
      std::cout << "Mutate Atom Tables" << std::endl;
      std::cout << "Before " << residue_spec_t(residue) <<std::endl;
      for(int i=0; i<nResidueAtoms; i++)
         std::cout << residue_atoms[i]->name << std::endl;
      std::cout << "To be replaced by: " << residue_spec_t(std_residue) << std::endl;
      for(int i=0; i<n_std_ResidueAtoms; i++)
         std::cout << std_residue_atoms[i]->name << std::endl;
   }

   // only touch the atoms with given alt conf, ignore the others.
   std::string to_residue_type = std_residue->GetResName();
   for(int i=0; i<nResidueAtoms; i++) {
      std::string atom_alt_conf(residue_atoms[i]->altLoc);
      if (atom_alt_conf == alt_conf) { 
         std::string residue_this_atom (residue_atoms[i]->name);
         if (coot::is_main_chain_p(residue_atoms[i])) {
            if (to_residue_type == "MSE") {
               residue_atoms[i]->Het = 1;
            } else {
               if (residue_atoms[i]->Het)
                  residue_atoms[i]->Het = 0; // from MSE to MET, say
            }
            if (to_residue_type == "PRO") {
               std::string atom_name(residue_atoms[i]->name);
               if (atom_name == " H  ")   // PDBv3 FIXME
                  residue->DeleteAtom(i);
            }
            if (to_residue_type == "GLY") {
               std::string atom_name(residue_atoms[i]->name);
               if (atom_name == " HA ")  // PDBv3 FIXME
                  residue->DeleteAtom(i);
            }
         } else {
            // don't delete OXT, but do delete other things
            std::string atom_name(residue_atoms[i]->GetAtomName());
            if (atom_name != " OXT")
               residue->DeleteAtom(i);
         }
      }
   }

   for(int i=0; i<n_std_ResidueAtoms; i++) {
      std::string std_residue_this_atom (std_residue_atoms[i]->name);
      if (! coot::is_main_chain_p(std_residue_atoms[i])) {
         if (is_from_shelx_ins_flag)
            std_residue_atoms[i]->occupancy = 11.0;
         std_residue_atoms[i]->tempFactor = b_factor;
         mmdb::Atom *copy_at = new mmdb::Atom;
         copy_at->Copy(std_residue_atoms[i]);
         residue->AddAtom(copy_at);
         if (use_old_seg_id) {
            strcpy(copy_at->segID, old_seg_id_for_residue_atoms.c_str());
         }
         if (alt_conf != "")
            strcpy(copy_at->altLoc, alt_conf.c_str());
      }
   }

   residue->SetResName(std_residue->GetResName());
   residue->TrimAtomTable();
}

// return state, 0 bad, 1 good
// 
int
coot::util::mutate(mmdb::Residue *res, mmdb::Residue *std_res_unoriented, const std::string &alt_conf,
                   short int shelx_flag, float b_factor) {

   // If you are looking here then you have mutated a residue with a pointer. You need
   // to use deep_copy_this_residue(residue_p) not residue_p.

   int istate = 0; 
   std::map<std::string, clipper::RTop_orth> rtops = get_ori_to_this_res(res);

   mmdb::PPAtom residue_atoms;
   int nResidueAtoms;
   std_res_unoriented->GetAtomTable(residue_atoms, nResidueAtoms);
   if (nResidueAtoms == 0) {
      std::cout << " something broken in atom residue selection in ";
      std::cout << "mutate, got 0 atoms" << std::endl;
   } else {
      for(int iat=0; iat<nResidueAtoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         clipper::Coord_orth co(at->x, at->y, at->z);
         std::map<std::string, clipper::RTop_orth>::const_iterator it = rtops.find(alt_conf);
         if (it != rtops.end()) { 
            clipper::Coord_orth rotted = co.transform(it->second);

            at->x = rotted.x();
            at->y = rotted.y();
            at->z = rotted.z();
         }
      }

      coot::util::mutate_internal(res, std_res_unoriented, alt_conf, shelx_flag, b_factor); // it's oriented now.
      istate = 1;
   }
   return istate;
}




// Here std_base is at some arbitary position when passed.
//
// use_old_style_naming means use pdb v2 atom names.
//
// You need to call mol->FinishStructEdit() after using this function.
//
// LSQ fit std_base onto residue, then change the atoms of residue to
// be those of std_base.
// 
void
coot::util::mutate_base(mmdb::Residue *residue, mmdb::Residue *std_base,
                        bool use_old_style_naming,
                        bool print_match_stats_flag,
                        float b_factor) {


   bool debug = false;
   
   std::vector<std::string> adenine;  // Pyrimidine
   adenine.push_back(" N9 ");
   adenine.push_back(" C8 ");
   adenine.push_back(" N7 ");
   adenine.push_back(" C5 ");
   adenine.push_back(" C4 ");
   // 
   adenine.push_back(" N1 ");
   adenine.push_back(" C2 ");
   adenine.push_back(" N3 ");
   adenine.push_back(" C6 ");
   adenine.push_back(" N6 ");

   std::vector<std::string> guanine; // Pyrimidine
   guanine.push_back(" N9 ");
   guanine.push_back(" C8 ");
   guanine.push_back(" N7 ");
   guanine.push_back(" C5 ");
   guanine.push_back(" C4 ");
   //
   guanine.push_back(" N1 ");
   guanine.push_back(" C2 ");
   guanine.push_back(" N3 ");
   guanine.push_back(" C6 ");
   guanine.push_back(" O6 ");
   guanine.push_back(" N2 "); // No matcher for this in adenine

   std::vector<std::string> thymine;  // Purine
   thymine.push_back(" N1 ");
   thymine.push_back(" C2 ");
   thymine.push_back(" N3 ");
   thymine.push_back(" C4 ");
   thymine.push_back(" C5 ");
   thymine.push_back(" C6 ");
   // 
   thymine.push_back(" O2 ");
   thymine.push_back(" O4 ");
   thymine.push_back(" C7 "); // was C5M
   thymine.push_back(" C5M");
   
   std::vector<std::string> cytosine;  // Purine
   cytosine.push_back(" N1 ");
   cytosine.push_back(" C2 ");
   cytosine.push_back(" N3 ");
   cytosine.push_back(" C4 ");
   cytosine.push_back(" C5 ");
   cytosine.push_back(" C6 ");
   // 
   cytosine.push_back(" O2 ");
   cytosine.push_back(" N4 ");
   
   std::vector<std::string> uracil;  // Purine
   uracil.push_back(" N1 ");
   uracil.push_back(" C2 ");
   uracil.push_back(" N3 ");
   uracil.push_back(" C4 ");
   uracil.push_back(" C5 ");
   uracil.push_back(" C6 ");
   // 
   uracil.push_back(" O2 ");
   uracil.push_back(" O4 ");
   

   // These next 2 are in match order, don't change it.
   std::vector<std::string> purine; // A and G
   purine.push_back(" N9 ");
   purine.push_back(" C4 ");
   purine.push_back(" C5 ");
   purine.push_back(" N7 ");
   purine.push_back(" C8 ");

   std::vector<std::string> pyrimidine; // T, U and C
   pyrimidine.push_back(" N1 ");
   pyrimidine.push_back(" C2 ");
   pyrimidine.push_back(" N3 ");
   pyrimidine.push_back(" C5 ");
   pyrimidine.push_back(" C6 ");
   pyrimidine.push_back(" C4 ");


   std::string old_seg_id_for_residue_atoms;
   bool use_old_seg_id = 0;
   try {
      old_seg_id_for_residue_atoms = coot::residue_atoms_segid(residue);
      use_old_seg_id = 1;
   }
   catch (const std::runtime_error &mess) {
   } 


   // We need to know whether we have purine or pyrimidine for both
   // the molecule base and the std_base.
   //
   // We need to get all (5 for pyrimidine, 6 for purine) the
   // coordinates for both bases.
   // 
   // If they are either or both are pyrimidine we match 5 atoms,
   // If they are both purine we match 6 atoms.


   // So what are the input base types?
   //
   // These for flags should be set to something after our test
   short int mol_base_is_pyrimidine = -1;
   short int mol_base_is_purine     = -1;
   short int std_base_is_pyrimidine = -1;
   short int std_base_is_purine     = -1;

   std::string mol_base_name = residue->GetResName();
   std::string std_base_name = std_base->GetResName();

   if (mol_base_name == "Ar" || mol_base_name == "Ad" ||
       mol_base_name == "Gr" || mol_base_name == "Gd" ||
       mol_base_name == "G"  || mol_base_name == "DG" ||
       mol_base_name == "A"  || mol_base_name == "DA" ) {
      mol_base_is_purine = 1;
      mol_base_is_pyrimidine = 0;
   }

   if (mol_base_name == "Cr" || mol_base_name == "Cd" ||
       mol_base_name == "Ur" || mol_base_name == "Ud" ||
       mol_base_name == "Tr" || mol_base_name == "Td" ||
       mol_base_name == "U"  || mol_base_name == "DT" ||
       mol_base_name == "Ur" || mol_base_name == "T"  ||
       mol_base_name == "DC" || mol_base_name == "C") {
      mol_base_is_pyrimidine = 1;
      mol_base_is_purine = 0;
   }

   if (std_base_name == "Ar" || std_base_name == "Ad" ||
       std_base_name == "Gr" || std_base_name == "Gd" ||
       std_base_name == "G"  || std_base_name == "DG" ||
       std_base_name == "A"  || std_base_name == "DA" ) {
      std_base_is_purine = 1;
      std_base_is_pyrimidine = 0;
   }

   if (std_base_name == "Cr" || std_base_name == "Cd" ||
       std_base_name == "Tr" || std_base_name == "Td" ||
       std_base_name == "Ur" || std_base_name == "Ud" ||
       std_base_name == "U"  || std_base_name == "DT" ||
       std_base_name == "T"  || std_base_name == "DC") { 
      std_base_is_pyrimidine = 1;
      std_base_is_purine = 0;
   }

   if ((mol_base_is_pyrimidine == -1) || (mol_base_is_purine == -1) || 
       (std_base_is_pyrimidine == -1) || (std_base_is_purine == -1) ) {

      std::cout << "ERROR:: mutate_base() unassigned type "
                << "mol_base_is_pyrimidine:" << " "
                << mol_base_is_pyrimidine << " "
                << "mol_base_is_purine: " << " "
                << mol_base_is_purine << " "
                << "std_base_is_pyrimidine: " << " "
                << std_base_is_pyrimidine << " "
                << "std_base_is_purine: " << " "
                << std_base_is_purine << " "
                << residue->GetResName() << " " << std_base->GetResName()
                << std::endl;

   } else {

      if (debug)
         std::cout << "DEBUG:: assigned types \n      "
                   << "mol_base_is_pyrimidine:"  << " "
                   << mol_base_is_pyrimidine     << "\n      "
                   << "mol_base_is_purine:    "  << " "
                   << mol_base_is_purine         << "\n      "
                   << "std_base_is_pyrimidine:" << " "
                   << std_base_is_pyrimidine     << "\n      "
                   << "std_base_is_purine:    " << " "
                   << std_base_is_purine         << "\n      res_name: "
                   << residue->GetResName()      << "    std_base name:"
                   << std_base->GetResName()
                   << std::endl;

      int n_match_atoms = 5;
      if (mol_base_is_pyrimidine && std_base_is_pyrimidine)
         n_match_atoms = 6;

      std::vector<std::string> moving_name_vector;
      std::vector<std::string> refrce_name_vector;

      if (std_base_is_purine)
         moving_name_vector = purine;
      else
         moving_name_vector = pyrimidine;

      if (mol_base_is_purine)
         refrce_name_vector = purine;
      else
         refrce_name_vector = pyrimidine;
      
      mmdb::PAtom *std_base_atoms;
      int n_std_base_atoms;

      mmdb::PAtom *mol_base_atoms;
      int n_mol_base_atoms;
      
      residue->GetAtomTable( mol_base_atoms, n_mol_base_atoms);
      std_base->GetAtomTable(std_base_atoms, n_std_base_atoms);

      std::vector<clipper::Coord_orth> refrce_atom_positions;
      std::vector<clipper::Coord_orth> moving_atom_positions;

      if (debug) { 
         for (unsigned int i=0; i<refrce_name_vector.size(); i++)
            std::cout << "ref base search atom :" << refrce_name_vector[i]
                      << ":" << std::endl;
         for (unsigned int i=0; i<moving_name_vector.size(); i++)
            std::cout << "mov base search atom :" << moving_name_vector[i]
                      << ":" << std::endl;
      }
         
      for (int j=0; j<n_match_atoms; j++) {
         for (int i=0; i<n_mol_base_atoms; i++) {
            std::string atom_name = mol_base_atoms[i]->name;
            if (refrce_name_vector[j] == atom_name) {
               refrce_atom_positions.push_back(clipper::Coord_orth(mol_base_atoms[i]->x,
                                                                   mol_base_atoms[i]->y,
                                                                   mol_base_atoms[i]->z));
               if (debug)
                  std::cout << "Found " << atom_name << " in reference " << std::endl;
            }
         }
      }

      for (int j=0; j<n_match_atoms; j++) {
         for (int i=0; i<n_std_base_atoms; i++) {
         std::string atom_name = std_base_atoms[i]->name;
            if (moving_name_vector[j] == atom_name) {
               moving_atom_positions.push_back(clipper::Coord_orth(std_base_atoms[i]->x,
                                                                   std_base_atoms[i]->y,
                                                                   std_base_atoms[i]->z));
               if (debug)
                  std::cout << "Found " << atom_name << " in moving (std) base " << std::endl;
            }
         }
      }

      if (int(refrce_atom_positions.size()) != n_match_atoms) {
         std::cout << "ERROR:: wrong number of reference atoms found! "
                   << refrce_atom_positions.size() << std::endl;
      } else {

         if (int(moving_atom_positions.size()) != n_match_atoms) {
            std::cout << "ERROR:: wrong number of moving atoms found! "
                   << moving_atom_positions.size() << std::endl;

         } else {

            clipper::RTop_orth rtop (moving_atom_positions, refrce_atom_positions);
      
            double sum_dist = 0.0;
            double sum_dist2 = 0.0;
            double mind =  999999999.9;
            double maxd = -999999999.9;
            double d;
            for (unsigned int i=0; i<refrce_atom_positions.size(); i++) {
               d = clipper::Coord_orth::length(refrce_atom_positions[i],
                                               clipper::Coord_orth(moving_atom_positions[i].transform(rtop)));
               sum_dist  += d;
               sum_dist2 += d*d;
               if (d>maxd)
                  maxd = d;
               if (d<mind)
                  mind = d;
            }
            double mean = sum_dist/double(moving_atom_positions.size());
            double var  = sum_dist2/double(moving_atom_positions.size()); // no mean*mean
            std::cout << "INFO:: " << moving_atom_positions.size() << " matched atoms had: \n"
                      << "   mean devi: " << mean << "\n"
                      << "    rms devi: " << sqrt(var) << "\n"
                      << "    max devi: " << maxd << "\n"
                      << "    min devi: " << mind << std::endl;

            mmdb::Atom *at;
            // We are going to delete the current atoms of the residue
            // and add the std_base ones.  First what *are* the atom
            // names we what to add or delete?
            // 
            std::vector<std::string> mol_base_atom_names;
            if (mol_base_name == "Ar" || mol_base_name == "Ad")
               mol_base_atom_names = adenine;
            if (mol_base_name == "Gr" || mol_base_name == "Gd")
               mol_base_atom_names = guanine;
            if (mol_base_name == "Cr" || mol_base_name == "Cd")
               mol_base_atom_names = cytosine;
            if (mol_base_name == "Tr" || mol_base_name == "Td")
               mol_base_atom_names = thymine;
            if (mol_base_name == "Ur" || mol_base_name == "Ud")
               mol_base_atom_names = uracil;

            // new names 
            if (mol_base_name == "A" || mol_base_name == "DA")
               mol_base_atom_names = adenine;
            if (mol_base_name == "G" || mol_base_name == "DG")
               mol_base_atom_names = guanine;
            if (mol_base_name == "C" || mol_base_name == "DC")
               mol_base_atom_names = cytosine;
            if (mol_base_name == "T" || mol_base_name == "DT")
               mol_base_atom_names = thymine;
            if (mol_base_name == "U" || mol_base_name == "DU")
               mol_base_atom_names = uracil;

            if (mol_base_atom_names.size() == 0) {
               std::cout << "ERROR:: muate_base(): ";
               std::cout << "failed to find mol_base_name for mol_base_atom_names\n";
            } else {
               
               std::vector<std::string> std_base_atom_names;

               if (std_base_name == "Ar" || std_base_name == "Ad")
                  std_base_atom_names = adenine;
               if (std_base_name == "Gr" || std_base_name == "Gd")
                  std_base_atom_names = guanine;
               if (std_base_name == "Cr" || std_base_name == "Cd")
                  std_base_atom_names = cytosine;
               if (std_base_name == "Tr" || std_base_name == "Td")
                  std_base_atom_names = thymine;
               if (std_base_name == "Ur" || std_base_name == "Ud")
                  std_base_atom_names = uracil;

               // new names
               if (std_base_name == "A" || std_base_name == "DA")
                  std_base_atom_names = adenine;
               if (std_base_name == "G" || std_base_name == "DG")
                  std_base_atom_names = guanine;
               if (std_base_name == "C" || std_base_name == "DC")
                  std_base_atom_names = cytosine;
               if (std_base_name == "T" || std_base_name == "DT")
                  std_base_atom_names = thymine;
               if (std_base_name == "U" || std_base_name == "DU")
                  std_base_atom_names = uracil;

            
               if (std_base_atom_names.size() == 0) {
                  std::cout << "ERROR:: muate_base(): ";
                  std::cout << "failed to find std_base_name for std_base_atom_names\n";
               } else {

                  // now find the atoms of the given residue and apply
                  // the transformation and add them to the residue;
                  
                  bool have_deleted = 0;
                  bool need_a_ter = 0;
                  for (unsigned int iat=0; iat<mol_base_atom_names.size(); iat++) {
                     for (int i=0; i<n_mol_base_atoms; i++) {
                        if (mol_base_atoms[i]) { 
                        
                           if (mol_base_atom_names[iat] == mol_base_atoms[i]->name) {

                              if (debug)
                                 std::cout << ".... Deleting Atom " << mol_base_atoms[i]->name
                                           << " i = " << i << std::endl;
                              
                              residue->DeleteAtom(i);
                              mol_base_atoms[i] = NULL;
                              have_deleted = 1;
                              break;
                           }
                        }
                     }
                  }
                   if (have_deleted)
                      residue->TrimAtomTable();


                  mmdb::PPAtom residue_atoms = 0;
                  int n_residue_atoms;
                  residue->GetAtomTable(residue_atoms, n_residue_atoms);
                  have_deleted = 0; // reset for Ter tests.

                  for (int i=0; i<n_residue_atoms; i++) { 
                     if (residue_atoms[i]->isTer()) {
                        if (debug)
                           std::cout << "..... atom " << i << " is a Ter" << std::endl;
                        residue->DeleteAtom(i);
                        need_a_ter = 1;
                        have_deleted = 1;
                     }
                  }
                   if (have_deleted)
                      residue->TrimAtomTable();
                  
                  for (unsigned int iat=0; iat<std_base_atom_names.size(); iat++) {
                     bool found = 0;
                     for (int i=0; i<n_std_base_atoms; i++) {
                        if (std_base_atom_names[iat] == std_base_atoms[i]->name) {
                           clipper::Coord_orth p(std_base_atoms[i]->x,
                                                 std_base_atoms[i]->y,
                                                 std_base_atoms[i]->z);
                           clipper::Coord_orth pt = p.transform(rtop);
                           std::string ele = std_base_atom_names[iat].substr(0,2);
                           at = new mmdb::Atom;
                           if (debug)
                              std::cout << ".... Adding Atom " << std_base_atoms[i]->name
                                        << std::endl;
                           at->SetCoordinates(pt.x(), pt.y(), pt.z(), 1.0, b_factor);
                           std::string new_atom_name = std_base_atoms[i]->name;
                           if (std_base_name == "Td")
                              if (new_atom_name == " C5M")
                                 if (! use_old_style_naming)
                                    new_atom_name = " C7 ";
                           at->SetAtomName(new_atom_name.c_str());
                           at->SetElementName(ele.c_str());
                           std::string new_alt_conf("");
                           // force it down the atom's throat :) [is there a better way?]
                           strncpy(at->altLoc, new_alt_conf.c_str(), 2);
                           residue->AddAtom(at);
                           if (use_old_seg_id)
                              strcpy(at->segID, old_seg_id_for_residue_atoms.c_str());
                           found = 1;
                           break;
                        }
                     }
                     if (! found) {
                        // std::cout << "... failed to find std_base_atom  "
                        // << std_base_atom_names[iat] << std::endl;
                     } 
                  }
            

                   std::string new_base_name =
                      coot::util::convert_base_name(std_base_name, use_old_style_naming);
                  residue->SetResName(new_base_name.c_str());
                  residue->TrimAtomTable();


                  if (need_a_ter) {
                     mmdb::Atom *at = new mmdb::Atom;
                     at->MakeTer();
                     residue->TrimAtomTable();
                  }
               }
            }
         }
      }
   }
}

std::string
coot::util::convert_base_name(const std::string &std_base_name, bool use_old_style_naming) {

   if (use_old_style_naming) { 
      return std_base_name;
   } else {
      if (std_base_name == "Cd")
         return "DC";
      if (std_base_name == "Ad")
         return "DA";
      if (std_base_name == "Gd")
         return "DG";
      if (std_base_name == "Td")
         return "DT";
      if (std_base_name == "Cr")
         return "C";
      if (std_base_name == "Ar")
         return "A";
      if (std_base_name == "Gr")
         return "G";
      if (std_base_name == "Ur")
         return "U";
      if (std_base_name == "Tr")
         return "T";
   }

   return std_base_name;
} 


int
coot::util::mutate(mmdb::Residue *residue_p, const std::string &new_res_name) {

   mmdb::Residue *std_residue = get_standard_residue_instance(new_res_name); // a deep copy
   if (std_residue) {
      std::string alt_conf("");
      return mutate(residue_p, std_residue, alt_conf, 0, 30.0);
   } else {
      std::cout << "ERROR:: when retriving standard residue for type \""
                << new_res_name << "\"" << std::endl;
      return -1;
   }

}

#include "atom-selection-container.hh"
#include "utils/coot-utils.hh"

 // a deep copy
mmdb::Residue *
coot::util::get_standard_residue_instance(const std::string &res_name) {


   if (standard_residues_map_cache.empty()) {
      // fill it
      std::string dir = package_data_dir();
      std::string standard_residues_file_name = append_dir_file(dir, "standard-residues.pdb");
      atom_selection_container_t asc = get_atom_selection(standard_residues_file_name, false, false, false);

      if (asc.read_success) {
         int imod = 1;
         mmdb::Model *model_p = asc.mol->GetModel(imod);
         if (model_p) {
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               int nres = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<nres; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  std::string res_name(residue_p->GetResName());
                  standard_residues_map_cache[res_name] = residue_p;
               }
            }
         }
      } else {
         std::cout << "ERROR:: failed to read " << standard_residues_file_name << std::endl;
      }
   }

   if (! standard_residues_map_cache.empty()) {
      std::map<std::string, mmdb::Residue *>::const_iterator it;
      it = standard_residues_map_cache.find(res_name);
      if (it != standard_residues_map_cache.end()) {
         mmdb::Residue *res = it->second;
         return deep_copy_this_residue(res);
      }
   }
   return 0; // failed, then
}
