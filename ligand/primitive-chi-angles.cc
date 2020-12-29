
#include <iostream>
#include <stdexcept>

#include "geometry/mol-utils.hh"
#include "coot-utils/coot-coord-utils.hh"

#include "primitive-chi-angles.hh"



std::vector<coot::alt_confed_chi_angles>
coot::primitive_chi_angles::get_chi_angles() {

   std::vector<coot::alt_confed_chi_angles> nv;
   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   residue->GetAtomTable(residue_atoms, n_residue_atoms);

   bool residue_has_alt_confs = 0;
   for (int i=0; i<n_residue_atoms; i++) {
      std::string alt_conf(residue_atoms[i]->altLoc);
      if (alt_conf != "") {
         residue_has_alt_confs = 1;
         break;
      }
   }

   std::vector<coot::atom_name_quad> atom_name_quad_list = get_atom_name_quads();
   if (atom_name_quad_list.size() == 0) {
      std::string mess = "Failed to find atom name quads for residue type ";
      mess += residue_name;
      throw std::runtime_error(mess);
   }
   if (residue_has_alt_confs == 0) {
      std::vector<coot::atom_index_quad> quads = get_quads(atom_name_quad_list, residue);
      // std::cout << "quads.size() " << quads.size() << std::endl;
      std::vector<std::pair<int, float> > v;
      for (unsigned int i_quad = 0; i_quad<quads.size(); i_quad++) {
         clipper::Coord_orth p1(atom_to_co(residue_atoms[quads[i_quad].index1]));
         clipper::Coord_orth p2(atom_to_co(residue_atoms[quads[i_quad].index2]));
         clipper::Coord_orth p3(atom_to_co(residue_atoms[quads[i_quad].index3]));
         clipper::Coord_orth p4(atom_to_co(residue_atoms[quads[i_quad].index4]));
         double tors =  clipper::Util::rad2d(clipper::Coord_orth::torsion(p1, p2, p3, p4));
         // std::cout << "pushing back plane residue torsion" << std::endl;
         v.push_back(std::pair<int,float>(i_quad+1, tors));
      }
      if (v.size() > 0) {
         nv.push_back(coot::alt_confed_chi_angles("", v));
      } else {
         if (false) // too verbose output
            std::cout << "INFO:: un-altconfed rotamer with no chis "
                      << coot::residue_spec_t(residue) << std::endl;
      }
   } else {
      // multiple quads from alt confed atoms:
      std::vector<coot::alt_confed_atom_index_quad> quads_vec =
         get_quads_using_altconfs(atom_name_quad_list, residue);
      // std::cout << "quads_vec.size() " << quads_vec.size() << std::endl;
//       for (unsigned int i_quad_set = 0; i_quad_set<quads_vec.size(); i_quad_set++) {
//          std::cout << "quads_vec[" << i_quad_set << "].size() " << quads_vec[i_quad_set].quad.size()
//                    << std::endl;
//       }
      for (unsigned int i_quad_set = 0; i_quad_set<quads_vec.size(); i_quad_set++) {
         std::vector<std::pair<int, float> > v;
         std::vector<coot::atom_index_quad> quad = quads_vec[i_quad_set].quad;
         for (unsigned int i_quad = 0; i_quad<quad.size(); i_quad++) {
            clipper::Coord_orth p1(atom_to_co(residue_atoms[quad[i_quad].index1]));
            clipper::Coord_orth p2(atom_to_co(residue_atoms[quad[i_quad].index2]));
            clipper::Coord_orth p3(atom_to_co(residue_atoms[quad[i_quad].index3]));
            clipper::Coord_orth p4(atom_to_co(residue_atoms[quad[i_quad].index4]));
            double tors =  clipper::Util::rad2d(clipper::Coord_orth::torsion(p1, p2, p3, p4));
            v.push_back(std::pair<int,float>(i_quad+1, tors));
         }
         if (v.size() > 0)
            nv.push_back(coot::alt_confed_chi_angles(quads_vec[i_quad_set].alt_conf, v));
         else
            std::cout << "not pushing back altconfed altconfed rotamer with no chis" << std::endl;
      }
   }

   if (0) {  // debugging 
      for (unsigned int inv=0; inv<nv.size(); inv++) {
	 std::cout << "DEBUG:: primitive_chi_angles:: nv index: " << inv << " chi_angles.size() "
		   << nv[inv].chi_angles.size() << std::endl;
	 for (unsigned int ichi=0; ichi<nv[inv].chi_angles.size(); ichi++) {
	    std::cout << "DEBUG:: primitive_chi_angles:: get_chi_angles "
		      << coot::residue_spec_t(residue)
		      << " chi set: "
		      << inv << " alt conf :" << nv[inv].alt_conf << ": "
		      << nv[inv].chi_angles[ichi].first
		      << " " << nv[inv].chi_angles[ichi].second << "\n" ;
	 }
      }
   }
   return nv;

}


void
coot::primitive_chi_angles::setup_chi_atom_quads() {

   add_chi_quad("VAL", " N  ", " CA ", " CB ", " CG1");

   add_chi_quad("TYR", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("TYR", " CA ", " CB ", " CG ", " CD1");

   add_chi_quad("TRP", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("TRP", " CA ", " CB ", " CG ", " CD1");

   add_chi_quad("THR", " N  ", " CA ", " CB ", " OG1");

   add_chi_quad("SER", " N  ", " CA ", " CB ", " OG ");

   add_chi_quad("PRO", " N  ", " CA ", " CB ", " CG ");
   //    add_chi_quad("PRO", " CA ", " CB ", " CG ", " CD ");  Hmmmm!

   add_chi_quad("PHE", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("PHE", " CA ", " CB ", " CG ", " CD1");
   
   add_chi_quad("MET", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("MET", " CA ", " CB ", " CG ", " SD ");
   add_chi_quad("MET", " CB ", " CG ", " SD ", " CE ");

   // Use the monomer dictionary bonding instead   
   add_chi_quad("MSE", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("MSE", " CA ", " CB ", " CG ", "SE  ");
   add_chi_quad("MSE", " CB ", " CG ", "SE  ", " CE ");

   add_chi_quad("LYS", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("LYS", " CA ", " CB ", " CG ", " CD ");
   add_chi_quad("LYS", " CB ", " CG ", " CD ", " CE ");
   add_chi_quad("LYS", " CG ", " CD ", " CE ", " NZ ");

   add_chi_quad("LEU", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("LEU", " CA ", " CB ", " CG ", " CD1");

   add_chi_quad("ILE", " N  ", " CA ", " CB ", " CG1");
   add_chi_quad("ILE", " CA ", " CB ", " CG1", " CD1");

   add_chi_quad("HIS", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("HIS", " CA ", " CB ", " CG ", " ND1");

   add_chi_quad("GLU", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("GLU", " CA ", " CB ", " CG ", " CD ");
   add_chi_quad("GLU", " CB ", " CG ", " CD ", " OE1");

   add_chi_quad("GLN", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("GLN", " CA ", " CB ", " CG ", " CD ");
   add_chi_quad("GLN", " CB ", " CG ", " CD ", " OE1");

   add_chi_quad("CYS", " N  ", " CA ", " CB ", " SG ");
   
   add_chi_quad("ASP", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("ASP", " CA ", " CB ", " CG ", " OD1");

   add_chi_quad("ASN", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("ASN", " CA ", " CB ", " CG ", " OD1");

   add_chi_quad("ARG", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("ARG", " CA ", " CB ", " CG ", " CD ");
   add_chi_quad("ARG", " CB ", " CG ", " CD ", " NE ");
   add_chi_quad("ARG", " CG ", " CD ", " NE ", " CZ ");

}

// Return chi angle pair for the chi angles in setup_chi_atom_pairs
// (above).  Not used by generic ligands (as far as I can see).
void
coot::primitive_chi_angles::add_chi_quad(const std::string &residue_type,
					 const std::string &atom_name_1,
					 const std::string &atom_name_2,
					 const std::string &atom_name_3,
					 const std::string &atom_name_4) {

   bool found_res = 0;
   for(unsigned int i=0; i< chi_angle_atoms_for_residue_type.size(); i++) {
      if (chi_angle_atoms_for_residue_type[i].residue_type == residue_type) {
         found_res = 1;
         chi_angle_atoms_for_residue_type[i].add_torsion_bond_by_name(atom_name_1, atom_name_2,
								      atom_name_3, atom_name_4);
	 break;
      }
   }

   // This is a departure from other chi_angle/monomer-utils code, at
   // startup of a primitive_chi_angles,
   // chi_angle_atoms_for_residue_type is empty so if we dont so this
   // block the we get the "Ooops ARG not found in
   // chi_angle_atoms_for_residue_type" message , but for chi_angles,
   // that does not happen.
   //
   if (! found_res) {
      coot::atom_name_quad quad(atom_name_1, atom_name_2, atom_name_3, atom_name_4);
      coot::residue_named_chi_angle_atom_name_quad_set_t named_quad(residue_type,quad);
      chi_angle_atoms_for_residue_type.push_back(named_quad);
   } 

//    if (found_res == 0) {
//       std::cout << "Oops, " << residue_type << " not found in chi_angle_atoms_for_residue_type"
//                 << std::endl;
//    }
}


void
coot::residue_named_chi_angle_atom_name_quad_set_t::add_torsion_bond_by_name(const std::string &atom_name_1,
                                                                             const std::string &atom_name_2,
                                                                             const std::string &atom_name_3,
                                                                             const std::string &atom_name_4) {

   name_quad.push_back(coot::atom_name_quad(atom_name_1,
                                            atom_name_2,
                                            atom_name_3,
                                            atom_name_4));
}



std::vector<coot::atom_name_quad>
coot::primitive_chi_angles::get_atom_name_quads() const {
   std::vector<coot::atom_name_quad> v;

   for (unsigned int ir=0; ir<chi_angle_atoms_for_residue_type.size(); ir++) {
      if (chi_angle_atoms_for_residue_type[ir].residue_type == residue_name) {
         v = chi_angle_atoms_for_residue_type[ir].name_quad;
         break;
      }
   }
   return v;
}

std::vector<coot::atom_index_quad>
coot::primitive_chi_angles::get_quads(const std::vector<coot::atom_name_quad> &atom_name_quads,
				      mmdb::Residue *residue) const {
   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   residue->GetAtomTable(residue_atoms, n_residue_atoms);
   return get_atom_index_quads(atom_name_quads, residue_atoms, n_residue_atoms);
}


std::vector<coot::atom_index_quad>
coot::primitive_chi_angles::get_atom_index_quads(const std::vector<coot::atom_name_quad> &atom_name_quads_in,
                                          const mmdb::PPAtom atoms, int nresatoms) const {

   std::vector<coot::atom_index_quad> v;
   for (unsigned int iquad=0; iquad<atom_name_quads_in.size(); iquad++) {
      int ifound = 0;
      for (int i1=0; i1<nresatoms; i1++) {
         std::string atom_name = atoms[i1]->name;
         if (atom_name == atom_name_quads_in[iquad].atom_name(0)) {
            for (int i2=0; i2<nresatoms; i2++) {
               std::string atom_name = atoms[i2]->name;
               if (atom_name == atom_name_quads_in[iquad].atom_name(1)) {
                  for (int i3=0; i3<nresatoms; i3++) {
                     std::string atom_name = atoms[i3]->name;
                     if (atom_name == atom_name_quads_in[iquad].atom_name(2)) {
                        for (int i4=0; i4<nresatoms; i4++) {
                           std::string atom_name = atoms[i4]->name;
                           if (atom_name == atom_name_quads_in[iquad].atom_name(3)) {
                              v.push_back(coot::atom_index_quad(i1, i2, i3, i4));
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   if (v.size() < atom_name_quads_in.size()) {
      bool write_status_of_missing_quads = false;
      if (write_status_of_missing_quads) {
         std::cout << "primitive chis: Failure to find correct atom quads in residue atoms\n" ;
         for (unsigned int iquad=0; iquad<atom_name_quads_in.size(); iquad++) {
            std::cout << "  quad needed: :"
                      << atom_name_quads_in[iquad].atom_name(0) << ":  :"
                      << atom_name_quads_in[iquad].atom_name(1) << ":  :"
                      << atom_name_quads_in[iquad].atom_name(2) << ":  :"
                      << atom_name_quads_in[iquad].atom_name(3) << ":\n";
         }
         for (unsigned int iv=0; iv<v.size(); iv++) {
            std::cout << "  found quad: "
                      << v[iv].index1 << "  "
                      << v[iv].index2 << "  "
                      << v[iv].index3 << "  "
                      << v[iv].index4 << "\n";
         }
      }
   } else {
      // std::cout << "found all quads in residue atoms\n" ;
   }
   if (0) { // debugging
      std::cout << "DEBUG:: get_atom_index_quads returns v with size " << v.size() << " from "
		<< "atom_name_quads_in size " << atom_name_quads_in.size() << std::endl;
   }
   return v;
}


std::vector<coot::alt_confed_atom_index_quad>
coot::primitive_chi_angles::get_quads_using_altconfs(const std::vector<coot::atom_name_quad> &atom_name_quads,
                                              mmdb::Residue *residue) const {

   std::vector<coot::alt_confed_atom_index_quad> alt_v;
   mmdb::PPAtom atoms = 0;
   int n_residue_atoms;
   residue->GetAtomTable(atoms, n_residue_atoms);

   std::vector<std::string> residue_alt_confs = coot::util::get_residue_alt_confs(residue);
   // remove "" from the residue_alt_confs:
   std::vector<std::string>::iterator it;
   for (it = residue_alt_confs.begin(); it != residue_alt_confs.end(); it++) {
      if (*it == "") {
         residue_alt_confs.erase(it);
         break;
      }
   }

   if (0) { 
      std::cout << "Checking residue " << coot::residue_spec_t(residue) <<  " against these alt confs: "
		<< std::endl;
      for (it = residue_alt_confs.begin(); it != residue_alt_confs.end(); it++) {
	 std::cout << *it << std::endl;
      }
   }

   for (unsigned int i_alt_conf=0; i_alt_conf<residue_alt_confs.size(); i_alt_conf++) {
      std::vector<coot::atom_index_quad> v;
      for (unsigned int iquad=0; iquad<atom_name_quads.size(); iquad++) {
         for (int i1=0; i1<n_residue_atoms; i1++) {
            std::string atom_name = atoms[i1]->name;
            std::string alt_conf_1 = atoms[i1]->altLoc;
            if (atom_name == atom_name_quads[iquad].atom_name(0)) {
               for (int i2=0; i2<n_residue_atoms; i2++) {
                  std::string atom_name = atoms[i2]->name;
                  std::string alt_conf_2 = atoms[i2]->altLoc;
                  if (atom_name == atom_name_quads[iquad].atom_name(1)) {
                     for (int i3=0; i3<n_residue_atoms; i3++) {
                        std::string atom_name = atoms[i3]->name;
                        std::string alt_conf_3 = atoms[i3]->altLoc;
                        if (atom_name == atom_name_quads[iquad].atom_name(2)) {
                           for (int i4=0; i4<n_residue_atoms; i4++) {
                              std::string atom_name = atoms[i4]->name;
                              std::string alt_conf_4 = atoms[i4]->altLoc;
                              if (atom_name == atom_name_quads[iquad].atom_name(3)) {
                                 if (alt_conf_4 == residue_alt_confs[i_alt_conf] || alt_conf_4 == "") {
                                    if (alt_conf_3 == residue_alt_confs[i_alt_conf] || alt_conf_3 == "") {
                                       if (alt_conf_2 == residue_alt_confs[i_alt_conf] || alt_conf_2 == "") {
                                          if (alt_conf_1 == residue_alt_confs[i_alt_conf] || alt_conf_1 == "") {
                                             v.push_back(coot::atom_index_quad(i1, i2, i3, i4));
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
      if (v.size()) {
         alt_v.push_back(coot::alt_confed_atom_index_quad(residue_alt_confs[i_alt_conf], v));
      }
   }

   if (0) { // debugging
      std::cout << " DEBUG:: returning " << alt_v.size() << " chi sets from "
		<< coot::residue_spec_t(residue) << std::endl;
      for (unsigned int i_alt_v=0; i_alt_v<alt_v.size(); i_alt_v++) {
	 for (unsigned int ichi=0; ichi<alt_v[i_alt_v].quad.size(); ichi++) {
	    std::cout << "DEBUG:: primitive_chi_angles:: quads_using_altconf "
		      << coot::residue_spec_t(residue)
		      << " chi set: "
		      << i_alt_v << " alt conf :" << alt_v[i_alt_v].alt_conf << ": "
		      << alt_v[i_alt_v].quad[ichi].index1 << " "
		      << alt_v[i_alt_v].quad[ichi].index2 << " "
		      << alt_v[i_alt_v].quad[ichi].index3 << " "
		      << alt_v[i_alt_v].quad[ichi].index4 << " "
		      << "\n" ;
	 }

      }
   }
   return alt_v;
}

clipper::Coord_orth
coot::primitive_chi_angles::atom_to_co(mmdb::Atom *at) const {
   return clipper::Coord_orth(at->x, at->y, at->z);
}
