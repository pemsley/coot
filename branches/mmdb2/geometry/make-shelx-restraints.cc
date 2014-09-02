
#include <fstream>
#include <iomanip>
#include "utils/coot-utils.hh"
#include "geometry/protein-geometry.hh"

void output_flats(const coot::dictionary_residue_restraints_t &rest,
		  std::ofstream &f) {

   for (unsigned int i=0; i<rest.plane_restraint.size(); i++) { 
      const coot::dict_plane_restraint_t &pr= rest.plane_restraint[i];
      f << "FLAT_" << rest.residue_info.comp_id;
      for (unsigned int iat=0; iat<pr.n_atoms(); iat++)
	 f << " " << coot::util::remove_whitespace(pr.atom_id(iat));
      f << "\n";
   }
} 

void output_chivs(const coot::dictionary_residue_restraints_t &rest,
		  std::ofstream &f) {

   if (rest.chiral_restraint.size() > 0) { 
      f << "CHIV_" << rest.residue_info.comp_id;
      for (unsigned int i=0; i<rest.chiral_restraint.size(); i++) { 
	 const coot::dict_chiral_restraint_t &r= rest.chiral_restraint[i];
	 f << " " << coot::util::remove_whitespace(r.atom_id_c_4c());
      }
      f << "\n";
   }
}

void output_dfixs(const coot::dictionary_residue_restraints_t &rest,
		  std::ofstream &f) {

   for (unsigned int i=0; i<rest.bond_restraint.size(); i++) { 
      const coot::dict_bond_restraint_t &r = rest.bond_restraint[i];
      f << "DFIX_" << rest.residue_info.comp_id;
      f << " " << std::setprecision(4) << r.value_dist();
      f << " " << coot::util::remove_whitespace(r.atom_id_1());
      f << " " << coot::util::remove_whitespace(r.atom_id_2());
      f << "\n";
   }
}

void output_dangs(const coot::dictionary_residue_restraints_t &rest,
		  std::ofstream &f) {

   // we have to convert angles to 1-3 distances.  Need to find bonds
   // of the angle atoms.
   
   for (unsigned int i=0; i<rest.angle_restraint.size(); i++) { 
      const coot::dict_angle_restraint_t &r = rest.angle_restraint[i];
      if (! rest.is_hydrogen(r.atom_id_1())) { 
	 if (! rest.is_hydrogen(r.atom_id_3())) { 
	    coot::dict_bond_restraint_t b1(r.atom_id_1(), r.atom_id_2(), "", 0.0, 0.0);
	    coot::dict_bond_restraint_t b2(r.atom_id_2(), r.atom_id_3(), "", 0.0, 0.0);
	    for (unsigned int ib1=0; ib1<rest.bond_restraint.size(); ib1++) {
	       if (rest.bond_restraint[ib1].matches_names(b1)) {
		  for (unsigned int ib2=0; ib2<rest.bond_restraint.size(); ib2++) {
		     if (rest.bond_restraint[ib2].matches_names(b2)) {
			double d1 = rest.bond_restraint[ib1].value_dist();
			double d2 = rest.bond_restraint[ib2].value_dist();
			double angle = clipper::Util::d2rad(r.angle());
			double d_sqrd = d1*d1 + d2*d2 - 2*d1*d2*cos(angle);
			if (d_sqrd > 0) {
			   double d = sqrt(d_sqrd);
			   f << "DANG_" << rest.residue_info.comp_id;
			   f << " " << std::setprecision(4) << d;
			   f << " " << coot::util::remove_whitespace(r.atom_id_1());
			   f << " " << coot::util::remove_whitespace(r.atom_id_3());
			   f << "\n";
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
}
   

int main(int argc, char **argv) {

   int status = 0;

   if (argc>1) {
      std::string mmcif_file_name(argv[1]);
      coot::protein_geometry geom;
      geom.set_verbose(false);
      int read_number = 0;
      geom.init_refmac_mon_lib(mmcif_file_name, read_number);
      std::vector<std::string> types = geom.monomer_types();

      if (types.size() > 0) {
	 
	 std::string file_name = "shelx-restraints.txt";
	 if (types.size() == 1)
	    file_name = "shelx-restraints-" + types[0] + ".txt";
	 std::ofstream f(file_name.c_str());

	 for (unsigned int i=0; i<types.size(); i++) { 
	    std::pair<bool, coot::dictionary_residue_restraints_t> r =
	       geom.get_monomer_restraints(types[0]);
	    if (r.first) {

	       // how can it not be?
	       const coot::dictionary_residue_restraints_t &rest = r.second;
	       // FLAT CHIV DFIX DANG

	       f << "\n";

	       output_flats(rest, f);
	       output_chivs(rest, f);
	       output_dfixs(rest, f);
	       output_dangs(rest, f);
	    }
	 }
      }
   }
   return status;
}
