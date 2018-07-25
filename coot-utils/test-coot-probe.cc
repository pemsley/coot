
#include <stdlib.h>
#include "atom-overlaps.hh"
#include "coot-coord-utils.hh"
#include "geometry/residue-and-atom-specs.hh"

int main(int argc, char **argv) {

   int status = 0;
   
   // do we want to probe a ligand or all atoms?
   mmdb::Residue *residue_p = 0; // set when specified on command line

   if (argc < 2) {
      std::cout << "Usage: test-coot-probe <filename>" << std::endl;
      exit(0);
   } else {

      coot::protein_geometry geom;
      geom.init_standard();
      
      // load all the codes in the input molecule that not already loaded 
      mmdb::Manager *mol = new mmdb::Manager;

      std::string file_name = argv[1];
      int read_status = mol->ReadCoorFile(file_name.c_str());
      coot::residue_spec_t spec;

      if (read_status == mmdb::Error_NoError) {

	 if (! spec.empty()) {
	    mmdb::Residue *residue_p = coot::util::get_residue(spec, mol);
	    if (residue_p) {
	       std::vector<mmdb::Residue *> neighbs = coot::residues_near_residue(residue_p, mol, 5);
	       try {
		  coot::atom_overlaps_container_t overlaps(residue_p, neighbs, mol, &geom, 0.5, 0.25);
		  coot::atom_overlaps_dots_container_t c = overlaps.contact_dots_for_ligand();
	       }
	       catch (const std::out_of_range &oor) {
		  std::cout << "ERROR:: " << oor.what() << std::endl;
	       }
	    } else {
	       std::cout << "Can't find residue" << spec << std::endl;
	    }
	 } else {
	    try {

	       int read_number = 40;
	       std::vector<std::string> rtv = coot::util::non_standard_residue_types_in_molecule(mol);
	       for (unsigned int i=0; i<rtv.size(); i++)
		  if (rtv[i] != "HOH")
		     geom.try_dynamic_add(rtv[i], read_number++);

	       // spike-length probe-radius
	       bool use_waters = true;
	       coot::atom_overlaps_container_t overlaps(mol, &geom, use_waters,
							0.5, 0.25);
	       double dot_density = 0.2;
	       coot::atom_overlaps_dots_container_t c = overlaps.all_atom_contact_dots(dot_density);
	    }
	    catch (const std::out_of_range &oor) {
	       std::cout << "ERROR:: " << oor.what() << std::endl;
	    }
	 }
      } else {
	 std::cout << "ERROR:: Failed to read " << file_name << std::endl;
      }
   }

   return status;

}
