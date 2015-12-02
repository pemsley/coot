
#include <clipper/ccp4/ccp4_map_io.h>

#include "../utils/coot-utils.hh"
#include "../coot-utils/residue-and-atom-specs.hh"
#include "../coot-utils/coot-coord-utils.hh"
#include "../coot-utils/coot-map-utils.hh"
#include "multi-peptide.hh"
#include "trace.hh"

int main(int argc, char **argv) {

   bool debug = false;

   std::string map_file_name = "trace-test.map";

   if (argc > 1) map_file_name = argv[1];

   if (coot::file_exists(map_file_name)) {

      try { 
	 clipper::CCP4MAPfile file;
	 clipper::Xmap<float> xmap;
	 file.open_read(map_file_name);
	 file.import_xmap(xmap);
	 file.close_read();

	 coot::trace t(xmap);

	 if (argc > 3) { 
	    std::string pdb_name = argv[2];
	    mmdb::Manager *mol = new mmdb::Manager;
	    std::cout << "Reading coordinate file: " << pdb_name.c_str() << "\n";
	    mmdb::ERROR_CODE err = mol->ReadCoorFile(pdb_name.c_str());
	    if (err) {
	       std::cout << "There was an error reading " << pdb_name.c_str() << ". \n";
	       std::cout << "ERROR " << err << " READ: "
			 << mmdb::GetErrorDescription(err) << std::endl;
	    } else {

	       coot::residue_spec_t spec("A", 425, "");
	       mmdb::Residue *r = coot::util::get_residue(spec, mol);

	       if (r) {

		  coot::protein_geometry geom;
		  geom.init_standard();
		  geom.remove_planar_peptide_restraint();
		  std::pair<float, float> mv = coot::util::mean_and_variance(xmap);
		  coot::minimol::fragment f =
		     coot::multi_build_N_terminal_ALA(r, "A", 20, 5000, geom, xmap, mv);
	       }
	    }
	 } else {

	    // not testing rama multi-build


	    // test from a pdb file
	    std::string test_pdb_file_name = "test-trace-template.pdb";

	    if (argc > 2) {
	       std::string fn = argv[2];
	       if (coot::file_exists(fn)) 
		  test_pdb_file_name = fn;
	    }

	 
	    // test with a null moll or flood mol
	    if (coot::file_exists(test_pdb_file_name)) {

	       mmdb::Manager *mol = new mmdb::Manager;
	       mmdb::ERROR_CODE err = mol->ReadCoorFile(test_pdb_file_name.c_str());
	       if (! err) {

		  if (coot::file_exists("test-scales")) {
		     t.optimize_weights(mol);

		     std::cout << "------------- Done optimize_weights() " << std::endl;

		  } else { 
		     std::cout << "----------------------------------------------------\n";
		     std::cout << "----------------------------------------------------\n";
		     std::cout << "running with test mol: " << std::endl;
		     std::cout << "----------------------------------------------------\n";
		     std::cout << "----------------------------------------------------\n";
		     t.test_model(mol);
		  }
	       }
	    } else {

	       t.action();

	    }
	 }
      }
      
      // problem reading the map, perhaps?
      // 
      catch (const clipper::Message_fatal &mess) {
	 std::cout << "ERROR:: " << mess.text() << std::endl;
      }
   }
   return 0;
}
