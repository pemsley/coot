
#include <mmdb2/mmdb_manager.h>
#include "coords/mmdb-extras.h"
#include "coords/mmdb.h"

#include "crankshaft.hh"

int main(int argc, char **argv) {

   int status = 0;

   if (argc > 1) {
      std::string pdb_file_name(argv[1]);
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, 1, 0);
      if (asc.read_success) {
	 coot::crankshaft cs(asc.mol);
	 // cs.test();
	 zo::rama_table_set zorts;
	 // coot::residue_spec_t rs("A", 51); // goodie
	 // coot::residue_spec_t rs("A", 40); // baddie
	 coot::residue_spec_t rs("A", 20); // baddie
	 std::vector<std::pair<float, float> > r = cs.spin_search(rs, zorts);
	 unsigned int n_samples = 30;
	 if (argc > 2) {
	    try {
	       // what's the std:: version of this?
	       n_samples = coot::util::string_to_int(argv[2]);
	    }
	    catch (const std::runtime_error &rte) {
	       std::cout << "WARNING:: " << rte.what() << std::endl;
	    }
	 }

	 if (false) {
	    if (r.size()) {
	       std::cout << "Residue " << rs << std::endl;
	       for (std::size_t i=0; i<r.size(); i++) {
		  std::cout << i << "   " << r[i].first << " " << r[i].second << std::endl;
	       }
	    }
	 }
	 if (true) {
	    // the moves the atoms
	    // cs.triple_spin_search(rs, zorts, true, n_samples);
	    cs.find_maxima(rs, zorts, n_samples);
	    // where did they move to?
	    asc.mol->WritePDBASCII("moved.pdb");
	 }
      } else {
	 std::cout << "read failure " << pdb_file_name << std::endl;
      }
   }

   return status;
};
