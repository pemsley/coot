
#include "glyco-torsions.hh"

int main(int argc, char **argv) {

   if (argc <= 2) {
      std::cout << "Usage: " << argv[0] << " template-pdb link-type new-comp-id (generate mode)\n";
      std::cout << "   or: " << argv[0] << " test link-type new-comp-id \n";
   } else { 

      if (std::string(argv[1]) == "test") {

	 // use the link-by-torsion reference files to make a new residue (i.e. to test
	 // the function)

	 if (argc > 3) {
	    std::string link_type = argv[2];
	    std::string new_residue_type = argv[3];
	    
	    std::string file_name = "link-by-torsion-to-pyranose-core-" + link_type + ".tab";
	    if (! coot::file_exists(file_name)) {
	       std::cout << "ERROR: file not found " << file_name << std::endl;
	    } else { 
	       coot::link_by_torsion_t l(file_name);
	       l.new_residue_type = new_residue_type;

	       // Get a base residue
	       CMMDBManager *mol = new CMMDBManager;
	       std::string pdb_file_name =
		  "pdb-templates/pyranose-pyranose-via-" + link_type + ".pdb";
	       
	       mol->ReadPDBASCII(pdb_file_name.c_str());
	       CResidue *base_residue_p = coot::util::get_first_residue(mol);
	       CResidue *r = l.make_residue(base_residue_p);
	       if (r) {
		  CMMDBManager *mol = coot::util::create_mmdbmanager_from_residue(r);
		  mol->WritePDBASCII("new-residue.pdb");
	       }
	    }
	 }

      } else {

	 // make the link-by-torsion reference files

	 std::string file_name = argv[1];
	 std::string link_type = argv[2]; // e.g. "ALPHA1-6";
	 std::string new_residue_type = "MAN";
      
	 CMMDBManager *mol = new CMMDBManager;
	 int status = mol->ReadPDBASCII(file_name.c_str());
	 if (status != Error_NoError) {
	    std::cout << "ERROR:: on reading " << file_name << std::endl;
	 } else {
	    std::pair<CResidue *, CResidue *> p = coot::link_by_torsion_t::get_residue_pair(mol);
	    if (p.first && p.second) {

	       // generate link torsions, write link torsions
	       coot::link_by_torsion_base_t to_core = coot::get_names_for_link_type(link_type);
	       if (! to_core.filled()) {
		  std::cout << "ERROR:: " << link_type << std::endl;
	       } else {
		  coot::link_by_torsion_t l_to_core(to_core, p.first, p.second);
		  if (l_to_core.filled()) {
		     std::string file_name = "link-by-torsion-to-pyranose-core-" + link_type + ".tab";
		     l_to_core.write(file_name);
		  }

		  // and the decorations on the core:
		  coot::link_by_torsion_base_t decor = coot::mannose_decorations();
		  coot::link_by_torsion_t l_decor(decor, p.first, p.second);
		  if (l_decor.filled()) {
		     std::string file_name =  new_residue_type + "decorations.tab";
		     l_decor.write(file_name);
		  }
	       }
	    }
	 }
      }
   }
   return 0;
}
   
