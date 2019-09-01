
#include <string>

#include "coords/mmdb-extras.h"
#include "coords/mmdb.h"

struct atom_vector_map_t {
   coot::atom_quad q1;
   coot::atom_quad q2;
   int res_no;
   std::string chain_id;
   std::string fn;
};

void print_tors(std::map<std::string, std::vector<atom_vector_map_t> > &atom_vector_map) {

   std::map<std::string, std::vector<atom_vector_map_t> >::const_iterator it;
   for(it=atom_vector_map.begin(); it!=atom_vector_map.end(); it++) {
      const std::string &key = it->first;
      const std::vector<atom_vector_map_t> &vv = it->second;
      for (std::size_t i=0; i<vv.size(); i++) {
	 const atom_vector_map_t &qp = vv[i];
	 double tors_1 = qp.q1.torsion();
	 double tors_2 = qp.q2.torsion();
	 std::cout << "pseudo-torsion " << " " << qp.fn << " " << qp.chain_id << " " << qp.res_no << " "
		   << key << "  N " << tors_1 << " C " << tors_2 << "\n";
      }
   }
}

#include "utils/coot-utils.hh"

int main(int argc, char **argv) {

   if (argc > 1) {

      std::string file_name(argv[1]);
      atom_selection_container_t asc = get_atom_selection(file_name, false, 0);
      if (asc.read_success) {

	 int imod = 1;
	 mmdb::Model *model_p = asc.mol->GetModel(imod);
	 if (model_p) {

	    // the lists of 4 atoms for each residue type
	    std::map<std::string, std::vector<atom_vector_map_t> > atom_vector_map;

	    int n_chains = model_p->GetNumberOfChains();
	    for (int ichain=0; ichain<n_chains; ichain++) {
	       mmdb::Chain *chain_p = model_p->GetChain(ichain);
	       int nres = chain_p->GetNumberOfResidues();
	       for (int ires=0; ires<nres; ires++) {
		  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
		  int n_atoms = residue_p->GetNumberOfAtoms();
		  std::vector<mmdb::Atom *> reference_atoms(4, nullptr);
		  for (int iat=0; iat<n_atoms; iat++) {
		     mmdb::Atom *at = residue_p->GetAtom(iat);
		     std::string atom_name(at->GetAtomName());
		     std::string altconf(at->altLoc);
		     if (altconf.empty()) {
			if (atom_name == " N  ") reference_atoms[0] = at;
			if (atom_name == " CA ") reference_atoms[1] = at;
			if (atom_name == " C  ") reference_atoms[2] = at;
			if (atom_name == " CB ") reference_atoms[3] = at;
		     }
		  }
		  if (reference_atoms[0] && reference_atoms[1] && reference_atoms[2] && reference_atoms[3]) {
		     std::string res_name(residue_p->GetResName());
		     // yes
		     coot::atom_quad q1(reference_atoms[0], reference_atoms[2],
					reference_atoms[1], reference_atoms[3]);
		     coot::atom_quad q2(reference_atoms[2], reference_atoms[0],
					reference_atoms[1], reference_atoms[3]);
		     atom_vector_map_t avm;
		     avm.q1 = q1;
		     avm.q2 = q2;
		     avm.fn = coot::util::file_name_non_directory(file_name);
		     avm.res_no = residue_p->GetSeqNum();
		     avm.chain_id = residue_p->GetChainID();
		     if (! avm.chain_id.empty())
			atom_vector_map[res_name].push_back(avm);
		  }
	       }
	    }

	    print_tors(atom_vector_map);

	 }
      }
   }
   return 0;
} 
