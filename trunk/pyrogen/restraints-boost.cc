
#include <GraphMol/GraphMol.h>

#include <boost/python.hpp>
using namespace boost::python;

#define HAVE_GSL
#include <ideal/simple-restraint.hh>
#include <coot-utils/coot-coord-utils.hh>
#include <lidia-core/rdkit-interface.hh>

#include "py-restraints.hh"


namespace coot {

   RDKit::ROMol *regularize(RDKit::ROMol &r);
   RDKit::ROMol *regularize_with_dict(RDKit::ROMol &r,
				      PyObject *py_restraints,
				      const std::string &comp_id);
   // This tries to get a residue from the dictionary using the
   // model_Cartn of the dict_atoms.  If that is not available, return
   // a molecule with atoms that do not have positions.
   RDKit::ROMol *rdkit_mol_chem_comp_pdbx(const std::string &chem_comp_dict_file_name,
					  const std::string &comp_id);
   // std::vector<RDKit::RWMol *> hydrogen_exchanges(const RDKit::ROMol &r);
   boost::python::list hydrogen_exchanges(const RDKit::ROMol &r);
}

BOOST_PYTHON_MODULE(libpyrogen_boost) {
   def("regularize",               coot::regularize,               return_value_policy<manage_new_object>());
   def("regularize_with_dict",     coot::regularize_with_dict,     return_value_policy<manage_new_object>());
   def("rdkit_mol_chem_comp_pdbx", coot::rdkit_mol_chem_comp_pdbx, return_value_policy<manage_new_object>());
   def("hydrogen_exchanges",       coot::hydrogen_exchanges);

   // to_python_converter<std::vector<std::string,class std::allocator<std::string> >, VecToList<std::string> >();
   // def("getListValue", getListValue);

   // class_<RDKit::RWMol, RDKit::RWMol*>("RDKit::RWMol");

}


RDKit::ROMol *
coot::regularize(RDKit::ROMol &mol_in) {

   RDKit::ROMol *m = new RDKit::ROMol(mol_in);
   return m;
}

RDKit::ROMol *
coot::regularize_with_dict(RDKit::ROMol &mol_in, PyObject *restraints_py, const std::string &comp_id) {

   coot::dictionary_residue_restraints_t dict_restraints = 
      monomer_restraints_from_python(restraints_py);
   RDKit::RWMol *m = new RDKit::RWMol(mol_in);
   CResidue *residue_p = make_residue(mol_in, 0, comp_id);
   if (! residue_p) {
      std::cout << "WARNING:: bad residue " << std::endl;
   } else {
      // deep copy residue and add to new molecule.
      CMMDBManager *cmmdbmanager = util::create_mmdbmanager_from_residue(residue_p);
      CResidue *new_residue_p = coot::util::get_first_residue(cmmdbmanager);
      PPCAtom residue_atoms = 0;
      int n_residue_atoms;
      new_residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      std::cout << "------------------ simple_refine() called from "
		<< "restraints-boost.cc:regularize_with_dict()" << std::endl;
      simple_refine(new_residue_p, cmmdbmanager, dict_restraints);
      std::cout << "------------------ simple_refine() finished" << std::endl;
      int iconf = 0;
      update_coords(m, iconf, new_residue_p);
   }
   return static_cast<RDKit::ROMol *>(m);
}


RDKit::ROMol *
coot::rdkit_mol_chem_comp_pdbx(const std::string &chem_comp_dict_file_name,
			       const std::string &comp_id) {

   RDKit::ROMol *mol = new RDKit::ROMol;
   coot::protein_geometry geom;
   geom.set_verbose(false);
   int read_number = 0;
   geom.init_refmac_mon_lib(chem_comp_dict_file_name, read_number);
   bool idealized = false;

   CResidue *r = geom.get_residue(comp_id, idealized);

   std::pair<bool, dictionary_residue_restraints_t> rest = geom.get_monomer_restraints(comp_id);
   if (rest.first) {
   
      if (r) {
	 // makes a 3d conformer

	 // 20140618: This was false - now try true so that it processes 110 (i.e. no more
	 //           valence error).
	 bool sanitize = true;
	 RDKit::RWMol mol_rw = coot::rdkit_mol(r, rest.second, "", sanitize);
	 RDKit::ROMol *m = new RDKit::ROMol(mol_rw);

	 // debug.  OK, so the bond orders are undelocalized here.
	 // debug_rdkit_molecule(&mol_rw);
      
	 return m;
      } else {
	 // makes a 2d conformer
	 RDKit::RWMol mol_rw = coot::rdkit_mol(rest.second);
	 RDKit::ROMol *m = new RDKit::ROMol(mol_rw);
	 bool canon_orient = false;
	 bool clear_confs = false;
	 int iconf = RDDepict::compute2DCoords(*m, NULL, canon_orient, clear_confs, 10, 20);
	 return m;
      }
   } 
   return mol;
} 



// std::vector<RDKit::RWMol *>
boost::python::list
coot::hydrogen_exchanges(const RDKit::ROMol &mol) {

   std::vector<RDKit::RWMol *> v;

   RDKit::ROMol *query_cooh = RDKit::SmartsToMol("[C^2](=O)O[H]");
   RDKit::ROMol *query_n = RDKit::SmartsToMol("[N^3;H2]");
   std::vector<RDKit::MatchVectType>  matches_cooh;
   std::vector<RDKit::MatchVectType>  matches_n;
   bool recursionPossible = true;
   bool useChirality = true;
   bool uniquify = true;
   int matched_cooh = RDKit::SubstructMatch(mol,*query_cooh,matches_cooh,uniquify,recursionPossible, useChirality);
   int matched_n    = RDKit::SubstructMatch(mol,*query_n,   matches_n,   uniquify,recursionPossible, useChirality);

   std::cout << "matches " << matches_cooh.size() << " " << matches_n.size() << std::endl;
   
   if (matches_cooh.size()) { 
      if (matches_n.size()) {
	 for (unsigned int imatch_cooh=0; imatch_cooh<matches_cooh.size(); imatch_cooh++) { 
	    for (unsigned int imatch_n=0; imatch_n<matches_n.size(); imatch_n++) {
	       std::cout << "matches_cooh: ";
	       for (unsigned int i=0; i<matches_cooh[imatch_cooh].size(); i++) { 
		  std::cout << " [" << matches_cooh[imatch_cooh][i].first 
			    << " "  << matches_cooh[imatch_cooh][i].second
			    << "]";
	       }
	       std::cout << std::endl;
	       std::cout << "matches_n: ";
	       for (unsigned int i=0; i<matches_n[imatch_n].size(); i++) { 
		  std::cout << " [" << matches_n[imatch_n][i].first 
			    << " "  << matches_n[imatch_n][i].second
			    << "]";
	       }
	       std::cout << std::endl;
	       
	       RDKit::RWMol *t_mol = new RDKit::RWMol(mol);
	       std::cout << "remove bond between atoms " << matches_cooh[imatch_cooh][2].second
			 << " and " << matches_cooh[imatch_cooh][3].second << std::endl;
	       std::cout << "add bond between " << matches_n[imatch_n][0].second
			 << " " << matches_cooh[imatch_cooh][3].second
			 << std::endl;

	       t_mol->removeBond(matches_cooh[imatch_cooh][2].second, matches_cooh[imatch_cooh][3].second);
	       t_mol->addBond(matches_n[imatch_n][0].second, matches_cooh[imatch_cooh][3].second);
	       RDKit::ATOM_SPTR at_o = (*t_mol)[matches_cooh[imatch_cooh][2].second];
	       RDKit::ATOM_SPTR at_n = (*t_mol)[matches_n[imatch_n][0].second];
	       std::cout << "at_o: atomic-number: " << at_o->getAtomicNum() << std::endl;
	       std::cout << "at_n: atomic-number: " << at_n->getAtomicNum() << std::endl;
	       at_o->setFormalCharge(at_o->getFormalCharge() - 1);
	       at_n->setFormalCharge(at_n->getFormalCharge() + 1);
	       v.push_back(t_mol);
	    }
	 }
      }
   }

   for (unsigned int i=0; i<v.size(); i++) { 
      std::cout << "   " << RDKit::MolToSmiles(*v[i]) << std::endl;
   }


   boost::python::list l;
   return l;

   // return v;
}
