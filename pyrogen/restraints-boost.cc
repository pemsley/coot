
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
   RDKit::ROMol *rdkit_mol_chem_comp_pdbx(const std::string &chem_comp_dict_file_name,
					  const std::string &comp_id,
					  const std::string &experimental_model_pdb_file_name);
}

BOOST_PYTHON_MODULE(restraints_boost) {
   def("regularize",  coot::regularize, return_value_policy<manage_new_object>());
   def("regularize_with_dict",  coot::regularize_with_dict, return_value_policy<manage_new_object>());
   def("rdkit_mol_chem_comp_pdbx", coot::rdkit_mol_chem_comp_pdbx, return_value_policy<manage_new_object>());
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
      CMMDBManager *cmmdbmanager = util::create_mmdbmanager_from_residue(NULL, residue_p);
      CResidue *new_residue_p = coot::util::get_first_residue(cmmdbmanager);
      PPCAtom residue_atoms = 0;
      int n_residue_atoms;
      new_residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      simple_refine(new_residue_p, cmmdbmanager, dict_restraints);
      int iconf = 0;
      update_coords(m, iconf, new_residue_p);
   }
   return static_cast<RDKit::ROMol *>(m);
}


RDKit::ROMol *
coot::rdkit_mol_chem_comp_pdbx(const std::string &chem_comp_dict_file_name,
			       const std::string &comp_id,
			       const std::string &experimental_model_pdb_file_name) {

   RDKit::ROMol *mol = new RDKit::ROMol;
   coot::protein_geometry geom;
   int read_number = 0;
   geom.init_refmac_mon_lib(chem_comp_dict_file_name, read_number);
   bool idealized = false;

   CResidue *r = geom.get_residue(comp_id, idealized);
   if (r) {
      RDKit::RWMol mol_rw = coot::rdkit_mol_sanitized(r, geom);
      RDKit::ROMol *m = new RDKit::ROMol(mol_rw);
      return m;
   }
   return mol;

} 
