
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
}

BOOST_PYTHON_MODULE(restraints_boost) {
   def("regularize",               coot::regularize,               return_value_policy<manage_new_object>());
   def("regularize_with_dict",     coot::regularize_with_dict,     return_value_policy<manage_new_object>());
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
			       const std::string &comp_id) {

   RDKit::ROMol *mol = new RDKit::ROMol;
   coot::protein_geometry geom;
   geom.set_verbose(false);
   int read_number = 0;
   geom.init_refmac_mon_lib(chem_comp_dict_file_name, read_number);
   bool idealized = false;

   CResidue *r = geom.get_residue(comp_id, idealized);

   if (r) {
      // makes a 3d conformer
      RDKit::RWMol mol_rw = coot::rdkit_mol_sanitized(r, geom);
      RDKit::ROMol *m = new RDKit::ROMol(mol_rw);

      // debug.  OK, so the bond orders are undelocalized here.
      debug_rdkit_molecule(&mol_rw);
      
      return m;
   } else {
      // makes a 2d conformer
      std::pair<bool, dictionary_residue_restraints_t> rest = geom.get_monomer_restraints(comp_id);
      if (rest.first) { 
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
