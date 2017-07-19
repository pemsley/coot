
//  None of these functions should depend on functions in src.
//
//  Consider putting these functions in a directory other than src.
//  pli perhaps?  Rehame pli to be more encompassing?
//

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#include <boost/python.hpp>

#include "utils/coot-utils.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "lidia-core/rdkit-interface.hh"

namespace coot {
   int import_rdkit_molecule(const RDKit::ROMol &m, int conf_id, const std::string &new_comp_id);
   boost::python::list extract_ligands_from_coords_file(const std::string &file_name);
   boost::python::object get_ligand_interactions(const std::string &file_name,
						 PyObject *ligand_spec);
}

//                      Functions that are importable from python

BOOST_PYTHON_MODULE(coot_extended) {

   def("extract_ligands_from_coords_file", coot::extract_ligands_from_coords_file);
   def("get_ligand_interactions", coot::get_ligand_interactions);

}

boost::python::list
coot::extract_ligands_from_coords_file(const std::string &file_name) {

   boost::python::list rdkit_mols_list;
   protein_geometry geom;

   if (coot::file_exists(file_name)) {
      mmdb::Manager *mol = new mmdb::Manager;
      mol->ReadCoorFile(file_name.c_str());
      std::vector<mmdb::Residue *> v = util::get_hetgroups(mol); // no waters

      std::cout << "Found " << v.size() << " hetgroups " << std::endl;
      if (v.size() > 0) {
	 int read_number = 0;
	 for (std::size_t i=0; i<v.size(); i++) {
	    std::string res_name = v[i]->GetResName();
	    int imol = 0;
	    if (geom.have_dictionary_for_residue_type(res_name, imol, read_number++)) { // autoloads
	       std::pair<bool, coot::dictionary_residue_restraints_t> rp =
		  geom.get_monomer_restraints(res_name, imol);
	       if (rp.first) {
		  try {
		     RDKit::RWMol rdkm = rdkit_mol(v[i], rp.second);
		     RDKit::ROMol *cm_p = new RDKit::ROMol(rdkm);
		     boost::shared_ptr<RDKit::ROMol> xx(cm_p);
		     // maybe I can append(xx) rather than needing this step:
		     boost::python::object obj(xx);
		     rdkit_mols_list.append(obj);
		  }
		  catch (const std::runtime_error &rte) {
		     std::cout << "WARNING:: " << rte.what() << std::endl;
		  }
		  catch (const std::exception &e) {
		     std::cout << "WARNING:: " << e.what() << std::endl;
		  }
	       }
	    }
	 }
      }
   }
   return rdkit_mols_list;
}

#include "pli/protein-ligand-interactions.hh"
#include "pli/pi-stacking.hh"
#include "coot-utils/reduce.hh"

// Pass also filename for the CCD for the ligand (or its neighbours)
// or (probably better) a directory.
//
boost::python::object
coot::get_ligand_interactions(const std::string &file_name,
			      PyObject *ligand_spec_py) {

   // this is how to convert (any) PyObject * to a boost::python::object
   // (consider if ref-counting is an issue).
   boost::python::object o(boost::python::handle<>(Py_False));
   float h_bond_dist_max = 3.6;
   float residues_near_radius = 5.0; // pass this

   if (PyList_Check(ligand_spec_py)) {
      if (PyObject_Length(ligand_spec_py) == 3) {
	 PyObject  *chain_id_py = PyList_GetItem(ligand_spec_py, 0);
	 PyObject     *resno_py = PyList_GetItem(ligand_spec_py, 1);
	 PyObject  *ins_code_py = PyList_GetItem(ligand_spec_py, 2);
	 if (PyInt_Check(resno_py)) {
	    int res_no = PyInt_AsLong(resno_py);
	    std::string chain_id = PyString_AsString(chain_id_py);
	    std::string ins_code  = PyString_AsString(ins_code_py);
	    residue_spec_t rs(chain_id, res_no, ins_code);
	    mmdb::Manager *mol = new mmdb::Manager;
	    mol->ReadCoorFile(file_name.c_str());
	    mmdb::Residue *residue_p = util::get_residue(rs, mol);
	    if (residue_p) {
	       // read in a dictionary
	       int imol = 0; // dummy
	       int read_number = 0;
	       protein_geometry geom;
	       geom.init_standard();
	       std::string rn = residue_p->GetResName();
	       geom.try_dynamic_add("MG", read_number++);
	       if (geom.have_dictionary_for_residue_type(rn, imol, read_number++)) { // autoloads
		  std::pair<bool, coot::dictionary_residue_restraints_t> rp =
		     geom.get_monomer_restraints(rn, imol);
		  if (rp.first) {
		     reduce r(mol, imol);
		     r.add_geometry(&geom);
		     r.add_hydrogen_atoms();

		     // protein-ligand interactions (various sorts of bonds) are calculated without
		     // needing rdkit.

		     // consider where a peptide is the ligand
		     std::vector<fle_ligand_bond_t> v =
			protein_ligand_interactions(residue_p, mol, &geom, h_bond_dist_max);

		     std::cout << "INFO:: found " << v.size() << " bonds/interactions " << std::endl;
		     for (std::size_t i=0; i<v.size(); i++)
			std::cout << "INFO::    " << v[i] << std::endl;

		     // now replace o:
		     PyObject *new_o_py = PyList_New(2);
		     PyList_SetItem(new_o_py, 1, Py_False); // replaced later maybe
		     PyObject *list_py = PyList_New(v.size());
		     for (std::size_t i=0; i<v.size(); i++) {
			PyObject *item_py = PyList_New(2);
			// transfer ligand-atom-spec, interacting-atom-spec, bond-type, bond-length
			PyList_SetItem(item_py, 0, PyInt_FromLong(v[i].bond_type));
			PyList_SetItem(item_py, 1, PyFloat_FromDouble(v[i].bond_length));
			PyList_SetItem(list_py, i, item_py);
		     }
		     PyList_SetItem(new_o_py, 0, list_py);

		     std::vector<mmdb::Residue *> neighb_residues =
			coot::residues_near_residue(residue_p, mol, residues_near_radius);

		     // pi-stacking interactions, using the CCD (or current refmac monomer
		     // library, I think) needs an rdkit molecule
		     //
		     try {
			RDKit::ROMol rdkm = rdkit_mol(residue_p, rp.second);
			pi_stacking_container_t pi_stack_info(rp.second, neighb_residues, residue_p, rdkm);
			std::cout << "INFO:: found " << pi_stack_info.stackings.size()
				  << " pi-stacking interactions" << std::endl;
			if (pi_stack_info.stackings.size() > 0) {
			   PyObject *pi_stack_info_py = PyList_New(pi_stack_info.stackings.size());
			   for (std::size_t j=0; j<pi_stack_info.stackings.size(); j++) {
			      PyObject *pi_stack_instance_py = PyList_New(1); // transfer more info
			      const pi_stacking_instance_t &stack_instance = pi_stack_info.stackings[j];
			      PyObject *pi_stack_py = PyList_New(1); // transfer more info
			      std::vector<std::string> lran = stack_instance.ligand_ring_atom_names;
			      PyObject *atom_name_list_py = PyList_New(lran.size());
			      for (std::size_t k=0; k<lran.size(); k++)
				 PyList_SetItem(atom_name_list_py, k, PyString_FromString(lran[k].c_str()));
			      PyList_SetItem(pi_stack_instance_py, 0, atom_name_list_py);
			      PyList_SetItem(pi_stack_info_py, j, pi_stack_instance_py);
			   }
			   PyList_SetItem(new_o_py, 1, pi_stack_info_py);
			}
		     }
		     catch (const std::runtime_error &rte) {
		     }

		     o = boost::python::object(boost::python::handle<>(new_o_py));

		  }
	       }
	    }
	 }
      }
   }
   return o;
}

#endif // MAKE_ENHANCED_LIGAND_TOOLS
