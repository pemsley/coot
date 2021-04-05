
//  None of these functions should depend on functions in src.
//
//  Consider putting these functions in a directory other than src.
//  pli perhaps?  Rehame pli to be more encompassing?
//

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#include <boost/python.hpp>

#include "utils/coot-utils.hh"
#include "geometry/mol-utils.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "lidia-core/rdkit-interface.hh"
#include "pli/specs.hh"
#include "pli/protein-ligand-interactions.hh"
#include "pli/pi-stacking.hh"

#include "coot-utils/reduce.hh"
#include "coot-utils/atom-overlaps.hh"



namespace coot {
   boost::python::list extract_ligands_from_coords_file(const std::string &file_name,
							   const std::string &ccd_dir_or_cif_file_name);
   boost::python::object get_ligand_interactions(const std::string &file_name,
						 PyObject *ligand_spec,
						 const std::string &dirname_or_cif_file_name);

   boost::python::object contact_dots_from_coordinates_file(const std::string &file_name,
							    bool add_hydrogens_flag);

   boost::python::object atom_overlaps_from_coordinates_file(const std::string &file_name,
							     bool add_hydrogens_flag);

   // helper function
   void get_ligand_interactions_read_cifs(const std::string &residue_name,
					  const std::string &ccd_dir_or_cif_file_name,
					  protein_geometry *geom_p);

   // helper function
   void extract_ligands_from_coords_file_try_read_cif(const std::string &ccd_dir_or_cif_file_name,
						      const std::string &res_name,
						      protein_geometry *geom_p,
						      int *cif_read_number_p);
}

//                      Functions that are importable from python

BOOST_PYTHON_MODULE(coot_extended) {

   def("extract_ligands_from_coords_file", coot::extract_ligands_from_coords_file);
   def("get_ligand_interactions", coot::get_ligand_interactions);
   def("contact_dots_from_coordinates_file", coot::contact_dots_from_coordinates_file);
   def("atom_overlaps_from_coordinates_file", coot::atom_overlaps_from_coordinates_file);

}


boost::python::list
coot::extract_ligands_from_coords_file(const std::string &file_name,
				       const std::string &ccd_dir_or_cif_file_name) {

   int rn = 42; // pass this and pass it back.

   boost::python::list rdkit_mols_list;
   bool debug = true;
   protein_geometry geom;
   geom.set_verbose(false);

   if (coot::file_exists(file_name)) {
      mmdb::Manager *mol = new mmdb::Manager;
      mol->ReadCoorFile(file_name.c_str());
      std::vector<mmdb::Residue *> v = util::get_hetgroups(mol); // no waters

      if (debug) {
	 std::cout << "Found " << v.size() << " hetgroups " << std::endl;
	 for (std::size_t i=0; i<v.size(); i++)
	    std::cout << " " << i << " " << residue_spec_t(v[i]) << " "
		      << v[i]->GetResName() << std::endl;
      }
      if (v.size() > 0) {
	 int read_number = 0;
	 for (std::size_t i=0; i<v.size(); i++) {
	    mmdb::Residue *residue_p = v[i];
	    std::string res_name = residue_p->GetResName();
	    int imol = 0;
	    extract_ligands_from_coords_file_try_read_cif(ccd_dir_or_cif_file_name, res_name, &geom, &rn);
	    rn++; // increment the cif read number (to prevent double-addition of bonds etc)
	    if (! geom.have_dictionary_for_residue_type(res_name, imol, read_number++)) { // autoloads
	       std::string message = "Missing dictionary for type " + res_name;
	       rdkit_mols_list.append(message);
	    } else {
	       std::pair<bool, dictionary_residue_restraints_t> rp =
		  geom.get_monomer_restraints(res_name, imol);
	       if (! rp.first) {
		  std::string message = "Missing dictionary for type " + res_name;
		  rdkit_mols_list.append(message);
	       } else {
		  try {
		     // rdkit_mol() will return a molecule with no atoms if
		     // there is only atoms with a non-blank alt conf
		     //
		     // so loop over the residue alt confs and make new molecules.
		     // I suppose the another/better way is to add conformers for
		     // every alt conf - not pass the alt conf as an arg
		     // Maybe later.
		     //
		     std::vector<std::string> ac = util::get_residue_alt_confs(residue_p);
		     for (std::size_t iac=0; iac<ac.size(); iac++) {
			RDKit::RWMol rdkm = rdkit_mol(residue_p, rp.second, ac[iac]);
			RDKit::ROMol *cm_p = new RDKit::ROMol(rdkm);
			boost::shared_ptr<RDKit::ROMol> xx(cm_p);
			// maybe I can append(xx) rather than needing this step:
			boost::python::object obj(xx);
			rdkit_mols_list.append(boost::python::object(obj));
			if (false)
			   std::cout << "added to list from residue "
				     << residue_spec_t(v[i]) << " " << res_name << std::endl;
		     }
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

// Pass also filename for the CCD for the ligand (or its neighbours)
// or (probably better) a directory.
//
// This needs to do a better job handling ligand residues with alt confs (5pb7)
//
boost::python::object
coot::get_ligand_interactions(const std::string &file_name,
			      PyObject *ligand_spec_py,
			      const std::string &dirname_or_cif_file_name) {

   // debug
   py_residue_spec_t debug = (ligand_spec_py);

   // this is how to convert (any) PyObject * to a boost::python::object
   // (consider if ref-counting is an issue).
   boost::python::object o(boost::python::handle<>(Py_False));
   float h_bond_dist_max = 3.6;
   float residues_near_radius = 5.0; // pass this
   std::string ligand_spec_py_str = PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyObject_Str(ligand_spec_py)));

   if (! PyList_Check(ligand_spec_py)) {
      // tell me what it was
      std::cout << "WARNING:: ligand spec was not a spec " << ligand_spec_py_str << std::endl;
   } else {
      Py_ssize_t sl = PyObject_Length(ligand_spec_py);
      if (sl != 3) {
	 std::cout << "WARNING:: ligand spec was not a triple list " << ligand_spec_py_str << std::endl;
      } else {
	 PyObject  *chain_id_py = PyList_GetItem(ligand_spec_py, 0);
	 PyObject     *resno_py = PyList_GetItem(ligand_spec_py, 1);
	 PyObject  *ins_code_py = PyList_GetItem(ligand_spec_py, 2);
	 if (! PyLong_Check(resno_py)) {
	    std::cout << "WARNING:: resno in spec is not Int in " << ligand_spec_py_str << std::endl;
	 } else {
	    int res_no = PyLong_AsLong(resno_py);
	    std::string chain_id = PyBytes_AS_STRING(PyUnicode_AsUTF8String(chain_id_py));
	    std::string ins_code = PyBytes_AS_STRING(PyUnicode_AsUTF8String(ins_code_py));
	    residue_spec_t rs(chain_id, res_no, ins_code);
	    mmdb::Manager *mol = new mmdb::Manager;
	    mol->ReadCoorFile(file_name.c_str());
	    mmdb::Residue *residue_p = util::get_residue(rs, mol);
	    if (! residue_p) {
	       std::cout << "WARNING:: ligand spec was not found " << residue_spec_t(residue_p) << std::endl;
	    } else {
	       // read in a dictionary
	       int imol = 0; // dummy
	       int read_number = 0;
	       protein_geometry geom;
	       geom.set_verbose(false);
	       geom.init_standard();
	       std::string rn = residue_p->GetResName();
	       get_ligand_interactions_read_cifs(rn, dirname_or_cif_file_name, &geom);
	       if (geom.have_dictionary_for_residue_type(rn, imol, read_number++)) { // autoloads
		  std::pair<bool, coot::dictionary_residue_restraints_t> rp =
		     geom.get_monomer_restraints(rn, imol);

		  // we want to add energy types to the dictionary (if it was a dictionary from the CCD).
		  //
		  if (rp.first) {
		     try {
			// rdkit_mol() will return a molecule with no atoms if
			// there is only atoms with a non-blank alt conf
			//
			// so loop over the residue alt confs and make new molecules.
			// I suppose the another/better way is to add conformers for
			// every alt conf - not pass the alt conf as an arg
			// Maybe later.
			//
			std::vector<std::string> ac = util::get_residue_alt_confs(residue_p);
			for (std::size_t iac=0; iac<ac.size(); iac++) {
			   RDKit::ROMol rdkm = rdkit_mol(residue_p, rp.second, ac[iac]);
			   // set_dictionary_atom_types(&rp.second);
			   set_dictionary_atom_types_from_mol(&rp.second, &rdkm);
			   geom.replace_monomer_restraints(rn, protein_geometry::IMOL_ENC_ANY, rp.second);

			   reduce r(mol, imol);
			   r.set_verbose_output(false); // don't tell me about missing atoms
			   r.add_geometry(&geom);
			   r.add_hydrogen_atoms();

			   // protein-ligand interactions (various sorts of bonds) are calculated without
			   // needing rdkit.

			   // consider where a peptide is the ligand
			   std::vector<fle_ligand_bond_t> v =
			      protein_ligand_interactions(residue_p, mol, &geom, h_bond_dist_max);

			   std::cout << "INFO:: found " << v.size() << " bonds/interactions " << std::endl;
			   for (std::size_t i=0; i<v.size(); i++)
			      std::cout << "INFO::  " << i << "  " << v[i] << std::endl;

			   // now replace o:
			   PyObject *new_o_py = PyList_New(2);
			   PyList_SetItem(new_o_py, 1, Py_False); // replaced later maybe
			   PyObject *list_py = PyList_New(v.size());
			   for (std::size_t i=0; i<v.size(); i++) {
			      PyObject *item_py = PyList_New(5);
			      // transfer ligand-atom-spec, interacting-atom-spec, bond-type, bond-length
			      py_atom_spec_t pas_1(v[i].ligand_atom_spec);
			      py_atom_spec_t pas_2(v[i].interacting_residue_atom_spec);
			      PyObject *is_water_py = Py_False;
			      if (v[i].is_H_bond_to_water)
				 is_water_py = Py_True;
			      PyList_SetItem(item_py, 0, PyLong_FromLong(v[i].bond_type));
			      PyList_SetItem(item_py, 1, PyFloat_FromDouble(v[i].bond_length));
			      PyList_SetItem(item_py, 2, pas_1.pyobject());
			      PyList_SetItem(item_py, 3, pas_2.pyobject());
			      PyList_SetItem(item_py, 4, is_water_py);
			      PyList_SetItem(list_py, i, item_py);
			   }
			   PyList_SetItem(new_o_py, 0, list_py);

			   // stackings
			
			   std::vector<mmdb::Residue *> neighb_residues =
			      coot::residues_near_residue(residue_p, mol, residues_near_radius);

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
				    PyList_SetItem(atom_name_list_py, k, PyUnicode_FromString(lran[k].c_str()));
				 PyList_SetItem(pi_stack_instance_py, 0, atom_name_list_py);
				 PyList_SetItem(pi_stack_info_py, j, pi_stack_instance_py);
			      }
			      PyList_SetItem(new_o_py, 1, pi_stack_info_py);
			   }

			   o = boost::python::object(boost::python::handle<>(new_o_py));

			}
		     }
		     catch (const std::runtime_error &rte) {
			std::cout << "WARNING:: " << rte.what() << std::endl;
		     }
		  }
	       }
	    }
	 }
      }
   }
   return o;
}

void
coot::get_ligand_interactions_read_cifs(const std::string &residue_name,
					const std::string &ccd_dir_or_cif_file_name,
					protein_geometry *geom_p) {

   if (file_exists(ccd_dir_or_cif_file_name)) {
      std::string cif_file_name = ccd_dir_or_cif_file_name;
      if (is_directory_p(ccd_dir_or_cif_file_name)) {
	 std::string dir_name = ccd_dir_or_cif_file_name;
	 //
	 // all files in this directory of just the one...?
	 // comma-separated list?
	 //
	 cif_file_name = util::append_dir_file(dir_name, residue_name + ".cif");
      }
      if (file_exists(cif_file_name)) {
	 int rn = 42;
	 geom_p->init_refmac_mon_lib(cif_file_name, rn++);
      }
   }
}

// helper function
void
coot::extract_ligands_from_coords_file_try_read_cif(const std::string &ccd_dir_or_cif_file_name,
						    const std::string &residue_name,
						    protein_geometry *geom_p,
						    int *cif_read_number_p) {

   // test that geom has restraints for residue_name first
   if (geom_p->have_dictionary_for_residue_type_no_dynamic_add(residue_name)) {
      // nothing
   } else {
      if (file_exists(ccd_dir_or_cif_file_name)) {
	 std::string cif_file_name = ccd_dir_or_cif_file_name;
	 if (is_directory_p(ccd_dir_or_cif_file_name)) {
	    std::string dir_name = ccd_dir_or_cif_file_name;

	    std::string cif_file_name = util::append_dir_file(dir_name, residue_name + ".cif");
	    if (file_exists(cif_file_name)) {
	       geom_p->init_refmac_mon_lib(cif_file_name, *cif_read_number_p);
	       *cif_read_number_p++;
	    }
	 } else {
	    geom_p->init_refmac_mon_lib(ccd_dir_or_cif_file_name, *cif_read_number_p);
	    *cif_read_number_p++;
	 }
      }
   }
}

// Input a pdb file that has Hydrogen atoms.
//
boost::python::object
coot::contact_dots_from_coordinates_file(const std::string &file_name, bool add_hydrogens_flag) {

   double dot_density = 0.15; // pass this?
   boost::python::object o(boost::python::handle<>(Py_False)); // return this

   protein_geometry geom;
   geom.set_verbose(false);
   geom.init_standard();

   if (coot::file_exists(file_name)) {
      mmdb::Manager *mol = new mmdb::Manager;
      mmdb::ERROR_CODE err = mol->ReadCoorFile(file_name.c_str());
      if (err) {
	 std::cout << "Error reading coordinates file " << file_name << std::endl;
      } else {
	 // Happy Path
	 try {

	    int read_number = 40;
	    std::vector<std::string> rtv = coot::util::non_standard_residue_types_in_molecule(mol);
	    for (unsigned int i=0; i<rtv.size(); i++)
	       if (rtv[i] != "HOH")
		  geom.try_dynamic_add(rtv[i], read_number++);

	    // spike-length probe-radius
	    bool ignore_waters = false;
	    coot::atom_overlaps_container_t overlaps(mol, &geom, ignore_waters, 0.5, 0.25);
	    coot::atom_overlaps_dots_container_t c = overlaps.all_atom_contact_dots(dot_density);
	    const std::map<std::string, std::vector<coot::atom_overlaps_dots_container_t::dot_t> > &dots = c.dots;
	    const coot::atom_overlaps_dots_container_t::spikes_t &clashes = c.clashes;
	    PyObject *atom_overlaps_py = PyList_New(2);
	    PyObject *dots_map_py = PyDict_New();
	    PyObject *clashes_py  = PyList_New(clashes.size());

	    std::map<std::string, std::vector<coot::atom_overlaps_dots_container_t::dot_t> >::const_iterator it;
	    for (it = dots.begin(); it!=dots.end(); it++) {
	       std::string col_key = it->first;
	       // std::cout << col_key << std::endl;
	       const std::vector<coot::atom_overlaps_dots_container_t::dot_t> &v = it->second;
	       PyObject *dots_list_py = PyList_New(v.size());
	       for (std::size_t ii=0; ii<v.size(); ii++) {
		  const coot::atom_overlaps_dots_container_t::dot_t &dot = v[ii];
		  // std::cout << "   " << dot.overlap << std::endl;
		  PyObject *dot_py = PyList_New(3);
		  PyObject *coords_py = PyList_New(3);
		  PyList_SetItem(coords_py, 0, PyFloat_FromDouble(dot.pos.x()));
		  PyList_SetItem(coords_py, 1, PyFloat_FromDouble(dot.pos.y()));
		  PyList_SetItem(coords_py, 2, PyFloat_FromDouble(dot.pos.z()));
		  PyList_SetItem(dot_py, 0, PyFloat_FromDouble(dot.overlap));
		  PyList_SetItem(dot_py, 1, coords_py);
		  PyList_SetItem(dot_py, 2, PyUnicode_FromString(dot.col.c_str()));
		  PyList_SetItem(dots_list_py, ii, dot_py);
	       }
	       PyDict_SetItemString(dots_map_py, col_key.c_str(), dots_list_py);
	    }

	    for (std::size_t ii=0; ii<clashes.size(); ii++) {
	       PyObject *o_py = PyList_New(3);
	       PyObject *p1_py = PyList_New(3);
	       PyObject *p2_py = PyList_New(3);
	       double dist = std::sqrt((clashes.positions[ii].first-clashes.positions[ii].second).lengthsq());
	       PyList_SetItem(p1_py, 0, PyFloat_FromDouble(clashes.positions[ii].first.x()));
	       PyList_SetItem(p1_py, 1, PyFloat_FromDouble(clashes.positions[ii].first.y()));
	       PyList_SetItem(p1_py, 2, PyFloat_FromDouble(clashes.positions[ii].first.z()));
	       PyList_SetItem(p2_py, 0, PyFloat_FromDouble(clashes.positions[ii].second.x()));
	       PyList_SetItem(p2_py, 1, PyFloat_FromDouble(clashes.positions[ii].second.y()));
	       PyList_SetItem(p2_py, 2, PyFloat_FromDouble(clashes.positions[ii].second.z()));
	       PyList_SetItem(o_py, 0, p1_py);
	       PyList_SetItem(o_py, 1, p2_py);
	       PyList_SetItem(o_py, 2, PyFloat_FromDouble(dist));
	       PyList_SetItem(clashes_py, ii, o_py);
	    }

	    PyObject *typed_clashes_py = PyList_New(2);
	    // PyList_SetItem(typed_clashes_py, 0, myPyString_FromString(clashes.type.c_str()));
	    PyList_SetItem(typed_clashes_py, 0, PyUnicode_FromString("clashes"));
	    PyList_SetItem(typed_clashes_py, 1, clashes_py);
	    PyList_SetItem(atom_overlaps_py, 0, dots_map_py);
	    PyList_SetItem(atom_overlaps_py, 1, typed_clashes_py);
	    o = boost::python::object(boost::python::handle<>(atom_overlaps_py));
	 }
	 catch (const std::out_of_range &oor) {
	    std::cout << "ERROR:: " << oor.what() << std::endl;
	 }
      }
   }
   return o;
}

boost::python::object
coot::atom_overlaps_from_coordinates_file(const std::string &file_name,
					  bool add_hydrogens_flag) {

   boost::python::object obj(boost::python::handle<>(Py_False)); // return this

   protein_geometry geom;
   geom.set_verbose(false);
   geom.init_standard();

   if (coot::file_exists(file_name)) {
      mmdb::Manager *mol = new mmdb::Manager;
      mmdb::ERROR_CODE err = mol->ReadCoorFile(file_name.c_str());
      if (err) {
	 std::cout << "Error reading coordinates file " << file_name << std::endl;
      } else {
	 // Happy Path
	 try {

	    bool proceed = true; // unless we can't find a dictionary
	    int read_number = 40;
	    std::vector<std::string> rtv = coot::util::non_standard_residue_types_in_molecule(mol);
	    for (unsigned int i=0; i<rtv.size(); i++)
	       if (rtv[i] != "HOH") {
		  int status = geom.try_dynamic_add(rtv[i], read_number++);
		  if (status == 0) {
		     // failure
		     std::cout << "Failed to add dictionary for " << rtv[i]  << std::endl;
		     proceed = false;
		  }
	       }

	    if (proceed) {

	       // OK, we set a non-False return value

	       // spike-length probe-radius (unused)
	       bool ignore_waters = false;
	       coot::atom_overlaps_container_t overlaps(mol, &geom, ignore_waters, 0.5, 0.25);
	       overlaps.make_all_atom_overlaps();
	       std::vector<coot::atom_overlap_t> olv = overlaps.overlaps;
	       // std::cout << "Found " << olv.size() << " atom overlaps" << std::endl;
	       PyObject *o_py = PyList_New(olv.size());
	       for (std::size_t ii=0; ii<olv.size(); ii++) {
		  const coot::atom_overlap_t &o = olv[ii];
		  if (false) // debug
		     std::cout << "Overlap " << ii << " "
			       << coot::atom_spec_t(o.atom_1) << " "
			       << coot::atom_spec_t(o.atom_2) << " overlap-vol "
			       << o.overlap_volume << " r_1 "
			       << o.r_1 << " r_2 " << o.r_2 << std::endl;
		  PyObject *item_dict_py = PyDict_New();
		  py_atom_spec_t spec_1(coot::atom_spec_t(o.atom_1));
		  py_atom_spec_t spec_2(coot::atom_spec_t(o.atom_2));
		  PyObject *r_1_py = PyFloat_FromDouble(o.r_1);
		  PyObject *r_2_py = PyFloat_FromDouble(o.r_2);
		  PyObject *ov_py  = PyFloat_FromDouble(o.overlap_volume);
		  PyDict_SetItemString(item_dict_py, "atom-1-spec", spec_1.pyobject());
		  PyDict_SetItemString(item_dict_py, "atom-2-spec", spec_2.pyobject());
		  PyDict_SetItemString(item_dict_py, "overlap-volume", ov_py);
		  PyDict_SetItemString(item_dict_py, "radius-1", r_1_py);
		  PyDict_SetItemString(item_dict_py, "radius-2", r_2_py);
		  PyList_SetItem(o_py, ii, item_dict_py);
	       }
	       obj = boost::python::object(boost::python::handle<>(o_py));
	    }
	 }

	 catch (const std::runtime_error &rte) {
	    std::cout << "WARNING:: " << rte.what() << std::endl;
	 }
	 catch (const std::exception &e) {
	    std::cout << "WARNING:: " << e.what() << std::endl;
	 }
      }
   }
   return obj;
}



#endif // MAKE_ENHANCED_LIGAND_TOOLS
