
/*  ----------------------------------------------------------------------- */
/*                  Pisa internal                                           */
/*  ----------------------------------------------------------------------- */
clipper::Coord_orth
make_complementary_dotted_surfaces(int imol_1, int imol_2,
                                   std::vector<coot::residue_spec_t> &r1,
                                   std::vector<coot::residue_spec_t> &r2);
#ifdef USE_GUILE
std::vector<coot::residue_spec_t>
residue_records_list_scm_to_residue_specs(SCM mol_1_residue_records,
                                          const std::string &chain_id);
SCM symbol_value_from_record(SCM record_1, const std::string &symbol);
#endif
#ifdef USE_PYTHON
std::vector<coot::residue_spec_t>
residue_records_list_py_to_residue_specs(PyObject *mol_1_residue_records,
                                         const std::string &chain_id);
//PyObject *symbol_value_from_record(PyObject *record_1, const std::string &symbol);
#endif

void
add_generic_object_bond(int imol1, int imol2,
                        const coot::atom_spec_t &atom_spec_1,
                        const coot::atom_spec_t &atom_spec_2,
                        int generic_object_number,
                        const std::string &colour);

void
pisa_interfaces_display_only(int imol_1, int imol_2, clipper::Coord_orth centre_pt);
std::string untangle_mmdb_chain_id_string(const std::string &mmdb_chain_id_in);

