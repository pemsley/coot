
bool residue_to_sdf_file(int imol, const char *chain_id, int res_no, const char *ins_code, 
			 const char *sdf_file_name, bool kekulize);

bool residue_to_mdl_file_for_mogul(int imol, const char *chain_id,
				   int res_no, const char *ins_code, 
				   const char *mdl_file_name);

bool show_feats(int imol, const char *chain_id, int resno, const char *ins_code);

// rdkit chemical features.
// now wrappes generate_meshes()
bool show_feats(int imol, mmdb::Residue *residue_p, const coot::protein_geometry &geom);

