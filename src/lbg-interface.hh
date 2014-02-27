
// 20140226: now we change things so that the interface functions always get generated,
//           and what the function does depends on MAKE_ENHANCED_LIGAND_TOOLS
// 
// #ifdef MAKE_ENHANCED_LIGAND_TOOLS
void residue_to_ligand_builder(int imol, const char *chain_id, int resno, const char *ins_code, double weight_for_3d_distances);
void smiles_to_ligand_builder(const char *smiles_string);
// #endif

