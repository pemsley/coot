
#include "protein_db/protein_db.h"  // Kevin's

mmdb::Manager *make_mol(const std::vector<ProteinDB::Chain> &chains, const std::string &chain_id, 
		       int first_resno, bool preserve_residue_names);
mmdb::Manager *make_mol(const ProteinDB::Chain &chain, const std::string &chain_id, int first_resno,
		       bool preserve_residue_names);
std::vector<mmdb::Residue *> 
add_chain_to_molecule(const ProteinDB::Chain &chain, const std::string &chain_id, 
		      int first_res_no, bool preserve_residue_names, mmdb::Manager *mol);

void add_cbs_and_os(std::vector<mmdb::Residue *> needs_cb_and_o, 
		    mmdb::Manager *mol);
