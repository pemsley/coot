
#include "protein_db/protein_db.h"  // Kevin's

CMMDBManager *make_mol(const std::vector<ProteinDB::Chain> &chains, const std::string &chain_id, 
		       int first_resno);
CMMDBManager *make_mol(const ProteinDB::Chain &chain, const std::string &chain_id, int first_resno);
std::vector<CResidue *> 
add_chain_to_molecule(const ProteinDB::Chain &chain, const std::string &chain_id, 
		      int first_res_no, CMMDBManager *mol);

void add_cbs_and_os(std::vector<CResidue *> needs_cb_and_o, 
		    CMMDBManager *mol);
