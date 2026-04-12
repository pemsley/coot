
#include <string>
#include <vector>

std::vector<int> make_ncs_maps(int imol_model, int imol_map);

//! \brief Make NCS-related copies of a ligand
//!
//! Given a ligand near the NCS master chain, use the NCS operators to generate
//! copies of the ligand positioned near each peer chain.
//! Returns a vector of molecule numbers for the newly created ligand copies.
std::vector<int> ncs_ligand(int imol_protein,
                            const std::string &ncs_master_chain_id,
                            int imol_ligand,
                            const std::string &chain_id_ligand,
                            int resno_ligand_start,
                            int resno_ligand_stop);
