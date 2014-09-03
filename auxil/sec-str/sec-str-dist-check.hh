
#include <string>
#include <mmdb2/mmdb_manager.h>

void distance_checks(mmdb::Model *model_p);
int SSE(mmdb::Residue *res);
CMMDBManager *get_mol(const std::string &filename);
float oxygen_check(mmdb::Residue *residue_p_1, mmdb::Residue *residue_p_2);




