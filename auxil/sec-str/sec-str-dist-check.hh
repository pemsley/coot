
#include <string>
#include "mmdb_manager.h"

void distance_checks(CModel *model_p);
int SSE(CResidue *res);
CMMDBManager *get_mol(const std::string &filename);
float oxygen_check(CResidue *residue_p_1, CResidue *residue_p_2);




