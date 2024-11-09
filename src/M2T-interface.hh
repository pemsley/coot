
#ifndef COOT_SRC_M2T_INTERFACE_HH
#define COOT_SRC_M2T_INTERFACE_HH

#include <string>

//! Update float parameter for MoleculesToTriangles molecular mesh
void M2T_updateFloatParameter(int imol, const std::string &param_name, float value);

//! Update int parameter for MoleculesToTriangles molecular mesh
void M2T_updateIntParameter(int imol, const std::string &param_name, int value);

#endif // COOT_SRC_M2T_INTERFACE_HH
