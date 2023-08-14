#ifndef LIGAND_BUILDER_UTILS_HPP
#define LIGAND_BUILDER_UTILS_HPP
#include <rdkit/GraphMol/RWMol.h>

namespace coot::ligand_editor {

void remove_non_polar_hydrogens(RDKit::RWMol* mol);

}

#endif //  LIGAND_BUILDER_UTILS_HPP