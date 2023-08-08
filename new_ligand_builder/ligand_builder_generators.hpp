#ifndef LIGAND_BUILDER_GENERATORS_HPP
#define LIGAND_BUILDER_GENERATORS_HPP
#include <string>

namespace coot::ligand_editor {


struct GeneratorRequest {
    enum class InputFormat: unsigned char {
        SMILES,
        MolFile
    } input_format;
    enum class Generator: unsigned char {
        Acedrg,
        Grade2
    } generator;
    std::string monomer_id;
};

} // namespace coot::ligand_editor

#endif // LIGAND_BUILDER_GENERATORS_HPP