#ifndef LIGAND_BUILDER_GENERATORS_HPP
#define LIGAND_BUILDER_GENERATORS_HPP
#include <string>
#include <gio/gio.h>

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
    std::string molecule_smiles;

    std::string get_filename() const;
};

inline GCancellable* global_generator_request_task_cancellable;

GCancellable* run_generator_request(GeneratorRequest request);

} // namespace coot::ligand_editor

#endif // LIGAND_BUILDER_GENERATORS_HPP