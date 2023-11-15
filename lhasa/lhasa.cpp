#include "lhasa.hpp"
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>
#include "glog_replacement.hpp"

std::unique_ptr<RDKit::RWMol> lhasa::rdkit_mol_from_smiles(std::string smiles) {
    std::unique_ptr<RDKit::RWMol> ret(RDKit::SmilesToMol(smiles));
    return ret;
}

std::string lhasa::rdkit_mol_to_smiles(RDKit::ROMol& mol) {
    auto ret = RDKit::MolToSmiles(mol);
    return ret;
}

void lhasa::append_from_smiles(CootLigandEditorCanvas& canvas, std::string smiles) {
    canvas.append_molecule(rdkit_mol_from_smiles(smiles));
}

coot::ligand_editor_canvas::ActiveTool lhasa::make_active_tool(emscripten::val tool) {
    using namespace coot::ligand_editor_canvas;

    // Just returns 'object'
    //std::string type_name = tool.typeOf().as<std::string>();

    std::string classname = tool["__proto__"]["constructor"]["name"].as<std::string>();
    // std::cout<< classname << '\n';

    // ActiveTool(ElementInsertion insertion) noexcept;
    // ActiveTool(BondModifier modifier) noexcept;
    // ActiveTool(TransformTool) noexcept;
    // ActiveTool(StructureInsertion insertion) noexcept;
    // ActiveTool(FlipTool) noexcept;
    if(classname == "LhasaDeleteTool") {
        return ActiveTool(DeleteTool());
    }
    if(classname == "LhasaChargeModifier") {
        return ActiveTool(ChargeModifier());
    }
    if(classname == "LhasaGeometryModifier") {
        return ActiveTool(GeometryModifier());
    }
    if(classname == "LhasaFormatTool") {
        return ActiveTool(FormatTool());
    }
    if(classname == "LhasaRemoveHydrogensTool") {
        return ActiveTool(RemoveHydrogensTool());
    }

    g_critical("%s does not correspond to any known tool type. Returning empty ActiveTool.", classname.c_str());
    return ActiveTool();

}