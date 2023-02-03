#ifndef COOT_LIGAND_EDITOR_CANVAS_TOOLS_HPP
#define COOT_LIGAND_EDITOR_CANVAS_TOOLS_HPP
#include "model.hpp"

namespace coot::ligand_editor_canvas {



class ActiveTool {
    public:
    enum class Variant: unsigned char {
        None,
        BondModifier,
        StructureInsertion,
        ElementInsertion,
        /// Stereo out
        GeometryModifier,
        DeleteHydrogens,
        Delete,
        Format,
        ChargeModifier
    };

    private:
    union {
        /// Valid for Variant::BondModifier
        BondModifier bond_modifier;
        /// Valid for Variant::ElementInsertion
        ElementInsertion element_insertion;
        /// Valid for Variant::StructureInsertion
        StructureInsertion structure_insertion;
    };
    Variant variant;

    public:
    ActiveTool() noexcept;

    Variant get_variant() const noexcept;

};

}

#endif // COOT_LIGAND_EDITOR_CANVAS_TOOLS_HPP