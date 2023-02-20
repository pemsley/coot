#ifndef COOT_LIGAND_EDITOR_CANVAS_TOOLS_HPP
#define COOT_LIGAND_EDITOR_CANVAS_TOOLS_HPP
#include "model.hpp"

namespace coot::ligand_editor_canvas {

namespace impl {
    /// This is a forward declaration of a structure defined at
    /// "core.hpp"
    struct WidgetCoreData;
    /// This is a forward declaration of a structure defined at
    /// "core.hpp"
    struct CootLigandEditorCanvasPriv;
}

class BondModifier {
    public:
    enum class BondModifierMode {
        Single,
        Double,
        Triple
    };
    private:

    BondModifierMode mode;
    public:

};

class ElementInsertion {
    public:
    enum class Element {
        C,
        N,
        O,
        S,
        P,
        H,
        F,
        Cl,
        Br,
        I,
        // todo: what is this? an arbitrary element?
        // [like eg. Selene for organoselene compounds]
        X
    };

    private:
    Element element;
    public:
    ElementInsertion(Element el) noexcept;

    Element get_element() const noexcept;
    const char* get_element_symbol() const noexcept;
};



class StructureInsertion {
    public:

    enum class Structure: unsigned int {
        CycloPropaneRing,
        CycloButaneRing,
        CycloPentaneRing,
        CycloHexaneRing,
        BenzeneRing,
        CycloHeptaneRing,
        CycloOctaneRing,
        // todo:
        // "env residues"
        // "key"

    };
    private:
    Structure structure;

    public:

};


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
    /// Non-owning pointer
    impl::WidgetCoreData* widget_data;

    /// Checks if the internal variant (the kind of the tool) matches what's expected (passed as argument).
    /// Throws an exception in case of a mismatch.
    void check_variant(Variant);

    public:
    ActiveTool() noexcept;
    ActiveTool(ElementInsertion insertion) noexcept;
    ActiveTool(BondModifier modifier) noexcept;

    Variant get_variant() const noexcept;
    /// Valid for Variant::ElementInsertion.
    /// Inserts currently chosen atom at the given coordinates.
    void insert_atom(int x, int y);

    /// Only meant to be invoked from within CootLigandEditorCanvas implementation
    ///
    /// Sets the pointer to access widget's data.
    void set_core_widget_data(impl::CootLigandEditorCanvasPriv* owning_widget) noexcept;
};

}

#endif // COOT_LIGAND_EDITOR_CANVAS_TOOLS_HPP