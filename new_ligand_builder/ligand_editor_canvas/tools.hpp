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
    std::optional<std::pair<unsigned int,unsigned int>> molecule_idx_and_first_atom_of_new_bond;
    bool is_in_drag;

    public:
    BondModifier(BondModifierMode) noexcept;

    CanvasMolecule::BondType get_target_bond_type() const noexcept;
    bool is_creating_bond() const noexcept;
    void begin_creating_bond(unsigned int molecule_idx, unsigned int atom_idx) noexcept;
    void finish_creating_bond() noexcept;
    std::optional<std::pair<unsigned int,unsigned int>> get_molecule_idx_and_first_atom_of_new_bond() const noexcept;

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
        I
    };

    private:
    std::variant<Element,unsigned int> element;
    public:
    ElementInsertion(const char* symbol);
    ElementInsertion(Element el) noexcept;

    unsigned int get_atomic_number() const noexcept;
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

    };
    private:
    Structure structure;

    public:
    StructureInsertion(Structure) noexcept;

    Structure get_structure() const noexcept;

};

class ChargeModifier {

};

class DeleteTool {

};

class TransformManager {
    public:
    enum class Mode {
        Rotation,
        Translation
    };
    private:

    struct RotationState {
        double last_absolute_angle;
        std::pair<int,int> original_rotation_pos;
        std::pair<int,int> current_rotation_pos;

        double get_current_absolute_angle(bool snap_to_angle) const;
        double get_current_angle_diff(bool snap_to_angle) const;
    };
    struct TranslationState {
        std::pair<int,int> prev_move_pos;
        std::pair<int,int> current_move_pos;

        std::pair<int,int> get_current_offset() const;
    };

    class IdleState {
        // empty
    };

    std::variant<RotationState, TranslationState, IdleState> state;
    std::optional<unsigned int> canvas_mol_idx;

    public:

    TransformManager() noexcept;

    bool is_active() const noexcept;

    void begin_transform(int x, int y, Mode mode) noexcept;
    void update_current_cursor_pos(int x, int y, bool snap) noexcept;

    void end_transform() noexcept;
    
    void set_canvas_molecule_index(unsigned int) noexcept;

    void apply_current_transform_state(impl::WidgetCoreData* widget_data, bool snapt_to_angle, bool about_to_end) const;
};

class MoveTool {

};

class RotateTool {
    
};

class GeometryModifier {

};

class FormatTool {

};

class FlipTool {
    FlipMode mode;
    public:
    FlipTool(FlipMode) noexcept;
    FlipMode get_mode() const noexcept;
};

class RemoveHydrogensTool {

};

class ActiveTool {
    public:
    enum class Variant: unsigned char {
        None,
        MoveTool,
        BondModifier,
        StructureInsertion,
        ElementInsertion,
        /// Stereo out
        GeometryModifier,
        Delete,
        Format,
        ChargeModifier,
        RotateTool,
        FlipTool,
        RemoveHydrogens
    };

    private:
    union {
        /// Valid for Variant::BondModifier
        BondModifier bond_modifier;
        /// Valid for Variant::ElementInsertion
        ElementInsertion element_insertion;
        /// Valid for Variant::StructureInsertion
        StructureInsertion structure_insertion;
        /// Valid for Variant::ChargeModifier
        ChargeModifier charge_modifier;
        /// Valid for Variant::Delete
        DeleteTool delete_tool;
        /// Valid for Variant::MoveTool
        MoveTool move_tool;
        /// Valid for Variant::GeometryModifier
        GeometryModifier geometry_modifier;
        /// Valid for Variant::Format
        FormatTool format_tool;
        /// Valid for Variant::RotateTool
        RotateTool rotate_tool;
        /// Valid for Variant::FlipTool
        FlipTool flip_tool;
        /// Valid for Variant::RemoveHydrogens
        RemoveHydrogensTool remove_hydrogens_tool;
    };
    Variant variant;
    /// Non-owning pointer
    impl::WidgetCoreData* widget_data;
    TransformManager transform_manager;

    /// Checks if the internal variant (the kind of the tool) matches what's expected (passed as argument).
    /// Throws an exception in case of a mismatch.
    void check_variant(Variant) const;

    /// Wraps `RDKit::MolOps::sanitizeMol()`.
    /// Used for checking validity and reversing kekulization.
    /// Does nothing when nonsensical molecules are allowed.
    void sanitize_molecule(RDKit::RWMol&) const;

    public:
    ActiveTool() noexcept;
    ActiveTool(ElementInsertion insertion) noexcept;
    ActiveTool(BondModifier modifier) noexcept;
    ActiveTool(DeleteTool) noexcept;
    ActiveTool(ChargeModifier) noexcept;
    ActiveTool(MoveTool) noexcept;
    ActiveTool(StructureInsertion insertion) noexcept;
    ActiveTool(GeometryModifier modifier) noexcept;
    ActiveTool(FormatTool) noexcept;
    ActiveTool(RotateTool) noexcept;
    ActiveTool(FlipTool) noexcept;
    ActiveTool(RemoveHydrogensTool) noexcept;

    Variant get_variant() const noexcept;
    /// Valid for Variant::ElementInsertion.
    /// Inserts currently chosen atom at the given coordinates.
    void insert_atom(int x, int y);
    /// Valid for Variant::BondModifier.
    /// Changes the bond found at the given coordinates.
    /// The kind of the bond depends upon current BondModifierMode.
    void alter_bond(int x, int y);
    bool is_creating_bond() const noexcept;
    /// Valid for Variant::BondModifier.
    void finish_creating_bond(int x, int y);
    /// Valid for Variant::ChargeModifier.
    /// Modifies the charge of bond found at the given coordinates.
    void alter_charge(int x, int y);
    /// Valid for Variant::Delete.
    /// Deletes whatever is found at the given coordinates
    void delete_at(int x, int y);
    /// Valid for Variant::StructureInsertion.
    /// Inserts currently chosen structure at the given coordinates.
    void insert_structure(int x, int y);
    /// Updates the current mouse coordinates
    /// allowing to compute adequate viewport rotation / translation
    void update_transform_cursor_pos(int x, int y, bool snap_to_angle) noexcept;
    /// Ends move or rotation
    void end_transform(bool snap_to_angle);
    /// Begins rotation or translation
    void begin_transform(int x, int y, TransformManager::Mode);
    /// Returns if the user is currently dragging their mouse
    /// to shift the viewport / rotate the molecule.
    bool is_in_transform() const noexcept;
    /// Valid for Variant::GeometryModifier.
    /// Changes geometry of the bond found at the given coordinates.
    void alter_geometry(int x, int y);
    /// Valid for Variant::Format.
    /// Forces re-computation of molecule geometry from scratch using RDKit.
    void format_at(int x, int y);
    /// Valid for Variant::FlipTool.
    /// Flip molecule under cursor around X/Y axis.
    void flip(int x, int y);
    /// Valid for Variant::BondModifier.
    std::optional<std::pair<unsigned int,unsigned int>> get_molecule_idx_and_first_atom_of_new_bond() const;

    /// Only meant to be invoked from within CootLigandEditorCanvas implementation
    ///
    /// Sets the pointer to access widget's data.
    void set_core_widget_data(impl::CootLigandEditorCanvasPriv* owning_widget) noexcept;
};

}

#endif // COOT_LIGAND_EDITOR_CANVAS_TOOLS_HPP