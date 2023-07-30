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

class Tool {

    public:

    /// Called when the click coordinates do not correspond to anything on canvas
    virtual void on_blank_space_click(impl::WidgetCoreData& widget_data, int x, int y);

    /// Returns true if `on_bond_click()` or `on_atom_click()` (respectively to what's been clicked) 
    /// should be called next (and then lastly `after_molecule_click()`)
    virtual bool on_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>&, CanvasMolecule&);

    virtual void on_bond_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>&, CanvasMolecule&, CanvasMolecule::Bond&);
    virtual void on_atom_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>&, CanvasMolecule&, CanvasMolecule::Atom&);

    /// Generic on-mouse-release event handler.
    /// No dedicated molecule/atom/bond handlers seem to be needed now.
    virtual void on_release(impl::WidgetCoreData& widget_data, int x, int y);

    virtual void after_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>&, CanvasMolecule&);

    /// Used to print tool-specific error messages should any handler throw an exception
    virtual std::string get_exception_message_prefix() const noexcept;

    /// Wraps `RDKit::MolOps::sanitizeMol()`.
    /// Used for checking validity and reversing kekulization.
    /// Does nothing when nonsensical molecules are allowed.
    static void sanitize_molecule(impl::WidgetCoreData&, RDKit::RWMol&);

    virtual ~Tool();
};

class BondModifier : public Tool {
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

    virtual bool on_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol) override;
    virtual void on_bond_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol, CanvasMolecule::Bond& bond) override;
    virtual void on_atom_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol, CanvasMolecule::Atom& atom) override;
    virtual std::string get_exception_message_prefix() const noexcept override;
    virtual void on_release(impl::WidgetCoreData& widget_data, int x, int y) override;
};

class ElementInsertion : public Tool {
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

    virtual bool on_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol) override;
    virtual void on_bond_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol, CanvasMolecule::Bond& bond) override;
    virtual void on_atom_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol, CanvasMolecule::Atom& atom) override;
    virtual void after_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol) override;
    virtual std::string get_exception_message_prefix() const noexcept override;
};



class StructureInsertion : public Tool {
    public:

    enum class Structure: unsigned int {
        CycloPropaneRing,
        CycloButaneRing,
        CycloPentaneRing,
        CycloHexaneRing,
        BenzeneRing,
        CycloHeptaneRing,
        CycloOctaneRing
    };
    private:
    Structure structure;

    static unsigned int append_carbon(RDKit::RWMol*, unsigned int target_idx, RDKit::Bond::BondType bond_type = RDKit::Bond::SINGLE);
    static unsigned int append_carbon_chain(RDKit::RWMol*, unsigned int chain_start_idx, std::size_t atom_count);
    void append_structure_to_atom(RDKit::RWMol*, unsigned int atom_idx) const;

    public:
    StructureInsertion(Structure) noexcept;

    Structure get_structure() const noexcept;

    virtual bool on_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol) override;
    virtual void on_bond_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol, CanvasMolecule::Bond& bond) override;
    virtual void on_atom_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol, CanvasMolecule::Atom& atom) override;
    virtual void after_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol) override;
    virtual void on_blank_space_click(impl::WidgetCoreData& widget_data, int x, int y) override;
    virtual std::string get_exception_message_prefix() const noexcept override;


};

class ChargeModifier : public Tool {
    public:

    virtual void on_atom_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol, CanvasMolecule::Atom& atom) override;
    virtual std::string get_exception_message_prefix() const noexcept override;

};

class DeleteTool : public Tool {
    public:

    virtual bool on_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol) override;
    virtual void on_bond_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol, CanvasMolecule::Bond& bond) override;
    virtual void on_atom_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol, CanvasMolecule::Atom& atom) override;
    virtual void after_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol) override;
    virtual std::string get_exception_message_prefix() const noexcept override;
};

class MoveTool {
    // This is an empty type.
    // It exists for the sake of maintaning API convention.
};

class RotateTool {
    // This is an empty type.
    // It exists for the sake of maintaning API convention.
};

class GeometryModifier : public Tool {
    public:

    virtual void on_bond_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol, CanvasMolecule::Bond& bond) override;
    virtual std::string get_exception_message_prefix() const noexcept override;
};

class FormatTool : public Tool {
    public:

    virtual bool on_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol) override;
    virtual std::string get_exception_message_prefix() const noexcept override;
};

class FlipTool : public Tool {
    FlipMode mode;
    public:
    FlipTool(FlipMode) noexcept;

    virtual bool on_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol) override;
    virtual std::string get_exception_message_prefix() const noexcept override;
};

class RemoveHydrogensTool : public Tool {
    public:

    virtual bool on_molecule_click(impl::WidgetCoreData& widget_data, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol) override;
    virtual std::string get_exception_message_prefix() const noexcept override;
};

class ActiveTool {
    public:
    enum class Mode: unsigned char {
        Rotate,
        Move,
        Tool
    };

    private:

    Mode mode;
    /// Non-owning pointer
    impl::WidgetCoreData* widget_data;
    std::unique_ptr<Tool> tool;
    TransformManager transform_manager;

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

    void on_click(int x, int y);
    void on_release(int x, int y);

    bool is_creating_bond() const noexcept;
    std::optional<std::pair<unsigned int,unsigned int>> get_molecule_idx_and_first_atom_of_new_bond() const noexcept;
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

    /// Only meant to be invoked from within CootLigandEditorCanvas implementation
    ///
    /// Sets the pointer to access widget's data.
    void set_core_widget_data(impl::CootLigandEditorCanvasPriv* owning_widget) noexcept;
};

}

#endif // COOT_LIGAND_EDITOR_CANVAS_TOOLS_HPP