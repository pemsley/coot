#ifndef COOT_LIGAND_EDITOR_CANVAS_CORE_HPP
#define COOT_LIGAND_EDITOR_CANVAS_CORE_HPP
#include "gtk/gtk.h"
#include <rdkit/GraphMol/RWMol.h>
#include <memory>
#include <vector>
#include "model.hpp"
#include "tools.hpp"

// GObject declaration 
G_BEGIN_DECLS   

#define COOT_LIGAND_EDITOR_CANVAS_TYPE (coot_ligand_editor_canvas_get_type ())
G_DECLARE_FINAL_TYPE  (CootLigandEditorCanvas, coot_ligand_editor_canvas, COOT, COOT_LIGAND_EDITOR_CANVAS, GtkWidget)

G_END_DECLS

namespace coot::ligand_editor_canvas::impl {

inline guint status_updated_signal;
inline guint scale_changed_signal;

/// This is here as a workaround.
/// 
/// GObject expects the first data field in the struct
/// to be an instance of the base class.
/// If this condition is not met, you end up with a SEGFAULT.
///
/// Thus CootLigandEditorCanvasPriv must first inherit from
/// a struct that satisfies this requirement.
/// This is how you can use (multiple) inheritance while
/// keeping the ABI happy.
struct CootLigandEditorCanvasPrivBase {
    GtkWidget parent;
};


struct StateSnapshot {

};

/// Used for widget's struct as a base class.
/// Useful for exposing inner state to the active tool.
struct WidgetCoreData {
    typedef std::vector<StateSnapshot> StateStack;
    typedef std::pair<CanvasMolecule::AtomOrBond,unsigned int> AtomOrBondWithMolIdx;
    typedef std::optional<AtomOrBondWithMolIdx> MaybeAtomOrBondWithMolIdx;

    protected:

    /// For Edit->Undo/Redo.
    /// To remember internal states
    std::unique_ptr<StateStack> state_stack;

    /// A snapshot preserving internal state
    /// from before the current edition began.
    /// nullptr if no edition is being done at the moment.
    std::unique_ptr<StateSnapshot> state_before_edition;

    public:
    /// molecules on the screen
    std::unique_ptr<std::vector<CanvasMolecule>> molecules;
    /// molecules (RDKit)
    std::unique_ptr<std::vector<std::shared_ptr<RDKit::RWMol>>> rdkit_molecules;

    float scale;

    public:

    /// Does Edit->Undo
    void undo_edition();

    /// Does Edit->Redo
    void redo_edition();

    /// Cancels the current edition
    /// and resets the state to the current snapshot.
    void rollback_current_edition();

    /// Snapshots the current state and opens new edition.
    ///
    /// This function must be called if one wishes
    /// to integrate any kind of state-altering operation
    /// with the Edit->Undo/Redo system
    void begin_edition();

    /// Completes the current edition.
    /// Moves the current snapshot to the state history stack.
    ///
    /// This function must be called after any kind 
    /// of state-altering operation if one wishes 
    /// to integrate it with the Edit->Undo/Redo system
    void finalize_edition();

    /// Goes over all molecules stored in the widget
    /// and calls CanvasMolecule::resolve_click(x,y) on each of them
    /// until an object matching the coordinates is found.
    /// The index number indicates which molecule the object comes from.
    /// If nothing matches the coordinates, nullopt is returned.
    MaybeAtomOrBondWithMolIdx resolve_click(int x, int y) const noexcept;

    /// Emits 'status-updated' signal.
    void update_status(const gchar* status_text) const noexcept;
};

/// This is the private struct for GObject
struct CootLigandEditorCanvasPriv : CootLigandEditorCanvasPrivBase, impl::WidgetCoreData {    

    std::unique_ptr<ActiveTool> active_tool;
};



}

#endif //#define COOT_LIGAND_EDITOR_CANVAS_CORE_HPP
