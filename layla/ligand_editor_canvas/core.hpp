/* layla/ligand_editor_canvas/core.hpp
 * 
 * Copyright 2023 by Global Phasing Ltd.
 * Author: Jakub Smulski
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifndef COOT_LIGAND_EDITOR_CANVAS_CORE_HPP
#define COOT_LIGAND_EDITOR_CANVAS_CORE_HPP
#include <gtk/gtk.h>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>
#include <memory>
#include <vector>
#include "model.hpp"
#include "pango/pango-layout.h"
#include "tools.hpp"

// GObject declaration 
G_BEGIN_DECLS   

#define COOT_LIGAND_EDITOR_CANVAS_TYPE (coot_ligand_editor_canvas_get_type ())
G_DECLARE_FINAL_TYPE  (CootLigandEditorCanvas, coot_ligand_editor_canvas, COOT, COOT_LIGAND_EDITOR_CANVAS, GtkWidget)

G_END_DECLS

namespace coot::ligand_editor_canvas::impl {

inline guint status_updated_signal;
inline guint scale_changed_signal;
inline guint smiles_changed_signal;
inline guint molecule_deleted_signal;

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
    std::unique_ptr<std::vector<CanvasMolecule>> molecules;
    std::unique_ptr<std::vector<std::shared_ptr<RDKit::RWMol>>> rdkit_molecules;

    StateSnapshot(const WidgetCoreData& core_data);
};

struct Renderer {
    cairo_t* cr;
    PangoLayout* pango_layout;
    /// Takes ownership of the pointers
    Renderer(cairo_t*,PangoLayout*);
    ~Renderer();
};

/// Used for widget's struct as a base class.
/// Useful for exposing inner state to the active tool.
struct WidgetCoreData {
    typedef std::vector<std::unique_ptr<StateSnapshot>> StateStack;
    typedef std::pair<CanvasMolecule::AtomOrBond,unsigned int> AtomOrBondWithMolIdx;
    typedef std::optional<AtomOrBondWithMolIdx> MaybeAtomOrBondWithMolIdx;

    /// Max length of edition history
    const static unsigned int MAX_STATE_STACK_LENGTH;
    /// Numbers of elements to be removed from the `state_stack`
    /// when its' maximum length gets exceeded
    const static unsigned int STATE_STACK_TRIM_BATCH_SIZE;

    protected:

    /// Current position in the state_stack, counting from the back.
    /// -1 if we're "fresh", that is, no action has been undone
    int state_stack_pos;
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
    /// Bond being currently created via click'n'drag
    std::optional<CurrentlyCreatedBond> currently_created_bond;

    float scale;

    bool allow_invalid_molecules;

    DisplayMode display_mode;

    void render(Renderer&);

    /// Does Edit->Undo
    void undo_edition();

    /// Does Edit->Redo
    void redo_edition();

    /// Cancels the current edition
    /// and resets the state to the current snapshot.
    void rollback_current_edition();

    /// Checks if we're currently inside an edition operation
    bool is_in_edition();

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

    void delete_molecule_with_idx(unsigned int idx) noexcept;

    /// Emits 'status-updated' signal.
    void update_status(const gchar* status_text) const noexcept;

    std::string build_smiles_string() const;
};

/// This is the private struct for GObject
struct CootLigandEditorCanvasPriv : CootLigandEditorCanvasPrivBase, impl::WidgetCoreData {    

    std::unique_ptr<ActiveTool> active_tool;
};



}

#endif //#define COOT_LIGAND_EDITOR_CANVAS_CORE_HPP
