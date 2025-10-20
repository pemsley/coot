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

#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>
// To check if Inchi support is enabled
#include <rdkit/RDGeneral/RDConfig.h>
#ifdef RDK_BUILD_INCHI_SUPPORT
#include <rdkit/GraphMol/inchi.h>
#else
#warning Your version of RDKit was built without InChI support. Molecule InChI key lookup will not be available.
#endif

#include <memory>
#include <vector>
#include <optional>
#include "render.hpp"
#include "model.hpp"
#include "tools.hpp"

#ifndef __EMSCRIPTEN__
    #include <gtk/gtk.h>
    #include <pango/pango-layout.h>

    // GObject declaration 
    G_BEGIN_DECLS   

    #define COOT_LIGAND_EDITOR_CANVAS_TYPE (coot_ligand_editor_canvas_get_type ())
    G_DECLARE_FINAL_TYPE  (CootLigandEditorCanvas, coot_ligand_editor_canvas, COOT, COOT_LIGAND_EDITOR_CANVAS, GtkWidget)

    G_END_DECLS

    #define _LIGAND_EDITOR_SIGNAL_EMIT(_instance_, _signal_name_) g_signal_emit((gpointer)(_instance_),impl::_signal_name_,0)
    #define _LIGAND_EDITOR_SIGNAL_EMIT_ARG(_instance_, _signal_name_, ...) g_signal_emit((gpointer)(_instance_),impl::_signal_name_,0,__VA_ARGS__)
#else // __EMSCRIPTEN__ defined
    #include "../../lhasa/glog_replacement.hpp"
    #include <sigc++-3.0/sigc++/sigc++.h>
    #define _LIGAND_EDITOR_SIGNAL_EMIT(_instance_, _signal_name_) _instance_->_signal_name_.emit()
    #define _LIGAND_EDITOR_SIGNAL_EMIT_ARG(_instance_, _signal_name_, ...) _instance_->_signal_name_.emit(__VA_ARGS__)
#endif


namespace coot::ligand_editor_canvas::impl {

#ifndef __EMSCRIPTEN__
inline guint status_updated_signal;
inline guint scale_changed_signal;
inline guint smiles_changed_signal;
inline guint molecule_deleted_signal;
inline guint qed_info_updated_signal;
#endif

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
    #ifndef __EMSCRIPTEN__
    GtkWidget parent;
    #else // __EMSCRIPTEN__ defined
    // Lhasa-specific includes/definitions
    #endif
};


struct StateSnapshot {
    std::unique_ptr<std::vector<std::optional<CanvasMolecule>>> molecules;
    std::unique_ptr<std::vector<std::optional<std::shared_ptr<RDKit::RWMol>>>> rdkit_molecules;

    StateSnapshot(const WidgetCoreData& core_data);
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
    std::unique_ptr<std::vector<std::optional<CanvasMolecule>>> molecules;
    /// molecules (RDKit)
    std::unique_ptr<std::vector<std::optional<std::shared_ptr<RDKit::RWMol>>>> rdkit_molecules;
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

    void delete_molecule_with_idx(unsigned int idx, bool integrate_with_edit_undo = true) noexcept;

    /// Emits 'status-updated' signal.
    void update_status(const char* status_text) const noexcept;

    /// Emits the 'smiles-changed' and 'qed-info-updated' signal
    void emit_mutation_signals() const noexcept;

    coot::ligand_editor_canvas::SmilesMap build_smiles() const;
    coot::ligand_editor_canvas::InchiKeyMap build_inchi_keys() const;

    unsigned int get_molecule_count_impl() const noexcept;
    /// Returns -1 if none
    int get_first_molecule_idx() const noexcept;

    /// Abstraction over gtk_widget_queue_draw
    void queue_redraw() const noexcept;
    /// Abstraction over gtk_widget_queue_resize
    void queue_resize() const noexcept;
};

/// This is the private struct for GObject
struct CootLigandEditorCanvasPriv : CootLigandEditorCanvasPrivBase, impl::WidgetCoreData {    

    std::unique_ptr<ActiveTool> active_tool;
};


} // namespace coot::ligand_editor_canvas::impl

#ifdef __EMSCRIPTEN__

//Forward declaration
namespace emscripten {
    struct val;
}

/// For Lhasa
struct CootLigandEditorCanvas : coot::ligand_editor_canvas::impl::CootLigandEditorCanvasPriv {

    
    sigc::signal<void(const char*)> status_updated_signal;
    sigc::signal<void(float)> scale_changed_signal;
    sigc::signal<void()> smiles_changed_signal;
    sigc::signal<void(int)> molecule_deleted_signal;
    sigc::signal<void(int, const coot::ligand_editor_canvas::CanvasMolecule::QEDInfo*)> qed_info_updated_signal;
    // Lhasa-only signals (for JS handlers):
    sigc::signal<void()> queue_redraw_signal;
    sigc::signal<void()> queue_resize_signal;

    public:

    enum class MeasurementDirection :unsigned char {
        HORIZONTAL = 0,
        VERTICAL = 1
    };

    struct SizingInfo {
        int requested_size;
    };

    // Implemented at 'ligand_editor_canvas.cpp'
    CootLigandEditorCanvas() noexcept;
    // Implemented at 'ligand_editor_canvas.cpp'
    ~CootLigandEditorCanvas() noexcept;

    void set_active_tool(std::unique_ptr<coot::ligand_editor_canvas::ActiveTool> active_tool);
    /// Returns the ID of new molecule or '-1' on error
    int append_molecule(std::shared_ptr<RDKit::RWMol> rdkit_mol) noexcept;
    void update_molecule_from_smiles(unsigned int molecule_idx, const std::string& smiles);
    void set_scale(float scale) noexcept;
    float get_scale() noexcept;
    void undo() noexcept;
    void redo() noexcept;
    // const RDKit::ROMol& get_rdkit_molecule(unsigned int index) noexcept;
    unsigned int get_molecule_count() noexcept;
    unsigned int get_idx_of_first_molecule() noexcept;
    unsigned int get_max_molecule_idx() noexcept;
    void set_allow_invalid_molecules(bool value) noexcept;
    bool get_allow_invalid_molecules() noexcept;
    coot::ligand_editor_canvas::DisplayMode get_display_mode() noexcept;
    void set_display_mode(coot::ligand_editor_canvas::DisplayMode value) noexcept;
    coot::ligand_editor_canvas::SmilesMap get_smiles() noexcept;
    coot::ligand_editor_canvas::InchiKeyMap get_inchi_keys() noexcept;
    std::string get_smiles_for_molecule(unsigned int molecule_idx) noexcept;
    std::string get_inchi_key_for_molecule(unsigned int molecule_idx) noexcept;
    std::string get_pickled_molecule(unsigned int molecule_idx) noexcept;
    std::string get_pickled_molecule_base64(unsigned int molecule_idx) noexcept;
    void clear_molecules() noexcept;

    /// For connecting javascript handlers to signals
    void connect(std::string signal_name, emscripten::val callback);

    // Implemented at 'ligand_editor_canvas.cpp'
    SizingInfo measure(MeasurementDirection orientation) const noexcept;
    // Implemented at 'ligand_editor_canvas.cpp'
    void on_hover(double x, double y, bool alt_pressed, bool control_pressed);
    // Implemented at 'ligand_editor_canvas.cpp'
    void on_scroll(double dx, double dy, bool control_pressed);
    // Implemented at 'ligand_editor_canvas.cpp'
    void on_left_click(double x, double y, bool alt_pressed, bool control_pressed, bool shift_pressed);
    // Implemented at 'ligand_editor_canvas.cpp'
    void on_left_click_released(double x, double y, bool alt_pressed, bool control_pressed, bool shift_pressed);
    // Implemented at 'ligand_editor_canvas.cpp'
    void on_right_click(double x, double y, bool alt_pressed, bool control_pressed, bool shift_pressed);
    // Implemented at 'ligand_editor_canvas.cpp'
    void on_right_click_released(double x, double y, bool alt_pressed, bool control_pressed, bool shift_pressed);


};
#endif


#endif //#define COOT_LIGAND_EDITOR_CANVAS_CORE_HPP
