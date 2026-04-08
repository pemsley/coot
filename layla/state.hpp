/* layla/state.hpp
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

#ifndef LAYLA_STATE_HPP
#define LAYLA_STATE_HPP
#include "inchi_key_database.hpp"
#include "ligand_editor_canvas.hpp"
#include "notifier.hpp"
#include <rdkit/GraphMol/RWMol.h>

#include "geometry/protein-geometry.hh"
#include "ligand_editor_canvas/core.hpp"
#include "ligand_editor_canvas/model.hpp"

namespace coot::layla {


/// Structure holding the state of the editor
class LaylaState {

    public:
    enum class ExportMode: unsigned char {
        PDF,
        PNG,
        SVG
    };

    struct InchiKeyLookupResult {
        std::string monomer_id;
        std::string chemical_name;
    };
    private:

    enum class UnsavedChangesDialogPurpose: unsigned char {
        CloseEditor,
        FileNew
    };

    /// The widget. Owned by glib/gtk
    CootLigandEditorCanvas* canvas;
    /// Owned by glib/gtk
    GtkWindow* main_window;
    /// Owned by glib/gtk
    GtkLabel* status_label;
    /// Owned by us.
    CootLaylaNotifier* notifier;

    bool unsaved_changes;
    std::optional<UnsavedChangesDialogPurpose> unsaved_changes_dialog_purpose;
    std::optional<unsigned int> current_filesave_molecule;
    std::optional<std::string> current_filesave_filename;

    /// Store of monomer library information
    std::unique_ptr<protein_geometry> monomer_library_info_store;
    std::unique_ptr<InchiKeyDatabase> inchi_key_database;

    /// Adds the molecule to the canvas.
    /// This function takes ownership of the molecule pointer.
    int append_molecule(RDKit::RWMol* molecule_ptr);

    void update_status(const char* new_status) noexcept;

    /// The core logic of saving files
    void save_file(unsigned int idx, const char* filename, GtkWindow* parent = nullptr) noexcept;
    void run_file_save_dialog(unsigned int molecule_idx) noexcept;

    public:
    LaylaState(CootLigandEditorCanvas* canvas_widget, GtkWindow* main_window, GtkLabel* status_label = nullptr) noexcept;
    ~LaylaState() noexcept;
    
    /// Returns a new reference to the event notifier object
    CootLaylaNotifier* get_notifier() const noexcept;

    ///Useful for signal handlers
    CootLigandEditorCanvas* get_canvas() const noexcept;

    /// Uses the `inchi_key_database` to look up th given inchi key
    std::optional<InchiKeyLookupResult> lookup_inchi_key(const std::string& inchi_key) const noexcept;

    /// Clears all molecules, file names, etc.
    void reset() noexcept;

    bool has_unsaved_changes() const noexcept;
    void unsaved_changes_dialog_accepted();

    // TOOLS
    void run_choose_element_dialog();
    // FILE
    void load_from_smiles();
    void file_new();
    void file_save();
    void file_save_as();
    void file_open();
    void file_export(ExportMode mode);
    void file_import_molecule();
    void file_fetch_molecule();
    void file_exit();
    // EDIT
    void edit_undo();
    void edit_redo();
    // Display
    void switch_display_mode(ligand_editor_canvas::DisplayMode mode);
};

/// Let this be the singleton used by the editor executable.
/// Could by used by Coot as well.
inline LaylaState* global_instance = nullptr;

/// Global GtkBuilder created from `layla.ui`.
/// Used for accessing widgets inside dialogs 
/// from within signal handlers.
inline GtkBuilder* global_layla_gtk_builder = nullptr;

void initialize_global_instance(CootLigandEditorCanvas* canvas, GtkWindow* win, GtkLabel* status_label = nullptr);

} // namespace coot::layla


#endif
