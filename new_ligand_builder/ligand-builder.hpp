#ifndef LIGAND_BUILDER_HPP
#define LIGAND_BUILDER_HPP
#include "ligand_editor_canvas.hpp"
#include <rdkit/GraphMol/RWMol.h>

namespace coot::ligand_editor {



/// Structure holding the state of the editor
class LigandBuilderState {

    public:
    enum class ExportMode: unsigned char {
        PDF,
        PNG,
        SVG
    };
    private:

    /// The widget. Owned by glib/gtk
    CootLigandEditorCanvas* canvas;
    /// Owned by glib/gtk
    GtkWindow* main_window;
    /// Owned by glib/gtk
    GtkLabel* status_label;

    std::optional<unsigned int> current_filesave_molecule;
    std::optional<std::string> current_filesave_filename;

    /// Adds the molecule to the canvas.
    /// This function takes ownership of the molecule pointer.
    void append_molecule(RDKit::RWMol* molecule_ptr);

    void update_status(const char* new_status) noexcept;

    /// The core logic of saving files
    void save_file(unsigned int idx, const char* filename, GtkWindow* parent = nullptr) noexcept;
    void run_file_save_dialog(unsigned int molecule_idx) noexcept;

    public:
    LigandBuilderState(CootLigandEditorCanvas* canvas_widget, GtkWindow* main_window, GtkLabel* status_label = nullptr) noexcept;
    
    // FILE
    void load_from_smiles();
    void file_new();
    void file_save();
    void file_save_as();
    void file_open();
    void file_export(ExportMode mode);
    void file_import_molecule();
    void file_fetch_molecule();
    // EDIT
    void edit_undo();
    void edit_redo();
    // Display
    // todo

};

/// Let this be the singleton used by the editor executable.
/// Could by used by Coot as well.
inline LigandBuilderState* global_instance;

void initialize_global_instance(CootLigandEditorCanvas* canvas, GtkWindow* win, GtkLabel* status_label = nullptr);

} // namespace coot::ligand_editor


#endif