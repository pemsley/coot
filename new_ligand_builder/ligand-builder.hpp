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

    /// the widget
    CootLigandEditorCanvas* canvas;
    GtkWindow* main_window;

    /// Updates the widget content with the molecule.
    /// This function takes ownership of the molecule pointer.
    void set_molecule(RDKit::RWMol* molecule_ptr);

    public:
    LigandBuilderState(CootLigandEditorCanvas* canvas_widget, GtkWindow* main_window) noexcept;
    
    // FILE
    void load_from_smiles();
    void file_new();
    void file_save();
    void file_open();
    void file_export(ExportMode mode);
    // EDIT
    void edit_undo();
    void edit_redo();
    // Display
    // todo

};

/// Let this be the singleton used by the editor executable
inline LigandBuilderState* global_instance;

void initialize_global_instance(CootLigandEditorCanvas* canvas, GtkWindow* win);

} // namespace coot::ligand_editor


#endif