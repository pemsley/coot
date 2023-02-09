#ifndef COOT_LIGAND_EDITOR_CANVAS_CORE_HPP
#define COOT_LIGAND_EDITOR_CANVAS_CORE_HPP
#include "gtk/gtk.h"
#include <rdkit/GraphMol/RWMol.h>
#include <memory>
#include <vector>
#include "tools.hpp"

// GObject declaration 
G_BEGIN_DECLS   

#define COOT_LIGAND_EDITOR_CANVAS_TYPE (coot_ligand_editor_canvas_get_type ())
G_DECLARE_FINAL_TYPE  (CootLigandEditorCanvas, coot_ligand_editor_canvas, COOT, COOT_LIGAND_EDITOR_CANVAS, GtkWidget)

G_END_DECLS

namespace coot::ligand_editor_canvas::impl {

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

/// Used for widget's struct as a base class.
/// Useful for exposing inner state to the active tool.
struct WidgetCoreData {
    /// molecules on the screen
    std::unique_ptr<std::vector<CanvasMolecule>> molecules;
    /// molecules (RDKit)
    std::unique_ptr<std::vector<std::shared_ptr<RDKit::RWMol>>> rdkit_molecules;
};

/// This is the private struct for GObject
struct CootLigandEditorCanvasPriv : CootLigandEditorCanvasPrivBase, impl::WidgetCoreData {    

    std::unique_ptr<ActiveTool> active_tool;
};



}

#endif //#define COOT_LIGAND_EDITOR_CANVAS_CORE_HPP
