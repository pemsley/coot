use crate::model::*;
use crate::widget::*;

use glib::gobject_ffi::*;
use glib::translate::*;
// #[no_mangle]
// pub extern "C" fn test(k: ActiveTool) {

// }

// The C Glib type
type CootLigandEditorCanvas =
    <imp::LigandEditorCanvas as glib::subclass::types::ObjectSubclass>::Instance;

#[no_mangle]
pub extern "C" fn coot_ligand_editor_canvas_new() -> *mut CootLigandEditorCanvas {
    let canvas = LigandEditorCanvas::new();
    let r = canvas.to_glib_full();
    r
}
