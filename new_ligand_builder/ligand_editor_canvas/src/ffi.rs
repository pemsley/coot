use crate::model::*;
use crate::widget::*;

use glib::gobject_ffi;
use glib::translate::*;
use glib::types::StaticType;
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

// Not sure if this is needed
#[no_mangle]
pub unsafe extern "C" fn coot_ligand_editor_canvas_free(zelf: *mut CootLigandEditorCanvas) {
    gobject_ffi::g_object_unref(zelf as *mut _);
}

// consider exporting
// is_TYPE

#[no_mangle]
pub unsafe extern "C" fn coot_ligand_editor_canvas_get_type() -> glib::ffi::GType {
    LigandEditorCanvas::static_type().into_glib()
}
