//use glib;

use glib::Object;

glib::wrapper! {
    pub struct LigandEditorCanvas(ObjectSubclass<imp::LigandEditorCanvas>)
        @extends gtk::Widget,
        @implements gtk::Accessible, gtk::Actionable, gtk::Buildable, gtk::ConstraintTarget;
}

impl LigandEditorCanvas {
    pub fn new() -> Self {
        Object::builder().build()
    }

    // pub fn with_label(label: &str) -> Self {
    //     Object::builder().property("label", label).build()
    // }
}
pub mod imp {
    use gtk::subclass::prelude::*;

    // Object holding the state
    #[derive(Default)]
    pub struct LigandEditorCanvas {
        
    }

    // The central trait for subclassing a GObject
    #[glib::object_subclass]
    impl ObjectSubclass for LigandEditorCanvas {
        const NAME: &'static str = "CootLigandEditorCanvas";
        type Type = super::LigandEditorCanvas;
        type ParentType = gtk::Widget;
    }

    // Trait shared by all GObjects
    impl ObjectImpl for LigandEditorCanvas {}

    // Trait shared by all widgets
    impl WidgetImpl for LigandEditorCanvas {}

}