#include "ligand_editor_canvas.hpp"

struct _CootLigandEditorCanvas {
    GtkWidget parent;

};

G_BEGIN_DECLS

G_DEFINE_TYPE(CootLigandEditorCanvas, coot_ligand_editor_canvas, GTK_TYPE_WIDGET)

// not sure what this is for or whether it is going to be needed at all
// struct _CootLigandEditorCanvasClass {
//     GObjectClass parent_class;
// };

void coot_ligand_editor_canvas_snapshot (GtkWidget *widget, GtkSnapshot *snapshot)
{
    CootLigandEditorCanvas* self = COOT_COOT_LIGAND_EDITOR_CANVAS(widget);
   
}

void coot_ligand_editor_canvas_measure
    (GtkWidget      *widget,
    GtkOrientation  orientation,
    int             for_size,
    int            *minimum_size,
    int            *natural_size,
    int            *minimum_baseline,
    int            *natural_baseline)
{
    CootLigandEditorCanvas* self = COOT_COOT_LIGAND_EDITOR_CANVAS(widget);

    switch (orientation)
    {
    case GTK_ORIENTATION_HORIZONTAL:{
    
        break;
    }
    case GTK_ORIENTATION_VERTICAL:{
    
        break;
    }
    default:
        break;
    }

}

static void on_left_click (
  GtkGestureClick* gesture_click,
  gint n_press,
  gdouble x,
  gdouble y,
  gpointer user_data
) {
    CootLigandEditorCanvas* self = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
    //gtk_gesture_set_state(GTK_GESTURE(gesture_click),GTK_EVENT_SEQUENCE_CLAIMED);
    
    //gtk_gesture_set_state(GTK_GESTURE(gesture_click),GTK_EVENT_SEQUENCE_NONE);
    
}

static void coot_ligand_editor_canvas_init(CootLigandEditorCanvas* self) {
    // This is the primary constructor

    GtkGesture* click_controller = gtk_gesture_click_new();
    // GtkEventController* hover_controller = gtk_event_controller_motion_new();

    // left mouse button
    gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(click_controller),GDK_BUTTON_PRIMARY);
    g_signal_connect(click_controller,"pressed",G_CALLBACK(on_left_click),self);

    // g_signal_connect(hover_controller,"motion",G_CALLBACK(on_hover),self);

    gtk_widget_add_controller(GTK_WIDGET(self),GTK_EVENT_CONTROLLER(click_controller));
    // gtk_widget_add_controller(GTK_WIDGET(self),GTK_EVENT_CONTROLLER(hover_controller));
}

static void coot_ligand_editor_canvas_dispose(GObject* _self) {

    CootLigandEditorCanvas* self = COOT_COOT_LIGAND_EDITOR_CANVAS(_self);

    G_OBJECT_CLASS(coot_ligand_editor_canvas_parent_class)->dispose(_self);
}

static void coot_ligand_editor_canvas_class_init(CootLigandEditorCanvasClass* klass) {
    // I think that this is a GObject class constructor that sets up the GObject class at runtime.
    GTK_WIDGET_CLASS(klass)->snapshot = coot_ligand_editor_canvas_snapshot;
    GTK_WIDGET_CLASS(klass)->measure = coot_ligand_editor_canvas_measure;
    G_OBJECT_CLASS(klass)->dispose = coot_ligand_editor_canvas_dispose;
    
}

CootLigandEditorCanvas* 
coot_ligand_editor_canvas_new()
{
    return COOT_COOT_LIGAND_EDITOR_CANVAS(g_object_new (COOT_LIGAND_EDITOR_CANVAS_TYPE, NULL));
}

G_END_DECLS