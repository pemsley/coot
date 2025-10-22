/* layla/ligand_editor_canvas.cpp
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

#include "ligand_editor_canvas.hpp"
#include "ligand_editor_canvas/core.hpp"
#include "ligand_editor_canvas/model.hpp"
#include "ligand_editor_canvas/tools.hpp"
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <exception>
#include <iterator>
#include <utility>
#include <algorithm>
#include <vector>
#include <memory>

#include "../utils/base64-encode-decode.hh"

// Prevents preprocessor substitution of `VERSION` in `MolPickler.h`
#ifndef RD_MOLPICKLE_H

#ifdef VERSION
#define __COOT_VERSION_VALUE VERSION
#undef VERSION
#endif

#include <rdkit/GraphMol/MolPickler.h>

#ifdef __COOT_VERSION_VALUE
#define VERSION __COOT_VERSION_VALUE
#undef __COOT_VERSION_VALUE
#endif

#endif //RD_MOLPICKLE_H

#ifndef __EMSCRIPTEN__
#include <cairo.h>
#include <pango/pango-font.h>
#include <pango/pangocairo.h>
#else
// Lhasa-specific includes
#endif

using namespace coot::ligand_editor_canvas;

// This can be safely omitted in Lhasa
#ifndef __EMSCRIPTEN__

/// Because of GObject's amazing macro system, 
/// I can't use a "typedef" to denote that
/// "_CootLigandEditorCanvas" is the same type as "coot::ligand_editor_canvas::impl::CootLigandEditorCanvasPriv".
///
/// A "using" statement doesn't compile too. I don't think there's anything I might be doing wrong.
///
/// An empty struct with inheritance should compile to the exactly same thing though.
struct _CootLigandEditorCanvas:  coot::ligand_editor_canvas::impl::CootLigandEditorCanvasPriv {
    friend void coot_ligand_editor_canvas_init_impl(CootLigandEditorCanvas* self);
    friend void coot_ligand_editor_canvas_dispose_impl(CootLigandEditorCanvas* self);
};
#endif

// Constructor

#ifndef __EMSCRIPTEN__
void coot_ligand_editor_canvas_init_impl(CootLigandEditorCanvas* self) {
#else // Lhasa code
CootLigandEditorCanvas::CootLigandEditorCanvas() noexcept {
    auto* self = this;
#endif
    g_debug("Instantiating CootLigandEditorCanvas.");
    self->active_tool = std::make_unique<ActiveTool>();
    self->active_tool->set_core_widget_data(static_cast<impl::CootLigandEditorCanvasPriv*>(self));
    self->molecules = std::make_unique<std::vector<std::optional<CanvasMolecule>>>();
    self->rdkit_molecules = std::make_unique<std::vector<std::optional<std::shared_ptr<RDKit::RWMol>>>>();
    self->currently_created_bond = std::nullopt;
    self->state_stack = std::make_unique<impl::WidgetCoreData::StateStack>();
    self->display_mode = DisplayMode::Standard;
    self->scale = 1.0;
    self->viewport_origin_offset = std::make_pair(0,0);
    self->allow_invalid_molecules = false;
    self->state_stack_pos = -1;
}

// Destructor

#ifndef __EMSCRIPTEN__
void coot_ligand_editor_canvas_dispose_impl(CootLigandEditorCanvas* self) {
#else // Lhasa code
CootLigandEditorCanvas::~CootLigandEditorCanvas() noexcept {
    auto* self = this;
#endif
    g_debug("De-instantiating CootLigandEditorCanvas.");
    self->molecules.reset(nullptr);
    self->active_tool.reset(nullptr);
    self->rdkit_molecules.reset(nullptr);
    self->state_stack.reset(nullptr);
}


#ifndef __EMSCRIPTEN__
G_BEGIN_DECLS

G_DEFINE_TYPE(CootLigandEditorCanvas, coot_ligand_editor_canvas, GTK_TYPE_WIDGET)

// not sure what this is for or whether it is going to be needed at all
// struct _CootLigandEditorCanvasClass {
//     GObjectClass parent_class;
// };

#endif

#ifndef __EMSCRIPTEN__
void coot_ligand_editor_canvas_snapshot (GtkWidget *widget, GtkSnapshot *snapshot)
{
    CootLigandEditorCanvas* self = COOT_COOT_LIGAND_EDITOR_CANVAS(widget);

    float w = (float) gtk_widget_get_size(widget,GTK_ORIENTATION_HORIZONTAL);
    float h = (float) gtk_widget_get_size(widget,GTK_ORIENTATION_VERTICAL);
    const graphene_rect_t background_rect = graphene_rect_t{{0,0},{w,h}};
    const GdkRGBA background_color = GdkRGBA{1.f,1.f,1.f,1.f};
    gtk_snapshot_append_color(snapshot, &background_color, &background_rect);
    PangoLayout* pango_layout = pango_layout_new(gtk_widget_get_pango_context(widget));
    cairo_t *cr = gtk_snapshot_append_cairo(snapshot, &background_rect);
    impl::Renderer ren(cr,pango_layout);
    self->render(ren);
}
#else
// Rendering in Lhasa is done via manually instantiating the renderer 
// and calling the 'render()' method
#endif

#ifndef __EMSCRIPTEN__
void coot_ligand_editor_canvas_measure(GtkWidget *widget, GtkOrientation orientation, int for_size, int *minimum_size, int *natural_size, int *minimum_baseline, int *natural_baseline)
{
    CootLigandEditorCanvas* self = COOT_COOT_LIGAND_EDITOR_CANVAS(widget);
#else
CootLigandEditorCanvas::SizingInfo CootLigandEditorCanvas::measure(CootLigandEditorCanvas::MeasurementDirection orientation) noexcept {
    auto* self = this;
    SizingInfo ret;
    int* minimum_size = &ret.requested_size;
    int* natural_size = &ret.requested_size;
#endif
    graphene_rect_t bounding_rect_for_all;
    graphene_rect_init(&bounding_rect_for_all, 0, 0, 0, 0);

    for(const auto& a: *self->molecules) {
        if(a.has_value()) {
            auto bounding_rect = a->get_on_screen_bounding_rect(self->viewport_origin_offset, self->scale);
            graphene_rect_union(&bounding_rect_for_all, &bounding_rect, &bounding_rect_for_all);
        }
    }
    switch (orientation) {
        #ifndef __EMSCRIPTEN__
        case GTK_ORIENTATION_HORIZONTAL:{
        #else
        case MeasurementDirection::HORIZONTAL:{
        #endif
            // For now:
            *natural_size = bounding_rect_for_all.size.width;
            *minimum_size = bounding_rect_for_all.size.width;
            break;
        }
        #ifndef __EMSCRIPTEN__
        case GTK_ORIENTATION_VERTICAL:{
        #else
        case MeasurementDirection::VERTICAL:{
        #endif
            // For now:
            *natural_size = bounding_rect_for_all.size.height;
            *minimum_size = bounding_rect_for_all.size.height;
            break;
        }
        default:
            break;
    }
    #ifdef __EMSCRIPTEN__
    // if(ret.requested_size == 0) {
    //     g_warning("FIXME: Overriding zeroed 'requested_size' with 600.");
    //     ret.requested_size = 600;
    // }
    return ret;
    #endif
    std::pair<int, int> viewport_offset = {0, 0};
    if(bounding_rect_for_all.origin.x < 0) {
        viewport_offset.first = bounding_rect_for_all.origin.x - impl::WidgetCoreData::VIEWPORT_OFFSET_PIXEL_MARGIN;
    }
    if(bounding_rect_for_all.origin.y < 0) {
        viewport_offset.second = bounding_rect_for_all.origin.y - impl::WidgetCoreData::VIEWPORT_OFFSET_PIXEL_MARGIN;
    }
    // It is controversial to make any mutations in a 'measure' method,
    // but this is the easiest way to implement panning via offsetting the viewport
    // without reworking the molecule mutation process.
    self->viewport_origin_offset = viewport_offset;
}

#ifndef __EMSCRIPTEN__
static void on_hover (
  GtkEventControllerMotion* controller,
  gdouble x,
  gdouble y,
  gpointer user_data
) {
    GdkEvent* event = gtk_event_controller_get_current_event(GTK_EVENT_CONTROLLER(controller));
    GdkModifierType modifiers = gdk_event_get_modifier_state(event);
    bool alt_pressed = modifiers & GDK_ALT_MASK;
    bool shift_pressed = modifiers & GDK_SHIFT_MASK;
    bool control_pressed = modifiers & GDK_CONTROL_MASK;

    CootLigandEditorCanvas* self = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
#else
void CootLigandEditorCanvas::on_hover(double x, double y, bool alt_pressed, bool control_pressed) {
    auto* self = this;
#endif

    // Clear all highlights first
    for(auto& molecule: *self->molecules) {
        if(molecule.has_value())    {
            molecule->clear_highlights();
        }
    }

    if(self->active_tool->is_in_transform()) {
        self->active_tool->update_transform_cursor_pos((int)x, (int)y, alt_pressed);
        self->queue_redraw();
        return;
    }

    // Update position of the second atom when creating a bond
    if(self->currently_created_bond.has_value()) {
        auto& new_bond = *self->currently_created_bond;
        new_bond.second_atom_x = x;
        new_bond.second_atom_y = y;
        // and set highlight for the first atom
        if(self->active_tool->is_creating_bond()) {
            auto [molecule_idx, atom_idx] = *self->active_tool->get_molecule_idx_and_first_atom_of_new_bond();
            auto& target = *(*self->molecules)[molecule_idx];
            target.add_atom_highlight(atom_idx, CanvasMolecule::HighlightType::Edition);
        }
    } else {
        self->active_tool->on_hover(alt_pressed, control_pressed, x, y);
    }

    // Highlights and snapping
    auto maybe_something_clicked = self->resolve_click(x, y);
    if(maybe_something_clicked.has_value()) {
        auto [bond_or_atom,molecule_idx] = *maybe_something_clicked;
        auto& target = *(*self->molecules)[molecule_idx];
        if(std::holds_alternative<CanvasMolecule::Atom>(bond_or_atom)) {
            auto atom = std::get<CanvasMolecule::Atom>(std::move(bond_or_atom));
            g_debug("Hovering on atom %u (%s)", atom.idx,atom.symbol.c_str());
            target.add_atom_highlight(atom.idx, CanvasMolecule::HighlightType::Hover);

            // Snapping to the target atom
            // when creating a bond
            if(self->currently_created_bond.has_value()) {
                auto& new_bond = *self->currently_created_bond;
                auto coords = target.get_on_screen_coords(atom.x, atom.y, self->viewport_origin_offset, self->scale);
                new_bond.second_atom_x = coords.first;
                new_bond.second_atom_y = coords.second;
            }
        } else { // a bond
            auto bond = std::get<CanvasMolecule::Bond>(std::move(bond_or_atom));
            g_debug("Hovering on bond between atoms %u and %u", bond.first_atom_idx, bond.second_atom_idx);
            target.add_bond_highlight(bond.first_atom_idx, bond.second_atom_idx, CanvasMolecule::HighlightType::Hover);
        }
    }
    self->queue_redraw();
}

#ifndef __EMSCRIPTEN__
static gboolean on_scroll(GtkEventControllerScroll* zoom_controller, gdouble dx, gdouble dy, gpointer user_data) {
    GdkEvent* event = gtk_event_controller_get_current_event(GTK_EVENT_CONTROLLER(zoom_controller));
    GdkModifierType modifiers = gdk_event_get_modifier_state(event);
    bool control_pressed = modifiers & GDK_CONTROL_MASK;

    CootLigandEditorCanvas* self = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
#else
void CootLigandEditorCanvas::on_scroll(double dx, double dy, bool control_pressed) {
    auto* self = this;
#endif
    if (control_pressed) {
        #ifndef __EMSCRIPTEN__
        self->scale *= (1.f - dy / 20.f);
        #else
        g_info("DeltaY = %f", dy);
        self->scale *= (1.f - dy / 800.f);
        #endif
        _LIGAND_EDITOR_SIGNAL_EMIT_ARG(self, scale_changed_signal,self->scale);
        self->queue_redraw();
        self->queue_resize();
        #ifndef __EMSCRIPTEN__
        return TRUE;
        #endif
    }
    #ifndef __EMSCRIPTEN__
    return FALSE;
    #endif
}

#ifndef __EMSCRIPTEN__
static void
on_left_click_released(
  GtkGestureClick* gesture_click,
  gint n_press,
  gdouble x,
  gdouble y,
  gpointer user_data
) {
    CootLigandEditorCanvas* self = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
    GdkEvent* event = gtk_event_controller_get_current_event(GTK_EVENT_CONTROLLER(gesture_click));
    GdkModifierType modifiers = gdk_event_get_modifier_state(event);
    bool alt_pressed = modifiers & GDK_ALT_MASK;
    bool shift_pressed = modifiers & GDK_SHIFT_MASK;
    bool control_pressed = modifiers & GDK_CONTROL_MASK;
#else
void CootLigandEditorCanvas::on_left_click_released(double x, double y, bool alt_pressed, bool control_pressed, bool shift_pressed) {
    auto* self = this;
#endif

    if(self->active_tool->is_in_transform()) {
        self->active_tool->end_transform(alt_pressed);
        return;
    }

    // `currently_created_bond` gets cleared here when appropriate
    self->active_tool->on_release(alt_pressed, control_pressed, x, y, false);
}

#ifndef __EMSCRIPTEN__
static void on_left_click(
  GtkGestureClick* gesture_click,
  gint n_press,
  gdouble x,
  gdouble y,
  gpointer user_data
) {
    CootLigandEditorCanvas* self = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
    GdkEvent* event = gtk_event_controller_get_current_event(GTK_EVENT_CONTROLLER(gesture_click));
    GdkModifierType modifiers = gdk_event_get_modifier_state(event);
    bool alt_pressed = modifiers & GDK_ALT_MASK;
    bool shift_pressed = modifiers & GDK_SHIFT_MASK;
    bool control_pressed = modifiers & GDK_CONTROL_MASK;
#else
void CootLigandEditorCanvas::on_left_click(double x, double y, bool alt_pressed, bool control_pressed, bool shift_pressed) {
    auto* self = this;
#endif

    if(!control_pressed) {
        if(shift_pressed) {
            self->active_tool->begin_transform(x, y, TransformManager::Mode::Rotation);
            return;
        } else if(alt_pressed) {
            self->active_tool->begin_transform(x, y, TransformManager::Mode::Translation);
            return;
        }
    }

    self->active_tool->on_click(alt_pressed, control_pressed, x, y, false);

    if(self->active_tool->is_creating_bond()) {
        CurrentlyCreatedBond new_bond;
        auto [mol_idx, atom_idx] = self->active_tool->get_molecule_idx_and_first_atom_of_new_bond().value();
        auto coords = *self->molecules->at(mol_idx)->get_on_screen_coords_of_atom(atom_idx, self->viewport_origin_offset, self->scale);
        new_bond.first_atom_x = coords.first;
        new_bond.first_atom_y = coords.second;
        new_bond.second_atom_x = coords.first;
        new_bond.second_atom_y = coords.second;
        self->currently_created_bond = new_bond;
    }
    //gtk_gesture_set_state(GTK_GESTURE(gesture_click),GTK_EVENT_SEQUENCE_CLAIMED);
}

#ifndef __EMSCRIPTEN__
static void
on_right_click_released(
  GtkGestureClick* gesture_click,
  gint n_press,
  gdouble x,
  gdouble y,
  gpointer user_data
) {
    CootLigandEditorCanvas* self = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
    GdkEvent* event = gtk_event_controller_get_current_event(GTK_EVENT_CONTROLLER(gesture_click));
    GdkModifierType modifiers = gdk_event_get_modifier_state(event);
    bool alt_pressed = modifiers & GDK_ALT_MASK;
    bool shift_pressed = modifiers & GDK_SHIFT_MASK;
    bool control_pressed = modifiers & GDK_CONTROL_MASK;
#else
void CootLigandEditorCanvas::on_right_click_released(double x, double y, bool alt_pressed, bool control_pressed, bool shift_pressed) {
    auto* self = this;
#endif

    self->active_tool->on_release(alt_pressed, control_pressed, x, y, true);
}

#ifndef __EMSCRIPTEN__
static void on_right_click(
  GtkGestureClick* gesture_click,
  gint n_press,
  gdouble x,
  gdouble y,
  gpointer user_data
) {
    CootLigandEditorCanvas* self = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
    GdkEvent* event = gtk_event_controller_get_current_event(GTK_EVENT_CONTROLLER(gesture_click));
    GdkModifierType modifiers = gdk_event_get_modifier_state(event);
    bool alt_pressed = modifiers & GDK_ALT_MASK;
    bool shift_pressed = modifiers & GDK_SHIFT_MASK;
    bool control_pressed = modifiers & GDK_CONTROL_MASK;
#else
void CootLigandEditorCanvas::on_right_click(double x, double y, bool alt_pressed, bool control_pressed, bool shift_pressed) {
    auto* self = this;
#endif

    self->active_tool->on_click(alt_pressed, control_pressed, x, y, true);
}


#ifndef __EMSCRIPTEN__
static void coot_ligand_editor_canvas_init(CootLigandEditorCanvas* self) {
    // This is the primary constructor
    
    // GObject doesn't run C++ constructors upon allocation
    // so we take care of this ourselves
    coot_ligand_editor_canvas_init_impl(self);
    GtkGesture* left_click_controller = gtk_gesture_click_new();
    GtkGesture* right_click_controller = gtk_gesture_click_new();
    GtkEventController* hover_controller = gtk_event_controller_motion_new();
    GtkEventController* zoom_controller = gtk_event_controller_scroll_new(GTK_EVENT_CONTROLLER_SCROLL_VERTICAL);

    // left mouse button
    gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(left_click_controller),GDK_BUTTON_PRIMARY);
    g_signal_connect(left_click_controller,"pressed",G_CALLBACK(on_left_click),self);
    g_signal_connect(left_click_controller,"released",G_CALLBACK(on_left_click_released),self);

    //right mouse button
    gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(right_click_controller),GDK_BUTTON_SECONDARY);
    g_signal_connect(right_click_controller,"pressed",G_CALLBACK(on_right_click),self);
    g_signal_connect(right_click_controller,"released",G_CALLBACK(on_right_click_released),self);

    g_signal_connect(hover_controller,"motion",G_CALLBACK(on_hover),self);

    g_signal_connect(zoom_controller, "scroll",G_CALLBACK(on_scroll), self);

    gtk_widget_add_controller(GTK_WIDGET(self),GTK_EVENT_CONTROLLER(left_click_controller));
    gtk_widget_add_controller(GTK_WIDGET(self),GTK_EVENT_CONTROLLER(right_click_controller));
    gtk_widget_add_controller(GTK_WIDGET(self),GTK_EVENT_CONTROLLER(hover_controller));
    gtk_widget_add_controller(GTK_WIDGET(self), GTK_EVENT_CONTROLLER(zoom_controller));
}
#endif


#ifndef __EMSCRIPTEN__
static void coot_ligand_editor_canvas_dispose(GObject* _self) {
    CootLigandEditorCanvas* self = COOT_COOT_LIGAND_EDITOR_CANVAS(_self);
    // GObject doesn't run C++ destructors
    // so we take care of this ourselves
    coot_ligand_editor_canvas_dispose_impl(self);
    G_OBJECT_CLASS(coot_ligand_editor_canvas_parent_class)->dispose(_self);
}
#endif

#ifndef __EMSCRIPTEN__
static void coot_ligand_editor_canvas_class_init(CootLigandEditorCanvasClass* klass) {
    // I think that this is a GObject class constructor that sets up the GObject class at runtime.
    impl::status_updated_signal = g_signal_new("status-updated",
        G_TYPE_FROM_CLASS (klass),
        (GSignalFlags) (G_SIGNAL_RUN_LAST | G_SIGNAL_NO_RECURSE | G_SIGNAL_NO_HOOKS),
        0 /* class offset.Subclass cannot override the class handler (default handler). */,
        NULL /* accumulator */,
        NULL /* accumulator data */,
        NULL /* C marshaller. g_cclosure_marshal_generic() will be used */,
        G_TYPE_NONE /* return_type */,
        1     /* n_params */,
        G_TYPE_STRING
    );
    impl::scale_changed_signal = g_signal_new("scale-changed",
        G_TYPE_FROM_CLASS (klass),
        (GSignalFlags) (G_SIGNAL_RUN_LAST | G_SIGNAL_NO_RECURSE | G_SIGNAL_NO_HOOKS),
        0 /* class offset.Subclass cannot override the class handler (default handler). */,
        NULL /* accumulator */,
        NULL /* accumulator data */,
        NULL /* C marshaller. g_cclosure_marshal_generic() will be used */,
        G_TYPE_NONE /* return_type */,
        1     /* n_params */,
        G_TYPE_FLOAT
    );
    impl::molecule_deleted_signal = g_signal_new("molecule-deleted",
        G_TYPE_FROM_CLASS (klass),
        (GSignalFlags) (G_SIGNAL_RUN_LAST | G_SIGNAL_NO_RECURSE | G_SIGNAL_NO_HOOKS),
        0 /* class offset.Subclass cannot override the class handler (default handler). */,
        NULL /* accumulator */,
        NULL /* accumulator data */,
        NULL /* C marshaller. g_cclosure_marshal_generic() will be used */,
        G_TYPE_NONE /* return_type */,
        1     /* n_params */,
        G_TYPE_INT
    );
    impl::smiles_changed_signal = g_signal_new("smiles-changed",
        G_TYPE_FROM_CLASS (klass),
        (GSignalFlags) (G_SIGNAL_RUN_LAST | G_SIGNAL_NO_RECURSE | G_SIGNAL_NO_HOOKS),
        0 /* class offset.Subclass cannot override the class handler (default handler). */,
        NULL /* accumulator */,
        NULL /* accumulator data */,
        NULL /* C marshaller. g_cclosure_marshal_generic() will be used */,
        G_TYPE_NONE /* return_type */,
        0     /* n_params */
    );
    impl::qed_info_updated_signal = g_signal_new("qed-info-updated",
        G_TYPE_FROM_CLASS (klass),
        (GSignalFlags) (G_SIGNAL_RUN_LAST | G_SIGNAL_NO_RECURSE | G_SIGNAL_NO_HOOKS),
        0 /* class offset.Subclass cannot override the class handler (default handler). */,
        NULL /* accumulator */,
        NULL /* accumulator data */,
        NULL /* C marshaller. g_cclosure_marshal_generic() will be used */,
        G_TYPE_NONE /* return_type */,
        2     /* n_params */,
        G_TYPE_INT,
        G_TYPE_POINTER
    );
    GTK_WIDGET_CLASS(klass)->snapshot = coot_ligand_editor_canvas_snapshot;
    GTK_WIDGET_CLASS(klass)->measure = coot_ligand_editor_canvas_measure;
    G_OBJECT_CLASS(klass)->dispose = coot_ligand_editor_canvas_dispose;
    
}
#endif

#ifndef __EMSCRIPTEN__
CootLigandEditorCanvas* 
coot_ligand_editor_canvas_new()
{
    return COOT_COOT_LIGAND_EDITOR_CANVAS(g_object_new (COOT_LIGAND_EDITOR_CANVAS_TYPE, NULL));
}
#else // Lhasa-specific code
CootLigandEditorCanvas* 
coot_ligand_editor_canvas_new()
{
    return new CootLigandEditorCanvas();
}
#endif

#ifndef __EMSCRIPTEN__
G_END_DECLS
#endif

void coot_ligand_editor_canvas_set_scale(CootLigandEditorCanvas* self, float display_scale) noexcept {
    self->scale = display_scale;
    _LIGAND_EDITOR_SIGNAL_EMIT_ARG(self, scale_changed_signal,self->scale);
    self->queue_redraw();
    self->queue_resize();
}

float coot_ligand_editor_canvas_get_scale(CootLigandEditorCanvas* self) noexcept {
    return self->scale;
}

void coot_ligand_editor_canvas_set_active_tool(CootLigandEditorCanvas* self, std::unique_ptr<ActiveTool>&& active_tool) {
    self->active_tool = std::move(active_tool);
    self->active_tool->set_core_widget_data(static_cast<impl::CootLigandEditorCanvasPriv*>(self));
    self->active_tool->on_load();
}

int coot_ligand_editor_canvas_append_molecule(CootLigandEditorCanvas* self, std::shared_ptr<RDKit::RWMol> rdkit_mol) noexcept {
    if(rdkit_mol->getNumAtoms() == 0) {
        self->update_status("Attempted to add an empty molecule!");
        g_warning("Attempted to add an empty molecule!");
        return -1;
    }
    try {
        g_debug("Appending new molecule to the widget...");
        // Might throw if the constructor fails.
        self->begin_edition();
        self->molecules->push_back(CanvasMolecule(rdkit_mol, self->allow_invalid_molecules));
        #ifndef __EMSCRIPTEN__
        self->molecules->back()->apply_canvas_translation(
            gtk_widget_get_size(GTK_WIDGET(self), GTK_ORIENTATION_HORIZONTAL) / 2.0, 
            gtk_widget_get_size(GTK_WIDGET(self), GTK_ORIENTATION_VERTICAL) / 2.0,
            self->scale
        );
        #else
        self->molecules->back()->apply_canvas_translation(
            self->measure(CootLigandEditorCanvas::MeasurementDirection::HORIZONTAL).requested_size / 2.0, 
            self->measure(CootLigandEditorCanvas::MeasurementDirection::VERTICAL).requested_size / 2.0,
            self->scale
        );
        #endif
        self->rdkit_molecules->push_back(std::move(rdkit_mol));
        self->finalize_edition();
        // Already called by finalize_edition()
        //self->queue_redraw();
        self->update_status("Molecule inserted.");
        return self->rdkit_molecules->size() - 1;
    }catch(const std::exception& e) {
        std::string msg = "2D representation could not be created: ";
        msg += e.what();
        msg += ". New molecule could not be added.";
        g_warning("coot_ligand_editor_canvas_append_molecule: %s\n%s",
            msg.c_str(), 
            self->allow_invalid_molecules ? "Invalid molecules are allowed." : "Invalid molecules are not allowed."
        );
        self->update_status(msg.c_str());
        self->rollback_current_edition();
        return -1;
    }
}

void coot_ligand_editor_canvas_update_molecule_from_smiles(CootLigandEditorCanvas* self, unsigned int molecule_idx, const char* smiles) {
    if(molecule_idx >= self->rdkit_molecules->size()) {
        return;
    }
    auto& target_mol_opt = (*self->rdkit_molecules)[molecule_idx];
    if(!target_mol_opt.has_value()) {
        return;
    }
    try{
        const auto* mol_ptr = RDKit::SmilesToMol(smiles, 0, !self->allow_invalid_molecules);
        if(mol_ptr) {
            self->begin_edition();
            *target_mol_opt->get() = std::move(*mol_ptr);
            auto& widget_mol = (*self->molecules)[molecule_idx];
            widget_mol->clear_cached_atom_coordinate_map();
            widget_mol->lower_from_rdkit(!self->allow_invalid_molecules);
            self->finalize_edition();
            self->update_status("Molecule updated from SMILES.");
            delete mol_ptr;
        }
    }catch(const std::exception& e) {
        std::string msg = "2D representation could not be created: ";
        msg += e.what();
        msg += ". Molecule could not be updated.";
        g_warning(
            "coot_ligand_editor_canvas_update_molecule_from_smiles: %s\n%s",msg.c_str(),
            self->allow_invalid_molecules ? "Invalid molecules are allowed." : "Invalid molecules are not allowed."
        );
        self->update_status(msg.c_str());
        self->rollback_current_edition();
    }
}

void coot_ligand_editor_canvas_undo_edition(CootLigandEditorCanvas* self) noexcept {
    self->undo_edition();
    self->queue_redraw();
    self->emit_mutation_signals();
}

void coot_ligand_editor_canvas_redo_edition(CootLigandEditorCanvas* self) noexcept {
    self->redo_edition();
    self->queue_redraw();
    self->emit_mutation_signals();
}

const RDKit::ROMol* coot_ligand_editor_canvas_get_rdkit_molecule(CootLigandEditorCanvas* self, unsigned int index) noexcept {
    if(self->rdkit_molecules->size() > index) {
        const auto& vec = *self->rdkit_molecules.get();
        const auto& ret_opt = vec[index];
        if(ret_opt.has_value()) {
            return ret_opt->get();
        }
    }
    return nullptr;
}

unsigned int coot_ligand_editor_canvas_get_molecule_count(CootLigandEditorCanvas* self) noexcept {
    return self->get_molecule_count_impl();
}

unsigned int coot_ligand_editor_canvas_get_idx_of_first_molecule(CootLigandEditorCanvas* self) noexcept {
    return self->get_first_molecule_idx();
}

unsigned int coot_ligand_editor_canvas_get_max_molecule_idx(CootLigandEditorCanvas* self) noexcept {
    return self->molecules->size() - 1;
}

void coot_ligand_editor_canvas_set_allow_invalid_molecules(CootLigandEditorCanvas* self, bool value) noexcept {
    g_debug("Invalid molecules allowed set to %s", value ? "true" : "false");
    self->allow_invalid_molecules = value;
}

bool coot_ligand_editor_canvas_get_allow_invalid_molecules(CootLigandEditorCanvas* self) noexcept {
    return self->allow_invalid_molecules;
}

DisplayMode coot_ligand_editor_canvas_get_display_mode(CootLigandEditorCanvas* self) noexcept {
    return self->display_mode;
}

void coot_ligand_editor_canvas_set_display_mode(CootLigandEditorCanvas* self, DisplayMode value) noexcept {
    self->display_mode = value;
    self->queue_redraw();
}

coot::ligand_editor_canvas::SmilesMap coot_ligand_editor_canvas_get_smiles(CootLigandEditorCanvas* self) noexcept {
    return self->build_smiles();
}

coot::ligand_editor_canvas::InchiKeyMap coot_ligand_editor_canvas_get_inchi_keys(CootLigandEditorCanvas* self) noexcept {
    return self->build_inchi_keys();
}

std::string coot_ligand_editor_canvas_get_smiles_for_molecule(CootLigandEditorCanvas* self, unsigned int molecule_idx) noexcept {
    if(molecule_idx < self->rdkit_molecules->size()) {
        const auto& mol_opt = (*self->rdkit_molecules)[molecule_idx];
        if(mol_opt.has_value()) {
            return RDKit::MolToSmiles(*mol_opt->get());
        }
    }
    return "";
}

std::string coot_ligand_editor_canvas_get_inchi_key_for_molecule(CootLigandEditorCanvas* self, unsigned int molecule_idx) noexcept {
    if(molecule_idx < self->rdkit_molecules->size()) {
        const auto& mol_opt = (*self->rdkit_molecules)[molecule_idx];
        if(mol_opt.has_value()) {
            #ifdef RDK_BUILD_INCHI_SUPPORT
            return RDKit::MolToInchiKey(*mol_opt->get());
            #else
            return "";
            #warning Your version of RDKit was built without InChI support. Molecule InChI key lookup will not be available.
            #endif
            
        }
    }
    return "";
}

std::string coot_ligand_editor_canvas_get_pickled_molecule(CootLigandEditorCanvas* self, unsigned int molecule_idx) noexcept {
    std::string ret = "";
    if(molecule_idx < self->rdkit_molecules->size()) {
        const auto& mol_opt = (*self->rdkit_molecules)[molecule_idx];
        if(mol_opt.has_value()) {
            RDKit::MolPickler::pickleMol(mol_opt->get(), ret);
        }
    }
    return ret;
}

std::string coot_ligand_editor_canvas_get_pickled_molecule_base64(CootLigandEditorCanvas* self, unsigned int molecule_idx) noexcept {
    std::string raw = coot_ligand_editor_canvas_get_pickled_molecule(self, molecule_idx);
    return moorhen_base64::base64_encode((const unsigned char*) raw.c_str(), raw.size());
}

#ifndef __EMSCRIPTEN__
void coot_ligand_editor_canvas_draw_on_cairo_surface(CootLigandEditorCanvas* self, cairo_t* cr) noexcept {
    PangoLayout* pango_layout = pango_cairo_create_layout(cr);
    PangoFontDescription* font_description = pango_font_description_new ();
    pango_font_description_set_family(font_description, "sans");

    pango_layout_set_font_description (pango_layout, font_description);
    impl::Renderer ren(cr, pango_layout);
    self->render(ren);

    pango_font_description_free(font_description);
}
#endif

void coot_ligand_editor_canvas_clear_molecules(CootLigandEditorCanvas* self) noexcept {
    self->begin_edition();
    unsigned int idx = 0;
    for(const auto& mol_opt: *self->molecules) {
        if(mol_opt.has_value()) {
            self->delete_molecule_with_idx(idx, false);
        }
        idx++;
    }
    self->rdkit_molecules->clear();
    self->molecules->clear();
    self->finalize_edition();
    self->update_status("Molecules cleared.");
    // Already called in finalize_edition()
    // self->queue_redraw();
    // self->emit_mutation_signals();
}
