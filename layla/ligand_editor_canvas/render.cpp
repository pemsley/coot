/* layla/ligand_editor_canvas/render.cpp
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

#include "render.hpp"
#include "model.hpp"
#ifdef __EMSCRIPTEN__
#include "../../lhasa/glog_replacement.hpp"
#else
#include <sstream>
#endif

using namespace coot::ligand_editor_canvas::impl;
using Atom = coot::ligand_editor_canvas::CanvasMolecule::Atom;
using BondGeometry = coot::ligand_editor_canvas::CanvasMolecule::BondGeometry;
using BondType = coot::ligand_editor_canvas::CanvasMolecule::BondType;
using DoubleBondDrawingDirection = coot::ligand_editor_canvas::CanvasMolecule::DoubleBondDrawingDirection;

const float MoleculeRenderContext::CENTERED_DOUBLE_BOND_LINE_SEPARATION = 0.2f;
// 10 degrees
const float MoleculeRenderContext::GEOMETRY_BOND_SPREAD_ANGLE = M_PI/18.f;
const float MoleculeRenderContext::WAVY_BOND_ARC_LENGTH = 0.25f;
const float MoleculeRenderContext::GEOMETRY_BOND_DASH_SEPARATION = 0.2f;

#ifndef __EMSCRIPTEN__
Renderer::Renderer(cairo_t* cr, PangoLayout* pango_layout) {
    this->cr = cr;
    this->pango_layout = pango_layout;
}
#else // __EMSCRIPTEN__ defined
// Lhasa-specific includes/definitions
Renderer::Renderer(emscripten::val text_measurement_function, Renderer::TextMeasurementCache& cache) 
 :Renderer(text_measurement_function) {
    this->tm_cache = &cache;
}

Renderer::Renderer(emscripten::val text_measurement_function) {
    this->text_measurement_function = text_measurement_function;
    this->tm_cache = nullptr;
    this->position.x = 0.f;
    this->position.y = 0.f;
    this->style.line_width = 1.0f;
    this->style.color.r = 0.f;
    this->style.color.g = 0.f;
    this->style.color.b = 0.f;
    this->style.color.a = 1.f;
}

bool Renderer::DrawingCommand::is_path() const {
    return std::holds_alternative<Renderer::Path>(this->content);
}

bool Renderer::DrawingCommand::is_text() const {
    return std::holds_alternative<Renderer::Text>(this->content);
}

const Renderer::Path& Renderer::DrawingCommand::as_path() const {
    return std::get<Renderer::Path>(this->content);
}

Renderer::Path& Renderer::DrawingCommand::as_path_mut() {
    return std::get<Renderer::Path>(this->content);
}

const Renderer::Text& Renderer::DrawingCommand::as_text() const {
    return std::get<Renderer::Text>(this->content);
}

const Renderer::Line& Renderer::PathElement::as_line() const {
    return std::get<Renderer::Line>(this->content);
}

const Renderer::Arc& Renderer::PathElement::as_arc() const {
    return std::get<Renderer::Arc>(this->content);
}

bool Renderer::PathElement::is_line() const {
    return std::holds_alternative<Renderer::Line>(this->content);
}

bool Renderer::PathElement::is_arc() const {
    return std::holds_alternative<Renderer::Arc>(this->content);
}

std::vector<Renderer::DrawingCommand> Renderer::get_commands() const {
    return this->drawing_commands;
}

const std::vector<Renderer::PathElement>& Renderer::Path::get_elements() const {
    return this->elements;
}

Renderer::Path Renderer::create_new_path() const {
    Path ret;
    ret.has_fill = false;
    ret.initial_point = this->position;
    ret.has_stroke = false;
    ret.closed = false;
    return ret;
}

Renderer::Path& Renderer::top_path() {
    if(this->drawing_commands.empty()) {
        this->drawing_commands.push_back({this->create_new_path()});
    }
    auto& last_el = this->drawing_commands.back();
    if(last_el.is_path()) {
        const auto& pth = last_el.as_path();
        if(!pth.closed) {
            return last_el.as_path_mut();
        }
    }
    this->drawing_commands.push_back({this->create_new_path()});
    return this->drawing_commands.back().as_path_mut();
}

Renderer::Path* Renderer::top_path_if_exists() {
    if(this->drawing_commands.empty()) {
        return nullptr;
    }
    auto& last_el = this->drawing_commands.back();
    if(last_el.is_path()) {
        return &last_el.as_path_mut();
    } else {
        return nullptr;
    }
}

#endif // Emscripten defined

Renderer::TextStyle::TextStyle() {
    this->positioning = TextPositioning::Normal;
    // We want to set them to "" to denote unspecified values by defult
    // this->weight = "normal";
    // this->size = "medium";
    this->color.r = 0.f;
    this->color.g = 0.f;
    this->color.b = 0.f;
    this->color.a = 1.f;
    this->specifies_color = false;
}

Renderer::TextSpan::TextSpan() {
    this->specifies_style = false;
    this->content = std::string();
}

Renderer::TextSpan::TextSpan(const std::string& caption) {
    this->specifies_style = false;
    this->content = caption;
}

Renderer::TextSpan::TextSpan(const std::vector<TextSpan>& subspans) {
    this->specifies_style = false;
    this->content = subspans;
}

Renderer::TextSpan::TextSpan(const Newline& nl) {
    this->specifies_style = false;
    this->content = nl;
}

bool Renderer::TextSpan::has_subspans() const {
    return std::holds_alternative<std::vector<TextSpan>>(this->content);
}

bool Renderer::TextSpan::is_newline() const {
    return std::holds_alternative<Newline>(this->content);
}

std::string& Renderer::TextSpan::as_caption() {
    return std::get<std::string>(this->content);
}

std::vector<Renderer::TextSpan>& Renderer::TextSpan::as_subspans() {
    return std::get<std::vector<TextSpan>>(this->content);
}

const std::string& Renderer::TextSpan::as_caption() const {
    return std::get<std::string>(this->content);
}

const std::vector<Renderer::TextSpan>& Renderer::TextSpan::as_subspans() const {
    return std::get<std::vector<TextSpan>>(this->content);
}

void Renderer::move_to(double x, double y) {
    #ifndef __EMSCRIPTEN__
    cairo_move_to(cr, x, y);
    #else // __EMSCRIPTEN__ defined
    this->position.x = x;
    this->position.y = y;
    // Should it always do this?
    // Cairo docs say so.
    this->new_sub_path();
    #endif
}

void Renderer::line_to(double x, double y) {
    #ifndef __EMSCRIPTEN__
    cairo_line_to(cr, x, y);
    #else // __EMSCRIPTEN__ defined
    Line line;
    line.start = this->position;
    line.end.x = x;
    line.end.y = y;
    graphene_point_t final_pos = line.end;
    this->top_path().elements.push_back({line});
    this->position = final_pos;
    #endif
}

void Renderer::arc(double x, double y, double radius, double angle_one, double angle_two) {
    #ifndef __EMSCRIPTEN__
    cairo_arc(cr, x, y, radius, angle_one, angle_two);
    #else // __EMSCRIPTEN__ defined
    Arc arc;
    arc.origin.x = x;
    arc.origin.y = y;
    arc.radius = radius;
    arc.angle_one = angle_one;
    arc.angle_two = angle_two;
    this->top_path().elements.push_back({arc});
    #endif
}

void Renderer::fill() {
    #ifndef __EMSCRIPTEN__
    cairo_fill(cr);
    #else // __EMSCRIPTEN__ defined
    auto* path = this->top_path_if_exists();
    if(!path) {
        g_warning("fill() called without a path.");
        return;
    }
    if(path->elements.empty()) {
        g_warning("fill() called with an empty path.");
        //return;
    }
    // Fill closes all opened subpaths I guess?
    path->closed = true;
    // Should I do anything else with it?

    path->has_fill = true;
    path->fill_color = this->style.color;
    #endif
}

void Renderer::stroke() {
    #ifndef __EMSCRIPTEN__
    cairo_stroke(cr);
    #else // __EMSCRIPTEN__ defined
    this->stroke_preserve();
    this->new_path();
    #endif
}

void Renderer::stroke_preserve() {
    #ifndef __EMSCRIPTEN__
    cairo_stroke_preserve(cr);
    #else // __EMSCRIPTEN__ defined
    auto* path = this->top_path_if_exists();
    if(!path) {
        g_warning("stroke() called without a path.");
        return;
    }
    if(path->elements.empty()) {
        g_warning("stroke() called with an empty path.");
        //return;
    }
    path->has_stroke = true;
    path->stroke_style = this->style;
    #endif
}

void Renderer::new_path() {
    #ifndef __EMSCRIPTEN__
    cairo_new_path(cr);
    #else // __EMSCRIPTEN__ defined
    auto* path = this->top_path_if_exists();
    if(path) {
        // if(!path->closed) {
        //     // Reset the path
        //     *path = this->create_new_path();
        // } else {
            // Actually creates another new path
            this->close_path_inner();
        // }
    }
    // Is this really what the function ought to do?
    #endif
}

#ifdef __EMSCRIPTEN__
void Renderer::close_path_inner() {
    // 2. Close the sub-path
    auto* path = this->top_path_if_exists();
    if(path) {
        if(!path->closed) {
            path->closed = true;
        } else {
            // Technically, just creating new path, means that the previous one is done with.
            this->drawing_commands.push_back({this->create_new_path()});
        }
    }

    // 3. 
    // To quote the docs:
    // "The behavior of cairo_close_path() is distinct from simply calling cairo_line_to() 
    // with the equivalent coordinate in the case of stroking. When a closed sub-path is stroked, 
    // there are no caps on the ends of the sub-path. Instead, there is a line join connecting 
    // the final and initial segments of the sub-path."
    //
    // Now, I don't exactly now what that means for me.
    // 
}
#endif

void Renderer::close_path() {
    #ifndef __EMSCRIPTEN__
    cairo_close_path(cr);
    #else // __EMSCRIPTEN__ defined
    // // 1. Add a line to the beginning of the path
    auto* path = this->top_path_if_exists();
    if(path) {
        if(!path->closed) {
            this->line_to(path->initial_point.x, path->initial_point.y);
        }
        this->close_path_inner();
    }
    #endif
}

void Renderer::new_sub_path() {
    #ifndef __EMSCRIPTEN__
    cairo_new_sub_path(cr);
    #else // __EMSCRIPTEN__ defined
    auto* path = this->top_path_if_exists();
    // I'm not sure if that's what the function is supposed to do.
    if(path) {
        if(!path->closed) {
            path->closed = true;
        }
    }
    this->drawing_commands.push_back({this->create_new_path()});

    #endif
}

void Renderer::set_source_rgb(double r, double g, double b) {
    #ifndef __EMSCRIPTEN__
    cairo_set_source_rgb(cr, r, g, b);
    #else // __EMSCRIPTEN__ defined
    this->style.color.r = r;
    this->style.color.g = g;
    this->style.color.b = b;
    this->style.color.a = 1.0f;
    #endif
}

void Renderer::set_source_rgba(double r, double g, double b, double a) {
    #ifndef __EMSCRIPTEN__
    cairo_set_source_rgba(cr, r, g, b, a);
    #else // __EMSCRIPTEN__ defined
    this->style.color.r = r;
    this->style.color.g = g;
    this->style.color.b = b;
    this->style.color.a = a;
    #endif
}

void Renderer::set_line_width(double width) {
    #ifndef __EMSCRIPTEN__
    cairo_set_line_width(cr, width);
    #else // __EMSCRIPTEN__ defined
    this->style.line_width = width;
    #endif
}

#ifndef __EMSCRIPTEN__
std::string Renderer::text_span_to_pango_markup(const TextSpan& span, const std::optional<TextStyle>& parent_style) const {
    std::string ret;
    auto should_specify_style = [&](){
        if(span.specifies_style) {

            if(parent_style.has_value()) {
                auto& pstyle = parent_style.value();
                //strange compiler error
                //return pstyle != span.style;
                return true;
            } else {
                return true;
            }
        }
        return false;
    };
    bool style_block = should_specify_style();
    if(style_block) {
        const auto& style = span.style;
        ret += "<span line_height=\"0.75\"";
        if(style.specifies_color) {
            std::stringstream html_color;
            html_color << "#";
            html_color << std::hex << std::setfill('0') << std::setw(2) << (int)(style.color.r * 255);
            html_color << std::hex << std::setfill('0') << std::setw(2) << (int)(style.color.g * 255);
            html_color << std::hex << std::setfill('0') << std::setw(2) << (int)(style.color.b * 255);
            html_color << std::hex << std::setfill('0') << std::setw(2) << (int)(style.color.a * 255);
            // std::cout<<html_color.str()<<' '<<style.color.r<<' '<<style.color.g<<' '<<style.color.b<<'\n';
            ret += "color=\"" + html_color.str() + "\" ";
        }
        if(!style.size.empty()) {
            ret += "size=\"" + style.size + "\" ";
        }
        if(!style.weight.empty()) {
            ret += "weight=\"" + style.weight + "\" ";
        }
        ret += ">";
        switch(style.positioning) {
            case TextPositioning::Sub: {
                ret += "<sub>";
                break;
            }
            case TextPositioning::Super: {
                // The string below begins with 
                // the invisible U+200B unicode character.
                // This is a workaround for what's likely 
                // a bug in pango font rendering engine.
                // Without it, the superscript is relative 
                // to the subscript (atom count)
                // instead of the atom's symbol
                ret += "​<sup>";
                break;
            }
            default:
            case TextPositioning::Normal: {
                // nothing
            }
        }
    }
    if(span.has_subspans()) {
        const auto& subspans = span.as_subspans();
        for(const auto& subspan: subspans) {
            std::optional<TextStyle> styleopt = subspan.specifies_style ? subspan.style : parent_style;
            ret += this->text_span_to_pango_markup(subspan, styleopt);
        }
    } else if(span.is_newline()) {
        ret += "\n";
    } else {
        const auto& caption = span.as_caption();
        // todo: escape characters!
        ret += caption;
    }
    if(style_block) {
        switch(span.style.positioning) {
            case TextPositioning::Sub: {
                ret += "</sub>";
                break;
            }
            case TextPositioning::Super: {
                ret += "</sup>";
                break;
            }
            default:
            case TextPositioning::Normal: {
                // nothing
            }
        }
        ret += "</span>";
    }
    return ret;
}
#else
std::optional<Renderer::TextSize> Renderer::TextMeasurementCache::lookup_span(const Renderer::TextSpan& text) const {
    return this->lookup_span(std::hash<TextSpan>{}(text));
}

std::optional<Renderer::TextSize> Renderer::TextMeasurementCache::lookup_span(Renderer::TextMeasurementCache::hash_t span_hash) const {
    auto it = this->cache.find(span_hash);
    if(it != this->cache.end()) {
        return it->second;
    }
    return std::nullopt;
}

void Renderer::TextMeasurementCache::add(Renderer::TextMeasurementCache::hash_t span_hash, Renderer::TextSize value) {
    this->cache.emplace(span_hash, value);
}

void Renderer::TextMeasurementCache::add(const Renderer::TextSpan& text, Renderer::TextSize value) {
    this->add(std::hash<TextSpan>{}(text), value);
}
std::size_t Renderer::TextMeasurementCache::size() const {
    return this->cache.size();
}

#endif

Renderer::TextSize Renderer::measure_text(const Renderer::TextSpan& text) {
    #ifndef __EMSCRIPTEN__
    std::string markup = this->text_span_to_pango_markup(text);
    pango_layout_set_markup(this->pango_layout, markup.c_str(), -1);
    TextSize ret;
    pango_layout_get_pixel_size(this->pango_layout, &ret.width, &ret.height);
    return ret;
    #else // __EMSCRIPTEN__ defined
    std::optional<TextMeasurementCache::hash_t> text_hash;
    if(this->tm_cache) {
        text_hash = std::hash<TextSpan>{}(text);
        auto cached_opt = this->tm_cache->lookup_span(text_hash.value());
        if(cached_opt.has_value()) {
            return cached_opt.value();
        }
    }
    // return {0,0};
    // The try..catch doesn't work for me.
    // try {
    // g_info("Measuring...");
    Renderer::Text wtext;
    wtext.origin.x = 0;
    wtext.origin.y = 0;
    wtext.spans.push_back(text);
    // g_info("Wrapper text has been built.");
    emscripten::val result = this->text_measurement_function(wtext);
    // g_info("Got result.");
    auto ret = result.as<TextSize>();
    // } catch(...) {
    //     return {0,0};
    // }
    if(this->tm_cache) {
        g_debug("About to hash TextSpan for caching.");
        this->tm_cache->add(text_hash.value(), ret);
        g_debug("TextSpan added to cache. Cache entries: %lu", this->tm_cache->size());
    }
    return ret;
    #endif
}

void Renderer::show_text(const Renderer::TextSpan& text_span) {
    #ifndef __EMSCRIPTEN__
    std::string markup = this->text_span_to_pango_markup(text_span);
    pango_layout_set_markup(this->pango_layout, markup.c_str(), -1);
    pango_cairo_show_layout(this->cr, this->pango_layout);
    #else // Lhasa
    Text text;
    if(text_span.has_subspans()) {
        // This is a deep copy. Yikes.
        text.spans = text_span.as_subspans();
        text.style = text_span.style;
    } else {
        text.spans.push_back(text_span);
        // Let's leave the style as it is, for now.
        // text.style
    }
    text.origin = this->position;
    this->drawing_commands.push_back(DrawingCommand{text});
    #endif
}

Renderer::~Renderer() {
    #ifndef __EMSCRIPTEN__
    g_object_unref(this->pango_layout);
    cairo_destroy(this->cr);
    #else // __EMSCRIPTEN__ defined
    // Lhasa-specific includes/definitions
    #endif
}

MoleculeRenderContext::MoleculeRenderContext(const CanvasMolecule& cm, Renderer& ren, DisplayMode mode, const std::pair<int, int>& viewport_offset, float canvas_scale) 
:canvas_molecule(cm), ren(ren), display_mode(mode) {
    scale_factor = CanvasMolecule::BASE_SCALE_FACTOR * canvas_scale;
    x_offset = scale_factor * canvas_molecule.x_canvas_translation - viewport_offset.first;
    y_offset = scale_factor * canvas_molecule.y_canvas_translation - viewport_offset.second;
}

MoleculeRenderContext::~MoleculeRenderContext() {

}

std::tuple<Renderer::TextSpan, bool, bool> MoleculeRenderContext::process_appendix(const std::string& symbol, const std::optional<Atom::Appendix>& appendix, const Renderer::TextStyle& inherited_style) {
    Renderer::TextSpan ret((std::vector<Renderer::TextSpan>()));
    /// Contains the symbol of the main atom
    Renderer::TextSpan symbol_span(symbol);
    bool reversed = false;
    bool vertical = false;
    if(!appendix.has_value()) {
        ret.as_subspans().push_back(symbol_span);
    } else {
        const auto& ap = appendix.value();
        //ret += "<span>";
        /// Contains symbols and indices of superatoms
        Renderer::TextSpan root_span((std::vector<Renderer::TextSpan>()));
        Renderer::TextSpan superatoms_symbol_span;
        auto make_index_span = [=](){
            Renderer::TextSpan index_span;
            index_span.specifies_style = true;
            index_span.style = inherited_style;
            index_span.style.positioning = Renderer::TextPositioning::Sub;
            return index_span;
        };
        auto index_span = make_index_span();

        for(auto i = ap.superatoms.begin(); i != ap.superatoms.end(); i++) {
            if(std::isdigit(*i)) {
                if(!superatoms_symbol_span.as_caption().empty()) {
                    root_span.as_subspans().push_back(std::move(superatoms_symbol_span));
                    superatoms_symbol_span = Renderer::TextSpan();
                }
                index_span.as_caption().push_back(*i);
            } else {
                if(!index_span.as_caption().empty()) {
                    root_span.as_subspans().push_back(std::move(index_span));
                    index_span = make_index_span();
                }
                superatoms_symbol_span.as_caption().push_back(*i);
            }
        }
        if(!superatoms_symbol_span.as_caption().empty()) {
            root_span.as_subspans().push_back(std::move(superatoms_symbol_span));
        }
        if(!index_span.as_caption().empty()) {
            root_span.as_subspans().push_back(std::move(index_span));
        }
        auto process_charge_span = [&](){
            if(ap.charge != 0) {
                Renderer::TextSpan charge_span;
                charge_span.specifies_style = true;
                charge_span.style = inherited_style;
                charge_span.style.positioning = Renderer::TextPositioning::Super;

                unsigned int charge_no_sign = std::abs(ap.charge);
                if(charge_no_sign > 1) {
                    charge_span.as_caption() += std::to_string(charge_no_sign);
                }
                charge_span.as_caption().push_back(ap.charge > 0 ? '+' : '-');
                ret.as_subspans().push_back(charge_span);
            }
        };
        if(ap.vertical) {
            vertical = true;
        }
        if (ap.reversed) {
            reversed = true;

            ret.as_subspans().push_back(root_span);
            if(vertical) {
                ret.as_subspans().push_back(Renderer::TextSpan(Renderer::TextSpan::Newline{}));
            }
            ret.as_subspans().push_back(symbol_span);
            process_charge_span();
        } else {
            ret.as_subspans().push_back(symbol_span);
            if(vertical) {
                process_charge_span();
                ret.as_subspans().push_back(Renderer::TextSpan(Renderer::TextSpan::Newline{}));
            }
            ret.as_subspans().push_back(root_span);
            if(!vertical) {
                process_charge_span();
            }
        }
        //ret += "</span>";
    }
    return std::make_tuple(ret, reversed, vertical);
}

std::pair<unsigned int,graphene_rect_t> MoleculeRenderContext::render_atom(const CanvasMolecule::Atom& atom, DisplayMode render_mode) {
    // pre-process text
    auto [r,g,b] = CanvasMolecule::atom_color_to_rgb(atom.color);
    // const std::string color_str = CanvasMolecule::atom_color_to_html(atom.color);
    
    // Span for the whole thing - includes symbol, appendix, index etc.
    Renderer::TextSpan atom_span((std::vector<Renderer::TextSpan>()));
    atom_span.specifies_style = true;
    atom_span.style.specifies_color = true;
    atom_span.style.color.r = r;
    atom_span.style.color.g = g;
    atom_span.style.color.b = b;
    #ifndef __EMSCRIPTEN__
    atom_span.style.size = render_mode != DisplayMode::AtomIndices ? "x-large" : "medium";
    #else
    atom_span.style.size = render_mode != DisplayMode::AtomIndices ? "medium" : "small";
    #endif
    atom_span.style.weight = atom.highlight != 0 ? "bold" : "normal";

    // Span for the atom symbol - solely.
    // This allows us to measure the size of the atom's symbol 
    // and then properly center/align the whole text
    Renderer::TextSpan raw_atom_span;
    raw_atom_span.specifies_style = true;
    raw_atom_span.style = atom_span.style;

    bool reversed = false;
    bool vertical = false;

    switch (render_mode) {
        case DisplayMode::AtomIndices: {
            raw_atom_span.as_caption() += atom.symbol;
            atom_span
                .as_subspans()
                .push_back(
                    Renderer::TextSpan(std::string(atom.symbol + ":" + std::to_string(atom.idx)))
                );
            break;
        }
        case DisplayMode::AtomNames: {
            if(atom.name.has_value()) {
                std::string atom_name = atom.name.value();
                raw_atom_span.as_caption() += atom_name;
                atom_span = raw_atom_span;
                break;
            } 
            // break;
            // We want to fall back to the standard case if the atom has no name.
        }
        default:
        case DisplayMode::Standard: {
            auto [appendix, p_reversed, p_vertical] = process_appendix(atom.symbol, atom.appendix, atom_span.style);
            reversed = p_reversed;
            vertical = p_vertical;
            raw_atom_span.as_caption() += atom.symbol;
            atom_span
                .as_subspans()
                .push_back(appendix);
            break;
        }
    }

    // Used to make the texts centered where they should be (appendix).
    // Measure the size of the "main" atom, without "appendix"
    Renderer::TextSize raw_size = ren.measure_text(raw_atom_span);
    // Measure full size of the text
    Renderer::TextSize size = ren.measure_text(atom_span);

    // g_info("Measurement results: raw %ix%i n %ix%i", raw_size.width, raw_size.height, size.width, size.height);

    #ifndef __EMSCRIPTEN__
    // todo: get rid of this '5' magic number - figure out what's wrong
    const int magic1 = 5;
    // For Lhasa
    const int magic2 = 0;
    // Magic number. This should be removed.
    // Workaround for pango giving us too high layout size.
    const float layout_too_high = 3.f;
    #else
    const int magic1 = 0;
    // Text appears too high. Manuall offset.
    const int magic2 = 15;
    const float layout_too_high = 0.f;
    #endif
    int layout_x_offset = reversed && !vertical ? size.width - raw_size.width / 2.f + magic1 : raw_size.width / 2.f;
    int layout_y_offset = vertical && reversed ? size.height - raw_size.height / 2.f : raw_size.height / 2.f;
    double origin_x = atom.x * scale_factor + x_offset - layout_x_offset;
    double origin_y = atom.y * scale_factor + y_offset - layout_y_offset;

    graphene_rect_t rect;
    rect.origin.x = origin_x;
    rect.origin.y = origin_y + layout_too_high;
    rect.size.width = size.width;
    rect.size.height = size.height - layout_too_high * 2.f;

    // g_info("Rect: x=%f, y=%f, width=%f, height=%f", rect.origin.x, rect.origin.y, rect.size.width, rect.size.height);

    // highlight
    process_atom_highlight(atom);
    // text
    ren.move_to(origin_x, origin_y + magic2);
    ren.show_text(atom_span);

    return std::make_pair(atom.idx, rect);
}

void MoleculeRenderContext::draw_atoms() {
    for(const auto& atom: canvas_molecule.atoms) {
        switch (display_mode) {
            case DisplayMode::AtomIndices: {
                atom_idx_to_canvas_rect.emplace(render_atom(atom,DisplayMode::AtomIndices));
                break;
            }
            case DisplayMode::AtomNames: {
                if(atom.name.has_value()) {
                    atom_idx_to_canvas_rect.emplace(render_atom(atom,DisplayMode::AtomNames));
                    break;
                }
                // We want to fall back to the standard case if the atom has no name
                // break;
            }
            default:
            case DisplayMode::Standard: {
                if(atom.symbol == "C") {
                    if(atom.appendix.has_value()) {
                        atom_idx_to_canvas_rect.emplace(render_atom(atom));
                    } else {
                        process_atom_highlight(atom);
                    }
                } else {
                    atom_idx_to_canvas_rect.emplace(render_atom(atom));
                }
                break;
            }
        }
    }
}

void MoleculeRenderContext::process_atom_highlight(const Atom& atom) {
    auto highlight = CanvasMolecule::determine_dominant_highlight(atom.highlight);
    if(highlight) {
        auto [r, g, b] = CanvasMolecule::hightlight_to_rgb(*highlight);
        ren.new_sub_path();
        // Go to the center of the atom
        //ren.move_to(atom.x * scale_factor + x_offset + CanvasMolecule::ATOM_HITBOX_RADIUS, atom.y * scale_factor + y_offset);
        ren.set_source_rgb(r, g, b);
        ren.arc(atom.x * scale_factor + x_offset, atom.y * scale_factor + y_offset, CanvasMolecule::ATOM_HITBOX_RADIUS, 0, M_PI * 2.0);
        ren.stroke_preserve();
        ren.set_source_rgba(r, g, b, 0.5);
        ren.fill();
    }
}

std::pair<graphene_point_t,graphene_point_t> MoleculeRenderContext::cropped_bond_coords(const graphene_point_t& first_atom, unsigned int first_atom_idx, const graphene_point_t& second_atom, unsigned int second_atom_idx) {
    // We pass in the bond vector so that it always points "away" from the point
    auto crop_line_against_rect = [](const graphene_rect_t& rect, float bond_vec_x, float bond_vec_y, const graphene_point_t& point){
        graphene_point_t ret = point;
        if(rect.origin.x < point.x && rect.origin.x + rect.size.width >= point.x
        && rect.origin.y < point.y && rect.origin.y + rect.size.height >= point.y) { // inside rectangle
            // We need to use the point of intersection formula
            // to find the intersection point between the bond and the edges.
            // We solve against two edges and then pick the solution lying the closest to 'point'.

            // Line from bond
            graphene_point_t point2 = point;
            point2.x += bond_vec_x;
            point2.y += bond_vec_y;
            const float bond_a_quot = (point.x - point2.x);
            const float bond_a = bond_a_quot == 0 ? -point.x : (point.y - point2.y) / bond_a_quot;
            const float bond_b = -1;
            const float bond_c = point.y - bond_a * point.x;

            // Pick line from the relevant vertical edge
            float vert_egde_c;
            const float vert_edge_a = 1;
            const float vert_edge_b = 0;
            if(bond_vec_x > 0) { // solve against right side
                vert_egde_c = -(rect.origin.x + rect.size.width);
            } else { // solve against left side
                vert_egde_c = -rect.origin.x;
            }

            // Pick line from the relevant horizontal edge
            float horiz_egde_c;
            const float horiz_edge_a = 0;
            const float horiz_edge_b = 1;
            if(bond_vec_y > 0) { // solve against the bottom
                horiz_egde_c = -(rect.origin.y + rect.size.height);
            } else { // solve against the top
                horiz_egde_c = -rect.origin.y;
            }
            graphene_point_t vert_edge_solution;
            vert_edge_solution.x = (bond_b * vert_egde_c - vert_edge_b * bond_c)/(bond_a * vert_edge_b - vert_edge_a * bond_b);
            vert_edge_solution.y = (vert_edge_a * bond_c - bond_a * vert_egde_c)/(bond_a * vert_edge_b - vert_edge_a * bond_b);
            graphene_point_t horiz_edge_solution;
            horiz_edge_solution.x = (bond_b * horiz_egde_c - horiz_edge_b * bond_c)/(bond_a * horiz_edge_b - horiz_edge_a * bond_b);
            horiz_edge_solution.y = (horiz_edge_a * bond_c - bond_a * horiz_egde_c)/(bond_a * horiz_edge_b - horiz_edge_a * bond_b);

            // if vert_edge_solution is farther from the point than horiz_edge_solution
            if(std::pow(vert_edge_solution.x - point.x, 2.f) + std::pow(vert_edge_solution.y - point.y, 2.f) 
                > std::pow(horiz_edge_solution.x - point.x, 2.f) + std::pow(horiz_edge_solution.y - point.y, 2.f)) {
                ret = horiz_edge_solution;
            } else {
                ret = vert_edge_solution;
            }
        } 
        return ret;
    };

    graphene_point_t a = first_atom;
    graphene_point_t b = second_atom;

    float bond_vec_x = second_atom.x - first_atom.x;
    float bond_vec_y = second_atom.y - first_atom.y;

    auto first_rect_iter = atom_idx_to_canvas_rect.find(first_atom_idx);
    if(first_rect_iter != atom_idx_to_canvas_rect.end()) {
        a = crop_line_against_rect(first_rect_iter->second, bond_vec_x, bond_vec_y, first_atom);
    }
    auto second_rect_iter = atom_idx_to_canvas_rect.find(second_atom_idx);
    if(second_rect_iter != atom_idx_to_canvas_rect.end()) {
        b = crop_line_against_rect(second_rect_iter->second, -bond_vec_x, -bond_vec_y, second_atom);
    }
    
    return std::make_pair(a,b);
}

void MoleculeRenderContext::draw_central_bond_line(const CanvasMolecule::Bond& bond) {
    graphene_point_t first_atom;
    first_atom.x = bond.first_atom_x * scale_factor + x_offset;
    first_atom.y = bond.first_atom_y * scale_factor + y_offset;

    graphene_point_t second_atom;
    second_atom.x = bond.second_atom_x * scale_factor + x_offset;
    second_atom.y = bond.second_atom_y * scale_factor + y_offset;

    auto [first,second] = cropped_bond_coords(first_atom,bond.first_atom_idx,second_atom,bond.second_atom_idx);
    
    ren.move_to(first.x, first.y);
    ren.line_to(second.x, second.y);
    ren.stroke();
}

void MoleculeRenderContext::draw_straight_wedge(const CanvasMolecule::Bond& bond, bool reversed) {
    graphene_point_t origin;
    origin.x = reversed ? bond.first_atom_x : bond.second_atom_x;
    origin.y = reversed ? bond.first_atom_y : bond.second_atom_y;
    auto origin_idx = reversed ? bond.first_atom_idx : bond.second_atom_idx;

    graphene_point_t target;
    target.x = reversed ? bond.second_atom_x : bond.first_atom_x;
    target.y = reversed ? bond.second_atom_y : bond.first_atom_y;
    auto target_idx = reversed ? bond.second_atom_idx : bond.first_atom_idx;

    origin.x *= scale_factor;
    origin.y *= scale_factor;
    target.x *= scale_factor;
    target.y *= scale_factor;

    origin.x += x_offset;
    origin.y += y_offset;
    target.x += x_offset;
    target.y += y_offset;

    auto [origin_cropped,target_cropped] = cropped_bond_coords(origin, origin_idx, target, target_idx);
    
    auto [pv_x,pv_y] = bond.get_perpendicular_versor();
    auto cropped_bond_len = std::sqrt(std::pow(target_cropped.x - origin_cropped.x, 2.f) + std::pow(target_cropped.y - origin_cropped.y, 2.f));
    auto v_x = pv_x * std::sin(GEOMETRY_BOND_SPREAD_ANGLE / 2.f) * cropped_bond_len;
    auto v_y = pv_y * std::sin(GEOMETRY_BOND_SPREAD_ANGLE / 2.f) * cropped_bond_len;

    ren.new_path();
    ren.move_to(origin_cropped.x, origin_cropped.y);
    ren.line_to(target_cropped.x + v_x, target_cropped.y + v_y);
    ren.stroke_preserve();

    ren.line_to(target_cropped.x - v_x, target_cropped.y - v_y);
    ren.stroke_preserve();

    ren.line_to(origin_cropped.x, origin_cropped.y);
    ren.stroke_preserve();

    ren.close_path();
    ren.fill();
}

void MoleculeRenderContext::draw_straight_dashed_bond(const CanvasMolecule::Bond& bond, bool reversed) {
    graphene_point_t origin;
    origin.x = reversed ? bond.first_atom_x : bond.second_atom_x;
    origin.y = reversed ? bond.first_atom_y : bond.second_atom_y;
    auto origin_idx = reversed ? bond.first_atom_idx : bond.second_atom_idx;

    graphene_point_t target;
    target.x = reversed ? bond.second_atom_x : bond.first_atom_x;
    target.y = reversed ? bond.second_atom_y : bond.first_atom_y;
    auto target_idx = reversed ? bond.second_atom_idx : bond.first_atom_idx;

    origin.x *= scale_factor;
    origin.y *= scale_factor;
    target.x *= scale_factor;
    target.y *= scale_factor;

    origin.x += x_offset;
    origin.y += y_offset;
    target.x += x_offset;
    target.y += y_offset;

    auto [current,target_cropped] = cropped_bond_coords(origin, origin_idx, target, target_idx);
    
    auto [pv_x,pv_y] = bond.get_perpendicular_versor();
    auto cropped_bond_len = std::sqrt(std::pow(target_cropped.x - current.x, 2.f) + std::pow(target_cropped.y - current.y, 2.f));

    float dashes = cropped_bond_len / (GEOMETRY_BOND_DASH_SEPARATION * scale_factor);
    unsigned int full_dashes = std::floor(dashes);

    float step_x = (target_cropped.x - current.x) / dashes;
    float step_y = (target_cropped.y - current.y) / dashes;

    auto v_x = pv_x * std::sin(GEOMETRY_BOND_SPREAD_ANGLE / 2.f) * cropped_bond_len;
    auto v_y = pv_y * std::sin(GEOMETRY_BOND_SPREAD_ANGLE / 2.f) * cropped_bond_len;
    
    for(unsigned int i = 0; i <= full_dashes; i++) {
        float spread_multiplier = (float) i / dashes;
        ren.move_to(current.x - v_x * spread_multiplier, current.y - v_y * spread_multiplier);
        ren.line_to(current.x + v_x * spread_multiplier, current.y + v_y * spread_multiplier);
        ren.stroke();
        current.x += step_x;
        current.y += step_y;
    }
}

void MoleculeRenderContext::draw_wavy_bond(const CanvasMolecule::Bond& bond) {
    graphene_point_t first_atom;
    first_atom.x = bond.first_atom_x * scale_factor + x_offset;
    first_atom.y = bond.first_atom_y * scale_factor + y_offset;

    graphene_point_t second_atom;
    second_atom.x = bond.second_atom_x * scale_factor + x_offset;
    second_atom.y = bond.second_atom_y * scale_factor + y_offset;

    auto [first,second] = cropped_bond_coords(first_atom,bond.first_atom_idx,second_atom,bond.second_atom_idx);
    auto full_vec_x = second.x - first.x;
    auto full_vec_y = second.y - first.y;

    const float wave_arc_radius = WAVY_BOND_ARC_LENGTH * scale_factor / 2.f;

    // The angle at which the bond points
    float base_angle = std::atan(full_vec_y / full_vec_x);
    float arcs_count = std::sqrt(std::pow(full_vec_x,2.f) + std::pow(full_vec_y,2.f)) / (WAVY_BOND_ARC_LENGTH * scale_factor);
    unsigned int rounded_arcs_count = std::floor(arcs_count);
    float step_x = full_vec_x / arcs_count;
    float step_y = full_vec_y / arcs_count;
    float current_x = first.x + step_x / 2.f;
    float current_y = first.y + step_y / 2.f;
    // It seems that for positive base_angle, 
    // 'true' is counter-clockwise
    // and 'false is clockwise.
    // For negative base_angle, it's the opposite.
    bool arc_direction = true;
    // for debugging
    // float l_angle_one = 0, l_angle_two = 0;

    // Two core angles for semi-circles
    float p1 = base_angle;
    float p2 = base_angle - M_PI;

    float angle_one, angle_two;

    for (unsigned int i = 0; i < rounded_arcs_count; i++) {
        float next_x = current_x + step_x;
        float next_y = current_y + step_y;
        
        if(arc_direction) {
            angle_one = p1;
            angle_two = p2;
        } else {
            angle_one = p2;
            angle_two = p1;
        }
        // l_angle_one = angle_one;
        // l_angle_two = angle_two;
        ren.new_sub_path();
        ren.arc(current_x, current_y, wave_arc_radius, angle_one, angle_two);
        ren.stroke();

        current_x = next_x;
        current_y = next_y;
        arc_direction = !arc_direction;
    }
    // Final part of the path. Truncated arc.
    float partial_arc_proportion = arcs_count - (float) rounded_arcs_count;
    float arccos_arg = 1.f - (partial_arc_proportion / WAVY_BOND_ARC_LENGTH / 2.f);
    
    // This is the angle for the final arc.
    float theta = std::acos(arccos_arg);

    // The magic behind 'step_x > 0'
    // is a bit of a mystery (derived empirically).
    // There's certain correlation with the sign of base_angle
    // but that's not the whole story.
    // What matters is that it allows for differentiating
    // between various cases, with angles from different quarters.
    float starting_angle = step_x > 0 ? p2 : p1;
    if(arc_direction) {
        if(step_x > 0) {
            angle_one = starting_angle - theta;
            angle_two = starting_angle;
        } else {
            angle_one = starting_angle;
            angle_two = starting_angle + theta;
        }
    } else {
        if(step_x > 0) {
            angle_one = starting_angle;
            angle_two = starting_angle + theta;
        } else {
            angle_one = starting_angle - theta;
            angle_two = starting_angle;
        }
    }

    // debugging stuff

    // std::string case_info;
    // if(theta > M_PI_2) {
    //     case_info += "T";
    // } else {
    //     case_info += "G";
    // }
    // if(base_angle > 0) {
    //     case_info += "A";
    // } else {
    //     case_info += "E";
    // }
    // if(arc_direction) {
    //     case_info += "K";
    // } else {
    //     case_info += "V";
    // }
    // case_info += angle_two - angle_one > 0 ? "O" : "U";
    // if(step_x > 0) {
    //     case_info += "Z";
    // } else {
    //     case_info += "V";
    // }
    // g_debug(
    //     "theta=%f, base_angle=%f a1=%f, a2=%f a2-a1=%f abs(a2-a1)=%f direction=%s p1=%f, p2=%f case_codename=%s",
    //     theta / M_PI * 180.f,
    //     base_angle / M_PI * 180.f,
    //     angle_one / M_PI * 180.f,
    //     angle_two / M_PI * 180.f,
    //     (angle_two - angle_one) / M_PI * 180.f,
    //     std::fabs(angle_two - angle_one) / M_PI * 180.f,
    //     arc_direction ? "true" : "false",
    //     l_angle_one / M_PI * 180.f,
    //     l_angle_two / M_PI * 180.f,
    //     case_info.c_str()
    // );
    ren.new_sub_path();
    ren.arc(current_x, current_y, wave_arc_radius, angle_one, angle_two);
    ren.stroke();
}

void MoleculeRenderContext::draw_side_bond_line(const CanvasMolecule::Bond& bond, bool addOrSub, std::optional<float> first_shortening_proportion, std::optional<float> second_shortening_proportion) {
    auto [pv_x,pv_y] = bond.get_perpendicular_versor();
    if (!addOrSub) { // change sign of the versor
        pv_x *= -1.f;
        pv_y *= -1.f;
    }
    // Convert the versor to a vector of the desired length
    pv_x *= CanvasMolecule::BOND_LINE_SEPARATION;
    pv_y *= CanvasMolecule::BOND_LINE_SEPARATION;

    auto [bond_vec_x, bond_vec_y] = bond.get_vector();
    auto first_x = bond.first_atom_x;
    auto second_x = bond.second_atom_x;
    auto first_y = bond.first_atom_y;
    auto second_y = bond.second_atom_y;

    if(first_shortening_proportion.has_value()) {
        first_x += first_shortening_proportion.value() * bond_vec_x;
        first_y += first_shortening_proportion.value() * bond_vec_y;
    }
    if(second_shortening_proportion.has_value()) {
        second_x -= second_shortening_proportion.value() * bond_vec_x;
        second_y -= second_shortening_proportion.value() * bond_vec_y;
    }

    // Points a and b represent the off-center bond before cropping.
    graphene_point_t a;
    a.x = (first_x + pv_x) * scale_factor + x_offset;
    a.y = (first_y + pv_y) * scale_factor + y_offset;
    graphene_point_t b;
    b.x = (second_x + pv_x) * scale_factor + x_offset;
    b.y = (second_y + pv_y) * scale_factor + y_offset;

    // We need to make sure that the off-center bond 
    // after cropping is not going to be longer than the center bond
    graphene_point_t first_atom_centered;
    first_atom_centered.x = first_x * scale_factor + x_offset;
    first_atom_centered.y = first_y * scale_factor + y_offset;

    graphene_point_t second_atom_centered;
    second_atom_centered.x = second_x * scale_factor + x_offset;
    second_atom_centered.y = second_y * scale_factor + y_offset;

    /// Centered bond cropped
    auto [first_c,second_c] = cropped_bond_coords(
        first_atom_centered,
        bond.first_atom_idx,
        second_atom_centered,
        bond.second_atom_idx
    );
    // Now we offset the center bond after cropping
    first_c.x += pv_x * scale_factor;
    first_c.y += pv_y * scale_factor;
    second_c.x += pv_x * scale_factor;
    second_c.y += pv_y * scale_factor;

    // Points a_cropped and b_cropped represent the off-center bond after cropping.
    auto [a_cropped,b_cropped] = cropped_bond_coords(a,bond.first_atom_idx,b,bond.second_atom_idx);

    // Now we need to make sure that the off-center bond 
    // after cropping is not going to be longer than the center bond
    if(bond_vec_x > 0) {
        // The beginning is shorter for the centered bond
        if(first_c.x > a_cropped.x) {
            a_cropped = first_c;
        }
        // The end is shorter for the centered bond
        if(second_c.x < b_cropped.x) {
            b_cropped = second_c;
        }
    } else {
        // The beginning is shorter for the centered bond
        if(first_c.x < a_cropped.x) {
            a_cropped = first_c;
        }
        // The end is shorter for the centered bond
        if(second_c.x > b_cropped.x) {
            b_cropped = second_c;
        }
    }
    if(bond_vec_y > 0) {
        // The beginning is shorter for the centered bond
        if(first_c.y > a_cropped.y) {
            a_cropped = first_c;
        }
        // The end is shorter for the centered bond
        if(second_c.y < b_cropped.y) {
            b_cropped = second_c;
        }
    } else {
        // The beginning is shorter for the centered bond
        if(first_c.y < a_cropped.y) {
            a_cropped = first_c;
        }
        // The end is shorter for the centered bond
        if(second_c.y > b_cropped.y) {
            b_cropped = second_c;
        }
    }
    ren.move_to(a_cropped.x, a_cropped.y);
    ren.line_to(b_cropped.x, b_cropped.y);
    ren.stroke();
}

void MoleculeRenderContext::draw_centered_double_bond(const CanvasMolecule::Bond& bond) {
    auto [pv_x,pv_y] = bond.get_perpendicular_versor();

    // Convert the versor to a vector of the desired length
    pv_x *= CENTERED_DOUBLE_BOND_LINE_SEPARATION / 2.f * scale_factor;
    pv_y *= CENTERED_DOUBLE_BOND_LINE_SEPARATION / 2.f * scale_factor;

    graphene_point_t first_atom;
    first_atom.x = bond.first_atom_x * scale_factor + x_offset;
    first_atom.y = bond.first_atom_y * scale_factor + y_offset;

    graphene_point_t second_atom;
    second_atom.x = bond.second_atom_x * scale_factor + x_offset;
    second_atom.y = bond.second_atom_y * scale_factor + y_offset;

    auto [first,second] = cropped_bond_coords(first_atom,bond.first_atom_idx,second_atom,bond.second_atom_idx);

    ren.move_to(first.x + pv_x, first.y + pv_y);
    ren.line_to(second.x + pv_x, second.y + pv_y);
    ren.stroke();

    ren.move_to(first.x - pv_x, first.y - pv_y);
    ren.line_to(second.x - pv_x, second.y - pv_y);
    ren.stroke();
}


void MoleculeRenderContext::draw_bonds() {
    for(const auto& bond: canvas_molecule.bonds) {
        auto highlight = CanvasMolecule::determine_dominant_highlight(bond->highlight);
        if(highlight) {
            ren.set_line_width(4.0);
            auto [r, g, b] = CanvasMolecule::hightlight_to_rgb(*highlight);
            ren.set_source_rgb(r, g, b);
        } else {
            ren.set_line_width(2.0);
            ren.set_source_rgb(0.0, 0.0, 0.0);
        }

        if(bond->geometry != BondGeometry::Flat && bond->type == BondType::Single) {
            switch (bond->geometry) {
                default:
                case BondGeometry::Unspecified:{
                    draw_wavy_bond(*bond.get());
                    break;
                }
                case BondGeometry::WedgeTowardsFirst:{
                    draw_straight_wedge(*bond.get(), true);
                    g_warning_once("todo: rendering bond geometry in rings");
                    break;
                }
                case BondGeometry::WedgeTowardsSecond:{
                    draw_straight_wedge(*bond.get(), false);
                    g_warning_once("todo: rendering bond geometry in rings");
                    break;
                }
                case BondGeometry::DashedTowardsFirst:{
                    draw_straight_dashed_bond(*bond.get(), true);
                    g_warning_once("todo: rendering bond geometry in rings");
                    break;
                }
                case BondGeometry::DashedTowardsSecond:{
                    draw_straight_dashed_bond(*bond.get(), false);
                    g_warning_once("todo: rendering bond geometry in rings");
                    break;
                }
            }
        } else {
            switch(bond->type) {
                case BondType::Double:{
                    DoubleBondDrawingDirection direction = 
                        bond->bond_drawing_direction.has_value() ? 
                        bond->bond_drawing_direction.value() 
                        : DoubleBondDrawingDirection::Primary;
                    bool direction_as_bool = true;

                    switch (direction) { 
                        case DoubleBondDrawingDirection::Secondary:{
                            direction_as_bool = false;
                            // no break here.
                        }
                        case DoubleBondDrawingDirection::Primary:{
                            draw_central_bond_line(*bond.get());
                            draw_side_bond_line(
                                *bond.get(),
                                direction_as_bool,
                                bond->first_shortening_proportion,
                                bond->second_shortening_proportion
                            );
                            break;
                        }
                        case DoubleBondDrawingDirection::Centered:{
                            draw_centered_double_bond(*bond.get());
                            break;
                        }
                    }
                    break;
                }
                case BondType::Triple:{
                    draw_central_bond_line(*bond.get());
                    g_warning_once("todo: Triple bonds might need truncating too.");
                    // "to the left"
                    draw_side_bond_line(*bond.get(), false, std::nullopt,std::nullopt);
                    // "to the right"
                    draw_side_bond_line(*bond.get(), true, std::nullopt,std::nullopt);
                    break;
                }
                default:
                case BondType::Single:{
                    draw_central_bond_line(*bond.get());
                    break;
                }
            }
        }
    }
}