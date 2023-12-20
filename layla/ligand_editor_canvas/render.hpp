/* layla/ligand_editor_canvas/render.hpp
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
#ifndef COOT_LIGAND_EDITOR_CANVAS_RENDER_HPP
#define COOT_LIGAND_EDITOR_CANVAS_RENDER_HPP

#ifndef __EMSCRIPTEN__
    #include <gtk/gtk.h>
    #include <pango/pango-layout.h>
#else // Lhasa-specific includes
    #include <graphene.h>
    #include <variant>
    #include <optional>
    #include <vector>
    #include <emscripten/val.h>
#endif
#include <map>
#include <memory>
#include <tuple>
#include "model.hpp"

namespace coot::ligand_editor_canvas::impl {


struct Renderer {
    struct Color {
        double r, g, b, a;
    };

    enum class TextPositioning: unsigned char {
        Normal,
        Sub,
        Super
    };

    struct TextSize {
        int width;
        int height;
    };

    struct TextStyle {
        TextPositioning positioning;
        /// Has to be compatible with both pango markup and HTML/CSS.
        ///
        /// Empty if unspecified
        std::string weight;
        /// Has to be compatible with both pango markup and HTML/CSS.
        ///
        /// Empty if unspecified
        std::string size;
        Color color;

        bool specifies_color;

        TextStyle();
    };

    class TextSpan {

        std::variant<std::string, std::vector<TextSpan>> content;

        public:
        TextStyle style;
        // Emscripten doesn't currently support std::optional
        bool specifies_style;

        bool has_subspans() const;
        std::string& as_caption();
        const std::string& as_caption() const;
        std::vector<TextSpan>& as_subspans();
        const std::vector<TextSpan>& as_subspans() const;
        TextSpan();
        TextSpan(const std::string&);
        TextSpan(const std::vector<TextSpan>&);
    };

    struct Text {
        TextStyle style;
        std::vector<TextSpan> spans;
        graphene_point_t origin;
    };

    #ifndef __EMSCRIPTEN__
    cairo_t* cr;
    PangoLayout* pango_layout;

    private:
    std::string text_span_to_pango_markup(const TextSpan&, const std::optional<TextStyle>& parent_style = std::nullopt) const;
    public:
    /// Takes ownership of the pointers
    Renderer(cairo_t*,PangoLayout*);
    Renderer(const Renderer&) = delete;

    #else // __EMSCRIPTEN__ defined
    //       Lhasa-specific includes/definitions
    struct Path;
    struct DrawingCommand;
    struct BrushStyle {
        Color color;
        double line_width;
    };
    private:

    BrushStyle style;
    graphene_point_t position;
    std::vector<DrawingCommand> drawing_commands;
    emscripten::val text_measurement_function;
    // WIP
    std::unique_ptr<std::vector<DrawingCommand>> currently_created_path;
    // WIP
    std::vector<std::vector<DrawingCommand>*> drawing_structure_stack;

    public:

    struct Line {
        graphene_point_t start, end;
        BrushStyle style;
    };

    struct Arc {
        graphene_point_t origin;
        double radius, angle_one, angle_two;
        bool has_fill;
        Color fill_color;
        bool has_stroke;
        BrushStyle stroke_style;
    };

    struct Path {
        std::vector<DrawingCommand> commands;
        Color fill_color;
        bool has_fill;
        // this needs work. Do we need a boolean here?
        BrushStyle stroke_style;
    };

    struct DrawingCommand {
        std::variant<Line, Arc, Path, Text> content;

        bool is_path();
        bool is_arc();
        bool is_line();
        bool is_text();

        const Path& as_path() const;
        const Arc& as_arc() const;
        const Line& as_line() const;
        const Text& as_text() const;
    };

    Renderer(emscripten::val text_measurement_function);

    std::vector<DrawingCommand> get_commands() const;
    #endif

    // PRIMITIVES

    void move_to(double x, double y);
    void line_to(double x, double y);
    void arc(double x, double y, double radius, double angle_one, double angle_two);
    void fill();
    void stroke();
    void stroke_preserve();

    // PATH

    void new_path();
    void close_path();
    void new_sub_path();

    // STYLE

    void set_source_rgb(double r, double g, double b);
    void set_source_rgba(double r, double g, double b, double a);
    void set_line_width(double width);

    // TEXT
    TextSize measure_text(const Renderer::TextSpan&);
    void show_text(const Renderer::TextSpan&);

    ~Renderer();
};

/// Encapsulates routines and state to render CanvasMolecule
class MoleculeRenderContext {

    const CanvasMolecule& canvas_molecule;
    Renderer& ren;
    DisplayMode display_mode;
    float scale_factor;
    float x_offset;
    float y_offset;
    // Used to truncate bonds not to cover atoms
    std::map<unsigned int,graphene_rect_t> atom_idx_to_canvas_rect;


    void process_atom_highlight(const CanvasMolecule::Atom& atom);
    /// Returns the appendix span + info if the appendix is reversed
    std::tuple<Renderer::TextSpan, bool> process_appendix(const std::string& symbol, const std::optional<CanvasMolecule::Atom::Appendix>& appendix, const Renderer::TextStyle& inherited_style);
    /// Returns a pair of atom index and bonding rect
    std::pair<unsigned int,graphene_rect_t> render_atom(const CanvasMolecule::Atom& atom, DisplayMode render_mode = DisplayMode::Standard);
    // Returns on-screen bond coordinates
    // cropped not to overlap with atoms' symbols and appendices.
    // Accepts on-screen coordinates of bond atoms.
    std::pair<graphene_point_t,graphene_point_t> cropped_bond_coords(const graphene_point_t& first_atom, unsigned int first_atom_idx, const graphene_point_t& second_atom, unsigned int second_atom_idx);

    void draw_central_bond_line(const CanvasMolecule::Bond& bond);
    void draw_straight_wedge(const CanvasMolecule::Bond& bond, bool reversed);
    void draw_straight_dashed_bond(const CanvasMolecule::Bond& bond, bool reversed);
    void draw_wavy_bond(const CanvasMolecule::Bond& bond);
    void draw_side_bond_line(const CanvasMolecule::Bond& bond, bool addOrSub, std::optional<float> first_shortening_proportion, std::optional<float> second_shortening_proportion);
    void draw_centered_double_bond(const CanvasMolecule::Bond& bond);

    static const float CENTERED_DOUBLE_BOND_LINE_SEPARATION;
    static const float GEOMETRY_BOND_SPREAD_ANGLE;
    static const float WAVY_BOND_ARC_LENGTH;
    static const float GEOMETRY_BOND_DASH_SEPARATION;

    public:
    MoleculeRenderContext(const CanvasMolecule& cm, Renderer& ren, DisplayMode mode);
    ~MoleculeRenderContext();

    void draw_atoms();
    void draw_bonds();
};

} //coot::ligand_editor_canvas::impl

#endif // COOT_LIGAND_EDITOR_CANVAS_RENDER_HPP