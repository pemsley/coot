/* layla/ligand_editor_canvas/model.hpp
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

#ifndef COOT_LIGAND_EDITOR_CANVAS_MODEL_HPP
#define COOT_LIGAND_EDITOR_CANVAS_MODEL_HPP
#include <memory>
#include <vector>
#include <variant>
#include <optional>
#include <rdkit/GraphMol/RWMol.h>
// Forward declaration of types defined at "render.hpp"
namespace coot::ligand_editor_canvas::impl {
    struct Renderer;
    class MoleculeRenderContext;
}

#ifndef __EMSCRIPTEN__
#include <gtk/gtk.h>
#else // __EMSCRIPTEN__ defined
// Lhasa-specific includes/definitions
#include <graphene.h>
#endif

namespace coot {
namespace ligand_editor_canvas {

struct CurrentlyCreatedBond {
    float first_atom_x;
    float first_atom_y;
    //unsigned int first_atom_idx;
    float second_atom_x;
    float second_atom_y;
};

enum class FlipMode: unsigned char {
    /// Along the X axis
    Horizontal,
    /// Along the Y axis
    Vertical
};

enum class DisplayMode: unsigned char {
    Standard,
    AtomIndices,
    AtomNames
};

const char* display_mode_to_string(DisplayMode mode) noexcept;
std::optional<DisplayMode> display_mode_from_string(const char*) noexcept;

/// Drawing-friendly representation of RDKit molecule
class CanvasMolecule {
    // Rendering is done via a separate class 
    // for the sake of code organization
    friend class impl::MoleculeRenderContext;
    public:
    enum class AtomColor: unsigned char {
        /// Carbon and hydrogens
        Black,
        /// For Chlorine and Flourine
        Green,
        /// For Nitrogen
        Blue,
        /// For Oxygen
        Red,
        /// For Sulphur
        Brown,
        /// For Bromine
        DarkRed,
        /// For Phosphorus
        Orange,
        /// For Iodine
        DarkBlue
        // are there more colors?
    };
    
    struct Atom {
        std::string symbol;
        std::optional<std::string> name;
        
        /// Appendix represents optional elements that appear after atom's symbol,
        /// e.g. charge, hydrogens
        struct Appendix {
            /// Ionization
            int charge;
            /// For CH4, remainder would be "H4".
            ///
            /// todo: Should I just change it to an int
            /// representing the number of hydrogens?
            std::string superatoms;
            /// NH2 VS H2N, CH3 VS H3C, etc.
            bool reversed;
            Appendix() noexcept;
        };
        /// Appendix is set when we draw groups.
        /// For SO2-, the appendix would be "O2-"
        std::optional<Appendix> appendix;
        AtomColor color;
        /// Position on canvas (x axis)
        float x;
        /// Position on canvas (y axis)
        float y;
        /// Corresponds to RDKit atom index
        unsigned int idx;
        bool highlighted;
    };
    enum class BondType: unsigned char {
        Single,
        Double,
        Triple
    };
    enum class BondGeometry: unsigned char {
        Flat,
        DashedTowardsFirst,
        DashedTowardsSecond,
        WedgeTowardsFirst,
        WedgeTowardsSecond,
        /// Wavy
        Unspecified
    };

    enum class DoubleBondDrawingDirection: unsigned char {
        Primary,
        Secondary,
        Centered
    };

    struct Bond {
        BondType type;
        BondGeometry geometry;
        /// Set for double bonds. 
        /// It's the proportion of the bond's original length
        /// by which the parallel segment of the bond 
        /// has to be shortend from the "first" side.
        std::optional<float> first_shortening_proportion;
        /// Set for double bonds. 
        /// It's the proportion of the bond's original length
        /// by which the parallel segment of the bond 
        /// has to be shortend from the "second" side.
        std::optional<float> second_shortening_proportion;
        /// For double bonds
        std::optional<DoubleBondDrawingDirection> bond_drawing_direction;
        float first_atom_x;
        float first_atom_y;
        unsigned int first_atom_idx;
        float second_atom_x;
        float second_atom_y;
        unsigned int second_atom_idx;
        bool highlighted;

        /// Returns an [x,y] pair of numbers
        std::pair<float,float> get_perpendicular_versor() const noexcept;
        /// Returns an [x,y] pair of numbers
        std::pair<float,float> get_versor() const noexcept;
        /// Returns an [x,y] pair of numbers
        std::pair<float,float> get_vector() const noexcept;

        float get_length() const noexcept;
    };
    typedef std::variant<CanvasMolecule::Atom,CanvasMolecule::Bond> AtomOrBond;
    typedef std::optional<AtomOrBond> MaybeAtomOrBond;
    private:

    static const float BOND_DISTANCE_BOUNDARY;
    static const float ATOM_HITBOX_RADIUS;
    static const float BASE_SCALE_FACTOR;
    static const float BOND_LINE_SEPARATION;

    static BondType bond_type_from_rdkit(RDKit::Bond::BondType);
    static AtomColor atom_color_from_rdkit(const RDKit::Atom *) noexcept;
    static std::tuple<float,float,float> atom_color_to_rgb(AtomColor) noexcept;
    static std::string atom_color_to_html(AtomColor) noexcept;

    std::shared_ptr<RDKit::RWMol> rdkit_molecule;
    std::vector<Atom> atoms;
    std::vector<std::shared_ptr<Bond>> bonds;

    /// X offset due to translation.
    /// Has to be multiplied by scale to get on-screen coordinates
    float x_canvas_translation;
    /// Y offset due to translation.
    /// Has to be multiplied by scale to get on-screen coordinates
    float y_canvas_translation;

    /// Scale used by the widget
    float canvas_scale;

    /// The top-left and bottom-right points, in between which the molecule lies.
    /// The coordinates are in "RDKit space".
    /// They have to be multiplied by scale and added to the offsets to get on-screen coordinates
    std::pair<RDGeom::Point2D,RDGeom::Point2D> bounding_atom_coords;

    /// Coordinate map built in the previous call
    /// to `compute_molecule_geometry()`
    /// stored for reference to maintain alignment 
    /// when recomputing molecule geometry
    std::optional<RDGeom::INT_POINT2D_MAP> cached_atom_coordinate_map;

    /// Cached bond map, computed while lowering from RDKit.
    /// Associates atom indices with lists of bonds (the pointers are shared with the `bonds` vector).
    /// Used for various lookups: while drawing, in the lowering process itself, etc.
    std::map<unsigned int,std::vector<std::shared_ptr<Bond>>> bond_map;


    /// Computes the scale used for drawing
    /// And interfacing with screen coordinates
    float get_scale() const noexcept;

    /// Uses RDDepict to get molecule depiction & geometry info
    /// 
    /// Part of the lowering process.
    RDGeom::INT_POINT2D_MAP compute_molecule_geometry() const;

    /// Builds the drawing-friendly 2D molecule representation
    /// based on geometry computed by RDKit.
    ///
    /// Part of the lowering process.
    void build_internal_molecule_representation(const RDGeom::INT_POINT2D_MAP& coordinate_map);

    /// Iterates over rings and sets the right alignment for double bonds and atom symbols inside of rings
    ///
    /// Part of the lowering process.
    void process_alignment_in_rings();

    /// Computes length proportions by which 
    /// the parallel segments of double bonds have to be shortened
    /// to main aesthetic proportions.
    ///
    /// Part of the lowering process.
    void shorten_double_bonds();

    public:

    CanvasMolecule(std::shared_ptr<RDKit::RWMol> rdkit_mol);

    /// Replaces the inner shared_ptr to the molecule
    /// from which the CanvasMolecule is lowered.
    ///
    /// Meant to be called when performing a deep-copy
    void update_source_molecule(std::shared_ptr<RDKit::RWMol> rdkit_mol);

    /// Clears the drawing-friendly 2D representation data
    /// and re-creates it from the internal RDKit::RWMol
    /// 
    /// If `sanitize_after` is true, the molecule will get sanitized
    /// after lowering
    void lower_from_rdkit(bool sanitize_after);

    /// Clears `cached_atom_coordinate_map`, 
    /// forcing the subsequent call to `compute_molecule_geometry()`
    /// to determine the shape of the molecule from scratch.
    ///
    /// This is what makes the "Format" tool work.
    void clear_cached_atom_coordinate_map();

    /// Updates the `cached_atom_coordinate_map` after an atom has been removed
    /// in such a way as to prevent the cached molecule geometry from being broken
    void update_cached_atom_coordinate_map_after_atom_removal(unsigned int removed_atom_idx);
    
    /// Sets the scale for drawing
    void set_canvas_scale(float scale);

    void apply_canvas_translation(int delta_x, int delta_y) noexcept;
    std::pair<float,float> get_on_screen_coords(float x, float y) const noexcept;
    std::optional<std::pair<float,float>> get_on_screen_coords_of_atom(unsigned int atom_idx) const noexcept;
    graphene_rect_t get_on_screen_bounding_rect() const noexcept;
    void perform_flip(FlipMode flip_mode);
    void rotate_by_angle(double radians);

    /// Draws the molecule using the renderer
    void draw(impl::Renderer& ren, DisplayMode display_mode) const noexcept;

    /// Checks if any object matches the click coordinates passed as arguments.
    /// Returns the thing that was clicked on (or nullopt if there's no match).
    MaybeAtomOrBond resolve_click(int x, int y) const noexcept;

    void highlight_atom(int atom_idx);
    void highlight_bond(int atom_a, int atom_b); 
    void clear_highlights();

    static RDKit::Bond::BondType bond_type_to_rdkit(BondType) noexcept;
    static BondGeometry bond_geometry_from_rdkit(RDKit::Bond::BondDir) noexcept;
    static BondGeometry cycle_bond_geometry(BondGeometry) noexcept;
    static RDKit::Bond::BondDir bond_geometry_to_rdkit(BondGeometry) noexcept;

};


} // namespace ligand_editor_canvas
} // namesapce coot

#endif // COOT_LIGAND_EDITOR_CANVAS_MODEL_HPP