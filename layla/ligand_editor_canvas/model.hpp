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

typedef std::map<unsigned int, std::string> SmilesMap;
typedef std::map<unsigned int, std::string> InchiKeyMap;

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

    typedef unsigned char highlight_t;
    enum class HighlightType: highlight_t {
        Hover = 1,
        Edition = 2,
        Error = 4,
        // A concept for the future
        Selection = 8
    };

   struct Atom {
        std::string symbol;
        std::optional<std::string> name;

        /// Appendix represents optional elements that appear after atom's symbol (i.e. superatoms),
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
            bool vertical;
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
        /// Highlight bitmask
        highlight_t highlight;
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
        /// Highlight bitmask
        highlight_t highlight;

        /// Returns an [x,y] pair of numbers
        std::pair<float,float> get_perpendicular_versor() const noexcept;
        /// Returns an [x,y] pair of numbers
        std::pair<float,float> get_versor() const noexcept;
        /// Returns an [x,y] pair of numbers
        std::pair<float,float> get_vector() const noexcept;

        float get_length() const noexcept;
    };

    struct QEDInfo {
        unsigned int number_of_hydrogen_bond_acceptors;
        unsigned int number_of_hydrogen_bond_donors;
        unsigned int number_of_rotatable_bonds;
        unsigned int number_of_aromatic_rings;
        unsigned int number_of_alerts;
        double molecular_weight;
        double alogp;         /// Hydrophobicity
        double molecular_polar_surface_area;
        double ads_mw;
        double ads_alogp;
        double ads_hba;
        double ads_hbd;
        double ads_psa;
        double ads_rotb;
        double ads_arom;
        double ads_alert;
        double qed_score;
    };

    typedef std::variant<CanvasMolecule::Atom,CanvasMolecule::Bond> AtomOrBond;
    typedef std::optional<AtomOrBond> MaybeAtomOrBond;
    static const float BASE_SCALE_FACTOR;

    private:

    static const float BOND_DISTANCE_BOUNDARY;
    static const float ATOM_HITBOX_RADIUS;
    static const float BOND_LINE_SEPARATION;
    static const float VERTICAL_SUPERATOM_ANGLE_THRESHOLD;

    static BondType bond_type_from_rdkit(RDKit::Bond::BondType);
    static AtomColor atom_color_from_rdkit(const RDKit::Atom *) noexcept;
    static std::tuple<float,float,float> atom_color_to_rgb(AtomColor) noexcept;
    static std::tuple<float,float,float> hightlight_to_rgb(HighlightType) noexcept;
    static std::string atom_color_to_html(AtomColor) noexcept;
    static std::optional<HighlightType> determine_dominant_highlight(highlight_t) noexcept;

    std::shared_ptr<RDKit::RWMol> rdkit_molecule;
    std::vector<Atom> atoms;
    std::vector<std::shared_ptr<Bond>> bonds;

    /// X offset due to translation.
    /// Has to be multiplied by scale and viewport offset must be added to get on-screen coordinates
    float x_canvas_translation;
    /// Y offset due to translation.
    /// Has to be multiplied by scale and viewport offset must be added to get on-screen coordinates
    float y_canvas_translation;

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

    /// QED info is updated while lowering from RDKit.
    std::optional<QEDInfo> qed_info;

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

    /// Updates the `qed_info` variable
    void update_qed_info();
    
    /// Manages error highlights
    /// Part of the lowering process.
    void process_problematic_areas(bool allow_invalid_molecules);

    public:

    CanvasMolecule(std::shared_ptr<RDKit::RWMol> rdkit_mol, bool allow_invalid_mol);

    /// Replaces the inner shared_ptr to the molecule
    /// from which the CanvasMolecule is lowered.
    ///
    /// Meant to be called when performing a deep-copy
    void update_source_molecule(std::shared_ptr<RDKit::RWMol> rdkit_mol);

    /// Clears the drawing-friendly 2D representation data
    /// and re-creates it from the internal RDKit::RWMol
    ///
    /// If `sanitize_after` is true, the molecule will get sanitized
    /// after lowering.
    /// QED gets recomputed and updated if `with_qed` is true
    void lower_from_rdkit(bool sanitize_after, bool with_qed = true);

    /// Clears `cached_atom_coordinate_map`,
    /// forcing the subsequent call to `compute_molecule_geometry()`
    /// to determine the shape of the molecule from scratch.
    ///
    /// This is what makes the "Format" tool work.
    void clear_cached_atom_coordinate_map();

    /// Updates the `cached_atom_coordinate_map` after an atom has been removed
    /// in such a way as to prevent the cached molecule geometry from being broken
    void update_cached_atom_coordinate_map_after_atom_removal(unsigned int removed_atom_idx);

    void apply_canvas_translation(int delta_x, int delta_y, float scale) noexcept;
    std::pair<float,float> get_on_screen_coords(float x, float y, const std::pair<int, int>& viewport_offset, float scale) const noexcept;
    std::optional<std::pair<float,float>> get_on_screen_coords_of_atom(unsigned int atom_idx, const std::pair<int, int>& viewport_offset, float scale) const noexcept;
    graphene_rect_t get_on_screen_bounding_rect(const std::pair<int, int>& viewport_offset, float scale) const noexcept;
    std::optional<QEDInfo> get_qed_info() const noexcept;
    void perform_flip(FlipMode flip_mode);
    void rotate_by_angle(double radians);

    /// Draws the molecule using the renderer
    void draw(impl::Renderer& ren, DisplayMode display_mode, const std::pair<int, int>& viewport_offset, float scale) const noexcept;

    /// Checks if any object matches the click coordinates passed as arguments.
    /// Returns the thing that was clicked on (or nullopt if there's no match).
    /// Works on canvas coordinates, i.e. those acquired after viewport offset had been accounted for
    MaybeAtomOrBond resolve_click(int x, int y, float scale) const noexcept;

    void add_atom_highlight(int atom_idx, HighlightType htype);
    void add_bond_highlight(unsigned int atom_a, unsigned int atom_b, HighlightType htype);
    void add_highlight_to_all_bonds(HighlightType htype);
    /// Clears the highlight flag of the given type (for both atoms and bonds)
    void clear_highlights(HighlightType htype = HighlightType::Hover);

    static RDKit::Bond::BondType bond_type_to_rdkit(BondType) noexcept;
    static BondGeometry bond_geometry_from_rdkit(RDKit::Bond::BondDir) noexcept;
    static BondGeometry cycle_bond_geometry(BondGeometry) noexcept;
    static RDKit::Bond::BondDir bond_geometry_to_rdkit(BondGeometry) noexcept;

};


} // namespace ligand_editor_canvas
} // namesapce coot

#endif // COOT_LIGAND_EDITOR_CANVAS_MODEL_HPP
