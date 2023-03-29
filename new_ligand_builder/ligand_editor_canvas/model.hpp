#ifndef COOT_LIGAND_EDITOR_CANVAS_MODEL_HPP
#define COOT_LIGAND_EDITOR_CANVAS_MODEL_HPP
#include <cstddef>
#include <gtk/gtk.h>
#include <memory>
#include <vector>
#include <variant>
#include <optional>
#include <rdkit/GraphMol/RWMol.h>

namespace coot {
namespace ligand_editor_canvas {


/// Drawing-friendly representation of RDKit molecule
class CanvasMolecule {
    public:
    enum class AtomColor: unsigned char {
        /// Carbon and hydrogens
        Black,
        /// For Chlorine
        Green,
        /// For Nitrogen
        Blue,
        /// For Oxygen
        Red
        // are there more colors?
    };
    
    struct Atom {
        std::string symbol;
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
        WedgeTowardsSecond
    };
    // todo: geometry support
    struct Bond {
        BondType type;
        BondGeometry geometry;
        float first_atom_x;
        float first_atom_y;
        unsigned int first_atom_idx;
        float second_atom_x;
        float second_atom_y;
        unsigned int second_atom_idx;
        bool highlighted;

        /// Returns an [x,y] pair of numbers
        std::pair<float,float> get_perpendicular_versor() const noexcept;
    };
    typedef std::variant<CanvasMolecule::Atom,CanvasMolecule::Bond> AtomOrBond;
    typedef std::optional<AtomOrBond> MaybeAtomOrBond;
    private:

    static const float ATOM_HITBOX_RADIUS;
    static const float BOND_DISTANCE_BOUNDARY;
    static const float BASE_SCALE_FACTOR;
    static const float BOND_LINE_SEPARATION;

    static BondType bond_type_from_rdkit(RDKit::Bond::BondType);
    static AtomColor atom_color_from_rdkit(const RDKit::Atom *) noexcept;
    static std::tuple<float,float,float> atom_color_to_rgb(AtomColor) noexcept;
    static std::string atom_color_to_html(AtomColor) noexcept;

    std::shared_ptr<RDKit::RWMol> rdkit_molecule;
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;

    /// X centering offset on canvas
    float x_canvas_size_adjustment;
    /// Y centering offset on canvas
    float y_canvas_size_adjustment;

    /// X offset due to translation
    int x_canvas_translation;
    /// Y offset due to translation
    int y_canvas_translation;

    /// Scale used by the widget
    float canvas_scale;


    /// Computes the scale used for drawing
    /// And interfacing with screen coordinates
    float get_scale() const noexcept;

    public:

    CanvasMolecule(std::shared_ptr<RDKit::RWMol> rdkit_mol);

    /// Replaces the inner shared_ptr to the molecule
    /// from which the CanvasMolecule is lowered.
    ///
    /// Meant to be called when performing a deep-copy
    void update_source_molecule(std::shared_ptr<RDKit::RWMol> rdkit_mol);

    /// Clears the drawing-friendly 2D representation data
    /// and re-creates it from the internal RDKit::RWMol
    void lower_from_rdkit();

    /// Changes the relative placement of the molecule on the screen
    void set_canvas_size_adjustment_from_bounds(const graphene_rect_t *bounds) noexcept;
    
    /// Sets the scale for drawing
    void set_canvas_scale(float scale);

    void apply_canvas_translation(int delta_x, int delta_y) noexcept;

    /// Draws the molecule on the widget.
    /// Be sure to call set_offset_from_bounds before drawing.
    void draw(GtkSnapshot* snapshot, PangoLayout* pango_layout, const graphene_rect_t *bounds) const noexcept;

    /// Checks if any object matches the click coordinates passed as arguments.
    /// Returns the thing that was clicked on (or nullopt if there's no match).
    MaybeAtomOrBond resolve_click(int x, int y) const noexcept;

    void highlight_atom(int atom_idx);
    void highlight_bond(int atom_a, int atom_b); 
    void clear_highlights();

    static RDKit::Bond::BondType bond_type_to_rdkit(BondType) noexcept;

};


} // namespace ligand_editor_canvas
} // namesapce coot

#endif // COOT_LIGAND_EDITOR_CANVAS_MODEL_HPP