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


/// For Edit Undo/Redo
class Operation {
    public:
    enum class OpType: unsigned char {
        ElementInsertion
    };
    class ElementInsertion {

    };
    private:
    union {
        ElementInsertion element_insertion;
    };
    OpType variant;
    public:
    
};

typedef std::vector<Operation> OperationStack;

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
    // todo: geometry support
    struct Bond {
        BondType type;
        float first_atom_x;
        float first_atom_y;
        unsigned int first_atom_idx;
        float second_atom_x;
        float second_atom_y;
        unsigned int second_atom_idx;
        bool highlighted;
    };

    private:

    typedef std::optional<std::variant<CanvasMolecule::Atom,CanvasMolecule::Bond>> MaybeAtomOrBond;
    
    /// Checks if any object matches the click coordinates passed as arguments.
    /// Returns the thing that was clicked on (or nullopt if there's no match).
    static MaybeAtomOrBond resolve_click(int x, int y);

    static BondType bond_type_from_rdkit(RDKit::Bond::BondType);

    std::shared_ptr<RDKit::RWMol> rdkit_molecule;
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;

    /// Clears the drawing-friendly 2D representation data
    /// and re-creates it from the internal RDKit::RWMol
    void lower_from_rdkit();

    public:

    CanvasMolecule(std::shared_ptr<RDKit::RWMol> rdkit_mol);

    void draw(GtkSnapshot* snapshot, PangoLayout* pango_layout, const graphene_rect_t *bounds) const noexcept;

};


} // namespace ligand_editor_canvas
} // namesapce coot

#endif // COOT_LIGAND_EDITOR_CANVAS_MODEL_HPP