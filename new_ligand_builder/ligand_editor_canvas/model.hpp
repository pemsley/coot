#ifndef COOT_LIGAND_EDITOR_CANVAS_MODEL_HPP
#define COOT_LIGAND_EDITOR_CANVAS_MODEL_HPP
#include <cstddef>
#include <gtk/gtk.h>
#include <memory>
#include <vector>
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
        int x;
        /// Position on canvas (y axis)
        int y;
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
        std::size_t first_atom_idx;
        std::size_t second_atom_idx;
        bool highlighted;
    };

    private:
    std::shared_ptr<RDKit::RWMol> rdkit_molecule;
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    public:

    void draw(GtkSnapshot* snapshot) const noexcept;

};


} // namespace ligand_editor_canvas
} // namesapce coot

#endif // COOT_LIGAND_EDITOR_CANVAS_MODEL_HPP