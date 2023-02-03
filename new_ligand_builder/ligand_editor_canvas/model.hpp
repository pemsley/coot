#ifndef COOT_LIGAND_EDITOR_CANVAS_MODEL_HPP
#define COOT_LIGAND_EDITOR_CANVAS_MODEL_HPP

namespace coot {
namespace ligand_editor_canvas {




class BondModifier {
    public:
    enum class BondModifierMode {
        Single,
        Double,
        Triple
    };
    private:

    BondModifierMode mode;
    public:

};

class ElementInsertion {
    public:
    enum class Element {
        C,
        N,
        O,
        S,
        P,
        H,
        F,
        Cl,
        Br,
        I,
        // todo: what is this? an arbitrary element?
        // [like eg. Selene for organoselene compounds]
        X
    };

    private:
    Element element;
    public:

};



class StructureInsertion {
    public:

    enum class Structure: unsigned int {
        CycloPropaneRing,
        CycloButaneRing,
        CycloPentaneRing,
        CycloHexaneRing,
        BenzeneRing,
        CycloHeptaneRing,
        CycloOctaneRing,
        // todo:
        // "env residues"
        // "key"

    };
    private:
    Structure structure;

    public:

};





/// Drawing-friendly representation of RDKit molecule
class CanvasMolecule {

    public:

};


}
}

#endif // COOT_LIGAND_EDITOR_CANVAS_MODEL_HPP