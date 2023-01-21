

pub mod bond_modifer {
    pub enum Mode {
        Single,
        Double,
        Triple
    }
}

pub mod element_insertion {
    pub enum Element {
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
    }
}

pub mod structure_insertion {
    pub enum Structure {
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
    }
}


#[repr(C)]
// #[export_name = "ligand_editor_ActiveTool"]
pub enum ActiveTool {
    BondModifier(bond_modifer::Mode),
    StructureInsertion,
    ElementInsertion,
    /// Stereo out
    GeometryModifier,
    DeleteHydrogens,
    Delete,
    Format,
    ChargeModifier
}

