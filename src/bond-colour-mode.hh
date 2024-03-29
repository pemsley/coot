#ifndef BOND_COLOUR_MODE_HH
#define BOND_COLOUR_MODE_HH


// This is a copy from molecule_class_info_t. They should be
// consistent, or better yet, include this file into
// molecule-class-info.h

namespace coot {
   enum { UNSET_TYPE = -1, NORMAL_BONDS=1, CA_BONDS=2,
	  COLOUR_BY_CHAIN_BONDS=3,
	  CA_BONDS_PLUS_LIGANDS=4, BONDS_NO_WATERS=5, BONDS_SEC_STRUCT_COLOUR=6,
	  BONDS_NO_HYDROGENS=15,
	  CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR=7,
	  CA_BONDS_PLUS_LIGANDS_B_FACTOR_COLOUR=14,
	  CA_BONDS_PLUS_LIGANDS_AND_SIDECHAINS=17,
	  COLOUR_BY_MOLECULE_BONDS=8,
	  COLOUR_BY_RAINBOW_BONDS=9,
	  COLOUR_BY_B_FACTOR_BONDS=10,
	  COLOUR_BY_OCCUPANCY_BONDS=11,
	  COLOUR_BY_USER_DEFINED_COLOURS____BONDS=12,
	  COLOUR_BY_USER_DEFINED_COLOURS_CA_BONDS=13 };
}


#endif // BOND_COLOUR_MODE_HH
