// -*-c++-*-

#ifndef ATOM_UTILS_HH
#define ATOM_UTILS_HH

#include <string>
#include <vector>

#include "Cartesian.h"

#include "mmdb_manager.h"
#include "mmdb-extras.h"
#include "mmdb-crystal.h"
#include "mmdb.h"

#include "Bond_lines.h" // needs mmdb-extras.h, mmdb-manager.h <string>, <vector>, mmdb-crystal.h


// atom_selection_container_t is in mmdb-extras.h


enum bond_type {ATOM_BONDS, CA_BONDS}; 


// 
void asc_to_graphics(atom_selection_container_t asc, 
		     std::string label, 
		     bond_type draw_bonds_style = ATOM_BONDS, 
		     float min_dist = 0.0,
		     float max_dist = 0.0); 


#endif // ATOM_UTILS_HH
