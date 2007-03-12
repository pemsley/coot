
#if defined _MSC_VER
#include <windows.h>
#endif
#include "atom-utils.h"

// The reason we go throught these hoops is that we don't want general
// library code (i.e. directories other than this and src) to be
// dependent on src include files.
//
// We want to be able to create a graphics object from anywhere in the
// code (for debugging purposes), so that means we have to include
// graphics-info.h (which is in src). 
//
// so in atom-utils.h we only use library includes, this file is the
// only place we depend on src/ include.
// 
#include "graphics-info.h"

// Not really a gl thing, but there is no where else for it.   Consider
// a new file for such things.
// 
void asc_to_graphics(atom_selection_container_t asc, 
		     std::string label, 
		     bond_type draw_bonds_style, //  = ATOM_BONDS, default arg
		     float min_dist, // default arg
		     float max_dist) {  // default arg


   graphics_info_t g; 
   int imol = g.n_molecules; 
   
   std::cout <<  "initializing molecule " << imol << std::endl; 

   g.molecules[imol].initialize_coordinate_things_on_read_molecule(label); 
   g.molecules[imol].atom_sel = asc; 

   switch (draw_bonds_style) { 

   case ATOM_BONDS:
      if (min_dist == 0.0) { 
	 if (max_dist == 0.0) { 
	    g.molecules[imol].makebonds(); 
	 } else { 
	    g.molecules[imol].makebonds(min_dist, max_dist); 
	 } 
      } else { 
	 if (max_dist == 0.0) { 
	    g.molecules[imol].makebonds(min_dist); // which is max_dist at the recieving end
	 }
      }
      break; 
      
   case CA_BONDS:
      // we don't have multiple arguments functions for make_ca_bonds
      // 
      g.molecules[imol].make_ca_bonds(min_dist, max_dist); 
      break; 

   default:
      std::cout << "bond type " << draw_bonds_style
		<< " unknown in asc_to_graphics" << std::endl;
   }

   std::cout << "DEBUG: molecule " << imol << " with title: "
	<< g.molecules[imol].name_ << std::endl; 
   
   g.n_molecules++; 
} 



