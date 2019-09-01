
// mmdb atom usage includes:
#include <string>
#include <vector>
#include <map> 

#include <mmdb2/mmdb_manager.h>

#include "coords/Cartesian.h"
#include "coords/mmdb-extras.h"

#include "graphical_skel.h"

// 
graphical_bonds_container
GraphicalSkel::make_graphical_bonds(const clipper::Xmap<float> &map,
				    const clipper::Xmap<int>   &l1,
				    coot::Cartesian centre_point,
				    float box_radius,
				    float cut_off) const {
      
   std::vector<std::vector<graphics_line_t> > cp_vec; // build them up in
					         // vector, convert on
					         // return

   
   clipper::Skeleton_fast<int,float>::Neighbours neigh( map );

   int toplevel = 0, level;
   clipper::Xmap_base::Map_reference_index iy;
   for ( iy = l1.first(); !iy.last(); iy.next() )
     toplevel = clipper::Util::max( toplevel, l1[iy] );

   // Get the limits of the box in which we want bones (i.e. grid)
   clipper::Coord_orth centre(centre_point.get_x(),
			      centre_point.get_y(),
			      centre_point.get_z());
   
   clipper::Coord_frac centre_f = centre.coord_frac(map.cell());

   clipper::Coord_frac box0(
    centre_f.u() - box_radius/map.cell().descr().a(),
    centre_f.v() - box_radius/map.cell().descr().b(),
    centre_f.w() - box_radius/map.cell().descr().c() );
  clipper::Coord_frac box1(
    centre_f.u() + box_radius/map.cell().descr().a(),
    centre_f.v() + box_radius/map.cell().descr().b(),
    centre_f.w() + box_radius/map.cell().descr().c() );

  clipper::Grid_map grid( box0.coord_grid(map.grid_sampling()),
			  box1.coord_grid(map.grid_sampling()));
  
  clipper::Coord_grid grid_bit(1,1,1);
   
  cp_vec.resize( toplevel + 1);
  clipper::Xmap_base::Map_reference_coord ix( map, grid.min() ), iu, iv, iw;
  for ( iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() )
    for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() )
      for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {

	   if (l1[iw] > 0) {
	      if (map[iw] > cut_off) {

		 clipper::Coord_grid gb0 = iw.coord() - grid_bit;
		 clipper::Coord_grid gb1 = iw.coord() + grid_bit;

		 clipper::Xmap_base::Map_reference_coord iyc = iw;

		 for ( int in = 0; in < neigh.size(); in++ ) {
		   iyc.set_coord( iw.coord() + neigh[in] );

		   if (l1[iyc] > 0) {
		       
		     if (map[iyc] > cut_off) { 
		       clipper::Coord_grid c_g_1 = iw.coord();
		       clipper::Coord_grid c_g_2 = iyc.coord(); 
			  
		       clipper::Coord_frac c_f_1 =
			 c_g_1.coord_frac(l1.grid_sampling());
		       clipper::Coord_frac c_f_2 =
			 c_g_2.coord_frac(l1.grid_sampling());
			  
		       clipper::Coord_orth c_o_1 = c_f_1.coord_orth(l1.cell());
		       clipper::Coord_orth c_o_2 = c_f_2.coord_orth(l1.cell()); 
			  
		       coot::Cartesian f(c_o_1.x(), c_o_1.y(), c_o_1.z());
		       coot::Cartesian s(c_o_2.x(), c_o_2.y(), c_o_2.z());
			  
		       coot::CartesianPair line(f, s); 
			  
		       level = clipper::Util::max( l1[iyc], l1[iw] ); // -1
		       graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
		       cp_vec[level].push_back(graphics_line_t(line, cc, false, false, -1, -1, -1));
		       //cp_vec[0].push_back(line);
		     }
		   }
		 }
	      }
	   }
      }
   
  // now convert the vector back to an allocated pointer
  // 
  std::cout << "skeletonization toplevel: " << toplevel << " \n";
  graphical_bonds_container gbonds;
  for ( int i = 1; i <= toplevel; i++ ) {
     // in this box, that is.
     // cout << "At skeleton level " << i << " there are " << cp_vec[i].size() << " points \n";
     gbonds.add_colour(cp_vec[i]);  
  }
  return gbonds;
}


// Old Style bonds in the clipper asymmetric unit
// 
graphical_bonds_container 
GraphicalSkel::make_graphical_bonds( const clipper::Xmap<float> &map,
				     const clipper::Xmap<int>   &l1 ) const {

   graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
   std::vector<graphics_line_t> cp_vec; 

   int n_lines = 0;
   float cut_off = 0.15;
   
   clipper::Skeleton_fast<int,float>::Neighbours neigh( map );


   clipper::Xmap_base::Map_reference_index ix;
   for (ix = l1.first(); !ix.last(); ix.next() )  { // iterator index.

      if (l1[ix] == 1) { 

	 if (map[ix] > cut_off) { 
	    
	    for (int ii=0; ii<neigh.size(); ii++) {
	       
	       clipper::Coord_grid c_g_2 = ix.coord() + neigh[ii];
	       
	       if (l1.get_data(c_g_2) == 1) {

		  if (map.get_data(c_g_2) > cut_off) { 

		     clipper::Coord_grid c_g_1 = ix.coord();
		  			
		     clipper::Coord_frac c_f_1 =
			c_g_1.coord_frac(l1.grid_sampling());
		     clipper::Coord_frac c_f_2 =
			c_g_2.coord_frac(l1.grid_sampling());
		  
		     clipper::Coord_orth c_o_1 = c_f_1.coord_orth(l1.cell());
		     clipper::Coord_orth c_o_2 = c_f_2.coord_orth(l1.cell()); 
		  
		     coot::Cartesian f(c_o_1.x(), c_o_1.y(), c_o_1.z());
		     coot::Cartesian s(c_o_2.x(), c_o_2.y(), c_o_2.z());
		  
		     coot::CartesianPair line(f, s); 

		     cp_vec.push_back(graphics_line_t(line, cc, false, false, -1, -1, -1));
		  }
	       }
	    }
	 }
      }
   }

   // now convert the vector back to an allocated pointer
   graphical_bonds_container a(cp_vec);
   return a;   
}




// we return the max level of (tip->core) of the skeleton. 
// 
int
GraphicalSkel::Pprune(const clipper::Xmap<float> &map,
		      clipper::Xmap<int> *l1,
		      float cut_off) {

   // Recall, that after Skeleton, l1 is a map that contains 0s
   // (representing non-skeleton) and 1 (representing skeleton).
   //
   // First, invert that, so that the skeleton is at -1 (we will
   // reinvert later).
   // 
   // The idea is that we progressively chop the tips/end-points so we
   // get some indication of points that are deeper in the skeleton
   // (therefore likely to be mainchain skeleton rather than
   // sidechain).  The more "distance" from a tip, the higher the
   // number of that gridpoint.
   // 
   // Note that loops don't have tips, so that they will have to be
   // treated differently.  (btw I believe that James' "prun" routine
   // does something similar).
   //
   // We manipulate the whole of the ASU.
   //

   // map_density_distribution(*l1);
   // map_format(*l1, 40); 

   // First, convert 1s to -1s.
   // 
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = l1->first(); !ix.last(); ix.next() )  { 
      if ( (*l1)[ix] == 1) { 
	 (*l1)[ix] = -1; 
      }
   }

   //
   int n_skelled = -1;
   int level = 0; 


   // convert this to a do loop, so that we don't have to pre-assign
   // n_skelled?
   // 
   while (n_skelled != 0) {
      level++; 
      n_skelled = Ptip_convert(map, l1, level, cut_off);
      std::cout << "n_skelled: " << n_skelled << " at level " << level << std::endl;
   }

   // now convert the stuff that didn't already get moved
   // (should be only loops (things without tips)).
   //


   for (ix = l1->first(); !ix.last(); ix.next() )  {
      if ( (*l1)[ix] == -1) { 
	 (*l1)[ix] = level; 
      }
   }

   return level;
} 


// Return the number of tips of the input skel map
//
int 
GraphicalSkel::N_tips(const clipper::Xmap<float> &map, 
		      const clipper::Xmap<int> &l1,
		      float cut_off) const { 


   int n_tip = 0; 

   int n_neigbs; 

   clipper::Skeleton_fast<int,float>::Neighbours neigh( map );
   clipper::Coord_grid c_g; 

   clipper::Xmap_base::Map_reference_index ix;
   for (ix = map.first(); !ix.last(); ix.next() )  {

      if ( l1[ix] > 0 ) { 

	 if ( map[ix] > cut_off ) { 
	    
	    n_neigbs = 0; 

	    for (int ii=0; ii<neigh.size(); ii++) {
	       
	       c_g = ix.coord() + neigh[ii];
	       
	       if (l1.get_data(c_g) > 0) {

		  if (map.get_data(c_g) > cut_off) { 

		     n_neigbs++;
		  }
	       }
	    }
	    
	    if (n_neigbs <= 2) { 

	       // was a tip then

	       n_tip++; 
	    }
	 } 
      }
   }

   return n_tip; 

}


// Now in the style of make_graphical_bonds (using Neighbours instead of a box)
//
// Return the number of points converted at this level
//  
int 
GraphicalSkel::Ptip_convert(const clipper::Xmap<float> &map,
			    clipper::Xmap<int> *l1, int level,
			    float cut_off) {

   clipper::Coord_grid c_g_2; 
   clipper::Skeleton_fast<int,float>::Neighbours neigh( map );
   int n_skel_neighbs = 0; 
   int n_converted = 0; 


   clipper::Xmap_base::Map_reference_index ix;
   for (ix = l1->first(); !ix.last(); ix.next() )  {

      if ( (*l1)[ix] == -1 ) { 
	 
	 if (map[ix] > cut_off) { 

	    n_skel_neighbs = 0; 
	    
	    for (int ii=0; ii<neigh.size(); ii++) {
	       
	       c_g_2 = ix.coord() + neigh[ii];
	       
	       if (l1->get_data(c_g_2) == -1) {
		  
		  if (map.get_data(c_g_2) > cut_off) {
		     
		     n_skel_neighbs++;
		     
		  }
	       }
	    }
      
	    if (n_skel_neighbs == 1) { 
	       
	       (*l1)[ix] = level;
	       n_converted++; 
	    }
	 }
      }
   }

   return n_converted; 

}

int 
GraphicalSkel::Ptip_convert_old(const clipper::Xmap<float> &map,
			    clipper::Xmap<int> *l1, int level,
			    float cut_off) {

   //
   int n_skel_neighbs;
   int n_converted = 0;
   int n_count_all = 0; 
   
   int n_hist[30];

   for (int ii=0; ii<30; ii++)
      n_hist[ii] = 0; 


   clipper::Coord_grid grid_bit(1,1,1);
   
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = l1->first(); !ix.last(); ix.next() )  {

      if ( (*l1)[ix] == -1) {

	 if ( map[ix] > cut_off ) { 

	    n_skel_neighbs = 0; 

	    //
	    clipper::Coord_grid g0 = ix.coord() - grid_bit;
	    clipper::Coord_grid g1 = ix.coord() + grid_bit;
      
	    // inner looping over a small box of density bounded by g0 and g1
	    // 
	    clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw;
      
	    i0 = clipper::Xmap_base::Map_reference_coord( *l1, g0 );
      
	    for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() ) {
	       for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() ) {
		  for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() ) {

		     n_count_all++;
		     
		     if ( (*l1)[iw] == -1 ) {

			if (map[iw] > cut_off) { 

			   n_skel_neighbs++;
			}
		     }
		  }
	       }
	    }

	    n_hist[n_skel_neighbs]++;
	    
	    if (n_skel_neighbs <= 2) { // was a tip (only one neighbour 
	       // plus itself).
	    
	       (*l1)[ix] = level;
	       // (*l1)[ix] = 0;
	       n_converted++; 
	    }
		  
	 }
      
      }
   }

//    for (int ii=0; ii<30; ii++) { 
//       cout << "skel neighbour hist: " << ii<< "   " <<  n_hist[ii] << endl;
//    }
   
   return n_converted; 

} 

