
#include "mmdb_manager.h"
#include "clipper/core/coords.h"

#include "coot-utils.hh"
#include "protein-geometry.hh"

namespace coot {

   class hole {
      CMMDBManager *mol;
      int radius_handle; // set by assign_vdw_radii;
      clipper::Coord_orth from_pt;
      clipper::Coord_orth to_pt;
      clipper::Coord_orth v_hat;
      
      void make_atom_selection(int selhnd, const clipper::Coord_orth &pt,
			       double radius_prev) const;
      std::pair<clipper::Coord_orth, double>
      optimize_point(const clipper::Coord_orth &pt, int selhnd);
      double sphere_size(const clipper::Coord_orth &pt, int selhnd) const;
      void assign_vdw_radii(const coot::protein_geometry &geom);
      void write_probe_path(const std::vector<std::pair<clipper::Coord_orth, double> > &probe_path) const;
      void write_probe_path_using_spheres(const std::vector<std::pair<clipper::Coord_orth, double> > &probe_path) const;
      std::pair<double, double> x_y(const double &r, const double &theta) const;
      colour_holder probe_size_to_colour(double ps) const {
	 colour_holder ch;
	 return ch;
      } 
      
   public:
      hole(CMMDBManager *mol, int selHnd,
	   const clipper::Coord_orth &from_pt,
	   const clipper::Coord_orth &to_pt,
	   const coot::protein_geometry &geom_in);
      
      void generate(); 
   };
} 
