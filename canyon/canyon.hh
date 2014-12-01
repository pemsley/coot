
#include <mmdb2/mmdb_manager.h>
#include <clipper/core/coords.h>

namespace coot {

   class trace_info_t {
   public:
      std::vector<clipper::Coord_orth> points;
      clipper::Array2d<double> beta_vec_hats;
      trace_info_t(const std::vector<clipper::Coord_orth> &points_in,
		   const clipper::Array2d<double> &beta_vec_hats_in) {
	 points = points_in;
	 beta_vec_hats = beta_vec_hats_in;
      }
      trace_info_t() {} // for dummy
   };

   class t_matrix : public clipper::Matrix<double> {
   public:
      t_matrix() {}
      t_matrix(unsigned int i, unsigned int j) : clipper::Matrix<double>(i,j) {
      }
      // The "special case" where XTX is a 3x3 matrix:
      // (using this makes it trivial to invert the result)
      clipper::Mat33<double> XTX() const {

	 clipper::Mat33<double> r;

	 for (unsigned int i=0; i<3; i++) {
	    for (unsigned int j=0; j<3; j++) {
	       double sum = 0;
	       for (unsigned int ii=0; ii<rows(); ii++) {
		  sum += (*this)(ii,i) * (*this)(ii,j);
	       }
	       r(j,i) = sum;
	    }
	 }
	 return r;
      }
   };

   class canyon {
      mmdb::Manager *mol;
      clipper::Coord_orth start_point;
      clipper::Coord_orth start_path_point;
      clipper::Coord_orth end_point;
      clipper::Coord_orth end_path_point; // into canyon from end point
      void update_atom_selection(const clipper::Coord_orth &pt);
      int selHnd; // get updated internally.
      clipper::Coord_orth add_probe_points(int iround, double frac,
					   clipper::Coord_orth &pt,
					   const trace_info_t &t1,
					   bool use_trace_info_for_path_uv,
					   bool make_surface_points);
      void output_grid() const;
      clipper::Coord_orth get_path_uv(const clipper::Coord_orth &start_path_uv,
				      const clipper::Coord_orth &end_path_uv, double frac) const;
      clipper::Coord_orth get_path_uv(double frac) const;
      clipper::Coord_orth get_path_uv(const trace_info_t &ti, double frac) const;
      void write_fit_path(const std::vector<clipper::Coord_orth> &probe_path_points,
			  const clipper::Array2d<double> &beta_vec_hats,
			  const std::string &file_name) const;

      enum { N_CANYON_STEPS = 120,
	     N_THETA_STEPS  = 50 };

      // These go together, for every N_CANYON_STEPS step, there is a
      // set of points in the plane and (probe_path_points) the origin
      // of the ring of points
      // 
      std::pair<bool, clipper::Coord_orth> surface_points[N_CANYON_STEPS][N_THETA_STEPS];
      clipper::Coord_orth  probe_path_points[N_CANYON_STEPS];

      clipper::Array2d<double>
      polynomial_path_fit(const std::vector<clipper::Coord_orth> &points) const;
      std::vector<double>
      polynomial_path_fit(const std::vector<clipper::Coord_orth> &points,
	 unsigned int points_vec_index) const;
      trace_info_t get_initial_trace(); 
      trace_info_t get_refined_trace(const trace_info_t &ti); 
      public:
      canyon(mmdb::Manager *mol);
      void set_start_point(     const clipper::Coord_orth &pt) { start_point      = pt; }
      void set_start_path_point(const clipper::Coord_orth &pt) { start_path_point = pt; }
      void set_end_point(       const clipper::Coord_orth &pt) { end_point        = pt; }
      void set_end_path_point(  const clipper::Coord_orth &pt) { end_path_point   = pt; }
      double tetrahedron_volume(const std::vector<clipper::Coord_orth> &positions) const;
      void trace();
      double trace_volume();
      
   };
}

