

#include <vector>
#include <utility>

#include "geometry/residue-and-atom-specs.hh"

namespace coot {

   class cablam_markup_t {
      // don't make invalid markups - check
      // the atoms before constructing this object
   public:
      clipper::Coord_orth O_prev_pos;
      clipper::Coord_orth O_this_pos;
      clipper::Coord_orth O_next_pos;
      clipper::Coord_orth CA_proj_point_prev;
      clipper::Coord_orth CA_proj_point_this;
      clipper::Coord_orth CA_proj_point_next;
      double score;

      cablam_markup_t() { score = -1;}
      cablam_markup_t(mmdb::Atom *O_prev_at,
                      mmdb::Atom *O_this_at,
                      mmdb::Atom *O_next_at,
                      mmdb::Atom *CA_prev_at,
                      mmdb::Atom *CA_this_at,
                      mmdb::Atom *CA_next_at,
                      mmdb::Atom *CA_next_next_at);
   };

   std::vector<cablam_markup_t>
   make_cablam_markups(const std::vector<std::pair<residue_spec_t, double> > &residues,
                       mmdb::Manager *mol);

   // parse this log file and call the above function for each cablam outlier residue
   std::vector<cablam_markup_t>
   make_cablam_markups(mmdb::Manager *mol, const std::string &cablam_output_file_name);
}
