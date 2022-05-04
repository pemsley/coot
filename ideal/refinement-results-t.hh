
#ifndef REFINEMENT_RESULTS_T_HH
#define REFINEMENT_RESULTS_T_HH

#include <vector>
#include "refinement-lights.hh"

#include <clipper/core/coords.h>
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

   // ---------------------------------------------------------------
   // ---------------------------------------------------------------
   //     class refinement_results_t, helper class for sending text
   //     results back to invoking function.  Returned by minimize()
   //     function.
   // ---------------------------------------------------------------
   // ---------------------------------------------------------------

   class refinement_results_for_rama_t {
   public:
      float ball_pos_x, ball_pos_y, ball_pos_z; // include file fun!
      atom_spec_t atom_spec_CA;
      float distortion;
      refinement_results_for_rama_t(mmdb::Atom *at_1,
                                    mmdb::Atom *at_2,
                                    mmdb::Atom *at_3,
                                    mmdb::Atom *at_4,
                                    mmdb::Atom *at_5,
                                    float distortion_in);
   };

   class refinement_results_nbc_baddie_t {
   public:
      atom_spec_t atom_spec_1;
      atom_spec_t atom_spec_2;
      clipper::Coord_orth atom_1_pos;
      clipper::Coord_orth atom_2_pos;
      clipper::Coord_orth mid_point;
      float score;
   };

   class refinement_results_t {
   public:

      bool found_restraints_flag; // if we found restraints or not.
      int progress; // GSL_CONTINUE, GSL_SUCCESS, GSL_ENOPROG (no progress)
      std::string info_text;
      int n_restraints;
      std::vector<refinement_lights_info_t> lights;
      bool refinement_results_contain_overall_nbc_score;
      bool refinement_results_contain_overall_rama_plot_score;
      float overall_nbc_score;
      // std::vector<std::pair<atom_spec_t, float> > sorted_nbc_baddies;
      std::vector<refinement_results_nbc_baddie_t> sorted_nbc_baddies;
      float overall_rama_plot_score;
      std::vector<std::pair<atom_spec_t, float> > sorted_rama_baddies;
      std::vector<refinement_results_for_rama_t> all_ramas;
      float overall_atom_pull_score;
      std::vector<std::pair<atom_spec_t, float> > sorted_atom_pulls; // all of them

      refinement_results_t(bool frf, int prog_in,
                           const std::vector<refinement_lights_info_t> &lights_in) {
         init();
         found_restraints_flag = frf;
         progress = prog_in;
         lights = lights_in;
     }
      refinement_results_t(bool frf, int prog_in, const std::string &info_in) {
         init();
         found_restraints_flag = frf;
         info_text = info_in;
         progress = prog_in;
      }
      refinement_results_t() {
         init();
      }
      explicit refinement_results_t(const std::string &s_in) {
         init();
         info_text = s_in;
         found_restraints_flag = false;
         progress = -1; // unset
      }
      void init() {
         info_text = "";
         found_restraints_flag = false;
         refinement_results_contain_overall_nbc_score = false;
         refinement_results_contain_overall_rama_plot_score = false;
         overall_nbc_score = 0.0;
         overall_rama_plot_score = 0.0;
         overall_atom_pull_score = 0.0;
         progress = 27; // GSL_ENOPROG
         n_restraints = -1; // unset
      }
      void show() const;
      std::pair<bool, std::string> hooray() const;
   };
}

#endif // REFINEMENT_RESULTS_T_HH
