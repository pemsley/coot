
#ifndef REFINEMENT_RESULTS_T_HH
#define REFINEMENT_RESULTS_T_HH

#include "refinement-lights.hh"

namespace coot {

   // ---------------------------------------------------------------
   // ---------------------------------------------------------------
   //     class refinement_results_t, helper class for sending text
   //     results back to invoking function.  Returned by minimize()
   //     function.
   // ---------------------------------------------------------------
   // ---------------------------------------------------------------
   
   class refinement_results_t { 
   public:
      bool found_restraints_flag; // if we found restraints or not.
      int progress; // GSL_CONTINUE, GSL_SUCCESS, GSL_ENOPROG (no progress)
      std::string info_text;
      std::vector<refinement_lights_info_t> lights;
      refinement_results_t(bool frf, int prog_in,
                           const std::vector<refinement_lights_info_t> &lights_in) {
                              found_restraints_flag = frf;
                              info_text = ""; // not used
                              progress = prog_in;
                              lights = lights_in;
                           }
      refinement_results_t(bool frf, int prog_in, const std::string &info_in) {
         found_restraints_flag = frf;
         info_text = info_in;
         progress = prog_in;
      }
      refinement_results_t() {
         info_text = "";
         found_restraints_flag = false;
      }
      refinement_results_t(const std::string &s_in) {
         info_text = s_in;
         found_restraints_flag = false;
         progress = -1; // unset
      }
      bool hooray() const;
   };

}

#endif // REFINEMENT_RESULTS_T_HH
