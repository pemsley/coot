
#include <optional>
namespace coot {

   // restraint types:
   //
   enum restraint_type_t {BOND_RESTRAINT=1, ANGLE_RESTRAINT=2, TORSION_RESTRAINT=4, PLANE_RESTRAINT=8,
                          NON_BONDED_CONTACT_RESTRAINT=16, CHIRAL_VOLUME_RESTRAINT=32, RAMACHANDRAN_RESTRAINT=64,
                          START_POS_RESTRAINT=128,
                          TARGET_POS_RESTRAINT=256, // restraint to make an atom be at a position
                          PARALLEL_PLANES_RESTRAINT=512,
                          GEMAN_MCCLURE_DISTANCE_RESTRAINT=1024,
                          TRANS_PEPTIDE_RESTRAINT=2048,
                          IMPROPER_DIHEDRAL_RESTRAINT=4096
   };

   class refinement_results_mini_stats_t {
   public:
      refinement_results_mini_stats_t() { is_set = false; nZ = std::nullopt; };
      // refinement_results_mini_stats_t(refinement_results_mini_stats_t &&) = default;
      // refinement_results_mini_stats_t(const refinement_results_mini_stats_t &) = default;
      refinement_results_mini_stats_t(restraint_type_t t, double dist, double target, double obs) :
         distortion(dist), target_value(target), observed_value(obs) {
            nZ = std::nullopt;
            is_set = true;
      }
      bool is_set;
      restraint_type_t type;
      double distortion;
      double target_value;
      double observed_value;
      std::optional<double> nZ;
   };

}
