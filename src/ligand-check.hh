

namespace coot { 
   class ligand_report_t {
   public:
      double mogul_percentile;
      double density_correlation_percentile;
      probe_clash_score_t pcs;
      ligand_report_t(double a, double b, probe_clash_score_t p) {
	 density_correlation_percentile = a;
	 mogul_percentile = b;
	 pcs = p;
      }
   };

}

