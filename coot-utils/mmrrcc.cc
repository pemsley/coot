

#include <clipper/ccp4/ccp4_map_io.h>
#include <clipper/core/atomsf.h>

#include "utils/coot-utils.hh"
#include "geometry/protein-geometry.hh"
#include "coot-map-utils.hh"
#include "atom-selection-container.hh"

int main(int argc, char **argv) {

   bool is_cryo_em  = true;

   if (argc > 3) {
      std::string pdb_file_name = argv[1];
      std::string chain_id      = argv[2];
      std::string map_file_name = argv[3];

      std::cout << "Getting atoms... " << std::endl;
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, true, true, false);
      if (asc.read_success) {

         std::cout << "pdb read success " << pdb_file_name << std::endl;

         clipper::CCP4MAPfile file;
         clipper::Xmap<float> xmap;
         std::cout << "# reading map" << std::endl;
         file.open_read(map_file_name);
         file.import_xmap(xmap);
         file.close_read();

         if (is_cryo_em)
            clipper::ScatteringFactors::selectScattteringFactorsType(clipper::SF_ELECTRON);

         bool exclude_NOC = true;
         float atom_mask_radius = 2.8;
         float NOC_mask_radius = 1.8;
         unsigned int n_residue_per_residue_range = 11;
         std::pair<std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>,
                   std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> > residue_stats =
            coot::util::map_to_model_correlation_stats_per_residue_run(asc.mol, chain_id, xmap,
                                                                       n_residue_per_residue_range, exclude_NOC,
                                                                       atom_mask_radius, NOC_mask_radius);
         std::cout << "INFO:: We got " << residue_stats.first.size()  << " residue all-atom correlations"   << std::endl;
         std::cout << "INFO:: We got " << residue_stats.second.size() << " residue side-chain correlations" << std::endl;

         std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>::const_iterator it;
         for (it=residue_stats.first.begin(); it!=residue_stats.first.end(); ++it) {
            const coot::residue_spec_t &rs(it->first);
            const coot::util::density_correlation_stats_info_t &stats(it->second);
            std::cout << "   all-atom-stats " << rs << " " << stats.correlation() << " from " << stats.n << " points ";
            std::cout << std::endl;
         }
         for (it=residue_stats.second.begin(); it!=residue_stats.second.end(); ++it) {
            const coot::residue_spec_t &rs(it->first);
            const coot::util::density_correlation_stats_info_t &stats(it->second);
            std::cout << "   side-chain-stats " << rs << " " << stats.correlation() << " from " << stats.n << " points ";
            std::cout << std::endl;
         }

         // calculate the delta between CC for MC and CC for SC. It will be negative in regions of out of register

         if (! residue_stats.second.empty()) {
            for (it=residue_stats.first.begin(); it!=residue_stats.first.end(); ++it) {
               const coot::residue_spec_t &rs_key(it->first);
               const coot::util::density_correlation_stats_info_t &stats_mc(it->second);
               std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>::const_iterator it_search =
                  residue_stats.second.find(rs_key);
               if (it_search != residue_stats.second.end()) {
                  const coot::util::density_correlation_stats_info_t &stats_sc(it_search->second);
                  double delta = stats_sc.correlation() - stats_mc.correlation();
                  std::cout << "   " << rs_key << " mc: " << stats_mc.correlation() << " sc: " << stats_sc.correlation()
                            << " sc-mc-delta: " << delta << std::endl;
               }
            }
         }
         
      } else {
         std::cout << "Failed to read " << pdb_file_name << std::endl;
      }
   } else {
      std::cout << "Usage: " << argv[0] << " pdb_file_name chain_id map_file_name"
                << std::endl;
   }

   return 0;
}

