
#include "molecules-container.hh"

std::pair<std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>,
          std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> >
molecules_container_t::mmrrcc_internal(const atom_selection_container_t &asc,
                                       const std::string &chain_id,
                                       unsigned int n_residue_per_residue_range,
                                       const clipper::Xmap<float> &xmap) const {

   bool exclude_NOC = true;
   float atom_mask_radius = 2.8;
   float NOC_mask_radius = 1.8;

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
      // std::cout << "mmrrcc:: all-atom-stats " << rs << " " << stats.correlation() << " from " << stats.n << " points ";
      // std::cout << std::endl;
   }
   for (it=residue_stats.second.begin(); it!=residue_stats.second.end(); ++it) {
      const coot::residue_spec_t &rs(it->first);
      const coot::util::density_correlation_stats_info_t &stats(it->second);
      // std::cout << "mmrrcc:: side-chain-stats " << rs << " " << stats.correlation() << " from " << stats.n << " points ";
      // std::cout << std::endl;
   }
   if (! residue_stats.second.empty()) {
      for (it=residue_stats.first.begin(); it!=residue_stats.first.end(); ++it) {
         const coot::residue_spec_t &rs_key(it->first);
         const coot::util::density_correlation_stats_info_t &stats_mc(it->second);
         std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>::const_iterator it_search =
            residue_stats.second.find(rs_key);
         if (it_search != residue_stats.second.end()) {
            const coot::util::density_correlation_stats_info_t &stats_sc(it_search->second);
            double delta = stats_sc.correlation() - stats_mc.correlation();
            if (false)
               std::cout << "   " << rs_key << " mc: " << stats_mc.correlation() << " sc: " << stats_sc.correlation()
                         << " sc-mc-delta: " << delta << std::endl;
         }
      }
   }

   return residue_stats;

}


std::pair<std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>,
          std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> >
molecules_container_t::get_mmrrcc(int imol, const std::string &chain_id, unsigned int n_residue_per_residue_range,
                                  int imol_map) const {

   std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> dummy;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_map)) {
         return mmrrcc_internal(molecules[imol].atom_sel, chain_id, n_residue_per_residue_range, molecules[imol_map].xmap);
      }
   }
   return std::make_pair(dummy, dummy);
}

std::pair<std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>,
          std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> >
molecules_container_t::mmrrcc(int imol, const std::string &chain_id,
                              int imol_map) const {

   unsigned int n_residue_per_residue_range = 11;
   std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> dummy;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_map)) {
         return mmrrcc_internal(molecules[imol].atom_sel, chain_id, n_residue_per_residue_range, molecules[imol_map].xmap);
      }
   }
   return std::make_pair(dummy, dummy);
}
