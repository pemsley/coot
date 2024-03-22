
#include "molecules-container.hh"

// This function does no normalisztion of the scales,
// presuming that they are pre-normalized.
// @return the index of the new map, or -1 on failure.
int
molecules_container_t::average_map(const std::string &imol_maps_string, std::vector<float> &scales) {

   int imol_new = -1;
   if (! scales.empty()) {
      std::vector<std::string> number_strings = coot::util::split_string(imol_maps_string, ":");
      std::vector<std::pair<int, float> > weights_and_maps;
      unsigned int scales_index = 0;
      for (const auto &item : number_strings) {
         int idx = coot::util::string_to_int(item);
         if (is_valid_map_molecule(idx)) {
            if (scales_index < scales.size()) {
               weights_and_maps.push_back(std::make_pair(idx, scales[scales_index]));
               scales_index++;
            }
         }
      }

      if (weights_and_maps.size() == scales.size()) {
         // now check that the grids are the same
         int idx = weights_and_maps[0].first;
         clipper::Grid_sampling gs_ref = molecules[idx].xmap.grid_sampling();
         unsigned int n_match = 0;
         for (unsigned int i=0; i<weights_and_maps.size(); i++) {
            idx = weights_and_maps[i].first;
            bool match = false;
            clipper::Grid_sampling gs_this = molecules[idx].xmap.grid_sampling();
            if (gs_ref.nu() == gs_this.nu())
               if (gs_ref.nv() == gs_this.nv())
                  if (gs_ref.nw() == gs_this.nw())
                     match = true;
            if (match) n_match++;
         }
         if (n_match == weights_and_maps.size()) {
            idx = weights_and_maps[0].first;
            std::string name = "Average map " + imol_maps_string;
            bool is_em_map = molecules[idx].is_EM_map();
            imol_new = molecules.size();
            std::vector<std::pair<clipper::Xmap<float>, float> > maps_and_weights(scales.size());
            // heavyweight stuff going on here
            for (unsigned int i=0; i<weights_and_maps.size(); i++) {
               idx = weights_and_maps[i].first;
               std::pair<clipper::Xmap<float>, float> p(molecules[idx].xmap, weights_and_maps[i].second);
               maps_and_weights[i] = p;
            }
            clipper::Xmap<float> xmap_new = coot::util::average_map(maps_and_weights);
            molecules.push_back(coot::molecule_t(name, imol_new, xmap_new, is_em_map));
         }
      }
   }
   return imol_new;
}

// This function does no normalisztion of the scales, presuming that they are pre-normalized.
// @return success status
bool
molecules_container_t::regen_map(int imol_map, const std::string &imol_maps_string, const std::vector<float> &scales) {

   bool status = false;
   // molecules.push_back(coot::molecule_t(name, imol_new, xmap_ref, is_em_map));

   if (is_valid_map_molecule(imol_map)) {
      if (! scales.empty()) {
         std::vector<std::string> number_strings = coot::util::split_string(imol_maps_string, ":");
         std::vector<std::pair<clipper::Xmap<float> *, float> > maps_and_scales_vec;
         unsigned int scales_index = 0;
         for (const auto &item : number_strings) {
            int idx = coot::util::string_to_int(item);
            if (is_valid_map_molecule(idx)) {
               if (scales_index < scales.size()) {
                  maps_and_scales_vec.push_back(std::make_pair(&molecules[idx].xmap, scales[scales_index]));
                  scales_index++;
               }
            }
         }

         if (maps_and_scales_vec.size() == scales.size()) {
            coot::util::regen_weighted_map(&molecules[imol_map].xmap, maps_and_scales_vec);
            status = true;
         }
      }
   }

   return status;
}
