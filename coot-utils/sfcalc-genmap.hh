
#ifndef SFCALC_GENMAP_HH
#define SFCALC_GENMAP_HH

#include <mmdb2/mmdb_manager.h>
#include <clipper/core/hkl_datatypes.h>
#include <clipper/core/xmap.h>

namespace coot {

   namespace util {

      class sfcalc_genmap_stats_t {
      public:
         class loc_table_t {
         public:
            class loc_table_item_t {
            public:
               float invresolsq;
               float scale;
               float lack_of_closure;
               loc_table_item_t(float i, float s, float loc) : invresolsq(i), scale(s), lack_of_closure(loc) {}
            };
            std::vector<loc_table_item_t> items;
            loc_table_t() {}
            void add(const loc_table_item_t &item) {
               items.push_back(item);
            }
            std::size_t size() const { return items.size(); }
         };
         float r_factor;
         float free_r_factor;
         float bulk_solvent_volume;
         float bulk_correction;
         unsigned int n_splines;
         loc_table_t loc_table;
         sfcalc_genmap_stats_t(float r_factor_in, float free_r_factor, float bulk_solvent_volume, float bulk_correction, unsigned int n_splines,
                               const loc_table_t &loc_table) : r_factor(r_factor_in), free_r_factor(free_r_factor), bulk_solvent_volume(bulk_solvent_volume),
                                                               bulk_correction(bulk_correction), n_splines(n_splines), loc_table(loc_table) {}
         sfcalc_genmap_stats_t() {
            r_factor = -1;
            free_r_factor = -1;
            bulk_solvent_volume = -1;
            bulk_correction = -1;
            n_splines = 0;
         }
      };

      void sfcalc_genmap(mmdb::Manager *mol,
                         const clipper::HKL_data<clipper::data32::F_sigF> &fobs,
                         const clipper::HKL_data<clipper::data32::Flag> &free,
                         clipper::Xmap<float> *xmap_p);

      sfcalc_genmap_stats_t sfcalc_genmaps_using_bulk_solvent(mmdb::Manager *mol,
                                                              const clipper::HKL_data<clipper::data32::F_sigF> &fobs,
                                                              const clipper::HKL_data<clipper::data32::Flag> &free,
                                                              const clipper::Cell &cell_for_fobs,
                                                              clipper::Xmap<float> *xmap_2fofc_p,
                                                              clipper::Xmap<float> *xmap_fofc_p);
   }
}

#endif // SFCALC_GENMAP_HH

