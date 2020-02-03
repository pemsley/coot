

namespace coot {
   namespace util {
      void sfcalc_genmap(mmdb::Manager *mol,
                         const clipper::HKL_data<clipper::data32::F_sigF> &fobs,
                         const clipper::HKL_data<clipper::data32::Flag> &free,
                         clipper::Xmap<float> *xmap_p);

   }
}
