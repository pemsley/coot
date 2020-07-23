
#ifndef COORDS_RAMACHANDRAN_CONTAINER_HH
#define COORDS_RAMACHANDRAN_CONTAINER_HH

#include "clipper/core/ramachandran.h"

class ramachandrans_container_t {

public:
   // make a heavyweight object   
   ramachandrans_container_t() {
#ifdef CLIPPER_HAS_TOP8000
      rama.init(            clipper::Ramachandran::All2);
      rama_gly.init(        clipper::Ramachandran::Gly2);
      rama_pro.init(        clipper::Ramachandran::Pro2);
      rama_ileval.init(     clipper::Ramachandran::IleVal2);
      rama_pre_pro.init(    clipper::Ramachandran::PrePro2);
      rama_non_gly_pro_pre_pro_ileval.init(clipper::Ramachandran::NoGPIVpreP2);
      // simple approximation
      rama_non_gly_pro.init(clipper::Ramachandran::NoGPIVpreP2);
#else
      rama.init(            clipper::Ramachandran::All5);
      rama_gly.init(        clipper::Ramachandran::Gly5);
      rama_pro.init(        clipper::Ramachandran::Pro5);
      rama_non_gly_pro.init(clipper::Ramachandran::NonGlyPro5);
#endif
   }
   clipper::Ramachandran rama;
   clipper::Ramachandran rama_gly;
   clipper::Ramachandran rama_pro;
   clipper::Ramachandran rama_non_gly_pro;
   clipper::Ramachandran rama_pre_pro;
#ifdef CLIPPER_HAS_TOP8000
   clipper::Ramachandran rama_ileval;
   clipper::Ramachandran rama_non_gly_pro_pre_pro_ileval;
#endif
};

#endif // COORDS_RAMACHANDRAN_CONTAINER_HH
