
#include "clipper/core/ramachandran.h"

class ramachandrans_container_t {

   clipper::Ramachandran rama;
   clipper::Ramachandran rama_gly;
   clipper::Ramachandran rama_pro;
   clipper::Ramachandran rama_non_gly_pro;
   clipper::Ramachandran rama_pre_pro;

public:
   // make a heavyweight object   
   ramachandrans_container_t() {
      rama.init(            clipper::Ramachandran::All5);
      rama_gly.init(        clipper::Ramachandran::Gly5);
      rama_pro.init(        clipper::Ramachandran::Pro5);
      rama_non_gly_pro.init(clipper::Ramachandran::NonGlyPro5);
   }
   
};

