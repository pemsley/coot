
#ifdef USE_PYTHON
#include "Python.h"
#endif // USE_PYTHON
#include "graphics-info.h"
#include "c-interface.h" // is_valid_model_molecule()
#include "cc-interface.hh"
#include "coot-utils/coot-map-utils.hh"
#include "coot-utils/emma.hh"


#ifdef USE_GUILE
SCM spherical_density_overlap(SCM i_scm, SCM j_scm) {

   SCM r = SCM_BOOL_F;
   int imol = scm_to_int(i_scm); // map molecule
   int imol_map = scm_to_int(j_scm); // map molecule

   if (is_valid_model_molecule(imol)) { 
      if (is_valid_map_molecule(imol_map)) {

	 graphics_info_t g;

	 const clipper::Xmap<float> &m = g.molecules[imol_map].xmap;
	 clipper::Coord_orth c(0,0,0); // (set-rotation-centre -15 -4 21)
	 coot::util::map_fragment_info_t mf(m, c, 50, true);

	 if (mf.xmap.is_null()) {
	    std::cout << "null map fragment xmap " << std::endl;
	 } else { 
	    clipper::CCP4MAPfile mapout;
	    mapout.open_write("map-fragment-at-origin.map");
	    mapout.export_xmap(mf.xmap);
	    mapout.close_write();
	       
	    mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
	    coot::util::emma sphd(mol, 5); // 5 is border
	    sphd.overlap_simple(mf.xmap);
	    // sphd.overlap(mf.xmap);
	 }
      }
   }

   return r;
}
#endif // USE_GUILE

