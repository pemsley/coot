#ifndef MAKE_A_DODEC_HH
#define MAKE_A_DODEC_HH

#include <vector>
#include "generic-vertex.hh"
#include "g_triangle.hh"

// spikey mode is 0 for smooth shaded and
// 1 for spikey/flat shaded

enum { DODEC_SPIKEY_MODE = 1, DODEC_SMOOTH_MODE = 2 };

std::pair<std::vector<vn_vertex>, std::vector<g_triangle> >
make_pentakis_dodec(int spikey_mode);


#endif // MAKE_A_DODEC_HH
