/*
 * src/make-a-dodec.hh
 *
 * Copyright 2020 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */
#ifndef MAKE_A_DODEC_HH
#define MAKE_A_DODEC_HH

#include <vector>
#include "generic-vertex.hh"
#include "coot-utils/g_triangle.hh"

// spikey mode is 0 for smooth shaded and
// 1 for spikey/flat shaded

enum { DODEC_SPIKEY_MODE = 1, DODEC_SMOOTH_MODE = 2 };

std::pair<std::vector<vn_vertex>, std::vector<g_triangle> >
make_pentakis_dodec(int spikey_mode);


#endif // MAKE_A_DODEC_HH
