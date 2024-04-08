/*
 * coords/ramachandran-container.cc
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
 */


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
