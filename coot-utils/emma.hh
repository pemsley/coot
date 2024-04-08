/*
 * coot-utils/emma.hh
 *
 * Copyright 2014 by Medical Research Council
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

#include <clipper/contrib/edcalc.h>
#include "atom-selection-container.hh"

namespace coot {

   namespace util {

      class emma {
	 void sfs_from_boxed_molecule(mmdb::Manager *mol, float border);
	 double f(double r) const;
	 double f(const clipper::HKL_info::HKL_reference_index &hri_model,
		  const clipper::HKL_info::HKL_reference_index &hri_map,
		  double va) const;
      public:
	 emma(mmdb::Manager *mol, float border) {
	    sfs_from_boxed_molecule(mol, border);
	 }
	 clipper::Spacegroup spacegroup;
	 clipper::Cell cell;
	 clipper::Resolution reso;
	 clipper::HKL_info hkl_info;
	 clipper::HKL_data<clipper::data32::F_phi> fc_from_model;
	 void overlap(const clipper::Xmap<float> &xmap) const;
	 void overlap_simple(const clipper::Xmap<float> &xmap) const;
	 void test() const;

	 // put atoms on grid, calc sfs then average by invresolsq
	 std::vector<std::pair<double, double> > spherically_averaged_FT();

      };


      std::vector<std::pair<double, double> >
      spherically_averaged_molecule(const atom_selection_container_t &asc,
				    float angstroms_per_bin);
   }

}
