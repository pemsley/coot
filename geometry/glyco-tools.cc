/* geometry/protein-geometry.cc
 * 
 * Copyright 2003, 2004, 2005, 2006 The University of York
 * Author: Paul Emsley
 * Copyright 2007, 2008, 2009, 2010, 2011, 2012 The University of Oxford
 * Copyright 2014, 2015, 2016 by Medical Research Council
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <algorithm>  // needed for sort? Yes.
#include <stdexcept>  // Thow execption.

#include "utils/win-compat.hh"
#include "mini-mol/atom-quads.hh"
#include "geometry/protein-geometry.hh"
#include "utils/coot-utils.hh"

#include <sys/types.h> // for stating
#include <sys/stat.h>

#if !defined _MSC_VER
#include <unistd.h>
#else
#define DATADIR "C:/coot/share"
#define PKGDATADIR DATADIR
#define S_ISDIR(m)  (((m) & S_IFMT) == S_IFDIR)
#define S_ISREG(m)  (((m) & S_IFMT) == S_IFREG)
#endif

#include "clipper/core/clipper_util.h"

#include "compat/coot-sysdep.h"

#include "lbg-graph.hh"


// use auto-load if not present
void
coot::protein_geometry::use_unimodal_ring_torsion_restraints(int imol, const std::string &res_name,
							     int mmcif_read_number) {

   bool minimal = false;
   int idx = get_monomer_restraints_index(res_name, imol, minimal);
   if (idx == -1) {
      try_dynamic_add(res_name, mmcif_read_number);
      idx = get_monomer_restraints_index(res_name, imol, minimal);
   }

   if (idx != -1) {
      // continue

      // (first is the imol)
      std::vector <dict_torsion_restraint_t> &torsion_restraints =
	 dict_res_restraints[idx].second.torsion_restraint;

      // We don't want to clear them all, just the ring torsions
      // orig_torsion_restraints.clear();
      std::vector<std::string> ring_atom_names;
      ring_atom_names.push_back(" C1 "); ring_atom_names.push_back(" C2 ");
      ring_atom_names.push_back(" C3 "); ring_atom_names.push_back(" C4 ");
      ring_atom_names.push_back(" C5 "); ring_atom_names.push_back(" O5 ");
      if (res_name == "XYP")
	 for (unsigned int i=0; i<ring_atom_names.size(); i++)
	    ring_atom_names[i][3] = 'B';

      if (false) {
	 for (unsigned int i=0; i<torsion_restraints.size(); i++)
	    std::cout << " torsion restraint: " << i << " " << torsion_restraints[i] << std::endl;
	 std::cout << "...............  pre-delete size: " << torsion_restraints.size()
		   << " for " << res_name << std::endl;
      }

      torsion_restraints.erase(std::remove_if(torsion_restraints.begin(),
					      torsion_restraints.end(),
					      restraint_eraser(ring_atom_names)), torsion_restraints.end());

      if (false)
	 std::cout << "............... post-delete size: " << torsion_restraints.size()
		   << " for " << res_name << std::endl;

      std::vector<atom_name_torsion_quad> quads = get_reference_monomodal_torsion_quads(res_name);
      for (unsigned int i=0; i<quads.size(); i++) {
	 const atom_name_torsion_quad &quad = quads[i];
	 dict_torsion_restraint_t tors(quad.id,
				       quad.atom_name(0), quad.atom_name(1), quad.atom_name(2), quad.atom_name(3),
				       quad.torsion, 4.0, 1);
	 torsion_restraints.push_back(tors);
      }
   }
}

// pass the atom names and the desired torsion value - sigma is not specified
// by the user.
void
coot::protein_geometry::use_unimodal_ring_torsion_restraints(int imol, const std::string &res_name,
							     const std::vector<coot::atom_name_torsion_quad> &tors_info_vec,
							     int mmcif_read_number) {

   bool minimal = false;
   int idx = get_monomer_restraints_index(res_name, imol, minimal);
   if (idx == -1) {
      try_dynamic_add(res_name, mmcif_read_number);
      idx = get_monomer_restraints_index(res_name, imol, minimal);
   }

   if (idx != -1) {
      // continue

      // (first is the imol)
      std::vector <dict_torsion_restraint_t> &torsion_restraints =
	 dict_res_restraints[idx].second.torsion_restraint;

      // We don't want to clear all torsion restraints, just the ring torsions
      //
      std::set<std::string> ring_atom_names;
      // we generate the ring atoms names from the atoms in the passed quads
      for (std::size_t i=0; i<tors_info_vec.size(); i++) {
	 const atom_name_torsion_quad &q = tors_info_vec[i];
	 ring_atom_names.insert(q.atom_name(0));
	 ring_atom_names.insert(q.atom_name(1));
	 ring_atom_names.insert(q.atom_name(2));
	 ring_atom_names.insert(q.atom_name(3));
      }

      if (false) {
	 std::cout << "...............  pre-delete size: " << torsion_restraints.size()
		   << " for " << res_name << std::endl;
      }

      torsion_restraints.erase(std::remove_if(torsion_restraints.begin(),
					      torsion_restraints.end(),
					      restraint_eraser(ring_atom_names)), torsion_restraints.end());

      if (false)
	 std::cout << "............... post-delete size: " << torsion_restraints.size()
		   << " for " << res_name << std::endl;

      for (unsigned int i=0; i<tors_info_vec.size(); i++) {
	 const atom_name_torsion_quad &quad = tors_info_vec[i];
	 // this could have a nicer interface - using the quad here.
	 dict_torsion_restraint_t tors(quad.id,
				       quad.atom_name(0), quad.atom_name(1), quad.atom_name(2), quad.atom_name(3),
				       quad.torsion, 4.0, 1);
	 torsion_restraints.push_back(tors);
      }
   }
}


std::vector<coot::atom_name_torsion_quad>
coot::protein_geometry::get_reference_monomodal_torsion_quads(const std::string &res_name) const {

   std::vector<coot::atom_name_torsion_quad> v;
   if (res_name == "MAN" || res_name == "GLC" || res_name == "GAL") {
      v.push_back(coot::atom_name_torsion_quad("var_4",  "O5", "C5", "C4", "C3", -61.35));
      v.push_back(coot::atom_name_torsion_quad("var_2",  "C1", "O5", "C5", "C4",  67.67));
      v.push_back(coot::atom_name_torsion_quad("var_11", "C5", "O5", "C1", "C2", -67.59));
      v.push_back(coot::atom_name_torsion_quad("var_9",  "C3", "C2", "C1", "O5",  61.20));
      v.push_back(coot::atom_name_torsion_quad("var_6",  "C5", "C4", "C3", "C2",  53.75));
   }
   if (res_name == "BMA") {
      v.push_back(coot::atom_name_torsion_quad("var_2",  "C1", "O5", "C5", "C4",  67.57));
      v.push_back(coot::atom_name_torsion_quad("var_4",  "O5", "C5", "C4", "C3", -61.31));
      v.push_back(coot::atom_name_torsion_quad("var_6",  "C5", "C4", "C3", "C2",  53.80));
      v.push_back(coot::atom_name_torsion_quad("var_9",  "C3", "C2", "C1", "O5",  61.28));
      v.push_back(coot::atom_name_torsion_quad("var_11", "C5", "O5", "C1", "C2", -67.52));
   }
   if (res_name == "NAG") {
      v.push_back(coot::atom_name_torsion_quad("var_2",  "C1", "O5", "C5", "C4",  61.18));
      v.push_back(coot::atom_name_torsion_quad("var_5",  "O5", "C5", "C4", "C3", -57.63));
      v.push_back(coot::atom_name_torsion_quad("var_7",  "C5", "C4", "C3", "C2",  57.01));
      v.push_back(coot::atom_name_torsion_quad("var_10", "C3", "C2", "C1", "O5",  57.60));
      v.push_back(coot::atom_name_torsion_quad("var_1",  "C5", "O5", "C1", "C2", -61.19));
   }
   if (res_name == "FUC") {
      v.push_back(coot::atom_name_torsion_quad("var_2",  "C1", "O5", "C5", "C4", -67.61));
      v.push_back(coot::atom_name_torsion_quad("var_5",  "C5", "C4", "C3", "C2", -53.87));
      v.push_back(coot::atom_name_torsion_quad("var_11", "C5", "O5", "C1", "C2",  67.595));
      v.push_back(coot::atom_name_torsion_quad("var_9",  "C3", "C2", "C1", "O5", -61.26));
      v.push_back(coot::atom_name_torsion_quad("var_4",  "O5", "C5", "C4", "C3",  61.32));
   }
   if (res_name == "XYP") {
      v.push_back(coot::atom_name_torsion_quad("var_4",  "O5B", "C5B", "C4B", "C3B", -61.35));
      v.push_back(coot::atom_name_torsion_quad("var_2",  "C1B", "O5B", "C5B", "C4B",  67.67));
      v.push_back(coot::atom_name_torsion_quad("var_11", "C5B", "O5B", "C1B", "C2B", -67.59));
      v.push_back(coot::atom_name_torsion_quad("var_9",  "C3B", "C2B", "C1B", "O5B",  61.20));
      v.push_back(coot::atom_name_torsion_quad("var_6",  "C5B", "C4B", "C3B", "C2B",  53.75));
   }

   return v;
}

