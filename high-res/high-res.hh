/*
 * high-res/high-res.hh
 *
 * Copyright 2007 by Medical Research Council
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
 * General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 */


#include "mini-mol/mini-mol.hh"

namespace coot { 

   class high_res {
      minimol::molecule globular_molecule;
      std::pair<clipper::Coord_orth, mmdb::Manager *> get_middle_pos(const minimol::molecule &mol) const;
      void fill_globular_protein(const minimol::molecule &mol,
				 const clipper::Coord_orth &target,
				 mmdb::Manager *mmdb_mol);
      void fill_globular_protein_by_fragments(const minimol::molecule &mol,
					      const clipper::Coord_orth &target,
					      mmdb::Manager *mmdb_mol);
      void make_trees();
      void mark_neighbours(int iatom, int igroup,
			   const std::string &atom_name,
			   const std::vector<std::vector<int> > &neighbours,
			   mmdb::PPAtom atom_selection, int uddhandle);
      minimol::molecule
      filter_on_groups(const std::vector<std::vector<int> > &groups,
		       mmdb::Manager *mol,
		       mmdb::PPAtom atom_selection,
		       int n_selected_atom) const;

      static bool fragment_sorter(const minimol::fragment &a,
				  const minimol::fragment &b);

   public:
      high_res(const minimol::molecule &m);
      // use fill_globular_protein_by_fragments:
      high_res(const minimol::molecule &m, int iflag);
      high_res(int i);
      // This one is for Kevin, who wanted to pass the program a
      // centre, not let it calculate it itself.
      high_res(const minimol::molecule &m,
	       const clipper::Coord_orth &given_centre);
      void buccafilter_neighbours(); // slim down the multiple copies in space
			  // and leave only one averaged trace
      void buccafilter(); // slim down the multiple copies in space
			   // and leave only first fragment of
			   // overlapping fragments (fragments sorted
			   // by length first).
      void add_os();
      void add_cbetas();
      void output_pdb(const std::string &s) const;
   }; 

}
