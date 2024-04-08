/*
 * mini-mol/mol-hierarchy.cc
 *
 * Copyright 2018 by Medical Research Council
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


	 for(unsigned int ifrag=0; ifrag<target_ca_coords.fragments.size(); ifrag++) {
	    for(int ires=target_ca_coords[ifrag].min_res_no(); ires<=target_ca_coords[ifrag].max_residue_number(); ires++) {
	       for (unsigned int iat=0; iat<target_ca_coords[ifrag][ires].atoms.size(); iat++) {
		  std::cout << " " << target_ca_coords[ifrag].fragment_id << " " << ires << " " << target_ca_coords[ifrag][ires]
			    << " " << target_ca_coords[ifrag][ires][iat].name
			    << " " << target_ca_coords[ifrag][ires][iat].pos.format() << std::endl;
	       }
	    }
	 }
