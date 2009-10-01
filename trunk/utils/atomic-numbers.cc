/* coot-utils/atom-numbers.cc
 * 
 * Copyright 2004 by The University of York
 * Author: Paul Emsley
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

#include "coot-utils.hh"

std::vector<std::pair<std::string, int> > 
coot::util::atomic_number_atom_list() {

   std::vector<std::pair<std::string, int> > atom_list;
   std::pair<std::string, int> a;
   
   a.first = "H" ; a.second =  1;
   atom_list.push_back(a);
   a.first = "He" ; a.second =  2;
   atom_list.push_back(a);
   a.first = "Li" ; a.second =  3;
   atom_list.push_back(a);
   a.first = "Be" ; a.second =  4;
   atom_list.push_back(a);
   a.first = "B" ; a.second =  5;
   atom_list.push_back(a);
   a.first = "C" ; a.second =  6;
   atom_list.push_back(a);
   a.first = "N" ; a.second =  7;
   atom_list.push_back(a);
   a.first = "O" ; a.second =  8;
   atom_list.push_back(a);
   a.first = "F" ; a.second =  9;
   atom_list.push_back(a);
   a.first = "Ne" ; a.second =  10;
   atom_list.push_back(a);
   a.first = "Na" ; a.second =  11;
   atom_list.push_back(a);
   a.first = "Mg" ; a.second =  12;
   atom_list.push_back(a);
   a.first = "Al" ; a.second =  13;
   atom_list.push_back(a);
   a.first = "Si" ; a.second =  14;
   atom_list.push_back(a);
   a.first = "P" ; a.second =  15;
   atom_list.push_back(a);
   a.first = "S" ; a.second =  16;
   atom_list.push_back(a);
   a.first = "Cl" ; a.second =  17;
   atom_list.push_back(a);
   a.first = "Ar" ; a.second =  18;
   atom_list.push_back(a);
   a.first = "K" ; a.second =  19;
   atom_list.push_back(a);
   a.first = "Ca" ; a.second =  20;
   atom_list.push_back(a);
   a.first = "Sc" ; a.second =  21;
   atom_list.push_back(a);
   a.first = "Ti" ; a.second =  22;
   atom_list.push_back(a);
   a.first = "V" ; a.second =  23;
   atom_list.push_back(a);
   a.first = "Cr" ; a.second =  24;
   atom_list.push_back(a);
   a.first = "Mn" ; a.second =  25;
   atom_list.push_back(a);
   a.first = "Fe" ; a.second =  26;
   atom_list.push_back(a);
   a.first = "Co" ; a.second =  27;
   atom_list.push_back(a);
   a.first = "Ni" ; a.second =  28;
   atom_list.push_back(a);
   a.first = "Cu" ; a.second =  29;
   atom_list.push_back(a);
   a.first = "Zn" ; a.second =  30;
   atom_list.push_back(a);
   a.first = "Ga" ; a.second =  31;
   atom_list.push_back(a);
   a.first = "Ge" ; a.second =  32;
   atom_list.push_back(a);
   a.first = "As" ; a.second =  33;
   atom_list.push_back(a);
   a.first = "Se" ; a.second =  34;
   atom_list.push_back(a);
   a.first = "Br" ; a.second =  35;
   atom_list.push_back(a);
   a.first = "Kr" ; a.second =  36;
   atom_list.push_back(a);
   a.first = "Rb" ; a.second =  37;
   atom_list.push_back(a);
   a.first = "Sr" ; a.second =  38;
   atom_list.push_back(a);
   a.first = "Y" ; a.second =  39;
   atom_list.push_back(a);
   a.first = "Zr" ; a.second =  40;
   atom_list.push_back(a);
   a.first = "Nb" ; a.second =  41;
   atom_list.push_back(a);
   a.first = "Mo" ; a.second =  42;
   atom_list.push_back(a);
   a.first = "Tc" ; a.second =  43;
   atom_list.push_back(a);
   a.first = "Ru" ; a.second =  44;
   atom_list.push_back(a);
   a.first = "Rh" ; a.second =  45;
   atom_list.push_back(a);
   a.first = "Pd" ; a.second =  46;
   atom_list.push_back(a);
   a.first = "Ag" ; a.second =  47;
   atom_list.push_back(a);
   a.first = "Cd" ; a.second =  48;
   atom_list.push_back(a);
   a.first = "In" ; a.second =  49;
   atom_list.push_back(a);
   a.first = "Sn" ; a.second =  50;
   atom_list.push_back(a);
   a.first = "Sb" ; a.second =  51;
   atom_list.push_back(a);
   a.first = "Te" ; a.second =  52;
   atom_list.push_back(a);
   a.first = "I" ; a.second =  53;
   atom_list.push_back(a);
   a.first = "Xe" ; a.second =  54;
   atom_list.push_back(a);
   a.first = "Cs" ; a.second =  55;
   atom_list.push_back(a);
   a.first = "Ba" ; a.second =  56;
   atom_list.push_back(a);
   a.first = "La" ; a.second =  57;
   atom_list.push_back(a);
   a.first = "Ce" ; a.second =  58;
   atom_list.push_back(a);
   a.first = "Pr" ; a.second =  59;
   atom_list.push_back(a);
   a.first = "Nd" ; a.second =  60;
   atom_list.push_back(a);
   a.first = "Pm" ; a.second =  61;
   atom_list.push_back(a);
   a.first = "Sm" ; a.second =  62;
   atom_list.push_back(a);
   a.first = "Eu" ; a.second =  63;
   atom_list.push_back(a);
   a.first = "Gd" ; a.second =  64;
   atom_list.push_back(a);
   a.first = "Tb" ; a.second =  65;
   atom_list.push_back(a);
   a.first = "Dy" ; a.second =  66;
   atom_list.push_back(a);
   a.first = "Ho" ; a.second =  67;
   atom_list.push_back(a);
   a.first = "Er" ; a.second =  68;
   atom_list.push_back(a);
   a.first = "Tm" ; a.second =  69;
   atom_list.push_back(a);
   a.first = "Yb" ; a.second =  70;
   atom_list.push_back(a);
   a.first = "Lu" ; a.second =  71;
   atom_list.push_back(a);
   a.first = "Hf" ; a.second =  72;
   atom_list.push_back(a);
   a.first = "Ta" ; a.second =  73;
   atom_list.push_back(a);
   a.first = "W" ; a.second =  74;
   atom_list.push_back(a);
   a.first = "Re" ; a.second =  75;
   atom_list.push_back(a);
   a.first = "Os" ; a.second =  76;
   atom_list.push_back(a);
   a.first = "Ir" ; a.second =  77;
   atom_list.push_back(a);
   a.first = "Pt" ; a.second =  78;
   atom_list.push_back(a);
   a.first = "Au" ; a.second =  79;
   atom_list.push_back(a);
   a.first = "Hg" ; a.second =  80;
   atom_list.push_back(a);
   a.first = "Tl" ; a.second =  81;
   atom_list.push_back(a);
   a.first = "Pb" ; a.second =  82;
   atom_list.push_back(a);
   a.first = "Bi" ; a.second =  83;
   atom_list.push_back(a);
   a.first = "Po" ; a.second =  84;
   atom_list.push_back(a);
   a.first = "At" ; a.second =  85;
   atom_list.push_back(a);
   a.first = "Rn" ; a.second =  86;
   atom_list.push_back(a);
   a.first = "Fr" ; a.second =  87;
   atom_list.push_back(a);
   a.first = "Ra" ; a.second =  88;
   atom_list.push_back(a);
   a.first = "Ac" ; a.second =  89;
   atom_list.push_back(a);
   a.first = "Th" ; a.second =  90;
   atom_list.push_back(a);
   a.first = "Pa" ; a.second =  91;
   atom_list.push_back(a);
   a.first = "U" ; a.second =  92;
   atom_list.push_back(a);
   a.first = "Np" ; a.second =  93;
   atom_list.push_back(a);
   a.first = "Pu" ; a.second =  94;
   atom_list.push_back(a);
   a.first = "Am" ; a.second =  95;
   atom_list.push_back(a);
   a.first = "Cm" ; a.second =  96;
   atom_list.push_back(a);
   a.first = "Bk" ; a.second =  97;
   atom_list.push_back(a);
   a.first = "Cf" ; a.second =  98;
   atom_list.push_back(a);

   return atom_list;

}

// return -1 on not found
int
coot::util::atomic_number(const std::string &atom_name, 
			  const std::vector<std::pair<std::string, int> > &atom_list) {

   int atom_number = -1;
   std::string up_atom_name = coot::util::upcase(atom_name);
   if (up_atom_name.length() == 2)
      if (up_atom_name[0] == ' ')
	 up_atom_name = up_atom_name[1];

   for (unsigned int i=0; i<atom_list.size(); i++) {
      if (coot::util::upcase(atom_list[i].first) == up_atom_name) {
	 atom_number = atom_list[i].second;
	 break;
      }
   }
   return atom_number;
}
