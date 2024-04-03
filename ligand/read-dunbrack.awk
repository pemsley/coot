 # ligand/read-dunbrack.awk
 #
 # Copyright 2007 by University of York
 # Author: Paul Emsley
 #
 # This file is part of Coot
 #
 # This program is free software; you can redistribute it and/or modify
 # it under the terms of the GNU Lesser General Public License as published
 # by the Free Software Foundation; either version 3 of the License, or (at
 # your option) any later version.
 #
 # This program is distributed in the hope that it will be useful, but
 # WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 # Lesser General Public License for more details.
 #
 # You should have received a copies of the GNU General Public License and
 # the GNU Lesser General Public License along with this program; if not,
 # write to the Free Software Foundation, Inc., 51 Franklin Street,
 # Fifth Floor, Boston, MA, 02110-1301, USA.
 # See http://www.gnu.org/licenses/
 #
 #

# usage: zcat bbind02.May.lib.gz | awk -f read-dunbrack.awk > chi-angles-autogen.cc
BEGIN { 
   print " ";
   print "#include \"chi-angles.hh\" ";
   print " ";
   print "void";
   print "coot::chi_angles::add_all_rotamers() { ";
   print " ";
}

/ARG|ASN|ASP|CYS|GLN|GLU|HIS|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL/ {
   restype = $1;
   r1 = $2;
   r2 = $3;
   r3 = $4;
   r4 = $5;
   n_r1 = $6;
   nr1234 = $7;
   p_r1234 = $8;
   sig_p_r1234 = $9;
   pr234_given_r1 = $10;
   sig_pr234_given_r1 = $11;
   chi1 = $12;
   sig_chi1 = $13;
   chi2 = $14;
   sig_chi2 = $15;
   chi3 = $16;
   sig_chi3 = $17;
   chi4 = $18;
   sig_chi4 = $19;

   if (restype == "PHE")
     if ((chi2+0) > 90.0)
       chi2 = chi2 - 180.0
   if (restype == "TYR")
     if ((chi2+0) > 90.0)
       chi2 = chi2 - 180.0

   if (sig_chi2 == "") { 
      chi2 = 0;
      sig_chi2 = -1;
   }
   if (sig_chi3 == "") { 
      chi3 = 0;
      sig_chi3 = -1;
   }
   if (sig_chi4 == "") { 
      chi4 = 0;
      sig_chi4 = -1;
   }
   print "      add_rotamer (std::string(" "\"" restype "\")" ", " r1 ", " r2 ", " r3 ", " r4 ", " n_r1 ", " nr1234 ", " p_r1234 ", " sig_p_r1234 ", " pr234_given_r1 ", " sig_pr234_given_r1 ", " chi1 ", " sig_chi1 ", " chi2 ", " sig_chi2 ", " chi3 ", " sig_chi3 ", " chi4 ", " sig_chi4 ");"
   # copy that MET rotamer to MSE too:
#    if (restype == "MET") 
#      print "      add_rotamer (std::string(" "\"" "MSE" "\")" ", " r1 ", " r2 ", " r3 ", " r4 ", " n_r1 ", " nr1234 ", " p_r1234 ", " sig_p_r1234 ", " pr234_given_r1 ", " sig_pr234_given_r1 ", " chi1 ", " sig_chi1 ", " chi2 ", " sig_chi2 ", " chi3 ", " sig_chi3 ", " chi4 ", " sig_chi4 ");"

      }

END {

   print " ";
   print "} ";

}
