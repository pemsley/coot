
 
#  Copyright 2007 by Paul Emsley
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or (at
#  your option) any later version.
 
#  This program is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
 
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
#  02110-1301, USA.


# run me with -F, 
#
# input file is Complete_rotamer_lib.csv e.g. in ~/rotamers on penelope

BEGIN { in_residue_table_flag = 0;
  print "#include \"chi-angles.hh\""
  print 
  print "void"
  print "coot::chi_angles::add_richardson_rotamers() {"
  print 
 }

{ 
  if (is_a_residue($2)) {
#     print $2 "is a resiude";
#     print "Residue", $2, "coming";
    this_residue=three_letter_code($2);
    in_residue_table_flag = 1;
    next; 
  }
}
 
in_residue_table_flag == 1 { 
   
#  print "Construct rot", this_residue;
#    for (i=1; i<=NF; i++) {
#       printf(":%s: ", $i);
#    }
#   printf("\n");

   seq_num = $1;
   tag = $2;
   name = $3;
   number = $4;
   percent = de_percent($5);
   alpha = de_percent($6);

#    print "de_percent on :" $6  ": gives :"alpha ":";

   beta = de_percent($7);
   other = de_percent($8);
   chi_1_mode = $10;
   chi_1_com  = $11;
   chi_2_mode = $13;
   chi_2_com  = $14;
   chi_3_mode = $16;
   chi_3_com  = $17;
   chi_4_mode = $19;
   chi_4_com  = $20;

   # print "chi_2_mode :" chi_2_mode ": chi_2_com :" chi_2_com ":"

#   printf("name: %s number: %s, percent: %s, alpha %s, beta: %s, other: %s, chi_1_mode: %s , chi_1_com: %s, chi_2_mode: %s, chi_2_com: %s\n", name, number, percent, alpha, beta, other, chi_1_mode, chi_1_com, chi_2_mode, chi_2_com);
   
   n_chis = get_n_chis(this_residue);
   # print "n_chis: ", n_chis 
   if (n_chis < 4) { 
     chi_4_mode = 0; 
     chi_4_com = 0;
   }
   if (n_chis < 3) { 
     chi_3_mode = 0; 
     chi_3_com = 0;
   }
   if (n_chis < 2) { 
     chi_2_mode = 0; 
     chi_2_com = 0;
   }

   # there are some additional so-called rotamers that don't have
   # chi_1_mode because they are not proper rotamers.
   if (chi_1_mode == "") {
     chi_1_mode = -5555;
   }
   if (chi_2_mode == "") {
     chi_2_mode = -5555;
   }
   if (chi_3_mode == "") {
     chi_3_mode = -5555;
   }
   if (chi_4_mode == "") {
     chi_4_mode = -5555;
   }

   if (length(name) >0) { 
     printf("       add_richardson_rotamer(\"%s\", \"%s\", %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);\n", this_residue, name, number, percent, alpha, beta, other, chi_1_mode, chi_1_com, chi_2_mode, chi_2_com, chi_3_mode, chi_3_com, chi_4_mode, chi_4_com);
     if (this_residue == "MET")
       printf("       add_richardson_rotamer(\"%s\", \"%s\", %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);\n", "MSE", name, number, percent, alpha, beta, other, chi_1_mode, chi_1_com, chi_2_mode, chi_2_com, chi_3_mode, chi_3_com, chi_4_mode, chi_4_com);
   }
}

NF == 0 { in_residue_table_flag = 0 }; 


function is_a_residue(res_str) { 

  if (res_str == "arginine") return 1; 
  if (res_str == "lysine") return 1; 
  if (res_str == "methionine") return 1; 
  if (res_str == "glutamate") return 1; 
  if (res_str == "glutamine") return 1; 
  if (res_str == "aspartate") return 1; 
  if (res_str == "asparagine") return 1; 
  if (res_str == "isoleucine") return 1; 
  if (res_str == "leucine") return 1; 
  if (res_str == "histidine") return 1; 
  if (res_str == "tryptophan") return 1; 
  if (res_str == "tyrosine") return 1; 
  if (res_str == "phenylalanine") return 1; 
  if (res_str == "proline") return 1; 
  if (res_str == "threonine") return 1; 
  if (res_str == "valine") return 1; 
  if (res_str == "serine") return 1; 
  if (res_str == "cysteine") 
    return 1; 
  else 
    return 0;
}

function get_n_chis(res_str) { 

  # print "get_n_chis for :" res_str ":"

  if (res_str == "ARG") return 4; 
  if (res_str == "LYS") return 4;
  if (res_str == "MET") return 3;
  if (res_str == "GLU") return 3; 
  if (res_str == "GLN") return 3; 
  if (res_str == "ASP") return 2; 
  if (res_str == "ASN") return 2; 
  if (res_str == "ILE") return 2; 
  if (res_str == "LEU") return 2; 
  if (res_str == "HIS") return 2; 
  if (res_str == "TRP") return 2; 
  if (res_str == "TYR") return 2; 
  if (res_str == "PHE") return 2; 
  if (res_str == "PRO") return 1;
  if (res_str == "THR") return 1; 
  if (res_str == "VAL") return 1; 
  if (res_str == "SER") return 1; 
  if (res_str == "CYS") 
    return 1; 
  else 
   return 0;
}

function three_letter_code(res_str) { 

  if (res_str == "arginine") return "ARG"; 
  if (res_str == "lysine")   return "LYS";
  if (res_str == "methionine") return "MET";
  if (res_str == "s-methionine") return "MSE";
  if (res_str == "glutamate") return "GLU"; 
  if (res_str == "glutamine") return "GLN"; 
  if (res_str == "aspartate") return "ASP"; 
  if (res_str == "asparagine") return "ASN"; 
  if (res_str == "isoleucine") return "ILE"; 
  if (res_str == "leucine") return "LEU"; 
  if (res_str == "histidine") return "HIS"; 
  if (res_str == "tryptophan") return "TRP"; 
  if (res_str == "tyrosine") return "TYR"; 
  if (res_str == "phenylalanine") return "PHE"; 
  if (res_str == "proline") return "PRO";
  if (res_str == "threonine") return "THR"; 
  if (res_str == "valine") return "VAL"; 
  if (res_str == "serine") return "SER"; 
  if (res_str == "cysteine") return "CYS";
  else 
    return "UNKNOWN";
}


function de_percent(str) { 
  
  l = length(str);
  if (l == 0) { 
    return str;
  } 
  last_char = substr(str, l);
  r = str;
  if (last_char == "%") {
    # print "returning chopped :" substr(str, 0, (l-1)) ":"
    r = substr(str, 0, (l-1));
  } else {
    # print "last_char was not % it was :" last_char ":"
  }
  return r;
}

END {
  print "}"
  print
}
