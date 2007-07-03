
# run me with -F, 
#
# input file is Complete_rotamer_lib.csv in ~/rotamers on penelope

BEGIN { in_residue_table_flag = 0; }

{ 
  if (is_a_residue($2)) {
#     print $2 "is a resiude";
#     print "Residue", $2, "coming";
    this_residue=$2;
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

#   printf("name: %s number: %s, percent: %s, alpha %s, beta: %s, other: %s, chi_1_mode: %s , chi_1_com: %s, chi_2_mode: %s, chi_2_com: %s\n", name, number, percent, alpha, beta, other, chi_1_mode, chi_1_com, chi_2_mode, chi_2_com);
   
   n_chis = get_n_chis(this_residue);
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

  if (res_str == "arginine") return 4; 
  if (res_str == "lysine")   return 4;
  if (res_str == "methionine") return 3;
  if (res_str == "glutamate") return 3; 
  if (res_str == "glutamine") return 3; 
  if (res_str == "aspartate") return 2; 
  if (res_str == "asparagine") return 2; 
  if (res_str == "isoleucine") return 1; 
  if (res_str == "leucine") return 2; 
  if (res_str == "histidine") return 2; 
  if (res_str == "tryptophan") return 2; 
  if (res_str == "tyrosine") return 2; 
  if (res_str == "phenylalanine") return 2; 
  if (res_str == "proline") return 1;
  if (res_str == "threonine") return 1; 
  if (res_str == "valine") return 1; 
  if (res_str == "serine") return 1; 
  if (res_str == "cysteine") 
    return 1; 
  else 
   return 0;
}

function de_percent(str) { 
  
  l = length(str);
  if (l == 0) { 
    return str;
  } 
  last_char = substr(str, l);
  r = str;
  if (last_char == "%") {
    r = substr(str, 0, l);
  }
  return r;

}