
#usage: awk -f angle-torsion-autogen.awk batched.dat > AngleInfo-angle-torsions-autogen.cc
#
# batched.data can be found at ~emsley/angles/batched.dat
#

BEGIN {
  print "";
  print "#include \"AngleInfo.h\"";
  print " ";
  print " ";
}

NF == 3 { 
  if ( $1 == "0" ) { 
    # a new function:
    print " ";
    print "void";
    torsion = $2 + 0;
    torsion = -torsion;  # kludge to fit to clipper convention
    v = torsion + 365;
    print "AngleInfo::from_batched_angle_torsions_bits_" v "() {";
    print " ";
  }

  print "     assign_angle_torsion("$1 "," torsion "," $3 ");";

  if ( $1 == "175" ) { 
    print "}";
  }
}

END {
}
