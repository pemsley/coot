
BEGIN { old_mat_no = -1;}

/^MTRIX/ { mat_no=$2+0; r1=$3; r2=$4; r3=$5; t1=$6;
  m[mat_no, mcount[mat_no]+0, 0] = r1;
  m[mat_no, mcount[mat_no]+0, 1] = r2;
  m[mat_no, mcount[mat_no]+0, 2] = r3;
  m[mat_no, mcount[mat_no]+0, 3] = t1;
  mcount[mat_no]++;
}

END { 
  print "(define imol 0)";
  # the chains are labelled "1", "2", "3", "4"

  n = 0;
  for (ichain=1; ichain<=4; ichain++) {
    incs_chain=1;
    for (imat_no=1; imat_no<=mat_no; imat_no++) { 
      # print "imat_no:",imat_no;
      print "(add-strict-ncs-matrix imol", "\"" ichain "\"","\"" ichain "-" incs_chain "\"";
      for (ncount=0; ncount<3; ncount++) { 
	  # print "ncount",ncount;
	  print "  ", m[imat_no, ncount, 0], m[imat_no, ncount, 1], m[imat_no, ncount, 2];
      }
      print "  ", m[imat_no, 0, 3], m[imat_no, 1, 3], m[imat_no, 2, 3], ")";
      incs_chain++;
    }
  }
}
