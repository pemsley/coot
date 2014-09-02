BEGIN {

  n_lines = 0;
  if (ARGC > 1) { 
    file_name = ARGV[1]
    new_file_name = make_fixed_format_file_name(file_name);
    if (file_exists(file_name)) { 
      while (getline lines[n_lines++] < file_name) {
      }
      process_lines(lines, n_lines, new_file_name)
    } else { 
      print "file not found:",file_name
    }
  } else {
    print "usage:", ARGV[0],"<cif-filename>"
  }
}

function file_exists(file_name) { 

  s = "test -f " file
  return (! system(s))
  
}

function make_fixed_format_file_name(file_name) { 
  
  l = length(file_name);
  stub = substr(file_name, 0, l-4);
  r = stub "-fixed-format.cif"
  return r;
} 

function process_lines(lines, n_lines, new_file_name) { 

  print "output file:", new_file_name;
  ok_to_print = 1;
  for (i=0; i<n_lines; i++) {
    if (lines[i] == "loop_")
      ok_to_print = check_printable(lines, n_lines, i);
    if (ok_to_print) 
      print lines[i] > new_file_name;
  }
} 

# return a boolean
# 
# look forward in lines from iline, count the number of lines not
# beginning with "_" or "loop_"
# 
# If n_data_lines > 0, then OK to print
# if 0 then no (return 0)
#
function check_printable(lines, n_lines, iline) { 

  for (ii=(iline+1); ii<n_lines; ii++) {
    # print "....... checking :" lines[ii] ":"
    if (lines[ii] == "loop_")
      return 0;
    if (substr(lines[ii], 0, 1) != "_") 
      return 1;
  }
  return 0;
}
