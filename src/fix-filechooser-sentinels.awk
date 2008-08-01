

/ .*_filechooserdialog. = / { 
  gsub("NULL\);", "NULL, NULL);");
  print $0;
  next;
}

{ print $0;}

