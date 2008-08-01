

/ .*_.*chooserdialog.* = gtk_file_chooser_dialog_new / { 
  gsub("NULL\);", "NULL, NULL);");
  print $0;
  next;
}

{ print $0;}

