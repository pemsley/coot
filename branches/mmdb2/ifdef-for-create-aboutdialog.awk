BEGIN { 
  marker = 0;
  s = "#if (((GTK_MAJOR_VERSION == 2) && (GTK_MINOR_VERSION > 5)) || GTK_MAJOR_VERSION > 2)";
}

/GtkWidget*/ { 
  if (marker == 2) { 
    marker = 3;
    print s;
  }
}

$0 == "}" { 

  if (marker == 1) 
    marker = 2;
  if (marker == 3) 
    marker = 4;
}

/create_sft_dialog .void./ { marker=1; } 



{ print $0;
  if (marker == 4) { 
    marker = 0;
    print "#endif /* GTK2 version */";
  }
}




  


