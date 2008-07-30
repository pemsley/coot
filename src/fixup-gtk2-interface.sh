
cp gtk2-interface.c gtk2-interface.c.copy
awk '
BEGIN { count = 0; }

count == 2 { count = 3; }
count == 1 { count = 2; }

$0 !~ "\"confirm_overwrite\"," { 
  print $0;
  if (count == 3) { 
    print "#endif /* (GTK_MINOR_VERSION > 9) */" ; 
    print ""; 
    count = 0; 
  }
} 

/\"confirm_overwrite\",/ { 
  print "#if (GTK_MINOR_VERSION > 9)" ; 
  print $0; 
  count = 1;
}' gtk2-interface.c.copy > gtk2-interface.c.copy-2

awk -f ../ifdef-for-create-aboutdialog.awk gtk2-interface.c.copy-2 > gtk2-interface.c

