
# use a better awk if we can find it
if [ -d $HOME/awk/bin ] ; then
   export PATH=$HOME/awk/bin:$PATH
fi

if [ -z "$1" ] ; then
   echo error in fixup-gtk2-interface.sh
   exit 2
fi

cp $1 $1.copy


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
}' $1.copy > gtk2-interface.c.copy-2

awk -f ../ifdef-for-create-aboutdialog.awk gtk2-interface.c.copy-2 > gtk2-interface.c.copy-3
awk -f  gtk2-interface-rot-trans-fixup.awk gtk2-interface.c.copy-3 > gtk2-interface.c.copy-4
# echo awk for sentinels
awk -f fix-filechooser-sentinels.awk gtk2-interface.c.copy-4 > gtk2-interface.c


