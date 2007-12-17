# BL says we need to make interface.h for gtk1 and gtk2 
# if we want to use filechooser. We assume we have working interface.h
# for gtk1

cp interface.h interface.h.gtk1

echo "#if (GTK_MAJOR_VERSION == 1)" >interface.h
cat interface.h.gtk1 >> interface.h
echo "#endif // GTK1" >> interface.h 
echo "#if (GTK_MAJOR_VERSION == 2)" >> interface.h
./add-c-extern.sh interface.h.gtk2
cat interface.h.gtk2 >> interface.h
echo "#endif // GTK2" >> interface.h

sed -e 's/#include "callbacks.h.gtk2"/#if (GTK_MAJOR_VERSION > 1)\n\n#include "callbacks.h"/' \
    -e 's/interface.h.gtk2/interface.h/' \
    -e 's/support.h.gtk2/support.h/' \
    -e 's/#include <unistd.h>/#ifndef _MSC_VER\n#include <unistd.h>\n#endif/' \
    gtk2-interface.c > gtk2-interface.tmp
echo '#endif // (GTK_MAJOR_VERSION > 1)' >> gtk2-interface.tmp

mv gtk2-interface.tmp gtk2-interface.c 

