
sed -e 's/#include "callbacks.h.gtk2"/#if (GTK_MAJOR_VERSION > 1)\n\n#include "callbacks.h"/' \
    -e 's/interface.h.gtk2/interface.h/' \
    -e 's/support.h.gtk2/support.h/' \
    -e 's/#include <unistd.h>/#ifndef _MSC_VER\n#include <unistd.h>\n#endif/' \
    gtk2-interface.c > gtk2-interface.tmp
echo '#endif // (GTK_MAJOR_VERSION > 1)' >> gtk2-interface.tmp

mv gtk2-interface.tmp gtk2-interface.c 


