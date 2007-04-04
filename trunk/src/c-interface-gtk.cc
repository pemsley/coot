

#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)

# else 

#include <iostream>
#include "c-interface.h"

void fileselection_sort_button_clicked_gtk1( GtkWidget *sort_button,
					     GtkWidget  *file_list) {

   std::cout << "---- sort button pressed!" << std::endl;
}

#endif // GTK2
