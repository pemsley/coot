
/  g_signal_connect ..gpointer. model_toolbar_rot_trans_toolbutton, ".*",/ { 

    print "#ifdef GTK_TYPE_MENU_TOOL_BUTTON";
    print $0;
    endif_next_null = 1;
    next;
} 

$1 == "NULL);" { 
    print $0;
    if (endif_next_null) { 
       print "#endif";
       endif_next_null = 0;
    }
    next; 
}

{ print $0; } 



