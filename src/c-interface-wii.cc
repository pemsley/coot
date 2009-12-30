
#include <iostream>
#include <stdio.h>

#include <gtk/gtk.h>
#include "c-interface.h"

#include "graphics-info.h" // which is where the static wiimote object
			   // is stored

#ifdef WII_INTERFACE

void cwiid_message_callback(cwiid_wiimote_t *wiimote,
			    int mesg_count,
			    union cwiid_mesg mesg_array[],
			    struct timespec *timestamp);


void setup_wii() {

   if (!g_thread_supported()) {
      g_thread_init(NULL);
   }
   gdk_threads_init();
   gdk_threads_enter();
   
   bdaddr_t bdaddr;
   // simulate BDADDR_ANY?
   for (int i=0; i<6; i++)
      bdaddr.b[i] = 0;
   struct acc_cal wm_cal, nc_cal;

   
   std::cout << "Trying to connect...." << std::endl;
   graphics_info_t::wiimote = cwiid_open(&bdaddr, CWIID_FLAG_MESG_IFC);
   if (! graphics_info_t::wiimote) {
      std::cout << "wiimote failed to connect" << std::endl;
   } else {
      std::cout << "==== wiimote connected! =============== " << std::endl;
      int calibrate = cwiid_get_acc_cal(graphics_info_t::wiimote,
					CWIID_EXT_NONE, &wm_cal);
      if (! calibrate) {
	 std::cout << "Unable to retrieve accelerometer calibration"
		   << std::endl;
      } else {
	 std::cout << "coot - got wiimote accel calibration" << std::endl;
      }

      cwiid_mesg_callback_t *callback = cwiid_message_callback;
      std::cout << "setting callback" << std::endl;
      int setting_callback_status =
	 cwiid_set_mesg_callback(graphics_info_t::wiimote,
				 callback);
      if (setting_callback_status) {
	 // fail
	 std::cout << "setting message callback failed" << std::endl;
	 cwiid_close(graphics_info_t::wiimote);
      } else {
	 std::cout << "message callback set OK!" << std::endl;
	 int flags = 0x7F;
	 int enable = cwiid_enable(graphics_info_t::wiimote, flags);
      }

      
      std::cout << "getting wiimote status" << std::endl;
      cwiid_request_status(graphics_info_t::wiimote);
      std::cout << "got wiimote status" << std::endl;
   }
   std::cout << "returning from setup" << std::endl;
   gdk_threads_leave();
}


void
cwiid_message_callback(cwiid_wiimote_t *wiimote,
		       int mesg_count,
		       union cwiid_mesg mesg_array[],
		       struct timespec *timestamp) {
   
   std::cout << "------- in cwiid_message_callback! " << std::endl;
   gdk_threads_enter();
   for (int i=0; i < mesg_count; i++) {
      switch (mesg_array[i].type) {
      case CWIID_MESG_STATUS:
	 std::cout << "   a status message " << std::endl;
	 break;
      case CWIID_MESG_BTN:
	 std::cout << "   a button message " << std::endl;
	 break;
      case CWIID_MESG_ACC:
	 std::cout << "   a accel message " << std::endl;
	 break;
      case CWIID_MESG_IR:
	 std::cout << "   an IR message " << std::endl;
	 break;
      case CWIID_MESG_NUNCHUK:
	 std::cout << "   a nunchuk_mesg message " << std::endl;
	 break;
      case CWIID_MESG_CLASSIC:
	 std::cout << "   a classic message " << std::endl;
	 break;
      case CWIID_MESG_ERROR:
	 std::cout << "   an error message " << std::endl;
	 break;
      case CWIID_MESG_UNKNOWN:
	 std::cout << "  an unknown message " << std::endl;
	 break;
      }
   }
   gdk_flush();
   gdk_threads_leave();
}




#endif // WII_INTERFACE
