
#include <Python.h>

#include <iostream>
#include <stdio.h>

#include <gtk/gtk.h>
#include "c-interface.h"

#include "graphics-info.h" // which is where the static wiimote object
			   // is stored

#include "globjects.h"

#include "compat/coot-sysdep.h"

#ifdef WII_INTERFACE_WIIUSE

#include "wiiuse.h"

#define MAX_WIIMOTES                            4

void handle_event(struct wiimote_t* wm) {
  
  gdk_threads_enter();

  printf("\n\n--- WII EVENT [id %i] ---\n", wm->unid);

  gint x=0, y=0;
  graphics_info_t g;
  GdkModifierType state;
  gdk_window_get_pointer(g.glarea->window, &x, &y, &state); 

  /* if a button is pressed, report it */
  if (IS_PRESSED(wm, WIIMOTE_BUTTON_A))           printf("A pressed\n");
  if (IS_PRESSED(wm, WIIMOTE_BUTTON_B))           printf("B pressed\n");
  if (IS_HELD(wm, WIIMOTE_BUTTON_LEFT)) {
    keypad_translate_xyz(1, -1);
    printf("LEFT held\n");
  }
  if (IS_HELD(wm, WIIMOTE_BUTTON_RIGHT)) {
    keypad_translate_xyz(1, 1);
    printf("RIGHT held\n");
  }
  if (IS_PRESSED(wm, WIIMOTE_BUTTON_UP)) {
    keypad_translate_xyz(2, 1);
    printf("UP pressed\n");
  }
  if (IS_PRESSED(wm, WIIMOTE_BUTTON_DOWN)) {
    keypad_translate_xyz(2, -1);
    printf("DOWN pressed\n");
  }
  if (IS_PRESSED(wm, WIIMOTE_BUTTON_LEFT)) {
    keypad_translate_xyz(1, -1);
    printf("LEFT pressed\n");
  }
  if (IS_PRESSED(wm, WIIMOTE_BUTTON_RIGHT)) {
    keypad_translate_xyz(1, 1);
    printf("RIGHT pressed\n");
  }
  if (IS_PRESSED(wm, WIIMOTE_BUTTON_MINUS))       printf("MINUS pressed\n");
  if (IS_PRESSED(wm, WIIMOTE_BUTTON_PLUS))        printf("PLUS pressed\n");
  if (IS_PRESSED(wm, WIIMOTE_BUTTON_ONE))         printf("ONE pressed\n");
  if (IS_PRESSED(wm, WIIMOTE_BUTTON_TWO))         printf("TWO pressed\n");
  if (IS_PRESSED(wm, WIIMOTE_BUTTON_HOME))        printf("HOME pressed\n");

  /*
   *      Pressing minus will tell the wiimote we are no longer interested in movement.
   *      This is useful because it saves battery power.
   */
  if (IS_JUST_PRESSED(wm, WIIMOTE_BUTTON_MINUS))
    wiiuse_motion_sensing(wm, 0);

  /*
   *      Pressing plus will tell the wiimote we are interested in movement.
   */
  if (IS_JUST_PRESSED(wm, WIIMOTE_BUTTON_PLUS))
    wiiuse_motion_sensing(wm, 1);

  /* if the accelerometer is turned on then print angles */
  if (WIIUSE_USING_ACC(wm)) {
    float roll_smooth = wm->orient.a_roll;
    float pitch_smooth = wm->orient.a_pitch;
    printf("wiimote roll  = %f [%f]\n", wm->orient.roll, wm->orient.a_roll);
    if (fabs(roll_smooth) > 7.) {
      rotate_y_scene(abs(int(roll_smooth)), roll_smooth/10.);
    }
    if (fabs(pitch_smooth) > 7.) {
      rotate_x_scene(abs(int(pitch_smooth)), pitch_smooth/10.);
      printf("BL DEBUG:: abs pitch %i\n", abs(int(pitch_smooth)));
    }
    printf("wiimote pitch = %f [%f]\n", wm->orient.pitch, wm->orient.a_pitch);
    printf("wiimote yaw   = %f\n", wm->orient.yaw);
  }



  /*
   *      Pressing B will toggle the rumble
   *
   *      if B is pressed but is not held, toggle the rumble
   */
  if (IS_JUST_PRESSED(wm, WIIMOTE_BUTTON_B))
    wiiuse_toggle_rumble(wm);

  if (IS_JUST_PRESSED(wm, WIIMOTE_BUTTON_UP))
    wiiuse_set_ir(wm, 1);
  if (IS_JUST_PRESSED(wm, WIIMOTE_BUTTON_DOWN))
    wiiuse_set_ir(wm, 0);

  gdk_flush(); // needed?
}

void handle_ctrl_status(struct wiimote_t* wm) {
	printf("\n\n--- CONTROLLER STATUS [wiimote id %i] ---\n", wm->unid);

	printf("attachment:      %i\n", wm->exp.type);
	printf("speaker:         %i\n", WIIUSE_USING_SPEAKER(wm));
	printf("ir:              %i\n", WIIUSE_USING_IR(wm));
	printf("leds:            %i %i %i %i\n", WIIUSE_IS_LED_SET(wm, 1), WIIUSE_IS_LED_SET(wm, 2), WIIUSE_IS_LED_SET(wm, 3), WIIUSE_IS_LED_SET(wm, 4));
	printf("battery:         %f %%\n", wm->battery_level);
}


int setup_wii() {
  
  int found, connected;
  graphics_info_t g;

  if (!g_thread_supported()) {
    g_thread_init(NULL);
  }
  gdk_threads_init();
  gdk_threads_enter();
   
  /*
   *      Initialize an array of wiimote objects.
   *
   *      The parameter is the number of wiimotes I want to create.
   */
  g.wiimotes =  wiiuse_init(MAX_WIIMOTES);
  found = wiiuse_find(g.wiimotes, MAX_WIIMOTES, 5);
  if (!found) {
    printf ("No wiimotes found.\n");
    return -1;
  }

  connected = wiiuse_connect(g.wiimotes, MAX_WIIMOTES);
  if (connected)
    printf("Connected to %i wiimotes (of %i found).\n", connected, found);
  else {
    printf("Failed to connect to any wiimote.\n");
    return -1;
  }
  /*
   *      Now set the LEDs and rumble for a second so it's easy
   *      to tell which wiimotes are connected (just like the wii does).
   */
  wiiuse_set_leds(g.wiimotes[0], WIIMOTE_LED_1);
  wiiuse_set_leds(g.wiimotes[1], WIIMOTE_LED_2);
  wiiuse_set_leds(g.wiimotes[2], WIIMOTE_LED_3);
  wiiuse_set_leds(g.wiimotes[3], WIIMOTE_LED_4);

  wiiuse_rumble(g.wiimotes[0], 1);
  wiiuse_rumble(g.wiimotes[1], 1);

  coot::usleep(200000);

  wiiuse_rumble(g.wiimotes[0], 0);
  wiiuse_rumble(g.wiimotes[1], 0);

  // battery warning, maybe should be a timeout
  for (int i = 0; i < found; ++i) {
    std::string txt = "";
    if (g.wiimotes[i]->battery_level < 0.1) {
      txt = "WARNING:: Batteries of Wiimote (";
      txt += coot::util::int_to_string(i + 1);
      txt += ") is low! (";
      txt += coot::util::float_to_string(g.wiimotes[i]->battery_level * 100.);
      txt += "%)";
    }
    if (txt != "")
      info_dialog(txt.c_str());
  }


  //while (1) {
  while (g.wiimotes) {
    if (wiiuse_poll(g.wiimotes, MAX_WIIMOTES)) {
      /*
       *      This happens if something happened on any wiimote.
       *      So go through each one and check if anything happened.
       */
      for (int i = 0; i < MAX_WIIMOTES; ++i) {
        switch (g.wiimotes[i]->event) {
        case WIIUSE_EVENT:
          /* a generic event occured */
          handle_event(g.wiimotes[i]);
          break;

        case WIIUSE_STATUS:
          /* a status event occured */
          handle_ctrl_status(g.wiimotes[i]);
          break;        

        case WIIUSE_DISCONNECT:
          printf("Disconnected wiimotes %i.\n", i);
          return 0;
        }
      }
    }
    while (gtk_events_pending ())
	  gtk_main_iteration ();
  }
  gdk_threads_leave();
  wiiuse_cleanup(g.wiimotes, MAX_WIIMOTES);
  return 0;

}

void stop_wii() {
  // stop all wii 
  graphics_info_t g;

  wiiuse_cleanup(g.wiimotes, MAX_WIIMOTES);
  g.wiimotes = NULL;
}

void wii_status() {
  graphics_info_t g;
  for (int i = 0; i < MAX_WIIMOTES; ++i) {
    wiiuse_status(g.wiimotes[i]);
  }
}

#endif // WII_INTERFACE_WIIUSE

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
