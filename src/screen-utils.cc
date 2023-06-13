
#include <string>
#include <iostream>
#include <glob.h>
#include <sys/stat.h>
#include <gtk/gtk.h>
#include "utils/coot-utils.hh"

void setup_application_icon(GtkWindow *window) {

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
      // 20220528-PE FIXME icons
#else

   std::string splash_screen_pixmap_dir = coot::package_data_dir();
   splash_screen_pixmap_dir += "/";
   splash_screen_pixmap_dir += "pixmaps";

   // over-ridden by user?
   char *s = getenv("COOT_PIXMAPS_DIR");
   if (s) {
      splash_screen_pixmap_dir = s;
   }

   // now add the application icon
   std::string app_icon_path = coot::util::append_dir_file(splash_screen_pixmap_dir, "coot-icon.png");

   struct stat buf;
   int status = stat(app_icon_path.c_str(), &buf);
   if (status == 0) { // icon file was found

      GdkPixbuf *icon_pixbuf =
	 gdk_pixbuf_new_from_file (app_icon_path.c_str(), NULL);
      if (icon_pixbuf) {
	 gtk_window_set_icon (GTK_WINDOW (window), icon_pixbuf);
	 g_object_unref(icon_pixbuf); // - what does this do?  I mean,
	 // why do I want to unref this pixbuf now? What does
	 // gtk_window_set_icon() do to the refcount? How do I know
	 // the refcount on a widget?
      }
   }

   // load svg/png files to antialias icons
   // maybe should go somewhere else.
   GtkIconSet* iconset;

   GtkIconFactory *iconfactory = gtk_icon_factory_new();
   GtkIconTheme *icon_theme = gtk_icon_theme_new();
   gtk_icon_theme_set_custom_theme(icon_theme, "coot");

   // GSList* stock_ids = gtk_stock_list_ids(); deprecated 
   GError *error = NULL;
   const char *stock_id;
   GdkPixbuf* pixbuf;

   glob_t myglob;
   int flags = 0;
   std::string glob_dir = splash_screen_pixmap_dir;
   std::string glob_file = glob_dir;
   glob_file += "/*.svg";
   glob(glob_file.c_str(), flags, 0, &myglob);
   glob_file = glob_dir;
   glob_file += "/*.png";
   glob(glob_file.c_str(), GLOB_APPEND, 0, &myglob);
   size_t count;
   char **p;
   for (p = myglob.gl_pathv, count = myglob.gl_pathc; count; p++, count--) {
      char *filename(*p);
      if (false) // remove the noise while I think about other things.
         std::cout << "setup_application_icon() filename " << filename << std::endl;
      pixbuf = gdk_pixbuf_new_from_file(filename, &error);
      if (error) {
	 g_print ("Error loading icon: %s\n", error->message);
	 g_error_free (error);
	 error = NULL;
      } else {
	 if (pixbuf) {
            // std::cout << "Replace gtk_icon_set_new_from_pixbuf()\n";
	    iconset = gtk_icon_set_new_from_pixbuf(pixbuf);
	    g_object_unref(pixbuf);
	    // may have to be adjusted for Windows!!
	    stock_id = coot::util::file_name_non_directory(filename).c_str();
	    std::string tmp = coot::util::file_name_non_directory(filename);
            stock_id = tmp.c_str();
            if (! tmp.empty()) {
               // std::cout << "factory adding icon " << stock_id << std::endl;
	       gtk_icon_factory_add(iconfactory, stock_id, iconset);
	       gtk_icon_factory_add_default(iconfactory);
	    }
	 }
      }
   }
   globfree(&myglob);

#endif

}
