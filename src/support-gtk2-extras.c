
#include <stdio.h>		/* for stderr on FC 3 */
#include "support.h"

#if (GTK_MAJOR_VERSION > 1) 
GList *pixmaps_directories;
gchar* find_pixmap_file(const gchar     *filename);
#endif 

#if (GTK_MAJOR_VERSION > 1) 
/* This is from GTK2's support.c - needed because create_pixbuf is
   used in the new aboutdialog */
/* This is an internally used function to create pixmaps. */
GdkPixbuf*
create_pixbuf                          (const gchar     *filename)
{
  gchar *pathname = NULL;
  GdkPixbuf *pixbuf;
  GError *error = NULL;

  if (!filename || !filename[0])
      return NULL;

  pathname = find_pixmap_file (filename);

  if (!pathname)
    {
      g_warning (_("Couldn't find pixmap file: %s"), filename);
      return NULL;
    }

  pixbuf = gdk_pixbuf_new_from_file (pathname, &error);
  if (!pixbuf)
    {
      fprintf (stderr, "Failed to load pixbuf file: %s: %s\n",
               pathname, error->message);
      g_error_free (error);
    }
  g_free (pathname);
  return pixbuf;
}
#endif /*  (GTK_MAJOR_VERSION > 1)  */

#if (GTK_MAJOR_VERSION > 1) 
/* This is an internally used function to find pixmap files. */
/* static */
gchar*
find_pixmap_file                       (const gchar     *filename)
{
  GList *elem;

  /* We step through each of the pixmaps directory to find it. */
  elem = pixmaps_directories;
  while (elem)
    {
      gchar *pathname = g_strdup_printf ("%s%s%s", (gchar*)elem->data,
                                         G_DIR_SEPARATOR_S, filename);
      if (g_file_test (pathname, G_FILE_TEST_EXISTS))
        return pathname;
      g_free (pathname);
      elem = elem->next;
    }
  return NULL;
}
#endif /*  (GTK_MAJOR_VERSION > 1)  */
