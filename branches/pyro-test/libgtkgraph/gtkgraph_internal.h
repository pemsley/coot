#ifndef GTK_GRAPH_INTERNAL_H
#define GTK_GRAPH_INTERNAL_H

/* This header defines function calls that are private to the internal working
of the library and should not generally be accessible to the user */

#include "gtkgraph.h"
#define TWOPI	(2.0 * M_PI)

/* Global Variables for Internal Consumption ONLY */
extern gint true_width;
extern gint true_height;
extern gint user_width;
extern gint user_height;
extern gint user_origin_x;
extern gint user_origin_y;
extern gint margin;

extern GdkPixmap *buffer;
extern GdkGC *BandWcontext;
extern GdkGC *BlueandWcontext;
extern GdkGC *Whitecontext;
extern GdkGC *Gridcontext;

/* Function prototypes in axis.c */
void gtk_graph_axis_scale_axis(GtkGraph *graph, gint user_width, gint user_height);
/* Function prototypes in trace.c */
void gtk_graph_trace_create_markers(GtkGraph *graph);
/* Function prototypes in polar.c */
void gtk_graph_polar_plot_axes(GtkGraph *graph);
void gtk_graph_polar_plot_traces(GtkGraph *graph);
/* Function prototypes in smith.c */
void gtk_graph_smith_plot_axes(GtkGraph *graph);
void gtk_graph_smith_plot_traces(GtkGraph *graph);
/* Function prototypes in polar_util.c */
gfloat wrap_angle(gfloat angle, GtkGraphPolarFormat format);
#endif
