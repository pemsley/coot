/*
 * File: gtkgraph.h
 * Auth: Andrew Hurrell
 */
#ifndef __GTK_GRAPH_H__
#define __GTK_GRAPH_H__

#include <gdk/gdk.h>
#include <gtk/gtk.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* Useful Global Variables */
extern GdkColor BLACK;
extern GdkColor WHITE;
extern GdkColor RED;
extern GdkColor GREEN;
extern GdkColor BLUE;
extern GdkColor LIGHT_GREY;
extern GdkColor MID_GREY;
extern GdkColor DARK_GREY;
extern GdkColor CYAN;
extern GdkColor MAGENTA;
extern GdkColor YELLOW;
extern GdkColor ORANGE;
extern GdkColor NAVY_BLUE;
extern GdkColor OLIVE_GREEN;
extern GdkColor PURPLE;
extern GdkColor BROWN;


/* 
 * --- Macros for conversion and type checking 
 */
#define GTK_GRAPH(obj)          	GTK_CHECK_CAST (obj, gtk_graph_get_type (), GtkGraph)
#define GTK_GRAPH_CLASS(klass)  GTK_CHECK_CLASS_CAST (klass, gtk_graph_get_type, GtkGraphClass)
#define GTK_IS_GRAPH(obj)       	GTK_CHECK_TYPE (obj, gtk_graph_get_type ())

/* --- enumerated types --- */

enum _GtkGraphAxisType
{
GTK_GRAPH_AXIS_INDEPENDANT	= 1 << 0,
GTK_GRAPH_AXIS_DEPENDANT	= 1 << 1
};
	

enum _GtkGraphType
{
XY    	= 1 << 0,
POLAR 	= 1 << 1,
SMITH 	= 1 << 2
};

enum _GtkGraphPolarType
{
DEGREES_SYMMETRIC     = 1 << 0,
DEGREES_ANTISYMMETRIC = 1 << 1,
RADIANS_SYMMETRIC     = 1 << 2,
RADIANS_ANTISYMMETRIC = 1 << 3
};

enum _GtkGraphNumberFormat
{
FLOATING_POINT	= 1 << 0,
ENGINEERING		= 1 << 1,
SCIENTIFIC		= 1 << 2
};

enum _GtkGraphCrossingType
{
GTK_GRAPH_AXISMAX   = 1 << 0,
GTK_GRAPH_AXISMIN   = 1 << 1,
GTK_GRAPH_USERVALUE = 1 << 2
};

enum _GtkGraphMarkerType
{
GTK_GRAPH_MARKER_NONE			= 1 << 0,
GTK_GRAPH_MARKER_SQUARE        	= 1 << 1,
GTK_GRAPH_MARKER_CIRCLE        	= 1 << 2,
GTK_GRAPH_MARKER_DIAMOND      	= 1 << 3,
GTK_GRAPH_MARKER_TRIANGLE      	= 1 << 4,
GTK_GRAPH_MARKER_PLUS           = 1 << 5,
GTK_GRAPH_MARKER_CROSS			= 1 << 6,
GTK_GRAPH_MARKER_STAR			= 1 << 7
};

enum _GtkGraphLineType
{
NO_LINE		= 1 << 0,
SOLID		= 1 << 1,
DASHED		= 1 << 2,
DASH_DOT	= 1 << 3,
DOTTED		= 1 << 4,
LONG_DASH	= 1 << 5
};
enum _GtkGraphAnnotationType
{
VERTICAL 	= 1 << 0,
HORIZONTAL 	= 1 << 1,
RADIAL 		= 1 << 2,
AZIMUTHAL	= 1 << 3,
VSWR		= 1 << 4,
Q		= 1 << 5
};

enum _GtkGraphPosition
{
GTK_GRAPH_NORTH         =   1 << 0,
GTK_GRAPH_NORTH_EAST    =   1 << 1,
GTK_GRAPH_EAST          =   1 << 2,
GTK_GRAPH_SOUTH_EAST    =   1 << 3,
GTK_GRAPH_SOUTH         =   1 << 4,
GTK_GRAPH_SOUTH_WEST    =   1 << 5,
GTK_GRAPH_WEST          =   1 << 6,
GTK_GRAPH_NORTH_WEST    =   1 << 7
};

/* --- Defining data structures. ---  */

typedef enum _GtkGraphAnnotationType GtkGraphAnnotationType;
typedef enum _GtkGraphType GtkGraphType;
typedef enum _GtkGraphPolarType GtkGraphPolarType;
typedef enum _GtkGraphNumberFormat GtkGraphNumberFormat;
typedef enum _GtkGraphCrossingType GtkGraphCrossingType;
typedef enum _GtkGraphMarkerType GtkGraphMarkerType;
typedef enum _GtkGraphLineType GtkGraphLineType;
typedef enum _GtkGraphPosition GtkGraphPosition;
typedef enum _GtkGraphAxisType GtkGraphAxisType;

struct _GtkGraphAxis
{
    gfloat maj_tick;
	gint n_maj_tick;
    gfloat min_tick;
	gint n_min_tick;
    gfloat axis_min;
    gfloat axis_max;
    GtkGraphCrossingType crossing_type; /* Should the other axis cross at max, min or user  defined value*/

    gfloat crossing_value;		/* The value at which the other axis cross*/
    GtkGraphNumberFormat format;
	gint precision;
    gchar *title;
	gboolean autoscale_limits;
	gboolean autoscale_tick;
 	gboolean grid_visible;
	
	gint pxls_per_maj_tick;
	gfloat scale_factor;
};

typedef struct _GtkGraphAxis GtkGraphAxis; 

struct _GtkGraphTraceFormat /* The formatting of individual traces*/
{
    GtkGraphMarkerType marker_type;
    gboolean marker_filled;
    GdkPixmap *marker;
    GdkGC *marker_gc;
	gint marker_size;
	
	GdkBitmap *mask;
    GdkGC *mask_gc;

    GtkGraphLineType line_type;
	GdkGC *line_gc;
    gboolean line_visible;
	
    gchar *legend_text;
};
typedef struct _GtkGraphTraceFormat GtkGraphTraceFormat;


struct _GtkGraphPolarFormat
{
	GtkGraphPolarType type;
	gfloat polar_start;
};
typedef struct _GtkGraphPolarFormat GtkGraphPolarFormat;

typedef struct _GtkGraphAnnotation GtkGraphAnnotation;
struct _GtkGraphAnnotation
{
GtkGraphAnnotationType type;
gfloat value;
gchar *text;
GtkGraphAnnotation *next;
}; 


typedef struct _GtkGraphTrace GtkGraphTrace;
struct _GtkGraphTrace
{
    gfloat *Xdata;
    gfloat *Ydata;
    gint   num_points;
    gfloat xmax;
    gfloat xmin;
    gfloat ymax;
    gfloat ymin;
    GtkGraphTraceFormat *format;
    GtkGraphTrace *next;
};


/*
 * GtkGraph:
 *
 * Contains all the structure of a GtkGraph
 */

typedef struct _GtkGraph GtkGraph;
struct _GtkGraph
{
GtkWidget drawing_area;		/* We need a windowed widget to act as placeholder   */
GtkGraphTrace *traces;
gint num_traces;
GtkGraphAnnotation *annotations;
gint num_annotations;	
GtkGraphAxis *dependant;
GtkGraphAxis *independant;
GtkGraphType graph_type;
GtkGraphPolarFormat polar_format;
gchar *title;
gchar *subtitle;
gint legend_visible;
GtkGraphPosition legend_position;
gfloat smith_Z0;
GtkGraph *next;
};



struct _GtkGraphClass
{
  GtkWidgetClass parent_class;
};
typedef struct _GtkGraphClass GtkGraphClass;



/* Function prototypes in gtkgraph.c */


guint gtk_graph_get_type (void);
GtkWidget *gtk_graph_new (GtkGraphType type);
void gtk_graph_set_title (GtkGraph *graph, gchar *title, gchar *subtitle);
void gtk_graph_legend_format(GtkGraph *graph, gboolean is_visible, GtkGraphPosition position);

/* Function prototypes in trace.c */

gint gtk_graph_trace_new(GtkGraph *graph);
void gtk_graph_trace_set_data(GtkGraph *graph, gint trace_id, gfloat *xd, gfloat *yd, gint n);
void gtk_graph_trace_format_line(GtkGraph *graph, gint trace_id, GtkGraphLineType type, gint width, GdkColor *line_color, gboolean visible);
void gtk_graph_trace_format_marker(GtkGraph *graph, gint trace_id, GtkGraphMarkerType type, gint marker_size, GdkColor *fg, GdkColor *bg, gboolean is_filled);
void gtk_graph_trace_format_title(GtkGraph *graph, gint trace_id, gchar *legend_text);

/* Function prototypes in axis.c */
void gtk_graph_axis_set_limits(GtkGraph *graph, GtkGraphAxisType axis, gfloat max, gfloat min);
void gtk_graph_axis_set_tick(GtkGraph *graph, GtkGraphAxisType axis, gfloat majtick, gfloat mintick);
void gtk_graph_axis_format(GtkGraph *graph, GtkGraphAxisType axis, GtkGraphNumberFormat number_format, gint precision, gchar *title);
void gtk_graph_axis_set_crossing(GtkGraph *graph, GtkGraphAxisType axis, GtkGraphCrossingType type, gfloat value);
void gtk_graph_axis_format_grid(GtkGraph *graph, GtkGraphAxisType axis, gboolean visible);

/* Function prototypes in annotations.c */
gint gtk_graph_annotation_new(GtkGraph *graph);
void gtk_graph_annotation_set_data(GtkGraph *graph, gint annotation_id, GtkGraphAnnotationType type, gfloat value, gchar *text);
void gtk_graph_plot_annotations(GtkGraph *graph);

/* Function prototypes in polar.c */
void gtk_graph_set_polar_format(GtkGraph *graph, GtkGraphPolarType type, gfloat polar_start);

/* Function prototypes in smith.c */
void gtk_graph_smith_set_Z0(GtkGraph *graph, gfloat Z0);
#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* __GTK_GRAPH_H__ */
