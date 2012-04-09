/*
 * File: gtkgraph.h
 * Auth: Andrew Hurrell
 */

/*! \mainpage GtkGraph: The Scientific Graphing Widget for GTK+
 *
 * \author Andrew Hurrell
 *
 * \section intro_sec Introduction
 *
 * GtkGraph is a simple to use widget for GTK+ V2.0 designed to allow easy presentation of
 * scientific data.  It has been designed with ease of use in mind and simple graphs can be
 * produced with only a few function calls.  Additional functions are provided to allow complete
 * control of the format in which data traces are displayed, and scaling of the graph axes.  A
 * selection of graph annotations are also available to enhance display or draw attention to
 * particular regions of the graph.
 *
 *Currently three types of graph are supported by GtkGraph:
 *
 *XY plots - where data pairs in the form (x co-ordinate, y co-ordinate) are plotted on a
 *conventional Cartesian XY Graph
 *
 *Polar plots - where data pairs in the form (radial distance, angle) are plotted on a circular
 *Polar Graph
 *
 *Smith charts - where complex electrical impedance or reflection coefficient data pairs are
 *plotted on a Smith Chart
 *
 * \section install_sec Installation
 *
 *GtkGraph uses the standard Gnu Build Tools therefore it should be sufficient to type
 *open a terminal window and type ./configure.  Once configure has done it's stuff then type
 *make, and you should then have the static library built
 *  
 */
 
 /*! \file gtkgraph.h
    \brief The main header file for GtkGraph
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
 
/*! Cast the object to a #GtkGraph */
#define GTK_GRAPH(obj) GTK_CHECK_CAST (obj, gtk_graph_get_type (), GtkGraph)
#define GTK_GRAPH_CLASS(klass) GTK_CHECK_CLASS_CAST (klass, gtk_graph_get_type, GtkGraphClass)
/*! Query whether the object is a #GtkGraph */
#define GTK_IS_GRAPH(obj) GTK_CHECK_TYPE (obj, gtk_graph_get_type ())

/* --- enumerated types --- */

/*!  Enumerated type to specify which type of graph is required */
typedef enum 
{
XY    	= 1 << 0,
POLAR 	= 1 << 1,
SMITH 	= 1 << 2
} GtkGraphType;

/*! \var GtkGraphType XY
 * An XY (Scatter) graph */
/*! \var GtkGraphType POLAR
 * A Polar plot */
 /*! \var GtkGraphType SMITH
 * A Smith's Chart */


/*! Enumerated type to specify which axis changes apply to */

typedef enum 
{
GTK_GRAPH_AXIS_INDEPENDANT = 1 << 0,
GTK_GRAPH_AXIS_DEPENDANT = 1 << 1
} GtkGraphAxisType;

/*! \var GtkGraphAxisType GTK_GRAPH_AXIS_INDEPENDANT
 * The independant axis of the graph; for XY graphs this is the X-axis,
 * for POLAR graphs this is the angle */
/*! \var GtkGraphAxisType GTK_GRAPH_AXIS_DEPENDANT
 * The dependant axis of the graph; for XY graphs this is the Y-axis,
 * for POLAR graphs this is the radial distance */


/*! Enumerated type to specify radial range of POLAR Graphs */
typedef enum 
{
DEGREES_SYMMETRIC     = 1 << 0,
DEGREES_ANTISYMMETRIC = 1 << 1,
RADIANS_SYMMETRIC     = 1 << 2,
RADIANS_ANTISYMMETRIC = 1 << 3
}GtkGraphPolarType;

/*! \var GtkGraphPolarType DEGREES_SYMMETRIC
 * Angle measured in degrees from -180 to 180 */
/*! \var GtkGraphPolarType DEGREES_ANTISYMMETRIC
 * Angle measured in degrees from 0 to 360 */
 /*! \var GtkGraphPolarType RADIANS_SYMMETRIC
 * Angle measured in radians from -PI to PI */
/*! \var GtkGraphPolarType RADIANS_ANTISYMMETRIC
 * Angle measured in radians from 0 to 2 PI */

/*! Enumerated type to specify how numbers should be displayed */
typedef enum 
{
FLOATING_POINT = 1 << 0,
ENGINEERING = 1 << 1,
SCIENTIFIC = 1 << 2
} GtkGraphNumberFormat;

/*! \var GtkGraphNumberFormat FLOATING_POINT
 * Standard Floating point */
/*! \var GtkGraphNumberFormat ENGINEERING
 * Engineering notation (i.e. 1.2E+3) but with exponent always as powers of 3 */
 /*! \var GtkGraphNumberFormat SCIENTIFIC
 * Scientific notation (i.e. 2.7E+4) but with exponent taking any integer value */


/*!  Enumerated type to specify where one axis should cross the other on an XY graph*/
typedef enum 
{
GTK_GRAPH_AXISMAX = 1 << 0,
GTK_GRAPH_AXISMIN = 1 << 1,
GTK_GRAPH_USERVALUE = 1 << 2
} GtkGraphCrossingType;
/*! \var GtkGraphCrossingType GTK_GRAPH_AXISMAX
 * The other axis should cross at the minumum value of this axis */
/*! \var GtkGraphCrossingType GTK_GRAPH_AXISMIN
 * The other axis should cross at the maximum value of this axis */
 /*! \var GtkGraphCrossingType GTK_GRAPH_USERVALUE
 * The other axis should cross at a user specified value on this axis */


/*! Enumerated type to specify which markers should be used to identify data points
 * on a trace*/
typedef enum 
{
GTK_GRAPH_MARKER_NONE = 1 << 0,
GTK_GRAPH_MARKER_SQUARE = 1 << 1,
GTK_GRAPH_MARKER_CIRCLE	= 1 << 2,
GTK_GRAPH_MARKER_DIAMOND = 1 << 3,
GTK_GRAPH_MARKER_TRIANGLE = 1 << 4,
GTK_GRAPH_MARKER_PLUS = 1 << 5,
GTK_GRAPH_MARKER_CROSS = 1 << 6,
GTK_GRAPH_MARKER_STAR = 1 << 7
} GtkGraphMarkerType;
/*! \var GtkGraphMarkerType GTK_GRAPH_MARKER_NONE
 * No marker */
/*! \var GtkGraphMarkerType GTK_GRAPH_MARKER_SQUARE
 * A Square */
 /*! \var GtkGraphMarkerType GTK_GRAPH_MARKER_CIRCLE
 * A Circle */
/*! \var GtkGraphMarkerType GTK_GRAPH_MARKER_DIAMOND
 * A Diamond */
/*! \var GtkGraphMarkerType GTK_GRAPH_MARKER_TRIANGLE
 * A Triangle */
 /*! \var GtkGraphMarkerType GTK_GRAPH_MARKER_PLUS
 * A Plus (+) */
/*! \var GtkGraphMarkerType GTK_GRAPH_MARKER_CROSS
 * A Cross (X) */
/*! \var GtkGraphMarkerType GTK_GRAPH_MARKER_STAR
 * A Star (*) */


/*!  Enumerated type to specify which line should be used to draw a trace*/
typedef enum 
{
NO_LINE = 1 << 0,
SOLID = 1 << 1,
DASHED = 1 << 2,
DASH_DOT = 1 << 3,
DOTTED = 1 << 4,
LONG_DASH = 1 << 5
} GtkGraphLineType;
/*! \var GtkGraphLineType NO_LINE
 * No line */
 /*! \var GtkGraphLineType SOLID
 * Solid line */
 /*! \var GtkGraphLineType DASHED
 * Dashed line */
 /*! \var GtkGraphLineType DASH_DOT
 * Alternate dash-dot line */
 /*! \var GtkGraphLineType DOTTED
 * Dotted line */
 /*! \var GtkGraphLineType LONG_DASH
 * Long Dashed line */


/*!  Enumerated type to specify the nature of any annotations drawn on a graph*/
typedef enum
{
VERTICAL = 1 << 0,
HORIZONTAL = 1 << 1,
RADIAL = 1 << 2,
AZIMUTHAL = 1 << 3,
VSWR = 1 << 4,
Q = 1 << 5
} GtkGraphAnnotationType;
/*! \var GtkGraphAnnotationType VERTICAL
 * Applies to XY Graphs: a Vertical Annotation line */
 /*! \var GtkGraphAnnotationType HORIZONTAL
 * Applies to XY Graphs: a Horizontal Annotation line */
  /*! \var GtkGraphAnnotationType RADIAL
 * Applies to Polar Plots: an Annotation line from the centre to the edge*/
  /*! \var GtkGraphAnnotationType AZIMUTHAL
 * Applies to Polar Plots: a Circular annotation line */
  /*! \var GtkGraphAnnotationType VSWR
 * Applies to Smith Charts: A Line of constant VSWR */
  /*! \var GtkGraphAnnotationType Q
 * Applies to Smith Charts: A Line of constant Q */
 
/*! Enumerated type to specify legend positions on a graph */ 
typedef enum
{
GTK_GRAPH_NORTH         =   1 << 0,/*! Top Edge Centred Horizontally */
GTK_GRAPH_NORTH_EAST    =   1 << 1,/*! Top Right Corner */
GTK_GRAPH_EAST          =   1 << 2,/*! Left Edge Centred Vertically*/
GTK_GRAPH_SOUTH_EAST    =   1 << 3,/*! Bottom Right Corner */
GTK_GRAPH_SOUTH         =   1 << 4,/*! Bottom Edge Centred Horizontally */
GTK_GRAPH_SOUTH_WEST    =   1 << 5,/*! Bottom Left Corner */
GTK_GRAPH_WEST          =   1 << 6,/*! Right Edge Centred Vertically*/
GTK_GRAPH_NORTH_WEST    =   1 << 7 /*! Top Left Corner */
} GtkGraphPosition;

/* --- Defining data structures. ---  */

/*!GtkGraphAxis describes the content and formatting of each axis.
 *  This is an opaque structure and  * its contents should not be written
 * to directly.  Instead please use the API functions gtk_graph_axis_set_limits()
 * gtk_graph_axis_format(), gtk_graph_axis_set_tick(), gtk_graph_axis_set_crossing()
 * or gtk_graph_axis_format_grid() to modify the contents of this structure
 */

typedef struct 
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
} GtkGraphAxis;


/*!GtkGraphTraceFormat structure contains all of the formatting information
 * for the formatting of individual traces.  This is an opaque structure and
 * its contents should not be written to directly.  Instead please use the
 * API functions such as gtk_graph_trace_format_line(), 
 * gtk_graph_trace_format_marker() or gtk_graph_trace_format_title() to modify
 * the contents of this structure
 */

typedef struct 
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
} GtkGraphTraceFormat;

/*!GtkGraphPolarFormat structure contains additional formatting information
 * for polar graphs.  This is an opaque structure and its contents should not
 * be written to directly.  Instead please use the API functions 
 * gtk_graph_set_polar_format() to modify the contents of this structure
 */

typedef struct 
{
GtkGraphPolarType type;
gfloat polar_start;
} GtkGraphPolarFormat;

/*! \struct GtkGraphAnnotation structure describes the formatting of each annotation 
 * added to a graph.  This is an opaque structure and its contents should not
 * be written to directly.  Instead please use the API function 
 * gtk_graph_annotation_set_data() to modify the contents of this structure
 */

typedef struct _GtkGraphAnnotation GtkGraphAnnotation;
struct _GtkGraphAnnotation
{
GtkGraphAnnotationType type;
gfloat value;
gchar *text;
GtkGraphAnnotation *next;
}; 

/*! \struct GtkGraphTrace structure contains the data describing each individual trace
 * of a graph.  This is an opaque structure and its contents should not
 * be written to directly.  Instead please use the API function 
 * gtk_graph_trace_set_data(), gtk_graph_trace_format_line(),
 * gtk_graph_trace_format_marker() and gtk_graph_trace_format_title()
 * to modify the contents of this structure*/

typedef struct _GtkGraphTrace GtkGraphTrace;
struct _GtkGraphTrace
{
gfloat *Xdata;
gfloat *Ydata;
gint num_points;
gfloat xmax;
gfloat xmin;
gfloat ymax;
gfloat ymin;
GtkGraphTraceFormat *format;
GtkGraphTrace *next;
} ;


/*! GtkGraph structure is the parent structure of a GtkGraph This is an opaque
 * structure and its contents should not be written to directly. A GtkGraph is 
 * created with gtk_graph_new() and its components are modified using the
 * many available API functions
*/

typedef struct 
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
} GtkGraph;



typedef struct 
{
GtkWidgetClass parent_class;
} GtkGraphClass;

/* Function prototypes in gtkgraph.c */

guint gtk_graph_get_type (void);

/*! \brief Create a new GtkGraph
* \param type the type of graph to be created
* \return A GtkWidget pointer to the graph just created
*/GtkWidget *gtk_graph_new (GtkGraphType type);

/*! \brief Set the title and/or subtitle of a GtkGraph
* \param graph the GtkGraph that is having its title/subtitle modified
* \param title a string containing the title text
* \param subtitle a string containing the subtitle text */
void gtk_graph_set_title (GtkGraph *graph, const gchar *title, const gchar *subtitle);

/*! \brief Format the legend of a GtkGraph
* \param graph the GtkGraph that is having its legend modified
* \param is_visible whether the legend is visible or not
* \param position the position of the legend on the GtkGraph */
void gtk_graph_legend_format(GtkGraph *graph, gboolean is_visible, GtkGraphPosition position);

/* Function prototypes in trace.c */

/*! \brief Create a new trace with a GtkGraph.
 *  \param graph the GtkGraph that will contain the newly created trace.
 *  \return an unique integer identifier of the newly created trace.
 *  \note The return value must be stored for by the user for later use (e.g. by gtk_graph_trace_set_data())*/
gint gtk_graph_trace_new(GtkGraph *graph);

/*! \brief Used to allocate a data set to a trace
 *  \param graph the GtkGraph that contains the trace being modified
 *  \param trace_id the unique integer identifier of the trace as returned by
 *  gtk_graph_trace_new()
 *  \param xd an array containing the \a x co-ordinates (must be n points long)
 *  \param yd an array containing the \a y co-ordinates (must be n points long)
 *  \param n the number of points in the x and y arrays */
void gtk_graph_trace_set_data(GtkGraph *graph, gint trace_id, gfloat *xd, gfloat *yd, gint n);

/*! \brief Change the formatting of the lines used to draw a particular trace
 *  \param graph the GtkGraph that contains the trace being modified
 *  \param trace_id the unique integer identifier of the trace as returned by
 *  gtk_graph_trace_new()
 *  \param type what type of line (see #GtkGraphLineType for details)
 *  \param width how wide is the line drawn
 *  \param line_color what color to draw the line
 *  \param visible whether the trace be drawn or not */
void gtk_graph_trace_format_line(GtkGraph *graph, gint trace_id, GtkGraphLineType type, gint width, GdkColor *line_color, gboolean visible);

/*! \brief Change the formatting of the marker used when drawing a particular trace
 *  \param graph the GtkGraph that contains the trace being modified
 *  \param trace_id the unique integer identifier of the trace as returned by
 *  gtk_graph_trace_new()
 *  \param type what type of marker (see #GtkGraphMarkerType for details) 
 *  \param marker_size size of marker (from 1-6)
 *  \param fg the foreground colour of the marker
 *  \param bg the background colour of the marker
 *  \param is_filled whether the marker should be filled or not  */ 
void gtk_graph_trace_format_marker(GtkGraph *graph, gint trace_id, GtkGraphMarkerType type, gint marker_size, GdkColor *fg, GdkColor *bg, gboolean is_filled);

/*! \brief Change the legend text for a particular trace
 *  \param graph the GtkGraph that contains the trace being modified
 *  \param trace_id the unique integer identifier of the trace as returned by
 *  gtk_graph_trace_new()
 *  \param legend_text a string containing the trace title as it will appear in the legend*/
void gtk_graph_trace_format_title(GtkGraph *graph, gint trace_id, const gchar *legend_text);

/* Function prototypes in axis.c */

/*! \brief set the maximum extents (limiting values) for one axis of a graph
 *  \param graph the GtkGraph that contains the axis being modified
 *  \param axis the axis being modified
 *  \param max the maximum value to display
 *  \param min the minimum value to display*/
void gtk_graph_axis_set_limits(GtkGraph *graph, GtkGraphAxisType axis, gfloat max, gfloat min);

/*! \brief set the tick increment for one axis of a graph
 *  \param graph the GtkGraph that contains the axis being modified
 *  \param axis the axis being modified
 *  \param majtick the spacing between major ticks
 *  \param mintick the spacing between minor ticks*/
void gtk_graph_axis_set_tick(GtkGraph *graph, GtkGraphAxisType axis, gfloat majtick, gfloat mintick);

/*! \brief set the format of numerical and axis labels for one axis of a graph
 *  \param graph the GtkGraph that contains the axis being modified
 *  \param axis the axis being modified
 *  \param number_format the format of numerical labels
 *  \param precision the number of decimal points to which numerical values should be displayed
 *  \param title a string containg the title to be displayed against that axis*/
void gtk_graph_axis_format(GtkGraph *graph, GtkGraphAxisType axis, GtkGraphNumberFormat number_format, gint precision, const gchar *title);

/*! \brief specifies how one axis should cross another on an XY graph
 *  \param graph the GtkGraph that contains the axis being modified
 *  \param axis the axis being modified
 *  \param type the #GtkGraphCrossingType specifying where the other axis should cross the axis currently being modified
 *  \param value if type was set to #GTK_GRAPH_USERVALUE this holds the value that the axis should cross at*/
void gtk_graph_axis_set_crossing(GtkGraph *graph, GtkGraphAxisType axis, GtkGraphCrossingType type, gfloat value);

/*! \brief specifies whether grid lines should extend outwards from the axis
 *  \param graph the GtkGraph that contains the axis being modified
 *  \param axis the axis being modified
 *  \param visible whether grid lines are visible for this axis or not*/
void gtk_graph_axis_format_grid(GtkGraph *graph, GtkGraphAxisType axis, gboolean visible);

/* Function prototypes in annotations.c */

/*! \brief Adds an annotation to a pre-existing GtkGraph graph.
 * \param graph The GtkGraph that you wish to add the annotation to
 * \return an unique integer identifying the new annotation.  
 * \note The identifier returned by this function must be
 * stored by the user if you wish to alter the properties of the annotation */
gint gtk_graph_annotation_new(GtkGraph *graph);

/*! \brief Determines where and what type of annotation should be displayed
 * \param graph  the GtkGraph that contains the annotation you wish to modify
 * \param annotation_id  the unique integer ID returned by gtk_graph_annotation_new()
 * \param type the type of annotation to display
 * \param value value at which to display the annotation
 * \param text any text that you want displayed next to the annotation */
void gtk_graph_annotation_set_data(GtkGraph *graph, gint annotation_id, GtkGraphAnnotationType type, gfloat value, gchar *text);

 
void gtk_graph_plot_annotations(GtkGraph *graph);

/* Function prototypes in polar.c */

/*! \brief Determines the formatting of Polar graphs
 * \param graph the #GtkGraph that contains the annotation you wish to modify
 * \param type the type of polar graph see #GtkGraphPolarType for more details
 * \param polar_start */
 void gtk_graph_set_polar_format(GtkGraph *graph, GtkGraphPolarType type, gfloat polar_start);

/* Function prototypes in smith.c */

/*! \brief All Smith's Charts are normalised to a particular impendance (Z0).  This function
 * allows the user to normalise the Smith's Chart to a value to something other than the
 * default value of 50 Ohms
 * \param graph the #GtkGraph containing the Smiths Chart that you wish to re-normalise
 * \param Z0 the new normalisation value*/
void gtk_graph_smith_set_Z0(GtkGraph *graph, gfloat Z0);
#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* __GTK_GRAPH_H__ */
