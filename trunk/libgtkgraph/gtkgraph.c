/*
 * File: GtkGraph.c
 * Auth: Andrew Hurrell  *
 * Simple gtk graphing widget
 */

#include <math.h>		/* for rint */
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <gtk/gtk.h>

#include "gtkgraph.h"
#include "gtkgraph_internal.h"

GdkColor BLACK 		= {0, 0x0000, 0x0000, 0x0000};
GdkColor WHITE 		= {1, 0xffff, 0xffff, 0xffff};
GdkColor RED   		= {2, 0xffff, 0x0000, 0x0000};
GdkColor GREEN 		= {3, 0x0000, 0xffff, 0x0000};
GdkColor BLUE  		= {4, 0x0000, 0x0000, 0xffff};
GdkColor LIGHT_GREY = {5, 0xbbbb, 0xbbbb, 0xbbbb};
GdkColor MID_GREY 	= {6, 0x7777, 0x7777, 0x7777};
GdkColor DARK_GREY 	= {7, 0x3333, 0x3333, 0x3333};
GdkColor CYAN 		= {8, 0x0000, 0xffff, 0xffff};
GdkColor MAGENTA 	= {9, 0xffff, 0x0000, 0xffff};
GdkColor YELLOW 	= {10, 0xffff, 0xffff, 0x0000};
GdkColor ORANGE		= {11, 0xffff, 0x7777, 0x0000};
GdkColor NAVY_BLUE	= {12, 0x0000, 0x0000, 0xaaaa};
GdkColor OLIVE_GREEN = {13, 0x3333, 0x7777, 0xaaaa};
GdkColor PURPLE		= {14, 0xaaaa, 0x0000, 0xffff};
GdkColor BROWN		= {15, 0x7777, 0x0000, 0x0000};

GdkPixmap *buffer = NULL;
GdkGC *BandWcontext = NULL;
GdkGC *BlueandWcontext = NULL;
GdkGC *Whitecontext = NULL;
GdkGC *Gridcontext = NULL;

gint true_width = 0;
gint true_height = 0;
gint user_width = 0;
gint user_height = 0;
gint user_origin_x = 0;
gint user_origin_y = 0;
gint margin = 10;

#define MAJ_TICK_LEN		8
#define MIN_TICK_LEN		3

/*  Class internal declarations:  */

static GtkWidgetClass *parent_class = NULL;

static void gtk_graph_class_init (GtkGraphClass *class);
static void gtk_graph_init (GtkGraph *graph);
static void gtk_graph_realize (GtkWidget *widget);
static void gtk_graph_draw (GtkWidget *widget, GdkRectangle *area);
static void gtk_graph_size_request (GtkWidget *widget, GtkRequisition *req);
static void gtk_graph_size_allocate (GtkWidget *widget, GtkAllocation *allocation);
static gint gtk_graph_expose (GtkWidget *widget, GdkEventExpose *event);
static void gtk_graph_destroy (GtkObject *object);
static void gtk_graph_create_pixmap (GtkGraph *graph);
static void gtk_graph_set_grid_clipping_rectangles(GtkGraph *graph);
static void gtk_graph_plot_axes (GtkGraph *graph);
static void gtk_graph_plot_axes_titles (GtkGraph *graph);
static void gtk_graph_plot_traces (GtkGraph *graph);
static void gtk_graph_plot_title (GtkGraph *graph);
static void gtk_graph_plot_legend(GtkGraph *graph);

/* static int rint(float f); - PE comments out 20111215 */

/**
 * gtk_graph_get_type:
 *
 * This function is needed by the class handler and should not normally be
 * used. Used to defined the GtkGraph class to GTK
 * 
 * Returns: a #guint.
 */
guint gtk_graph_get_type (void)
{
  static guint graph_type = 0;

    if (!graph_type)    /* --- If not created yet --- */
		{

        /* --- Create a graph_info object --- */
        GtkTypeInfo graph_info =
            {
			"GtkGraph",
			sizeof (GtkGraph),
			sizeof (GtkGraphClass),
			(GtkClassInitFunc) gtk_graph_class_init,
			(GtkObjectInitFunc) gtk_graph_init,
			NULL,
			NULL,
			NULL
            };
  
        /* --- Tell GTK about it - get a unique identifying key --- */
        graph_type = gtk_type_unique (gtk_widget_get_type (), &graph_info);
    }
    return graph_type;
}


/*
 * gtk_graph_class_init
 *
 * Override any methods for the graph class that are needed for
 * the graph class to behave properly.  Here, the functions that
 * cause painting to occur are overridden.
 *
 * class - object definition class.
 */
static void gtk_graph_class_init (GtkGraphClass *class)
{
    GtkObjectClass *object_class;
    GtkWidgetClass *widget_class;

    /* --- Get the widget class --- */
    object_class = (GtkObjectClass *) class;
    widget_class = (GtkWidgetClass *) class;
    parent_class = gtk_type_class (gtk_widget_get_type ());

    /* --- Override object destroy --- */
    object_class->destroy = gtk_graph_destroy;

    /* --- Override these methods --- */
    widget_class->realize = gtk_graph_realize;
    widget_class->size_request = gtk_graph_size_request;
    widget_class->size_allocate = gtk_graph_size_allocate;
    widget_class->expose_event = gtk_graph_expose;
    widget_class->size_allocate = gtk_graph_size_allocate;
}

/* gtk_graph_init: Called each time a graph item gets created. This initializes fields in our structure. */
static void gtk_graph_init (GtkGraph *graph)
{
    gdk_color_alloc(gdk_colormap_get_system(), &WHITE);
    gdk_color_alloc(gdk_colormap_get_system(), &RED);
    gdk_color_alloc(gdk_colormap_get_system(), &GREEN);
    gdk_color_alloc(gdk_colormap_get_system(), &BLUE);
    gdk_color_alloc(gdk_colormap_get_system(), &LIGHT_GREY);
	gdk_color_alloc(gdk_colormap_get_system(), &MID_GREY);
	gdk_color_alloc(gdk_colormap_get_system(), &DARK_GREY);
	gdk_color_alloc(gdk_colormap_get_system(), &CYAN);
	gdk_color_alloc(gdk_colormap_get_system(), &MAGENTA);
	gdk_color_alloc(gdk_colormap_get_system(), &YELLOW);
	gdk_color_alloc(gdk_colormap_get_system(), &ORANGE);
	gdk_color_alloc(gdk_colormap_get_system(), &NAVY_BLUE);
	gdk_color_alloc(gdk_colormap_get_system(), &OLIVE_GREEN);
	gdk_color_alloc(gdk_colormap_get_system(), &PURPLE);
	gdk_color_alloc(gdk_colormap_get_system(), &BROWN);

}

/**
 * gtk_graph_new:
 * @type:  the type of #GtkGraph that you wish to create
 *
 * Create a new GtkGraph of type @type.  Currently this can be either
 * XY, POLAR or SMITH
 * 
 * Returns: a pointer of the form #GtkWidget to the newly created #GtkGraph.
 */

GtkWidget* gtk_graph_new (GtkGraphType type)
{
  GtkWidget *widget;
  GtkGraph *graph;
	
  widget = gtk_type_new (gtk_graph_get_type ());
  graph = GTK_GRAPH(widget);
  graph->graph_type = type;
  graph->traces = NULL;
  graph->num_traces = 0;
  graph->annotations = NULL;
  graph->num_annotations = 0;
  graph->smith_Z0 = 50.0;

  if (type == POLAR)
	{
	graph->polar_format.type = DEGREES_SYMMETRIC;
	graph->polar_format.polar_start = 60;
	}

  graph->dependant = (GtkGraphAxis *) malloc (sizeof(GtkGraphAxis));
  graph->dependant->autoscale_limits = TRUE;
  graph->dependant->autoscale_tick = TRUE;
  graph->dependant->crossing_type = GTK_GRAPH_AXISMIN;
  graph->dependant->crossing_value = 0;
  graph->dependant->precision = 0;
  graph->dependant->format = FLOATING_POINT;
  graph->dependant->grid_visible = FALSE;
  graph->dependant->title = NULL;

  graph->independant = (GtkGraphAxis *) malloc (sizeof(GtkGraphAxis));
  graph->independant->crossing_value = 0;  
  graph->independant->crossing_type = GTK_GRAPH_AXISMIN;
  graph->independant->autoscale_limits = TRUE; 
  graph->independant->autoscale_tick = TRUE; 
  graph->independant->precision = 0;
  graph->independant->format = FLOATING_POINT; 
  graph->independant->grid_visible = FALSE;
  graph->independant->title = NULL;

  graph->title = NULL;
  graph->subtitle = NULL;
  graph->legend_visible = TRUE;
  graph->legend_position = GTK_GRAPH_NORTH_WEST;

  return GTK_WIDGET(graph);
}


/* gtk_graph_realize: Associate the widget with an x-window.*/
static void gtk_graph_realize (GtkWidget *widget)
{
GtkGraph *darea;
GdkWindowAttr attributes;
gint attributes_mask;

  /* --- Check for failures --- */
g_return_if_fail (widget != NULL);
g_return_if_fail (GTK_IS_GRAPH (widget));

darea = GTK_GRAPH (widget);
GTK_WIDGET_SET_FLAGS (widget, GTK_REALIZED);

/* --- attributes to create the window --- */
attributes.window_type = GDK_WINDOW_CHILD;
attributes.x = widget->allocation.x;
attributes.y = widget->allocation.y;
attributes.width = widget->allocation.width;
attributes.height = widget->allocation.height;
attributes.wclass = GDK_INPUT_OUTPUT;
attributes.visual = gtk_widget_get_visual (widget);
attributes.colormap = gtk_widget_get_colormap (widget);
attributes.event_mask = gtk_widget_get_events (widget) | GDK_EXPOSURE_MASK;

/* --- We're passing in x, y, visual and colormap values --- */
attributes_mask = GDK_WA_X | GDK_WA_Y | GDK_WA_VISUAL | GDK_WA_COLORMAP;

/* --- Create the window --- */
widget->window = gdk_window_new (gtk_widget_get_parent_window (widget), &attributes, attributes_mask);
gdk_window_set_user_data (widget->window, darea);

widget->style = gtk_style_attach (widget->style, widget->window);
gtk_style_set_background (widget->style, widget->window, GTK_STATE_NORMAL);

if (BlueandWcontext)
	gdk_gc_unref (BlueandWcontext);
BlueandWcontext = gdk_gc_new(widget->window);
gdk_gc_set_foreground(BlueandWcontext, &BLUE);
gdk_gc_set_background(BlueandWcontext, &WHITE);
    
if (BandWcontext)
  	gdk_gc_unref (BandWcontext);
BandWcontext = gdk_gc_new(widget->window);
gdk_gc_set_background(BandWcontext, &WHITE);
gdk_gc_set_foreground(BandWcontext, &BLACK);

if (Whitecontext)
	gdk_gc_unref (Whitecontext);
Whitecontext = gdk_gc_new(widget->window);
gdk_gc_set_foreground(Whitecontext, &WHITE);
gdk_gc_set_background(Whitecontext, &WHITE);
  
if (Gridcontext)
	gdk_gc_unref (Gridcontext);
Gridcontext = gdk_gc_new(widget->window);
gdk_gc_set_foreground(Gridcontext, &LIGHT_GREY);
gdk_gc_set_background(Gridcontext, &WHITE);

gtk_graph_create_pixmap (darea);
}

/* gtk_graph_draw: Draw the widget.*/
/* The philosophy for drawing a cartesian graph is as follows
*	
* 1) Clear the drawing area workspace.  This is accomplished by drawing a filled white
* rectangle over the entire drawing area, and then drawing an unfilled blue rectangle on
* top of it as a border
* 2) Plot the title.  Since either or both of  title and subtitle may be null it makes
* little sense to reserve part of the drawing area for them.  Instead draw them first, and
* if the height of the text they contain is non-zero then within gtk_graph_plot_title()
* adjust the user_origin_y variable to account for the space used.
* 3) Plot the axis labels.  The Y-axis is easiest to deal with.  The label for this axis
* will always be beneath the subtitle, and thus if height of the Y-axis label text is
* non-zero user_origin_y is adjusted accordingly.
* 
*/

static void gtk_graph_draw (GtkWidget *widget, GdkRectangle *area)
{
    GtkGraph *graph = GTK_GRAPH (widget);

    g_return_if_fail (widget != NULL);  /* --- Check for obvious problems --- */
    g_return_if_fail (GTK_IS_GRAPH (widget));
	g_return_if_fail (GTK_WIDGET_REALIZED(widget));
               
    // Clear the graph first 
  	
 	gdk_draw_rectangle(buffer, Whitecontext, TRUE, 0, 0, true_width-1, true_height-1);
   	gdk_draw_rectangle(buffer, BlueandWcontext, FALSE, 0, 0, true_width-1, true_height-1);

 	gtk_graph_plot_title (graph);   /* Title First - since this may affect the screen area available
									* to draw on an hence it will adjust user_origin_y and user_height */
	
	gtk_graph_plot_axes_titles(graph); /* user_origin_y and user_height can also be changed here */
  	
	/* all user_(origin_x, origin_y, height and width) should now be fixed and we can now
	*  scale the axes and set the clipping rectange for the grid */

	gtk_graph_axis_scale_axis (graph, user_width, user_height);  /* Scale the axes - this also determines the screen scaling factors */
	gtk_graph_set_grid_clipping_rectangles(graph);
	
	switch (graph->graph_type)
		{
		case POLAR:
			gtk_graph_polar_plot_axes (graph);   /* Now plot the axes to screen, */
			gtk_graph_polar_plot_traces (graph);
			break;
		case SMITH:
			gtk_graph_smith_plot_axes(graph);
			gtk_graph_smith_plot_traces(graph);
			break;
		default:
			
			gtk_graph_plot_axes (graph);   /* Now plot the axes to screen, */
			gtk_graph_plot_traces (graph);  /* then the data traces themselves */
			break;
		}
	gtk_graph_plot_annotations(graph);
	gtk_graph_plot_legend (graph); /* and finally the legend to describe the trace */
		
    /* Once all relevant information has been transfered into the buffer then copy to screen */
    gdk_draw_pixmap(widget->window, BandWcontext, buffer,0, 0, 0, 0, true_width, true_height);
}

/* gtk_graph_size_request: How big should the widget be?  */
static void gtk_graph_size_request (GtkWidget *widget, GtkRequisition *req)
{
    GtkGraph *graph = GTK_GRAPH (widget);

    g_return_if_fail (widget != NULL);  /* --- Check for obvious problems --- */
    g_return_if_fail (GTK_IS_GRAPH (widget));

	if (graph->graph_type == XY)
		{
		req->width = 550;
		req->height = 400;
		}
	else
		{
		req->width = 380;
		req->height = 480;
		}
}

/* gtk_graph_expose: The graph widget has been exposed and needs to be painted.*/
static void gtk_graph_size_allocate (GtkWidget *widget, GtkAllocation *allocation)
{
gint min = allocation->width;
g_return_if_fail (widget != NULL);
g_return_if_fail (GTK_IS_GRAPH (widget));
g_return_if_fail (allocation != NULL);

GtkGraph *graph = GTK_GRAPH(widget);

if (graph->graph_type != XY) // Needs to be a circular graph, hence drawn within a square.
	{
	min = allocation->width;
	if (allocation->height < allocation->width)
		min = allocation ->height;
	allocation->height = min;
	allocation->width = min;
	}

if (GTK_WIDGET_REALIZED (widget)) // Then we are moving / resizing an existing widget
 	gdk_window_move_resize (widget->window, allocation->x, allocation->y, allocation->width, allocation->height);

widget->allocation = *allocation;
gtk_graph_create_pixmap (GTK_GRAPH (widget));
} 
 
 
static gint gtk_graph_expose (GtkWidget *widget, GdkEventExpose *event)
{
    /* --- Do error checking --- */
    g_return_val_if_fail (widget != NULL, FALSE);
    g_return_val_if_fail (GTK_IS_GRAPH (widget), FALSE);
    g_return_val_if_fail (event != NULL, FALSE);

    if (event->count > 0)
        return (FALSE);

    gtk_graph_create_pixmap (GTK_GRAPH(widget));
    gtk_graph_draw (widget, NULL);    /* --- Draw the graph --- */
    return (FALSE);
}

static void gtk_graph_destroy (GtkObject *object)
{
    GtkGraph *graph;

    /* --- Do error checking --- */
    g_return_if_fail (object != NULL);
    g_return_if_fail (GTK_IS_GRAPH (object));

    graph = GTK_GRAPH (object);    /* --- Convert to graph object --- */
    g_free (graph->dependant);
    g_free (graph->independant);
    graph->dependant = NULL;
    graph->independant = NULL;
	

    /* --- Call parent destroy --- */
    GTK_OBJECT_CLASS (parent_class)->destroy (object);
}


static void gtk_graph_create_pixmap (GtkGraph *graph)
{
GtkWidget *widget;

g_return_if_fail (graph != NULL);
g_return_if_fail (GTK_IS_GRAPH (graph));

  if (GTK_WIDGET_REALIZED (graph))
    {
    widget = GTK_WIDGET (graph);

    if (buffer)
	   gdk_pixmap_unref (buffer);

    buffer = gdk_pixmap_new (widget->window, widget->allocation.width, widget->allocation.height, -1);
	true_width = widget->allocation.width;
    true_height = widget->allocation.height;

	if (graph->graph_type == XY)
		{
		user_origin_x = margin;
		user_origin_y = margin / 2.0;
		user_height = true_height - user_origin_y - margin;
		user_width = true_width - user_origin_x - 2.0 * margin;
		}
	else
		{	
		user_origin_x = 2 * margin;
		user_origin_y = 3 * margin;
		user_height = true_height - 5.0 * margin;
		user_width = true_width - 4.0 * margin;
		}

	}
}

static void gtk_graph_set_grid_clipping_rectangles(GtkGraph *graph)
{
GdkRectangle rect;
GtkGraphTrace *tmp;
gint n;

rect.x = user_origin_x;
rect.y = user_origin_y;
rect.width = user_width+1;
rect.height = user_height+1;
gdk_gc_set_clip_origin(Gridcontext, 0, 0);
gdk_gc_set_clip_rectangle(Gridcontext, &rect);

/* These settings are intial guesses.  The exact values will be known once the axes
*  have been scaled */
tmp = graph->traces;
for (n = 0 ; n < graph->num_traces ; n++)
	{
	gdk_gc_set_clip_origin(tmp->format->line_gc, 0, 0);		/* Set the clipping for each trace */
	gdk_gc_set_clip_rectangle(tmp->format->line_gc, &rect);
	gdk_gc_set_clip_origin(tmp->format->marker_gc, 0, 0);		/* Set the clipping for each trace */
	gdk_gc_set_clip_rectangle(tmp->format->marker_gc, &rect);
	tmp = tmp->next;	/* and then move onto the next trace */
	}
	
}
/* gtk_graph_plot_axes: Draw the graph axes on the pixmap.*/
static void gtk_graph_plot_axes (GtkGraph *graph)
{
gint x_baseline, y_baseline, x_coord, y_coord;
gint i, j, text_width, text_height;
gchar *tbuffer, *format_string = NULL;
gfloat min_tick_value;

PangoFontDescription *fontdesc = NULL;
PangoLayout *layout = NULL;

g_return_if_fail (graph != NULL);
g_return_if_fail (GTK_IS_GRAPH (graph));
g_return_if_fail (GTK_WIDGET_REALIZED (graph));
	
fontdesc = pango_font_description_from_string("Sans 10");
layout = gtk_widget_create_pango_layout(GTK_WIDGET(graph), NULL);
pango_layout_set_font_description(layout, fontdesc);	

/* Draw the X axis */

/* Check at which value on the Y-axis wants we should have the X-axis cross */

switch (graph->dependant->crossing_type)
    {
    case GTK_GRAPH_AXISMAX:
        x_baseline = user_origin_y;
        break;
    case GTK_GRAPH_AXISMIN:
        x_baseline = user_origin_y + graph->dependant->n_maj_tick * graph->dependant->pxls_per_maj_tick;
        break;
    case GTK_GRAPH_USERVALUE:
    default:
        if (graph->dependant->crossing_value < graph->dependant->axis_min)
            x_baseline = user_origin_y + graph->dependant->n_maj_tick * graph->dependant->pxls_per_maj_tick;
        else if (graph->dependant->crossing_value > graph->dependant->axis_max)
	        x_baseline = user_origin_y;
        else
        	x_baseline = (int) rint((graph->dependant->axis_max - graph->dependant->crossing_value) * graph->dependant->scale_factor) + user_origin_y;
        break;
    }        	

for (i = 0 ; i <= graph->independant->n_maj_tick ; i++)
	{
	x_coord = user_origin_x -1 + i * graph->independant->pxls_per_maj_tick;
	switch(graph->independant->format)
		{
		case ENGINEERING:
			format_string = g_strdup_printf("%%.%de", graph->independant->precision);
			break;
		case SCIENTIFIC:
			format_string = g_strdup_printf("%%.%de", graph->independant->precision);
			break;
		case FLOATING_POINT:
		default:
			format_string = g_strdup_printf("%%.%df", graph->independant->precision);
			break;
		}
	tbuffer = g_strdup_printf(format_string, graph->independant->axis_min + i*graph->independant->maj_tick );
	g_free(format_string);
	if (graph->independant->grid_visible)
		gdk_draw_line(buffer, Gridcontext, x_coord, user_origin_y, x_coord, user_origin_y + graph->dependant->pxls_per_maj_tick * graph->dependant->n_maj_tick);
	gdk_draw_line(buffer, BandWcontext, x_coord, x_baseline, x_coord, x_baseline + MAJ_TICK_LEN);
	pango_layout_set_text(layout, tbuffer, -1);
	g_free(tbuffer);
	pango_layout_get_pixel_size(layout, &text_width, &text_height);
	gdk_draw_layout(buffer, BandWcontext, x_coord - text_width / 2, x_baseline + MAJ_TICK_LEN + 2, layout);
	for (j = 1 ; j <= graph->independant->n_min_tick ; j++)
		{
		if (i != graph->independant->n_maj_tick)
			{
			min_tick_value = (gfloat) j * (gfloat) graph->independant->pxls_per_maj_tick / (gfloat) (graph->independant->n_min_tick );
			x_coord = user_origin_x -1+ i * graph->independant->pxls_per_maj_tick + (gint) min_tick_value;
			gdk_draw_line(buffer, BandWcontext, x_coord, x_baseline, x_coord, x_baseline + MIN_TICK_LEN);
			}
		}
	}

/* Now Draw the Y axis */
/* Check at which value on the X-axis wants we should have the Y-axis cross */

switch (graph->independant->crossing_type)
    {
    case GTK_GRAPH_AXISMAX:
       		y_baseline = user_origin_x + user_width;
        break;
    case GTK_GRAPH_AXISMIN:
        	y_baseline = user_origin_x-1;
        break;
    case GTK_GRAPH_USERVALUE:
    default:
        if (graph->independant->crossing_value < graph->independant->axis_min)
        	y_baseline = user_origin_x-1;
       	else if (graph->independant->crossing_value > graph->independant->axis_max)
       		y_baseline = user_origin_x + user_width;
  		else
  			y_baseline = (int) rint((graph->independant->crossing_value - graph->independant->axis_min) * graph->independant->scale_factor) + user_origin_x-1;
        break;
    }        	

for (i = 0 ; i <= graph->dependant->n_maj_tick ; i++)
	{
	y_coord = user_origin_y + i * graph->dependant->pxls_per_maj_tick;
	switch(graph->dependant->format)
		{
		case ENGINEERING:
			format_string = g_strdup_printf("%%.%de", graph->dependant->precision);
			break;
		case SCIENTIFIC:
			format_string = g_strdup_printf("%%.%de", graph->dependant->precision);
			break;
		case FLOATING_POINT:
		default:
			format_string = g_strdup_printf( "%%.%df", graph->dependant->precision);
			break;
		}
	tbuffer = g_strdup_printf(format_string, graph->dependant->axis_max - i*graph->dependant->maj_tick );
	g_free(format_string);
	if (graph->dependant->grid_visible)
		gdk_draw_line(buffer, Gridcontext, y_baseline-1, y_coord, y_baseline + graph->independant->pxls_per_maj_tick * graph->independant->n_maj_tick, y_coord);
	gdk_draw_line(buffer, BandWcontext, y_baseline-1, y_coord, y_baseline - MAJ_TICK_LEN, y_coord);
	pango_layout_set_text(layout, tbuffer, -1);
	g_free(tbuffer);
	pango_layout_get_pixel_size(layout, &text_width, &text_height);
	gdk_draw_layout(buffer, BandWcontext, y_baseline - MAJ_TICK_LEN - 2 - text_width, y_coord - text_height / 2, layout);
	for (j = 1 ; j <= graph->dependant->n_min_tick ; j++)
		{
		if (i != graph->dependant->n_maj_tick)
			{
			min_tick_value = (gfloat) j * (gfloat) graph->dependant->pxls_per_maj_tick / (gfloat) (graph->dependant->n_min_tick + 1);
			y_coord = user_origin_y + i * graph->dependant->pxls_per_maj_tick + min_tick_value;
			gdk_draw_line(buffer, BandWcontext, y_baseline-1, y_coord, y_baseline - MIN_TICK_LEN, y_coord);
			}
		}
	}

gdk_draw_line (buffer, BandWcontext, user_origin_x-1, x_baseline, user_origin_x-1  + graph->independant->pxls_per_maj_tick * graph->independant->n_maj_tick, x_baseline);
gdk_draw_line (buffer, BandWcontext, y_baseline, user_origin_y, y_baseline, user_origin_y + graph->dependant->pxls_per_maj_tick * graph->dependant->n_maj_tick);
pango_font_description_free(fontdesc);
g_object_unref(layout);
}
/* gtk_graph_plot_axes_titles */
static void gtk_graph_plot_axes_titles (GtkGraph *graph)
{
gint text_width, text_height, y_range, label_len;
PangoFontDescription *fontdesc = NULL;
PangoLayout *layout = NULL;
gint i;
gchar *label_buffer = NULL, *format_string = NULL;
gfloat global_X_max = -1E99, global_Y_max= -1E99, global_X_min=1E99, global_Y_min=1E99;
GtkGraphTrace *t = graph->traces;
g_return_if_fail (graph != NULL);
g_return_if_fail (GTK_IS_GRAPH (graph));
g_return_if_fail (GTK_WIDGET_REALIZED (graph));
	
fontdesc = pango_font_description_from_string("Sans 10");
layout = gtk_widget_create_pango_layout(GTK_WIDGET(graph), NULL);
pango_layout_set_font_description(layout, fontdesc);		

/* Dichotomy: we need axis_max and axis_min to determine values of crossing points and
*  positions of maximum label length but we can't determine them till after we know the
*  user_width / user_height from this function.  Hence must locally find axis_max/min */

for (i = 0 ; i < graph->num_traces ; i++)
	{
	if (t->xmax > global_X_max)
		global_X_max = t->xmax;
	if (t->xmin < global_X_min)
		global_X_min = t->xmin;
	if (t->ymax > global_Y_max)
		global_Y_max = t->ymax;
	if (t->ymin < global_Y_min)
		global_Y_min = t->ymin;
	t = t->next;
	}
	
if (graph->dependant->title != NULL)
	{
	pango_layout_set_text(layout, graph->dependant->title, -1);
	pango_layout_get_pixel_size(layout, &text_width, &text_height);
	gdk_draw_layout(buffer, BandWcontext, margin, user_origin_y, layout);
	user_origin_y += (text_height + margin / 2.0);
	}
/* Nothing should affect user_origin_y from here, so thus calculate user_height from it */
	
user_height = true_height - user_origin_y - margin / 2.0;

if (graph->independant->title != NULL)
	{
	pango_layout_set_text(layout, graph->independant->title, -1);
	pango_layout_get_pixel_size(layout, &text_width, &text_height);
	gdk_draw_layout(buffer, BandWcontext, user_origin_x + user_width / 2 - text_width / 2, user_origin_y + user_height - text_height, layout);
	user_height -= text_height;		
	}

/* Calculate the height of the X-axis label text Since this may need to be subtracted
	from the user_height variable if the X-axis crosses the Y-axis somewhere the minumum */
label_buffer = g_strdup_printf("1234567890E+1");
pango_layout_set_text(layout, label_buffer, -1);
pango_layout_get_pixel_size(layout, &text_width, &text_height);
g_free(label_buffer);

/* If the crossing_value is within a specified value (10% of the total range) of the axis_max
	minumum then we need to reduce user_height to account for the text height.  If the crossing_value is more than 5%
	greater than the minumum, then we have no need to worry as the axis labels will be within the area of the main graph*/

y_range = global_Y_max - global_Y_min;
if ((graph->dependant->crossing_type == GTK_GRAPH_AXISMIN) || ((graph->dependant->crossing_type == GTK_GRAPH_USERVALUE) && (graph->dependant->crossing_value < global_Y_min + 0.1 * y_range)))
	user_height -= (text_height + MAJ_TICK_LEN);

if (graph->graph_type != XY)	/* Polar Plots and Smiths Charts also have additional labels */
	user_height -= text_height;	/*outside the unit circle that needs to be accounted for */

/* Thus far we have considered factors affecting y-origin and height only, but the length of
text labels on the Y-axis affects the position of the x-origin and width, so now conduct
the same sort of exercise as above. The maximum width of the Y-axis label text will occur
at the axis maximum (or possibly the minimum due to the extra length of a minus sign */
switch(graph->dependant->format)
	{
	case FLOATING_POINT:
		format_string = g_strdup_printf( "%%.%df", graph->dependant->precision);
		break;
	case ENGINEERING:
		format_string = g_strdup_printf("%%.%de", graph->dependant->precision);
		break;
	case SCIENTIFIC:
		format_string = g_strdup_printf("%%.%de", graph->dependant->precision);
		break;
	}
label_buffer = g_strdup_printf(format_string, global_Y_max);
pango_layout_set_text(layout, label_buffer, -1);
pango_layout_get_pixel_size(layout, &text_width, &text_height);
g_free(label_buffer);
label_len = text_width;
	
label_buffer = g_strdup_printf(format_string, global_Y_min);

pango_layout_set_text(layout, label_buffer, -1);
pango_layout_get_pixel_size(layout, &text_width, &text_height);	
g_free(label_buffer);	
g_free(format_string);	

user_origin_x += (MAX(text_width, label_len) + MAJ_TICK_LEN);
user_width -= (MAX(text_width, label_len) + MAJ_TICK_LEN);
	
pango_font_description_free(fontdesc);
g_object_unref(layout);	
}

/* gtk_graph_plot_traces: Draw the graph traces on the pixmap.*/
static void gtk_graph_plot_traces (GtkGraph *graph)
{
gint i = 0, n;
GtkGraphTrace *tmp;
GdkPoint *pts = NULL;
GdkRectangle clip_rect;

g_return_if_fail (graph != NULL);
g_return_if_fail (GTK_IS_GRAPH (graph));
g_return_if_fail (GTK_WIDGET_REALIZED (graph));

tmp = graph->traces;
for (n = 0 ; n < graph->num_traces ; n++)
	{
	/* Make sure that there is some data in the trace */
	if (tmp->Xdata == NULL || tmp->Ydata == NULL) 
		continue;

	/* Assign the storage for the co-ordinates of each data point */

	pts = (GdkPoint *) g_malloc ((tmp->num_points) * sizeof(GdkPoint));

	/* Set a clipping rectangle for each trace */
	
	clip_rect.x = user_origin_x;
	clip_rect.y = user_origin_y;
	clip_rect.width = user_origin_x + graph->independant->pxls_per_maj_tick * graph->independant->n_maj_tick;
	clip_rect.height = user_origin_y + graph->dependant->pxls_per_maj_tick * graph->dependant->n_maj_tick;
	gdk_gc_set_clip_origin(tmp->format->line_gc, 0, 0);
	gdk_gc_set_clip_rectangle(tmp->format->line_gc, &clip_rect);

	/* Calculate the positions of the co-ordinates */
	
	for (i = 0 ; i < tmp->num_points ; i++)
		{
		pts[i].x = (int) rint((tmp->Xdata[i] - graph->independant->axis_min ) * graph->independant->scale_factor) + user_origin_x;
		pts[i].y = (int) rint((graph->dependant->axis_max - tmp->Ydata[i] ) * graph->dependant->scale_factor) + user_origin_y;
		if (tmp->Ydata[i] < graph->dependant->axis_min)  /* Implement Crude bottom end clipping */
			pts[i].y = (int) rint((graph->dependant->axis_max - graph->dependant->axis_min) * graph->dependant->scale_factor) + user_origin_y;
			
		}

	gdk_draw_lines (buffer,tmp->format->line_gc, pts, tmp->num_points);/* Draw the lines */	

	for (i = 0 ; i < tmp->num_points ; i++)/* and then draw the markers */	
        if (tmp->format->marker_type != GTK_GRAPH_MARKER_NONE)
			{
			gdk_gc_set_clip_origin(tmp->format->marker_gc, pts[i].x - tmp->format->marker_size, pts[i].y - tmp->format->marker_size);
            gdk_draw_pixmap(buffer, tmp->format->marker_gc, tmp->format->marker, 0, 0, pts[i].x - tmp->format->marker_size, pts[i].y - tmp->format->marker_size, -1, -1);
			}

	/* Free up all storage after use */
	g_free(pts);
	
	/* and then move onto the next trace */
	tmp = tmp->next;
	}
}

static void gtk_graph_plot_title (GtkGraph *graph)
{
gint text_width, text_height;	
PangoFontDescription *fontdesc = NULL;
PangoLayout *layout = NULL;

g_return_if_fail (graph != NULL);
g_return_if_fail (GTK_IS_GRAPH (graph));
g_return_if_fail (GTK_WIDGET_REALIZED (graph));

if (graph->title == NULL)
	return;

fontdesc = pango_font_description_from_string("Sans 11");
layout = gtk_widget_create_pango_layout(GTK_WIDGET(graph), graph->title);
pango_layout_set_font_description(layout, fontdesc);
pango_layout_get_pixel_size(layout, &text_width, &text_height);

gdk_draw_layout(buffer, BandWcontext, true_width / 2.0 - text_width / 2.0, user_origin_y , layout);
pango_font_description_free(fontdesc);

user_origin_y += text_height;
user_height -= text_height;

if (graph->subtitle == NULL)
	return;
	
fontdesc = pango_font_description_from_string("Sans 10");
pango_layout_set_text(layout, graph->subtitle, -1);
pango_layout_set_font_description(layout, fontdesc);	
pango_layout_get_pixel_size(layout, &text_width, &text_height);
gdk_draw_layout(buffer, BandWcontext, true_width / 2.0 - text_width / 2.0, user_origin_y, layout);

user_origin_y += text_height;
user_height -= text_height;

pango_font_description_free(fontdesc);
g_object_unref(layout);
}

/**
 * gtk_graph_set_title:
 * @graph:  the #GtkGraph that you wish to set titles for
 * @title:	the title to be added to @graph
 * @subtitle:	the subtitle to be added to @graph, and displayed underneath the title
 *
 * Adds a title and/or subtitle the #GtkGraph @graph.  Either or both may be NULL.
 */


void gtk_graph_set_title (GtkGraph *graph, const gchar *title, const gchar *subtitle)
{
g_return_if_fail (graph != NULL);
g_return_if_fail (GTK_IS_GRAPH (graph));

if (graph->title != NULL)
    free(graph->title);
if (graph->subtitle != NULL)
    free(graph->subtitle);

graph->title = g_strdup(title);
graph->subtitle = g_strdup(subtitle);
}
  

static void gtk_graph_plot_legend (GtkGraph *graph)                  
{
gint text_height, text_width;	
PangoFontDescription *fontdesc = NULL;
PangoLayout *layout = NULL;

gint i, ypos = 0, max_len = 0, len = 0;
gint legend_x = 0, legend_y = 0, legend_width = 0, legend_height = 0, legend_margin;
gchar *text_buffer = NULL;
GtkGraphTrace *tmp;

g_return_if_fail (graph != NULL);
g_return_if_fail (GTK_IS_GRAPH (graph));
g_return_if_fail (GTK_WIDGET_REALIZED (graph));
if (graph->legend_visible == FALSE)
	return;

fontdesc = pango_font_description_from_string("Sans 8");
layout = gtk_widget_create_pango_layout(GTK_WIDGET(graph), NULL);
pango_layout_set_font_description(layout, fontdesc);

legend_margin = 4;
	
tmp = graph->traces;
for (i = 0 ; i < graph->num_traces ; i++)
	{
	if (tmp->format->legend_text == NULL)
		{
		pango_layout_set_text(layout, "Trace 00", -1);
		pango_layout_get_pixel_size(layout, &len, &text_height);
		}
	else
		{
		pango_layout_set_text(layout, tmp->format->legend_text, -1);
		pango_layout_get_pixel_size(layout, &text_width, &text_height);
		}
	if (text_width > max_len)
		max_len = text_width;
	}

if (graph->legend_position != GTK_GRAPH_NORTH && graph->legend_position != GTK_GRAPH_SOUTH)
    {
	legend_width = max_len + 40;
	legend_height = graph->num_traces * (text_height + legend_margin);
    }
else
	{
	}
    
switch (graph->legend_position)
    {
    case GTK_GRAPH_NORTH:
        break;
    case GTK_GRAPH_NORTH_EAST:
        legend_x = user_origin_x + user_width - legend_width - legend_margin;
        legend_y = user_origin_y + legend_margin;
        break;
    case GTK_GRAPH_EAST:
        legend_x = user_origin_x + user_width - legend_width - legend_margin;
        legend_y = user_origin_y + user_height / 2 - legend_height / 2;
        break;
    case GTK_GRAPH_SOUTH_EAST:
        legend_x = user_origin_x + user_width - legend_width - legend_margin;
		if (graph->graph_type == XY)	
			legend_y = user_origin_y + graph->dependant->n_maj_tick * graph->dependant->pxls_per_maj_tick - legend_height;
		else
			legend_y = user_origin_y + 2.0 * graph->dependant->n_maj_tick * graph->dependant->pxls_per_maj_tick - legend_height;			
        break;
    case GTK_GRAPH_SOUTH:
        legend_x = true_width / 2 - legend_width;
		if (graph->graph_type == XY)
			legend_y = user_origin_y + graph->dependant->n_maj_tick * graph->dependant->pxls_per_maj_tick - legend_height;
		else
			legend_y = user_origin_y + 2.0 * graph->dependant->n_maj_tick * graph->dependant->pxls_per_maj_tick - legend_height;
        break;      
    case GTK_GRAPH_SOUTH_WEST:
        legend_x = user_origin_x + legend_margin;
		if (graph->graph_type == XY)
			legend_y = user_origin_y + graph->dependant->n_maj_tick * graph->dependant->pxls_per_maj_tick - legend_height - legend_margin;
		else
			legend_y = user_origin_y + 2.0 * graph->dependant->n_maj_tick * graph->dependant->pxls_per_maj_tick - legend_height;
        break;
    case GTK_GRAPH_WEST:
        legend_x = user_origin_x + legend_margin;
        legend_y = user_origin_y + user_height / 2 - legend_height/2;
        break;
    case GTK_GRAPH_NORTH_WEST:
        legend_x = user_origin_x + legend_margin;
        legend_y = user_origin_y + legend_margin;
        break;
    }

gdk_draw_rectangle(buffer, Whitecontext, TRUE, legend_x, legend_y, legend_width, legend_height);
gdk_draw_rectangle(buffer, BandWcontext, FALSE, legend_x, legend_y, legend_width, legend_height);

tmp = graph->traces;
for (i = 0 ; i < graph->num_traces ; i++)
	{
    ypos =	(legend_margin + text_height) * (2 * i + 1) / 2  ;

	gdk_draw_line(buffer, tmp->format->line_gc, legend_x + legend_margin, legend_y + ypos, legend_x + 25, legend_y + ypos);
	if (tmp->format->marker_type != GTK_GRAPH_MARKER_NONE)
		{
		gdk_gc_set_clip_origin(tmp->format->marker_gc, legend_x + 15 - tmp->format->marker_size, legend_y + ypos - tmp->format->marker_size);
        gdk_draw_pixmap(buffer, tmp->format->marker_gc, tmp->format->marker, 0, 0, legend_x + 15 - tmp->format->marker_size, legend_y + ypos - tmp->format->marker_size, -1, -1);
		}
	
	if  (text_buffer != NULL);
		free(text_buffer);

	if (tmp->format->legend_text == NULL)
	  text_buffer = g_strdup_printf("Trace %d", i);
	else
	  text_buffer = g_strdup_printf(tmp->format->legend_text, "%s"); /* add a format arg --PE 20111215 */

	pango_layout_set_text(layout, text_buffer, -1);
	pango_layout_get_pixel_size(layout, &len, &text_height);
	gdk_draw_layout(buffer, BandWcontext, legend_x + 35, legend_y + ypos - text_height /2, layout);	
	
	tmp=tmp->next;
	}
pango_font_description_free(fontdesc);
g_object_unref(layout);
}

/**
 * gtk_graph_legend_format:
 * @graph:  the #GtkGraph that contains the legend you wish to format
 * @is_visible:	set legend visibility to TRUE or FALSE
 * @position:
 *
 * Format the legend within the #GtkGraph @graph
 */

void gtk_graph_legend_format(GtkGraph *graph, gboolean is_visible, GtkGraphPosition position)
{
g_return_if_fail (graph != NULL);
g_return_if_fail (GTK_IS_GRAPH (graph));

graph->legend_visible = is_visible;
graph->legend_position = position;
}

/* libm.a version of rint() doesn't seem to link with dev-cpp so here's my own */

/* static int rint(float f) */
/* { */
/* int temp_i; */
/* float temp_f; */

/* temp_i = (int) f; */
/* temp_f = f - (float) temp_i; */

/* if (temp_f >= 0.5) */
/*     return temp_i + 1; */
/* else */
/*     return temp_i; */
/* } */
