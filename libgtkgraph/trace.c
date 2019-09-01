#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "gtkgraph.h"
#include "gtkgraph_internal.h"


/* Declaration of all local functions */

GtkGraphTrace *gtk_graph_trace_allocate(void);

/* Externally referenceable functions */

GtkGraphTrace *gtk_graph_trace_allocate(void)
{
GtkGraphTrace *tmp;

tmp = (GtkGraphTrace *) g_malloc(sizeof(GtkGraphTrace));

if (tmp == (GtkGraphTrace *) NULL)
    {
        (void) fprintf(stderr,"malloc failed at: %s\n","gtk_graph_trace_allocate");
        return ((GtkGraphTrace *) NULL);
    }
tmp->Xdata = NULL;
tmp->Ydata = NULL;
tmp->num_points = 0;
tmp->xmax = 0;
tmp->xmin = 0;
tmp->ymax = 0;
tmp->ymin = 0;
tmp->format = (GtkGraphTraceFormat *) g_malloc (sizeof(GtkGraphTraceFormat));
tmp->next = NULL;

tmp->format->marker_type = GTK_GRAPH_MARKER_NONE;
tmp->format->marker = NULL;
tmp->format->marker_gc = NULL;
tmp->format->mask = NULL;
tmp->format->mask_gc = NULL;
tmp->format->marker_filled = FALSE;
tmp->format->marker_size = 5;
	
tmp->format->line_type = SOLID;
tmp->format->line_visible = TRUE;
tmp->format->line_gc = NULL;
	
tmp->format->legend_text = NULL;

return(tmp);
}

/**
 * gtk_graph_trace_new:
 * @graph:  the #GtkGraph that you wish to add a new trace to
 *
 * Creates a new data trace within the #GtkGraph @graph.  The integer returned
 * by this function must be stored within the program since it is used whenever
 * this trace needs to be referenced again (e.g when setting data with
 * #gtk_graph_trace_set_data
 * 
 * Returns: a unique integer identifier for the new trace
 */
gint gtk_graph_trace_new(GtkGraph *graph)
{

int i;
GtkGraphTrace *new_trace, *tmp;
GtkWidget *widget = NULL;

g_return_val_if_fail (graph != NULL, FALSE);
g_return_val_if_fail (GTK_IS_GRAPH (graph), FALSE);
g_return_val_if_fail (GTK_WIDGET_REALIZED(GTK_WIDGET(graph)), FALSE);

widget = GTK_WIDGET(graph);
new_trace = gtk_graph_trace_allocate();
	
if (new_trace->format->line_gc)
	gdk_gc_unref (new_trace->format->line_gc);
new_trace->format->line_gc = gdk_gc_new(widget->window);

if (new_trace->format->marker_gc)
	gdk_gc_unref (new_trace->format->marker_gc);
new_trace->format->marker_gc = gdk_gc_new(widget->window);

if (new_trace->format->mask_gc)
	gdk_gc_unref (new_trace->format->mask_gc);
new_trace->format->mask_gc = gdk_gc_new(widget->window);

gdk_gc_set_line_attributes(new_trace->format->line_gc, 0, GDK_LINE_SOLID, GDK_CAP_BUTT, GDK_JOIN_ROUND);

if (graph->traces == NULL)
	graph->traces = new_trace;
else
	{
	for (i = 0, tmp = graph->traces; tmp->next != NULL ; tmp = tmp->next , i++)  // Should advance us to last item in list
        ;
	tmp->next = new_trace;
	}
graph->num_traces += 1;
return graph->num_traces - 1;
}
/**
 * gtk_graph_trace_set_data:
 * @graph:  the #GtkGraph containing the trace to be modified
 * @trace_id: 	the unique identifier of the trace
 * @xd:	an array containing the x co-ordinates of each data point
 * @yd: an array containing the y co-ordinates of each data point
 * @n:	the number of data points
 *
 * Assigns to trace @trace_id the co-ordinate pairs stored within
 * the @xd and @yd data arrays.  Note @xd and @yd must be @n data
 * points long
 */
void gtk_graph_trace_set_data(GtkGraph *graph, gint trace_id, gfloat *xd, gfloat *yd, gint n)
{
GtkGraphTrace *t;
gint i;
gfloat xmax_store = -1.0E99, xmin_store = 1.0E99, ymax_store = -1.0E99, ymin_store = 1.0E99;

g_return_if_fail (graph != NULL);
g_return_if_fail (GTK_IS_GRAPH (graph));
g_return_if_fail (graph->traces != NULL);
g_return_if_fail (trace_id < graph->num_traces);

t = graph->traces;
for (i = 0 ; i < trace_id ; i++)
	t = t->next;

if (t->Xdata != NULL)
    g_free(t->Xdata);
if (t->Ydata != NULL)
    g_free(t->Ydata);
    
t->Xdata = (gfloat *) malloc (n * sizeof(gfloat));
t->Ydata = (gfloat *) malloc (n * sizeof(gfloat));

for (i = 0 ; i < n ; i++)
    {
    t->Xdata[i] = xd[i];
    t->Ydata[i] = yd[i];

    if (yd[i] > ymax_store)
      ymax_store = yd[i];

    if (yd[i] < ymin_store)
      ymin_store = yd[i];

    if (xd[i] > xmax_store)
      xmax_store = xd[i];

    if (xd[i] < xmin_store )
      xmin_store = xd[i];

    }
t->ymax = ymax_store;
t->ymin = ymin_store;
t->xmax = xmax_store;
t->xmin = xmin_store;
t->num_points = n;

}
void gtk_graph_trace_format_marker(GtkGraph *graph, gint trace_id, GtkGraphMarkerType type, gint marker_size, GdkColor *fg, GdkColor *bg, gboolean is_filled)
{
GtkGraphTrace *t = graph->traces;;
GtkWidget *w = GTK_WIDGET(graph);
GtkGraphTraceFormat *current;

GdkColor zero;
GdkColor one;
GdkPixmap *target;
GdkGC *target_gc;
GdkPoint corner[4];

gint i, size = marker_size * 2 + 1;
	
g_return_if_fail (graph != NULL);
g_return_if_fail (GTK_IS_GRAPH (graph));
g_return_if_fail (graph->traces != NULL);
g_return_if_fail (trace_id < graph->num_traces);

for (i = 0 ; i < trace_id ; i++) /* Find the required trace */
	t = t->next;

current = t->format;
current->marker_type = type;
current->marker_filled = is_filled;
current->marker_size = marker_size;

/* Free up existing marker and allocate storage for a new one */
if (current->marker != NULL)
	gdk_pixmap_unref(current->marker);
current->marker = gdk_pixmap_new(w->window, size, size, -1);

if (current->mask != NULL)
	gdk_bitmap_unref(current->mask);
current->mask = gdk_pixmap_new(w->window, size, size, 1);

/* Free up existing Graphics Context for the marker and allocate a new one */
if (current->marker_gc != NULL)
	gdk_gc_unref(current->marker_gc);
current->marker_gc = gdk_gc_new(w->window);

if (current->mask_gc != NULL)
	gdk_gc_unref(current->mask_gc);
current->mask_gc = gdk_gc_new(current->mask);

gdk_color_black (gdk_colormap_get_system (), &zero);

if (zero.pixel != 0)
	{
    gdk_color_white (gdk_colormap_get_system(), &zero);
    gdk_color_black (gdk_colormap_get_system(), &one);
	}
else
	gdk_color_white (gdk_colormap_get_system(), &one);

/* Clear the mask then reset the correct bg color*/

gdk_gc_set_background (current->mask_gc, &zero);
gdk_gc_set_foreground (current->mask_gc, &zero);
gdk_draw_rectangle(current->mask, current->mask_gc, TRUE, 0, 0, size, size);    
gdk_gc_set_foreground (current->mask_gc, &one);

/* Now clear the marker and then reset the correct fg color*/
gdk_gc_set_background (current->marker_gc, bg);
gdk_gc_set_foreground (current->marker_gc, &WHITE);
gdk_draw_rectangle(current->marker, current->marker_gc, TRUE, 0, 0, size, size);    
gdk_gc_set_foreground (current->marker_gc, fg);

/* Now draw the required markers */
/* We run through the loop twice, once to draw them marker and once for the mask */

for (i = 0 ; i <=1  ; i++)
    {
    if (i == 0)
		{
        target = current->marker;
		target_gc = current->marker_gc;
		}
    else
		{
        target = current->mask;
		target_gc = current->mask_gc;
		}

	switch (current->marker_type)
	  {
	  case GTK_GRAPH_MARKER_NONE:
	    break;
	    
	  case GTK_GRAPH_MARKER_SQUARE:
	    if (is_filled)
	      if (i == 0)
		{
		  gdk_gc_set_foreground(target_gc, bg);
		  gdk_draw_rectangle(target, target_gc, TRUE, 0, 0, size, size);
		  gdk_gc_set_foreground(target_gc, fg);
		  gdk_draw_rectangle(target, target_gc, FALSE, 0, 0, size-1 , size-1 );
		}
	      else {
		gdk_draw_rectangle(target, target_gc, TRUE, 0, 0, size, size);
	      }
	    else
	      gdk_draw_rectangle(target, target_gc, FALSE, 0, 0, size - 1, size - 1); 
	    break;

	  case GTK_GRAPH_MARKER_CIRCLE:
			if (is_filled)
				if (i == 0)
					{
					gdk_gc_set_foreground(target_gc, bg);
					gdk_draw_arc(target, target_gc, TRUE, 0, 0, size, size, 0, 23040);
					gdk_gc_set_foreground(target_gc, fg);
				    gdk_draw_arc(target, target_gc, FALSE, 0, 0, size-1, size-1, 0, 23040);
					}
				else
				   gdk_draw_rectangle(target, target_gc, TRUE, 0, 0, size, size);
			else
			   gdk_draw_arc(target, target_gc, FALSE, 0, 0, size-1, size-1, 0, 23040);
			break;

	  case GTK_GRAPH_MARKER_DIAMOND:
			corner[0].x = marker_size; corner[0].y = 0;
			corner[1].x = size - 1; corner[1].y = marker_size;
			corner[2].x = marker_size; corner[2].y = size - 1;
			corner[3].x = 0; corner[3].y = marker_size;
			if (is_filled)
				if (i == 0)
					{
					gdk_gc_set_foreground(target_gc, bg);
					gdk_draw_polygon(target, target_gc, TRUE, corner, 4);
					gdk_gc_set_foreground(target_gc, fg);
					corner[1].x = size - 2;
					gdk_draw_polygon(target, target_gc, FALSE, corner, 4);
					corner[1].x = size - 1;
					}
				else
					gdk_draw_polygon(target, target_gc, TRUE, corner, 4);
			else
				gdk_draw_polygon(target, target_gc, FALSE, corner, 4);
			break;
    	case GTK_GRAPH_MARKER_TRIANGLE:
			corner[0].x = marker_size; corner[0].y = 0;
			corner[1].x = 1; corner[1].y = size - 3;
			corner[2].x = size - 2; corner[2].y = size - 3;
			gdk_draw_polygon(target, target_gc, is_filled, corner, 3);
	   	    break;
    	case GTK_GRAPH_MARKER_PLUS:
			if (is_filled)
				if (i == 0)
					{
					gdk_gc_set_foreground(target_gc, bg);
					gdk_draw_rectangle(target, target_gc, TRUE, 0, 0, size, size);
					gdk_gc_set_foreground(target_gc, fg);
					gdk_draw_line(target, target_gc, marker_size, 0, marker_size, size - 1);
					gdk_draw_line(target, target_gc, 0, marker_size, size - 1, marker_size);
					}
				else
					gdk_draw_rectangle(target, target_gc, TRUE, 0, 0, size, size);
			else
				{
				gdk_draw_line(target, target_gc, marker_size, 0, marker_size, size - 1);
				gdk_draw_line(target, target_gc, 0, marker_size, size - 1, marker_size);
				}
			break;
    	case GTK_GRAPH_MARKER_CROSS:
   			if (is_filled)
				if (i == 0)
					{
					gdk_gc_set_foreground(target_gc, bg);
					gdk_draw_rectangle(target, target_gc, TRUE, 0, 0, size, size);
					gdk_gc_set_foreground(target_gc, fg);
					gdk_draw_line(target, target_gc, 0, 0, size - 1, size - 1);
					gdk_draw_line(target, target_gc, size - 1, 0, 0, size - 1);
					}
				else
					gdk_draw_rectangle(target, target_gc, TRUE, 0, 0, size, size);
			else
				{
				gdk_draw_line(target, target_gc, 0, 0, size - 1, size - 1);
				gdk_draw_line(target, target_gc, size - 1, 0, 0, size - 1);
				}
			break;
	   case GTK_GRAPH_MARKER_STAR :
   			if (is_filled)
				if (i == 0)
					{
					gdk_gc_set_foreground(target_gc, bg);
					gdk_draw_rectangle(target, target_gc, TRUE, 0, 0, size, size);
					gdk_gc_set_foreground(target_gc, fg);
					gdk_draw_line(target, target_gc, 0, 0, size - 1, size - 1);
					gdk_draw_line(target, target_gc, size - 1, 0, 0, size - 1);
					gdk_draw_line(target, target_gc, marker_size, 0, marker_size, size - 1);
					gdk_draw_line(target, target_gc, 0, marker_size, size - 1, marker_size);
					}
				else
					gdk_draw_rectangle(target, target_gc, TRUE, 0, 0, size, size);
			else
				{
				gdk_draw_line(target, target_gc, 0, 0, size - 1, size - 1);
				gdk_draw_line(target, target_gc, size - 1, 0, 0, size - 1);
				gdk_draw_line(target, target_gc, marker_size, 0, marker_size, size - 1);
				gdk_draw_line(target, target_gc, 0, marker_size, size - 1, marker_size);
				}
    		break;
	   }
    }
gdk_gc_set_clip_mask(current->marker_gc, current->mask);
}

/*! \fn gtk_graph_trace_format_line
 */

void gtk_graph_trace_format_line(GtkGraph *graph, gint trace_id, GtkGraphLineType type, gint width, GdkColor *line_color, gboolean visible )
{
GtkGraphTrace *t;
gint i;
gint8 dash_list[4];
g_return_if_fail (graph != NULL);
g_return_if_fail (GTK_IS_GRAPH (graph));
g_return_if_fail (graph->traces != NULL);
g_return_if_fail (trace_id < graph->num_traces);

t = graph->traces;
for (i = 0 ; i < trace_id ; i++)
	t = t->next;

t->format->line_visible = visible;

gdk_gc_set_foreground(t->format->line_gc, line_color);
gdk_gc_set_background(t->format->line_gc, &WHITE);
switch (type)
	{
	case NO_LINE:
		gdk_gc_set_line_attributes(t->format->line_gc, width, GDK_LINE_SOLID, GDK_CAP_BUTT, GDK_JOIN_ROUND);
		t->format->line_visible = FALSE;
		break;
	case SOLID:
		gdk_gc_set_line_attributes(t->format->line_gc, width, GDK_LINE_SOLID, GDK_CAP_BUTT, GDK_JOIN_BEVEL);
		t->format->line_visible = TRUE;
		break;
	case DASHED:
		gdk_gc_set_line_attributes(t->format->line_gc, width, GDK_LINE_ON_OFF_DASH, GDK_CAP_BUTT, GDK_JOIN_ROUND);
		dash_list[0] = 6;
		dash_list[1] = 6;
		gdk_gc_set_dashes(t->format->line_gc, 0, dash_list, 2);
		t->format->line_visible = TRUE;
		break;
	case DASH_DOT:
		gdk_gc_set_line_attributes(t->format->line_gc, width, GDK_LINE_ON_OFF_DASH, GDK_CAP_BUTT, GDK_JOIN_ROUND);
		dash_list[0] = 3;
		dash_list[1] = 2;
		dash_list[2] = 1;
		dash_list[3] = 2;
		gdk_gc_set_dashes(t->format->line_gc, 0, dash_list, 4);
		t->format->line_visible = TRUE;
		break;
	case DOTTED:
		gdk_gc_set_line_attributes(t->format->line_gc, width, GDK_LINE_ON_OFF_DASH, GDK_CAP_BUTT, GDK_JOIN_ROUND);
		dash_list[0] = 2;
		dash_list[1] = 3;
		gdk_gc_set_dashes(t->format->line_gc, 0, dash_list, 2);
		t->format->line_visible = TRUE;
		break;
	case LONG_DASH:
		gdk_gc_set_line_attributes(t->format->line_gc, width, GDK_LINE_ON_OFF_DASH, GDK_CAP_BUTT, GDK_JOIN_ROUND);
		dash_list[0] = 10;
		dash_list[1] = 6;
		gdk_gc_set_dashes(t->format->line_gc, 0, dash_list, 2);
		t->format->line_visible = TRUE;
		break;
	}

}
void gtk_graph_trace_format_title(GtkGraph *graph, gint trace_id, const gchar *legend_text)
{
GtkGraphTrace *t;
int i;
g_return_if_fail (graph != NULL);
g_return_if_fail (GTK_IS_GRAPH (graph));
g_return_if_fail (graph->traces != NULL);
g_return_if_fail (trace_id < graph->num_traces);

t = graph->traces;
for (i = 0 ; i < trace_id ; i++)
	t = t->next;

if (t->format->legend_text != NULL)
	g_free (t->format->legend_text);
t->format->legend_text = g_strdup(legend_text);
}
