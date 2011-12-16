#include <gtk/gtk.h>
#include <math.h>
#include <complex.h>
#include "gtkgraph.h"
#include "gtkgraph_internal.h"

static gint centre_x;
static gint centre_y;

/* Important to note that the Smith Chart actually plots reflection coefficient and thus */
/* should range with |R| between the limits of +1 and -1. Although the grid lines appear */
/* polar the chart is actually plotted in an X/Y co-ordinate system                      */

/* Smith.c: Routines for plotting Smith charts */
void gtk_graph_smith_plot_axes(GtkGraph *graph)
{
gint i ;
gint x, y, w, h;
gfloat ri[] = {0.2, 0.5, 1.0, 2.0, 5.0}, rxl[]={0.20, 0.5, 1, 2, 5}, coeff1, r, store, arc_angle;
gfloat yend, xend; // The x and y co-ordinates of the end of the reactance lines
gchar label[10];
gint text_width, text_height;
gint xpos = 0, ypos = 0;		// X and Y co-ordinate locations of radial labels

PangoFontDescription *fontdesc = NULL;
PangoLayout *layout = NULL;

g_return_if_fail (graph != NULL);
g_return_if_fail (GTK_IS_GRAPH (graph));
g_return_if_fail (GTK_WIDGET_REALIZED (graph));
g_return_if_fail (graph->graph_type == SMITH);

fontdesc = pango_font_description_from_string("Sans 10");
layout = gtk_widget_create_pango_layout(GTK_WIDGET(graph), NULL);
pango_layout_set_font_description(layout, fontdesc);

/* Plot the concentric grid */

centre_x = user_origin_x + user_width / 2;
centre_y = user_origin_y + user_width / 2;
gdk_draw_arc(buffer, Gridcontext, FALSE, centre_x - graph->dependant->scale_factor, centre_y - graph->dependant->scale_factor, 2.0 * graph->dependant->scale_factor, 2.0 * graph->dependant->scale_factor, 0, 23040);

for (i = 0 ; i < 5 ; i++)
	{
/* Draw the circles of Constant Resistance */
	coeff1 = ri[i] / (1.0+ri[i]);
	r = 1.0 / (ri[i] + 1.0);
	x = centre_x + (gint)((coeff1 - r) * graph->dependant->scale_factor);
	y = centre_y - (gint) (r * graph->dependant->scale_factor);
	w = 2.0 * r * graph->dependant->scale_factor;
	h = 2.0 * r * graph->dependant->scale_factor;
	gdk_draw_arc(buffer, Gridcontext, FALSE, x, y, w, h, 0, 23039);
		
	g_snprintf(label, 10, "%.0f", ri[i] * graph->smith_Z0);
    pango_layout_set_text(layout, label, -1);
    pango_layout_get_pixel_size(layout, &text_width, &text_height);
	xpos = x - text_width / 2;
	ypos = centre_y ;
	gdk_draw_layout(buffer, BandWcontext, xpos, ypos, layout);

/* Draw the circles of Constant Reactance */

	r = 1.0 / rxl[i];
	x = centre_x + (gint) ((1 - r) * graph->dependant->scale_factor);
	w = 2.0 * r * graph->dependant->scale_factor;
	h = 2.0 * r * graph->dependant->scale_factor;

/* Now Label the Constant Reactance curves */
	g_snprintf(label, 10, "+j%0.f", rxl[i]*graph->smith_Z0);
    pango_layout_set_text(layout, label, -1);
    pango_layout_get_pixel_size(layout, &text_width, &text_height);
	store = 2.0 * r/(r*r + 1);
	yend = store * (gfloat)(graph->dependant->scale_factor);
	xend = sqrt(1 - store*store) * (gfloat)(graph->dependant->scale_factor);
	
	switch(i)
		{
		case 0:
		case 1:
			g_snprintf(label, 10, "+j%0.f", rxl[i]*graph->smith_Z0);
			pango_layout_set_text(layout, label, -1);
			pango_layout_get_pixel_size(layout, &text_width, &text_height);
			ypos = centre_y - (gint)yend;
			xpos = centre_x - (gint)xend;
			gdk_draw_layout(buffer, BandWcontext, xpos - text_width, ypos - text_height, layout);
			g_snprintf(label, 10, "-j%0.f", rxl[i]*graph->smith_Z0);
			pango_layout_set_text(layout, label, -1);
			ypos = centre_y + (gint)yend;
			gdk_draw_layout(buffer, BandWcontext, xpos - text_width, ypos, layout);
			arc_angle = atan((gfloat)(centre_x + graph->dependant->scale_factor - xpos)/(gfloat)((centre_y + r * graph->dependant->scale_factor)-ypos)) *180.0 / 3.141592654;
			gdk_draw_arc(buffer, Gridcontext, FALSE, x, centre_y, w, h, 5760, (gint)(arc_angle * 64.0)); // Capacitive Circles
			gdk_draw_arc(buffer, Gridcontext, FALSE, x, centre_y - (gint) ((2 *r) * graph->dependant->scale_factor), w, h, 17280 - (gint) (arc_angle*64.0), (gint) (arc_angle*64.0));	// Inductive Circles
			break;
		case 2:
			g_snprintf(label, 10, "+j%0.f", rxl[i]*graph->smith_Z0);
			pango_layout_set_text(layout, label, -1);
			pango_layout_get_pixel_size(layout, &text_width, &text_height);
			ypos = centre_y - (gint)yend;
			xpos = centre_x - (gint)xend;
			gdk_draw_layout(buffer, BandWcontext, xpos - text_width/2, ypos - text_height, layout);
			g_snprintf(label, 10, "-j%0.f", rxl[i]*graph->smith_Z0);
			pango_layout_set_text(layout, label, -1);
			ypos = centre_y + yend;
			gdk_draw_layout(buffer, BandWcontext, xpos - text_width/2, ypos, layout);
			arc_angle = atan((gfloat)(centre_x + graph->dependant->scale_factor - xpos)/(gfloat)(ypos - (centre_y + r * graph->dependant->scale_factor))) *180.0 / 3.141592654;
			gdk_draw_arc(buffer, Gridcontext, FALSE, x, centre_y, w, h, 5760, (gint) (arc_angle*64.0)); // Capacitive Circles
			gdk_draw_arc(buffer, Gridcontext, FALSE, x, centre_y - (gint) ((2 *r) * graph->dependant->scale_factor), w, h, 17280 - (gint) (arc_angle*64.0), (gint) (arc_angle*64.0));	// Inductive Circles
			break;
		case 3:
		case 4:
			g_snprintf(label, 10, "+j%0.f", rxl[i]*graph->smith_Z0);
			pango_layout_set_text(layout, label, -1);
			pango_layout_get_pixel_size(layout, &text_width, &text_height);
			ypos = centre_y - (gint)yend;
			xpos = centre_x + (gint)xend;
			gdk_draw_layout(buffer, BandWcontext, xpos, ypos - text_height, layout);
			g_snprintf(label, 10, "-j%0.f", rxl[i]*graph->smith_Z0);
			pango_layout_set_text(layout, label, -1);
			ypos = centre_y + (gint)yend;
			gdk_draw_layout(buffer, BandWcontext, xpos, ypos, layout);
			arc_angle = atan((gfloat)(centre_x + graph->dependant->scale_factor - xpos)/(gfloat)(ypos - (centre_y + r * graph->dependant->scale_factor))) *180.0 / 3.141592654;
			gdk_draw_arc(buffer, Gridcontext, FALSE, x, centre_y, w, h, 5760, 11520 - (gint) (arc_angle*64.0)); // Capacitive Circles
			gdk_draw_arc(buffer, Gridcontext, FALSE, x, centre_y - (gint) ((2 *r) * graph->dependant->scale_factor), w, h, 5760 + (gint) (arc_angle*64.0), 11520 - (gint) (arc_angle*64.0));	// Inductive Circles
			break;
		}
	}
	
/* Draw the Real Axis */	
gdk_draw_line(buffer, Gridcontext, centre_x - user_width/2, centre_y, centre_x + user_width/2, centre_y);	
g_snprintf(label, 10, "+inf");
pango_layout_set_text(layout, label, -1);
pango_layout_get_pixel_size(layout, &text_width, &text_height);
xpos = centre_x + user_width/2.0 - text_width/2 ;
ypos = centre_y ;
gdk_draw_layout(buffer, BandWcontext, xpos, ypos, layout);
	
pango_font_description_free(fontdesc);
g_object_unref(layout);	
}


void gtk_graph_smith_plot_traces(GtkGraph *graph)
{
gint i, n;
GtkGraphTrace *tmp;
gfloat CA, CB, CC, CD;

GdkPoint *pts = NULL;

g_return_if_fail (graph != NULL);
g_return_if_fail (GTK_IS_GRAPH (graph));
g_return_if_fail (GTK_WIDGET_REALIZED (graph));
g_return_if_fail (graph->graph_type == SMITH);
	
tmp = graph->traces;
for (n = 0 ; n < graph->num_traces ; n++)
	{
	/* Make sure that there is some data in the trace */
	if (tmp->Xdata == NULL || tmp->Ydata == NULL) 
		continue;

	/* Assign the storage for the co-ordinates of each data point */

	pts = (GdkPoint *) g_malloc ((tmp->num_points) * sizeof(GdkPoint));

	/* The Xdata array contains the resistance whilst the Ydata array has the reactance */
	
	CA = tmp->Xdata[0] - graph->smith_Z0;
	CB = tmp->Ydata[0];
	CC = tmp->Xdata[0] + graph->smith_Z0;
	CD = tmp->Ydata[0];
		
	pts[0].x = centre_x + (CA*CC + CB*CD)/(CC*CC + CD*CD) * graph->dependant->scale_factor;
	pts[0].y = centre_y - (CB*CC - CA*CD)/(CC*CC + CD*CD) * graph->dependant->scale_factor;
	
 	for (i = 1 ; i < tmp->num_points ; i++)
      	{
		CA = tmp->Xdata[i] - graph->smith_Z0;
		CB = tmp->Ydata[i];
		CC = tmp->Xdata[i] + graph->smith_Z0;
		CD = tmp->Ydata[i];
		
		pts[i].x = centre_x + (CA*CC + CB*CD)/(CC*CC + CD*CD) * graph->dependant->scale_factor;
		pts[i].y = centre_y - (CB*CC - CA*CD)/(CC*CC + CD*CD) * graph->dependant->scale_factor;
		}

	gdk_draw_lines (buffer, tmp->format->line_gc, pts, tmp->num_points);/* Draw the lines */	
	for (i = 1 ; i < tmp->num_points ; i++)/* and then draw the markers */	
        if (tmp->format->marker_type != GTK_GRAPH_MARKER_NONE)
			{
			gdk_gc_set_clip_origin(tmp->format->marker_gc, pts[i].x-5, pts[i].y-5);
            gdk_draw_pixmap(buffer, tmp->format->marker_gc, tmp->format->marker, 0, 0, pts[i].x-5, pts[i].y-5, -1, -1);
			}

	g_free(pts);// Free up all storage after use 
	tmp = tmp->next;// and then move onto the next trace 
	
	}
	
}

/**
 * gtk_graph_smith_set_Z0:
 * @graph:  the #GtkGraph containing the Smith Chart to be modified
 * @Z0: the unique identifier of the trace
 *
 * Sets the normalising value for Smith's Charts.  All Smith's Charts
 * are normalised to a particular value - A #GtkGraph defaults to a
 * normalisation of 50 Ohms.  This function allows the user to select
 * any @Z0 that they choose
 * 
 */
void gtk_graph_smith_set_Z0(GtkGraph *graph, gfloat Z0)
{
g_return_if_fail (graph != NULL);
g_return_if_fail (GTK_IS_GRAPH (graph));
g_return_if_fail (GTK_WIDGET_REALIZED (graph));

graph->smith_Z0 = Z0;
}
