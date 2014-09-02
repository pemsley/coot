#include <math.h>
#include <stdio.h>

#include "gtkgraph.h"
#include "gtkgraph_internal.h"

static gint centre_x;
static gint centre_y;
static gint radius;

#define TWOPI	(2.0 * M_PI)

/* Polar.c: Routines for plotting Polar plots */
void gtk_graph_polar_plot_axes(GtkGraph *graph)
{
 gint i, j;

 gint xpos = 0, ypos = 0;		// X and Y co-ordinate locations of radial labels
 gint x, y, w, h;
 gint cleft, cright, ctop, cbot, chalfleft, chalfright, chalftop, chalfbot;
 gfloat angle, new_angle, r;

 gchar label[10], format[10];
 gint text_width, text_height;
 gint max_store, min_tick_value;
 PangoFontDescription *fontdesc = NULL;
 PangoLayout *layout = NULL;

g_return_if_fail (graph != NULL);
g_return_if_fail (GTK_IS_GRAPH (graph));
g_return_if_fail (GTK_WIDGET_REALIZED (graph));
g_return_if_fail (graph->graph_type == POLAR);
	
 fontdesc = pango_font_description_from_string("Sans 10");
 layout = gtk_widget_create_pango_layout(GTK_WIDGET(graph), NULL);
 pango_layout_set_font_description(layout, fontdesc);
  
 /* Plot the concentric grid */

 centre_x = user_origin_x + user_width / 2;
 centre_y = user_origin_y + user_height / 2;
 radius = graph->dependant->n_maj_tick * graph->dependant->pxls_per_maj_tick;

 for (i = 0 ; i <= graph->dependant->n_maj_tick ; i++)
	{
	x = centre_x - i * graph->dependant->pxls_per_maj_tick;
	y = centre_y - i * graph->dependant->pxls_per_maj_tick;
	w = 2.0 * i * graph->dependant->pxls_per_maj_tick;
	h = 2.0 * i * graph->dependant->pxls_per_maj_tick;
	gdk_draw_arc(buffer, Gridcontext, FALSE, x, y, w, h, 0, 23040);
	}
/* Determine the position of the radial lines and plot them */

 cleft = centre_x - radius;
 cright = centre_x + radius;
 ctop = centre_y - radius;
 cbot = centre_y + radius;
 r = (gfloat)(radius) / sqrt(2.0);
 chalfleft = centre_x - (gint)(r);
 chalfright = centre_x + (gint)(r);
 chalftop = centre_y - (gint)(r);
 chalfbot = centre_y + (gint)(r);
	
 gdk_draw_line(buffer, Gridcontext, cleft, centre_y, cright, centre_y);
 gdk_draw_line(buffer, Gridcontext, chalfleft, chalftop, chalfright, chalfbot);
 gdk_draw_line(buffer, Gridcontext, chalfright, chalftop, chalfleft, chalfbot);
 gdk_draw_line(buffer, BandWcontext, centre_x, cbot, centre_x, centre_y);
 gdk_draw_line(buffer, Gridcontext, centre_x, ctop, centre_x, centre_y);

/* Finally label the radial lines */

switch(graph->polar_format.type)
    {
    case RADIANS_SYMMETRIC:
    case RADIANS_ANTISYMMETRIC:
        for (i = 0 ; i < 8 ; i++)
            {
            angle = M_PI / 4 * i - graph->polar_format.polar_start;
            new_angle = wrap_angle(angle, graph->polar_format);
            g_snprintf(label, 10, "%.3f", new_angle);
        	pango_layout_set_text(layout, label, -1);
        	pango_layout_get_pixel_size(layout, &text_width, &text_height);
			switch(i)
				{
				case 0:
					xpos = centre_x - text_width / 2;
					ypos = ctop - text_height - 2;
					break;
				case 1:
					xpos = chalfright + 2;
					ypos = chalftop - text_height - 2;
					break;
				case 2:
					xpos = cright + 4;
					ypos = centre_y - text_height / 2;
					break;
				case 3:
					xpos = chalfright + 2;
					ypos = chalfbot + 2;
					break;				
				case 4:
					xpos = centre_x - text_width / 2;
					ypos = cbot + 2;
					break;
				case 5:
					xpos = chalfleft - text_width - 2;
					ypos = chalfbot + 2;
					break;
				case 6:
					xpos = cleft - text_width - 2;
					ypos = centre_y - text_height / 2;
					break;
				case 7:
					xpos = chalfleft - text_width - 2;
					ypos = chalftop - text_height - 2;
					break;
				}
           	gdk_draw_layout(buffer, BandWcontext, xpos, ypos, layout);
			}
        break;
    case DEGREES_SYMMETRIC:
    case DEGREES_ANTISYMMETRIC:
        for (i = 0 ; i < 8 ; i++)
            {
            angle = 45.0 * i - graph->polar_format.polar_start;
            new_angle = wrap_angle(angle, graph->polar_format);
            g_snprintf(label, 10, "%.0f", new_angle);
            pango_layout_set_text(layout, label, -1);
           	pango_layout_get_pixel_size(layout, &text_width, &text_height);
			switch(i)
				{
				case 0:
					xpos = centre_x - text_width / 2;
					ypos = ctop - text_height - 2;
					break;
				case 1:
					xpos = chalfright + 2;
					ypos = chalftop - text_height - 2;
					break;
				case 2:
					xpos = cright + 4;
					ypos = centre_y - text_height / 2;
					break;
				case 3:
					xpos = chalfright + 2;
					ypos = chalfbot + 2;
					break;				
				case 4:
					xpos = centre_x - text_width / 2;
					ypos = cbot + 2;
					break;
				case 5:
					xpos = chalfleft - text_width - 2;
					ypos = chalfbot + 2;
					break;
				case 6:
					xpos = cleft - text_width - 2;
					ypos = centre_y - text_height / 2;
					break;
				case 7:
					xpos = chalfleft - text_width - 2;
					ypos = chalftop - text_height - 2;
					break;
				}
           	gdk_draw_layout(buffer, BandWcontext, xpos, ypos, layout);
			}
        break;
    }

/* Determine the Scaling for and Plot the Y-Axis */

 max_store = (int) rint(-floor(log10(graph->dependant->maj_tick)));	/* simplest axis labels */

 if (fabs(max_store) >= 5)
	g_snprintf(format, 10, "%%.1E");
 else if (max_store <= 0)
	g_snprintf(format, 10, "%%.0f");
 else
	g_snprintf(format, 10, "%%.%df", max_store);

 for (i = 0; i <= graph->dependant->n_maj_tick ; i++)
	{
	gdk_draw_line(buffer, BandWcontext, centre_x, centre_y + i * graph->dependant->pxls_per_maj_tick, centre_x - 6, centre_y + i * graph->dependant->pxls_per_maj_tick);
	g_snprintf(label, 10, format, graph->dependant->axis_min + i * graph->dependant->maj_tick);/* simplest axis labels */
	pango_layout_set_text(layout, label, -1);
	pango_layout_get_pixel_size(layout, &text_width, &text_height); 
	gdk_draw_layout(buffer, BandWcontext, centre_x - 6 - text_width, centre_y + i * graph->dependant->pxls_per_maj_tick - text_height/2, layout);
	for (j = 1 ; j < graph->dependant->n_min_tick ; j++)
		if (i < graph->dependant->n_maj_tick)
			{
			min_tick_value = (gfloat) j * (gfloat) graph->dependant->pxls_per_maj_tick / (gfloat) (graph->dependant->n_min_tick);
			y = centre_y + i * graph->dependant->pxls_per_maj_tick + (gint) min_tick_value;
			gdk_draw_line(buffer, BandWcontext, centre_x , y, centre_x - 3, y);
			}
	}

if (graph->dependant->title != NULL)
	{
	pango_layout_set_text(layout, graph->dependant->title, -1);
	pango_layout_get_pixel_size(layout, &text_width, &text_height);
	gdk_draw_layout(buffer, BandWcontext, centre_x + margin/2, cbot - text_height, layout);
	}
	
 pango_font_description_free(fontdesc);
 g_object_unref(layout);	
}



void gtk_graph_polar_plot_traces(GtkGraph *graph)
{
gint i, n;
GtkGraphTrace *tmp;
gfloat theta, new_angle, r;
GdkPoint *pts = NULL;

g_return_if_fail (graph != NULL);
g_return_if_fail (GTK_IS_GRAPH (graph));
g_return_if_fail (GTK_WIDGET_REALIZED (graph));
g_return_if_fail (graph->graph_type == POLAR);
	
tmp = graph->traces;
for (n = 0 ; n < graph->num_traces ; n++)
	{
	/* Make sure that there is some data in the trace */
	if (tmp->Xdata == NULL || tmp->Ydata == NULL) 
		continue;

	/* Assign the storage for the co-ordinates of each data point */

	pts = (GdkPoint *) g_malloc ((tmp->num_points) * sizeof(GdkPoint));

// Calculate the co-ordinates of the data points 

	new_angle = wrap_angle(tmp->Xdata[0], graph->polar_format);
	if (graph->polar_format.type >> 2) // i.e. one of the radian options 
		theta = new_angle + graph->polar_format.polar_start;
	else							// otherwise we must be in degrees 
      	theta = (new_angle + graph->polar_format.polar_start) * M_PI / 180.0;

/* Deal with the Goldilocks problem */
	if (tmp->Ydata[0] < graph->dependant->axis_min) /*this one's too small */
		r = 0;
	else if (tmp->Ydata[0] > graph->dependant->axis_max) /* this one's too large */
		r = (graph->dependant->axis_max - graph->dependant->axis_min) * graph->dependant->scale_factor;
	else /* this one's just right */
		r = (tmp->Ydata[0] - graph->dependant->axis_min) * graph->dependant->scale_factor; 
	
    pts[0].x = centre_x + (gint) rint(r * sin(theta));
	pts[0].y = centre_y - (gint) rint(r * cos(theta));

 	for (i = 1 ; i < tmp->num_points ; i++)
      	{
		new_angle = wrap_angle(tmp->Xdata[i], graph->polar_format);
		if (graph->polar_format.type >> 2) // i.e. one of the radian options 
			theta = new_angle + graph->polar_format.polar_start;
		else							// otherwise we must be in degrees 
			theta = (new_angle + graph->polar_format.polar_start) * M_PI / 180.0;

/* Deal with the Goldilocks problem */
		if (tmp->Ydata[i] < graph->dependant->axis_min) /*this one's too small */
			r = 0;
		else if (tmp->Ydata[i] > graph->dependant->axis_max) /* this one's too large */
			r = (graph->dependant->axis_max - graph->dependant->axis_min) * graph->dependant->scale_factor;
		else /* this one's just right */
			r = (tmp->Ydata[i] - graph->dependant->axis_min) * graph->dependant->scale_factor;
	
		pts[i].x = centre_x + (int) rint(r * sin(theta));
		pts[i].y = centre_y - (int) rint(r * cos(theta));
		}

	gdk_draw_lines (buffer, tmp->format->line_gc, pts, tmp->num_points);/* Draw the lines */	

	for (i = 1 ; i < tmp->num_points ; i++)/* and then draw the markers */	
        if (tmp->format->marker_type != GTK_GRAPH_MARKER_NONE)
			{
			gdk_gc_set_clip_origin(tmp->format->marker_gc, pts[i].x - tmp->format->marker_size, pts[i].y - tmp->format->marker_size);
            gdk_draw_pixmap(buffer, tmp->format->marker_gc, tmp->format->marker, 0, 0, pts[i].x - tmp->format->marker_size, pts[i].y - tmp->format->marker_size, -1, -1);
			}

	g_free(pts);// Free up all storage after use 
	tmp = tmp->next;// and then move onto the next trace 
	
	}

}


void gtk_graph_set_polar_format(GtkGraph *graph, GtkGraphPolarType type, gfloat polar_start)
{
g_return_if_fail (graph != NULL);
g_return_if_fail (GTK_IS_GRAPH (graph));
g_return_if_fail (graph->graph_type == POLAR);

graph->polar_format.type = type;
graph->polar_format.polar_start = polar_start;
}
