#include <stdio.h>
#include <math.h>
#include "libgtkgraph/gtkgraph.h"
#include "libgtkgraph/gtkgraph_internal.h"
/* Declaration of all local functions */

static GtkGraphAnnotation *gtk_graph_annotation_allocate(void);

/* Externally referenceable functions */

static GtkGraphAnnotation *gtk_graph_annotation_allocate(void)
{
GtkGraphAnnotation *tmp;

tmp = (GtkGraphAnnotation *) g_malloc(sizeof(GtkGraphAnnotation));

if (tmp == (GtkGraphAnnotation *) NULL)
    {
        (void) fprintf(stderr,"malloc failed at: %s\n","gtk_graph_trace_allocate");
        return ((GtkGraphAnnotation *) NULL);
    }
tmp->type = VERTICAL;
tmp->value = 0;
tmp->text = NULL;
tmp->next = NULL;

return(tmp);
}

/**
 * gtk_graph_annotation_new:
 * @graph:  The #GtkGraph that you wish to add the annotation to
 *	
 * Adds an annotation to a pre-existing #GtkGraph @graph. 
 * 
 * Returns: an unique integer identifying the new annotation.  This must be
 * stored if you wish to alter the properties of the annotation
 */
gint gtk_graph_annotation_new(GtkGraph *graph)
{

int i;
GtkGraphAnnotation *new_annotation, *tmp;

g_return_val_if_fail (graph != NULL, FALSE);
g_return_val_if_fail (GTK_IS_GRAPH (graph), FALSE);

new_annotation = gtk_graph_annotation_allocate();

if (graph->annotations == NULL)
	graph->annotations = new_annotation;
else
	{
	for (i = 0, tmp = graph->annotations; tmp->next != NULL ; tmp = tmp->next , i++)  // Should advance us to last item in list
        ;
	tmp->next = new_annotation;
	}
graph->num_annotations += 1;
return graph->num_annotations - 1;
}

/**
 * gtk_graph_annotation_set_data:
 * @graph:  the #GtkGraph that contains the annotation you wish to modify
 * @annotation_id:  the unique integer ID returned by gtk_graph_annotation_new()
 * @type:	
 * @value:	value at which to display the annotation
 * @text:  any text that you want displayed next to the
 *
 * The function description goes here. You can use @par1 to refer to parameters
 * so that they are highlighted in the output. You can also use %constant
 * for constants, function_name2() for functions and #GtkWidget for links to
 * other declarations (which may be documented elsewhere).
 */
void gtk_graph_annotation_set_data(GtkGraph *graph, gint annotation_id, GtkGraphAnnotationType type, gfloat value, gchar *text)
{
GtkGraphAnnotation *t;
gint i;
g_return_if_fail (graph != NULL);
g_return_if_fail (GTK_IS_GRAPH (graph));
if (graph->annotations == NULL)
	return;
g_return_if_fail (annotation_id < graph->num_annotations);

t = graph->annotations;
for (i = 0 ; i < annotation_id ; i++)
	t = t->next;

t->type = type;
t->value = value;
t->text = g_strdup(text);

}
void gtk_graph_plot_annotations(GtkGraph *graph)
{
gint n;
gint plotting_value, text_width, text_height;
GtkGraphAnnotation *tmp;
PangoFontDescription *fontdesc = NULL;
PangoLayout *layout = NULL;
gint centre_x = user_width / 2.0 + user_origin_x, xpos;
gint centre_y = user_height / 2.0 + user_origin_y, ypos;	
gfloat theta, new_angle, r;
gchar *text_buffer = NULL;
//gint horizontal_position = user_origin_x;
gint vertical_position = user_origin_y;

g_return_if_fail (graph != NULL);
g_return_if_fail (GTK_IS_GRAPH (graph));
g_return_if_fail (GTK_WIDGET_REALIZED(GTK_WIDGET(graph)));

if (graph->annotations == NULL)
	return;

fontdesc = pango_font_description_from_string("Sans 9");
layout = gtk_widget_create_pango_layout(GTK_WIDGET(graph), NULL);
pango_layout_set_font_description(layout, fontdesc);		
	
tmp = graph->annotations;
for (n = 0 ; n < graph->num_annotations ; n++)
	{
	switch(tmp->type)
		{
		case HORIZONTAL:
			if (graph->graph_type != XY)	/* Skip on if the graph is not the correct type */
				continue;

/* Ensure that we are within the plotable range and then draw the line*/

			plotting_value = (gint) rint((graph->dependant->axis_max - tmp->value ) * graph->dependant->scale_factor) + user_origin_y;
			if (tmp->value < graph->dependant->axis_min)
				plotting_value = (gint) rint((graph->dependant->axis_max - graph->dependant->axis_min) * graph->dependant->scale_factor) + user_origin_y;
			if (tmp->value > graph->dependant->axis_max)
				plotting_value = user_origin_y;
			gdk_draw_line(buffer, BlueandWcontext, user_origin_x, plotting_value, (graph->independant->n_maj_tick * graph->independant->pxls_per_maj_tick ) + user_origin_x, plotting_value);

/* Set text and then plot above line if value is below half way and vice versa */
			
			if (text_buffer)
				g_free(text_buffer);
			if (tmp->text)
				text_buffer = g_strdup_printf( "%s\nY: %.2f", tmp->text, tmp->value);
			else
				text_buffer = g_strdup_printf("Y: %.2f", tmp->value);
			pango_layout_set_text(layout, text_buffer, -1);
			pango_layout_get_pixel_size(layout, &text_width, &text_height);

			if (tmp->value <= (graph->dependant->axis_max - graph->dependant->axis_min) / 2.0)
				gdk_draw_layout(buffer, BandWcontext, (graph->independant->n_maj_tick * graph->independant->pxls_per_maj_tick )/2 + user_origin_x, plotting_value -  text_height - 2, layout);
			else
				gdk_draw_layout(buffer, BandWcontext, (graph->independant->n_maj_tick * graph->independant->pxls_per_maj_tick )/2 + user_origin_x, plotting_value + 2, layout);
			break;
		case VERTICAL:
			if (graph->graph_type != XY)	/* Skip on if the graph is not the correct type */
				continue;

/* Ensure that we are within the plotable range and then draw the line*/
			
			plotting_value = (gint) rint((tmp->value - graph->independant->axis_min ) * graph->independant->scale_factor) + user_origin_x;
			if (tmp->value < graph->independant->axis_min)
				plotting_value = user_origin_x;
			if (tmp->value > graph->independant->axis_max)
				plotting_value = user_origin_x  + graph->independant->pxls_per_maj_tick * graph->independant->n_maj_tick;
			gdk_draw_line(buffer, BlueandWcontext, plotting_value, user_origin_y, plotting_value, user_origin_y + graph->dependant->pxls_per_maj_tick * graph->dependant->n_maj_tick);

/* Set text and then plot to the left if line on the right of the middle and vice versa */

			if (text_buffer)
				g_free(text_buffer);
			if (tmp->text)
				text_buffer = g_strdup_printf("%s\nX: %.2f", tmp->text, tmp->value);
			else
				text_buffer = g_strdup_printf("X: %.2f", tmp->value);
			pango_layout_set_text(layout, text_buffer, -1);
			pango_layout_get_pixel_size(layout, &text_width, &text_height);
			if (tmp->value <= (graph->independant->axis_max - graph->independant->axis_min) / 2.0)
				gdk_draw_layout(buffer, BandWcontext, plotting_value +2, vertical_position, layout);
			else
				gdk_draw_layout(buffer, BandWcontext, plotting_value - text_width - 2, vertical_position, layout);
			vertical_position += text_height;
			break;
		case RADIAL:
			if (graph->graph_type != POLAR)
				continue;
			new_angle = wrap_angle(tmp->value, graph->polar_format);
			if (graph->polar_format.type >> 2) // i.e. one of the radian options 
				theta = new_angle + graph->polar_format.polar_start;
			else							// otherwise we must be in degrees 
				theta = (new_angle + graph->polar_format.polar_start) * M_PI / 180.0;
			xpos = centre_x + (gint) (user_width / 2.0 * sin(theta));
			ypos = centre_y - (gint) (user_height / 2.0 * cos(theta));
			gdk_draw_line(buffer, BlueandWcontext, centre_x, centre_y, xpos, ypos);
			if (text_buffer)
				g_free(text_buffer);
			if (tmp->text)
				text_buffer = g_strdup_printf("%s\nTheta: %.2f", tmp->text, tmp->value);
			else
				text_buffer = g_strdup_printf("Theta: %.1f", tmp->value);
			pango_layout_set_text(layout, text_buffer, -1);
			pango_layout_get_pixel_size(layout, &text_width, &text_height);
			if (xpos > centre_x)
				{
				if (ypos >= centre_y)
					gdk_draw_layout(buffer, BandWcontext, xpos, ypos, layout);
				else
					gdk_draw_layout(buffer, BandWcontext, xpos, ypos-text_height - 2, layout);
				}
			else
				{
				if (ypos >= centre_y)
					gdk_draw_layout(buffer, BandWcontext, xpos - text_width - 2, ypos, layout);
				else
					gdk_draw_layout(buffer, BandWcontext, xpos - text_width - 2, ypos - text_height - 2, layout);
				}
			break;
		case AZIMUTHAL:
			if (graph->graph_type != POLAR)
				continue;
			plotting_value = (gint) rint(((tmp->value - graph->dependant->axis_min ) * graph->dependant->scale_factor));
			if (tmp->value < graph->dependant->axis_min)
				plotting_value = 0;
			if (tmp->value > graph->dependant->axis_max)
				plotting_value = (gint) rint((graph->dependant->axis_max - graph->dependant->axis_min) * graph->dependant->scale_factor);
			gdk_draw_arc(buffer, BlueandWcontext, FALSE, centre_x - plotting_value, centre_y - plotting_value, 2.0 * plotting_value, 2.0 * plotting_value, 0, 23040);
			if (text_buffer)
				g_free(text_buffer);
			if (tmp->text)
				text_buffer = g_strdup_printf( "%s\nRadius: %.1f", tmp->text, tmp->value);
			else
				text_buffer = g_strdup_printf( "Radius: %.1f", tmp->value);
			pango_layout_set_text(layout, text_buffer, -1);
			pango_layout_get_pixel_size(layout, &text_width, &text_height);
			gdk_draw_layout(buffer, BandWcontext, centre_x - plotting_value / 1.414 - text_width - 2, centre_y - plotting_value / 1.414 - text_height, layout);			
			break;
		case VSWR:
			if (graph->graph_type != SMITH)
				continue;
			plotting_value = fabs((tmp->value - 1.0)/(tmp->value + 1.0))*graph->dependant->scale_factor;
			gdk_draw_arc(buffer, BlueandWcontext, FALSE, centre_x - plotting_value, centre_y - plotting_value, 2.0 * plotting_value, 2.0 * plotting_value, 0, 23040);
			if (text_buffer)
				g_free(text_buffer);
			if (tmp->text)
				text_buffer = g_strdup_printf( "%s\nVSWR: %.1f", tmp->text, tmp->value);
			else
				text_buffer = g_strdup_printf( "VSWR: %.1f", tmp->value);				
			pango_layout_set_text(layout, text_buffer, -1);
			pango_layout_get_pixel_size(layout, &text_width, &text_height);
			break;
		case Q:
			if (graph->graph_type != SMITH)
				continue;
			plotting_value = (1.0 / tmp->value) * graph->dependant->scale_factor;
			r = sqrt(1.0 / (tmp->value * tmp->value) + 1.0) * graph->dependant->scale_factor;
			theta = atan(1.0 / tmp->value) * 180.0 / M_PI * 64;
			gdk_draw_arc(buffer, BlueandWcontext, FALSE, centre_x - r, centre_y + plotting_value - r, 2.0 * r, 2.0 * r, theta, 11520 - 2.0 * theta);
			gdk_draw_arc(buffer, BlueandWcontext, FALSE, centre_x - r, centre_y - plotting_value - r, 2.0 * r, 2.0 * r, 11520 + theta, 11520 - 2.0 * theta);
			if (text_buffer)
				g_free(text_buffer);
			if (tmp->text)
				text_buffer = g_strdup_printf( "%s\nQ: %.1f", tmp->text, tmp->value);
			else
				text_buffer = g_strdup_printf( "Q: %.1f", tmp->value);				
			pango_layout_set_text(layout, text_buffer, -1);
			pango_layout_get_pixel_size(layout, &text_width, &text_height);
			gdk_draw_layout(buffer, BandWcontext, centre_x - text_width / 2, centre_y - plotting_value + r, layout);			
			break;
			
		}
	/* and then move onto the next trace */
	tmp = tmp->next;		
	}
pango_font_description_free(fontdesc);
g_object_unref(layout);	
}
