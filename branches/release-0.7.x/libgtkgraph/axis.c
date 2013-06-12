#include <gtk/gtk.h>
#include <math.h>
#include <stdio.h>

#include "gtkgraph.h"
#include "gtkgraph_internal.h"

/* Local Function Definitions */

static double nicenum(double x, int round);

void gtk_graph_axis_set_tick(GtkGraph *graph, GtkGraphAxisType axis, gfloat majtick, gfloat mintick)
{
GtkGraphAxis *target;
g_return_if_fail (GTK_IS_GRAPH (graph));
	
switch(axis)
	{
	case GTK_GRAPH_AXIS_INDEPENDANT:
		target = graph->independant;
		break;
	case GTK_GRAPH_AXIS_DEPENDANT:
		target = graph->dependant;
		break;
	default:
		return;
		break;
	}

target->maj_tick = majtick;

if (mintick)
	target->min_tick = mintick;
else
	target->min_tick = target->maj_tick;

target->autoscale_tick = FALSE;
}

/* Set_Limits: Set upper and lower limits, as well as tick increments */
void gtk_graph_axis_set_limits(GtkGraph *graph, GtkGraphAxisType axis, gfloat max, gfloat min)
{
GtkGraphAxis *target;
g_return_if_fail (GTK_IS_GRAPH (graph));
	
switch(axis)
	{
	case GTK_GRAPH_AXIS_INDEPENDANT:
		target = graph->independant;
		break;
	case GTK_GRAPH_AXIS_DEPENDANT:
		target = graph->dependant;
		break;
	default:
		return;
		break;
	}

target->axis_max = max;
target->axis_min = min;

target->autoscale_limits = FALSE;
}
/* Allow user to set the value at which the other axis crosses the specified axis */
void gtk_graph_axis_set_crossing(GtkGraph *graph, GtkGraphAxisType axis, GtkGraphCrossingType type, gfloat crossing_value)
{
GtkGraphAxis *target;
g_return_if_fail (GTK_IS_GRAPH (graph));
	
target = graph->independant;				/* Assume it's the independant axis */
if (axis == GTK_GRAPH_AXIS_DEPENDANT)		/* and correct things it its not */
	target = graph->dependant;

g_return_if_fail (target != NULL);

target->crossing_type = type;
target->crossing_value = 0;
if (type == GTK_GRAPH_USERVALUE)
    target->crossing_value = crossing_value;    
}
/* Scale_Axis: A simple Axis Scaler */
void gtk_graph_axis_scale_axis(GtkGraph *graph, gint user_width, gint user_height)
{
GtkGraphTrace *t;
gfloat global_X_max = -1E99, global_Y_max= -1E99, global_X_min=1E99, global_Y_min=1E99;
gint i;

/* --- Do error checking --- */
g_return_if_fail (graph != NULL);
g_return_if_fail (GTK_IS_GRAPH (graph));

t = graph->traces;

if (graph->independant->autoscale_limits || graph->dependant->autoscale_limits)
	{
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
	}

switch (graph->graph_type)
	{
	case XY:
		if (graph->independant->autoscale_limits)		/* we're autoscaling the independant axis limits*/
			{
			graph->independant->axis_max = nicenum(global_X_max, 0);
			graph->independant->axis_min = nicenum(global_X_min, 0);
			}
		if (graph->independant->autoscale_tick)		/* we're autoscaling the independant axis ticks*/
			{
			graph->independant->n_maj_tick = 10;
			graph->independant->n_min_tick = 5;
			graph->independant->maj_tick = (graph->independant->axis_max - graph->independant->axis_min) / graph->independant->n_maj_tick;
			graph->independant->min_tick = graph->independant->maj_tick/ graph->independant->n_min_tick;	
			}
		else
			{
			graph->independant->n_maj_tick = ceil((graph->independant->axis_max - graph->independant->axis_min) / (float) graph->independant->maj_tick);
			graph->independant->n_min_tick = graph->independant->maj_tick / graph->independant->min_tick;
			}
		if (graph->dependant->autoscale_limits)		/* we're autoscaling the dependant axis limits*/
			{
			graph->dependant->axis_max = nicenum(global_Y_max, 0);
			graph->dependant->axis_min = nicenum(global_Y_min, 0);
			}
		if (graph->dependant->autoscale_tick)		/* we're autoscaling the dependant axis ticks*/
			{
			graph->dependant->n_maj_tick = 10;
			graph->dependant->n_min_tick = 5;
			graph->dependant->maj_tick = (graph->dependant->axis_max - graph->dependant->axis_min) / graph->dependant->n_maj_tick;
			graph->dependant->min_tick = graph->dependant->maj_tick/ graph->dependant->n_min_tick;	
			}
		else
			{
			graph->dependant->n_maj_tick = ceil((graph->dependant->axis_max - graph->dependant->axis_min) / (float)graph->dependant->maj_tick);
			graph->dependant->n_min_tick = graph->dependant->maj_tick / graph->dependant->min_tick;
			}

		graph->independant->pxls_per_maj_tick = user_width / graph->independant->n_maj_tick;
		graph->independant->scale_factor = (gfloat) (graph->independant->pxls_per_maj_tick) / graph->independant->maj_tick;
		graph->dependant->pxls_per_maj_tick = user_height / graph->dependant->n_maj_tick;
		graph->dependant->scale_factor = (gfloat) (graph->dependant->pxls_per_maj_tick) / graph->dependant->maj_tick;
		break;
	case POLAR:
		if (graph->dependant->autoscale_limits)		/* we're autoscaling the dependant axis limits*/
			{
			graph->dependant->axis_max = nicenum(global_Y_max, 0);
			graph->dependant->axis_min = nicenum(global_Y_min, 0);
			}
		if (graph->dependant->autoscale_tick)		/* we're autoscaling the dependant axis ticks*/
			{
			graph->dependant->n_maj_tick = 10;
			graph->dependant->n_min_tick = 5;
			graph->dependant->maj_tick = (graph->dependant->axis_max - graph->dependant->axis_min) / graph->dependant->n_maj_tick;
			graph->dependant->min_tick = graph->dependant->maj_tick/ graph->dependant->n_min_tick;	
			}
		else
			{
			graph->dependant->n_maj_tick = ceil((graph->dependant->axis_max - graph->dependant->axis_min) / (float)graph->dependant->maj_tick);
			graph->dependant->n_min_tick = graph->dependant->maj_tick / graph->dependant->min_tick;
			}

		switch(graph->polar_format.type)
			{
			case DEGREES_SYMMETRIC:
				graph->independant->axis_min = -180.0;
				graph->independant->axis_max = 179.9;
				break;
			case DEGREES_ANTISYMMETRIC:
				graph->independant->axis_min = 0.0;
				graph->independant->axis_max = 359.9;
				break;
			 case RADIANS_SYMMETRIC:
				graph->independant->axis_min = -3.141592654;
				graph->independant->axis_max = 3.141592654;
				break;
			 case RADIANS_ANTISYMMETRIC:
				graph->independant->axis_min = 0.0;
				graph->independant->axis_max = 2*3.141592654;
				break;
			 }
		graph->dependant->pxls_per_maj_tick = user_height / (2.0 * graph->dependant->n_maj_tick);
		graph->dependant->scale_factor = (gfloat) (graph->dependant->pxls_per_maj_tick) / graph->dependant->maj_tick;	
		break;
	case SMITH:
		graph->dependant->n_maj_tick = 1.0;
		graph->dependant->maj_tick = 1.0;
		graph->dependant->pxls_per_maj_tick = user_height / (2.0 * graph->dependant->n_maj_tick);
		graph->dependant->scale_factor = (gfloat) (graph->dependant->pxls_per_maj_tick) / graph->dependant->maj_tick;			
		break;
	}
}



/* Format_Axis: Specifiy Axis Labels and Titles */
void gtk_graph_axis_format(GtkGraph *graph, GtkGraphAxisType axis, GtkGraphNumberFormat number_format, gint precision, const gchar *title)
{
GtkGraphAxis *target;
g_return_if_fail (GTK_IS_GRAPH (graph));
	
target = graph->independant;				/* Assume it's the independant axis */
if (axis == GTK_GRAPH_AXIS_DEPENDANT)		/* and correct things it its not */
	target = graph->dependant;

g_return_if_fail (target != NULL);					/* --- Do error checking --- */
target->format = number_format;
target->precision = precision;

if (target->title != NULL)
	g_free(target->title);
target->title = g_strdup(title);
}

/* Format_Grid: Specify whether Grid is visible [and how to draw it (to be implemented)] */

void gtk_graph_axis_format_grid(GtkGraph *graph, GtkGraphAxisType axis, gboolean visible)
{
GtkGraphAxis *target;
g_return_if_fail (GTK_IS_GRAPH (graph));
	
target = graph->independant;						/* Assume it's the independant axis */
if (axis == GTK_GRAPH_AXIS_DEPENDANT)		/* and correct things it its not */
	target = graph->dependant;

g_return_if_fail (target != NULL);   				/* --- Do error checking --- */
target->grid_visible = visible;
}

/* Nicenum: provides nice numerical limits when scaling an axis */

#define expt(a, n) pow(a, (double)(n))

static double nicenum(double x, int round)
{
    int exp;				/* exponent of x */
    double f;				/* fractional part of x */
    double nf;				/* nice, rounded fraction */
	double sign = (x >= 0) ? 1.0 : -1.0;
	double ax = x / sign;	/* Absolute value of x*/

    exp = floor(log10(ax));
    f = ax/expt(10., exp);		/* between 1 and 10 */
    if (round)
		{
		if (f<1.5)
			nf = 1.;
		else if (f<3.)
			nf = 2.;
		else if (f<7.)
			nf = 5.;
		else
			nf = 10.;
		}
    else
		{
		if (f<=1.)
			nf = 1.;
		else if (f<=2.)
			nf = 2.;
		else if (f<=5.)
			nf = 5.;
		else
			nf = 10.;
		}
   return sign*nf*expt(10., exp);
}

#undef expt
