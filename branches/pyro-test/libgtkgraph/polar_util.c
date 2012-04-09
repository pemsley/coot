#include <math.h>
#include "gtkgraph.h"
#include "gtkgraph_internal.h"

gfloat wrap_angle(gfloat angle, GtkGraphPolarFormat format)
{
gint num_cycles;
gfloat new_angle;

switch (format.type)
	{
	case DEGREES_SYMMETRIC:
      	if (angle < -180.0 || angle >= 180)
			{
            num_cycles = (gint) floor((angle + 180) / 360.0);
            new_angle = angle - num_cycles * 360.0;
            return new_angle;
            }
		break;
	case DEGREES_ANTISYMMETRIC:
       	if (angle < 0.0 || angle >= 360)
          	{
            num_cycles = (gint) floor(angle / 360.0);
            new_angle = angle - num_cycles * 360.0;
            return new_angle;
            }
		break;
	case RADIANS_SYMMETRIC:

		if (angle < -M_PI || angle >= M_PI)
          	{
            num_cycles = (gint) ((angle + M_PI) / TWOPI);
            new_angle = angle - (gfloat) num_cycles * TWOPI;
            return new_angle;
            }
		break;
	case RADIANS_ANTISYMMETRIC:
      	if (angle < 0.0 || angle >= TWOPI)
         	{
            num_cycles = (gint) floor(angle / TWOPI);
            new_angle = angle - (gfloat) num_cycles * TWOPI;
            return new_angle;
            }
		break;

	default:
		return angle;
		break;
	}
return angle;
}
