/*
 * (c) Copyright 1993, 1994, Silicon Graphics, Inc.
 * ALL RIGHTS RESERVED
 * Permission to use, copy, modify, and distribute this software for
 * any purpose and without fee is hereby granted, provided that the above
 * copyright notice appear in all copies and that both the copyright notice
 * and this permission notice appear in supporting documentation, and that
 * the name of Silicon Graphics, Inc. not be used in advertising
 * or publicity pertaining to distribution of the software without specific,
 * written prior permission.
 *
 * THE MATERIAL EMBODIED ON THIS SOFTWARE IS PROVIDED TO YOU "AS-IS"
 * AND WITHOUT WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR OTHERWISE,
 * INCLUDING WITHOUT LIMITATION, ANY WARRANTY OF MERCHANTABILITY OR
 * FITNESS FOR A PARTICULAR PURPOSE.  IN NO EVENT SHALL SILICON
 * GRAPHICS, INC.  BE LIABLE TO YOU OR ANYONE ELSE FOR ANY DIRECT,
 * SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY
 * KIND, OR ANY DAMAGES WHATSOEVER, INCLUDING WITHOUT LIMITATION,
 * LOSS OF PROFIT, LOSS OF USE, SAVINGS OR REVENUE, OR THE CLAIMS OF
 * THIRD PARTIES, WHETHER OR NOT SILICON GRAPHICS, INC.  HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH LOSS, HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE
 * POSSESSION, USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 * US Government Users Restricted Rights
 * Use, duplication, or disclosure by the Government is subject to
 * restrictions set forth in FAR 52.227.19(c)(2) or subparagraph
 * (c)(1)(ii) of the Rights in Technical Data and Computer Software
 * clause at DFARS 252.227-7013 and/or in similar or successor
 * clauses in the FAR or the DOD or NASA FAR Supplement.
 * Unpublished-- rights reserved under the copyright laws of the
 * United States.  Contractor/manufacturer is Silicon Graphics,
 * Inc., 2011 N.  Shoreline Blvd., Mountain View, CA 94039-7311.
 *
 * OpenGL(TM) is a trademark of Silicon Graphics, Inc.
 */
/*
 * Trackball code:
 *
 * Implementation of a virtual trackball.
 * Implemented by Gavin Bell, lots of ideas from Thant Tessman and
 *   the August '88 issue of Siggraph's "Computer Graphics," pp. 121-129.
 *
 * Original code from:
 * David M. Ciemiewicz, Mark Grossman, Henry Moreton, and Paul Haeberli
 *
 * Much mucking with by:
 * Gavin Bell
 * GLM conversion: Paul Emsley
 */

#if defined (USE_PYTHON)
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#include "python-3-interface.hh"
#endif

#include "graphics-info.h"

// static
float
graphics_info_t::trackball_project_to_sphere(float r, float x, float y) {

    float oor2 = 0.70710678118654752440;
    float d = sqrt(x*x + y*y);
    if (d < r * oor2) {    /* Inside sphere */
        return sqrt(r*r - d*d);
    } else {
       /* On hyperbola */
       float t = r * oor2;
       return t*t/d;
    }
}

// static
glm::quat
graphics_info_t::trackball_to_quaternion(float p1x, float p1y, float p2x, float p2y, float trackball_size) {

   if (p1x == p2x && p1y == p2y) {
        /* Zero rotation */
        return glm::quat(1,0,0,0);
    }

    float d_mult = trackball_size - 0.3; /* or some such */

    /*
     * First, figure out z-coordinates for projection of P1 and P2 to
     * deformed sphere
     */
    glm::vec3 p1(p1x, p1y, trackball_project_to_sphere(trackball_size, p1x, p1y));
    glm::vec3 p2(p2x, p2y, trackball_project_to_sphere(trackball_size, p2x, p2y));

    /*
     *  Now, we want the cross product of P1 and P2
     */
    // vcross(p2,p1,a);
    glm::vec3 a = glm::normalize(glm::cross(p2,p1));

    /*
     *  Figure out how much to rotate around that axis.
     */
    glm::vec3 d(p1 - p2);

    float t = glm::length(d) * d_mult / (trackball_size);

    /*
     * Avoid problems with out-of-control values...
     */
    if (t >  1.0f) t =  1.0f;
    if (t < -1.0f) t = -1.0f;

    /* how much to rotate about axis */
    float phi = 2.0f * asin(t);

    glm::vec3 qaaa(a * sinf(phi/2.0f));
    glm::quat q(cosf(phi/2.0f), qaaa.x, qaaa.y, qaaa.z);
    return q;
}
