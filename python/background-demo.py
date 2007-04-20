# Makefile.am
# 
# Copyright 2005, 2006, 2007 The University of York
# Copyright 2005, 2006, 2007 Bernhard Lohkamp
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA

#;; A silly (but somewhat amusing) demo that changes the background colour
#;; over a range of colours.
#;;
#;; The scene is spun too, for extra eye-candy (so have a molecule
#;; loaded before you run this).
#;;
#;; What's the point?  It's to show that we can display the result of
#;; an arbitarily complex computation, i.e, we have a real extension
#;; language, not just a list of commands (like some other molecular
#;; graphics programs).

# n-steps is divided into 3 sets:
n_steps = 300.0

def background_demo():

    def rgb_colours(nlocal):

	r_step = 1 / n_steps
        red = 0.0
        green =0.0
        blue = 0.0

        if (n_local < (n_steps / 3)):
            red = 3 * r_step * n_local
            green = 0.0
        elif (n_local < ((2 * n_steps) / 3)):
            red = 1 - ((n_local - (n_steps / 3)) * 3 * r_step) 
            green = ((n_local - (n_steps / 3)) * 3 * r_step)
        else:
            red = 0.0
            green = 1 - (n_local - ((n_steps / 3) * 2)) * 3 * r_step

	return red,green,blue


    n_local = 1.0
    while (n_local<=n_steps):
	      r_col = rgb_colours(n_local)
# 	      print "colours ",r_col
	      set_background_colour(r_col[0],r_col[1],r_col[2])
	      rotate_y_scene(10,0.1) # n-frames frame-interval(degrees)
	      n_local = n_local + 1.0

background_demo()


