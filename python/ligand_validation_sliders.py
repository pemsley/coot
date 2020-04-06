#  Copyright 2016 by Medical Research Council
#  Author: Paul Emsley
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or (at
#  your option) any later version.
#  
#  This program is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
#  02110-1301, USA

import gtk
import pango
import math
import string
import types

class ligand_validation_metrics_gui:

    def __init__(self, ligand_validation_specs):

        window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        title = "Ligand Validation Report for " # ...
        title += ligand_validation_specs.title

        window.set_title(title)
        self.vbox = gtk.VBox(False, 0)
        self.vbox.set_border_width(5)
        h_sep = gtk.HSeparator()
	self.bar_length = 300
	self.bar_height = 10
	self.abs_bar_height = int(self.bar_height * 1.6)
	self.x_bar_offset = 160
	self.y_bar_offset =  20
	self.y_pixels_bar_spacing = 30
	self.abs_bar_width = 6
        da = gtk.DrawingArea()
        da.set_size_request(560,280)
        close_button = gtk.Button("  Close  ")
        hbox = gtk.HBox(False, 0)
        self.vbox.pack_start(da,         False, 6)
        self.vbox.pack_start(h_sep,      False, 6)
        hbox.pack_start(close_button,    True, 6)
        self.vbox.pack_end(hbox,         False, 6)
        window.add(self.vbox)
        window.show_all()
        da.connect("expose-event", self.on_drawing_area_expose)
        close_button.connect("clicked", lambda a : window.destroy())
        self.pangolayout = da.create_pango_layout("")
        style = da.get_style()
	gc = style.fg_gc[gtk.STATE_NORMAL]
	self.setup_colour_bar_buff()
        self.lvs = ligand_validation_specs

    def setup_colour_bar_buff(self):

	c = 3*self.bar_length*self.bar_height*['\0']
	for j in range(self.bar_length):

	    f_j = float(j)/float(self.bar_length)
	    # we need g to go 255:255:0
	    r = int(255*(1-math.pow(f_j, 5)))
	    # we need g to go 0:255:0
	    g_1 = f_j                #  0 : 0.5 : 1
	    g_2 = 2 * (g_1 - 0.5)    # -1 : 0.  : 1
	    g_3 = g_2 * g_2          #  1 : 0.  : 1
	    g_4 = 1 - g_3            #  0 : 1.  : 0
	    g = int(240*g_4)
	    b = int(255*math.pow(f_j, 0.2))
	    for i in range(self.bar_height):
		idx = 3*(self.bar_length*i + j)
		c[idx  ] = chr(r)
		c[idx+1] = chr(g)
		c[idx+2] = chr(b)
	self.colour_bar_buff = string.join(c, '')

    def draw_sliders(self, ligand_validation_specs, da, gc):

        for idx in range(len(ligand_validation_specs.metrics)):
            m = ligand_validation_specs.metrics[idx]
            # print 'render metric', m.metric_name, m.percentile_abs
            self.draw_slider(m.metric_name, 0, m.percentile_abs, m.percentile_rel, '{0:.3f}'.format(m.value), idx, da, gc)

	values = [m.percentile_abs for m in ligand_validation_specs.metrics]
	self.draw_connecting_lines(values, da, gc)

        # perhaps we should check that the ligand metrics were (each good)?
        return len(ligand_validation_specs.metrics)

                 
    # return True if the bar for the absolute percentile was drawn,
    # otherwise return False
    #
    def draw_slider(self, name, x_for_rj, percentile_abs_str, percentile_rel_str, value_str, slider_no, da, gc):

	y = self.y_bar_offset + self.y_pixels_bar_spacing*(slider_no+1)

	# colour bar
	self.draw_rgb_image(da, gc, self.x_bar_offset, y)
	
	local_gc = gc;
	local_gc.set_foreground(local_gc.get_colormap().alloc_color("#888888"))
	
	da.window.draw_rectangle(local_gc, False, self.x_bar_offset, y, self.bar_length, self.bar_height)
	local_gc.set_foreground(local_gc.get_colormap().alloc_color("#000000"))

	# Metric text
	pangolayout = da.create_pango_layout("")
	pangolayout.set_justify(1)
	pangolayout.set_alignment(pango.ALIGN_RIGHT)
        pangolayout.set_text(name)
	pangolayout.context_changed()

        da.window.draw_layout(gc, 4+x_for_rj, y-3, pangolayout)

	# Values text
	if isinstance(value_str, bytes):
	    x_for_value = self.x_bar_offset + self.bar_length + 12
	    pangolayout.set_text(value_str)
	    # print "Drawing Value value", x_for_value, y, value_str
	    da.window.draw_layout(gc, x_for_value, y-6, pangolayout)

	# Bars for percentile scores
	try:
	    if percentile_rel_str != False:
		rel = float(percentile_rel_str)
		self.arcs_for_rel(rel, y, da, gc)
	    if percentile_abs_str != False:
		abs = float(percentile_abs_str)
		self.bar_for_abs(abs, y, da, gc)
		return True
		
	    rel = float(rel_str)
	except TypeError as e:
	    print(e)

	# hopefully we don't get here.
	return False

    # pass abs values
    def draw_connecting_lines(self, pc_ranks, da, gc):

	if len(pc_ranks) > 1:
	    for slider_no in range(len(pc_ranks)-1):
		y_min_1 = 20 + 30*(slider_no+1)
		y_min_2 = 20 + 30*(slider_no+2)
		abs_1 = float(pc_ranks[slider_no])
		abs_2 = float(pc_ranks[slider_no+1])
		
		x_1 = int(self.x_bar_offset + 0.01 * abs_1 * self.bar_length)
		y_1 = int(y_min_1 + self.bar_height + 1)
		
		x_2 = int(self.x_bar_offset + 0.01 * abs_2 * self.bar_length)
		y_2 = int(y_min_2 - self.bar_height * 0.25)
		
		da.window.draw_line(gc, x_1, y_1, x_2, y_2)
	
    def on_drawing_area_expose(self, drawing_area, event):
        style = drawing_area.get_style()
	gc = style.fg_gc[gtk.STATE_NORMAL]
        n_sliders = self.draw_sliders(self.lvs, drawing_area, gc)
	self.draw_top_labels(drawing_area, gc)
	self.draw_bottom_labels(drawing_area, gc, n_sliders) # Worse, Better
        self.draw_key(drawing_area, gc, n_sliders) # percentile box descriptions

    def draw_rgb_image(self, da, gc, x, y):
	da.window.draw_rgb_image(gc, x, y, self.bar_length, self.bar_height,
				 gtk.gdk.RGB_DITHER_NONE, self.colour_bar_buff, self.bar_length*3)
        
    def bar_for_abs(self, abs_percent, y_min, da, gc):

	x = int(self.x_bar_offset + 0.01 * abs_percent * self.bar_length - 0.5 * self.abs_bar_width)
	y = int(y_min - self.bar_height * 0.25)
	da.window.draw_rectangle(gc, True, x, y, self.abs_bar_width, self.abs_bar_height)
	
    def arcs_for_rel(self, rel_percent, y_min, da, gc):
	x = int(self.x_bar_offset + 0.01 * rel_percent * self.bar_length - 0.5 * self.abs_bar_width)
	y = int(y_min - self.bar_height * 0.25)
	da.window.draw_rectangle(gc, False, x, y, self.abs_bar_width, self.abs_bar_height)

    def draw_top_labels(self, da, gc):
	pangolayout = da.create_pango_layout("")

        pangolayout = da.create_pango_layout("")
        font_desc = pango.FontDescription('Sans 13')
	pangolayout.set_font_description(font_desc)
	
        pangolayout.set_text('Metric')
	y_level = 15
        da.window.draw_layout(gc, 80, y_level, pangolayout)
        pangolayout.set_text('Percentile Ranks')
        da.window.draw_layout(gc, 245, y_level, pangolayout)
        pangolayout.set_text('Value')
        da.window.draw_layout(gc, 470, y_level, pangolayout)

    # Worse, Better
    def draw_bottom_labels(self, da, gc, n_sliders): 
	x_worse  = self.x_bar_offset
	y_worse  = self.y_bar_offset + self.y_pixels_bar_spacing * n_sliders + 18
	x_better = self.x_bar_offset + self.bar_length - 32
	y_better = y_worse

        pl_wb = da.create_pango_layout("")
        pl_wb.set_text('Worse')
        font_desc = pango.FontDescription('Sans 8')
	font_desc.set_style(pango.STYLE_ITALIC)
	pl_wb.set_font_description(font_desc)
        da.window.draw_layout(gc, x_worse,  y_worse,  pl_wb)
	
        pl_wb.set_text('Better')
        da.window.draw_layout(gc, x_better, y_better, pl_wb)
	
    
    # percentile box descriptions
    def draw_key(self, da, gc, n_sliders):
	x_key_box_abs =  self.x_bar_offset
	y_key_box_abs = self.y_bar_offset + self.y_pixels_bar_spacing * (n_sliders + 1) + 10

	x_key_box_rel = x_key_box_abs
	y_key_box_rel = y_key_box_abs + 20

	x_key_1 =  x_key_box_abs + 10
	y_key_1  = y_key_box_abs 
	
	x_key_2 =  x_key_1
	y_key_2  = y_key_1 + 20

	da.window.draw_rectangle(gc, True, x_key_box_abs, y_key_box_abs,
				 self.abs_bar_width, self.abs_bar_height)
	
        pl = da.create_pango_layout("")
        font_desc = pango.FontDescription('Sans 9')
	pl.set_font_description(font_desc)
	
        pl.set_text('Percentile relative to all x-ray structures')
        da.window.draw_layout(gc, x_key_1, y_key_1, pl)
	
	#da.window.draw_rectangle(gc, False, x_key_box_rel, y_key_box_rel,
	#			 self.abs_bar_width, self.abs_bar_height)
        # pl.set_text('Percentile relative to x-ray structures of similar resolution')
        # da.window.draw_layout(gc, x_key_2, y_key_2, pl)

class ligand_stats:

    class stat_t:
        # percentile is in 0-100 range, can be strings, 
        def __init__(self, metric_name, pc_metric_value_abs, pc_metric_value_rel, value):
            self.metric_name      = metric_name
            self.percentile_abs = pc_metric_value_abs
            self.percentile_rel = pc_metric_value_rel
            self.value       = value

    def __init__(self):
        self.title = 'Test-ligand'
        self.metrics = []

    def add_metric(self, metric_name, pc_value_abs, pc_value_rel, value):
        s = self.stat_t(metric_name, pc_value_abs, pc_value_rel, value)
        self.metrics.append(s)

    # ligand stats are not very resolution dependent, so a
    # ligand_stats metric can be generated without
    # resolution-dependent (relative) percentiles
    #
    def add_metric(self, metric_name, pc_value_abs, value):
        s = self.stat_t(metric_name, pc_value_abs, False, value)
        self.metrics.append(s)

# we want to be able to call this function from scheme.  So pass a list of stats
#
def ligand_validation_metrics_gui_list_wrapper(stats_list):

    ls = ligand_stats()
    for metric in stats_list:
	metric_name = metric[0]
	percentile  = metric[1]
	value       = metric[2]
	ls.add_metric(metric_name, percentile, value)
    ligand_validation_metrics_gui(ls)

# if __name__ == '__main__':

#     ls = ligand_stats()
#     ls.add_metric("metric-name-1", 33, 45, '0.78')
#     ls.add_metric("metric-name-2", 30, 66, '0.93')
#     ls.add_metric("metric-name-3", 53, 45, '22.78')
#     ls.add_metric("metric-name-4", 23, 25, '1.78')
#     ligand_validation_metrics_gui(ls)

#     gtk.main()
