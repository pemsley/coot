#
# Copyright (C) 2000-2005 by Yasushi Saito (yasushi.saito@gmail.com)
# 
# Pychart is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2, or (at your option) any
# later version.
#
# Pychart is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
from pychart import *
theme.get_options()
theme.output_format="png"
theme.scale_factor=4.2
theme.default_font_size=6
theme.reinitialize()

can = canvas.default_canvas()

# We have 10 sample points total.  The first value in each tuple is
# the X value, and subsequent values are Y values for different lines.
#
data = chart_data.read_csv("burn-up.tab", delim=" ")

# The format attribute specifies the text to be drawn at each tick mark.
# Here, texts are rotated -60 degrees ("/a-60"), left-aligned ("/hL"),
# and numbers are printed as integers ("%d").
#

x_day_range = 180

x_tick_interval = 2
if (x_day_range > 40):
   x_tick_interval = 5
if (x_day_range > 80):
   x_tick_interval = 10
if (x_day_range > 120):
   x_tick_interval = 20

xaxis = axis.X(tic_interval = x_tick_interval, label="Days (since pre-release start)")
yaxis = axis.Y(tic_interval = 20, label="Dev Points")

# Define the drawing area. "y_range=(0,None)" tells that the Y minimum
# is 0, but the Y maximum is to be computed automatically. Without
# y_ranges, Pychart will pick the minimum Y value among the samples,
# i.e., 20, as the base value of Y axis.
ar = area.T(x_axis=xaxis, y_axis=yaxis, x_range=(0,x_day_range), y_range=(0,240))

# The first plot extracts Y values from the 2nd column
# ("ycol=1") of DATA ("data=data"). X values are takes from the first
# column, which is the default.
plot = line_plot.T(label="Done", data=data, ycol=1)
# plot2 = line_plot.T(label="Total", data=data, ycol=2, tick_mark=tick_mark.square)
plot2 = line_plot.T(label="Total Scope for 0.6", data=data, ycol=2)

ar.add_plot(plot, plot2)

# The call to ar.draw() usually comes at the end of a program.  It
# draws the axes, the plots, and the legend (if any).

ar.draw()


yloc = ar.loc[1] + ar.size[1] + 50
ytip = ar.loc[1] + ar.size[1]
ybot = ar.loc[1]

yloc = 20
ybot = 0
theme.default_font_size=4


def describeEvent(days, label, off):
    x1 = ar.x_pos(days)
    # can.line(line_style.black_dash1, x1, ybot, x1, ytip)
    # can.line(line_style.black_dash1, x1, ybot, x1, 2)
    tb = text_box.T(text=label, loc=(x1+off, yloc), shadow=(1,-1,fill_style.gray70))
    tb.add_arrow((x1, 0))
    tb.draw()
    
describeEvent(5, "Christmas\nholidays start", -10)

xpscale = 0.9

tb = text_box.T(loc=(21*xpscale, 43), text="Boston", shadow=(1,-1,fill_style.gray70))
# tb.add_arrow((ar.x_pos(data[33][0]), ar.y_pos(data[33][1])), "cb")
tb.add_arrow((ar.x_pos(data[33][0]), ar.y_pos(data[33][1])), "c")
tb.draw()

tb = text_box.T(loc=(40*xpscale, 73), text="Papers & Happy", shadow=(1,-1,fill_style.gray70))
tb.add_arrow((ar.x_pos(data[74][0]), ar.y_pos(data[74][1])), "c")
tb.draw()

tb = text_box.T(loc=(80*xpscale, 72), text="CCP4 Dev", shadow=(1,-1,fill_style.gray70))
tb.add_arrow((ar.x_pos(data[115][0]), ar.y_pos(data[115][1])), "c")
tb.draw()

tb = text_box.T(loc=(85*xpscale, 22), text="Dutch Meeting", shadow=(1,-1,fill_style.gray70))
tb.add_arrow((ar.x_pos(data[122][0]), ar.y_pos(data[122][1])), "tc")
tb.draw()

