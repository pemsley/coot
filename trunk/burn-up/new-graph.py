#
# Copyright (C) 2000-2005 by Yasushi Saito (yasushi.saito@gmail.com)
# Copyright 2009 by The University of Oxford
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
theme.scale_factor=5
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

x_day_range = 355

x_tick_interval = 2
if (x_day_range > 40):
   x_tick_interval = 5
if (x_day_range > 80):
   x_tick_interval = 10
if (x_day_range > 120):
   x_tick_interval = 20
if (x_day_range > 220):
   x_tick_interval = 25
if (x_day_range > 280):
   x_tick_interval = 50

xaxis = axis.X(tic_interval = x_tick_interval, label="Days (since pre-release start)")
yaxis = axis.Y(tic_interval = 20, label="Dev Points")

# Define the drawing area. "y_range=(0,None)" tells that the Y minimum
# is 0, but the Y maximum is to be computed automatically. Without
# y_ranges, Pychart will pick the minimum Y value among the samples,
# i.e., 20, as the base value of Y axis.
ar = area.T(x_axis=xaxis, y_axis=yaxis, x_range=(0,x_day_range), y_range=(0,280))

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
    
describeEvent(5, "Christmas\nholidays start", -20)

xpscale = 0.6

# crisis:
# crisis_pt=59 280 days?
# crisis_pt=57
crisis_pt=50
can.line(line_style.black_dash1, crisis_pt, ybot, crisis_pt, ytip)
tb = text_box.T(text="Crisis!", loc=(25, 110), shadow=(1,-1,fill_style.gray70), bottom_fudge=2)
tb.add_arrow((crisis_pt, ytip))
tb.draw()
    
tb = text_box.T(loc=(18*xpscale, 43), text="Boston", shadow=(1,-1,fill_style.gray70), bottom_fudge=2)
# tb.add_arrow((ar.x_pos(data[33][0]), ar.y_pos(data[33][1])), "cb")
tb.add_arrow((ar.x_pos(data[33][0]), ar.y_pos(data[33][1])), "c")
tb.draw()

tb = text_box.T(loc=(20*xpscale, 90), text="Papers & Happy", shadow=(1,-1,fill_style.gray70), bottom_fudge=3)
tb.add_arrow((ar.x_pos(data[74][0]), ar.y_pos(data[74][1])), "c")
tb.draw()

tb = text_box.T(loc=(68*xpscale, 60), text="CCP4 Dev", shadow=(1,-1,fill_style.gray70), bottom_fudge=2)
tb.add_arrow((ar.x_pos(data[115][0]), ar.y_pos(data[115][1])), "c")
tb.draw()

tb = text_box.T(loc=(50*xpscale, 8), text="Dutch Meeting", shadow=(1,-1,fill_style.gray70), bottom_fudge=3)
tb.add_arrow((ar.x_pos(data[122][0]), ar.y_pos(data[122][1])), "tc")
tb.draw()

tb = text_box.T(loc=(100*xpscale, 20), text="NIH & MAMCM", shadow=(1,-1,fill_style.gray70), bottom_fudge=3)
tb.add_arrow((ar.x_pos(data[160][0]), ar.y_pos(data[160][1])), "tc")
tb.draw()

tb = text_box.T(loc=(130*xpscale, 31), text="Paper Submitted", shadow=(1,-1,fill_style.gray70), bottom_fudge=3)

tb.add_arrow((ar.x_pos(data[206][0]), ar.y_pos(data[206][1])), "tc")
tb.draw()

#tb = text_box.T(loc=(170*xpscale, 45), text="Havana", shadow=(1,-1,fill_style.gray70), bottom_fudge=3, bg_style=fill_style.gray90)
tb = text_box.T(loc=(155*xpscale, 42), text="Havana", shadow=(1,-1,fill_style.gray70), bottom_fudge=3)
tb.add_arrow((ar.x_pos(data[262][0]), ar.y_pos(data[262][1])), "tc")
tb.draw()

tb = text_box.T(loc=(168*xpscale, 54), text="Boston", shadow=(1,-1,fill_style.gray70), bottom_fudge=3)
tb.add_arrow((ar.x_pos(data[343][0]), ar.y_pos(data[343][1])), "tc")
tb.draw()

tb = text_box.T(loc=(110*xpscale, 110), text="User Feedback", shadow=(1,-1,fill_style.gray70), bottom_fudge=3)
tb.add_arrow((ar.x_pos(data[366][0]), ar.y_pos(data[366][1])), "c")
tb.draw()

tb = text_box.T(loc=(186*xpscale, 63), text="HEC-2009", shadow=(1,-1,fill_style.gray70), bottom_fudge=3)
tb.add_arrow((ar.x_pos(data[414][0]), ar.y_pos(data[414][1])), "tc")
tb.draw()

tb = text_box.T(loc=(201*xpscale, 72), text="CSH-2009", shadow=(1,-1,fill_style.gray70), bottom_fudge=3)
tb.add_arrow((ar.x_pos(data[440][0]), ar.y_pos(data[440][1])), "tc")
tb.draw()


# take-home:
#   Max sustainable rate: 2.2 dev-pts/day
#  Real (long-term) rate: 0.7 dev-pts/day
# Unmassaged scope creep: 0.8 dev-pts/day
#  -> Oh dear.
#
# Lessons:
#
# for 100-day release and same scope creep
#    -> impossible
# for 100-day release cycle and no scope creep
#    -> spec for 70 dev points
# for 100-day release cycle and 50% reduction in scope creep
#    -> spec for 30 dev points
#
#
# given that, 160 points (even cut down to 120 pts) was laughable.


