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
from datetime import date
from datetime import timedelta

debug_mode = False
# debug_mode = True

def describeEvent(days, label, off):
    x1 = ar.x_pos(days)
    # can.line(line_style.black_dash1, x1, ybot, x1, ytip)
    # can.line(line_style.black_dash1, x1, ybot, x1, 2)
    tb = text_box.T(text=label, loc=(x1+off, yloc), shadow=(1,-1,fill_style.gray70))
    tb.add_arrow((x1, 0))
    tb.draw()

def db_print(arg, val):
    global debug_mode
    if debug_mode:
        print arg, val
    else:
        pass
    
def predict_release(data):

    cg_prediction = predict_release_by_crossing_graphs(data)
    db_print('%% crossing-graph prediction returned ', cg_prediction)
    if cg_prediction:
       today = cg_prediction[2]
       db_print("%%cg_prediction: ", cg_prediction)
       dpb_prediction = predict_release_from_dev_points(data)
       weight = get_weight(data)
       pr_1 = cg_prediction[0]
       pr_2 = dpb_prediction[0]
       db_print("%% combining pr_1 and pr_2 using weight ", [pr_1, pr_2, weight])
       pred_x_combined = pr_1 * weight + pr_2 * (1.0 - weight)
       pred_y_combined = cg_prediction[1] * weight + dpb_prediction[1] * (1.0 - weight)
       db_print("%% pred_x_combined ", pred_x_combined)
       db_print("%% pred_y_combined ", pred_y_combined)
       return [pred_x_combined, pred_y_combined, today]
    else:
       return False

def get_weight(data):
    # return the weight, 0-1, how much emphasis should be put on the
    # crossing-graphs-based prediction
    # What fraction of scope have we done as of today?
    last = len(data) -1
    c1 = data[last][1] # amount "done" as of today
    c2 = data[last][2] # scope esitmate as of today
    db_print("weight from division of ", [c1, c2])
    weight = c1/c2;
    return weight

def predict_release_from_dev_points(data):
    # scope fraction prediction, that is: how much of the first scope will be the 
    # scope at completion? (Should be 1.0 of course).
    scope_fp = 1.0
    last = len(data) -1
    c2 = data[last][2] # latest scope
    c0 = data[0][2]    #  first scope
    # pred_x = c0 * 4/3.14159
    # pred_y = c0 * 2
    pred_x = c0 * 3.14/2
    pred_y = c0 * scope_fp
    db_print("%% predict_release_from_dev_points: first,last scope: ", [c0, c2])
    db_print("%% predict_release_from_dev_points: day,dev_points ", [pred_x, pred_y])
    return [pred_x, pred_y]

# return either False or a list of [X_pred, Y_pred, today]
#
def predict_release_by_crossing_graphs(data):
    db_print('data', data)
    last = len(data) -1
    if (last > 1):
        c1 = data[0][1] # first how-much-done
        c2 = data[0][2] # first scope
        X_today = data[last][0] # how many days have passed
        db_print('c1 (first how-much done):', c1)
        db_print('c2 (first scope):', c2)
        db_print('X_today:', X_today)
        if X_today > 0:
            db_print("data_stuff: ", (data[last][1], c1, X_today))
            m1 = (data[last][1] - c1)/X_today
            m2 = (data[last][2] - c2)/X_today
            m_diff = m1 - m2
            db_print('m1  ("done" gradient):', m1) # 'done' gradient
            db_print('m2 ("scope" gradient):', m2) # scope gradient
            db_print('m_diff:', m_diff)
            if m_diff > 0:
                X_pred = (c2-c1)/(m1-m2) 
                Y_pred_1 = m1 * X_pred + c1
                Y_pred_2 = m2 * X_pred + c2
                return [X_pred+1, Y_pred_1, X_today]
            else:
                return False
        else:
            return False
    else:
        return False

def x_ticks(x_day_range):
   x_tick_interval = 2
   if (x_day_range >= 40):
      x_tick_interval = 5
   if (x_day_range > 80):
      x_tick_interval = 10
   if (x_day_range > 120):
      x_tick_interval = 20
   if (x_day_range > 220):
      x_tick_interval = 25
   if (x_day_range > 280):
      x_tick_interval = 50
   return x_tick_interval

def draw_prediction_on_canvas(pr, can):

   if pr:
      X_pred = pr[0]
      Y_pred = pr[1]
      X_today = pr[2]
      t_delta = timedelta(days=X_pred-X_today)
      now = date.today()
      pred_date = now + t_delta
      pred_date_str = pred_date.strftime("%A %d %b %Y")
      
      s = "/4/oProjected Release: Day "
      int_X_pred = int(X_pred+1)
      s += str(int_X_pred)
      s += ",\n"
      s += pred_date_str
      # can.show(70,90, s)
      # can.show(100,90, s)
      can.show(28, 105, s)

      #   print "X_pred: :", X_pred, "Y_pred: :", Y_pred
      can.rectangle(line_style.default, fill_style.default,
                    x1=ar.x_pos(int_X_pred-0.3), x2=ar.x_pos(int_X_pred+0.3),
                    y1=ar.y_pos(Y_pred-0.5), y2=ar.y_pos(Y_pred+0.5))
      can.line(line_style.black_dash1, ar.x_pos(int_X_pred), 0,
               ar.x_pos(int_X_pred), ar.y_pos(Y_pred))

# Annotations:

def annotation_box(ar, box_text, loc_x, loc_y, done_or_scope, data_index, arrow_position):
    db_print("in annotation_box: data is ", data)
    tb = text_box.T(loc=(loc_x, loc_y), text=box_text,
                    shadow=(1,-1,fill_style.gray70), bottom_fudge=2)
    tb.add_arrow((ar.x_pos(data[data_index][0]),
                  ar.y_pos(data[data_index][done_or_scope])),
                 arrow_position)
    tb.draw()



#############################################################################

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
x_label = "Number of Days (since starting release development)"


# The format attribute specifies the text to be drawn at each tick mark.
# Here, texts are rotated -60 degrees ("/a-60"), left-aligned ("/hL"),
# and numbers are printed as integers ("%d").
#

x_day_range = 100

x_tick_interval = x_ticks(x_day_range)

xaxis = axis.X(tic_interval = x_tick_interval, label=x_label)
yaxis = axis.Y(tic_interval = 20, label="Dev Points")

pr = predict_release(data)

if debug_mode == False:

    # Define the drawing area. "y_range=(0,None)" tells that the Y minimum
    # is 0, but the Y maximum is to be computed automatically. Without
    # y_ranges, Pychart will pick the minimum Y value among the samples,
    # i.e., 20, as the base value of Y axis.
    ar = area.T(x_axis=xaxis, y_axis=yaxis, x_range=(0,x_day_range), y_range=(0,220))

    # The first plot extracts Y values from the 2nd column
    # ("ycol=1") of DATA ("data=data"). X values are takes from the first
    # column, which is the default.
    plot = line_plot.T(label="Done", data=data, ycol=1)
    # plot2 = line_plot.T(label="Total", data=data, ycol=2, tick_mark=tick_mark.square)
    plot2 = line_plot.T(label="Total Scope for 0.8.9.2", data=data, ycol=2)

    ar.add_plot(plot, plot2)
    draw_prediction_on_canvas(pr, can)
    # The call to ar.draw() usually comes at the end of a program.  It
    # draws the axes, the plots, and the legend (if any).
    #
    ar.draw()
        
    yloc = ar.loc[1] + ar.size[1] + 50
    ytip = ar.loc[1] + ar.size[1]
    ybot = ar.loc[1]
 
    yloc = 20
    ybot = 0
    theme.default_font_size=4

    # annotation_box(ar, "The Great CSHL Purge",       80, 90,   2, 12, "lt")


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


