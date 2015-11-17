
import os
import coot
import math
import cairo
import cairoplot

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

# import matplotlib new

def pathology_plots(mtz, fp, sigfp):
    try: 
        print "getting data for", mtz
        data = coot.pathology_data(mtz, fp, sigfp)
        # print len(data), data
        prefix,tail = os.path.splitext(mtz)
        for i in range(4):
            png_file_name = prefix + "-" + str(i) + ".png"
            print "making plot", png_file_name
            yt = 'I'
            xt = 'Resolution'
            if i == 1:
                yt = 'I/SIGI'
            if i == 2:
                yt = 'SIGI'
            if i == 3:
                yt = 'I/SIGI'
                xt = "I"
            cairoplot.scatter_plot(png_file_name,
                                   data = data[i+1],
                                   # series_colors = ['black'],
                                   series_colors = [[0.2, 0.2, 0.2]],
                                   width = 600, height = 600, border = 20,
                                   discrete = True,
                                   x_title = xt,
                                   y_title = yt,
                                   dots = True,
                                   axis = True)

    except TypeError as e:
                print "caught TypeError:", e


def new_pathology_plots(mtz, fp, sigfp):
    data = coot.pathology_data(mtz, fp, sigfp)
    try: 
        prefix,tail = os.path.splitext(mtz)
        labels = [("Resolution", "FP"),
                  ('Resolution', 'FP/SIGFP'),
                  ('FP', 'SIGFP'),
                  ('FP', 'FP/SIGFP')]
        z = zip(range(4), labels)
        invresolsq_max = data[0]
        for ii in range(4):
            i = ii+1 # data and labels are indexed differently
            l_0 = labels[ii][0].replace('/', '_')
            l_1 = labels[ii][1].replace('/', '_')
            png_file_name = prefix + "-" + l_0 + "-vs-" + l_1 + ".png"
            n_x_labels = 5
            x_labels = ["low", '.', "medium", '.', "high"]

            if l_0 == "Resolution":
                for j in range(n_x_labels):
                    f = invresolsq_max * j / n_x_labels
                    if f == 0:
                        x_label = "Inf"
                    else:
                        x_label = str(round(math.sqrt(1/f), 3));
                    print j, x_label
                    x_labels[j] = x_label
            else:
                x_labels = None
            
            print "making plot", png_file_name

            fig1 = plt.figure()
            xdata = [x[0] for x in data[i+1]]
            ydata = [y[1] for y in data[i+1]]
            # l = plt.scatter(xdata, ydata, 'r-')
            s = 0.1
            l = plt.scatter(xdata, ydata, s)
            plt.savefig(png_file_name)
            

    except TypeError as e:
                print "caught TypeError:", e

