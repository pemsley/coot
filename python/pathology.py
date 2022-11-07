
import os
import coot
import math

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt


def cairoplot_pathology_plots(mtz, fp, sigfp):
    import cairo
    import cairoplot
    try:
        print("getting data for", mtz)
        data = coot.pathology_data(mtz, fp, sigfp)
        # print len(data), data
        prefix,tail = os.path.splitext(mtz)
        for i in range(4):
            png_file_name = prefix + "-" + str(i) + ".png"
            print( "making plot", png_file_name)
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
                print( "caught TypeError:", e)

# these can be I, SIGI also
#
def pathology_plots(mtz, fp, sigfp):

# these can be I, SIGI also
#
def plots(mtz, fp, sigfp):

    # meta data
    intensity_data = False
    if fp[0] == 'I':
        intensity_data = True

    data = pathology.PathologyData(mtz, fp, sigfp)
    # try:
    if True:
        prefix,tail = os.path.splitext(mtz)
        fp_lab = "FP"
        sigfp_lab = "SIGFP"
        if (intensity_data):
            fp_lab = "I"
            sigfp_lab = "SIGI"

        labels = [("Resolution", fp_lab),
                  ('Resolution', fp_lab + '/' + sigfp_lab),
                  (fp_lab, sigfp_lab),
                  (fp_lab, fp_lab + '/' + sigfp_lab)]
        z = zip(range(4), labels)
        if True:
            ii = 0 # label indexing hack
            l_0 = labels[ii][0].replace('/', '_')
            l_1 = labels[ii][1].replace('/', '_')
            png_file_name = prefix + "-" + l_0 + "-vs-" + l_1 + ".png"

            n_x_labels = 5
            """
            if l_0 == "Resolution":
                for j in range(n_x_labels):
                    f = data.invresolsq_max() * j / n_x_labels
                    if f == 0:
                        x_label = "Inf"
                    else:
                        x_label = str(round(math.sqrt(1/f), 3));
                    # print j, x_label
                    x_labels[j] = x_label
            else:
                x_labels = None
            """

            fig1 = plt.figure()
            xdata = [x[0] for x in data.fp_vs_reso()]
            ydata = [y[1] for y in data.fp_vs_reso()]

            ymin = 0
            if (intensity_data):
                ymin = min(ydata)

            xmax = max(xdata)
            ymax = max(ydata)
            plt.xlim(0, xmax)
            plt.ylim(ymin, ymax)

            s = 8
            # If this is F data the the y axis min should be 0.
            # The x-axis min should be 0.
            print("making plot", png_file_name)
            l = plt.scatter(xdata, ydata, s, alpha=0.2)
            plt.savefig(png_file_name)


    # except TypeError as e:
    #            print "caught TypeError:", e
