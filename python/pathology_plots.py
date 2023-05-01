
import os
import coot
import math

import matplotlib
matplotlib.use('agg')
matplotlib.rc('font', family='Arial')
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
            print("making plot", png_file_name)
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
                print("caught TypeError:", e)

# these can be I, SIGI also
#
def plots(mtz, fp, sigfp):

    # meta data
    intensity_data = False
    if fp[0] == 'I':
        intensity_data = True
        
    data = coot.pathology_data(mtz, fp, sigfp)
    print('len data is', len(data))
    print('data[0] is', data[0])

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

    for icount in range(len(labels)):

        label_pair = labels[icount];
        
        ii = 0 # label indexing hack
        l_0 = label_pair[0].replace('/', '_')
        l_1 = label_pair[1].replace('/', '_')
        png_file_name = prefix + "-" + l_0 + "-vs-" + l_1 + ".png"

        fig1 = plt.figure()
        ax = plt.gca()
        plt.title(label_pair[1] + " vs. " + label_pair[0])
        ax.set_xlabel(label_pair[0])
        ax.set_ylabel(label_pair[1])
        xdata = [x[0] for x in data[icount+1]]
        ydata = [y[1] for y in data[icount+1]]

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
        plt.close()
        
