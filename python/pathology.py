
import os
import coot
import matplotlib

def pathology_plots(mtz, fp, sigfp):
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

            

    except TypeError as e:
                print "caught TypeError:", e

