
from __future__ import print_function

import math
import numpy

import os
n_term_columns = int(os.popen('tput cols', 'r').read().split()[0])
numpy.set_printoptions(precision=4, linewidth=(n_term_columns-2))

def read_data(file_name, add_specs=False):

    data = []
    f = open(file_name)
    if f:
        l = True
        res_count = 0
        while l:
            l = f.readline();
            if len(l) > 40:
                if l[:5] == 'Sums:':
                    words = l.split()
                    # print(words)
                    try:
                        w_08 = int(words[8])
                        w_09 = float(words[9])
                        w_12 = int(words[12])
                        w_13 = float(words[13])
                        w_16 = int(words[16])
                        w_17 = float(words[17])
                        w_20 = int(words[20])
                        w_21 = float(words[21])
                        w_24 = int(words[24])
                        w_25 = float(words[25])
                        w_28 = int(words[28])
                        w_29 = float(words[29])
                        # print(w_28, w_29)
                        # t = ([w_08, w_09], [w_12, w_13], [w_16, w_17], [w_20, w_21], [w_24, w_25], [w_28, w_29])
                        # we don't want sums for waters or metals, so add a sanity
                        # check here: Also do a check on counts would be useful
                        if w_09 > 0.1: # needs optimization
                            if w_08 > 100: # needs optimization
                                t = [w_09, w_13, w_17, w_21, w_25, w_29]
                                if add_specs:
                                    # what is concat with space?
                                    spec = words[1] + ' ' + words[2] + ' ' + words[3] + ' ' + words[4] + ' ' + words[5]
                                    p = (spec,t)
                                    data.append(p)
                                else:
                                    data.append(t)

                        # next
                        res_count = res_count + 1
                    except TypeError as e:
                        print('problem', e)

    f.close()
    return data

# return a tuple of the means, the SigmaInverse and the determinant of Sigma
#
def find_stats(data_in, n_var):

    # trim the data to the first n_var variables
    data = [line[:n_var] for line in data_in ]

    sums     = [0] * n_var
    p        = [0] * n_var
    products = [p] * n_var
    means    = [0] * n_var

    n_data = len(data)
    print('n_data: ', n_data)
    for line in data:
        for idx in range(n_var):
            sums[idx] += line[idx]

    for idx in range(n_var):
        m = sums[idx]/n_data
        means[idx] = m

    for id in range(len(data)):
        for i in range(n_var):
            data[id][i] -= means[i]

    d = numpy.array(data)
    dt = d.transpose()

    pp = numpy.dot(dt, d)

    ppn = pp * (1.0/n_data)

    for i in range(n_var):
        # print('Type index ', i, ' has mean ', means[i], ' sd ', math.sqrt(ppn[i][i]))
        print('Type index {} has mean {:5.4f} sd {:5.4f} '.format(i, means[i], math.sqrt(ppn[i][i])))

    ppm = numpy.matrix(ppn)

    if False:

        print('pp')
        print(pp)
        print(pp.shape)

        print('ppn')
        print(ppn)
        print(ppn.shape)
    
        print('ppm')
        print(ppm)
        print(ppm.shape)

    det = numpy.linalg.det(ppm)
    print("determinant {}".format(det))

    ppi = numpy.linalg.inv(ppm)

    if False:
        print('ppi')
        print(ppi)
        print(ppi.shape)

        # test the dot product of the matrix and its inverse
        p_pi = numpy.dot(ppn, ppi)

        print('p_pi')
        print(p_pi)
        print(p_pi.shape)

    return (means, ppi, det)

def test_probability_model(means, SigmaInverse, Sigma_det):

    n_vars = len(means)

    means = means[:n_var]

    test_x = [0.41682, 0.22097, 0.22779, 0.0374, 0.073089 ]
    test_x = test_x[:n_var]

    P = probability_model(test_x, means, SigmaInverse, Sigma_det)

    return P

def triple_mult(x_in, InvSigma):

    x = numpy.array(x_in)
    xt = x.transpose()

    part_1 = xt.dot(InvSigma)
    r = part_1.dot(x)

    if False:
        print('triple_mult(): part_1 {}'.format(part_1))
        print('triple_mult(): x {}     '.format(x))
        print('triple_mult(): r {}     '.format(r))

    return r

def test_trip_mult():

    test_x = [0.41682, 0.2209]

    data = [[1,2], [3,4]]

    m = numpy.array(data)

    m_inv = numpy.linalg.inv(m)

    det = numpy.linalg.det(m)

    p = triple_mult(test_x, m_inv, det)

    print('p: {}'.format(p))

def probability_model(x_in, means, InvSigma, det_Sigma):

    # the size of the x_in will be chopped to match that of means.
    # 

    n_var = len(means)
    normalizing = 1.0/pow(2 * math.pi, 0.5 * n_var)

    # print("Here with x_in", x_in)
    # print("Here with means", means)

    x = x_in[:n_var]
    for i in range(n_var):
        x[i] -= means[i]

    xt_S_x = triple_mult(x, InvSigma)

    coeff = -0.5 * xt_S_x

    if False:
        print('xt_S_x {}'.format(xt_S_x))
        print('coeff {} '.format(coeff))

    p = normalizing / math.sqrt(det_Sigma) * math.exp(coeff)

    return p

def get_probabilities(our_residue_fractions, means, InvSigma, det_Sigma):

    for residue_spec,fractions in our_residue_fractions:
        x = fractions
        p = probability_model(x, means, InvSigma, det_Sigma)
        print(residue_spec, "probability-density ", p)


# data = read_data('out-1000')
data = read_data('out.orig')

n_var = 5

(means, SigmaInverse, Sigma_det) = find_stats(data, n_var)

# P = test_probability_model(means, SigmaInverse, Sigma_det)
# print("P: {}".format(P))

# test that model against "our" protein model
#
our_fractions = read_data("A-chain.tab", add_specs=True)

r = get_probabilities(our_fractions, means, SigmaInverse, Sigma_det)

