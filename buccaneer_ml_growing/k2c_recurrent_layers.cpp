/**
k2c_recurrent_layers.c
This file is part of keras2c
Copyright 2020 Rory Conlin
Licensed under MIT License
https://github.com/f0uriest/keras2c
 */

#include <math.h>
#include <stdio.h>
#include "k2c_include.h"


/**
 * Cell for the LSTM layer.
 * "units" is the dimension of the output space
 *
 * :param state: array[2*units] recurrent state.
 * :param input: array of input data.
 * :param kernel: kernel tensor.
 * :param recurrent_kernel: recurrent kernel tensor
 * :param bias: bias tensor.
 * :param fwork: array[8*units] working storage.
 * :param recurrent_activation: activation function to apply to internal state.
 * :param output_activation: activation function to apply to output.
 */
void k2c_lstmcell(float * state, const float * input, const k2c_tensor* kernel,
                  const k2c_tensor* recurrent_kernel, const k2c_tensor* bias, float * fwork,
                  k2c_activationType *recurrent_activation,
                  k2c_activationType *output_activation) {


    const size_t units = recurrent_kernel->shape[1];
    const size_t in_width = kernel->shape[0]/4;

    float *h_tm1 = &state[0];  // previous memory state
    float *c_tm1 = &state[units];  // previous carry state
    const size_t outrows = 1;
    const float * const Wi = &kernel->array[0];
    const float * const Wf = &kernel->array[in_width*units];
    const float * const Wc = &kernel->array[2*in_width*units];
    const float * const Wo = &kernel->array[3*in_width*units];
    const float * const Ui = &recurrent_kernel->array[0];
    const float * const Uf = &recurrent_kernel->array[units*units];
    const float * const Uc = &recurrent_kernel->array[2*units*units];
    const float * const Uo = &recurrent_kernel->array[3*units*units];
    const float * const bi = &bias->array[0];
    const float * const bf = &bias->array[units];
    const float * const bc = &bias->array[2*units];
    const float * const bo = &bias->array[3*units];
    float *xi = &fwork[0];
    float *xf = &fwork[units];
    float *xc = &fwork[2*units];
    float *xo = &fwork[3*units];
    float *yi = &fwork[4*units];
    float *yf = &fwork[5*units];
    float *yc = &fwork[6*units];
    float *yo = &fwork[7*units];


    //xi = input*Wi + bi;
    k2c_affine_matmul(xi, input, Wi, bi, outrows, units, in_width);
    //xf = input*Wf + bf;
    k2c_affine_matmul(xf, input, Wf, bf, outrows, units, in_width);
    //xc = input*Wc + bc;
    k2c_affine_matmul(xc, input, Wc, bc, outrows, units, in_width);
    //xo = input*Wo + bo;
    k2c_affine_matmul(xo, input, Wo, bo, outrows, units, in_width);

    // yi = recurrent_activation(xi + h_tm1*Ui);
    k2c_affine_matmul(yi, h_tm1, Ui, xi, outrows, units, units);
    recurrent_activation(yi, units);

    // yf = recurrent_activation(xf + h_tm1*Uf);
    k2c_affine_matmul(yf, h_tm1, Uf, xf, outrows, units, units);
    recurrent_activation(yf, units);

    // yc = yf.*c_tm1 + yi.*output_activation(xc + h_tm1*Uc);
    k2c_affine_matmul(yc, h_tm1, Uc, xc, outrows, units, units);
    output_activation(yc, units);
    for (size_t i=0; i < units; ++i) {
        yc[i] = yf[i]*c_tm1[i] + yi[i]*yc[i];
    }

    // yo = recurrent_activation(xo + h_tm1*Uo);
    k2c_affine_matmul(yo, h_tm1, Uo, xo, outrows, units, units);
    recurrent_activation(yo, units);

    // h = yo.*output_activation(yc);
    // state = [h;yc];
    for (size_t i=0; i < units; ++i) {
        state[units+i] = yc[i];
    }

    output_activation(yc, units);

    for (size_t i=0; i < units; ++i) {
        state[i] = yo[i]*yc[i];
    }

}


/**
 * Long Short-Term Memory (LSTM) layer.
 * "units" is the dimension of the output space
 *
 * :param output: output tensor.
 * :param input: input tensor.
 * :param state: array[2*units] recurrent state.
 * :param kernel: kernel tensor.
 * :param recurrent_kernel: recurrent kernel tensor
 * :param bias: bias tensor.
 * :param fwork: array[8*units] working storage.
 * :param go_backwards: whether to process input sequences forwards (1) or backwards (0).
 * :param return_sequences: whether to return the last output in the output sequence (0), or the full sequence (1).
 * :param recurrent_activation: activation function to apply to internal state.
 * :param output_activation: activation function to apply to output.
 */
void k2c_lstm(k2c_tensor* output, const k2c_tensor* input, float * state,
              const k2c_tensor* kernel, const k2c_tensor* recurrent_kernel,
              const k2c_tensor* bias, float * fwork, const int go_backwards,
              const int return_sequences, k2c_activationType *recurrent_activation,
              k2c_activationType *output_activation) {


    const size_t in_height = input->shape[0];
    const size_t in_width = input->shape[1];
    const size_t units = recurrent_kernel->shape[1];
    if (go_backwards) {
        for (int i=in_height-1; i>-1; --i) {
            k2c_lstmcell(state, &input->array[i*in_width], kernel, recurrent_kernel,
                         bias, fwork, recurrent_activation, output_activation);
            if (return_sequences) {
                for (size_t j=0; j<units; ++j) {
                    output->array[(in_height-1-i)*units+j] = state[j];
                }
            }
        }
    }
    else {
        for (size_t i=0; i < in_height; ++i) {
            k2c_lstmcell(state, &input->array[i*in_width], kernel, recurrent_kernel,
                         bias, fwork, recurrent_activation, output_activation);
            if (return_sequences) {
                for (size_t j=0; j<units; ++j) {
                    output->array[i*units+j] = state[j];
                }
            }
        }
    }
    if (!return_sequences) {
        for (size_t i=0; i < units; ++i) {
            output->array[i] = state[i];
        }
    }
}


/**
 * Cell for the RNN layer.
 * "units" is the dimension of the output space
 *
 * :param state: array[units] recurrent state.
 * :param input: array of input data.
 * :param kernel: kernel tensor.
 * :param recurrent_kernel: recurrent kernel tensor
 * :param bias: bias tensor.
 * :param fwork: array[2*units] working storage.
 * :param output_activation: activation function to apply to output.
 */
void k2c_simpleRNNcell(float * state, const float * input, const k2c_tensor* kernel,
                       const k2c_tensor* recurrent_kernel, const k2c_tensor* bias,
                       float * fwork, k2c_activationType *output_activation) {

    const size_t units = recurrent_kernel->shape[1];
    const size_t in_width = kernel->shape[0];

    const size_t outrows = 1;
    float *h1 = &fwork[0];
    float *h2 = &fwork[units];
    // h1 = input*kernel+bias
    k2c_affine_matmul(h1,input,kernel->array,bias->array,outrows,units,in_width);

    // h2 = state*recurrent_kernel + h1
    k2c_affine_matmul(h2,state,recurrent_kernel->array,h1,outrows,units,units);
    output_activation(h2,units);

    for (size_t i=0; i<units; ++i) {
        state[i] = h2[i];
    }
}


/**
 * Fully-connected RNN where the output is to be fed back to input.
 * "units" is the dimension of the output space
 *
 * :param output: output tensor.
 * :param input: input tensor.
 * :param state: array[units] recurrent state.
 * :param kernel: kernel tensor.
 * :param recurrent_kernel: recurrent kernel tensor
 * :param bias: bias tensor.
 * :param fwork: array[2*units] working storage.
 * :param go_backwards: whether to process input sequences forwards (1) or backwards (0).
 * :param return_sequences: whether to return the last output in the output sequence (0), or the full sequence (1).
 * :param output_activation: activation function to apply to output.
 */
void k2c_simpleRNN(k2c_tensor* output, const k2c_tensor* input, float * state,
                   const k2c_tensor* kernel, const k2c_tensor* recurrent_kernel,
                   const k2c_tensor* bias, float * fwork, const int go_backwards,
                   const int return_sequences, k2c_activationType *output_activation) {

    const size_t in_width = input->shape[1];
    const size_t in_height = input->shape[0];
    const size_t units = recurrent_kernel->shape[1];

    if (go_backwards) {
        for (int i=in_height-1; i>-1; --i) {
            k2c_simpleRNNcell(state,&input->array[i*in_width],kernel,recurrent_kernel,bias,
                              fwork, output_activation);
            if (return_sequences) {
                for (size_t j=0; j<units; ++j) {
                    output->array[(in_height-1-i)*units+j] = state[j];
                }
            }
        }
    }
    else {
        for (size_t i=0; i<in_height; ++i) {
            k2c_simpleRNNcell(state,&input->array[i*in_width],kernel,recurrent_kernel,bias,
                              fwork, output_activation);
            if (return_sequences) {
                for (size_t j=0; j<units; ++j) {
                    output->array[i*units+j] = state[j];
                }
            }
        }
    }
    if (!return_sequences) {
        for (size_t i=0; i < units; ++i) {
            output->array[i] = state[i];
        }
    }
}


/**
 * Cell for the GRU layer.
 * "units" is the dimension of the output space
 *
 * :param state: array[units] recurrent state.
 * :param input: array of input data.
 * :param kernel: kernel tensor.
 * :param recurrent_kernel: recurrent kernel tensor
 * :param bias: bias tensor.
 * :param fwork: array[6*units] working storage.
 * :param reset_after: whether to apply the reset gate before (0) or after (1) the matrix multiplication.
 * :param recurrent_activation: activation function to apply to internal state.
 * :param output_activation: activation function to apply to output.
 */
void k2c_grucell(float * state, const float * input, const k2c_tensor* kernel,
                 const k2c_tensor* recurrent_kernel, const k2c_tensor* bias, float * fwork,
                 const int reset_after, k2c_activationType *recurrent_activation,
                 k2c_activationType *output_activation) {

    const size_t units = recurrent_kernel->shape[1];
    const size_t in_width = kernel->shape[0]/3;

    float *h_tm1 = &state[0];
    const size_t outrows = 1;
    const float * const Wz = &kernel->array[0];
    const float * const Wr = &kernel->array[in_width*units];
    const float * const Wh = &kernel->array[2*in_width*units];
    const float * const Uz = &recurrent_kernel->array[0];
    const float * const Ur = &recurrent_kernel->array[units*units];
    const float * const Uh = &recurrent_kernel->array[2*units*units];
    const float * const bz = &bias->array[0];
    const float * const br = &bias->array[units];
    const float * const bh = &bias->array[2*units];
    const float * const rbz = &bias->array[3*units];
    const float * const rbr = &bias->array[4*units];
    const float * const rbh = &bias->array[5*units];
    float *xz = &fwork[0];
    float *xr = &fwork[units];
    float *xh = &fwork[2*units];
    float *yz = &fwork[3*units];
    float *yr = &fwork[4*units];
    float *yh = &fwork[5*units];

    //     x_z = input*kernel_z + input_bias_z
    k2c_affine_matmul(xz, input, Wz, bz, outrows, units, in_width);
    //    x_r = input@kernel_r + input_bias_r
    k2c_affine_matmul(xr, input, Wr, br, outrows, units, in_width);
    //    x_h = input@kernel_h + input_bias_h
    k2c_affine_matmul(xh, input, Wh, bh, outrows, units, in_width);

    //   recurrent_z = h_tm1@recurrent_kernel_z
    k2c_affine_matmul(yz, h_tm1, Uz, rbz, outrows, units, units);
    //    recurrent_r = h_tm1@recurrent_kernel_r
    k2c_affine_matmul(yr, h_tm1, Ur, rbr, outrows, units, units);

    //    z = np.tanh(x_z + recurrent_z)
    //    r = np.tanh(x_r + recurrent_r)
    for (size_t i=0; i<units; ++i) {
        yz[i] = xz[i] + yz[i];
        yr[i] = xr[i] + yr[i];
    }
    recurrent_activation(yz, units);
    recurrent_activation(yr, units);

    //    reset gate applied after/before matrix multiplication
    if (reset_after) {
        //        recurrent_h = h_tm1*recurrent_kernel_h + recurrent_bias_h
        k2c_affine_matmul(yh, h_tm1, Uh, rbh, outrows, units, units);
        //        recurrent_h = r .* recurrent_h
        for (size_t i=0; i<units; ++i) {
            yh[i] = yr[i] * yh[i];
        }
    }
    else {
        //        recurrent_h = (r .* h_tm1)*recurrent_kernel_h
        for (size_t i=0; i<units; ++i) {
            yh[i] = yr[i]*h_tm1[i];
        }
        k2c_matmul(xz, yh, Uh, outrows, units, units); //reuse xz as new yh
        for (size_t i=0; i<units; ++i) {
            yh[i] = xz[i];
        }
    }
    //    hh = np.tanh(x_h + recurrent_h)
    for (size_t i=0; i<units; ++i) {
        xr[i] = xh[i] + yh[i];  // reuse xr = hh
    }
    output_activation(xr, units);

    //    h = z .* h_tm1 + (1 - z) .* hh
    for (size_t i=0; i<units; ++i) {
        state[i] = yz[i] * h_tm1[i] + (1.0f-yz[i])*xr[i];
    }
}


/**
 * Gated Recurrent Unit.
 * "units" is the dimension of the output space
 *
 * :param output: output tensor.
 * :param input: input tensor.
 * :param state: array[units] recurrent state.
 * :param kernel: kernel tensor.
 * :param recurrent_kernel: recurrent kernel tensor
 * :param bias: bias tensor.
 * :param fwork: array[6*units] working storage.
 * :param reset_after: whether to apply the reset gate before (0) or after (1) the matrix multiplication.
 * :param go_backwards: whether to process input sequences forwards (1) or backwards (0).
 * :param return_sequences: whether to return the last output in the output sequence (0), or the full sequence (1).
 * :param recurrent_activation: activation function to apply to internal state.
 * :param output_activation: activation function to apply to output.
 */
void k2c_gru(k2c_tensor* output, const k2c_tensor* input, float * state,
             const k2c_tensor* kernel, const k2c_tensor* recurrent_kernel,
             const k2c_tensor* bias, float * fwork, const int reset_after,
             const int go_backwards, const int return_sequences,
             k2c_activationType *recurrent_activation,
             k2c_activationType *output_activation) {


    const size_t in_width = input->shape[1];
    const size_t in_height = input->shape[0];
    const size_t units = recurrent_kernel->shape[1];

    if (go_backwards) {
        for (int i=in_height-1; i>-1; --i) {
            k2c_grucell(state, &input->array[i*in_width], kernel, recurrent_kernel, bias,
                        fwork, reset_after, recurrent_activation, output_activation);
            if (return_sequences) {
                for (size_t j=0; j<units; ++j) {
                    output->array[(in_height-1-i)*units+j] = state[j];
                }
            }
        }
    }
    else {
        for (size_t i=0; i<in_height; ++i) {
            k2c_grucell(state, &input->array[i*in_width], kernel, recurrent_kernel, bias,
                        fwork, reset_after, recurrent_activation, output_activation);
            if (return_sequences) {
                for (size_t j=0; j<units; ++j) {
                    output->array[i*units+j] = state[j];
                }
            }
        }
    }

    if (!return_sequences) {
        for (size_t i=0; i<units; ++i) {
            output->array[i] = state[i];
        }
    }
}
