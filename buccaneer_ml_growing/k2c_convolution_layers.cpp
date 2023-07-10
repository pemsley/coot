/**
k2c_convolution_layers.c
This file is part of keras2c
Copyright 2020 Rory Conlin
Licensed under MIT License
https://github.com/f0uriest/keras2c
 */


#include <math.h>
#include <stdio.h>
#include <string.h>
#include "k2c_include.h"


/**
 * 1D (temporal) Padding.
 *
 * :param output: tensor to store padded output data.
 * :param input: tensor to pad.
 * :param fill: value to fill in padded areas.
 * :param pad: array[2] of how many rows to pad. Order is {before dim 1, after dim 1}.
 */
void k2c_pad1d(k2c_tensor* output, const k2c_tensor* input, const float fill,
               const size_t * pad) {

    const size_t in_width = input->shape[1];
    const size_t pad_top = pad[0];

    // set output array to fill value
    if (fabs(fill) < 1e-6) {
        // fill is ~zero, use memset
        memset(output->array,0,output->numel*sizeof(output->array[0]));
    }
    else {
        for(size_t i=0; i<output->numel; ++i) {
            output->array[i] = fill;
        }
    }

    // memcpy the old array in the right place
    const size_t offset = pad_top*in_width;
    memcpy(&output->array[offset],&input->array[0],
           input->numel*sizeof(input->array[0]));
}


/**
 * 2D (spatial) Padding.
 *
 * :param output: tensor to store padded output data.
 * :param input: tensor to pad.
 * :param fill: value to fill in padded areas.
 * :param pad: array[4] of how many rows/cols to pad. Order is {before dim 1, after dim 1, before dim 2, after dim 2}.
 */
void k2c_pad2d(k2c_tensor* output, const k2c_tensor* input, const float fill,
               const size_t * pad) {

    const size_t in_height = input->shape[0];
    const size_t in_width = input->shape[1];
    const size_t in_channels = input->shape[2];
    const size_t pad_top = pad[0];
    const size_t pad_left = pad[2];
    const size_t pad_right = pad[3];

    // set output array to fill value
    if (fabs(fill) < 1e-6) {
        // fill is ~zero, use memset
        memset(output->array,0,output->numel*sizeof(output->array[0]));
    }
    else {
        for(size_t i=0; i<output->numel; ++i) {
            output->array[i] = fill;
        }
    }
    // memcpy the old array in the middle
    size_t offset = in_channels*(pad_left+pad_right+in_width)*pad_top +
                    in_channels*pad_left;
    const size_t num = in_channels*in_width;
    const size_t step = num+in_channels*(pad_left+pad_right);
    for (size_t i=0; i<in_height; ++i) {
        memcpy(&output->array[offset],
               &input->array[i*num],
               num*sizeof(input->array[0]));
        offset += step;
    }
}


/**
 * 3D (spatial or spatio-temporal) Padding.
 *
 * :param output: tensor to store padded output data.
 * :param input: tensor to pad.
 * :param fill: value to fill in padded areas.
 * :param pad: array[6] of how many rows/cols to pad. Order is {before dim 1, after dim 1, before dim 2, after dim 2, before dim 3, after dim 3}.
 */
void k2c_pad3d(k2c_tensor* output, const k2c_tensor* input, const float fill,
               const size_t * pad) {

    const size_t dim1 = input->shape[0];
    const size_t dim2 = input->shape[1];
    const size_t dim3 = input->shape[2];
    const size_t outdim1 = dim1 + pad[0] + pad[1];
    const size_t outdim2 = dim2 + pad[2] + pad[3];
    const size_t outdim3 = dim3 + pad[4] + pad[5];
    const size_t in_channels = input->shape[3];

    // set output array to fill value
    if (fabs(fill) < 1e-6) {
        // fill is ~zero, use memset
        memset(output->array,0,output->numel*sizeof(output->array[0]));
    }
    else {
        for(size_t i=0; i<output->numel; ++i) {
            output->array[i] = fill;
        }
    }
    // memcpy the old array in the middle
    const size_t offset1 = in_channels*(outdim2*outdim3)*pad[0] + in_channels*outdim3*pad[2] + in_channels*pad[4];
    const size_t num = in_channels*dim3;
    const size_t outstep2 = num+in_channels*(pad[4]+pad[5]);
    const size_t outstep1 = outdim2*outdim3*in_channels;
    const size_t instep1 = dim2*dim3*in_channels;
    const size_t instep2 = dim3*in_channels;

    for (size_t i=0; i<dim1; ++i) {
        for (size_t j=0; j<dim2; ++j) {
            memcpy(&output->array[offset1+i*outstep1 + j*outstep2],
                   &input->array[i*instep1+j*instep2],
                   num*sizeof(input->array[0]));
        }
    }
}


/**
 * 1D (temporal) Convolution.
 * Assumes a "channels last" structure.
 *
 * :param output: output tensor.
 * :param input: input tensor.
 * :param kernel: kernel tensor.
 * :param bias: bias tensor.
 * :param stride: stride length of the convolution.
 * :param dilation: dilation rate to use for dilated convolution.
 * :param activation: activation function to apply to output.
 */
void k2c_conv1d(k2c_tensor* output, const k2c_tensor* input, const k2c_tensor* kernel,
                const k2c_tensor* bias, const size_t stride, const size_t dilation,
                k2c_activationType *activation) {

    memset(output->array,0,output->numel*sizeof(output->array[0]));

    const size_t out_times = output->shape[0];
    const size_t out_channels = output->shape[1];
    const size_t in_channels = input->shape[1];

    for (size_t x0=0; x0 < out_times; ++x0) {
        for (size_t z=0; z < kernel->shape[0]; ++z) {
            for (size_t q=0; q < in_channels; ++q) {
                for (size_t k=0; k < out_channels; ++k) {
                    output->array[x0*out_channels + k] +=
                        kernel->array[z*(kernel->shape[2]*kernel->shape[1]) +
                                                                            q*(kernel->shape[2]) + k]*
                        input->array[(x0*stride + dilation*z)*in_channels + q];
                }
            }
        }
    }
    k2c_bias_add(output,bias);
    activation(output->array,output->numel);
}


/**
 * 2D (spatial) Convolution.
 * Assumes a "channels last" structure.
 *
 * :param output: output tensor.
 * :param input: input tensor.
 * :param kernel: kernel tensor.
 * :param bias: bias tensor.
 * :param stride: array[2] of stride length of the convolution. Order is {stride dim 1, stride dim 2}.
 * :param dilation: array[2] dilation rate to use for dilated convolution. Order is {dilation dim 1, dilation dim 2}.
 * :param activation: activation function to apply to output.
 */
void k2c_conv2d(k2c_tensor* output, const k2c_tensor* input, const k2c_tensor* kernel,
                const k2c_tensor* bias, const size_t * stride, const size_t * dilation,
                k2c_activationType *activation) {

    memset(output->array,0,output->numel*sizeof(output->array[0]));

    const size_t out_rows = output->shape[0];
    const size_t out_cols = output->shape[1];
    const size_t out_channels = output->shape[2];
    const size_t in_channels = input->shape[2];

    for (size_t x0=0; x0 < out_rows; ++x0) {
        for (size_t x1=0; x1 < out_cols; ++x1) {
            for (size_t z0=0; z0 < kernel->shape[0]; ++z0) {
                for (size_t z1=0; z1 < kernel->shape[1]; ++z1) {
                    for (size_t q=0; q < in_channels; ++q) {
                        for (size_t k=0; k < out_channels; ++k) {
                            output->array[x0*(output->shape[2]*output->shape[1])
                                          + x1*(output->shape[2]) + k] +=
                                              kernel->array[z0*(kernel->shape[3]*kernel->shape[2]*kernel->shape[1])
                                                            + z1*(kernel->shape[3]*kernel->shape[2])
                                                            + q*(kernel->shape[3]) + k]*
                                              input->array[(x0*stride[0]
                                                            + dilation[0]*z0)*(input->shape[2]*input->shape[1])
                                                           + (x1*stride[1] + dilation[1]*z1)*(input->shape[2]) + q];
                        }
                    }
                }
            }
        }
    }
    k2c_bias_add(output,bias);
    activation(output->array,output->numel);
}


/**
 * 3D (spatial or spatio-temporal) Convolution.
 * Assumes a "channels last" structure.
 *
 * :param output: output tensor.
 * :param input: input tensor.
 * :param kernel: kernel tensor.
 * :param bias: bias tensor.
 * :param stride: array[3] of stride length of the convolution. Order is {stride dim 1, stride dim 2, stride dim 3}.
 * :param dilation: array[3] dilation rate to use for dilated convolution. Order is {dilation dim 1, dilation dim 2, dilation dim 3}.
 * :param activation: activation function to apply to output.
 */
void k2c_conv3d(k2c_tensor* output, const k2c_tensor* input, const k2c_tensor* kernel,
                const k2c_tensor* bias, const size_t * stride, const size_t * dilation,
                k2c_activationType *activation) {

    memset(output->array,0,output->numel*sizeof(output->array[0]));
    const size_t dim1 = output->shape[0];
    const size_t dim2 = output->shape[1];
    const size_t dim3 = output->shape[2];
    const size_t out_channels = output->shape[3];
    const size_t in_channels = input->shape[3];

    for (size_t x0=0; x0 < dim1; ++x0) {
        for (size_t x1=0; x1 < dim2; ++x1) {
            for (size_t x2=0; x2<dim3; ++x2) {
                for (size_t z0=0; z0 < kernel->shape[0]; ++z0) {
                    for (size_t z1=0; z1 < kernel->shape[1]; ++z1) {
                        for (size_t z2=0; z2 < kernel->shape[2]; ++z2) {
                            for (size_t q=0; q < in_channels; ++q) {
                                for (size_t k=0; k < out_channels; ++k) {
                                    output->array[x0*(output->shape[3]*output->shape[2]
                                                      *output->shape[1])
                                                  + x1*(output->shape[3]*output->shape[2])
                                                  + x2*(output->shape[3]) + k] +=
                                                      kernel->array[z0*(kernel->shape[4]*kernel->shape[3]
                                                                        *kernel->shape[2]*kernel->shape[1])
                                                                    + z1*(kernel->shape[4]*kernel->shape[3]
                                                                          *kernel->shape[2])
                                                                    + z2*(kernel->shape[4]*kernel->shape[3])
                                                                    + q*(kernel->shape[4]) + k]
                                                      *input->array[(x0*stride[0] + dilation[0]*z0)
                                                                    *(input->shape[3]*input->shape[2]*input->shape[1])
                                                                    + (x1*stride[1] + dilation[1]*z1)
                                                                    *(input->shape[3]*input->shape[2])
                                                                    + (x2*stride[2] + dilation[2]*z2)
                                                                    *(input->shape[3]) + q];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    k2c_bias_add(output,bias);
    activation(output->array,output->numel);
}


/**
 * 1D (temporal) Cropping.
 *
 * :param output: tensor to store cropped output data.
 * :param input: tensor to crop.
 * :param pad: array[2] of how many rows to crop. Order is {before dim 1, after dim 1}.
 */
void k2c_crop1d(k2c_tensor* output, const k2c_tensor* input, const size_t * crop) {

    const size_t offset = crop[0]*input->shape[1];
    memcpy(&output->array[0],&input->array[offset],
           output->numel*sizeof(output->array[0]));
}


/**
 * 2D (spatial) Cropping.
 *
 * :param output: tensor to store cropped output data.
 * :param input: tensor to crop.
 * :param pad: array[4] of how many rows/cols to crop. Order is {before dim 1, after dim 1, before dim 2, after dim 2}.
 */
void k2c_crop2d(k2c_tensor* output, const k2c_tensor* input, const size_t * crop) {

    const size_t out_height = output->shape[0];
    const size_t in_width = input->shape[1];
    const size_t in_channels = input->shape[2];
    const size_t crop_top = crop[0];
    const size_t crop_left = crop[2];
    const size_t crop_right = crop[3];

    size_t offset = in_channels*in_width*crop_top + in_channels*crop_left;
    const size_t num = in_channels*(in_width-crop_left-crop_right);
    for (size_t i=0; i<out_height; ++i) {
        memcpy(&output->array[i*num],&input->array[offset],num*sizeof(input->array[0]));
        offset += in_width*in_channels;
    }
}


/**
 * 3D (spatial or spatio-temporal) Cropping.
 *
 * :param output: tensor to store cropped output data.
 * :param input: tensor to crop.
 * :param pad: array[6] of how many rows/cols to crop. Order is {before dim 1, after dim 1, before dim 2, after dim 2, before dim 3, after dim 3}.
 */
void k2c_crop3d(k2c_tensor* output, const k2c_tensor* input, const size_t * crop) {

    const size_t dim1 = input->shape[0];
    const size_t dim2 = input->shape[1];
    const size_t dim3 = input->shape[2];
    const size_t outdim1 = dim1 - crop[0] - crop[1];
    const size_t outdim2 = dim2 - crop[2] - crop[3];
    const size_t outdim3 = dim3 - crop[4] - crop[5];
    const size_t in_channels = input->shape[3];

    const size_t offset1 = in_channels*(dim2*dim3)*crop[0] +
                           in_channels*dim3*crop[2] + in_channels*crop[4];
    const size_t num = in_channels*outdim3;
    const size_t instep2 = num+in_channels*(crop[4]+crop[5]);
    const size_t instep1 = dim2*dim3*in_channels;
    const size_t outstep1 = outdim2*outdim3*in_channels;
    const size_t outstep2 = outdim3*in_channels;

    for (size_t i=0; i<outdim1; ++i) {
        for (size_t j=0; j<outdim2; ++j) {
            memcpy(&output->array[i*outstep1 + j*outstep2],
                   &input->array[offset1+i*instep1+j*instep2],
                   num*sizeof(input->array[0]));
        }
    }
}


/**
 * 1D (temporal) Upsampling.
 * Repeats each temporal step size times along the time axis.
 *
 * :param output: output tensor.
 * :param input: input tensor.
 * :param size: Upsampling factor.
 */
void k2c_upsampling1d(k2c_tensor* output, const k2c_tensor* input, const size_t size) {

    const size_t in_height = input->shape[0];
    const size_t in_width = input->shape[1];

    for (size_t i=0; i<in_height; ++i) {
        for (size_t j=0; j<size; ++j) {
            for (size_t k=0; k<in_width; ++k) {
                output->array[(size*i+j)*in_width + k] = input->array[i*in_width+k];
            }
        }
    }
}


/**
 * 2D (spatial) Upsampling.
 * Repeats the rows and columns of the data by size[0] and size[1] respectively.
 *
 * :param output: output tensor.
 * :param input: input tensor.
 * :param size: array[2] of upsampling factors. Order is {upsampling dim 1, upsampling dim 2}.
 */
void k2c_upsampling2d(k2c_tensor* output, const k2c_tensor* input, const size_t * size) {

    const size_t out_height = output->shape[0];
    const size_t out_width = output->shape[1];
    const size_t channels = input->shape[2];

    for (size_t i=0; i<out_height; ++i) {
        for (size_t j=0; j<out_width; ++j) {
            const size_t insub[K2C_MAX_NDIM] = {i/size[0],j/size[1],0};
            const size_t outsub[K2C_MAX_NDIM] = {i,j,0};
            memcpy(&output->array[k2c_sub2idx(outsub,output->shape,output->ndim)],
                   &input->array[k2c_sub2idx(insub,input->shape,input->ndim)],
                   channels*sizeof(input->array[0]));
        }
    }
}


/**
 * 2D (spatial) Upsampling.
 * Repeats the 1st, 2nd and 3rd dimensions of the data by size[0], size[1] and size[2] respectively.
 *
 * :param output: output tensor.
 * :param input: input tensor.
 * :param size: array[3] of upsampling factors. Order is {upsampling dim 1, upsampling dim 2, upsampling dim 3}.
 */
void k2c_upsampling3d(k2c_tensor* output, const k2c_tensor* input, const size_t * size) {

    const size_t dim1 = output->shape[0];
    const size_t dim2 = output->shape[1];
    const size_t dim3 = output->shape[2];
    const size_t channels = input->shape[3];

    for (size_t i=0; i<dim1; ++i) {
        for (size_t j=0; j<dim2; ++j) {
            for (size_t k=0; k<dim3; ++k) {
                const size_t insub[K2C_MAX_NDIM] = {i/size[0],j/size[1],k/size[2],0};
                const size_t outsub[K2C_MAX_NDIM] = {i,j,k,0};
                memcpy(&output->array[k2c_sub2idx(outsub,output->shape,output->ndim)],
                       &input->array[k2c_sub2idx(insub,input->shape,input->ndim)],
                       channels*sizeof(input->array[0]));
            }
        }
    }
}
