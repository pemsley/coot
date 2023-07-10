/**
k2c_merge_layers.c
This file is part of keras2c
Copyright 2020 Rory Conlin
Licensed under MIT License
https://github.com/f0uriest/keras2c
 */

#include <string.h>
#include <stdarg.h>
#include "k2c_include.h"

/**
 * Element-wise sum of several tensors.
 *
 * :param output: output tensor.
 * :param num_tensors: number of tensors being summed.
 * :param ...: variadic. Tensors to be summed.
 */
void k2c_add(k2c_tensor* output, const size_t num_tensors,...) {

    va_list args;
    const k2c_tensor *arrptr;
    va_start (args, num_tensors);
    memset(output->array, 0, output->numel*sizeof(output->array[0]));

    for (size_t i = 0; i < num_tensors; ++i) {
        arrptr = va_arg(args, k2c_tensor*);
        for (size_t j=0; j<output->numel; ++j) {
            output->array[j] += arrptr->array[j];
        }
    }
    va_end (args);
}


/**
 * Element-wise difference of two tensors.
 *
 * :param output: output tensor.
 * :param num_tensors: number of tensors being summed. Not used but kept for a consistent API with other merge layers.
 * :param tensor1: first input tensor.
 * :param tensor2: second input tensor.
 */
void k2c_subtract(k2c_tensor* output, const size_t num_tensors,
                  const k2c_tensor* tensor1, const k2c_tensor* tensor2) {

    for (size_t i=0; i<output->numel; ++i) {
        output->array[i] = tensor1->array[i]-
                           tensor2->array[i];
    }
}


/**
 * Element-wise product of several tensors.
 *
 * :param output: output tensor.
 * :param num_tensors: number of tensors being multiplied.
 * :param ...: variadic. Tensors to be multiplied.
 */
void k2c_multiply(k2c_tensor* output, const size_t num_tensors,...) {

    va_list args;
    const k2c_tensor *arrptr;
    va_start (args, num_tensors);

    for (size_t i=0; i<output->numel; ++i) {
        output->array[i] = 1.0f;
    }

    for (size_t i = 0; i < num_tensors; ++i) {
        arrptr = va_arg(args, k2c_tensor*);
        for (size_t j=0; j<output->numel; ++j) {
            output->array[j] *= arrptr->array[j];
        }
    }
    va_end (args);
}


/**
 * Element-wise average of several tensors.
 *
 * :param output: output tensor.
 * :param num_tensors: number of tensors being averaged.
 * :param ...: variadic. Tensors to be averaged.
 */
void k2c_average(k2c_tensor* output, const size_t num_tensors,...) {

    va_list args;
    const k2c_tensor *arrptr;
    const float num_tensors_inv = 1.0f/num_tensors;

    va_start (args, num_tensors);
    memset(output->array, 0, output->numel*sizeof(output->array[0]));
    for (size_t i = 0; i < num_tensors; ++i) {
        arrptr = va_arg(args, k2c_tensor*);
        for (size_t j=0; j<output->numel; ++j) {
            output->array[j] += arrptr->array[j]*num_tensors_inv;
        }
    }
    va_end (args);
}


/**
 * Element-wise maximum of several tensors.
 *
 * :param output: output tensor.
 * :param num_tensors: number of tensors over which to take max.
 * :param ...: variadic. Tensors to take the max of.
 */
void k2c_max(k2c_tensor* output, const size_t num_tensors,...) {

    va_list args;
    const k2c_tensor *arrptr;
    va_start (args, num_tensors);
    arrptr = va_arg(args, k2c_tensor*);

    for (size_t i=0; i<output->numel; ++i) {
        output->array[i] = arrptr->array[i];
    }

    for (size_t i = 0; i < num_tensors-1; ++i) {
        arrptr = va_arg(args, k2c_tensor*);
        for (size_t j=0; j<output->numel; ++j) {
            if (output->array[j]<arrptr->array[j]) {
                output->array[j] = arrptr->array[j];
            }
        }
    }
    va_end (args);
}


/**
 * Element-wise minimum of several tensors.
 *
 * :param output: output tensor.
 * :param num_tensors: number of tensors over which to take min.
 * :param ...: variadic. Tensors to take the min of.
 */
void k2c_min(k2c_tensor* output, const size_t num_tensors,...) {

    va_list args;
    const k2c_tensor *arrptr;
    va_start (args, num_tensors);
    arrptr = va_arg(args, k2c_tensor*);

    for (size_t i=0; i<output->numel; ++i) {
        output->array[i] = arrptr->array[i];
    }

    for (size_t i = 0; i < num_tensors-1; ++i) {
        arrptr = va_arg(args, k2c_tensor*);
        for (size_t j=0; j<output->numel; ++j) {
            if (output->array[j]>arrptr->array[j]) {
                output->array[j] = arrptr->array[j];
            }
        }
    }
    va_end (args);
}


/**
 * Concatenation of several tensors.
 *
 * :param output: output tensor.
 * :param axis: axis along which to concatenate.
 * :param num_tensors: number of tensors being concatenated.
 * :param ...: variadic. Tensors to concatenate.
 */
void k2c_concatenate(k2c_tensor* output, const size_t axis, const size_t num_tensors,...) {

    va_list args;
    const k2c_tensor* arrptr;
    size_t  offset = 0;
    size_t outidx;
    size_t insub[K2C_MAX_NDIM], outsub[K2C_MAX_NDIM];
    va_start (args, num_tensors);

    for (size_t i=0; i<num_tensors; ++i) {
        arrptr = va_arg(args, k2c_tensor*);
        for (size_t j=0; j<arrptr->numel; ++j) {
            k2c_idx2sub(j,insub,arrptr->shape,arrptr->ndim);
            memcpy(outsub,insub,K2C_MAX_NDIM*sizeof(size_t));
            outsub[axis] += offset;
            outidx = k2c_sub2idx(outsub,output->shape, output->ndim);
            output->array[outidx] = arrptr->array[j];
        }
        offset += arrptr->shape[axis];
    }
    va_end (args);
}
