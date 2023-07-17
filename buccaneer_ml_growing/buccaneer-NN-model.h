#ifndef NNSELECT_H_
#define NNSELECT_H_
#include "k2c_tensor_include.h"
class nnmodel {
public:
	nnmodel(k2c_tensor *masking_input_input, k2c_tensor *dense_output,
			float *lstm_output_array, float *lstm_kernel_array,
			float *lstm_recurrent_kernel_array, float *lstm_bias_array,
			float *lstm_1_output_array, float *lstm_1_kernel_array,
			float *lstm_1_recurrent_kernel_array, float *lstm_1_bias_array,
			float *lstm_2_output_array, float *lstm_2_kernel_array,
			float *lstm_2_recurrent_kernel_array, float *lstm_2_bias_array,
			float *lstm_3_output_array, float *lstm_3_kernel_array,
			float *lstm_3_recurrent_kernel_array, float *lstm_3_bias_array,
			float *lstm_4_output_array, float *lstm_4_kernel_array,
			float *lstm_4_recurrent_kernel_array, float *lstm_4_bias_array,
			float *dense_kernel_array, float *dense_bias_array, size_t FrgNum);

	static void nnmodel_initialize(float **lstm_output_array,
			float **lstm_kernel_array, float **lstm_recurrent_kernel_array,
			float **lstm_bias_array, float **lstm_1_output_array,
			float **lstm_1_kernel_array, float **lstm_1_recurrent_kernel_array,
			float **lstm_1_bias_array, float **lstm_2_output_array,
			float **lstm_2_kernel_array, float **lstm_2_recurrent_kernel_array,
			float **lstm_2_bias_array, float **lstm_3_output_array,
			float **lstm_3_kernel_array, float **lstm_3_recurrent_kernel_array,
			float **lstm_3_bias_array, float **lstm_4_output_array,
			float **lstm_4_kernel_array, float **lstm_4_recurrent_kernel_array,
			float **lstm_4_bias_array, float **dense_kernel_array,
			float **dense_bias_array, const char *filename);
	static void nnmodel_terminate(float *lstm_output_array,
			float *lstm_kernel_array, float *lstm_recurrent_kernel_array,
			float *lstm_bias_array, float *lstm_1_output_array,
			float *lstm_1_kernel_array, float *lstm_1_recurrent_kernel_array,
			float *lstm_1_bias_array, float *lstm_2_output_array,
			float *lstm_2_kernel_array, float *lstm_2_recurrent_kernel_array,
			float *lstm_2_bias_array, float *lstm_3_output_array,
			float *lstm_3_kernel_array, float *lstm_3_recurrent_kernel_array,
			float *lstm_3_bias_array, float *lstm_4_output_array,
			float *lstm_4_kernel_array, float *lstm_4_recurrent_kernel_array,
			float *lstm_4_bias_array, float *dense_kernel_array,
			float *dense_bias_array);
};

#endif /* nnmodel */
