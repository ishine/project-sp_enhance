/*
 * File: conv.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 24-Oct-2018 17:37:29
 */

#ifndef CONV_H
#define CONV_H

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
 

/* Function Declarations */
extern void b_conv(const omlsa_float32_t A_data[], omlsa_float32_t C_data[], int C_size[1]);
extern void c_conv(const omlsa_float32_t A_data[], omlsa_float32_t C_data[], int C_size[1]);
extern void conv(const omlsa_float32_t A[3], const omlsa_float32_t B[257], omlsa_float32_t C[259]);

#endif

/*
 * File trailer for conv.h
 *
 * [EOF]
 */
