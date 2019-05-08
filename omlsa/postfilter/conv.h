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
extern void b_conv(const double A_data[], double C_data[], int C_size[1]);
extern void c_conv(const double A_data[], double C_data[], int C_size[1]);
extern void conv(const double A[3], const double B[257], double C[259]);

#endif

/*
 * File trailer for conv.h
 *
 * [EOF]
 */
