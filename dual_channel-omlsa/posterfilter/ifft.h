/*
 * File: ifft.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 09-Nov-2018 14:14:28
 */

#ifndef IFFT_H
#define IFFT_H

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include"omlsa_def.h"
typedef  struct {
   omlsa_float32_t re;
   omlsa_float32_t im;
}creal_T;
 typedef unsigned char boolean_T;
#define true 1
#define false 0
/* Function Declarations */
extern void ifft(const creal_T x[512], creal_T y[512]);

#endif

/*
 * File trailer for ifft.h
 *
 * [EOF]
 */
