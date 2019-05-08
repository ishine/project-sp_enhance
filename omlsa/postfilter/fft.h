/*
 * File: fft.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 09-Nov-2018 14:14:28
 */

#ifndef FFT_H
#define FFT_H

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
 typedef unsigned char boolean_T;
typedef  struct {
   double re;
   double im;
}creal_T;

#define true 1
#define false 0

/* Function Declarations */
extern void fft(const double x[512], creal_T y[512]);

#endif

/*
 * File trailer for fft.h
 *
 * [EOF]
 */
