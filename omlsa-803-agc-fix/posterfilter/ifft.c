/*
 * File: ifft.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 09-Nov-2018 14:14:28
 */

/* Include Files */
#include "ifft.h"
#include "omlsa_tab.h"
#include "omlsa_tab.h"
/* Function Definitions */

/*
 * Arguments    : const creal_T x[512]
 *                creal_T y[512]
 * Return Type  : void
 */


typedef unsigned char boolean_T;
void ifft(const creal_T x[512], creal_T y[512])
{
 

  int ix;
  int ju;
  int iy;
  int i;
  boolean_T tst;
  omlsa_float32_t temp_re;
  omlsa_float32_t temp_im;
  int iheight;
  int istart;
  int j;
  omlsa_float32_t twid_re;

  

  omlsa_float32_t twid_im;

  int ihi;
  ix = 0;
  ju = 0;
  iy = 0;
  for (i = 0; i < 511; i++) {
    y[iy] = x[ix];
    iy = 512;
    tst = true;
    while (tst) {
      iy >>= 1;
      ju ^= iy;
      tst = ((ju & iy) == 0);
    }

    iy = ju;
    ix++;
  }

  y[iy] = x[ix];
  for (i = 0; i <= 511; i += 2) {
    temp_re = y[i + 1].re;
    temp_im = y[i + 1].im;
    y[i + 1].re = y[i].re - y[i + 1].re;
    y[i + 1].im = y[i].im - y[i + 1].im;
    y[i].re += temp_re;
    y[i].im += temp_im;
  }

  iy = 2;
  ix = 4;
  ju = 128;
  iheight = 509;
  while (ju > 0) {
    for (i = 0; i < iheight; i += ix) {
      temp_re = y[i + iy].re;
      temp_im = y[i + iy].im;
      y[i + iy].re = y[i].re - temp_re;
      y[i + iy].im = y[i].im - temp_im;
      y[i].re += temp_re;
      y[i].im += temp_im;
    }

    istart = 1;
    for (j = ju; j < 256; j += ju) {
      twid_re = dv3[j];
      twid_im = dv4[j];
      i = istart;
      ihi = istart + iheight;
      while (i < ihi) {
        temp_re = twid_re * y[i + iy].re - twid_im * y[i + iy].im;
        temp_im = twid_re * y[i + iy].im + twid_im * y[i + iy].re;
        y[i + iy].re = y[i].re - temp_re;
        y[i + iy].im = y[i].im - temp_im;
        y[i].re += temp_re;
        y[i].im += temp_im;
        i += ix;
      }

      istart++;
    }

    //ju /= 2;
    ju >>=1;
    
    iy = ix;
    ix <<= 1;
    iheight -= iy;
  }

  for (iy = 0; iy < 512; iy++) {
    y[iy].re *= 0.001953125f;
    y[iy].im *= 0.001953125f;
  }
}

/*
 * File trailer for ifft.c
 *
 * [EOF]
 */
