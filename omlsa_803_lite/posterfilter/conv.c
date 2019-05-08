/*
 * File: conv.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 24-Oct-2018 17:37:29
 */

/* Include Files */
#include"omlsa_def.h"
#include "conv.h"

/* Function Definitions */

/*
 * Arguments    : const omlsa_float32_t A_data[]
 *                omlsa_float32_t C_data[]
 *                int C_size[1]
 * Return Type  : void
 */
void b_conv(const omlsa_float32_t A_data[], omlsa_float32_t C_data[], int C_size[1])
{
  int jC;
  int jA2;
  omlsa_float32_t s;
  int k;
  static const omlsa_float32_t B[3] = { 0.24999999999999997f, 0.5f, 0.24999999999999997f };

  C_size[0] = 259;
  for (jC = 0; jC < 259; jC++) {
    if (257 < jC + 1) {
      jA2 = 256;
    } else {
      jA2 = jC;
    }

    s = 0.0;
    if (3 < jC + 2) {
      k = jC - 2;
    } else {
      k = 0;
    }

    while (k + 1 <= jA2 + 1) {
      s += A_data[k] * B[jC - k];
      k++;
    }

    C_data[jC] = s;
  }
}

/*
 * Arguments    : const omlsa_float32_t A_data[]
 *                omlsa_float32_t C_data[]
 *                int C_size[1]
 * Return Type  : void
 */
void c_conv(const omlsa_float32_t A_data[], omlsa_float32_t C_data[], int C_size[1])
{
  int jC;
  int jA2;
  omlsa_float32_t s;
  int k;
  static const omlsa_float32_t B[31] = { 0.000600459987399049f, 0.0023787646090222894f,
    0.0052665746155454614f, 0.0091529130879203884f, 0.013888430218137428f,
    0.019291142738590943f, 0.025153427436995991f, 0.031249999999999997f,
    0.037346572563004006f, 0.04320885726140905f, 0.048611569781862561f,
    0.053347086912079608f, 0.057233425384454542f, 0.060121235390977711f,
    0.061899540012600951f, 0.0625f, 0.061899540012600951f, 0.060121235390977711f,
    0.057233425384454542f, 0.053347086912079608f, 0.048611569781862561f,
    0.04320885726140905f, 0.037346572563004006f, 0.031249999999999997f,
    0.025153427436995991f, 0.019291142738590943f, 0.013888430218137428f,
    0.0091529130879203884f, 0.0052665746155454614f, 0.0023787646090222894f,
    0.000600459987399049f };

  C_size[0] = 287;
  for (jC = 0; jC < 287; jC++) {
    if (257 < jC + 1) {
      jA2 = 256;
    } else {
      jA2 = jC;
    }

    s = 0.0;
    if (31 < jC + 2) {
      k = jC - 30;
    } else {
      k = 0;
    }

    while (k + 1 <= jA2 + 1) {
      s += A_data[k] * B[jC - k];
      k++;
    }

    C_data[jC] = s;
  }
}

/*
 * Arguments    : const omlsa_float32_t A[3]
 *                const omlsa_float32_t B[257]
 *                omlsa_float32_t C[259]
 * Return Type  : void
 */
void conv(const omlsa_float32_t A[3], const omlsa_float32_t B[257], omlsa_float32_t C[259])
{
  int jC;
  int jA2;
  omlsa_float32_t s;
  int k;
  for (jC = 0; jC < 259; jC++) {
    if (3 < jC + 1) {
      jA2 = 2;
    } else {
      jA2 = jC;
    }

    s = 0.0;
    if (257 < jC + 2) {
      k = jC - 256;
    } else {
      k = 0;
    }

    while (k + 1 <= jA2 + 1) {
      s += A[k] * B[jC - k];
      k++;
    }

    C[jC] = s;
  }
}

/*
 * File trailer for conv.c
 *
 * [EOF]
 */
