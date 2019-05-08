#include"omlsa_def.h"
int FrameCnt = 1;

double pcm[FRAME_LEN];
double Y_real[FRAME_LEN ];
double Y_imag[FRAME_LEN ];
double winy[FRAME_LEN];

double Y_2[FRAME_LEN21];
double S_win[3];
double Sf[FRAME_LEN21+3];
double S[FRAME_LEN21];
double Sy[FRAME_LEN21];
double Sft[FRAME_LEN21];
double St[FRAME_LEN21];
double Smin[FRAME_LEN21];
double Smint[FRAME_LEN21];
double Smact[FRAME_LEN21];
double Smactt[FRAME_LEN21];
double I[FRAME_LEN21];
double Conv_I[FRAME_LEN21+3];
double Conv_Y[FRAME_LEN21+3];

double phat[FRAME_LEN21];
double qhat[FRAME_LEN21];
double Sr_Y2[FRAME_LEN21];
double Sr_S[FRAME_LEN21];

double alpha_dt[FRAME_LEN21];
double lambda_d[FRAME_LEN21];
double lambda_dav[FRAME_LEN21];
double lambda_dt[FRAME_LEN21];
double lambda_d_long[FRAME_LEN21];
double lambda_d_global[FRAME_LEN21];

double gamma[FRAME_LEN21];
double eta[FRAME_LEN21];
double v[FRAME_LEN21];
double eta_2term[FRAME_LEN21];

double S_his0[FRAME_LEN21][8];
double S_his[FRAME_LEN21][8];
double xi[FRAME_LEN21];
double xi_local[FRAME_LEN21+3], xi_local_dB[FRAME_LEN21],  P_local[FRAME_LEN21];
double xi_global[FRAME_LEN21+32],xi_global_dB[FRAME_LEN21],  P_global[FRAME_LEN21]; 
double xi_frame, ex_xi_frame, xi_frame_dB,P_frame;
double out_buf[FRAME_LEN];

// win
double Window[FRAME_LEN];
double Win2[FRAME_LEN];
double Cwin= 0;

typedef  struct {
   double re;
   double im;
}creal_T;
  creal_T Y_COMPLEX[512];
    creal_T X_COMPLEX[512];
creal_T x_g_c[FFT_LEN];
void ifft(const creal_T x[512], creal_T y[512]);
 
double alpha_dt_long[FRAME_LEN21];
double q[FRAME_LEN21];
double PH1[FRAME_LEN21];
double GH1[FRAME_LEN21];
double GH0[FRAME_LEN21];
double G[FRAME_LEN21];
double X_G_real[FRAME_LEN];
double X_G_imag[FRAME_LEN];
double x_ifft[FFT_LEN];

extern double HannWin[FRAME_LEN];