

//#include "omlsa_typedefs.h"
  

#include<math.h>


#define Fs          16000
#define FRAME_LEN   512
#define FRAME_SHIFT  (FRAME_LEN*3/4)
#define FRAME_LEN21  (FRAME_LEN/2 +1 )
#define FRAME_LEN41  (FRAME_LEN/4)
#define FFT_LEN FRAME_LEN
#define  M  9 //(log10((int)FFT_LEN) / log10(2))

#define M_PI 3.1415926

#define W 1
#define alpha_eta 0.95
#define alpha_s   0.9

#define  delta_s  1.67 	//	% 2.4)  Local minimum factor
#define  Bmin     1.66      //
#define  delta_Y2 4.6 	//	% 2.4)  Local minimum factor
#define  delta_yt 3 
#define  alpha_d  0.85
#define  alpha_d_long   0.99 
#define  beta   0.7 
#define  xi_min_dB  -10 
#define  xi_max_dB  -5 
#define  xi_p_min_dB  0
#define  xi_p_max_dB  10
#define  P_min  0.005
#define max(a,b) ((a>b)? a:b)
#define min(a,b) ((a<b)? a:b)
#define  qmax  0.9980

#define  eta_min_dB -18.0	
#define  eta_min   pow(10.0,eta_min_dB/10.0) //   10^(eta_min_dB/10)
#define  Gmin   0.1259 //sqrt(eta_min)	//   % Gain floor


#define k2_local 16// (int)(round(500.0/Fs*FRAME_LEN+1)) //% 500hz
#define k3_local 112//(int)(round(3500.0/Fs*FRAME_LEN+1))// %3500hz

void HannWinMake(double*win, int len);
void D_M_ANC_WIN_FFT(double* pX , int len, double*pReaOut,  double*pImgOut);
void HammingWin75Overlap(int len);
int IFFT_CAL(double* input_real, double*input_imag, double*fft_nextReal );
//void fft(const double x[512], creal_T_local y[512]);
//void ifft(creal_T_local y[512]);

void b_conv(const double A_data[], double C_data[], int C_size[1]);
void c_conv(const double A_data[], double C_data[], int C_size[1]);

double expint(double x);
void FftWin();