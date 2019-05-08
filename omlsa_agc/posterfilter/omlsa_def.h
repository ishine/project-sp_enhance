

//#include "omlsa_typedefs.h"
  

#include<math.h>

typedef float omlsa_float32_t;
typedef double omlsa_float64_t;
#define Fs          16000
#define FRAME_LEN   512
#define FRAME_SHIFT  (FRAME_LEN*3/4)
#define FRAME_LEN21  (FRAME_LEN/2 +1 )
#define FRAME_LEN41  (FRAME_LEN/4)
#define FFT_LEN FRAME_LEN
#define  M  9 //(log10((int)FFT_LEN) / log10(2))

#define M_PI 3.1415926f

#define W 1
#define alpha_eta 0.95f
#define alpha_s   0.9f

#define  delta_s  1.67f 	//	% 2.4)  Local minimum factor
#define  Bmin     1.66f      //
#define  delta_Y2 4.6f 	//	% 2.4)  Local minimum factor
#define  delta_yt 3 
#define  alpha_d  0.85f
#define  alpha_d_long   0.99f 
#define  beta   0.7f 
#define  xi_min_dB  -10 
#define  xi_max_dB  -5 
#define  xi_p_min_dB  0
#define  xi_p_max_dB  10
#define  P_min  0.005f
#define max_local(a,b) ((a>b)? a:b)
#define min_local(a,b) ((a<b)? a:b)
#define abs_local(a)  ((a>0)?a:(-a))
#define  qmax  0.9980f

#define  eta_min_dB -18.0f	
#define  eta_min   (float)pow(10.0f,eta_min_dB/10.0f) //   10^(eta_min_dB/10)
#define  Gmin   0.1259f //sqrt(eta_min)	//   % Gain floor


#define k2_local 16// (int)(round(500.0/Fs*FRAME_LEN+1)) //% 500hz
#define k3_local 112//(int)(round(3500.0/Fs*FRAME_LEN+1))// %3500hz

void HannWinMake(omlsa_float32_t*win, int len);
void D_M_ANC_WIN_FFT(omlsa_float32_t* pX , int len, omlsa_float32_t*pReaOut,  omlsa_float32_t*pImgOut);
void HammingWin75Overlap(int len);
int IFFT_CAL(omlsa_float32_t* input_real, omlsa_float32_t*input_imag, omlsa_float32_t*fft_nextReal );
//void fft(const omlsa_float32_t x[512], creal_T_local y[512]);
//void ifft(creal_T_local y[512]);

void b_conv(const omlsa_float32_t A_data[], omlsa_float32_t C_data[], int C_size[1]);
void c_conv(const omlsa_float32_t A_data[], omlsa_float32_t C_data[], int C_size[1]);

omlsa_float32_t expint(omlsa_float64_t x);
void FftWin();