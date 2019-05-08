#include"omlsa_def.h"
int FrameCnt = 1;


typedef  struct {
   omlsa_float32_t re;
   omlsa_float32_t im;
}creal_T;
//all   30832
 short out_buf[FRAME_LEN];                   //1024    
//omlsa_float32_t pcm[FRAME_LEN];              
short  pcm[FRAME_LEN];                       // 1024

omlsa_float32_t S[FRAME_LEN21];				 //1028
omlsa_float32_t Smin[FRAME_LEN21];			 //1028
omlsa_float32_t Smact[FRAME_LEN21];		     //1028
omlsa_float32_t St[FRAME_LEN21];			 //1028
omlsa_float32_t Smint[FRAME_LEN21];			 //1028
omlsa_float32_t Smactt[FRAME_LEN21];		 //1028
omlsa_float32_t Sy[FRAME_LEN21];			 //1028
omlsa_float32_t lambda_dav[FRAME_LEN21];	 //1028
omlsa_float32_t lambda_d_long[FRAME_LEN21];	 //1028
omlsa_float32_t lambda_d[FRAME_LEN21];	     //1028
omlsa_float32_t eta_2term[FRAME_LEN21];		 //1028					  
omlsa_float32_t S_his0[FRAME_LEN21][8];		 //1028*8
omlsa_float32_t S_his[FRAME_LEN21][8];		 //1028*8
omlsa_float32_t xi[FRAME_LEN21];			 //1028     
 
omlsa_float32_t xi_frame, ex_xi_frame, xi_frame_dB,P_frame;


// win
 
omlsa_float32_t Cwin= 0;
   
void ifft(const creal_T x[512], creal_T y[512]);	