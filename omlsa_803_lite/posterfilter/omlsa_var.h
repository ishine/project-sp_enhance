#include"omlsa_def.h"
int FrameCnt = 1;


typedef  struct {
   omlsa_float32_t re;
   omlsa_float32_t im;
}creal_T;

 short out_buf[FRAME_LEN];         //A  4k   

short pcm[FRAME_LEN];             //A  

omlsa_float32_t Sf[FRAME_LEN21+3];			//B
omlsa_float32_t S[FRAME_LEN21];				//B
omlsa_float32_t Smin[FRAME_LEN21];			//B
omlsa_float32_t Smact[FRAME_LEN21];			//B
omlsa_float32_t St[FRAME_LEN21];			//B
omlsa_float32_t Smint[FRAME_LEN21];			//B
omlsa_float32_t Smactt[FRAME_LEN21];		//B
omlsa_float32_t S_his0[FRAME_LEN21][4];		//8*B
omlsa_float32_t S_his[FRAME_LEN21][4];		//8*B 


omlsa_float32_t Sy[FRAME_LEN21];			//B

omlsa_float32_t lambda_dav[FRAME_LEN21];	//B
omlsa_float32_t lambda_dt[FRAME_LEN21];		//B
omlsa_float32_t lambda_d_long[FRAME_LEN21];	//B
omlsa_float32_t lambda_d[FRAME_LEN21];	    //B   
omlsa_float32_t eta_2term[FRAME_LEN21];		//B    21k									    
omlsa_float32_t xi[FRAME_LEN21];			//B    38k

 
omlsa_float32_t xi_frame, ex_xi_frame, xi_frame_dB,P_frame;
omlsa_float32_t alpha_dt ;		//B
omlsa_float32_t alpha_dt_long ; //B

// win
omlsa_float32_t Window[FRAME_LEN];                 //A C
omlsa_float32_t Cwin= 0;
  





void ifft(const creal_T x[512], creal_T y[512]);	