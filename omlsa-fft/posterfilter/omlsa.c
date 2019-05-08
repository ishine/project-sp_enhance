#include"omlsa_def.h"
#include"omlsa_var.h"
/*
   omlsa.c 2018/11/1

  Speech enhancement with om-lsa and mcra algorithm

  paper: (1) Speech enhancement for non-stationary noise enviroments. By Israel Cohen   2001
         (2) Noise Spectrum Estimation in Adverse Environments: Improved Minima Controlled Recursive Averaging. By Israel Cohen   2003
  
*/


/*
  (1) win & fft, smoothing process

  (2) noise power estimate 

  (3) speech absent probability

  (4) gamma, eta, v, GH1
*/


omlsa_float32_t Window[FRAME_LEN] = 
{  
0.005614f, 0.005616f, 0.005624f, 0.005636f, 0.005653f, 0.005675f, 0.005702f, 0.005733f, 0.005770f, 0.005811f, 0.005858f, 0.005909f, 0.005965f, 0.006025f, 0.006091f, 0.006161f, 
0.006236f, 0.006316f, 0.006401f, 0.006491f, 0.006585f, 0.006684f, 0.006788f, 0.006896f, 0.007009f, 0.007127f, 0.007249f, 0.007376f, 0.007508f, 0.007644f, 0.007785f, 0.007931f, 
0.008080f, 0.008235f, 0.008394f, 0.008557f, 0.008725f, 0.008897f, 0.009074f, 0.009255f, 0.009440f, 0.009630f, 0.009823f, 0.010022f, 0.010224f, 0.010430f, 0.010641f, 0.010856f, 
0.011075f, 0.011297f, 0.011524f, 0.011755f, 0.011990f, 0.012229f, 0.012472f, 0.012718f, 0.012968f, 0.013223f, 0.013480f, 0.013742f, 0.014007f, 0.014276f, 0.014548f, 0.014824f, 
0.015103f, 0.015386f, 0.015672f, 0.015962f, 0.016255f, 0.016551f, 0.016850f, 0.017153f, 0.017459f, 0.017767f, 0.018079f, 0.018394f, 0.018712f, 0.019032f, 0.019356f, 0.019682f, 
0.020011f, 0.020343f, 0.020678f, 0.021015f, 0.021354f, 0.021696f, 0.022041f, 0.022388f, 0.022737f, 0.023089f, 0.023442f, 0.023798f, 0.024157f, 0.024517f, 0.024879f, 0.025243f, 
0.025609f, 0.025977f, 0.026347f, 0.026718f, 0.027092f, 0.027467f, 0.027843f, 0.028221f, 0.028600f, 0.028981f, 0.029363f, 0.029747f, 0.030131f, 0.030517f, 0.030904f, 0.031292f, 
0.031681f, 0.032071f, 0.032462f, 0.032853f, 0.033246f, 0.033639f, 0.034033f, 0.034427f, 0.034822f, 0.035217f, 0.035613f, 0.036009f, 0.036405f, 0.036802f, 0.037199f, 0.037596f, 
0.037993f, 0.038389f, 0.038786f, 0.039183f, 0.039579f, 0.039976f, 0.040372f, 0.040767f, 0.041162f, 0.041557f, 0.041951f, 0.042344f, 0.042737f, 0.043129f, 0.043520f, 0.043911f, 
0.044300f, 0.044689f, 0.045076f, 0.045463f, 0.045848f, 0.046232f, 0.046615f, 0.046996f, 0.047376f, 0.047755f, 0.048132f, 0.048508f, 0.048882f, 0.049254f, 0.049625f, 0.049994f, 
0.050361f, 0.050726f, 0.051089f, 0.051450f, 0.051809f, 0.052166f, 0.052521f, 0.052874f, 0.053224f, 0.053573f, 0.053918f, 0.054262f, 0.054602f, 0.054941f, 0.055277f, 0.055610f, 
0.055940f, 0.056268f, 0.056593f, 0.056915f, 0.057234f, 0.057550f, 0.057864f, 0.058174f, 0.058481f, 0.058785f, 0.059086f, 0.059384f, 0.059679f, 0.059970f, 0.060258f, 0.060542f, 
0.060823f, 0.061101f, 0.061375f, 0.061646f, 0.061913f, 0.062176f, 0.062436f, 0.062692f, 0.062944f, 0.063192f, 0.063437f, 0.063677f, 0.063914f, 0.064147f, 0.064376f, 0.064601f, 
0.064822f, 0.065039f, 0.065251f, 0.065460f, 0.065664f, 0.065865f, 0.066061f, 0.066252f, 0.066440f, 0.066623f, 0.066802f, 0.066976f, 0.067146f, 0.067312f, 0.067473f, 0.067630f, 
0.067782f, 0.067929f, 0.068072f, 0.068211f, 0.068345f, 0.068474f, 0.068599f, 0.068719f, 0.068835f, 0.068945f, 0.069051f, 0.069153f, 0.069249f, 0.069341f, 0.069428f, 0.069511f, 
0.069588f, 0.069661f, 0.069729f, 0.069792f, 0.069851f, 0.069904f, 0.069953f, 0.069997f, 0.070036f, 0.070070f, 0.070099f, 0.070123f, 0.070143f, 0.070158f, 0.070167f, 0.070172f, 
0.070172f, 0.070167f, 0.070158f, 0.070143f, 0.070123f, 0.070099f, 0.070070f, 0.070036f, 0.069997f, 0.069953f, 0.069904f, 0.069851f, 0.069792f, 0.069729f, 0.069661f, 0.069588f, 
0.069511f, 0.069428f, 0.069341f, 0.069249f, 0.069153f, 0.069051f, 0.068945f, 0.068835f, 0.068719f, 0.068599f, 0.068474f, 0.068345f, 0.068211f, 0.068072f, 0.067929f, 0.067782f, 
0.067630f, 0.067473f, 0.067312f, 0.067146f, 0.066976f, 0.066802f, 0.066623f, 0.066440f, 0.066252f, 0.066061f, 0.065865f, 0.065664f, 0.065460f, 0.065251f, 0.065039f, 0.064822f, 
0.064601f, 0.064376f, 0.064147f, 0.063914f, 0.063677f, 0.063437f, 0.063192f, 0.062944f, 0.062692f, 0.062436f, 0.062176f, 0.061913f, 0.061646f, 0.061375f, 0.061101f, 0.060823f, 
0.060542f, 0.060258f, 0.059970f, 0.059679f, 0.059384f, 0.059086f, 0.058785f, 0.058481f, 0.058174f, 0.057864f, 0.057550f, 0.057234f, 0.056915f, 0.056593f, 0.056268f, 0.055940f, 
0.055610f, 0.055277f, 0.054941f, 0.054602f, 0.054262f, 0.053918f, 0.053573f, 0.053224f, 0.052874f, 0.052521f, 0.052166f, 0.051809f, 0.051450f, 0.051089f, 0.050726f, 0.050361f, 
0.049994f, 0.049625f, 0.049254f, 0.048882f, 0.048508f, 0.048132f, 0.047755f, 0.047376f, 0.046996f, 0.046615f, 0.046232f, 0.045848f, 0.045463f, 0.045076f, 0.044689f, 0.044300f, 
0.043911f, 0.043520f, 0.043129f, 0.042737f, 0.042344f, 0.041951f, 0.041557f, 0.041162f, 0.040767f, 0.040372f, 0.039976f, 0.039579f, 0.039183f, 0.038786f, 0.038389f, 0.037993f, 
0.037596f, 0.037199f, 0.036802f, 0.036405f, 0.036009f, 0.035613f, 0.035217f, 0.034822f, 0.034427f, 0.034033f, 0.033639f, 0.033246f, 0.032853f, 0.032462f, 0.032071f, 0.031681f, 
0.031292f, 0.030904f, 0.030517f, 0.030131f, 0.029747f, 0.029363f, 0.028981f, 0.028600f, 0.028221f, 0.027843f, 0.027467f, 0.027092f, 0.026719f, 0.026347f, 0.025977f, 0.025609f, 
0.025243f, 0.024879f, 0.024517f, 0.024157f, 0.023798f, 0.023442f, 0.023089f, 0.022737f, 0.022388f, 0.022041f, 0.021696f, 0.021354f, 0.021015f, 0.020678f, 0.020343f, 0.020011f, 
0.019682f, 0.019356f, 0.019032f, 0.018712f, 0.018394f, 0.018079f, 0.017767f, 0.017459f, 0.017153f, 0.016850f, 0.016551f, 0.016255f, 0.015962f, 0.015672f, 0.015386f, 0.015103f, 
0.014824f, 0.014548f, 0.014276f, 0.014007f, 0.013742f, 0.013480f, 0.013223f, 0.012968f, 0.012718f, 0.012472f, 0.012229f, 0.011990f, 0.011755f, 0.011524f, 0.011298f, 0.011075f, 
0.010856f, 0.010641f, 0.010430f, 0.010224f, 0.010022f, 0.009823f, 0.009630f, 0.009440f, 0.009255f, 0.009074f, 0.008897f, 0.008725f, 0.008557f, 0.008394f, 0.008235f, 0.008080f, 
0.007931f, 0.007785f, 0.007644f, 0.007508f, 0.007376f, 0.007249f, 0.007127f, 0.007009f, 0.006896f, 0.006788f, 0.006684f, 0.006585f, 0.006491f, 0.006401f, 0.006316f, 0.006236f, 
0.006161f, 0.006091f, 0.006025f, 0.005965f, 0.005909f, 0.005858f, 0.005811f, 0.005770f, 0.005733f, 0.005702f, 0.005675f, 0.005653f, 0.005636f, 0.005624f, 0.005616f, 0.005614f, 
};

void fft(const omlsa_float32_t x[512], creal_T  y[512]);
 

extern void rdft(int n, int isgn, float *a);
extern void cdft(int n, int isgn, float *a);
int MinSeg = 0;


double winy_double[FFT_LEN*2];
void rdft(int n, int isgn, float *a);
void PostFilterInit(){
	int i;	 
	  
	//HammingWin75Overlap(FRAME_LEN);

	Cwin = 11.313709f;
	 
	for(i=0;i< FRAME_LEN21; i++) eta_2term[i] = 1; 
}
   int I_empty= 1;  omlsa_float32_t m_P_local;
   omlsa_float32_t *xi_local_ptr, *xi_global_ptr;omlsa_float32_t xi_peak_dB;
    	 	
#define REV_SUCESS 0
    
      omlsa_float32_t winy[FRAME_LEN*2];  //  A
	     
      omlsa_float32_t winy2[FRAME_LEN*2];  //  A
      int winy_fix[FRAME_LEN*2];  //  A

      extern void rfftlocal_fast(float* fftinput ,float* fftoutput );
	  extern void rifftlocal_fast(float* fftinput ,float* fftoutput );
	  extern void armfft(float *fftin );
	  extern void fftarm(float* fftin);
	  extern void offt(float *fftin );

	   creal_T ri_part[FFT_LEN];   
	   
	   creal_T ri_part2[FFT_LEN];   

	   float winyout[FFT_LEN*2];

int PostFilterProcess(short*pIn, int Inlen, short*pOut, int*Outlen){

	int i,rev ,size; omlsa_float32_t *Sf_ptr, *Conv_I_ptr;


   //  omlsa_float32_t lambda_d_global[FRAME_LEN21];	  
     omlsa_float32_t gamma[FRAME_LEN21];			 
     omlsa_float32_t v[FRAME_LEN21];				 
     creal_T Y_COMPLEX[FFT_LEN];       
	 omlsa_float32_t Y_2[FRAME_LEN21];	
	
     rev = REV_SUCESS;

// Temp
    

     	for(i=0;i< Inlen; i++){	
	    	pcm[i+FRAME_SHIFT] = (omlsa_float32_t)pIn[i] ;		 	
	     }
 
 
	 for(i=0;i<FRAME_LEN;i++){
#if 1
		winy[i] = pcm[i]*Window[i];
		winy[i+FRAME_LEN] =0 ;

		winy2[i] = pcm[i]*Window[i];
		winy2[i+FRAME_LEN] =0 ;
#else 
	    winy[ 2* i] = pcm[ i] *Window[i];
	    winy[ 2*i+1] = 0; 
#endif
	
	  } 

	  fft(winy,Y_COMPLEX);

	 //fftarm(winy  );
	// offt(winy);
	// armfft(winy);
   


	 rdft(512,1,winy);


     rfftlocal_fast(winy2,   (omlsa_float32_t*)winyout);

     ri_part[0].re =  winy[0];
	 ri_part[0].im =  0;

     ri_part[256].re =  winy[1];
	 ri_part[256].im =  0;

	 for(i=1;i<256;i++){
		 ri_part[i].re = winy[2*i];
		 ri_part[i].im = winy[2*i+1];	 
	 }
	 ////////
     ri_part2[0].re =  winyout[0];
	 ri_part2[0].im =  0;

     ri_part2[256].re =  winyout[1];
	 ri_part2[256].im =  0;

	 for(i=1;i<256;i++){
		 ri_part2[i].re = winyout[2*i];
		 ri_part2[i].im = winyout[2*i+1];	 
	 }
	  
     rdft(512,-1,winy);
	 rifftlocal_fast(winyout,winy2  );

	 for(i=1;i<512;i++){
		 winy[i] /= 256.0f; 
	 }
   
	 for(i=0;i<384;i++){
	   pcm[i ] = pcm[i+128];  
	 }
	  
		 
 
		 
		FrameCnt++;

		return rev;
}
 
  
//W0=win2(1:Mno);
//for k=Mno:Mno:M-1
//    swin2=lnshift(win2,k);
//    W0=W0+swin2(1:Mno);
//end
//W0=mean(W0)^0.5;
//win=win/W0;
//Cwin=sum(win.^2)^0.5;
//win=win/Cwin;
 

void HammingWin75Overlap(int len) {

   int i ; 
   omlsa_float32_t Win2[FRAME_LEN];  
   omlsa_float32_t W0[FRAME_LEN41]; 
   omlsa_float32_t W0_avg;//, Cwin;
 
    for(i=0;i<len;i++){
        Window[i] = 0.54f-0.46f*cosf(2*M_PI*(i)/(len-1));
		Win2[i] = Window[i]*Window[i];		 
	}
   
     W0_avg = 0;
     for(i=0;i<FRAME_LEN41;i++){
       W0[i] = Win2[i] + Win2[i+len/4] + Win2[i+len/2] + Win2[i+len*3/4];
	   W0_avg +=  W0[i] ;
     }

	 W0_avg /= FRAME_LEN41;

	 W0_avg= sqrtf(W0_avg);


     Cwin = 0;
    for(i=0;i<FRAME_LEN;i++){
		Window[i]/= W0_avg;
		Cwin+= Window[i]*Window[i];
	}

	Cwin = sqrtf(Cwin);
  
    for(i=0;i<FRAME_LEN;i++){
		Window[i]/= Cwin;	 
	}

   
}





 