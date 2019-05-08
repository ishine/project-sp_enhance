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

 omlsa_float32_t omlsa_shared_buf[1040];
 
 int om_shared_cnt;
 
omlsa_float32_t* om_mem_alloc( int size ){
    omlsa_float32_t* rev;
    rev = &omlsa_shared_buf[om_shared_cnt];
    om_shared_cnt += size;

	 //printf("om_shared_cnt is %d\n ",  om_shared_cnt);

	return rev;
}

void om_mem_free(int size){
    om_shared_cnt -=size;   
}


void fft(const omlsa_float32_t x[512], creal_T  y[512]);
 
int MinSeg = 0;

void PostFilterInit(){
	int i;	 
	  
	//HammingWin75Overlap(FRAME_LEN);

	Cwin = 128;//11.313709f;
	 
	for(i=0;i< FRAME_LEN21; i++) eta_2term[i] = 1; 
}
   int I_empty= 1;  omlsa_float32_t m_P_local;
   omlsa_float32_t *xi_local_ptr, *xi_global_ptr;omlsa_float32_t xi_peak_dB;
    	 	

   void fft_opt_512(float* fftinput ,float* fftoutput);
   void ifft_opt_512(float* fftinput ,float* fftoutput);

#define REV_SUCESS 0
    
int PostFilterProcess(short*pIn, int Inlen, short*pOut, int*Outlen){

	int i,rev ,size; omlsa_float32_t *Sf_ptr, *Conv_I_ptr;
omlsa_float32_t *Sf;
 omlsa_float32_t *P_global ;
 omlsa_float32_t *xi_global;
  omlsa_float32_t *P_local ;
  omlsa_float32_t *xi_local ;

	   omlsa_float32_t* Conv_Y_ptr;
	   
      omlsa_float32_t gamma , eta, v;			 
  				 
     creal_T Y_COMPLEX[FRAME_LEN21];       
	 omlsa_float32_t Y_2[FRAME_LEN21];	
	
     rev = REV_SUCESS;

// Temp
    

     	for(i=0;i< Inlen; i++){	
	    	pcm[i+FRAME_SHIFT] = (omlsa_float32_t)pIn[i] ;		 	
	     }
 

    {
     //omlsa_float32_t winy[FRAME_LEN];  //  A
	 omlsa_float32_t*winy;

	 winy = om_mem_alloc(FRAME_LEN);
 
     for(i=0;i<FRAME_LEN/2;i++){
     	winy[i] = pcm[i]*Window[i];
    } 
    for(i=FRAME_LEN/2;i<FRAME_LEN;i++){
	    winy[i] = pcm[i]*Window[FRAME_LEN-i-1];
    }  
	   fft_opt_512(winy, (omlsa_float32_t*)Y_COMPLEX );


       om_mem_free(FRAME_LEN);
    }


	 Y_2[0]   = Y_COMPLEX[0].re * Y_COMPLEX[0].re ;
	 Y_2[256] = Y_COMPLEX[0].im * Y_COMPLEX[0].im ;

	for(i=1;i< FRAME_LEN21-1; i++){	
		Y_2[i] = Y_COMPLEX[i].re * Y_COMPLEX[i].re +  Y_COMPLEX[i].im * Y_COMPLEX[i].im;
	}

	/*  if FrameCnt==1 lambda_d=Y_2;  end*/
	if( FrameCnt==1){
		for(i=0; i<FRAME_LEN21; i++) 
			lambda_dav[i] = lambda_d[i] = Y_2[i] ;
	}
    
      //omlsa_float32_t Sf[FRAME_LEN21+3];
     Sf = om_mem_alloc(FRAME_LEN21+3);
	b_conv(Y_2, Sf, &size);  
	Sf_ptr = &Sf[1];      
 
	if(FrameCnt==1){
		for(i=0; i<FRAME_LEN21; i++){		  
			Sy[i] = Y_2[i];	 S[i]  = Sf_ptr[i];  St[i] = Sf_ptr[i];
		}
	}
	else{
		for(i=0; i<FRAME_LEN21; i++){	
			S[i] = alpha_s * S[i] + (1-alpha_s)*Sf_ptr[i];
		}
	}
  	om_mem_free(FRAME_LEN21+3);


     if(FrameCnt<15) { 
		 for(i=0; i<FRAME_LEN21; i++) {
			 Smin[i] = S[i]; 
			 Smact[i] = S[i];
		 }
	 }
     else{
		 for(i=0; i<FRAME_LEN21; i++)  {
			 Smin[i] = min_local(Smin[i], S[i]); 
			 Smact[i] =  min_local(Smact[i], S[i]); 
		 }
	 }
  
   
   {
				
	//omlsa_float32_t I[FRAME_LEN21];
	//omlsa_float32_t Conv_I[FRAME_LEN21+3];			 
    //omlsa_float32_t Conv_Y[FRAME_LEN21+3];	

    omlsa_float32_t*Conv_I= om_mem_alloc( FRAME_LEN21+3);  //[FRAME_LEN21+3];	//1  Conv_I   FRAME_LEN21+3
    omlsa_float32_t*I = om_mem_alloc( FRAME_LEN21);                             //2  I       2*FRAME_LEN21+3

	for(i=0; i<FRAME_LEN21; i++) {
		if(( Y_2[i]<(delta_Y2*Bmin*Smin[i])) && (S[i] < delta_s*Bmin*Smin[i]) ) {
			I[i]= 1; 
		}
		else {
			I[i]= 0;
		}
	}	

	b_conv(I, Conv_I, &size);   


     om_mem_free(FRAME_LEN21);                                 //2 I   FRAME_LEN21+3

	Conv_I_ptr = &Conv_I[1];

	{

		omlsa_float32_t*Conv_Y= om_mem_alloc( FRAME_LEN21+3); //  2  Conv_Y   2*FRAME_LEN21+6
	    //omlsa_float32_t Y2_I[FRAME_LEN21];
         omlsa_float32_t * Y2_I = om_mem_alloc( FRAME_LEN21); // 3 Y2_I    3*FRAME_LEN21+6
     

		for(i=0; i<FRAME_LEN21; i++){
			Y2_I[i] =  Y_2[i]*I[i];
		}
 
		b_conv(Y2_I, Conv_Y, &size);  

        om_mem_free(FRAME_LEN21);                               // 3 Y2_I   2*FRAME_LEN21+6

		Conv_Y_ptr = &Conv_Y[1];
		 
 
	}

	}
 
     if(FrameCnt<15){
	     for(i=0;i<FRAME_LEN21;i++){
  		     St[i] = S[i];
			 Smint[i] = St[i];
			 Smactt[i] = St[i];
		 }		  
	 }
	 else{
		 for(i=0;i<FRAME_LEN21;i++){

	       if( Conv_I_ptr[i])
	          St[i]=alpha_s*St[i]+(1-alpha_s)*Conv_Y_ptr[i]/Conv_I_ptr[i];
		   else 
			    St[i]=alpha_s*St[i]+(1-alpha_s)*St[i];


		   
            Smint[i]=min_local(Smint[i],St[i]);
            Smactt[i]=min_local(Smactt[i],St[i]);
		 }
	 }

	 om_mem_free(2*FRAME_LEN21+6);

{
   ex_xi_frame = xi_frame;
 xi_frame = 0;
         
		for(i=0;i<FRAME_LEN21;i++){ 

			 omlsa_float32_t phat , qhat, Sr_Y2,Sr_S   ;


             gamma  = Y_2[i] / max_local(lambda_d[i], 1e-10f); 	    
             eta  = alpha_eta * eta_2term[i] + (1-alpha_eta) * max_local(gamma -1,0);
             eta  = max_local(eta ,eta_min);
             v   = gamma * eta / (1+eta );   	
 
             xi[i] = beta * xi[i] + (1-beta) * eta ; 
   
              if(i>=2){
                  xi_frame +=  xi[i];
             }			 

			  Sr_Y2  = Y_2[i] / Bmin / max_local(Smint[i], 1e-10f);
              Sr_S  = S[i]/Bmin/max_local(Smint[i],1e-10f); 
  
			 if(Sr_Y2 >1 && Sr_Y2 <delta_yt && Sr_S <delta_s){
				  qhat= (delta_yt-Sr_Y2 )/(delta_yt-1);
				  phat =1 /(1+ (qhat/(1-qhat))*(1+eta )*(float)expf(-v ));
			 }
			 else if(Sr_Y2 >=delta_yt || Sr_S >=delta_s){
				  phat  = 1;
			 }
			 else{
				  phat  = 0;
			 }

		     alpha_dt  = alpha_d + (1- alpha_d) * phat ;  
             lambda_dav[i] = alpha_dt * lambda_dav[i] + (1- alpha_dt ) * Y_2[i];

			  if(FrameCnt<15){
				  lambda_d_long[i]=lambda_dav[i];
			  }
			  else{
				  alpha_dt_long =alpha_d_long +(1-alpha_d_long)*phat ;   
                  lambda_d_long[i]=alpha_dt_long *lambda_d_long[i]+(1-alpha_dt_long) *Y_2[i]; 
			  }
		 }
 

 }
     MinSeg = MinSeg +1;

	 if(MinSeg==15){
		 
		 MinSeg = 0;

		 if(FrameCnt==15){
			 for(i=0;i<FRAME_LEN21;i++){
				   S_his0[i][0]
				  =S_his0[i][1]
				  =S_his0[i][2]
				  =S_his0[i][3]
				 // =S_his0[i][4]
				//  =S_his0[i][5]
				  = S[i];

				   S_his[i][0]
				  =S_his[i][1]
				  =S_his[i][2]
				  =S_his[i][3]
				 // =S_his[i][4]
				 // =S_his[i][5]
				   = St[i];
			 }
		 }
		 else{
			  for(i=0;i<FRAME_LEN21;i++){

			   Smin[i] = min_local(S_his0[i][3], Smact[i]);
               //   Smin[i] = min_local(S_his0[i][5], Smact[i]);
  				//  Smin[i] = min_local(S_his0[i][4], Smin[i]);					  
  				//  Smin[i] = min_local(S_his0[i][3], Smin[i]);
  				  Smin[i] = min_local(S_his0[i][2], Smin[i]);
   			      Smin[i] = min_local(S_his0[i][1], Smin[i]);
                  Smin[i] = min_local(S_his0[i][0], Smin[i]);
				   
				//  S_his0[i][5] = S_his0[i][4];
  				//  S_his0[i][4] = S_his0[i][3];			 
				  S_his0[i][3] = S_his0[i][2];
  				  S_his0[i][2] = S_his0[i][1];
				  S_his0[i][1] = S_his0[i][0];
				  S_his0[i][0] =  Smact[i];
				   
                  Smin[i] = min_local(Smin[i], Smact[i]);
				  Smact[i] = S[i];
 

                    Smint[i] = min_local(S_his[i][3], Smact[i]);
                //  Smint[i] = min_local(S_his[i][5], Smact[i]);
				//  Smint[i] = min_local(S_his[i][4], Smint[i]);
				//  Smint[i] = min_local(S_his[i][3], Smint[i]);
  				  Smint[i] = min_local(S_his[i][2], Smint[i]);
   			      Smint[i] = min_local(S_his[i][1], Smint[i]);
                  Smint[i] = min_local(S_his[i][0], Smint[i]);

		  
				//  S_his[i][5] = S_his[i][4];			 
				//  S_his[i][4] = S_his[i][3];
				  S_his[i][3] = S_his[i][2];
  				  S_his[i][2] = S_his[i][1];
				  S_his[i][1] = S_his[i][0];
				  S_his[i][0] =  Smactt[i];
				  
                  Smint[i] = min_local(Smint[i], Smactt[i]);
				  Smactt[i] = St[i];
			  }

			   Smactt[0] =  Smactt[0];
		 }         
	 }

 
     for(i=0;i<FRAME_LEN21;i++){
		 lambda_d[i] = 1.4685f*lambda_dav[i];		 
	 }


     { 		
   
    // omlsa_float32_t P_local[FRAME_LEN21];    		    
	// omlsa_float32_t P_global[FRAME_LEN21]; 
	// omlsa_float32_t xi_local[FRAME_LEN21+3];              
	// omlsa_float32_t xi_global[FRAME_LEN21+32];
      P_local = om_mem_alloc(FRAME_LEN21);
      xi_local = om_mem_alloc(FRAME_LEN21+3);

 
       xi_frame /=(FRAME_LEN21-2);
 
       b_conv(xi,xi_local,&size);   xi_local_ptr = &xi_local[1];
      
 
      for(i=0;i<FRAME_LEN21;i++){
        omlsa_float32_t xi_local_dB; 


		  if(xi_local[i]>0){
			  xi_local_dB = 10*(omlsa_float32_t)log10f(xi_local_ptr[i]);
		  }
		  else{
			  xi_local_dB  = -100;
		  }

			if(xi_local_dB <=xi_min_dB){
				P_local[i] = P_min ;
			}
			else if((xi_local_dB >xi_min_dB) && (xi_local_dB <xi_max_dB)){
				P_local[i] =P_min+(xi_local_dB -xi_min_dB)/(xi_max_dB-xi_min_dB)*(1-P_min); 
			}
			else {
				P_local[i] = 1;
			}
	  }
      
	   om_mem_free(FRAME_LEN21+3);


        P_global = om_mem_alloc(FRAME_LEN21);
	    xi_global = om_mem_alloc(FRAME_LEN21+32);

        c_conv(xi,xi_global,&size ); xi_global_ptr = &xi_global[15];

	    for(i=0;i<FRAME_LEN21;i++){
             omlsa_float32_t   xi_global_dB ;


		  if(xi_global[i]>0){
			  xi_global_dB  = 10*(omlsa_float32_t)log10f(xi_global_ptr[i]);
		  }
		  else{
			  xi_global_dB  = -100;
		  }

		  if(xi_global_dB <=xi_min_dB){
		  	P_global[i] = P_min ;
		  }
		  else if(xi_global_dB >xi_min_dB && xi_global_dB <xi_max_dB){
		  	P_global[i] = P_min+(xi_global_dB  -xi_min_dB)/(xi_max_dB-xi_min_dB)*(1-P_min); 
		  }
		  else {
		  	P_global[i] = 1;
		  }
	  }
      om_mem_free(FRAME_LEN21+32);
	  
	  if(xi_frame>0){
  	    xi_frame_dB = 10*(omlsa_float32_t)log10f(xi_frame);
	  }
	  else {
	    xi_frame_dB = -100;
	  }
  
     
        m_P_local = 0;
        for(i=2;i<k2_local+k3_local-3;i++){          
		   m_P_local += P_local[i];
		}
		  
		m_P_local =  m_P_local/(k2_local+k3_local-3-2);
 
        if(m_P_local<0.25f){
			for(i=k2_local;i<k3_local;i++){
				P_local[i]=P_min;
			}
		}
 

         if ((m_P_local<0.5f) && (FrameCnt>120)){
             for(i=8;i<FRAME_LEN21-8;i++){			 
	 			 if(  lambda_d_long[i] > 2.5f*(lambda_d_long[i] + lambda_d_long[i-2])  ){
					    P_local[i+6] = P_local[i+7]= P_local[i+8] = P_min;
				 }
			 }
		 }
   
       if( xi_frame_dB<=xi_min_dB){
            P_frame=P_min; 
	   }   
	   else if(xi_frame >= ex_xi_frame) {
		    xi_peak_dB=min_local(max_local(xi_frame_dB,xi_p_min_dB),xi_p_max_dB);  
		    P_frame=1;    
	   }
	   else if(xi_frame_dB>=xi_peak_dB+xi_max_dB ){
		    P_frame=1;
	   }
	   else if(xi_frame_dB<=xi_peak_dB+xi_min_dB ){
            P_frame=P_min;     
	   }
	   else{
            P_frame=P_min+(xi_frame_dB-xi_min_dB-xi_peak_dB)/(xi_max_dB-xi_min_dB)*(1-P_min);
	   }
  
	

        for(i=0;i<FRAME_LEN21;i++){	
			omlsa_float32_t q, PH1, GH0,GH1,G;
 	        omlsa_float32_t minval,lambda_d_global;	
           q = 1 - P_frame*P_global[i]*P_local[i];
           q = min_local(qmax,q);

           if((i<=2)||(i>=FRAME_LEN21-3)){
                lambda_d_global =  lambda_d[i];
			}
  
           	minval = min_local(lambda_d[i] , lambda_d[i-3] );
           lambda_d_global  = min_local(minval , lambda_d[i+3] );		   
 
           gamma  = Y_2[i] / max_local(lambda_d[i], 1e-10f);
           eta = alpha_eta*eta_2term[i] + (1-alpha_eta)*max_local(gamma-1,0);
           eta = max_local(eta,eta_min);
           v   = gamma*eta/(1+eta);    

           if(q<0.9f){
               //    PH1  = 1/( 1 +(q /(1 - q)) * (1+eta) * (omlsa_float32_t)expf(-v )  );
			         PH1  = (1 - q)/(1 - q + q * (1+eta) * (omlsa_float32_t)expf(-v ));
		   }
		   else{
			   PH1  = 0;
		   }        

           if(v>5){
  	          GH1= eta/(1+eta);
           }
           else if(v<=5 && v>0){
                GH1= (eta/(eta+1))*(omlsa_float32_t)expf(0.5*expint(v )); // need
           }	
           else {
                GH1= 1;
           }

		    GH1*= 0.5;

		      Sy[i]=0.8f*Sy[i]+0.2f*Y_2[i];   
 
			  GH0= Gmin * (omlsa_float32_t)sqrt(lambda_d_global/(Sy[i]+ 1e-10f));         	
			  G= powf(GH1 ,PH1 ) * (omlsa_float32_t)powf(GH0 ,(1-PH1 ));  

			 // G*=0.5;

              eta_2term[i]=GH1 *GH1 *gamma;



			if(i<3 ||i==(FRAME_LEN21-1)){
			   Y_COMPLEX[i].re = 0;
  			   Y_COMPLEX[i].im = 0;
			}
			else {
 			   Y_COMPLEX[i].re = Y_COMPLEX[i].re *G;
			   Y_COMPLEX[i].im = Y_COMPLEX[i].im *G;			
			}
		}

         om_mem_free(FRAME_LEN21*2); 

        }

         Y_COMPLEX[0].re = 0.0f; 
        Y_COMPLEX[0].im = 0.0f;
        Y_COMPLEX[256].re = 0.0f; 
        Y_COMPLEX[256].im = 0.0f;
		 
        {
		    omlsa_float32_t x_g_c[FFT_LEN];
	 

			ifft_opt_512((omlsa_float32_t*)Y_COMPLEX,(omlsa_float32_t*)x_g_c );
  
		 
		    for(i=0;i<FRAME_LEN/2;i++){	 
 	             x_g_c[i] =  x_g_c[i]  *Window[i]* Cwin ;
            }
            for(i=FRAME_LEN/2;i<FRAME_LEN;i++){	 
 	             x_g_c[i] =  x_g_c[i]  *Window[FRAME_LEN-i-1]* Cwin  ;
            }
   
		    for(i=0;i<FRAME_LEN;i++){ 

				float temp;

				temp = out_buf[i] + x_g_c[i] ;

				out_buf[i] =(short)((temp>32767)? 32767: ((temp<-32768)? -32768:temp));

				if(i<FRAME_LEN/4)
				    pOut[i] =  (short)out_buf[i];
		    }
        }
		  	
	 
	  
		for(i=0;i<FRAME_LEN;i++){

			 if(i<FRAME_SHIFT){
			    out_buf[i] = out_buf[i+FRAME_LEN/4];
			    pcm[i ] = pcm[i +	FRAME_LEN41]; //   [[...] . ]  <- [.[...]]
			 }
			else {
			    out_buf[i] = 0;
			}
		}
 
		 
		FrameCnt++;

		return rev;
}
 
   




 