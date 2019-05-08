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

   

void fft(const omlsa_float32_t x[512], creal_T  y[512]);
 
int MinSeg = 0;

void PostFilterInit(){
	int i;	 
	  
	HammingWin75Overlap(FRAME_LEN);
	 
	for(i=0;i< FRAME_LEN21; i++) eta_2term[i] = 1; 
}
   int I_empty= 1;  omlsa_float32_t m_P_local;
   omlsa_float32_t *xi_local_ptr, *xi_global_ptr;omlsa_float32_t xi_peak_dB;
    	 	
#define REV_SUCESS 0


   //Ram struct to 
  

int PostFilterProcess(short*pIn, int Inlen, short*pOut, int*Outlen){

	int i,rev ; int size; omlsa_float32_t *Sf_ptr, *Conv_I_ptr;
     omlsa_float32_t lambda_d_global[FRAME_LEN21];	 //  B 
     omlsa_float32_t gamma[FRAME_LEN21];				 //  B  
     omlsa_float32_t v[FRAME_LEN21];					 //  B  T  
	 omlsa_float32_t Y_2[FRAME_LEN21];	
     creal_T Y_COMPLEX[FFT_LEN];  //2A T
	
     rev = REV_SUCESS;

// Temp
   
	// Temp var

    if(FrameCnt==1){
     	for(i=0;i< Inlen; i++){	
	    	pcm[i] = (omlsa_float32_t)pIn[i]/32768.0f;		 	
	     }
	}
	else{
     	for(i=0;i< Inlen; i++){	
	    	pcm[i+FRAME_SHIFT] = (omlsa_float32_t)pIn[i]/32768.0f;		 	
	     }
	}


    {
     omlsa_float32_t winy[FRAME_LEN];  //  A
	 for(i=0;i<FRAME_LEN;i++){
		winy[i] = pcm[i]*Window[i];
	  } 

	  fft(winy,Y_COMPLEX);
    }
	for(i=0;i< FRAME_LEN21; i++){	
		Y_2[i] = Y_COMPLEX[i].re * Y_COMPLEX[i].re +  Y_COMPLEX[i].im * Y_COMPLEX[i].im;
	}

	/*  if FrameCnt==1 lambda_d=Y_2;  end*/
	if( FrameCnt==1){
		for(i=0; i<FRAME_LEN21; i++) 
			lambda_dav[i] = lambda_d[i] = Y_2[i] ;
	}
    
	for(i=0; i<FRAME_LEN21; i++){ 
		gamma[i] = Y_2[i] / max_local(lambda_d[i], 1e-10f); 	    
		eta[i] = alpha_eta * eta_2term[i] + (1-alpha_eta) * max_local(gamma[i]-1,0);
		eta[i]=max_local(eta[i],eta_min);
        v[i]  = gamma[i] * eta[i] / (1+eta[i]);    
	}
 
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
 
	 for(i=0; i<FRAME_LEN21; i++) {
		 Sft[i] = St[i];
	 }
   
   {
				
	 omlsa_float32_t I[FRAME_LEN21];
	 omlsa_float32_t Conv_I[FRAME_LEN21+3];			 
     omlsa_float32_t Conv_Y[FRAME_LEN21+3];	


	for(i=0; i<FRAME_LEN21; i++) {
		if(( Y_2[i]<(delta_Y2*Bmin*Smin[i])) && (S[i] < delta_s*Bmin*Smin[i]) ) {
			I[i]= 1; 
		}
		else {
			I[i]= 0;
		}
	}	

	b_conv(I, Conv_I, &size);  
	Conv_I_ptr = &Conv_I[1];

	{
	    omlsa_float32_t Y2_I[FRAME_LEN21];
        omlsa_float32_t* Conv_Y_ptr;

		for(i=0; i<FRAME_LEN21; i++){
			Y2_I[i] =  Y_2[i]*I[i];
		}
		b_conv(Y2_I, Conv_Y, &size);  
		Conv_Y_ptr = &Conv_Y[1];
		 
		for(i=0; i<FRAME_LEN21; i++){
			if( Conv_I_ptr[i])
			   Sft[i] = Conv_Y_ptr[i]/Conv_I_ptr[i];
		}
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
		    St[i]=alpha_s*St[i]+(1-alpha_s)*Sft[i];
            Smint[i]=min_local(Smint[i],St[i]);
            Smactt[i]=min_local(Smactt[i],St[i]);
		 }
	 }

{
	 omlsa_float32_t phat[FRAME_LEN21];	
	 omlsa_float32_t qhat[FRAME_LEN21];				 
     omlsa_float32_t Sr_Y2[FRAME_LEN21];			 
     omlsa_float32_t Sr_S[FRAME_LEN21];	
  
        for(i=0;i<FRAME_LEN21;i++){         
		 	 Sr_Y2[i] = Y_2[i] / Bmin / max_local(Smint[i], 1e-10f);
			 Sr_S[i] = S[i]/Bmin/max_local(Smint[i],1e-10f); 
		}
         
		for(i=0;i<FRAME_LEN21;i++){ 
  
			 if(Sr_Y2[i]>1 && Sr_Y2[i]<delta_yt && Sr_S[i]<delta_s){
				  qhat[i]= (delta_yt-Sr_Y2[i])/(delta_yt-1);
				  phat[i]=1 /(1+ (qhat[i]/(1-qhat[i]))*(1+eta[i])*(float)exp(-v[i]));
			 }
			 else if(Sr_Y2[i]>=delta_yt || Sr_S[i]>=delta_s){
				  phat[i] = 1;
			 }
			 else{
				  qhat[i] = 1; phat[i] = 0;
			 }
		 }
 
         for(i=0;i<FRAME_LEN21;i++){ 
			 alpha_dt[i] = alpha_d + (1- alpha_d) * phat[i];  
             lambda_dav[i] = alpha_dt[i]* lambda_dav[i] + (1- alpha_dt[i]) * Y_2[i];
		 }
 
     if(FrameCnt<15){
		 for(i=0;i<FRAME_LEN21;i++)
           lambda_d_long[i]=lambda_dav[i];
	 }
	 else{ 
		 for(i=0;i<FRAME_LEN21;i++){
			 alpha_dt_long[i] =alpha_d_long +(1-alpha_d_long)*phat[i];
             lambda_d_long[i]=alpha_dt_long[i] *lambda_d_long[i]+(1-alpha_dt_long[i]) *Y_2[i];
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
				  =S_his0[i][4]
				  =S_his0[i][5]
				  =S_his0[i][6]
				  =S_his0[i][7] = S[i];

				   S_his[i][0]
				  =S_his[i][1]
				  =S_his[i][2]
				  =S_his[i][3]
				  =S_his[i][4]
				  =S_his[i][5]
				  =S_his[i][6]
				  =S_his[i][7] = St[i];
			 }
		 }
		 else{
			  for(i=0;i<FRAME_LEN21;i++){

				  Smin[i] = min_local(S_his0[i][6], Smact[i]);
  				  Smin[i] = min_local(S_his0[i][5], Smin[i]);
   			      Smin[i] = min_local(S_his0[i][4], Smin[i]);
                  Smin[i] = min_local(S_his0[i][3], Smin[i]);
  				  Smin[i] = min_local(S_his0[i][2], Smin[i]);
   			      Smin[i] = min_local(S_his0[i][1], Smin[i]);
                  Smin[i] = min_local(S_his0[i][0], Smin[i]);

				  S_his0[i][7] = S_his0[i][6];
  				  S_his0[i][6] = S_his0[i][5];
				  S_his0[i][5] = S_his0[i][4];
				  S_his0[i][4] = S_his0[i][3];
				  S_his0[i][3] = S_his0[i][2];
  				  S_his0[i][2] = S_his0[i][1];
				  S_his0[i][1] = S_his0[i][0];
				  S_his0[i][0] =  Smact[i];
				   
                  Smin[i] = min_local(Smin[i], Smact[i]);
				  Smact[i] = S[i];

				  Smint[i] = min_local(S_his[i][6], Smactt[i]);
  				  Smint[i] = min_local(S_his[i][5], Smint[i]);
   			      Smint[i] = min_local(S_his[i][4], Smint[i]);
                  Smint[i] = min_local(S_his[i][3], Smint[i]);
  				  Smint[i] = min_local(S_his[i][2], Smint[i]);
   			      Smint[i] = min_local(S_his[i][1], Smint[i]);
                  Smint[i] = min_local(S_his[i][0], Smint[i]);

				  S_his[i][7] = S_his[i][6];
  				  S_his[i][6] = S_his[i][5];
				  S_his[i][5] = S_his[i][4];
				  S_his[i][4] = S_his[i][3];
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

    //  lambda_d1=1.4685*lambda_d ;
     for(i=0;i<FRAME_LEN21;i++){
		 lambda_d[i] = 1.4685f*lambda_dav[i];
		 lambda_d_global[i] = lambda_d[i];
	 }


{ 		
     omlsa_float32_t q[FRAME_LEN21];  //B
	 { 
		 
     omlsa_float32_t P_local[FRAME_LEN21];    		    
	 omlsa_float32_t P_global[FRAME_LEN21]; 
	 omlsa_float32_t xi_local[FRAME_LEN21+3];       
     omlsa_float32_t xi_local_dB[FRAME_LEN21];      
	 omlsa_float32_t xi_global[FRAME_LEN21+32];     
	 omlsa_float32_t xi_global_dB[FRAME_LEN21]; 
	  		 //B T
         ex_xi_frame = xi_frame;
		 xi_frame = 0;
     
	  for(i=0;i<FRAME_LEN21;i++){
		  xi[i] = beta * xi[i] + (1-beta) * eta[i]; 
         
		 if(i>=2){
          xi_frame +=  xi[i];
		}
	  }

     xi_frame /=(FRAME_LEN21-2);
 
       b_conv(xi,xi_local,&size);   xi_local_ptr = &xi_local[1];
       c_conv(xi,xi_global,&size ); xi_global_ptr = &xi_global[15];
 
      for(i=0;i<FRAME_LEN21;i++){
		  if(xi_local[i]>0){
			  xi_local_dB[i] = 10*(omlsa_float32_t)log10(xi_local_ptr[i]);
		  }
		  else{
			  xi_local_dB[i] = -100;
		  }
		  if(xi_global[i]>0){
			  xi_global_dB[i] = 10*(omlsa_float32_t)log10(xi_global_ptr[i]);
		  }
		  else{
			  xi_global_dB[i] = -100;
		  }
	  }
	  if(xi_frame>0){
  	    xi_frame_dB = 10*(omlsa_float32_t)log10(xi_frame);
	  }
	  else {
	    xi_frame_dB = -100;
	  }
 
        for(i=0;i<FRAME_LEN21;i++){
			if(xi_local_dB[i]<=xi_min_dB){
				P_local[i] = P_min ;
			}
			else if((xi_local_dB[i]>xi_min_dB) && (xi_local_dB[i]<xi_max_dB)){
				P_local[i] =P_min+(xi_local_dB[i] -xi_min_dB)/(xi_max_dB-xi_min_dB)*(1-P_min); 
			}
			else {
				P_local[i] = 1;
			}
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
 
        for(i=0;i<FRAME_LEN21;i++){
			if(xi_global_dB[i]<=xi_min_dB){
				P_global[i] = P_min ;
			}
			else if(xi_global_dB[i]>xi_min_dB && xi_global_dB[i]<xi_max_dB){
				P_global[i] =P_min+(xi_global_dB[i] -xi_min_dB)/(xi_max_dB-xi_min_dB)*(1-P_min); 
			}
			else {
				P_global[i] = 1;
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

		  q[i] = 1 - P_frame*P_global[i]*P_local[i];
		  q[i] = min_local(qmax,q[i]);
	  }
    }
 
      for(i=0;i<FRAME_LEN21;i++){
          gamma[i] = Y_2[i] / max_local(lambda_d_global[i], 1e-10f);
		  eta[i] = alpha_eta*eta_2term[i] + (1-alpha_eta)*max_local(gamma[i]-1,0);
		  eta[i] = max_local(eta[i],eta_min);
		     v[i]= gamma[i]*eta[i]/(1+eta[i]);    
	  }

{

     omlsa_float32_t PH1[FRAME_LEN21];			 //B T
     omlsa_float32_t GH1[FRAME_LEN21];			 //B T
     omlsa_float32_t GH0[FRAME_LEN21];			 //B T
     omlsa_float32_t G[FRAME_LEN21];				 //B T
				
 
        for(i=0;i<FRAME_LEN21;i++){
           if(q[i]<0.9f){
               PH1[i] = 1/( 1 +(q[i] /(1 - q[i])) * (1+eta[i]) * (omlsa_float32_t)exp(-v[i])  );
		   }
		   else{
			   PH1[i] = 0;
		   }            
		}
 
       for(i=0;i<FRAME_LEN21;i++){
		   if(v[i]>5){
			   GH1[i] = eta[i]/(1+eta[i]);
		   }
		   else if(v[i]<=5 && v[i]>0){
                GH1[i] = (eta[i]/(eta[i]+1))*(omlsa_float32_t)exp(0.5*expint(v[i])); // need
		   }	
		   else {
		        GH1[i] = 1;
		   }
	   }
 
        for(i=3;i<FRAME_LEN21-3;i++){
			omlsa_float32_t minval;
  
			minval = min_local(lambda_d[i] , lambda_d[i-3] );
            lambda_d_global[i] = min_local(minval , lambda_d[i+3] );
		}
        
		for(i=0;i<FRAME_LEN21;i++){
			Sy[i]=0.8f*Sy[i]+0.2f*Y_2[i];    
		}
		
         for(i=0;i<FRAME_LEN21;i++){
            GH0[i]= Gmin * (omlsa_float32_t)sqrt(lambda_d_global[i]/(Sy[i]+ 1e-10f));         			
		 }

		  for(i=0;i<FRAME_LEN21;i++){      
		    G[i]= powf(GH1[i],PH1[i]) * (omlsa_float32_t)powf(GH0[i] ,(1-PH1[i]));  
		  }
       
         for(i=0;i<FRAME_LEN21;i++){
             eta_2term[i]=GH1[i]*GH1[i]*gamma[i];
		 }
  

        for(i=0;i<FRAME_LEN21;i++){
		
			if(i<3 ||i==(FRAME_LEN21-1)){
			   Y_COMPLEX[i].re = 0;
  			   Y_COMPLEX[i].im = 0;
			}
			else {
 			   Y_COMPLEX[i].re = Y_COMPLEX[i].re * G[i];
			   Y_COMPLEX[i].im = Y_COMPLEX[i].im * G[i];			
			}
		}  


        for(i=FRAME_LEN21;i<FRAME_LEN ;i++){
			Y_COMPLEX[i].re = Y_COMPLEX[FRAME_LEN-i].re;
			Y_COMPLEX[i].im = -Y_COMPLEX[FRAME_LEN-i].im;
		}
		 
        {
		    creal_T x_g_c[FFT_LEN];
	
		    ifft(Y_COMPLEX, x_g_c  );
  
		   for(i=0;i<FRAME_LEN;i++){	 
		    	x_g_c[i].re =  x_g_c[i].re *Window[i]* Cwin*Cwin;
		    }
  
		    for(i=0;i<FRAME_LEN;i++){ 
				out_buf[i] = out_buf[i] + x_g_c[i].re;
		    }
        }
		} 	
		}
		for(i=0;i<FRAME_LEN/4;i++){
			omlsa_float32_t temp;

			temp = (omlsa_float32_t)( out_buf[i]* 32768 );

			temp = (temp>32767)? 32767: ((temp<-32768)? -32768:temp);

			pOut[i] =  (short)temp;
		}

		for(i=0;i<FRAME_LEN;i++){

			 if(i<FRAME_SHIFT)
			   out_buf[i] = out_buf[i+FRAME_LEN/4];
			else 
			    out_buf[i] = 0;
		}
		 
	    for(i=0;i< FRAME_SHIFT; i++){	
		    pcm[i ] = pcm[i +	FRAME_LEN41]; //   [[...] . ]  <- [.[...]]
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





 