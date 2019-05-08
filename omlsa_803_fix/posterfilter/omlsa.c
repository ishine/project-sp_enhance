#include"omlsa_def.h"
#include"omlsa_var.h"
#include"omlsa_tab.h"

#include<stdio.h>

void gpio_on(){};
  void gpio_on1(){};
  void gpio_on2(){};
  void gpio_off(){};
  void gpio_off1(){};
  void gpio_off2(){};

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
 
void fft_opt_512(float* fftinput ,float* fftoutput);

int MinSeg = 0;


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


omlsa_float32_t Cwin_2 ;

void PostFilterInit(){
	int i;	 
	   

	Cwin = 11.313709f;
	Cwin_2 = 128;
	 
	for(i=0;i< FRAME_LEN21; i++) eta_2term[i] = 1; 
}
   int I_empty= 1;  omlsa_float32_t m_P_local;
   omlsa_float32_t *xi_local_ptr, *xi_global_ptr;omlsa_float32_t xi_peak_dB;
    	 	
#define REV_SUCESS 0

omlsa_float32_t gamma,eta;        //[FRAME_LEN21];			 
omlsa_float32_t v;                //[FRAME_LEN21];				 
creal_T Y_COMPLEX[FRAME_LEN21];       //     
omlsa_float32_t Y_2[FRAME_LEN21];	
     	  
void ifft_opt_512(float* fftinput ,float* fftoutput);


#define fft_512_opt 1

int PostFilterProcess(short*pIn, int Inlen, short*pOut, int*Outlen){

	int i,rev ,size; omlsa_float32_t *Sf_ptr, *Conv_I_ptr;

	 omlsa_float32_t *winy;  //  A
     omlsa_float32_t *I; 
     omlsa_float32_t *Conv_I;
     omlsa_float32_t *Conv_Y;
     omlsa_float32_t *P_local ;    	
     omlsa_float32_t *P_global ; 
     omlsa_float32_t *xi_local ;     
     omlsa_float32_t *xi_global ; 
     omlsa_float32_t *lambda_d_global;	
	 omlsa_float32_t *Sf;
	 omlsa_float32_t *Sft;
     omlsa_float32_t *Y2_I;
     omlsa_float32_t* Conv_Y_ptr;
	 float *x_g_c;
	  
     rev = REV_SUCESS;
	 
	  gpio_on();
      gpio_on1();

     	for(i=0;i< Inlen; i++){		    
			pcm[i+FRAME_SHIFT] = pIn[i] ;		 		 	
	     }
  
	   winy = om_mem_alloc(512); // alloc  1
	   
	   {
		   omlsa_float32_t* p0, *p1,*p2;
	    
	   for(i=0;i<FRAME_LEN;i++){
		  winy[i] = (omlsa_float32_t)pcm[i]*Window[i];// / 32768.0f;		   		   
	   } 
	   
	   }

  
	   fft_opt_512(winy, (omlsa_float32_t*)Y_COMPLEX );
 
 
    
 	 gpio_off1();
	   // start section 1
      
   om_mem_free(512);	 // free 1  

 gpio_on2();
  
	 Y_2[0]   = Y_COMPLEX[0].re * Y_COMPLEX[0].re ;
	 Y_2[256] = Y_COMPLEX[0].im * Y_COMPLEX[0].im ;

	for(i=1;i< FRAME_LEN21-1; i++){	
		Y_2[i] = Y_COMPLEX[i].re * Y_COMPLEX[i].re +  Y_COMPLEX[i].im * Y_COMPLEX[i].im;
	}
  
	if( FrameCnt==1){
		for(i=0; i<FRAME_LEN21; i++) 
			lambda_dav[i] = lambda_d[i] = Y_2[i] ;
	}
      
 //0.138

    //omlsa_float32_t Sf[FRAME_LEN21+3];			//B
    Sf = om_mem_alloc(FRAME_LEN21+3);  // alloc Sf
	b_conv(Y_2, Sf, &size);  
	Sf_ptr = &Sf[1];      
 
	if(FrameCnt==1){
		for(i=0; i<FRAME_LEN21; i++){		  
		//	Sy[i] = Y_2[i];	 S[i] = Sf_ptr[i];  St[i] = Sf_ptr[i];
			Sy[i] = Y_2[i];	 S[i] = *Sf_ptr;  St[i] = *Sf_ptr++;
		}
	}
	else{
		for(i=0; i<FRAME_LEN21; i++){	
		//	S[i] = alpha_s * S[i] + alpha_s_C*Sf_ptr[i];
			S[i] = alpha_s * S[i] + alpha_s_C*   *Sf_ptr++;
		}
	}

	om_mem_free(FRAME_LEN21+3); // free Sf
  	
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
  
	I = om_mem_alloc(FRAME_LEN21);        // alloc I
	Conv_I = om_mem_alloc(FRAME_LEN21+3); // alloc Conv_I
	Conv_Y = om_mem_alloc(FRAME_LEN21+3); // alloc Conv_Y
 
	for(i=0; i<FRAME_LEN21; i++) {
		if(( Y_2[i]<(delta_Y2_S*Smin[i])) && (S[i] < delta_s_S*Smin[i]) ) {
			I[i]= 1.0f; 
		}
		else {
			I[i]= 0.0f;
		}
	}	

	b_conv(I, Conv_I, &size);  
	
  
		Y2_I = om_mem_alloc(FRAME_LEN21); // alloc Y2_I

		for(i=0; i<FRAME_LEN21; i++){
			Y2_I[i] =  Y_2[i]*I[i];
		}
		b_conv(Y2_I, Conv_Y, &size);  
 
		Conv_Y_ptr = &Conv_Y[1];
        Conv_I_ptr = &Conv_I[1];
  //46us

      //om_mem_free(FRAME_LEN21); // free 4
  
     if(FrameCnt<15){
	     for(i=0;i<FRAME_LEN21;i++){
  		     St[i] = S[i];
			 Smint[i] = St[i];
			 Smactt[i] = St[i];
		 }		  
	 }
	 else{
		 for(i=0;i<FRAME_LEN21;i++){
             omlsa_float32_t Sft;

			 if(Conv_I_ptr[i]){
                St[i]=alpha_s*St[i]+alpha_s_C*Conv_Y_ptr[i]/Conv_I_ptr[i];
			 }
            else {
				St[i]=alpha_s*St[i]+alpha_s_C*St[i];
			}
			 
            Smint[i] = min_local(Smint[i],St[i]);
            Smactt[i]= min_local(Smactt[i],St[i]);
		 }
	 }

        om_mem_free(FRAME_LEN21*4 + 6);  // free  
	 
      
	 ex_xi_frame = xi_frame;     
	 xi_frame = 0;
   
//0.176ms    
	for(i=0;i<FRAME_LEN21;i++){ 

		omlsa_float32_t phat , qhat, Sr_Y2,Sr_S , alpha_dt, alpha_dt_long,  Sr_f_B;
		  
	    gamma = Y_2[i] / max_local(lambda_d[i], 1e-10f); 	    
	    eta = alpha_eta * eta_2term[i] + alpha_eta_C * max_local(gamma-1.0f,0.0f);
	    eta = max_local(eta,eta_min);
        v = gamma * eta / (1.0f+eta);   

  
        xi[i] = beta * xi[i] + (1-beta) * eta;   
        if(i>=2){
            xi_frame +=  xi[i];
		}
 
 
		  Sr_f_B = Bmin_Inv / max_local(Smint[i], 1e-10f);
  
		  Sr_Y2  = Y_2[i] *Sr_f_B;
          Sr_S  = S[i]*Sr_f_B; 
 
		 if(Sr_Y2 >1.0f && Sr_Y2 <delta_yt && Sr_S <delta_s){
			  
 
               qhat= (delta_yt-Sr_Y2 )*0.5f;
			   phat =(1.0f-qhat) /((1.0f-qhat) + (qhat)*(1.0f+eta)*(float)expf(-v));
 
		 }
		 else if(Sr_Y2 >=delta_yt || Sr_S >=delta_s){
			  phat  = 1.0f;
		 }
		 else{
			  phat  = 0.0f;
		 }

	     alpha_dt = alpha_d + alpha_d_C * phat ;  
         lambda_dav[i] = alpha_dt* lambda_dav[i] + (1.0f- alpha_dt) * Y_2[i];

		  if(FrameCnt<15){
			  lambda_d_long[i]=lambda_dav[i];
		  }
		  else{ 
			  alpha_dt_long =alpha_d_long +alpha_d_long_C*phat ;    
              lambda_d_long[i]=alpha_dt_long *lambda_d_long[i]+(1-alpha_dt_long) *Y_2[i]; 
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
	  
     for(i=0;i<FRAME_LEN21;i++){
		 lambda_d[i] = 1.4685f*lambda_dav[i];	 
	 }

	 // end of section 1

 gpio_off2();
 gpio_on1();

	 // start of section 2
     { 		 
	   P_local   = om_mem_alloc(FRAME_LEN21);	 
	   xi_local  = om_mem_alloc(FRAME_LEN21+3);
	  
  
       xi_frame /=(FRAME_LEN21-2);
 
       b_conv(xi,xi_local,&size);   xi_local_ptr = &xi_local[1];
       
 
      for(i=0;i<FRAME_LEN21;i++){
          omlsa_float32_t xi_local_dB  ; 
		 
		  if(xi_local[i]>0.0f){ 
			  xi_local_dB =  4.3429447587231051f*(omlsa_float32_t)logf(xi_local_ptr[i]); 
		  }
		  else{
			  xi_local_dB  = -100.0f;
		  }

			if(xi_local_dB <=xi_min_dB){
				P_local[i] = P_min ;
			}
			else if((xi_local_dB >xi_min_dB) && (xi_local_dB <xi_max_dB)){ 
				P_local[i] =P_min+(xi_local_dB -xi_min_dB)  * 0.201005f;// /(xi_max_dB-xi_min_dB)*(1-P_min);  
			}
			else {
				P_local[i] = 1;
			}
	  }

      om_mem_free(FRAME_LEN21+3); 
 
      xi_global = om_mem_alloc(FRAME_LEN21+32);
      P_global  = om_mem_alloc(FRAME_LEN21);

	  c_conv(xi,xi_global,&size ); xi_global_ptr = &xi_global[15];

	  for(i=0;i<FRAME_LEN21;i++){
          omlsa_float32_t  xi_global_dB;

		  if(xi_global[i]>0){  
			  xi_global_dB  = 4.3429447587231051f*(omlsa_float32_t)logf(xi_global_ptr[i]);			   
		  }
		  else{
			  xi_global_dB  = -100;
		  }

		  if(xi_global_dB <=xi_min_dB){
		  	P_global[i] = P_min ;
		  }
		  else if(xi_global_dB >xi_min_dB && xi_global_dB <xi_max_dB){
 
		  	P_global[i] = P_min+(xi_global_dB  -xi_min_dB) * 0.201005f;
 
		  }
		  else {
		  	P_global[i] = 1;
		  }
	  }
 
       om_mem_free(FRAME_LEN21+32); // set free  xi_local xi_global
	  
	  if(xi_frame>0){
 
  	    xi_frame_dB =  4.3429447587231051f*(omlsa_float32_t)logf(xi_frame);
 
	  }
	  else {
	    xi_frame_dB = -100.0f;
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
		    P_frame=1.0f;    
	   }
	   else if(xi_frame_dB>=xi_peak_dB+xi_max_dB ){
		    P_frame=1.0f;
	   }
	   else if(xi_frame_dB<=xi_peak_dB+xi_min_dB ){
            P_frame=P_min;     
	   }
	   else{
 
            P_frame=P_min+(xi_frame_dB-xi_min_dB-xi_peak_dB)* 0.201005f;
 
	   }

	   // end of section 2
   gpio_off1();
   
   
    gpio_on2();


	   // start of secton 3
 
        for(i=0;i<FRAME_LEN21;i++){	
			omlsa_float32_t q, PH1, GH0,GH1,G,  lambda_d_global;


			if((i<=2)||(i>=(FRAME_LEN21-3))){
			    lambda_d_global = lambda_d[i];
			}
			else{
			    omlsa_float32_t minval;
			    minval = min_local(lambda_d[i] , lambda_d[i-3] );
                lambda_d_global  = min_local(minval , lambda_d[i+3] );
			}


            gamma = Y_2[i] / max_local(lambda_d[i], 1e-10f);
            eta = alpha_eta*eta_2term[i] + 0.05f*max_local(gamma-1.0f,0.0f);
            eta = max_local(eta,eta_min);
            v   = gamma*eta/(1.0f+eta);    
            
  
           q = 1.0f - P_frame*  (*P_global++)  * (*P_local++);
 
           q = min_local(qmax,q);
 
           if(q<0.9f){ 
               PH1  = (1.0f - q)/( (1.0f - q) +  q  * (1.0f+eta) * (omlsa_float32_t)expf(-v)  );
		   }
		   else{
			   PH1  = 0;
		   }        

           if(v>5.0f){
  	          GH1= eta/(1.0f+eta);
           }
           else if(v<=5.0f && v>0.0f){ 
              GH1= (eta/(eta+1.0f))*(omlsa_float32_t)expf(0.5f*expint(v)); // need
           }	
           else {
              GH1= 1.0f;
           }

		   Sy[i]=0.8f*Sy[i]+0.2f*Y_2[i];   
 
		   GH0= Gmin * (omlsa_float32_t)sqrtf(lambda_d_global /(Sy[i]+ 1e-10f));         	
 
		   G= powf(GH1 ,PH1 ) * (omlsa_float32_t)powf(GH0 ,(1-PH1 ));  

           eta_2term[i]=GH1 *GH1 *gamma;

			if(i<3 ||i==(FRAME_LEN21-1)){
			   Y_COMPLEX[i].re = 0.0f;
  			   Y_COMPLEX[i].im = 0.0f;
			}
			else {
 			   Y_COMPLEX[i].re = Y_COMPLEX[i].re *G;
			   Y_COMPLEX[i].im = Y_COMPLEX[i].im *G;			
			}
		}
  
		   om_mem_free(2*FRAME_LEN21);


 gpio_off2();

		   // enf of section 3
 
        Y_COMPLEX[0].re = 0.0f; 
        Y_COMPLEX[0].im = 0.0f;
        Y_COMPLEX[256].re = 0.0f; 
        Y_COMPLEX[256].im = 0.0f;
  
        {
			 
		    x_g_c = ( float *)om_mem_alloc(FFT_LEN );
	   gpio_on2();//0.622ms
		    
	  
	        ifft_opt_512((omlsa_float32_t*)Y_COMPLEX,(omlsa_float32_t*)x_g_c );
			 
   
       gpio_off2();
	
	// 88 us

 
{
            omlsa_float32_t *p0, *p1;
            short *p2;
            int temp;
			p0 = (omlsa_float32_t*) x_g_c;
            p1 = &Window[0];
			p2 = out_buf;

		    for(i=0;i<64;i++){ 
				 
				temp = (int)(*p0++ *  (*p1++) * Cwin_2) + *p2; 
				//p0++; 
				*p2 =(short)((temp>32767)? 32767: ((temp<-32768)? -32768:temp));				
				*pOut++ =  (short)*p2++;

				temp = (int)(*p0++ *  (*p1++) * Cwin_2) + *p2; 
				//p0++; 
				*p2 =(short)((temp>32767)? 32767: ((temp<-32768)? -32768:temp));				
				*pOut++ =  (short)*p2++;				
		    }	

			  for(i=0;i<192;i++){ 
				
				temp = (int)(*p0++ * (*p1++)*Cwin_2) + *p2; 
			//	p0++;
				*p2++ =(short)((temp>32767)? 32767: ((temp<-32768)? -32768:temp));

				temp = (int)(*p0++ * (*p1++)*Cwin_2) + *p2; 
			//	p0++;
				*p2++ =(short)((temp>32767)? 32767: ((temp<-32768)? -32768:temp));

		    }					
}
	 
 	
			om_mem_free(FFT_LEN );
        }
		
		}
	  
 
	   // 36 us
 
        {
			short *p0, *p1, *p2, *p3;

			p0 = out_buf;
			p1 = &out_buf[FRAME_LEN41];
			p2 = pcm;
			p3 = &pcm[FRAME_LEN41];

		for(i=0;i<96;i++){		 
			*p0++ = *p1++; 	*p0++ = *p1++; *p0++ = *p1++; *p0++ = *p1++;
			*p2++ = *p3++; 	*p2++ = *p3++; *p2++ = *p3++; *p2++ = *p3++;				 
		}
          
		  p0 = &out_buf[FRAME_SHIFT];

	    for(i= 0;i<32; i++){         
			 *p0++ = 0; *p0++ = 0; *p0++ = 0; *p0++ = 0;			 			 
		}

		}
	  
	  
 
		 gpio_off();
		FrameCnt++;

		return rev;
}
 
 




 