#include"omlsa_def.h"
#include"omlsa_var.h"
#include"omlsa_tab.h"

#include<stdio.h>




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


 omlsa_float32_t omlsa_shared_buf[1300];
 
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

void PostFilterInit(){
	int i;	 
	   

	Cwin = 11.313709f;
	 
	for(i=0;i< FRAME_LEN21; i++) eta_2term[i] = 1; 
}
   int I_empty= 1;  omlsa_float32_t m_P_local;
   omlsa_float32_t *xi_local_ptr, *xi_global_ptr;omlsa_float32_t xi_peak_dB;
    	 	
#define REV_SUCESS 0
    
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
	 creal_T *x_g_c;

    omlsa_float32_t gamma ;			 // omlsa_float32_t gamma[FRAME_LEN21];			 
     omlsa_float32_t v[FRAME_LEN21];				 
     creal_T Y_COMPLEX[FFT_LEN];       
	 omlsa_float32_t Y_2[FRAME_LEN21];	
	
     rev = REV_SUCESS;


     	for(i=0;i< Inlen; i++){		    
			pcm[i+FRAME_SHIFT] = pIn[i] ;		 		 	
	     }
  
	   winy = om_mem_alloc(512); // alloc  1
	    
	   for(i=0;i<FRAME_LEN;i++){
		  winy[i] = (omlsa_float32_t)pcm[i]*Window[i];// / 32768.0f;
	   } 
	   


	   fft(winy,Y_COMPLEX);	   
	   om_mem_free(512);	 // free 1  
    
	   // start section 1


 
    {

#if 1		
    omlsa_float32_t* pf1,*pf2, *pf3;
 
	pf1 = &Y_COMPLEX[0].re;  // pf2 = &Y_COMPLEX[0].im;
	pf3 = Y_2;

	for(i=0;i< FRAME_LEN21; i++){	
	 	*pf3      = (*pf1) * (*pf1++);
	 	(*pf3++) += (*pf1) * (*pf1++); 		 
	}

#else 
	for(i=0;i< FRAME_LEN21; i++){	
		Y_2[i] = Y_COMPLEX[i].re * Y_COMPLEX[i].re +  Y_COMPLEX[i].im * Y_COMPLEX[i].im;
	}
#endif
    }

 
	if( FrameCnt==1){
		for(i=0; i<FRAME_LEN21; i++) 
			lambda_dav[i] = lambda_d[i] = Y_2[i] ;
	}
    
//	for(i=0; i<FRAME_LEN21; i++){ 	
//		gamma = Y_2[i] / max_local(lambda_d[i], 1e-10f); 	    
//		eta[i] = alpha_eta * eta_2term[i] + alpha_eta_C * max_local(gamma-1,0.0f);
//		eta[i] = max_local(eta[i],eta_min);
//        v[i] = gamma * eta[i] / (1+eta[i]);   	 
//	}
 

    //omlsa_float32_t Sf[FRAME_LEN21+3];			//B
    Sf = om_mem_alloc(FRAME_LEN21+3);  // alloc 1
	b_conv(Y_2, Sf, &size);  
	Sf_ptr = &Sf[1];      
 
	if(FrameCnt==1){
		for(i=0; i<FRAME_LEN21; i++){		  
			Sy[i] = Y_2[i];	 S[i] = Sf_ptr[i];  St[i] = Sf_ptr[i];
		}
	}
	else{
		for(i=0; i<FRAME_LEN21; i++){	
			S[i] = alpha_s * S[i] + alpha_s_C*Sf_ptr[i];
		}
	}

	om_mem_free(FRAME_LEN21+3); // free 1
  	
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

	 //omlsa_float32_t Sft[FRAME_LEN21];			//B
	 Sft = om_mem_alloc(FRAME_LEN21);  // alloc 1
 
	 for(i=0; i<FRAME_LEN21; i++) {
		 Sft[i] = St[i];
	 }
   
     	 
	I = om_mem_alloc(FRAME_LEN21);
	Conv_I = om_mem_alloc(FRAME_LEN21+3); // alloc 2
	Conv_Y = om_mem_alloc(FRAME_LEN21+3); // alloc 3
 
	for(i=0; i<FRAME_LEN21; i++) {
		if(( Y_2[i]<(delta_Y2_S*Smin[i])) && (S[i] < delta_s_S*Smin[i]) ) {
			I[i]= 1.0f; 
		}
		else {
			I[i]= 0.0f;
		}
	}	

	b_conv(I, Conv_I, &size);  
	Conv_I_ptr = &Conv_I[1];
  
		Y2_I = om_mem_alloc(FRAME_LEN21); // alloc 4

		for(i=0; i<FRAME_LEN21; i++){
			Y2_I[i] =  Y_2[i]*I[i];
		}
		b_conv(Y2_I, Conv_Y, &size);  
		Conv_Y_ptr = &Conv_Y[1];
		 
		for(i=0; i<FRAME_LEN21; i++){
			if( Conv_I_ptr[i])
			   Sft[i] = Conv_Y_ptr[i]/Conv_I_ptr[i];
		}
 

	     om_mem_free(FRAME_LEN21); // free 4

	 

	 om_mem_free(FRAME_LEN21*3 + 6);  // free 3 2
 
     if(FrameCnt<15){
	     for(i=0;i<FRAME_LEN21;i++){
  		     St[i] = S[i];
			 Smint[i] = St[i];
			 Smactt[i] = St[i];
		 }		  
	 }
	 else{
		 for(i=0;i<FRAME_LEN21;i++){
		    St[i]=alpha_s*St[i]+alpha_s_C*Sft[i];
            Smint[i] = min_local(Smint[i],St[i]);
            Smactt[i]= min_local(Smactt[i],St[i]);
		 }
	 }
	  
	 om_mem_free(FRAME_LEN21); // free 1
	  
	for(i=0;i<FRAME_LEN21;i++){ 

		 omlsa_float32_t phat , qhat, Sr_Y2,Sr_S , alpha_dt, alpha_dt_long,  Sr_f_B;
		 omlsa_float32_t  eta,v;

         gamma = Y_2[i] / max_local(lambda_d[i], 1e-10f); 	    
         eta = alpha_eta * eta_2term[i] + alpha_eta_C * max_local(gamma-1,0.0f);
         eta = max_local(eta,eta_min);
          v = gamma * eta / (1+eta);   	 
 
#if OPT_SEC1
		  Sr_f_B = Bmin_Inv / max_local(Smint[i], 1e-10f);
  
		  Sr_Y2  = Y_2[i] *Sr_f_B;
          Sr_S  = S[i]*Sr_f_B; 
#else 
		  Sr_Y2  = Y_2[i] / Bmin / max_local(Smint[i], 1e-10f);
          Sr_S  = S[i]/Bmin/max_local(Smint[i],1e-10f); 
#endif
		 if(Sr_Y2 >1 && Sr_Y2 <delta_yt && Sr_S <delta_s){
			  qhat= (delta_yt-Sr_Y2 )/(delta_yt-1.0f);
			 
#if OPT_SEC1
			   phat =(1.0f-qhat) /((1.0f-qhat) + (qhat)*(1.0f+eta )*(float)expf(-v ));
#else 
			  phat =1.0f /(1.0f + (qhat/(1.0f-qhat))*(1.0f+eta[i])*(float)expf(-v[i]));
#endif
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

#if OPT_SEC1
			  alpha_dt_long =alpha_d_long +alpha_d_long_C*phat ;   
#else 
			  alpha_dt_long =alpha_d_long +(1-alpha_d_long)*phat ;   
#endif
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




	 // start of section 2
     { 		 
	   P_local   = om_mem_alloc(FRAME_LEN21);
	   P_global  = om_mem_alloc(FRAME_LEN21);
	   xi_local  = om_mem_alloc(FRAME_LEN21+3);
	   xi_global = om_mem_alloc(FRAME_LEN21+32);
 
       ex_xi_frame = xi_frame;
	   xi_frame = 0;
     

#if OPT_SEC2
	   xi[0] = beta * xi[0] + beta_C * eta[0];         
	   xi[1] = beta * xi[1] + beta_C * eta[1];  

	  for(i=2;i<FRAME_LEN21;i++){ 
		  xi[i] = beta * xi[i] + beta_C * eta[i];           	 
          xi_frame +=  xi[i];
		
	  }
#else 
	  for(i=0;i<FRAME_LEN21;i++){ 
		  xi[i] = beta * xi[i] + (1-beta) * eta[i];           
		 if(i>=2){
           xi_frame +=  xi[i];
		}
	  }
#endif
	   
       xi_frame /=(FRAME_LEN21-2);
 
       b_conv(xi,xi_local,&size);   xi_local_ptr = &xi_local[1];
       c_conv(xi,xi_global,&size ); xi_global_ptr = &xi_global[15];
 
      for(i=0;i<FRAME_LEN21;i++){
          omlsa_float32_t xi_global_dB ,xi_local_dB; 
		 
		  if(xi_local[i]>0){
#if OPT_SEC2
			  xi_local_dB =  4.3429447587231051f*(omlsa_float32_t)logf(xi_local_ptr[i]);
#else 
			  xi_local_dB = 10*(omlsa_float32_t)log10(xi_local_ptr[i]);
#endif
		  }
		  else{
			  xi_local_dB  = -100.0f;
		  }

			if(xi_local_dB <=xi_min_dB){
				P_local[i] = P_min ;
			}
			else if((xi_local_dB >xi_min_dB) && (xi_local_dB <xi_max_dB)){
#if OPT_SEC2
				P_local[i] =P_min+(xi_local_dB -xi_min_dB)  * 0.201005f;// /(xi_max_dB-xi_min_dB)*(1-P_min); 
#else 
                P_local[i] =P_min+(xi_local_dB -xi_min_dB)/(xi_max_dB-xi_min_dB)*(1-P_min); 
#endif
			}
			else {
				P_local[i] = 1;
			}

		  if(xi_global[i]>0){
#if OPT_SEC2
			  xi_global_dB  = 4.3429447587231051f*(omlsa_float32_t)logf(xi_local_ptr[i]);
#else 
			   xi_global_dB  = 10*(omlsa_float32_t)log10(xi_global_ptr[i]);
#endif
		  }
		  else{
			  xi_global_dB  = -100;
		  }

		  if(xi_global_dB <=xi_min_dB){
		  	P_global[i] = P_min ;
		  }
		  else if(xi_global_dB >xi_min_dB && xi_global_dB <xi_max_dB){
#if OPT_SEC2
		  	P_global[i] = P_min+(xi_global_dB  -xi_min_dB) * 0.201005f;
#else 
            P_global[i] = P_min+(xi_global_dB  -xi_min_dB)/(xi_max_dB-xi_min_dB)*(1-P_min); 
#endif
		  }
		  else {
		  	P_global[i] = 1;
		  }
	  }


	  om_mem_free(FRAME_LEN21*2+35); // set free  xi_local xi_global

	  
	  if(xi_frame>0){
#if OPT_SEC2
  	    xi_frame_dB =  4.3429447587231051f*(omlsa_float32_t)logf(xi_frame);
#else 
  	    xi_frame_dB = 10*(omlsa_float32_t)log10(xi_frame);

#endif
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
#if OPT_SEC2
            P_frame=P_min+(xi_frame_dB-xi_min_dB-xi_peak_dB)* 0.201005f;
#else 
            P_frame=P_min+(xi_frame_dB-xi_min_dB-xi_peak_dB)/(xi_max_dB-xi_min_dB)*(1-P_min);
#endif
	   }

	   // end of section 2



	   // start of secton 3
 

      {
 
		 
		lambda_d_global = om_mem_alloc(FRAME_LEN21);
 
		
 
 		lambda_d_global[0]=  lambda_d[0];
		lambda_d_global[1]=  lambda_d[1];
		lambda_d_global[2]=  lambda_d[2];
		lambda_d_global[FRAME_LEN21-3]=  lambda_d[FRAME_LEN21-3];
		lambda_d_global[FRAME_LEN21-2]=  lambda_d[FRAME_LEN21-2];
		lambda_d_global[FRAME_LEN21-1]=  lambda_d[FRAME_LEN21-1];

        for(i=3;i<FRAME_LEN21-3;i++){
			omlsa_float32_t minval;
  
			minval = min_local(lambda_d[i] , lambda_d[i-3] );
            lambda_d_global[i] = min_local(minval , lambda_d[i+3] );
		}
		 
        for(i=0;i<FRAME_LEN21;i++){	
			omlsa_float32_t q, PH1, GH0,GH1,G,v;//,eta;
 		
           q = 1.0f - P_frame*P_global[i]*P_local[i];
           q = min_local(qmax,q);

          gamma = Y_2[i] / max_local(lambda_d[i], 1e-10f);

#if 0
          eta = alpha_eta*eta_2term[i] + 0.05f*max_local(gamma-1.0f,0.0f);
          eta = max_local(eta,eta_min);
          v[i]   = gamma*eta/(1.0f+eta);  
#else 
          eta[i] = alpha_eta*eta_2term[i] + 0.05f*max_local(gamma-1.0f,0.0f);
          eta[i] = max_local(eta[i],eta_min);
          v    = gamma*eta[i]/(1.0f+eta[i]); 
#endif


           if(q<0.9f){
 
               PH1  = (1.0f - q)/( (1.0f - q) +  q  * (1.0f+eta[i]) * (omlsa_float32_t)expf(-v )  );
 
		   }
		   else{
			   PH1  = 0;
		   }        


#if 0
           if(v[i]>5){
  	          GH1= eta/(1.0f+eta);
           }
           else if(v[i]<=5.0f && v[i]>0.0f){ 
                GH1= (eta/(eta+1.0f))*(omlsa_float32_t)expf(0.5f*expint(v[i])); // need
           }	
#else 
           if(v[i]>5){
  	          GH1= eta[i]/(1.0f+eta[i]);
           }
           else if(v[i]<=5.0f && v[i]>0.0f){ 
                GH1= (eta[i]/(eta[i]+1.0f))*(omlsa_float32_t)expf(0.5f*expint(v )); // need
           }	
#endif
           else {
                GH1= 1.0f;
           }

		      Sy[i]=0.8f*Sy[i]+0.2f*Y_2[i];   
 
 
			  GH0= Gmin * (omlsa_float32_t)sqrtf(lambda_d_global[i]/(Sy[i]+ 1e-10f));         	
 
			  G= powf(GH1 ,PH1 ) * (omlsa_float32_t)powf(GH0 ,(1-PH1 ));  
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
		    
        }
           
		   om_mem_free(3*FRAME_LEN21);


		   // enf of section 3

        for(i=FRAME_LEN21;i<FRAME_LEN ;i++){
			Y_COMPLEX[i].re = Y_COMPLEX[FRAME_LEN-i].re;
			Y_COMPLEX[i].im = -Y_COMPLEX[FRAME_LEN-i].im;
		}
		 
        {
			 
		    x_g_c = ( creal_T *)om_mem_alloc(FFT_LEN*2);
	  
		    ifft(Y_COMPLEX, x_g_c  );
  
		  // for(i=0;i<FRAME_LEN;i++){	 
		  //  	x_g_c[i].re =  x_g_c[i].re *Window[i]*Cwin_2 ;//*32768;
		  //  }
  
#if 1
			{
            omlsa_float32_t *p0, *p1;
            short *p2;

			p0 = (omlsa_float32_t*) &x_g_c[0].re;
            p1 = &Window[0];
			p2 = out_buf;

		    for(i=0;i<128;i++){ 
				int temp;
				temp = (int)(*p0++ *  (*p1++) * 128) + *p2; 
				p0++; 
				*p2 =(short)((temp>32767)? 32767: ((temp<-32768)? -32768:temp));				
				*pOut++ =  (short)*p2++;
		    }	

			  for(i=0;i<384;i++){ 
				int temp;
				temp = (int)(*p0++ * (*p1++)*128) + *p2; 
				p0++;
				*p2++ =(short)((temp>32767)? 32767: ((temp<-32768)? -32768:temp));
								 				   
		    }	
			}

	 
	#else 
		    for(i=0;i<FRAME_LEN;i++){ 

				int temp;

				temp = (int)( x_g_c[i].re *Window[i]*Cwin_2) + out_buf[i];//*32768;
 
				out_buf[i] =(short)((temp>32767)? 32767: ((temp<-32768)? -32768:temp));
				
				if(i<128)
				    pOut[i] =  (short)out_buf[i];
		    }
    #endif	
			om_mem_free(FFT_LEN*2);
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





 