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

  (3) speech absent probability   xi -> xi_ xi_global, xi_frame ->q   

  (4) gamma, eta, v, GH1
*/
 
void fft(const double x[512], creal_T  y[512]);
 
int MinSeg = 0;

void PostFilterInit(){
	int i;	 
	HannWinMake(S_win, 3);

	HammingWin75Overlap(FRAME_LEN);
	 FftWin();
	for(i=0;i< FRAME_LEN21; i++) eta_2term[i] = 1; 
}
   int I_empty= 1;  double m_P_local;
     double *xi_local_ptr, *xi_global_ptr;double xi_peak_dB;

int PostFilterProcess(short*pIn, int Inlen, short*pOut, int*Outlen){

	int i, jj; int size;

	double *Sf_ptr, *Conv_I_ptr;

 
    if(FrameCnt==1){
     	for(i=0;i< Inlen; i++){	
	    	pcm[i] = (double)pIn[i]/32768.0;		 	
	     }
	}
	else{
     	for(i=0;i< Inlen; i++){	
	    	pcm[i+FRAME_SHIFT] = (double)pIn[i]/32768.0;		 	
	     }
	}

	for(i=0;i<FRAME_LEN;i++){

		winy[i] = pcm[i]*Window[i];
	}


	D_M_ANC_WIN_FFT(pcm,  FRAME_LEN, Y_real,Y_imag  );


	fft(winy,Y_COMPLEX);
	

#if 1
 
	//  |Y|^2
	for(i=0;i< FRAME_LEN21; i++){	
	//	Y_2[i] = Y_real[i]*Y_real[i] + Y_imag[i]*Y_imag[i];		
		Y_2[i] = Y_COMPLEX[i].re * Y_COMPLEX[i].re +  Y_COMPLEX[i].im * Y_COMPLEX[i].im;

	}

	/*  if FrameCnt==1 lambda_d=Y_2;  end*/
	if( FrameCnt==1){
		for(i=0; i<FRAME_LEN21; i++) 
			lambda_dav[i] = lambda_d[i] = Y_2[i] ;
	}

	//    gmma=Y_2./max(lambda_d,1e-10);
     //       eta=alpha_eta*eta_2term+(1-alpha_eta)*max(gamma-1,0);
     //   eta=max(eta,eta_min);
    //    v  = gmma .* eta ./ (1+eta);     
	for(i=0; i<FRAME_LEN21; i++){ 
		gamma[i] = Y_2[i] / max(lambda_d[i], 1e-10); 	    
		eta[i] = alpha_eta * eta_2term[i] + (1-alpha_eta) * max(gamma[i]-1,0);
		eta[i]=max(eta[i],eta_min);
        v[i]  = gamma[i] * eta[i] / (1+eta[i]);    
	}
 
     // Sf=conv(b,Y_2);  % smooth over frequency
     // Sf=Sf(w+1:M21+w);
	b_conv(Y_2, Sf, &size);  
	Sf_ptr = &Sf[1];      

      //  if FrameCnt==1     % new version omlsa3
      //         Sy=Y_2;
      //         S=Sf;
      //         St=Sf;
      //         lambda_d=Y_2;      
      //  else
      //         S=alpha_s*S+(1-alpha_s)*Sf;     % smooth over time      
      //  end    
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


 //       if l<14+l_fnz     % new version omlsa3
 //           Smin=S;
 //           SMact=S;
 //       else
 //           Smin=min(Smin,S);
 //           SMact=min(SMact,S);
 //       end

 			
     if(FrameCnt<15) { 
		 for(i=0; i<FRAME_LEN21; i++) {
			 Smin[i] = S[i]; 
			 Smact[i] = S[i];
		 }
	 }
     else{
		 for(i=0; i<FRAME_LEN21; i++)  {
			 Smin[i] = min(Smin[i], S[i]); 
			 Smact[i] =  min(Smact[i], S[i]); 
		 }
	 }

	 // Sft = St;
	 for(i=0; i<FRAME_LEN21; i++) Sft[i] = St[i];

 
   //     I_f=double(Ya2<delta_y*Bmin.*Smin & S<delta_s*Bmin.*Smin);
  //      conv_I=conv(b,I_f);
  //      conv_I=conv_I(w+1:M21+w);
  //      Sft=St;
  //      idx=find(conv_I);
  //      if ~isempty(idx)       
  //         conv_Y=conv(b,I_f.*Ya2);
  //         conv_Y=conv_Y(w+1:M21+w);
  //         Sft(idx)=conv_Y(idx)./conv_I(idx);
  //      end
	

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
		double Y2_I[FRAME_LEN21];
        double* Conv_Y_ptr;
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


 //   if FrameCnt<14+1     % new version omlsa3
 //       St=S;
 //       Smint=St;
 //       SMactt=St;
 //   else
 //       St=alpha_s*St+(1-alpha_s)*Sft;
 //       Smint=min(Smint,St);
 //       SMactt=min(SMactt,St);
 //   end
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
            Smint[i]=min(Smint[i],St[i]);
            Smactt[i]=min(Smactt[i],St[i]);
		 }
	 }

 //      qhat=ones(M21,1);    
 //      phat=zeros(M21,1);
 //      Sr_Y2 = Y_2./Bmin./max(Smint,1e-10);   % like the Sr for this      
 //      Sr_S = S./Bmin./max(Smint,1e-10);      % new version

 
        for(i=0;i<FRAME_LEN21;i++){         
		 	 Sr_Y2[i] = Y_2[i] / Bmin / max(Smint[i], 1e-10);
			 Sr_S[i] = S[i]/Bmin/max(Smint[i],1e-10); 
		}
//       idx=find(Sr_Y2>1 & Sr_Y2<delta_yt & Sr_S<delta_s);
//       qhat(idx)=(delta_yt-Sr_Y2(idx))/(delta_yt-1);
//       phat(idx)=1./(1+   qhat(idx)./(1-qhat(idx)).*(1+eta(idx)).*exp(-v(idx)));
//       phat(Sr_Y2>=delta_yt | Sr_S>=delta_s)=1; % both '>' , set 1
         for(i=0;i<FRAME_LEN21;i++){ 

			 if(i==40)
				 i = i;

			 if(Sr_Y2[i]>1 && Sr_Y2[i]<delta_yt & Sr_S[i]<delta_s){
				 qhat[i]= (delta_yt-Sr_Y2[i])/(delta_yt-1);
				// phat[i]=1 /((1+qhat[i])/(1-qhat[i])*(1+eta[i])*exp(-v[i]));
				  phat[i]=1 /(1+ (qhat[i]/(1-qhat[i]))*(1+eta[i])*exp(-v[i]));
			 }
			 else if(Sr_Y2[i]>=delta_yt || Sr_S[i]>=delta_s){
				 phat[i] = 1;
			 }
			 else{
				 qhat[i] = 1; phat[i] = 0;
			 }
		 }

 //     alpha_dt = alpha_d + (1-alpha_d) * phat;  
 //     lambda_d = alpha_dt .* lambda_d + (1- alpha_dt) .* Y_2;

         for(i=0;i<FRAME_LEN21;i++){ 
			 alpha_dt[i] = alpha_d + (1- alpha_d) * phat[i];  
             lambda_dav[i] = alpha_dt[i]* lambda_dav[i] + (1- alpha_dt[i]) * Y_2[i];
		 }

 //     if FrameCnt<15     % new version omlsa3
 //          lambda_d_long=lambda_d;
 //      else
 //          alpha_dt_long=alpha_d_long+(1-alpha_d_long)*phat;
 //          lambda_d_long=alpha_dt_long.*lambda_d_long+(1-alpha_dt_long).*Y_2[i];
 //     end   
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

 //    % update Smin Stmp S_his
 //        if l==Vwin-1+l_fnz    % new version omlsa3
 //            SW=repmat(S,1,Nwin);
 //            SWt=repmat(St,1,Nwin);
 //        else
 
 //            SWt=[SWt(:,2:Nwin) SMactt];
 //            Smint=min(SWt,[],2);
 //            SMactt=St;
 //        end             
 //   end
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

				  Smin[i] = min(S_his0[i][6], Smact[i]);
  				  Smin[i] = min(S_his0[i][5], Smin[i]);
   			      Smin[i] = min(S_his0[i][4], Smin[i]);
                  Smin[i] = min(S_his0[i][3], Smin[i]);
  				  Smin[i] = min(S_his0[i][2], Smin[i]);
   			      Smin[i] = min(S_his0[i][1], Smin[i]);
                  Smin[i] = min(S_his0[i][0], Smin[i]);

				  S_his0[i][7] = S_his0[i][6];
  				  S_his0[i][6] = S_his0[i][5];
				  S_his0[i][5] = S_his0[i][4];
				  S_his0[i][4] = S_his0[i][3];
				  S_his0[i][3] = S_his0[i][2];
  				  S_his0[i][2] = S_his0[i][1];
				  S_his0[i][1] = S_his0[i][0];
				  S_his0[i][0] =  Smact[i];
				   
                  Smin[i] = min(Smin[i], Smact[i]);
				  Smact[i] = S[i];

				  Smint[i] = min(S_his[i][6], Smactt[i]);
  				  Smint[i] = min(S_his[i][5], Smint[i]);
   			      Smint[i] = min(S_his[i][4], Smint[i]);
                  Smint[i] = min(S_his[i][3], Smint[i]);
  				  Smint[i] = min(S_his[i][2], Smint[i]);
   			      Smint[i] = min(S_his[i][1], Smint[i]);
                  Smint[i] = min(S_his[i][0], Smint[i]);

				  S_his[i][7] = S_his[i][6];
  				  S_his[i][6] = S_his[i][5];
				  S_his[i][5] = S_his[i][4];
				  S_his[i][4] = S_his[i][3];
				  S_his[i][3] = S_his[i][2];
  				  S_his[i][2] = S_his[i][1];
				  S_his[i][1] = S_his[i][0];
				  S_his[i][0] =  Smactt[i];
				  
                  Smint[i] = min(Smint[i], Smactt[i]);
				  Smactt[i] = St[i];
			  }

			   Smactt[0] =  Smactt[0];
		 }         
	 }

    //  lambda_d1=1.4685*lambda_d ;
     for(i=0;i<FRAME_LEN21;i++){
		 lambda_d[i] = 1.4685*lambda_dav[i];
		 lambda_d_global[i] = lambda_d[i];
	 }


	 // estimate speech absense probability 'q'  
 //   xi = beta .* xi + (1-beta) .* eta;       % eq. (23)           
 //   xi_local=conv(xi,b_xi_local);            % eq. (24)      
 //   xi_local=xi_local(w_xi_local+1:M21+w_xi_local);
 //   xi_global=conv(xi,b_xi_global);          % eq. (24)      
 //   xi_global=xi_global(w_xi_global+1:M21+w_xi_global);
    


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
 
 //   if xi_local>0,  xi_local_dB=10*log10(xi_local)  ; else xi_local_dB=-100;  end
 //   if xi_global>0, xi_global_dB=10*log10(xi_global); else xi_global_dB=-100; end
 //   if xi_frame>0,  xi_frame_dB=10*log10(xi_frame)  ; else xi_frame_dB=-100;  end

      for(i=0;i<FRAME_LEN21;i++){
		  if(xi_local[i]>0){
			  xi_local_dB[i] = 10*log10(xi_local_ptr[i]);
		  }
		  else{
			  xi_local_dB[i] = -100;
		  }
		  if(xi_global[i]>0){
			  xi_global_dB[i] = 10*log10(xi_global_ptr[i]);
		  }
		  else{
			  xi_global_dB[i] = -100;
		  }
	  }
	  if(xi_frame>0){
  	    xi_frame_dB = 10*log10(xi_frame);
	  }
	  else {
	    xi_frame_dB = -100;
	  }

//     P_local = ones(M21,1);
//     P_local(xi_local_dB<=xi_min_dB) = P_min;% use 0.005 instead of 0, costrain the value
//     idx=find(xi_local_dB>xi_min_dB & xi_local_dB<xi_max_dB);
//     P_local(idx) = P_min+(xi_local_dB(idx)-xi_min_dB)/(xi_max_dB-xi_min_dB)*(1-P_min); 
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
//        m_P_local=mean(P_local(3:(k2_local+k3_local-3)));    % average probability of speech presence 500 ~3500hz

        

        m_P_local = 0;
        for(i=2;i<k2_local+k3_local-3;i++){    //        
		   m_P_local += P_local[i];
		}
		  
		m_P_local =  m_P_local/(k2_local+k3_local-3-2);
//    
//        if m_P_local<0.25
//            P_local(k2_local:k3_local)=P_min;    % reset P_local (frequency>500Hz) for low probability of speech presence
//        end         
//    
        if(m_P_local<0.25){
			for(i=k2_local;i<k3_local;i++){
				P_local[i]=P_min;
			}
		}

	if(FrameCnt==200)
			FrameCnt= FrameCnt;

//       if tone_flag               % new version
//           if (m_P_local<0.5) && (FrameCnt>120)
//               idx=find( alpha_dt_long(8:(M21-8)) > 2.5*(alpha_dt_long(10:(M21-6))+alpha_dt_long(6:(M21-10))) );
//               P_local([idx+6;idx+7;idx+8])=P_min;   % remove interfering tonals
//           end
//       end 

         if ((m_P_local<0.5) && (FrameCnt>120)){
             for(i=8;i<FRAME_LEN21-8;i++){			 
				//  idx=find( lambda_dav_long(8:(M21-8)) > 2.5*(lambda_dav_long(10:(M21-6))+lambda_dav_long(6:(M21-10))) ){
				 if(  lambda_d_long[i] > 2.5*(lambda_d_long[i] + lambda_d_long[i-2])  ){
					    P_local[i+6] = P_local[i+7]= P_local[i+8] = P_min;
				 }
			 }
		 }
 
//     % P_global
//     P_global = ones(M21,1);
//     P_global(xi_global_dB <= xi_min_dB) = P_min;
//     idx = find(xi_global_dB > xi_min_dB & xi_global_dB < xi_max_dB);
//     P_global(idx) = P_min + (xi_global_dB(idx) - xi_min_dB)/ (xi_max_dB - xi_min_dB) * (1-P_min);
//     
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
		 
//     
//     % P_frame  Fig.3  eq.26 & eq.27
//     if xi_frame_dB<=xi_min_dB     
//         P_frame=P_min;        
//     elseif xi_frame >= ex_xi_frame             
//         xi_peak_dB=min(max(xi_frame_dB,xi_p_min_dB),xi_p_max_dB);            
//         P_frame=1;        
//     elseif xi_frame_dB>=xi_peak_dB+xi_max_dB % eq.(27)  u = 1  xi_frame > xi_peak * xi_max
//         P_frame=1;        
//     elseif xi_frame_dB<=xi_peak_dB+xi_min_dB  % eq.(27) u = 0  xi_frame < xi_peak * xi_min         
//         P_frame=P_min;        
//     else
//          %  P_frame= (xi_frame_dB-xi_min_dB-xi_peak_dB)/(xi_max_dB-xi_min_dB);             
//          P_frame=P_min+(xi_frame_dB-xi_min_dB-xi_peak_dB)/(xi_max_dB-xi_min_dB)*(1-P_min);         
//     end  
       
       if( xi_frame_dB<=xi_min_dB){
            P_frame=P_min; 
	   }   
	   else if(xi_frame >= ex_xi_frame) {
		    xi_peak_dB=min(max(xi_frame_dB,xi_p_min_dB),xi_p_max_dB);  
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

	

         
//      
//      
//    q = 1 - P_frame .* P_global .* P_local;     
//    q = min(qmax, q);
      for(i=0;i<FRAME_LEN21;i++){

		  q[i] = 1 - P_frame*P_global[i]*P_local[i];
		  q[i] = min(qmax,q[i]);
	  }
//     
//    % (2) est speech present probability 'P' and gain
//
//        gmma=X_F_2./max(lambda_d1,1e-10);                       % postier SNR eq.10 in paper
//        eta=alpha_eta*eta_2term +(1-alpha_eta)*max(gmma-1,0); % proior SNR  eq.28  
//        eta=max(eta,eta_min);
//        v=gmma.*eta./(1+eta);                                % eq.10       
//        
      for(i=0;i<FRAME_LEN21;i++){
          gamma[i] = Y_2[i] / max(lambda_d_global[i], 1e-10);
		  eta[i] = alpha_eta*eta_2term[i] + (1-alpha_eta)*max(gamma[i]-1,0);
		  eta[i] = max(eta[i],eta_min);
		  v[i]=gamma[i]*eta[i]/(1+eta[i]);    
	  }



//     
//     PH1=zeros(M21,1);
//     idx=find(q<0.9);
//     PH1(idx)=1./(1+  q(idx)./(1-q(idx)).*(1+eta(idx)).*exp(-v(idx)));
        for(i=0;i<FRAME_LEN21;i++){
           if(q[i]<0.9){
               PH1[i] = 1/( 1 +(q[i] /(1 - q[i])) * (1+eta[i]) * exp(-v[i])  );
		   }
		   else{
			   PH1[i] = 0;
		   }            
		}
//    
//    % Gain 
//     GH1=ones(M21,1);
//     idx=find(v>5);
//     GH1(idx) = (eta(idx)./(1+eta(idx)));
//     idx=find(v<=5 & v>0);
//     GH1(idx)=eta(idx)./(1+eta(idx)).*exp(0.5*expint(v(idx)));
       for(i=0;i<FRAME_LEN21;i++){
		   if(v[i]>5){
			   GH1[i] = eta[i]/(1+eta[i]);
		   }
		   else if(v[i]<=5 && v[i]>0){
                GH1[i] = (eta[i]/(eta[i]+1))*exp(0.5*expint(v[i])); // need
		   }	
		   else {
		        GH1[i] = 1;
		   }


		   GH1[i]*=0.5;
	   }

//  
//    % G = GH1.^P .* (Gmin.^(1-P));
//         
//     if tone_flag   % new version     
//         lambda_d_global=lambda_d1;   % new version         
//         lambda_d_global(4:M21-3)
//       =min([lambda_d_global(4:M21-3),  //   3~M21-4
//             lambda_d_global(1:M21-6),  //   0~M21-7
//             lambda_d_global(7:M21  )]  //   6~M21-1
 //            ,[],2);   % new version         
//         Sy=0.8*Sy+0.2*X_F_2;    % new version
//         GH0=Gmin*(lambda_d_global./(Sy+1e-10)).^0.5;   % new version omlsa3        
//     else   % new version            
//         GH0=Gmin;   %#ok<UNRCH> % new version   
//     end   % new version
//        
//     G=GH1.^PH1.*GH0.^(1-PH1);     
//     eta_2term=GH1.^2.*gmma;
//       
        for(i=3;i<FRAME_LEN21-3;i++){
			double minval;

		 

			minval = min(lambda_d[i] , lambda_d[i-3] );
            lambda_d_global[i] = min(minval , lambda_d[i+3] );
		}
        
		for(i=0;i<FRAME_LEN21;i++){
			Sy[i]=0.8*Sy[i]+0.2*Y_2[i];    
		}
		
         for(i=0;i<FRAME_LEN21;i++){
            GH0[i]= Gmin * sqrt(lambda_d_global[i]/(Sy[i]+ 1e-10));         			
		 }

		  for(i=0;i<FRAME_LEN21;i++){      
		    G[i]= pow(GH1[i],PH1[i]) * pow(GH0[i] ,(1-PH1[i]));  
		  }
       
         for(i=0;i<FRAME_LEN21;i++){
             eta_2term[i]=GH1[i]*GH1[i]*gamma[i];
		 }
     

 //   X_G(1:M21) = X_F(1:M21).* G;        
 //   X_G(M21+1:FFT_LEN) = conj(X_G(M21-1:-1:2));
 //  
 

        
        for(i=0;i<FRAME_LEN21;i++){
			X_G_real[i] = Y_real[i] * G[i];
			X_G_imag[i] = Y_imag[i] * G[i];
		}  
        for(i=FRAME_LEN21;i<FRAME_LEN-1;i++){
			X_G_real[i] = X_G_real[FRAME_LEN-i];
			X_G_imag[i] = -X_G_imag[FRAME_LEN-i];
		}

	

		
#endif

        
     //   X=[zeros(3,1); G(4:M21-1).*Y(4:M21-1); 0];
     //   X(M21+1:M)=conj(X(M21-1:-1:2)); %extend the anti-symmetric range of the spectum
     //   x=Cwin^2*win.*real(ifft(X));


        for(i=0;i<FRAME_LEN21;i++){
		
			if(i<3 ||i==(FRAME_LEN21-1)){
			   X_COMPLEX[i].re = 0;
  			   X_COMPLEX[i].im = 0;
			}
			else {
 			   X_COMPLEX[i].re = Y_COMPLEX[i].re * G[i];
			   X_COMPLEX[i].im = Y_COMPLEX[i].im * G[i];			
			}
		}  
        for(i=FRAME_LEN21;i<FRAME_LEN ;i++){
			X_COMPLEX[i].re = X_COMPLEX[FRAME_LEN-i].re;
			X_COMPLEX[i].im = -X_COMPLEX[FRAME_LEN-i].im;
		}
		 
		ifft(X_COMPLEX, x_g_c  );


        for(i=0;i<FRAME_LEN;i++){
		//	x_ifft[i]=x_ifft[i]*Window[i];
			x_ifft[i]=  x_g_c[i].re   *Window[i]* Cwin*Cwin;
		}



		for(i=0;i<FRAME_LEN;i++){
			out_buf[i] = out_buf[i] + x_ifft[i];
		}

		for(i=0;i<FRAME_LEN/4;i++){
			int temp;

			temp = ( out_buf[i]* 32768 );

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
   double W0[FRAME_LEN41]; 
   double W0_avg;//, Cwin;
 
    for(i=0;i<len;i++){
        Window[i] = 0.54-0.46*cos(2*M_PI*(i )/(len-1));
		Win2[i] = Window[i]*Window[i];		 
	}
   
     W0_avg = 0;
     for(i=0;i<FRAME_LEN41;i++){
       W0[i] = Win2[i] + Win2[i+len/4] + Win2[i+len/2] + Win2[i+len*3/4];
	   W0_avg +=  W0[i] ;
     }

	 W0_avg /= FRAME_LEN41;

	 W0_avg= sqrt(W0_avg);


     Cwin = 0;
    for(i=0;i<FRAME_LEN;i++){
		Window[i]/= W0_avg;
		Cwin+= Window[i]*Window[i];
	}

	Cwin = sqrt(Cwin);
  
    for(i=0;i<FRAME_LEN;i++){
		Window[i]/= Cwin;	 
	}

   
}





 