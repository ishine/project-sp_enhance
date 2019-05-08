#include<math.h>
#include"omlsa_def.h"
 
 
static int FFT_CAL(double* input, double*fftOutReal, double*fftOutImg);
  int IFFT_CAL(double* input, double*fft_nextReal );

static double W_FFT[2*FFT_LEN];
static double IW_FFT[2*FFT_LEN];
double HannWin[FRAME_LEN];
extern double Window[FRAME_LEN];

void HannWinMake(double*win, int len){
	int i;
	for(i=0;i< len;i++){
		win[i] =  0.5 *(1 - (cos(2*M_PI*(i)/len)));
	}
}

void HanningWin(int len){
	int i;

	for(i=0;i<len;i++){	
		HannWin[i] = 0.5 *(1 - (cos(2*M_PI*(i)/len)));
	}
}
 
void FftWin(){
	int n;

    for(n=0;n<FFT_LEN ;n++)
	{
		IW_FFT[n] = W_FFT[n]  =     cos( n *2*   M_PI / FFT_LEN);  //
		W_FFT[n + FFT_LEN  ] =  sin( n * 2* M_PI / FFT_LEN); 	  // 
		IW_FFT[n + FFT_LEN] =  - W_FFT[n + FFT_LEN  ] ;		
	}
}

void D_M_ANC_WIN_FFT(double* pX , int len, double*pReaOut,  double*pImgOut){

	int i;

   double FftInputRealX[FFT_LEN*2]= {0};
 
	// hanning window
	for(i=0; i < len; i++){
	   FftInputRealX[i]  = pX[i] * Window[i];  	 	  
	}	  

	// fft operation
	FFT_CAL(FftInputRealX, pReaOut,pImgOut);	 
}
 
  int IFFT_CAL(double* input_real, double*input_imag, double*fft_nextReal )
{
	int k,temp,n, temp1,j,m,i;
	double a,b, c;
	int  fft_list[(FFT_LEN)][M]; 
	short  fft_list1[(FFT_LEN)]; 
 	double fft_input[2*(FFT_LEN)] = {0};
	double fft0[2*(FFT_LEN)] = {0};
	double fft_temp[2*FFT_LEN] = {0};
 
	for (k= 0; k < (FFT_LEN); k++)
	{
			 temp = k;
			 fft_list1[k] = 0;
		
			 for (n= 0; n < M; n++)
			 {
				fft_list[k][n] = (temp&0x1);
				temp>>=1;
				fft_list1[k] += fft_list[k][n] << (M-1-n);
			 }
	}
	 
		for (n = 0; n < (FFT_LEN); n ++)
		{
 
			fft0[n] = (double)input_real[fft_list1[n]] ; 
			fft0[n + FFT_LEN] = (double)input_imag[fft_list1[n]] ; 
		}
		 

		 for( n = 0; n< (FFT_LEN); n++)
		 {
			fft_input[n] =  fft0[n];
			fft_input[n+ FFT_LEN] =  fft0[n+ FFT_LEN];
		 }
		 


		for(i = 1; i < (FFT_LEN) ; i <<=1)
		{	 
			int m; m = (FFT_LEN)/(i*2);		

			for (k = 0; k < (FFT_LEN); k+=(i*2))  
			{
				for (n = k; n< k+ i  ;n++)  
				{
				a = ( fft_input[n + i  + (FFT_LEN)] * IW_FFT[(FFT_LEN)+ m*(n - k)]);  		
				b  =  ( fft_input[n +  i  ] * IW_FFT[ m*(n - k)] );		 
				b  = a + b;  
				c =   fft_input[n]  ;
				fft_temp[n]	= c + b;  fft_temp[n + i ] = c - b;

				a = ( fft_input[n + i + (FFT_LEN)] * IW_FFT[ m *(n - k)]);		 
				b =  ( fft_input[n + i] * IW_FFT[(FFT_LEN) + m * (n - k)]);	 
				b = a - b;  
				c =  fft_input[n + (FFT_LEN)]  ;	 
				fft_temp[n + (FFT_LEN)] = c + b;
		 		fft_temp[n + i + (FFT_LEN)] = c - b;
				}
			}
             for (n = 0; n< 2*(FFT_LEN); n++)
             {
				fft_input[ n ] = fft_temp[n];
             }
		}


	  for (n = 0; n<   FFT_LEN  ; n++){
		fft_nextReal[n] = 	(fft_temp[n] / FFT_LEN );
		//fft_nextImg[n]  =(int)	(fft_temp[n+FFT_LEN] / FFT_LEN );
	  }

		return 0;
} 

static int FFT_CAL(double* input, double*fftOutReal, double*fftOutImg)
{
	int k,temp,n, temp1,j,m,i;
	double a,b, c;
	int  fft_list[(FFT_LEN)][M]; 
	short  fft_list1[(FFT_LEN)]; 
 	double fft_input[2*(FFT_LEN)] = {0};
	double fft0[2*(FFT_LEN)] = {0};
	double fft_temp[2*FFT_LEN] = {0};

	// �����λ����
	for (k= 0; k < (FFT_LEN); k++)
	{
			 temp = k;
			 fft_list1[k] = 0;
		
			 for (n= 0; n < M; n++)
			 {
				fft_list[k][n] = (temp&0x1);
				temp>>=1;
				fft_list1[k] += fft_list[k][n] << (M-1-n);
			 }
	}
	 
		for (n = 0; n < (FFT_LEN); n ++)
		{
 
			fft0[n] = (double)input[fft_list1[n]] ; 
			fft0[n + FFT_LEN] = (double)input[fft_list1[n] + FFT_LEN] ; 
		}
		 

		 for( n = 0; n< (FFT_LEN); n++)
		 {
			fft_input[n] =  fft0[n];
			fft_input[n+ FFT_LEN] =  fft0[n+ FFT_LEN];
		 }
		  
		for(i = 1; i < (FFT_LEN) ; i <<=1)
		{	 
			int m; m = (FFT_LEN)/(i*2);		

			for (k = 0; k < (FFT_LEN); k+=(i*2))  
			{
				for (n = k; n< k+ i  ;n++)  
				{
				a = ( fft_input[n + i  + (FFT_LEN)] * W_FFT[(FFT_LEN)+ m*(n - k)]);  		
				b  =  ( fft_input[n +  i  ] * W_FFT[ m*(n - k)] );		 
				b  = a + b;   
				c =   fft_input[n]  ;
				fft_temp[n]	= c + b;  fft_temp[n + i ] = c - b;

				a = ( fft_input[n + i + (FFT_LEN)] * W_FFT[ m *(n - k)]);		 
				b =  ( fft_input[n + i] * W_FFT[(FFT_LEN) + m * (n - k)]);	 
				b = a - b;   
				c =  fft_input[n + (FFT_LEN)]  ;	 
				fft_temp[n + (FFT_LEN)] = c + b;
		 		fft_temp[n + i + (FFT_LEN)] = c - b;
				}
			}
             for (n = 0; n< 2*(FFT_LEN); n++)
             {
				fft_input[ n ] = fft_temp[n];
             }
		}
		 
	  for (n = 0; n < FFT_LEN; n++){
		fftOutReal[n] =(double)	fft_temp[n];
		fftOutImg[n]  =(double)	fft_temp[n+ FFT_LEN];
	  }

   	 return 0;
} 
