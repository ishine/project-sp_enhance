
#include"posterfilter\libomlsa_postfilter.h"
#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#define INLEN   512
#define INLEN41 128
#include <math.h>

#include "agc\libWnDigitalAGC.h"

short pcminput[INLEN41];
short pcminputR[INLEN41];
int Inlen;

short pcmoutput[INLEN41];

int outlen;

int maincnt = 0;

//extern float Window[512]; 


short postertestdata[128] = {

 0xffff,0x30fb,0x5a82,0x7641,0x7fff,0x7641,0x5a83,0x30fb,0x0,0xcf05,0xa57e,0x89c0,0x8001,0x89bf,0xa57f,0xcf05,0x0,0x30fc,0x5a82,0x7641,0x7fff,0x7641,0x5a82,0x30fc,0x0,0xcf04,0xa57e,0x89c0,0x8001,0x89bf,0xa57e,0xcf05,0xffff,0x30fb,0x5a81,0x7641,0x7fff,0x7641,0x5a82,0x30fc,0x0,0xcf05,0xa57e,0x89bf,0x8001,0x89c0,0xa57e,0xcf05,0x0,0x30fc,0x5a82,0x7641,0x7fff,0x7640,0x5a82,0x30fc,0xffff,0xcf05,0xa57f,0x89bf,0x8001,0x89bf,0xa57e,0xcf05,0x0,0x30fc,0x5a82,0x7641,0x7fff,0x7641,0x5a81,0x30fc,0x0,0xcf05,0xa57e,0x89bf,0x8001,0x89bf,0xa57e,0xcf04,0x0,0x30fb,0x5a82,0x7641,0x7fff,0x7641,0x5a82,0x30fb,0x0,0xcf04,0xa57e,0x89bf,0x8001,0x89bf,0xa57e,0xcf05,0xffff,0x30fc,0x5a81,0x7641,0x7fff,0x7640,0x5a81,0x30fb,0x0,0xcf04,0xa57e,0x89c0,0x8001,0x89c0,0xa57e,0xcf05,0xffff,0x30fc,0x5a81,0x7641,0x7fff,0x7640,0x5a81,0x30fb,0x0,0xcf04,0xa57e,0x89c0,0x8001,0x89c0,0xa57e,0xcf05,
};
 

#define loge10val  0.43429447587231051f
void main(){
 
      FILE*fin,*finR, *fout,  *ftest;
     char namein[]= {".wav"};
     char nameout[]={"omlsa_lite_fixagc.wav"};
	 char name[100]  ={"../voice3/t192l"};// {"../ysy/2018-11-16T09-53-45_10uF,4dBm"};  
	 char nameR[100]  ={"../voice3/t192r"};
     char nameo[100];
	  
	 float x,y,z, w; 

	 int i;
	 	 ftest = fopen("test.bin" ,"wb");
	   


	 printf("Enter the name of input file: Notice!!! if the input file is 'test.wav', just enter 'test'\n");

	// scanf("%s", name);
	  
     memcpy(nameo,name, sizeof(name));

     strcat(name,namein); 
     fin = fopen(name ,"rb");

     strcat(nameR,namein); 
     finR = fopen(nameR ,"rb");

	 if(fin==NULL){
		 printf("No such file, check again\n");
		 system("pause");

		 return;
	 }
    
     strcat(nameo,nameout)  ;
     fout = fopen(nameo,"wb");

	 	 
	 fread(pcminput, 0x2c, 1, fin); 
	// fwrite(pcminput,0x2c,1,fout);

	 fread(pcminputR, 0x2c, 1, finR); 
	 

	 PostFilterInit();
	   WnAgcInit(30, 11);
	  printf("Start processing!\n");

     while(feof(fin)==NULL){
		 int i;
	  
		     fread(pcminput, INLEN41*2, 1, fin); 
		     fread(pcminputR, INLEN41*2, 1, finR);    
	     //  memcpy(pcminputR,postertestdata, 160);

            // PostFilterProcess(pcminput, INLEN41,pcmoutput, &outlen);
		    Dual_Channel_PostFilterProcess(pcminput,pcminputR ,INLEN41,pcmoutput, &outlen);
		//    WnAgcProcess(pcminput, 128, pcmoutput, &outlen  ); 


	         fwrite(pcmoutput,INLEN41*2,1,fout);

		 maincnt = 1;
	 }


	 printf("End of process!\n");

	// system("pause");

}