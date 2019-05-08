
#include"postfilter\libomlsa_postfilter.h"
#include<string.h>
#include<stdio.h>
#define INLEN   512
#define INLEN41 128


short pcminput[INLEN41];

int Inlen;

short pcmoutput[INLEN41];

int outlen;

int maincnt = 0;


extern double expint(double x);
double val;
void main(){
 
    FILE*fin, *fout;
     char namein[]= {".wav"};
     char nameout[]={"omlsaout.wav"};
     char name[100] = {"../voice/omtest2"};  

     char nameo[100];

	 val = expint(0.3318);

     memcpy(nameo,name, sizeof(name));

     strcat(name,namein); 
     fin = fopen(name ,"rb");
    
     strcat(nameo,nameout)  ;
     fout = fopen(nameo,"wb");

	 	 
	 fread(pcminput, 0x2c, 1, fin); 
	 fwrite(pcminput,0x2c,1,fout);

	 PostFilterInit();

     while(feof(fin)==NULL){

		 if(!maincnt) {		 	
			 short tempbuf[512];

			 fread(tempbuf, INLEN*2, 1, fin); 
             PostFilterProcess(tempbuf, INLEN,pcmoutput, &outlen);
		 }
		 else{
		     fread(pcminput, INLEN41*2, 1, fin); 
             PostFilterProcess(pcminput, INLEN41,pcmoutput, &outlen);
		 }

	     fwrite(pcmoutput,INLEN41*2,1,fout);

		 maincnt = 1;
	 }


}