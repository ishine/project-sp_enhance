#include"AGC\libWnDigitalAGC.h"
#include"posterfilter\libomlsa_postfilter.h"
#include<string.h>
#include<stdio.h>
#define INLEN   512
#define INLEN41 128


short pcminput[INLEN41];

int Inlen;

short pcmoutput[INLEN41];

int outlen;

int maincnt = 0;


// for agc 
int Agc_ready =0;
int Revagc,Outlen ;
short InputAgc[640];
short OutAgc[320];


void main(){
 
    FILE*fin, *fout,*foutagc;
     char namein[]= {".wav"};
     char nameout[]={"omlsa_opt.wav"};
	 char nameagcout[]={"omlsa_opt_agc.wav"};
     char name[100] = {"../voice/2018-11-16T13-14-31_10uF"};  

     char nameo[100];
	 char nameoagc[100];
	 WnAgcInit(15, 0); // agc init
	  
     memcpy(nameo,name, sizeof(name));
     memcpy(nameoagc,name, sizeof(name));
     strcat(name,namein); 
     fin = fopen(name ,"rb");
    
     strcat(nameo,nameout)  ;
     fout = fopen(nameo,"wb");

     strcat(nameoagc,nameagcout)  ;
     foutagc = fopen(nameoagc,"wb");


	 	 
	 fread(pcminput, 0x2c, 1, fin); 
	 fwrite(pcminput,0x2c,1,fout);
	 fwrite(pcminput,0x2c,1,foutagc);
	 PostFilterInit();

     while(feof(fin)==NULL){  
		fread(pcminput, INLEN41*2, 1, fin); 
        PostFilterProcess(pcminput, INLEN41,pcmoutput, &outlen);		 
	    fwrite(pcmoutput,INLEN41*2,1,fout);

		memcpy(&InputAgc[maincnt*128], pcmoutput ,  INLEN41*2);

		maincnt += 1;

		if(maincnt==5){ // 640 
			
   		    maincnt = 0;

		    Revagc = WnAgcProcess(InputAgc, 320, OutAgc, &Outlen);
			fwrite(OutAgc,320*2,1,foutagc);

			Revagc = WnAgcProcess(&InputAgc[320], 320, OutAgc, &Outlen);
			fwrite(OutAgc,320*2,1,foutagc);

		}

	 }


	 fclose(fin);
	 fclose(fout);
	 fclose(foutagc);

}