
#include<string.h>
#include<stdio.h>
#include<stdlib.h>
 #include "signal_processing_library.h"
 #include "gain_control.h"

short Outagc[320];
  short Input160[320];
  int len;
  int audio_frame;
  int frame_size;
  unsigned char EncOutput[160];
  int max_packet;
    unsigned char saturationWarning;
	int micLevelIn = 0;
	int micLevelOut = 0;
	int inMicLevel= 0;//  = micLevelOut;
	int outMicLevel = 0;

	int nRet = 0;
	int nAgcRet;
	short *pData    = NULL;
	short *pOutData = NULL;
	void *agcHandle = NULL;



	int frameSize;

	extern void* WnAgcInit();

	void* agcHandle ;

unsigned char wavheader[44];

void main(){


	FILE* finput, *foutput;
     char namein[]= {".wav"};
     char nameout[]={"_agc_new.wav"};
	 char name[100]  ={"../voice/omtest3omlsa_lite_agc_803"};// {"../ysy/2018-11-16T09-53-45_10uF,4dBm"};  

     char nameo[100];

     memcpy(nameo,name, sizeof(name));

     strcat(name,namein); 
     finput = fopen(name ,"rb");



	 if(finput==NULL){
		 printf("No such file, check again\n");
		 system("pause");

		 return;
	 }
    
     strcat(nameo,nameout)  ;
     foutput = fopen(nameo,"wb");

		fread( wavheader ,0x2c,1,finput);
		fwrite(wavheader,0x2c,1,foutput);


	agcHandle = 	WnAgcInit();


	while(!(feof(finput))){
	
	fread( Input160 ,640,1,finput);
	   //memcpy(Input160, alltab, 320   );
	nAgcRet = WebRtcAgc_Process(agcHandle, Input160, 0, 320, Outagc,NULL, inMicLevel, &outMicLevel, 0, &saturationWarning);
	 

 printf("outMicLevel is %d \n", outMicLevel);

	fwrite(Outagc,640,1,foutput);
	}
	fclose(foutput);

	 system("pause");
}