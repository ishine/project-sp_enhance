
#include"posterfilter\libomlsa_postfilter.h"
#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#define INLEN   512
#define INLEN41 128


short pcminput[INLEN41];

int Inlen;

short pcmoutput[INLEN41];

int outlen;

int maincnt = 0;

extern float Window[512]; 
 
void main(){
 
    FILE*fin, *fout,  *ftest;
     char namein[]= {".wav"};
     char nameout[]={"omlsaout_opt.wav"};
	 char name[100]  ={"../voice/omtest2"};// {"../ysy/2018-11-16T09-53-45_10uF,4dBm"};  

     char nameo[100];


	 printf("Enter the name of input file: Notice!!! if the input file is 'test.wav', just enter 'test'\n");

	// scanf("%s", name);
	  
     memcpy(nameo,name, sizeof(name));

     strcat(name,namein); 
     fin = fopen(name ,"rb");
	   
	 ftest = fopen("test.bin" ,"wb");


	 if(fin==NULL){
		 printf("No such file, check again\n");
		 system("pause");

		 return;
	 }
    
     strcat(nameo,nameout)  ;
     fout = fopen(nameo,"wb");

	 	 
	 fread(pcminput, 0x2c, 1, fin); 
	 fwrite(pcminput,0x2c,1,fout);

	 PostFilterInit();

	  printf("Start processing!\n");

     while(feof(fin)==NULL){
		 int i;
	 
		 {
		     fread(pcminput, INLEN41*2, 1, fin); 
             PostFilterProcess(pcminput, INLEN41,pcmoutput, &outlen);
		 }

		// if(maincnt==0)
		// for(i=0; i < 512 ; i++){
		// 
		//   fprintf(ftest, "%8ff, " , Window[i] );
		//
		//
		//   if( !((i+1)%16))
		//	   fprintf(ftest, "\n"  );
		//
		// 
		// }




	     fwrite(pcmoutput,INLEN41*2,1,fout);

		 maincnt = 1;
	 }


	 printf("End of process!\n");

	 system("pause");

}