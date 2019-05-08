
#include"libomlsa_postfilter.h"
#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#define INLEN   320
#define INLEN41 80


short pcminput[INLEN41];

int Inlen;

short pcmoutput[INLEN41];

int outlen;

int maincnt = 0;


 
void main(){
 
    FILE*fin, *fout;
     char namein[]= {".wav"};
     char nameout[]={"omlsaout_varlen.wav"};
     char name[100]={"../mi/1"};  

     char nameo[100];


//	 printf("Enter the name of input file: Notice!!! if the input file is 'test.wav', just enter 'test'\n");

//	 scanf("%s", name);
	  
     memcpy(nameo,name, sizeof(name));

     strcat(name,namein); 
     fin = fopen(name ,"rb");

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

		// if(!maincnt) {		 	
		//	 short tempbuf[512];
		//
		//	 fread(tempbuf, INLEN*2, 1, fin); 
        //     PostFilterProcess(tempbuf, INLEN,pcmoutput, &outlen);
		// }
		// else
		 {
		     fread(pcminput, INLEN41*2, 1, fin); 
             PostFilterProcess(pcminput, INLEN41,pcmoutput, &outlen);
		 }

	     fwrite(pcmoutput,INLEN41*2,1,fout);

		 maincnt = 1;
	 }


	 printf("End of process!\n");

	 system("pause");

}