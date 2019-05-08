
#include"libomlsa_postfilter.h"
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

int posterfilterinit = 0;


short postertestdata[128] = {

 0xffff,0x30fb,0x5a82,0x7641,0x7fff,0x7641,0x5a83,0x30fb,0x0,0xcf05,0xa57e,0x89c0,0x8001,0x89bf,0xa57f,0xcf05,0x0,0x30fc,0x5a82,0x7641,0x7fff,0x7641,0x5a82,0x30fc,0x0,0xcf04,0xa57e,0x89c0,0x8001,0x89bf,0xa57e,0xcf05,0xffff,0x30fb,0x5a81,0x7641,0x7fff,0x7641,0x5a82,0x30fc,0x0,0xcf05,0xa57e,0x89bf,0x8001,0x89c0,0xa57e,0xcf05,0x0,0x30fc,0x5a82,0x7641,0x7fff,0x7640,0x5a82,0x30fc,0xffff,0xcf05,0xa57f,0x89bf,0x8001,0x89bf,0xa57e,0xcf05,0x0,0x30fc,0x5a82,0x7641,0x7fff,0x7641,0x5a81,0x30fc,0x0,0xcf05,0xa57e,0x89bf,0x8001,0x89bf,0xa57e,0xcf04,0x0,0x30fb,0x5a82,0x7641,0x7fff,0x7641,0x5a82,0x30fb,0x0,0xcf04,0xa57e,0x89bf,0x8001,0x89bf,0xa57e,0xcf05,0xffff,0x30fc,0x5a81,0x7641,0x7fff,0x7640,0x5a81,0x30fb,0x0,0xcf04,0xa57e,0x89c0,0x8001,0x89c0,0xa57e,0xcf05,0xffff,0x30fc,0x5a81,0x7641,0x7fff,0x7640,0x5a81,0x30fb,0x0,0xcf04,0xa57e,0x89c0,0x8001,0x89c0,0xa57e,0xcf05,
};
 
 
void posterfiltertest(){
	
	  if(posterfilterinit==0){
	      PostFilterInit();
		}
		
		posterfilterinit = 1;
		 
		memcpy( pcminput, postertestdata,  256 );
   
        PostFilterProcess(pcminput, INLEN41,pcmoutput, &outlen);
	  
}