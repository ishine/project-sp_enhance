 

#include "HighPass300hzFilter.h"
#include <math.h>
#include<string.h>
 
typedef   double real64_T;

#define FLOAFVER 1

#define ORDER 3


 #define MWSPT_NSEC 7
typedef double real64_T;
 
const static int NL[MWSPT_NSEC] = { 1,3,1,3,1,2,1 };
const static real64_T NUM[MWSPT_NSEC][3] = {
  {
     0.9616076266275,                 0,                 0 
  },
  {
                   1,                -2,                 1 
  },
  {
     0.9100023323578,                 0,                 0 
  },
  {
                   1,                -2,                 1 
  },
  {
     0.9443110625264,                 0,                 0 
  },
  {
                   1,                -1,                 0 
  },
  {
                   1,                 0,                 0 
  }
};
const static int DL[MWSPT_NSEC] = { 1,3,1,3,1,2,1 };
const static real64_T DEN[MWSPT_NSEC][3] = {
  {
                   1,                 0,                 0 
  },
  {
                   1,   -1.916526647419,   0.9299038590912 
  },
  {
                   1,                 0,                 0 
  },
  {
                   1,   -1.813675007231,      0.8263343222 
  },
  {
                   1,                 0,                 0 
  },
  {
                   1,  -0.8886221250528,                 0 
  },
  {
                   1,                 0,                 0 
  }
};


  static  float InterValue400[MWSPT_NSEC+1][ORDER] = {0};
  
  

  void TwoPointFilter400hz(int index, float* in , float *out)
  {
	  int j;
 
	   for(j = ORDER-1; j>0; j--)
	  	 out[j] = out[j-1];   
 
	   out[0] = NUM[index][0] * in[0] ;

	   for(j = 1; j < ORDER; j++ )
	   {
		   out[0]+=  NUM[index][j] * in[j] -  DEN[index][j] *  out[j];
	   } 

  }
   
   
  void WnHighPassFilterSingle400hz(short x, short*y )
  {
      long long temp ;
	  
	  int i, j; 
	  
	  for(j = 2; j>0; j--)
	  {
	  
	  InterValue400[0][j] = InterValue400[0][j-1];		   
	  	   
	  }
	  
	  InterValue400[0][0] = x;
 
	  //new sample for output	  
	  for( i = 0; i < MWSPT_NSEC-1; i++ )
	  {
		  TwoPointFilter400hz(i, InterValue400[i], InterValue400[i+1]);
	      
	  }
  
	  temp = (long long )InterValue400[MWSPT_NSEC - 1][0]  ;
  
	  
	  temp = (temp > 32767) ? 32767  : (temp < -32767 ? -32767 :  temp);
	  *y =  (short  )temp;
  }

 
  static int mutecnt400;    
    
	 
   void WnHighPassFilter400hz(short* In, short *Out , int length, int mutelen, int VoiceStartFlag)
  {
    int i ;

	if(VoiceStartFlag){
	   mutecnt400 = 0;
	}	
	 
	for(i = 0; i < length; i++)
	{ 
		WnHighPassFilterSingle400hz(In[i],  &Out[i] );	
	}

	/**/
	if(mutecnt400 < mutelen){
		mutecnt400+=length;
		memset(Out,0, length*2);
	}
  }

    