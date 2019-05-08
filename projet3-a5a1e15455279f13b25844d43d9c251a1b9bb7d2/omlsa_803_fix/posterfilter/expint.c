/*
   Copy from matlab expint.m
   yang.yang@nanosic.com
*/
#include"omlsa_def.h"

#include<float.h>
#include<math.h>
#define Inf   FLT_MAX //DBL_MAX
omlsa_float32_t polyval(omlsa_float32_t x, omlsa_float32_t*coeff,int len);
omlsa_float32_t  eps(omlsa_float32_t val);

int j_max = 0;

static omlsa_float32_t  j_f_tab[50] = {
0.0f,	1.000000f, 0.500000f, 0.333333f, 0.250000f, 0.200000f, 0.166667f, 0.142857f, 0.125000f, 0.111111f, 0.100000f, 0.090909f, 0.083333f, 0.076923f, 0.071429f, 0.066667f, 
0.062500f, 0.058824f, 0.055556f, 0.052632f, 0.050000f, 0.047619f, 0.045455f, 0.043478f, 0.041667f, 0.040000f, 0.038462f, 0.037037f, 0.035714f, 0.034483f, 0.033333f, 0.032258f, 
0.031250f, 0.030303f, 0.029412f, 0.028571f, 0.027778f, 0.027027f, 0.026316f, 0.025641f, 0.025000f, 0.024390f, 0.023810f, 0.023256f, 0.022727f, 0.022222f, 0.021739f, 0.021277f, 
0.020833f, 0.020408f, 
};

omlsa_float32_t expint(omlsa_float32_t x){
 

   static omlsa_float32_t p[9] =  {
	                   -3.602693626336023e-09f, -4.819538452140960e-07f, -2.569498322115933e-05f, 
                       -6.973790859534190e-04f, -1.019573529845792e-02f, -7.811863559248197e-02f,
                       -3.012432892762715e-01f, -7.773807325735529e-01f,  8.267661952366478e+00f};

    static int n; omlsa_float32_t am2,am1,bm2,bm1; omlsa_float32_t alpha,a , b,beta_l;
     omlsa_float32_t polyv,y;
     omlsa_float32_t xk, yk, pterm, term,f,oldf ; 
     int j;  omlsa_float32_t j_f , b_f;
 
     polyv = polyval(x, p, 9); 
   if(0<=polyv){
       omlsa_float32_t egamma = -0.57721566490153286061f;
       
       xk = x;
       yk =  egamma - logf(xk);
	   j = 1;
       pterm = xk;
       term = xk;
 
     while(fabsf(term)> eps(yk)){
         yk = yk+term;


         j=j+1;

#if 1
		 if(j >=50 ){
		    j_f = 1.0f/(omlsa_float32_t)j;
		 }
		 else{
		 
			 j_f = j_f_tab[j ];
		 }
		  
         pterm = -xk*pterm*j_f;
         term = pterm*j_f;
#else 
         pterm = -xk*pterm/j;
         term = pterm/j;
#endif
     }

       y = yk;
   }
  
  
   if(polyv<0.0f){
       n = 1; 
	   xk = x;
	   am2 = 0.0f; bm2 = 0.0f;  am1 = 1.0f;  
       bm1 = xk;  f = am1 / bm1;
       oldf = (omlsa_float32_t)Inf;
     

       if(fabsf(f-oldf)> (2.2204e-14f *fabsf(f))){  //  eps(class(f)) = eps('omlsa_float32_t') 浮点相对精度
         
#if 1 
		    // alpha = n-1+(j/2);           
             a = am1 +   am2;
             b = bm1 +   bm2;
			  
			 b_f = 1.0f/b;
           
             am2 = am1*b_f;
             bm2 = bm1*b_f;
             am1 = a*b_f;
             bm1 = 1;
             f = am1;
            
                          // calculate 
            // calculate the coefficients for j odd
             
             beta_l = xk;
             a = beta_l * am1 + am2;
             b = beta_l + 1.0f;

			 b_f = 1.0f/b;
             am2 = am1 *b_f;
             bm2 = b_f;
             am1 = a*b_f;
             
             oldf = f;
             f = am1;
             
#else 
             alpha = n-1+(j/2); 
          
             a = am1 + alpha * am2;
             b = bm1 + alpha * bm2;
           
             am2 = am1 / b;
             bm2 = bm1 / b;
             am1 = a / b;
             bm1 = 1;
             f = am1;
             j = j+1;
                          // calculate 
            // calculate the coefficients for j odd
             alpha = (j-1)/2;
             beta_l = xk;
             a = beta_l * am1 + alpha * am2;
             b = beta_l * bm1 + alpha * bm2;
             am2 = am1 / b;
             bm2 = bm1 / b;
             am1 = a / b;
             bm1 = 1;
             oldf = f;
             f = am1;
             j = j+1;

#endif
       }

       y = expf(-xk) * f;// - 1i*pi*(((xk)<0)&(imag(xk)==0)); 
	    
   }


//
//y = reshape(y,siz);
    return (omlsa_float32_t)y;
}

/*
   polyval from matlab function 'polyval'.

   y = x^8 + x^7 + x^6+x^5+x^4+x^3+x^2+x^1 + 1;
*/

omlsa_float32_t polyval(omlsa_float32_t x, omlsa_float32_t*coeff,int len){
    
    omlsa_float32_t temp1, temp2, rev;
	omlsa_float32_t x_2 =  x*x;
    temp1 = x;	
    temp2 = x*x;

	rev = coeff[8];

    rev += temp1 * coeff[7] + temp2 *coeff[6];  temp1 *=x_2; temp2 *=x_2;
    rev += temp1 * coeff[5] + temp2 *coeff[4];  temp1 *=x_2; temp2 *=x_2;
    rev += temp1 * coeff[3] + temp2 *coeff[2];  temp1 *=x_2; temp2 *=x_2;
    rev += temp1 * coeff[1] + temp2 *coeff[0];
    return rev;
}


#define eps1  2.2204e-16f
omlsa_float32_t  eps(omlsa_float32_t val){
    omlsa_float32_t rev;
    omlsa_float32_t  temp;
    temp = 1.0f;
    val = fabs(val);

    rev = eps1;

    if(val>1.0f){         
         while(temp> val){

#if 1
             temp*=2.0f;
             rev*=2.0f;
#else 
             temp*=2;
             rev*=2;
#endif
         }
    }
    else{
          while( temp>val){
#if 1
             temp*=0.5f;
             rev*=0.5f;
#else 
             temp/=2;
             rev/=2;
#endif
         }        
    }

 
	return rev;

}