/* ----------------------------------------------------------------------    
* Copyright (C) 2010-2013 ARM Limited. All rights reserved.    
*    
* $Date:        16. October 2013  
* $Revision: 	V1.4.2  
*    
* Project: 	    CMSIS DSP Library    
* Title:	    arm_cfft_f32.c   
*    
* Description:	Combined Radix Decimation in Frequency CFFT Floating point processing function
*    
* Target Processor: Cortex-M4/Cortex-M3/Cortex-M0
*  
* Redistribution and use in source and binary forms, with or without 
* modification, are permitted provided that the following conditions
* are met:
*   - Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   - Redistributions in binary form must reproduce the above copyright
*     notice, this list of conditions and the following disclaimer in
*     the documentation and/or other materials provided with the 
*     distribution.
*   - Neither the name of ARM LIMITED nor the names of its contributors
*     may be used to endorse or promote products derived from this
*     software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
* FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
* COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
* BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.   
* -------------------------------------------------------------------- */


 #include "./include/arm_math.h"
 #include "./include/arm_common_tables.h"

 
 

extern void arm_radix8_butterfly_f32(
  float32_t * pSrc,
  uint16_t fftLen,
  const float32_t * pCoef,
  uint16_t twidCoefModifier);

extern void arm_bitreversal_32(
		uint32_t * pSrc,
		const uint16_t bitRevLen,
		const uint16_t * pBitRevTable);
 
 
 
void arm_cfft_radix8by4_f32( arm_cfft_instance_f32 * S, float32_t * p1) 
{
   uint32_t    L  = S->fftLen >> 1;
   float32_t * pCol1, *pCol2, *pCol3, *pCol4, *pEnd1, *pEnd2, *pEnd3, *pEnd4;
   //const float32_t *tw2, *tw3, *tw4, *tw4_pls;
   const float32_t *tw21, *tw31, *tw41, *tw4_pls1;
   const float32_t *tw21_I, *tw31_I, *tw41_I, *tw4_pls1_I;
   float32_t * p2 = p1 + L;
   float32_t * p3 = p2 + L;
   float32_t * p4 = p3 + L;
   float32_t t2[4], t3[4], t4[4], twR, twI;
   float32_t p1ap3_0, p1sp3_0, p1ap3_1, p1sp3_1;
   float32_t m0, m1, m2, m3;
   uint32_t l, twMod2, twMod3, twMod4;

   float32_t twR__pls, twI__pls;
   uint32_t  twMod21, twMod31, twMod41;

   pCol1 = p1;         // points to real values by default
   pCol2 = p2;
   pCol3 = p3;
   pCol4 = p4;
   pEnd1 = p2 - 1;     // points to imaginary values by default
   pEnd2 = p3 - 1;
   pEnd3 = p4 - 1;
   pEnd4 = pEnd3 + L;
   
  // tw2 = tw3 = tw4 = (float32_t *) S->pTwiddle; // 1024 cosine wav   256 * 4
  // tw4_pls = (float32_t *) &S->pTwiddle[124];  // add for ram cost saving

   tw21 = tw31 = tw41 = (float32_t *) twiddleCoef_256_opt; // 1024 cosine wav   256 * 4
   tw4_pls1 = (float32_t *) &twiddleCoef_256_opt[62];  // add for ram cost saving    
  
   tw21_I = tw31_I = tw41_I = (float32_t *) &twiddleCoef_256_opt[64]; // 1024 cosine wav   256 * 4
   tw4_pls1_I = (float32_t *) &twiddleCoef_256_opt[2];  // add for ram cost saving    
   
   L >>= 1;

   // do four dot Fourier transform

   twMod2 = 2;
   twMod3 = 4;
   twMod4 = 6;
   twMod21 = 1;
   twMod31 = 2;
   twMod41 = 3;

   // TOP
   p1ap3_0 = p1[0] + p3[0];
   p1sp3_0 = p1[0] - p3[0];
   p1ap3_1 = p1[1] + p3[1];
   p1sp3_1 = p1[1] - p3[1];

   // col 2
   t2[0] = p1sp3_0 + p2[1] - p4[1];
   t2[1] = p1sp3_1 - p2[0] + p4[0];
   // col 3
   t3[0] = p1ap3_0 - p2[0] - p4[0];
   t3[1] = p1ap3_1 - p2[1] - p4[1];
   // col 4
   t4[0] = p1sp3_0 - p2[1] + p4[1];
   t4[1] = p1sp3_1 + p2[0] - p4[0];
   // col 1
   *p1++ = p1ap3_0 + p2[0] + p4[0];
   *p1++ = p1ap3_1 + p2[1] + p4[1];

   // Twiddle factors are ones
   *p2++ = t2[0];
   *p2++ = t2[1];
   *p3++ = t3[0];
   *p3++ = t3[1];
   *p4++ = t4[0];
   *p4++ = t4[1];
   
  // tw2 += twMod2;
  // tw3 += twMod3;
  // tw4 += twMod4;

   tw21 += twMod21;
   tw31 += twMod31;
   tw41 += twMod41;
 
   tw21_I -= twMod21;
   tw31_I -= twMod31;
   tw41_I -= twMod41;

   for (l = (L - 2) >> 1; l > 0; l-- ) 
   {

      // TOP
      p1ap3_0 = p1[0] + p3[0];
      p1sp3_0 = p1[0] - p3[0];
      p1ap3_1 = p1[1] + p3[1];
      p1sp3_1 = p1[1] - p3[1];
      // col 2
      t2[0] = p1sp3_0 + p2[1] - p4[1];
      t2[1] = p1sp3_1 - p2[0] + p4[0];
      // col 3
      t3[0] = p1ap3_0 - p2[0] - p4[0];
      t3[1] = p1ap3_1 - p2[1] - p4[1];
      // col 4
      t4[0] = p1sp3_0 - p2[1] + p4[1];
      t4[1] = p1sp3_1 + p2[0] - p4[0];
      // col 1 - top
      *p1++ = p1ap3_0 + p2[0] + p4[0];
      *p1++ = p1ap3_1 + p2[1] + p4[1];

      // BOTTOM
      p1ap3_1 = pEnd1[-1] + pEnd3[-1];
      p1sp3_1 = pEnd1[-1] - pEnd3[-1];
      p1ap3_0 = pEnd1[0] + pEnd3[0];
      p1sp3_0 = pEnd1[0] - pEnd3[0];
      // col 2
      t2[2] = pEnd2[0]  - pEnd4[0] + p1sp3_1;
      t2[3] = pEnd1[0] - pEnd3[0] - pEnd2[-1] + pEnd4[-1];
      // col 3
      t3[2] = p1ap3_1 - pEnd2[-1] - pEnd4[-1];
      t3[3] = p1ap3_0 - pEnd2[0]  - pEnd4[0];
      // col 4
      t4[2] = pEnd2[0]  - pEnd4[0]  - p1sp3_1;
      t4[3] = pEnd4[-1] - pEnd2[-1] - p1sp3_0;
      // col 1 - Bottom
      *pEnd1-- = p1ap3_0 + pEnd2[0]  + pEnd4[0];
      *pEnd1-- = p1ap3_1 + pEnd2[-1] + pEnd4[-1];

      // COL 2
      // read twiddle factors
   //   twR = *tw2++;
   //   twI = *tw2++;
	  
	  twR  = *tw21++;
      twI  = *tw21_I -- ;

      // multiply by twiddle factors
      //  let    Z1 = a + i(b),   Z2 = c + i(d)
      //   =>  Z1 * Z2  =  (a*c - b*d) + i(b*c + a*d)
      // Top
      m0 = t2[0] * twR;
      m1 = t2[1] * twI;
      m2 = t2[1] * twR;
      m3 = t2[0] * twI;
      
      *p2++ = m0 + m1;
      *p2++ = m2 - m3;
      // use vertical symmetry col 2
      // 0.9997 - 0.0245i  <==>  0.0245 - 0.9997i
      // Bottom
      m0 = t2[3] * twI;
      m1 = t2[2] * twR;
      m2 = t2[2] * twI;
      m3 = t2[3] * twR;
      
      *pEnd2-- = m0 - m1;
      *pEnd2-- = m2 + m3;

      // COL 3
  //    twR = tw3[0];
  //    twI = tw3[1];
  //    tw3 += twMod3;

	  twR  = *tw31 ;
      twI  = *tw31_I   ;

	  tw31 += 2;
	  tw31_I -= 2;

      // Top
      m0 = t3[0] * twR;
      m1 = t3[1] * twI;
      m2 = t3[1] * twR;
      m3 = t3[0] * twI;
      
      *p3++ = m0 + m1;
      *p3++ = m2 - m3;
      // use vertical symmetry col 3
      // 0.9988 - 0.0491i  <==>  -0.9988 - 0.0491i
      // Bottom
      m0 = -t3[3] * twR;
      m1 = t3[2] * twI;
      m2 = t3[2] * twR;
      m3 = t3[3] * twI;
      
      *pEnd3-- = m0 - m1;
      *pEnd3-- = m3 - m2;
      
      // COL 4

  
	   if(l>10){
        //    twR  = tw4[0];
        //    twI  = tw4[1];    
		//    tw4 += twMod4;
	        twR  = *tw41  ;
            twI  = *tw41_I;
			tw41 +=3; 
			tw41_I -=3;
	   } 
	   else {
        //    twR  = - tw4_pls[0];
		//	twI  = tw4_pls[1];   
        //    tw4_pls -= twMod4; 

            twR   = - *tw4_pls1;
			twI   = *tw4_pls1_I ;   
            tw4_pls1 -= 3; 
			tw4_pls1_I += 3;
	   }
	    
      // Top
      m0 = t4[0] * twR;
      m1 = t4[1] * twI;
      m2 = t4[1] * twR;
      m3 = t4[0] * twI;
      
      *p4++ = m0 + m1;
      *p4++ = m2 - m3;
      // use vertical symmetry col 4
      // 0.9973 - 0.0736i  <==>  -0.0736 + 0.9973i
      // Bottom
      m0 = t4[3] * twI;
      m1 = t4[2] * twR;
      m2 = t4[2] * twI;
      m3 = t4[3] * twR;
      
      *pEnd4-- = m0 - m1;
      *pEnd4-- = m2 + m3;
   }

   //MIDDLE
   // Twiddle factors are 
   //  1.0000  0.7071-0.7071i  -1.0000i  -0.7071-0.7071i
   p1ap3_0 = p1[0] + p3[0];
   p1sp3_0 = p1[0] - p3[0];
   p1ap3_1 = p1[1] + p3[1];
   p1sp3_1 = p1[1] - p3[1];

   // col 2
   t2[0] = p1sp3_0 + p2[1] - p4[1];
   t2[1] = p1sp3_1 - p2[0] + p4[0];
   // col 3
   t3[0] = p1ap3_0 - p2[0] - p4[0];
   t3[1] = p1ap3_1 - p2[1] - p4[1];
   // col 4
   t4[0] = p1sp3_0 - p2[1] + p4[1];
   t4[1] = p1sp3_1 + p2[0] - p4[0];
   // col 1 - Top
   *p1++ = p1ap3_0 + p2[0] + p4[0];
   *p1++ = p1ap3_1 + p2[1] + p4[1];
   
   // COL 2
  twR = *tw21 ;
  twI = *tw21_I  ;
        
 


   m0 = t2[0] * twR;
   m1 = t2[1] * twI;
   m2 = t2[1] * twR;
   m3 = t2[0] * twI;
   
   *p2++ = m0 + m1;
   *p2++ = m2 - m3;
      // COL 3
   twR  = *tw31 ;
   twI  = *tw31_I   ;
   
   m0 = t3[0] * twR;
   m1 = t3[1] * twI;
   m2 = t3[1] * twR;
   m3 = t3[0] * twI;
   
   *p3++ = m0 + m1;
   *p3++ = m2 - m3;
   // COL 4

 
   twR = - *tw4_pls1;
   twI   = *tw4_pls1_I ;   
 
   
   m0 = t4[0] * twR;
   m1 = t4[1] * twI;
   m2 = t4[1] * twR;
   m3 = t4[0] * twI;
   
   *p4++ = m0 + m1;
   *p4++ = m2 - m3;

   // first col
   arm_radix8_butterfly_f32( pCol1, L, (float32_t *) S->pTwiddle, 4u);
   // second col
   arm_radix8_butterfly_f32( pCol2, L, (float32_t *) S->pTwiddle, 4u);
   // third col
   arm_radix8_butterfly_f32( pCol3, L, (float32_t *) S->pTwiddle, 4u);
   // fourth col
   arm_radix8_butterfly_f32( pCol4, L, (float32_t *) S->pTwiddle, 4u);

}

 
/**
* @addtogroup ComplexFFT   
* @{   
*/

/**   
* @details   
* @brief       Processing function for the floating-point complex FFT.
* @param[in]      *S    points to an instance of the floating-point CFFT structure.  
* @param[in, out] *p1   points to the complex data buffer of size <code>2*fftLen</code>. Processing occurs in-place.  
* @param[in]     ifftFlag       flag that selects forward (ifftFlag=0) or inverse (ifftFlag=1) transform.  
* @param[in]     bitReverseFlag flag that enables (bitReverseFlag=1) or disables (bitReverseFlag=0) bit reversal of output.  
* @return none.  
*/
extern
void arm_bitreversal_f32(
float32_t * pSrc,
uint16_t fftSize,
uint16_t bitRevFactor,
uint16_t * pBitRevTab);
extern
 void biterverse(uint32_t*In1,  uint16_t In2,  uint16_t *In3);

void arm_cfft_f32( 
   const arm_cfft_instance_f32 * S, 
   float32_t * p1,
   uint8_t ifftFlag,
   uint8_t bitReverseFlag)
{

   uint32_t  L = S->fftLen, l;
   float32_t invL, * pSrc;

  if(ifftFlag == 1u)
  {
	  /*  Conjugate input data  */
	  pSrc = p1 + 1;
	  for(l=0; l<L; l++) {
		  *pSrc = -*pSrc;
		   pSrc += 2;
	  }
  }

 
	     arm_cfft_radix8by4_f32  ( (arm_cfft_instance_f32 *) S, p1);
  
	 	biterverse((uint32_t*)p1,S->bitRevLength,(uint16_t*)S->pBitRevTable);

  if(ifftFlag == 1u)
  {
	  invL = 1.0f/(float32_t)L;
	  /*  Conjugate and scale output data */
	  pSrc = p1;
	  for(l=0; l<L; l++) {
  		 *pSrc++ *=   invL ;
		 *pSrc  = -(*pSrc) * invL;
          pSrc++;
	  }
  }
}

