/* ----------------------------------------------------------------------    
* Copyright (C) 2010-2013 ARM Limited. All rights reserved.    
*    
* $Date:        16. October 2013  
* $Revision: 	V1.4.2  
*    
* Project: 	    CMSIS DSP Library    
* Title:	    arm_cfft_radix8_f32.c    
*    
* Description:	Radix-8 Decimation in Frequency CFFT & CIFFT Floating point processing function        
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
 
/**    
* @ingroup groupTransforms    
*/

/**    
* @defgroup Radix8_CFFT_CIFFT Radix-8 Complex FFT Functions    
*    
* \par    
* Complex Fast Fourier Transform(CFFT) and Complex Inverse Fast Fourier Transform(CIFFT) is an efficient algorithm to compute Discrete Fourier Transform(DFT) and Inverse Discrete Fourier Transform(IDFT).    
* Computational complexity of CFFT reduces drastically when compared to DFT.    
* \par    
* This set of functions implements CFFT/CIFFT    
* for floating-point data types.  The functions operates on in-place buffer which uses same buffer for input and output.    
* Complex input is stored in input buffer in an interleaved fashion.    
*    
* \par    
* The functions operate on blocks of input and output data and each call to the function processes    
* <code>2*fftLen</code> samples through the transform.  <code>pSrc</code>  points to In-place arrays containing <code>2*fftLen</code> values.    
* \par   
* The <code>pSrc</code> points to the array of in-place buffer of size <code>2*fftLen</code> and inputs and outputs are stored in an interleaved fashion as shown below.    
* <pre> {real[0], imag[0], real[1], imag[1],..} </pre>    
*    
* \par Lengths supported by the transform:   
* \par    
* Internally, the function utilize a Radix-8 decimation in frequency(DIF) algorithm    
* and the size of the FFT supported are of the lengths [ 64, 512, 4096].   
*     
*    
* \par Algorithm:    
*    
* <b>Complex Fast Fourier Transform:</b>    
* \par     
* Input real and imaginary data:    
* <pre>    
* x(n) = xa + j * ya    
* x(n+N/4 ) = xb + j * yb    
* x(n+N/2 ) = xc + j * yc    
* x(n+3N 4) = xd + j * yd    
* </pre>    
* where N is length of FFT    
* \par    
* Output real and imaginary data:    
* <pre>    
* X(4r) = xa'+ j * ya'    
* X(4r+1) = xb'+ j * yb'    
* X(4r+2) = xc'+ j * yc'    
* X(4r+3) = xd'+ j * yd'    
* </pre>    
* \par    
* Twiddle factors for Radix-8 FFT:    
* <pre>    
* Wn = co1 + j * (- si1)    
* W2n = co2 + j * (- si2)    
* W3n = co3 + j * (- si3)    
* </pre>    
*    
* \par    
* \image html CFFT.gif "Radix-8 Decimation-in Frequency Complex Fast Fourier Transform"    
*    
* \par    
* Output from Radix-8 CFFT Results in Digit reversal order. Interchange middle two branches of every butterfly results in Bit reversed output.    
* \par    
* <b> Butterfly CFFT equations:</b>    
* <pre>    
* xa' = xa + xb + xc + xd    
* ya' = ya + yb + yc + yd    
* xc' = (xa+yb-xc-yd)* co1 + (ya-xb-yc+xd)* (si1)    
* yc' = (ya-xb-yc+xd)* co1 - (xa+yb-xc-yd)* (si1)    
* xb' = (xa-xb+xc-xd)* co2 + (ya-yb+yc-yd)* (si2)    
* yb' = (ya-yb+yc-yd)* co2 - (xa-xb+xc-xd)* (si2)    
* xd' = (xa-yb-xc+yd)* co3 + (ya+xb-yc-xd)* (si3)    
* yd' = (ya+xb-yc-xd)* co3 - (xa-yb-xc+yd)* (si3)    
* </pre>    
*    
* \par    
* where <code>fftLen</code> length of CFFT/CIFFT; <code>ifftFlag</code> Flag for selection of CFFT or CIFFT(Set ifftFlag to calculate CIFFT otherwise calculates CFFT);    
* <code>bitReverseFlag</code> Flag for selection of output order(Set bitReverseFlag to output in normal order otherwise output in bit reversed order);     
* <code>pTwiddle</code>points to array of twiddle coefficients; <code>pBitRevTable</code> points to the array of bit reversal table.    
* <code>twidCoefModifier</code> modifier for twiddle factor table which supports all FFT lengths with same table;     
* <code>pBitRevTable</code> modifier for bit reversal table which supports all FFT lengths with same table.    
* <code>onebyfftLen</code> value of 1/fftLen to calculate CIFFT;    
*   
* \par Fixed-Point Behavior    
* Care must be taken when using the fixed-point versions of the CFFT/CIFFT function.    
* Refer to the function specific documentation below for usage guidelines.    
*/


/*    
* @brief  Core function for the floating-point CFFT butterfly process.   
* @param[in, out] *pSrc            points to the in-place buffer of floating-point data type.   
* @param[in]      fftLen           length of the FFT.   
* @param[in]      *pCoef           points to the twiddle coefficient buffer.   
* @param[in]      twidCoefModifier twiddle coefficient modifier that supports different size FFTs with the same twiddle factor table.   
* @return none.   
*/
typedef float float32_t;
typedef int int32_t;
typedef unsigned int uint32_t;
typedef unsigned short uint16_t;
typedef short int16_t;
typedef unsigned char uint8_t;
typedef char int8_t;
 

/**    
* @} end of Radix8_CFFT_CIFFT group    
*/
#if 0

extern float32_t twiddleCoef_256_opt[65];

void arm_radix8_butterfly_f32 (
float32_t * pSrc,
uint16_t fftLen,
const float32_t * pCoef,
uint16_t twidCoefModifier)
{
   uint32_t ia1, ia2, ia3, ia4, ia5, ia6, ia7;
   uint32_t i1, i2, i3, i4, i5, i6, i7, i8;
   uint32_t id;
   uint32_t n1, n2, j;
   
   float32_t r1, r2, r3, r4, r5, r6, r7, r8;
   float32_t t1, t2;
   float32_t s1, s2, s3, s4, s5, s6, s7, s8;
   float32_t p1, p2, p3, p4;
   float32_t co2, co3, co4, co5, co6, co7, co8;
   float32_t si2, si3, si4, si5, si6, si7, si8;

   float32_t co2_pls, co3_pls, co4_pls, co5_pls, co6_pls, co7_pls, co8_pls;
   float32_t si2_pls, si3_pls, si4_pls, si5_pls, si6_pls, si7_pls, si8_pls;

   float32_t * pCoef8, *pCoef16, *pCoef24, *pCoef32,*pCoef40,*pCoef48 , *pCoef56;
   float32_t * pCoefR, *pCoefI ;



   const float32_t C81 = 0.70710678118f;
   
   pCoefR =  &twiddleCoef_256_opt[0];
   pCoefI =  &twiddleCoef_256_opt[0];

   n2 = fftLen;
   
   do 
   {
      n1 = n2;
      n2 = n2 >> 3;
      i1 = 0;
      
      do
      {
         i2 = i1 + n2;
         i3 = i2 + n2;
         i4 = i3 + n2;
         i5 = i4 + n2;
         i6 = i5 + n2;
         i7 = i6 + n2;
         i8 = i7 + n2;
         r1 = pSrc[2 * i1] + pSrc[2 * i5];
         r5 = pSrc[2 * i1] - pSrc[2 * i5];
         r2 = pSrc[2 * i2] + pSrc[2 * i6];
         r6 = pSrc[2 * i2] - pSrc[2 * i6];
         r3 = pSrc[2 * i3] + pSrc[2 * i7];
         r7 = pSrc[2 * i3] - pSrc[2 * i7];
         r4 = pSrc[2 * i4] + pSrc[2 * i8];
         r8 = pSrc[2 * i4] - pSrc[2 * i8];
         t1 = r1 - r3;
         r1 = r1 + r3;
         r3 = r2 - r4;
         r2 = r2 + r4;
         pSrc[2 * i1] = r1 + r2;   
         pSrc[2 * i5] = r1 - r2;
         r1 = pSrc[2 * i1 + 1] + pSrc[2 * i5 + 1];
         s5 = pSrc[2 * i1 + 1] - pSrc[2 * i5 + 1];
         r2 = pSrc[2 * i2 + 1] + pSrc[2 * i6 + 1];
         s6 = pSrc[2 * i2 + 1] - pSrc[2 * i6 + 1];
         s3 = pSrc[2 * i3 + 1] + pSrc[2 * i7 + 1];
         s7 = pSrc[2 * i3 + 1] - pSrc[2 * i7 + 1];
         r4 = pSrc[2 * i4 + 1] + pSrc[2 * i8 + 1];
         s8 = pSrc[2 * i4 + 1] - pSrc[2 * i8 + 1];
         t2 = r1 - s3;
         r1 = r1 + s3;
         s3 = r2 - r4;
         r2 = r2 + r4;
         pSrc[2 * i1 + 1] = r1 + r2;
         pSrc[2 * i5 + 1] = r1 - r2;
         pSrc[2 * i3]     = t1 + s3;
         pSrc[2 * i7]     = t1 - s3;
         pSrc[2 * i3 + 1] = t2 - r3;
         pSrc[2 * i7 + 1] = t2 + r3;
         r1 = (r6 - r8) * C81;
         r6 = (r6 + r8) * C81;
         r2 = (s6 - s8) * C81;
         s6 = (s6 + s8) * C81;
         t1 = r5 - r1;
         r5 = r5 + r1;
         r8 = r7 - r6;
         r7 = r7 + r6;
         t2 = s5 - r2;
         s5 = s5 + r2;
         s8 = s7 - s6;
         s7 = s7 + s6;
         pSrc[2 * i2]     = r5 + s7;
         pSrc[2 * i8]     = r5 - s7;
         pSrc[2 * i6]     = t1 + s8;
         pSrc[2 * i4]     = t1 - s8;
         pSrc[2 * i2 + 1] = s5 - r7;
         pSrc[2 * i8 + 1] = s5 + r7;
         pSrc[2 * i6 + 1] = t2 - r8;
         pSrc[2 * i4 + 1] = t2 + r8;
         
         i1 += n1;
      } while(i1 < fftLen);
      
      if(n2 < 8)
         break;
      
      ia1 = 0;
      j = 1;


	  
      
      do
      {      
         /*  index calculation for the coefficients */
         id  = ia1 + twidCoefModifier;
         ia1 = id;
         ia2 = ia1 + id;
         ia3 = ia2 + id;
         ia4 = ia3 + id;
         ia5 = ia4 + id;
         ia6 = ia5 + id;
         ia7 = ia6 + id;
 
		 switch(j){
		              
		 case 1:      co2 = pCoefR[4]; co3 = pCoefR[8];co4 = pCoefR[12];co5 = pCoefR[16];
			          co6 = pCoefR[20]; co7 = pCoefR[24];co8 = pCoefR[28]; 
 
					  si2=  pCoefR[60]; si3=  pCoefR[56];si4=  pCoefR[52];si5=  pCoefR[48];
			          si6=  pCoefR[44]; si7=  pCoefR[48];si8=  pCoefR[36]; 	
  
			 break;
		 case 2:      co2 = pCoefR[8]; co3 = pCoefR[16];co4 = pCoefR[24];co5 = pCoefR[32];
			          co6 = pCoefR[40]; co7 = pCoefR[48];co8 = pCoefR[56]; 
		  
					  si2=  pCoefR[56]; si3=  pCoefR[48];si4=  pCoefR[40];si5=  pCoefR[32];
			          si6=  pCoefR[24]; si7=  pCoefR[16];si8=   pCoefR[8]; 	
  
			 break;
		 case 3: 

					  co2 = pCoefR[12]; co3 = pCoefR[24];co4 = pCoefR[36];co5 = pCoefR[48];
			          co6 = pCoefR[60]; co7 = -pCoefR[56];co8 = - pCoefR[44]; 
					    
					  si2=  pCoefR[64-12]; si3=  pCoefR[64-24];si4=  pCoefR[64-36];si5=  pCoefR[64-48];
			          si6=  pCoefR[64-60]; si7=  pCoefR[64-56];si8=    pCoefR[64-44]; 	
 
			 break;
			          
		 case 4:      co2 = pCoefR[16]; co3 = pCoefR[32];co4 = pCoefR[48];co5 = pCoefR[64];
			          co6 =- pCoefR[48]; co7 = -pCoefR[32];co8 = - pCoefR[16]; 
 
					  si2=  pCoefI[64-16]; si3=  pCoefI[64-32];si4=  pCoefI[64-48];si5=  pCoefI[0];
			          si6=   pCoefI[64-48]; si7=   pCoefI[64-32];si8=    pCoefI[64-16]; 	
 

			 break;
			          
		 case 5:      co2 = pCoefR[20]; co3 = pCoefR[40];co4 = pCoefR[60];co5 =- pCoefR[48];
			          co6 =- pCoefR[28]; co7 = -pCoefR[8];co8 = - pCoefR[12]; 
					   
					  si2=  pCoefI[44]; si3=  pCoefI[24];si4=  pCoefI[4];si5=  pCoefI[16];
			          si6=  pCoefI[36]; si7=  pCoefI[56];si8=  - pCoefI[52]; 	
 
			 break;
			         
		 case 6:      co2 = pCoefR[24]; co3 = pCoefR[48];co4 =- pCoefR[56];co5 =- pCoefR[32];
			          co6 =- pCoefR[8]; co7 = -pCoefR[16];co8 = - pCoefR[40]; 

					  si2=  pCoefI[40]; si3=  pCoefI[16];si4=  pCoefI[8];si5=  pCoefI[32];
			          si6=  pCoefI[56]; si7= - pCoefI[48];si8=  - pCoefI[24]; 	
  

			 break;
			         
		 case 7:     
					  co2 = pCoefR[28]; co3 = pCoefR[56];co4 =- pCoefR[44];co5 =- pCoefR[16];
			          co6 =- pCoefR[12]; co7 = -pCoefR[40];co8 =  pCoefR[60]; 


					  si2=  pCoefI[36]; si3=  pCoefI[8];si4=  pCoefI[20];si5=   pCoefI[48];
			          si6=  -pCoefI[52]; si7=  -pCoefI[24];si8= - pCoefI[4]; 	
 
 

			 break;		 
		 default: break;
		 }
		     
	  
         
         
         i1 = j;
         
         do
         {
            /*  index calculation for the input */
            i2 = i1 + n2;
            i3 = i2 + n2;
            i4 = i3 + n2;
            i5 = i4 + n2;
            i6 = i5 + n2;
            i7 = i6 + n2;
            i8 = i7 + n2;
            r1 = pSrc[2 * i1] + pSrc[2 * i5];
            r5 = pSrc[2 * i1] - pSrc[2 * i5];
            r2 = pSrc[2 * i2] + pSrc[2 * i6];
            r6 = pSrc[2 * i2] - pSrc[2 * i6];
            r3 = pSrc[2 * i3] + pSrc[2 * i7];
            r7 = pSrc[2 * i3] - pSrc[2 * i7];
            r4 = pSrc[2 * i4] + pSrc[2 * i8];
            r8 = pSrc[2 * i4] - pSrc[2 * i8];
            t1 = r1 - r3;
            r1 = r1 + r3;
            r3 = r2 - r4;
            r2 = r2 + r4;
            pSrc[2 * i1] = r1 + r2;
            r2 = r1 - r2;
            s1 = pSrc[2 * i1 + 1] + pSrc[2 * i5 + 1];
            s5 = pSrc[2 * i1 + 1] - pSrc[2 * i5 + 1];
            s2 = pSrc[2 * i2 + 1] + pSrc[2 * i6 + 1];
            s6 = pSrc[2 * i2 + 1] - pSrc[2 * i6 + 1];
            s3 = pSrc[2 * i3 + 1] + pSrc[2 * i7 + 1];
            s7 = pSrc[2 * i3 + 1] - pSrc[2 * i7 + 1];
            s4 = pSrc[2 * i4 + 1] + pSrc[2 * i8 + 1];
            s8 = pSrc[2 * i4 + 1] - pSrc[2 * i8 + 1];
            t2 = s1 - s3;
            s1 = s1 + s3;
            s3 = s2 - s4;
            s2 = s2 + s4;
            r1 = t1 + s3;
            t1 = t1 - s3;
            pSrc[2 * i1 + 1] = s1 + s2;
            s2 = s1 - s2;
            s1 = t2 - r3;
            t2 = t2 + r3;
            p1 = co5 * r2;
            p2 = si5 * s2;
            p3 = co5 * s2;
            p4 = si5 * r2;
            pSrc[2 * i5]     = p1 + p2;
            pSrc[2 * i5 + 1] = p3 - p4;
            p1 = co3 * r1;
            p2 = si3 * s1;
            p3 = co3 * s1;
            p4 = si3 * r1;
            pSrc[2 * i3]     = p1 + p2;
            pSrc[2 * i3 + 1] = p3 - p4;
            p1 = co7 * t1;
            p2 = si7 * t2;
            p3 = co7 * t2;
            p4 = si7 * t1;
            pSrc[2 * i7]     = p1 + p2;
            pSrc[2 * i7 + 1] = p3 - p4;
            r1 = (r6 - r8) * C81;
            r6 = (r6 + r8) * C81;
            s1 = (s6 - s8) * C81;
            s6 = (s6 + s8) * C81;
            t1 = r5 - r1;
            r5 = r5 + r1;
            r8 = r7 - r6;
            r7 = r7 + r6;
            t2 = s5 - s1;
            s5 = s5 + s1;
            s8 = s7 - s6;
            s7 = s7 + s6;
            r1 = r5 + s7;
            r5 = r5 - s7;
            r6 = t1 + s8;
            t1 = t1 - s8;
            s1 = s5 - r7;
            s5 = s5 + r7;
            s6 = t2 - r8;
            t2 = t2 + r8;
            p1 = co2 * r1;
            p2 = si2 * s1;
            p3 = co2 * s1;
            p4 = si2 * r1;
            pSrc[2 * i2]     = p1 + p2;
            pSrc[2 * i2 + 1] = p3 - p4;
            p1 = co8 * r5;
            p2 = si8 * s5;
            p3 = co8 * s5;
            p4 = si8 * r5;
            pSrc[2 * i8]     = p1 + p2;
            pSrc[2 * i8 + 1] = p3 - p4;
            p1 = co6 * r6;
            p2 = si6 * s6;
            p3 = co6 * s6;
            p4 = si6 * r6;
            pSrc[2 * i6]     = p1 + p2;
            pSrc[2 * i6 + 1] = p3 - p4;
            p1 = co4 * t1;
            p2 = si4 * t2;
            p3 = co4 * t2;
            p4 = si4 * t1;
            pSrc[2 * i4]     = p1 + p2;
            pSrc[2 * i4 + 1] = p3 - p4;
            
            i1 += n1;
         } while(i1 < fftLen);
         
         j++;
      } while(j < n2);
      
      twidCoefModifier <<= 3;
   } while(n2 > 7);   
}
#else 

#define OPT 1

extern float32_t twiddleCoef_256_opt[65];

void arm_radix8_butterfly_f32 (
float32_t * pSrc,
uint16_t fftLen,
const float32_t * pCoef,
uint16_t twidCoefModifier)
{
   uint32_t ia1, ia2, ia3, ia4, ia5, ia6, ia7;
   uint32_t i1, i2, i3, i4, i5, i6, i7, i8;
   uint32_t id;
   uint32_t n1, n2, j;
   
   float32_t r1, r2, r3, r4, r5, r6, r7, r8;
   float32_t t1, t2;
   float32_t s1, s2, s3, s4, s5, s6, s7, s8;
   float32_t p1, p2, p3, p4;
   float32_t co2, co3, co4, co5, co6, co7, co8;
   float32_t si2, si3, si4, si5, si6, si7, si8;

   float32_t co2_pls, co3_pls, co4_pls, co5_pls, co6_pls, co7_pls, co8_pls;
   float32_t si2_pls, si3_pls, si4_pls, si5_pls, si6_pls, si7_pls, si8_pls;

   float32_t * pCoef8, *pCoef16, *pCoef24, *pCoef32,*pCoef40,*pCoef48 , *pCoef56;
   float32_t * pCoefR, *pCoefI ;



   const float32_t C81 = 0.70710678118f;
   
   pCoefR =  &twiddleCoef_256_opt[0];
   pCoefI =  &twiddleCoef_256_opt[0];

   n2 = fftLen;
   
   do 
   {
      n1 = n2;
      n2 = n2 >> 3;
      i1 = 0;
      
      do
      {
         i2 = i1 + n2;
         i3 = i2 + n2;
         i4 = i3 + n2;
         i5 = i4 + n2;
         i6 = i5 + n2;
         i7 = i6 + n2;
         i8 = i7 + n2;
         r1 = pSrc[2 * i1] + pSrc[2 * i5];
         r5 = pSrc[2 * i1] - pSrc[2 * i5];
         r2 = pSrc[2 * i2] + pSrc[2 * i6];
         r6 = pSrc[2 * i2] - pSrc[2 * i6];
         r3 = pSrc[2 * i3] + pSrc[2 * i7];
         r7 = pSrc[2 * i3] - pSrc[2 * i7];
         r4 = pSrc[2 * i4] + pSrc[2 * i8];
         r8 = pSrc[2 * i4] - pSrc[2 * i8];
         t1 = r1 - r3;
         r1 = r1 + r3;
         r3 = r2 - r4;
         r2 = r2 + r4;
         pSrc[2 * i1] = r1 + r2;   
         pSrc[2 * i5] = r1 - r2;
         r1 = pSrc[2 * i1 + 1] + pSrc[2 * i5 + 1];
         s5 = pSrc[2 * i1 + 1] - pSrc[2 * i5 + 1];
         r2 = pSrc[2 * i2 + 1] + pSrc[2 * i6 + 1];
         s6 = pSrc[2 * i2 + 1] - pSrc[2 * i6 + 1];
         s3 = pSrc[2 * i3 + 1] + pSrc[2 * i7 + 1];
         s7 = pSrc[2 * i3 + 1] - pSrc[2 * i7 + 1];
         r4 = pSrc[2 * i4 + 1] + pSrc[2 * i8 + 1];
         s8 = pSrc[2 * i4 + 1] - pSrc[2 * i8 + 1];
         t2 = r1 - s3;
         r1 = r1 + s3;
         s3 = r2 - r4;
         r2 = r2 + r4;
         pSrc[2 * i1 + 1] = r1 + r2;
         pSrc[2 * i5 + 1] = r1 - r2;
         pSrc[2 * i3]     = t1 + s3;
         pSrc[2 * i7]     = t1 - s3;
         pSrc[2 * i3 + 1] = t2 - r3;
         pSrc[2 * i7 + 1] = t2 + r3;
         r1 = (r6 - r8) * C81;
         r6 = (r6 + r8) * C81;
         r2 = (s6 - s8) * C81;
         s6 = (s6 + s8) * C81;
         t1 = r5 - r1;
         r5 = r5 + r1;
         r8 = r7 - r6;
         r7 = r7 + r6;
         t2 = s5 - r2;
         s5 = s5 + r2;
         s8 = s7 - s6;
         s7 = s7 + s6;
         pSrc[2 * i2]     = r5 + s7;
         pSrc[2 * i8]     = r5 - s7;
         pSrc[2 * i6]     = t1 + s8;
         pSrc[2 * i4]     = t1 - s8;
         pSrc[2 * i2 + 1] = s5 - r7;
         pSrc[2 * i8 + 1] = s5 + r7;
         pSrc[2 * i6 + 1] = t2 - r8;
         pSrc[2 * i4 + 1] = t2 + r8;
         
         i1 += n1;
      } while(i1 < fftLen);
      
      if(n2 < 8)
         break;
      
      ia1 = 0;
      j = 1;


	  
      
      do
      {      
         /*  index calculation for the coefficients */
         id  = ia1 + twidCoefModifier;
         ia1 = id;
         ia2 = ia1 + id;
         ia3 = ia2 + id;
         ia4 = ia3 + id;
         ia5 = ia4 + id;
         ia6 = ia5 + id;
         ia7 = ia6 + id;
 
		 switch(j){
		              
		 case 1:    
					  co2 = pCoefR[4]; co3 = pCoefR[8];co4 = pCoefR[12];co5 = pCoefR[16];
			          co6 = pCoefR[20]; co7 = pCoefR[24];co8 = pCoefR[28]; 
 
 					  si2=  pCoefR[60]; si3=  pCoefR[64-8];si4=  pCoefR[64-12];si5=  pCoefR[64-16];
			          si6=  pCoefR[44]; si7=  pCoefR[64-24];si8=  pCoefR[64-28]; 	

 

			 break;
		 case 2:     
			        
					  co2 = pCoefR[8]; co3 = pCoefR[16];co4 = pCoefR[24];co5 = pCoefR[32];
			          co6 = pCoefR[40]; co7 = pCoefR[48];co8 = pCoefR[56]; 
  
 
				 	  si2=  pCoefR[64-8]; si3=  pCoefR[64-16];si4=  pCoefR[64-24];si5=  pCoefR[64-32];
			          si6=  pCoefR[64-40]; si7=  pCoefR[64-48];si8=   pCoefR[64-56]; 	
 		 
			 break;
		 case 3:    
					  co2 = pCoefR[12]; co3 = pCoefR[24];co4 = pCoefR[36];co5 = pCoefR[48];
			          co6 = pCoefR[60]; co7 = -pCoefR[56];co8 = - pCoefR[44]; 
 
					  si2=  pCoefR[64-12]; si3=  pCoefR[64-24];si4=  pCoefR[64-36];si5=  pCoefR[64-48];
			          si6=  pCoefR[64-60]; si7=  pCoefR[64-56];si8=    pCoefR[64-44]; 	
 
			 break;
			          
		 case 4:   
					  co2 = pCoefR[16]; co3 = pCoefR[32];co4 = pCoefR[48];co5 = pCoefR[64];
			          co6 =- pCoefR[48]; co7 = -pCoefR[32];co8 = - pCoefR[16]; 

 
					  si2=  pCoefR[64-16]; si3=  pCoefR[64-32];si4=  pCoefR[64-48];si5=  pCoefR[0];
			          si6=   pCoefR[64-48]; si7=   pCoefR[64-32];si8=    pCoefR[64-16]; 	
 

			 break;
			          
		 case 5:   
					  co2 = pCoefR[20]; co3 = pCoefR[40];co4 = pCoefR[60];co5 =- pCoefR[48];
			          co6 =- pCoefR[28]; co7 = -pCoefR[8];co8 = - pCoefR[12]; 

 
			 		  si2=  pCoefR[64-20]; si3=  pCoefR[64-40];si4=  pCoefR[64-60];si5=  pCoefR[64-48];
			          si6=  pCoefR[64-28]; si7=  pCoefR[64-8];si8=  - pCoefR[64-12]; 	
 			 

			 break;
			         
		 case 6:       co2 = pCoefR[24]; co3 = pCoefR[48];co4 =- pCoefR[56];co5 =- pCoefR[32];
			          co6 =- pCoefR[8]; co7 = -pCoefR[16];co8 = - pCoefR[40]; 
  
 
					  si2=  pCoefR[64-24]; si3=  pCoefR[64-48];si4=  pCoefR[64-56];si5=  pCoefR[64-32];
			          si6=  pCoefR[64-8]; si7= - pCoefR[64-16];si8=  - pCoefR[64-40]; 	
 

			 break;
			         
		 case 7:     
					  co2 = pCoefR[28]; co3 = pCoefR[56];co4 =- pCoefR[44];co5 =- pCoefR[16];
			          co6 =- pCoefR[12]; co7 = -pCoefR[40];co8 =  pCoefR[60]; 
 

					  si2=  pCoefR[64-28]; si3=  pCoefR[64-56];si4=  pCoefR[64-44];si5=   pCoefR[64-16];
			          si6=  -pCoefR[64-12]; si7=  -pCoefR[64-40];si8= - pCoefR[64-60]; 	
 

					  // 57	113	89	33	25	81	137 

			 break;		 
		 default: break;
		 }
		     
	  
 
         
         
         i1 = j;
         
         do
         {
            /*  index calculation for the input */
            i2 = i1 + n2;
            i3 = i2 + n2;
            i4 = i3 + n2;
            i5 = i4 + n2;
            i6 = i5 + n2;
            i7 = i6 + n2;
            i8 = i7 + n2;
            r1 = pSrc[2 * i1] + pSrc[2 * i5];
            r5 = pSrc[2 * i1] - pSrc[2 * i5];
            r2 = pSrc[2 * i2] + pSrc[2 * i6];
            r6 = pSrc[2 * i2] - pSrc[2 * i6];
            r3 = pSrc[2 * i3] + pSrc[2 * i7];
            r7 = pSrc[2 * i3] - pSrc[2 * i7];
            r4 = pSrc[2 * i4] + pSrc[2 * i8];
            r8 = pSrc[2 * i4] - pSrc[2 * i8];
            t1 = r1 - r3;
            r1 = r1 + r3;
            r3 = r2 - r4;
            r2 = r2 + r4;
            pSrc[2 * i1] = r1 + r2;
            r2 = r1 - r2;
            s1 = pSrc[2 * i1 + 1] + pSrc[2 * i5 + 1];
            s5 = pSrc[2 * i1 + 1] - pSrc[2 * i5 + 1];
            s2 = pSrc[2 * i2 + 1] + pSrc[2 * i6 + 1];
            s6 = pSrc[2 * i2 + 1] - pSrc[2 * i6 + 1];
            s3 = pSrc[2 * i3 + 1] + pSrc[2 * i7 + 1];
            s7 = pSrc[2 * i3 + 1] - pSrc[2 * i7 + 1];
            s4 = pSrc[2 * i4 + 1] + pSrc[2 * i8 + 1];
            s8 = pSrc[2 * i4 + 1] - pSrc[2 * i8 + 1];
            t2 = s1 - s3;
            s1 = s1 + s3;
            s3 = s2 - s4;
            s2 = s2 + s4;
            r1 = t1 + s3;
            t1 = t1 - s3;
            pSrc[2 * i1 + 1] = s1 + s2;
            s2 = s1 - s2;
            s1 = t2 - r3;
            t2 = t2 + r3;
            p1 = co5 * r2;
            p2 = si5 * s2;
            p3 = co5 * s2;
            p4 = si5 * r2;
            pSrc[2 * i5]     = p1 + p2;
            pSrc[2 * i5 + 1] = p3 - p4;
            p1 = co3 * r1;
            p2 = si3 * s1;
            p3 = co3 * s1;
            p4 = si3 * r1;
            pSrc[2 * i3]     = p1 + p2;
            pSrc[2 * i3 + 1] = p3 - p4;
            p1 = co7 * t1;
            p2 = si7 * t2;
            p3 = co7 * t2;
            p4 = si7 * t1;
            pSrc[2 * i7]     = p1 + p2;
            pSrc[2 * i7 + 1] = p3 - p4;
            r1 = (r6 - r8) * C81;
            r6 = (r6 + r8) * C81;
            s1 = (s6 - s8) * C81;
            s6 = (s6 + s8) * C81;
            t1 = r5 - r1;
            r5 = r5 + r1;
            r8 = r7 - r6;
            r7 = r7 + r6;
            t2 = s5 - s1;
            s5 = s5 + s1;
            s8 = s7 - s6;
            s7 = s7 + s6;
            r1 = r5 + s7;
            r5 = r5 - s7;
            r6 = t1 + s8;
            t1 = t1 - s8;
            s1 = s5 - r7;
            s5 = s5 + r7;
            s6 = t2 - r8;
            t2 = t2 + r8;
            p1 = co2 * r1;
            p2 = si2 * s1;
            p3 = co2 * s1;
            p4 = si2 * r1;
            pSrc[2 * i2]     = p1 + p2;
            pSrc[2 * i2 + 1] = p3 - p4;
            p1 = co8 * r5;
            p2 = si8 * s5;
            p3 = co8 * s5;
            p4 = si8 * r5;
            pSrc[2 * i8]     = p1 + p2;
            pSrc[2 * i8 + 1] = p3 - p4;
            p1 = co6 * r6;
            p2 = si6 * s6;
            p3 = co6 * s6;
            p4 = si6 * r6;
            pSrc[2 * i6]     = p1 + p2;
            pSrc[2 * i6 + 1] = p3 - p4;
            p1 = co4 * t1;
            p2 = si4 * t2;
            p3 = co4 * t2;
            p4 = si4 * t1;
            pSrc[2 * i4]     = p1 + p2;
            pSrc[2 * i4 + 1] = p3 - p4;
            
            i1 += n1;
         } while(i1 < fftLen);
         
         j++;
      } while(j < n2);
      
      twidCoefModifier <<= 3;
   } while(n2 > 7);   
}
#endif