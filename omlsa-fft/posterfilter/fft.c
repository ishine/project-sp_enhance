/*
 * File: fft.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 09-Nov-2018 14:14:28
 */

/* Include Files */
 
#include "fft.h"

/* Function Definitions */

/*
 * Arguments    : const omlsa_float32_t x[512]
 *                creal_T y[512]
 * Return Type  : void
 */
void fft(const omlsa_float32_t x[512], creal_T y[512])
{
  int ix;
  int ju;
  int iy;
  int i;
  boolean_T tst;
  omlsa_float32_t temp_re;
  omlsa_float32_t temp_im;
  int iheight;
  int istart;
  int j;
  omlsa_float32_t twid_re;
  static const omlsa_float32_t dv1[257] = { 1.0f, 0.9999247018391445f, 0.99969881869620425f,
    0.99932238458834954f, 0.99879545620517241f, 0.99811811290014918f,
    0.99729045667869021f, 0.996312612182778f, 0.99518472667219693f,
    0.99390697000235606f, 0.99247953459871f, 0.99090263542778f, 0.989176509964781f,
    0.98730141815785843f, 0.98527764238894122f, 0.98310548743121629f,
    0.98078528040323043f, 0.97831737071962765f, 0.97570213003852857f,
    0.97293995220556018f, 0.970031253194544f, 0.96697647104485207f,
    0.96377606579543984f, 0.96043051941556579f, 0.95694033573220882f,
    0.95330604035419386f, 0.94952818059303667f, 0.94560732538052128f,
    0.94154406518302081f, 0.937339011912575f, 0.932992798834739f,
    0.92850608047321559f, 0.92387953251128674f, 0.91911385169005777f,
    0.91420975570353069f, 0.90916798309052238f, 0.90398929312344334f,
    0.89867446569395382f, 0.89322430119551532f, 0.88763962040285393f,
    0.881921264348355f, 0.8760700941954066f, 0.87008699110871146f,
    0.8639728561215867f, 0.85772861000027212f, 0.8513551931052652f,
    0.84485356524970712f, 0.83822470555483808f, 0.83146961230254524f,
    0.82458930278502529f, 0.81758481315158371f, 0.81045719825259477f,
    0.80320753148064494f, 0.79583690460888357f, 0.78834642762660634f,
    0.78073722857209449f, 0.773010453362737f, 0.765167265622459f,
    0.75720884650648457f, 0.74913639452345937f, 0.74095112535495922f,
    0.73265427167241282f, 0.724247082951467f, 0.71573082528381859f,
    0.70710678118654757f, 0.69837624940897292f, 0.68954054473706683f,
    0.680600997795453f, 0.67155895484701833f, 0.66241577759017178f,
    0.65317284295377676f, 0.64383154288979139f, 0.63439328416364549f,
    0.62485948814238634f, 0.61523159058062682f, 0.60551104140432555f,
    0.59569930449243336f, 0.58579785745643886f, 0.57580819141784534f,
    0.56573181078361312f, 0.55557023301960218f, 0.54532498842204646f,
    0.53499761988709715f, 0.524589682678469f, 0.51410274419322166f,
    0.50353838372571758f, 0.49289819222978404f, 0.48218377207912272f,
    0.47139673682599764f, 0.46053871095824f, 0.44961132965460654f,
    0.43861623853852766f, 0.42755509343028208f, 0.41642956009763715f,
    0.40524131400498986f, 0.3939920400610481f, 0.38268343236508978f,
    0.37131719395183749f, 0.35989503653498811f, 0.34841868024943456f,
    0.33688985339222005f, 0.32531029216226293f, 0.31368174039889152f,
    0.30200594931922808f, 0.29028467725446233f, 0.27851968938505306f,
    0.26671275747489837f, 0.25486565960451457f, 0.24298017990326387f,
    0.23105810828067111f, 0.2191012401568698f, 0.20711137619221856f,
    0.19509032201612825f, 0.18303988795514095f, 0.17096188876030122f,
    0.15885814333386145f, 0.14673047445536175f, 0.13458070850712617f,
    0.1224106751992162f, 0.11022220729388306f, 0.0980171403295606f,
    0.0857973123444399f, 0.073564563599667426f, 0.061320736302208578f,
    0.049067674327418015f, 0.036807222941358832f, 0.024541228522912288f,
    0.012271538285719925f, 0.0f, -0.012271538285719925f, -0.024541228522912288f,
    -0.036807222941358832f, -0.049067674327418015f, -0.061320736302208578f,
    -0.073564563599667426f, -0.0857973123444399f, -0.0980171403295606f,
    -0.11022220729388306f, -0.1224106751992162f, -0.13458070850712617f,
    -0.14673047445536175f, -0.15885814333386145f, -0.17096188876030122f,
    -0.18303988795514095f, -0.19509032201612825f, -0.20711137619221856f,
    -0.2191012401568698f, -0.23105810828067111f, -0.24298017990326387f,
    -0.25486565960451457f, -0.26671275747489837f, -0.27851968938505306f,
    -0.29028467725446233f, -0.30200594931922808f, -0.31368174039889152f,
    -0.32531029216226293f, -0.33688985339222005f, -0.34841868024943456f,
    -0.35989503653498811f, -0.37131719395183749f, -0.38268343236508978f,
    -0.3939920400610481f, -0.40524131400498986f, -0.41642956009763715f,
    -0.42755509343028208f, -0.43861623853852766f, -0.44961132965460654f,
    -0.46053871095824f, -0.47139673682599764f, -0.48218377207912272f,
    -0.49289819222978404f, -0.50353838372571758f, -0.51410274419322166f,
    -0.524589682678469f, -0.53499761988709715f, -0.54532498842204646f,
    -0.55557023301960218f, -0.56573181078361312f, -0.57580819141784534f,
    -0.58579785745643886f, -0.59569930449243336f, -0.60551104140432555f,
    -0.61523159058062682f, -0.62485948814238634f, -0.63439328416364549f,
    -0.64383154288979139f, -0.65317284295377676f, -0.66241577759017178f,
    -0.67155895484701833f, -0.680600997795453f, -0.68954054473706683f,
    -0.69837624940897292f, -0.70710678118654757f, -0.71573082528381859f,
    -0.724247082951467f, -0.73265427167241282f, -0.74095112535495922f,
    -0.74913639452345937f, -0.75720884650648457f, -0.765167265622459f,
    -0.773010453362737f, -0.78073722857209449f, -0.78834642762660634f,
    -0.79583690460888357f, -0.80320753148064494f, -0.81045719825259477f,
    -0.81758481315158371f, -0.82458930278502529f, -0.83146961230254524f,
    -0.83822470555483808f, -0.84485356524970712f, -0.8513551931052652f,
    -0.85772861000027212f, -0.8639728561215867f, -0.87008699110871146f,
    -0.8760700941954066f, -0.881921264348355f, -0.88763962040285393f,
    -0.89322430119551532f, -0.89867446569395382f, -0.90398929312344334f,
    -0.90916798309052238f, -0.91420975570353069f, -0.91911385169005777f,
    -0.92387953251128674f, -0.92850608047321559f, -0.932992798834739f,
    -0.937339011912575f, -0.94154406518302081f, -0.94560732538052128f,
    -0.94952818059303667f, -0.95330604035419386f, -0.95694033573220882f,
    -0.96043051941556579f, -0.96377606579543984f, -0.96697647104485207f,
    -0.970031253194544f, -0.97293995220556018f, -0.97570213003852857f,
    -0.97831737071962765f, -0.98078528040323043f, -0.98310548743121629f,
    -0.98527764238894122f, -0.98730141815785843f, -0.989176509964781f,
    -0.99090263542778f, -0.99247953459871f, -0.99390697000235606f,
    -0.99518472667219693f, -0.996312612182778f, -0.99729045667869021f,
    -0.99811811290014918f, -0.99879545620517241f, -0.99932238458834954f,
    -0.99969881869620425f, -0.9999247018391445f, -1.0f };

  omlsa_float32_t twid_im;
  static const omlsa_float32_t dv2[257] = { 0.0f, -0.012271538285719925f,
    -0.024541228522912288f, -0.036807222941358832f, -0.049067674327418015f,
    -0.061320736302208578f, -0.073564563599667426f, -0.0857973123444399f,
    -0.0980171403295606f, -0.11022220729388306f, -0.1224106751992162f,
    -0.13458070850712617f, -0.14673047445536175f, -0.15885814333386145f,
    -0.17096188876030122f, -0.18303988795514095f, -0.19509032201612825f,
    -0.20711137619221856f, -0.2191012401568698f, -0.23105810828067111f,
    -0.24298017990326387f, -0.25486565960451457f, -0.26671275747489837f,
    -0.27851968938505306f, -0.29028467725446233f, -0.30200594931922808f,
    -0.31368174039889152f, -0.32531029216226293f, -0.33688985339222005f,
    -0.34841868024943456f, -0.35989503653498811f, -0.37131719395183749f,
    -0.38268343236508978f, -0.3939920400610481f, -0.40524131400498986f,
    -0.41642956009763715f, -0.42755509343028208f, -0.43861623853852766f,
    -0.44961132965460654f, -0.46053871095824f, -0.47139673682599764f,
    -0.48218377207912272f, -0.49289819222978404f, -0.50353838372571758f,
    -0.51410274419322166f, -0.524589682678469f, -0.53499761988709715f,
    -0.54532498842204646f, -0.55557023301960218f, -0.56573181078361312f,
    -0.57580819141784534f, -0.58579785745643886f, -0.59569930449243336f,
    -0.60551104140432555f, -0.61523159058062682f, -0.62485948814238634f,
    -0.63439328416364549f, -0.64383154288979139f, -0.65317284295377676f,
    -0.66241577759017178f, -0.67155895484701833f, -0.680600997795453f,
    -0.68954054473706683f, -0.69837624940897292f, -0.70710678118654757f,
    -0.71573082528381859f, -0.724247082951467f, -0.73265427167241282f,
    -0.74095112535495922f, -0.74913639452345937f, -0.75720884650648457f,
    -0.765167265622459f, -0.773010453362737f, -0.78073722857209449f,
    -0.78834642762660634f, -0.79583690460888357f, -0.80320753148064494f,
    -0.81045719825259477f, -0.81758481315158371f, -0.82458930278502529f,
    -0.83146961230254524f, -0.83822470555483808f, -0.84485356524970712f,
    -0.8513551931052652f, -0.85772861000027212f, -0.8639728561215867f,
    -0.87008699110871146f, -0.8760700941954066f, -0.881921264348355f,
    -0.88763962040285393f, -0.89322430119551532f, -0.89867446569395382f,
    -0.90398929312344334f, -0.90916798309052238f, -0.91420975570353069f,
    -0.91911385169005777f, -0.92387953251128674f, -0.92850608047321559f,
    -0.932992798834739f, -0.937339011912575f, -0.94154406518302081f,
    -0.94560732538052128f, -0.94952818059303667f, -0.95330604035419386f,
    -0.95694033573220882f, -0.96043051941556579f, -0.96377606579543984f,
    -0.96697647104485207f, -0.970031253194544f, -0.97293995220556018f,
    -0.97570213003852857f, -0.97831737071962765f, -0.98078528040323043f,
    -0.98310548743121629f, -0.98527764238894122f, -0.98730141815785843f,
    -0.989176509964781f, -0.99090263542778f, -0.99247953459871f,
    -0.99390697000235606f, -0.99518472667219693f, -0.996312612182778f,
    -0.99729045667869021f, -0.99811811290014918f, -0.99879545620517241f,
    -0.99932238458834954f, -0.99969881869620425f, -0.9999247018391445f, -1.0f,
    -0.9999247018391445f, -0.99969881869620425f, -0.99932238458834954f,
    -0.99879545620517241f, -0.99811811290014918f, -0.99729045667869021f,
    -0.996312612182778f, -0.99518472667219693f, -0.99390697000235606f,
    -0.99247953459871f, -0.99090263542778f, -0.989176509964781f,
    -0.98730141815785843f, -0.98527764238894122f, -0.98310548743121629f,
    -0.98078528040323043f, -0.97831737071962765f, -0.97570213003852857f,
    -0.97293995220556018f, -0.970031253194544f, -0.96697647104485207f,
    -0.96377606579543984f, -0.96043051941556579f, -0.95694033573220882f,
    -0.95330604035419386f, -0.94952818059303667f, -0.94560732538052128f,
    -0.94154406518302081f, -0.937339011912575f, -0.932992798834739f,
    -0.92850608047321559f, -0.92387953251128674f, -0.91911385169005777f,
    -0.91420975570353069f, -0.90916798309052238f, -0.90398929312344334f,
    -0.89867446569395382f, -0.89322430119551532f, -0.88763962040285393f,
    -0.881921264348355f, -0.8760700941954066f, -0.87008699110871146f,
    -0.8639728561215867f, -0.85772861000027212f, -0.8513551931052652f,
    -0.84485356524970712f, -0.83822470555483808f, -0.83146961230254524f,
    -0.82458930278502529f, -0.81758481315158371f, -0.81045719825259477f,
    -0.80320753148064494f, -0.79583690460888357f, -0.78834642762660634f,
    -0.78073722857209449f, -0.773010453362737f, -0.765167265622459f,
    -0.75720884650648457f, -0.74913639452345937f, -0.74095112535495922f,
    -0.73265427167241282f, -0.724247082951467f, -0.71573082528381859f,
    -0.70710678118654757f, -0.69837624940897292f, -0.68954054473706683f,
    -0.680600997795453f, -0.67155895484701833f, -0.66241577759017178f,
    -0.65317284295377676f, -0.64383154288979139f, -0.63439328416364549f,
    -0.62485948814238634f, -0.61523159058062682f, -0.60551104140432555f,
    -0.59569930449243336f, -0.58579785745643886f, -0.57580819141784534f,
    -0.56573181078361312f, -0.55557023301960218f, -0.54532498842204646f,
    -0.53499761988709715f, -0.524589682678469f, -0.51410274419322166f,
    -0.50353838372571758f, -0.49289819222978404f, -0.48218377207912272f,
    -0.47139673682599764f, -0.46053871095824f, -0.44961132965460654f,
    -0.43861623853852766f, -0.42755509343028208f, -0.41642956009763715f,
    -0.40524131400498986f, -0.3939920400610481f, -0.38268343236508978f,
    -0.37131719395183749f, -0.35989503653498811f, -0.34841868024943456f,
    -0.33688985339222005f, -0.32531029216226293f, -0.31368174039889152f,
    -0.30200594931922808f, -0.29028467725446233f, -0.27851968938505306f,
    -0.26671275747489837f, -0.25486565960451457f, -0.24298017990326387f,
    -0.23105810828067111f, -0.2191012401568698f, -0.20711137619221856f,
    -0.19509032201612825f, -0.18303988795514095f, -0.17096188876030122f,
    -0.15885814333386145f, -0.14673047445536175f, -0.13458070850712617f,
    -0.1224106751992162f, -0.11022220729388306f, -0.0980171403295606f,
    -0.0857973123444399f, -0.073564563599667426f, -0.061320736302208578f,
    -0.049067674327418015f, -0.036807222941358832f, -0.024541228522912288f,
    -0.012271538285719925f, -0.0f };

  int ihi;
  ix = 0;
  ju = 0;
  iy = 0;
  for (i = 0; i < 511; i++) {
    y[iy].re = x[ix];
    y[iy].im = 0.0;
    iy = 512;
    tst = true;
    while (tst) {
      iy >>= 1;
      ju ^= iy;
      tst = ((ju & iy) == 0);
    }

    iy = ju;
    ix++;
  }

  y[iy].re = x[ix];
  y[iy].im = 0.0;
  for (i = 0; i <= 511; i += 2) {
    temp_re = y[i + 1].re;
    temp_im = y[i + 1].im;
    y[i + 1].re = y[i].re - y[i + 1].re;
    y[i + 1].im = y[i].im - y[i + 1].im;
    y[i].re += temp_re;
    y[i].im += temp_im;
  }

  iy = 2;
  ix = 4;
  ju = 128;
  iheight = 509;
  while (ju > 0) {
    for (i = 0; i < iheight; i += ix) {
      temp_re = y[i + iy].re;
      temp_im = y[i + iy].im;
      y[i + iy].re = y[i].re - temp_re;
      y[i + iy].im = y[i].im - temp_im;
      y[i].re += temp_re;
      y[i].im += temp_im;
    }

    istart = 1;
    for (j = ju; j < 256; j += ju) {
      twid_re = dv1[j];
      twid_im = dv2[j];
      i = istart;
      ihi = istart + iheight;
      while (i < ihi) {
        temp_re = twid_re * y[i + iy].re - twid_im * y[i + iy].im;
        temp_im = twid_re * y[i + iy].im + twid_im * y[i + iy].re;
        y[i + iy].re = y[i].re - temp_re;
        y[i + iy].im = y[i].im - temp_im;
        y[i].re += temp_re;
        y[i].im += temp_im;
        i += ix;
      }

      istart++;
    }

    ju /= 2;
    iy = ix;
    ix <<= 1;
    iheight -= iy;
  }
}

/*
 * File trailer for fft.c
 *
 * [EOF]
 */

#include "../TransformFunctions/Include/arm_math.h"
 
 extern  arm_cfft_instance_f32 arm_cfft_sR_f32_len512;

//void arm_cfft_f32( 
//  const arm_cfft_instance_f32 * S, 
//  float32_t * p1,
//  uint8_t ifftFlag,
//  uint8_t bitReverseFlag);

 arm_cfft_radix2_instance_q31 S;

 int fft_init = 0;
 
void fftarm(float * fftin){
  
 
	//arm_cfft_f32(&arm_cfft_sR_f32_len512, fftin,  0,   0);
 
}
void fftarm_fixq31(int * fftin){
	if(fft_init==0){	
//		arm_cfft_radix2_init_q31(&S, 512,0,1 );
	}

//	arm_cfft_radix2_q31(&S, fftin );

}

void rdft(int n, int isgn, float *a);


void offt(float *fftin ){

	int j;

	rdft(512, 1, fftin);
	rdft(512, -1, fftin);
   
	for (j = 0; j <= 512 - 1; j++) {
                fftin[j] *= 2.0 / 512;
            }
}

extern
void arm_rfft_fast_f32(
arm_rfft_fast_instance_f32 * S,
float32_t * p, float32_t * pOut,
uint8_t ifftFlag);

int armfft_init = 0;
 arm_rfft_fast_instance_f32   Sf32;

 float fftoutputram[1024];
 float fftoutputram2[1024];
  

int fftinit_flag;
arm_rfft_fast_instance_f32   Sf31;
void rfftlocal_fast(float* fftinput ,float* fftoutput ){
	int i;
	if(!fftinit_flag){
 	arm_rfft_fast_init_f32(&Sf31,512);
		fftinit_flag = 1;
  }

	
	arm_rfft_fast_f32(&Sf31, fftinput, fftoutput,0);

//	 for(i=0;i< 512;i++){ 	 		 
// 		 fftoutputram[i] = (float)fftoutput[ i];
//   }	

//    arm_rfft_fast_f32(&Sf31, fftoutput, fftinput,1);
//
// for(i=0;i< 512;i++){ 	 		 
// 		 fftoutputram[i] = (float)fftinput[2*i];
//   }	
}

void rifftlocal_fast(float* fftinput ,float* fftoutput ){
	 arm_rfft_fast_f32(&Sf31, fftinput, fftoutput,1);
}