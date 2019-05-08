 #include "signal_processing_library.h"
 #include "gain_control.h"

#include<stdio.h>
//   enum
//   {
//       kAgcModeUnchanged,
//       kAgcModeAdaptiveAnalog,
//       kAgcModeAdaptiveDigital,
//       kAgcModeFixedDigital
//   };



	void *agcHandleIn  ;
    int minLevel = 0;
	int maxLevel = 255;
	int agcMode  =    kAgcModeFixedDigital  ;//     kAgcModeUnchanged,  kAgcModeAdaptiveAnalog, kAgcModeAdaptiveDigital,  kAgcModeFixedDigital
	int inMicLevel= 0;//  = micLevelOut;
	int outMicLevel = 0;
	unsigned char saturationWarning;


	WebRtcAgc_config_t agcConfig;

void WnAgcInit(int MaxGainDb, int TargetLevelDb){
	 
    WebRtcAgc_Create(&agcHandleIn);

	WebRtcAgc_Init(agcHandleIn, minLevel, maxLevel, agcMode, 16000);

	agcConfig.compressionGaindB = MaxGainDb;
	agcConfig.limiterEnable     = 1;
	agcConfig.targetLevelDbfs   = TargetLevelDb;

	WebRtcAgc_set_config(agcHandleIn, agcConfig);
}
/*
    WnAgcProcess
 */
int WnAgcProcess(short*Input,int Inputlen, short*Output, int*Outputlen){

	 int nAgcRet;

	 nAgcRet = WebRtcAgc_Process(agcHandleIn, Input, 0, Inputlen, Output,0, inMicLevel, &outMicLevel, 0, &saturationWarning);
     //nAgcRet = WebRtcAgc_Process2(agcHandleIn, Input, 0, Inputlen, Output,0, inMicLevel, &outMicLevel, 0, &saturationWarning);

	 printf("outMicLevel is %d  \n ", outMicLevel);

	 *Outputlen = Inputlen;

	 return nAgcRet;
}
