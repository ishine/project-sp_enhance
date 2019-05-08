 #include "signal_processing_library.h"
 #include "gain_control.h"

	void *agcHandleIn  ;
    int minLevel = 0;
	int maxLevel = 255;
	int agcMode  = kAgcModeFixedDigital;// 3;      kAgcModeAdaptiveDigital,   kAgcModeFixedDigital
  
	int inMicLevel= 0;//  = micLevelOut;
	int outMicLevel = 0;
	int saturationWarning;


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

	 *Outputlen = Inputlen;

	 return nAgcRet;
}
