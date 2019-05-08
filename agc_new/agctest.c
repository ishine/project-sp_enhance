 #include "signal_processing_library.h"
 #include "gain_control.h"



//enum
//{
//    kAgcModeUnchanged,
//    kAgcModeAdaptiveAnalog,
//    kAgcModeAdaptiveDigital,
//    kAgcModeFixedDigital
//};

	void *agcHandleIn  ;
		int minLevel = 0;
	int maxLevel = 255;
	int agcMode  = kAgcModeAdaptiveAnalog;  // kAgcModeAdaptiveAnalog

	
	WebRtcAgc_config_t agcConfig;

void* WnAgcInit(){
	 
    WebRtcAgc_Create(&agcHandleIn);


	WebRtcAgc_Init(agcHandleIn, minLevel, maxLevel, agcMode, 16000);

	agcConfig.compressionGaindB = 30;
	agcConfig.limiterEnable     = 1;
	agcConfig.targetLevelDbfs   = 3;
	WebRtcAgc_set_config(agcHandleIn, agcConfig);


	return agcHandleIn;
}