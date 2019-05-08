
extern void PostFilterInit();

extern int PostFilterProcess(short*pIn, int Inlen, short*pOut, int*Outlen);

extern int Dual_Channel_PostFilterProcess(short*pIn, short*pInR, int Inlen, short*pOut, int*Outlen);