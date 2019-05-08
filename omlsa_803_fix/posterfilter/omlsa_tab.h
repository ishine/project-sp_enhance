
#include "omlsa_def.h"
#ifdef TAB_FLASH 
  const omlsa_float32_t  dv1[257];
  const omlsa_float32_t  dv2[257];
  const omlsa_float32_t  Window[FRAME_LEN];
  const omlsa_float32_t  dv3[257];
  const omlsa_float32_t  dv4[257];
#else 
  omlsa_float32_t   dv1[257];
  omlsa_float32_t   dv2[257];
  omlsa_float32_t   Window[FRAME_LEN];
  omlsa_float32_t   dv3[257];
  omlsa_float32_t   dv4[257];
#endif
 