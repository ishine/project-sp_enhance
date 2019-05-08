
#include "./include/arm_math.h"
#include "./include/arm_common_tables.h"

  float32_t twiddleCoef_256_opt[65] = {
    1.000000000f,  // 0
    0.999698819f,  
    0.998795456f,  
    0.997290457f,  
    0.995184727f,  
    0.992479535f,  
    0.989176510f,  
    0.985277642f,  
    0.980785280f,  
    0.975702130f,  
    0.970031253f,  
    0.963776066f,  
    0.956940336f,  
    0.949528181f,  
    0.941544065f,  
    0.932992799f,  
    0.923879533f,  
    0.914209756f,  
    0.903989293f,  
    0.893224301f,  
    0.881921264f,  
    0.870086991f,  
    0.857728610f,  
    0.844853565f,  
    0.831469612f,  
    0.817584813f,  
    0.803207531f,  
    0.788346428f,  
    0.773010453f,  
    0.757208847f,  
    0.740951125f,  
    0.724247083f,  
    0.707106781f,  // 32
    0.689540545f,  
    0.671558955f,  
    0.653172843f,  
    0.634393284f,  
    0.615231591f,  
    0.595699304f,  
    0.575808191f,  
    0.555570233f,  
    0.534997620f,  
    0.514102744f,  
    0.492898192f,  
    0.471396737f,  
    0.449611330f,  
    0.427555093f,  
    0.405241314f,  
    0.382683432f,  
    0.359895037f,  
    0.336889853f,  
    0.313681740f,  
    0.290284677f,  
    0.266712757f,  
    0.242980180f,  
    0.219101240f,  
    0.195090322f,  
    0.170961889f,  
    0.146730474f,  
    0.122410675f,  
    0.098017140f,  
    0.073564564f,  
    0.049067674f,  
    0.024541229f,  
    0.000000000f,  // 64 
	 
};

  uint16_t armBitRevIndexTable256[ARMBITREVINDEXTABLE_256_TABLE_LENGTH] = 
{
   //8x4, size 440
   8,512, 16,1024, 24,1536, 32,64, 40,576, 48,1088, 56,1600, 64,128, 72,640, 
   80,1152, 88,1664, 96,192, 104,704, 112,1216, 120,1728, 128,256, 136,768, 
   144,1280, 152,1792, 160,320, 168,832, 176,1344, 184,1856, 192,384, 
   200,896, 208,1408, 216,1920, 224,448, 232,960, 240,1472, 248,1984, 
   256,512, 264,520, 272,1032, 280,1544, 288,640, 296,584, 304,1096, 312,1608, 
   320,768, 328,648, 336,1160, 344,1672, 352,896, 360,712, 368,1224, 376,1736, 
   384,520, 392,776, 400,1288, 408,1800, 416,648, 424,840, 432,1352, 440,1864, 
   448,776, 456,904, 464,1416, 472,1928, 480,904, 488,968, 496,1480, 504,1992, 
   520,528, 512,1024, 528,1040, 536,1552, 544,1152, 552,592, 560,1104, 
   568,1616, 576,1280, 584,656, 592,1168, 600,1680, 608,1408, 616,720, 
   624,1232, 632,1744, 640,1032, 648,784, 656,1296, 664,1808, 672,1160, 
   680,848, 688,1360, 696,1872, 704,1288, 712,912, 720,1424, 728,1936, 
   736,1416, 744,976, 752,1488, 760,2000, 768,1536, 776,1552, 784,1048, 
   792,1560, 800,1664, 808,1680, 816,1112, 824,1624, 832,1792, 840,1808, 
   848,1176, 856,1688, 864,1920, 872,1936, 880,1240, 888,1752, 896,1544, 
   904,1560, 912,1304, 920,1816, 928,1672, 936,1688, 944,1368, 952,1880, 
   960,1800, 968,1816, 976,1432, 984,1944, 992,1928, 1000,1944, 1008,1496, 
   1016,2008, 1032,1152, 1040,1056, 1048,1568, 1064,1408, 1072,1120, 
   1080,1632, 1088,1536, 1096,1160, 1104,1184, 1112,1696, 1120,1552, 
   1128,1416, 1136,1248, 1144,1760, 1160,1664, 1168,1312, 1176,1824, 
   1184,1544, 1192,1920, 1200,1376, 1208,1888, 1216,1568, 1224,1672, 
   1232,1440, 1240,1952, 1248,1560, 1256,1928, 1264,1504, 1272,2016, 
   1288,1312, 1296,1408, 1304,1576, 1320,1424, 1328,1416, 1336,1640, 
   1344,1792, 1352,1824, 1360,1920, 1368,1704, 1376,1800, 1384,1432, 
   1392,1928, 1400,1768, 1416,1680, 1432,1832, 1440,1576, 1448,1936, 
   1456,1832, 1464,1896, 1472,1808, 1480,1688, 1488,1936, 1496,1960, 
   1504,1816, 1512,1944, 1520,1944, 1528,2024, 1560,1584, 1592,1648, 
   1600,1792, 1608,1920, 1616,1800, 1624,1712, 1632,1808, 1640,1936, 
   1648,1816, 1656,1776, 1672,1696, 1688,1840, 1704,1952, 1712,1928, 
   1720,1904, 1728,1824, 1736,1952, 1744,1832, 1752,1968, 1760,1840, 
   1768,1960, 1776,1944, 1784,2032, 1864,1872, 1848,1944, 1872,1888, 
   1880,1904, 1888,1984, 1896,2000, 1912,2032, 1904,2016, 1976,2032,
   1960,1968, 2008,2032, 1992,2016, 2024,2032
};



  float32_t twiddleCoef_rfft_512_opt[129] = {
    0.000000000f,  // 256
    0.012271538f,  // 255
    0.024541229f,  // 
    0.036807223f,  
    0.049067674f,  
    0.061320736f,  
    0.073564564f,  
    0.085797312f,  
    0.098017140f,  
    0.110222207f,  
    0.122410675f,  
    0.134580709f,  
    0.146730474f,  
    0.158858143f,  
    0.170961889f,  
    0.183039888f,  
    0.195090322f,  
    0.207111376f,  
    0.219101240f,  
    0.231058108f,  
    0.242980180f,  
    0.254865660f,  
    0.266712757f,  
    0.278519689f,  
    0.290284677f,  
    0.302005949f,  
    0.313681740f,  
    0.325310292f,  
    0.336889853f,  
    0.348418680f,  
    0.359895037f,  
    0.371317194f,  
    0.382683432f,  
    0.393992040f,  
    0.405241314f,  
    0.416429560f,  
    0.427555093f,  
    0.438616239f,  
    0.449611330f,  
    0.460538711f,  
    0.471396737f,  
    0.482183772f,  
    0.492898192f,  
    0.503538384f,  
    0.514102744f,  
    0.524589683f,  
    0.534997620f,  
    0.545324988f,  
    0.555570233f,  
    0.565731811f,  
    0.575808191f,  
    0.585797857f,  
    0.595699304f,  
    0.605511041f,  
    0.615231591f,  
    0.624859488f,  
    0.634393284f,  
    0.643831543f,  
    0.653172843f,  
    0.662415778f,  
    0.671558955f,  
    0.680600998f,  
    0.689540545f,  
    0.698376249f,  
    0.707106781f,  
    0.715730825f,  
    0.724247083f,  
    0.732654272f,  
    0.740951125f,  
    0.749136395f,  
    0.757208847f,  
    0.765167266f,  
    0.773010453f,  
    0.780737229f,  
    0.788346428f,  
    0.795836905f,  
    0.803207531f,  
    0.810457198f,  
    0.817584813f,  
    0.824589303f,  
    0.831469612f,  
    0.838224706f,  
    0.844853565f,  
    0.851355193f,  
    0.857728610f,  
    0.863972856f,  
    0.870086991f,  
    0.876070094f,  
    0.881921264f,  
    0.887639620f,  
    0.893224301f,  
    0.898674466f,  
    0.903989293f,  
    0.909167983f,  
    0.914209756f,  
    0.919113852f,  
    0.923879533f,  
    0.928506080f,  
    0.932992799f,  
    0.937339012f,  
    0.941544065f,  
    0.945607325f,  
    0.949528181f,  
    0.953306040f,  
    0.956940336f,  
    0.960430519f,  
    0.963776066f,  
    0.966976471f,  
    0.970031253f,  
    0.972939952f,  
    0.975702130f,  
    0.978317371f,  
    0.980785280f,  
    0.983105487f,  
    0.985277642f,  
    0.987301418f,  
    0.989176510f,  
    0.990902635f,  
    0.992479535f,  
    0.993906970f,  
    0.995184727f,  
    0.996312612f,  
    0.997290457f,  
    0.998118113f,  
    0.998795456f,  
    0.999322385f,  
    0.999698819f,  
    0.999924702f,   
    1.000000000f,  // 0
   
};
 