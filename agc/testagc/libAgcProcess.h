/*
 * libAgcProcess.h
 *
 *  Created on: 2017-9-1
 *      Author: Clarence
 */

#ifndef LIBAGCPROCESS_H_
#define LIBAGCPROCESS_H_


/*
    WnAgcInit, Init function of agc.
 */
extern
void  WnAgcInit(int MaxGainDb,      // Max change of the original vale, 20 is suggested
		        int TargetLevelDb); // Target Db level,  10 mean -10Db

/*
  WnAgcProcess, agc process function
  return value: 0 - Normal  -1 - Error
 */

extern
int WnAgcProcess(short*Input,    // Input data
		         int Inputlen,   // Input data length (sample)
		         short*Output,   // Output data
		         int*Outputlen); // Output data length (sample), always the same as input length



/*Example & Test data:

  short InputAgc[320];
  int Outlen;
  int Inlen = 320;
 short OutAgc[320];

  int Revagc;


void AgcTest(){

    WnAgcInit(20, 11);


    while(1){

        memcpy(InputAgc, alltab, 640   );
	    // nAgcRet = WebRtcAgc_Process(agcHandle, Input, NULL, 160, Outagc,NULL, inMicLevel, &outMicLevel, 0, &saturationWarning);

	     Revagc = WnAgcProcess(InputAgc, Inlen, OutAgc, &Outlen);}

    }


const
short alltab[320] = {
    10034,     4383,     1309,     1162,      148,      105,     -388,     -411,     -624,     -622,     -638,
    -609,     -474,     -434,     -190,     -163,      150,      129,      476,      368,      715,
     481,      811,      409,      707,       79,      356,     -617,     -367,    -2148,    -2293,
  -11863,    10035,     4384,     1309,     1162,      149,      104,     -388,     -412,     -624,
    -622,     -637,     -608,     -474,     -433,     -189,     -161,      151,      129,      474,
     368,      717,      481,      809,      407,      709,       81,      357,     -617,     -369,
   -2149,    -2293,   -11863,    10034,     4383,     1309,     1162,      150,      106,     -386,
    -410,     -625,     -622,     -636,     -608,     -476,     -435,     -190,     -162,      150,
     131,      476,      367,      717,      481,      810,      407,      709,       80,      357,
    -616,     -368,    -2148,    -2294,   -11863,    10034,     4384,     1309,     1163,      149,
     106,     -387,     -412,     -626,     -622,     -638,     -610,     -474,     -433,     -191,
    -162,      150,      129,      476,      368,      716,      483,      811,      408,      709,
      80,      358,     -618,     -369,    -2150,    -2292,   -11863,    10034,     4384,     1309,
    1161,      150,      105,     -388,     -411,     -624,     -622,     -636,     -610,     -474,
    -433,     -190,     -162,      149,      129,      475,      368,      715,      480,      810,
     407,      708,       80,      356,     -618,     -367,    -2150,    -2294,   -11863,
     10034,     4383,     1309,     1162,      148,      105,     -388,     -411,     -624,     -622,     -638,
        -609,     -474,     -434,     -190,     -163,      150,      129,      476,      368,      715,
         481,      811,      409,      707,       79,      356,     -617,     -367,    -2148,    -2293,
      -11863,    10035,     4384,     1309,     1162,      149,      104,     -388,     -412,     -624,
        -622,     -637,     -608,     -474,     -433,     -189,     -161,      151,      129,      474,
         368,      717,      481,      809,      407,      709,       81,      357,     -617,     -369,
       -2149,    -2293,   -11863,    10034,     4383,     1309,     1162,      150,      106,     -386,
        -410,     -625,     -622,     -636,     -608,     -476,     -435,     -190,     -162,      150,
         131,      476,      367,      717,      481,      810,      407,      709,       80,      357,
        -616,     -368,    -2148,    -2294,   -11863,    10034,     4384,     1309,     1163,      149,
         106,     -387,     -412,     -626,     -622,     -638,     -610,     -474,     -433,     -191,
        -162,      150,      129,      476,      368,      716,      483,      811,      408,      709,
          80,      358,     -618,     -369,    -2150,    -2292,   -11863,    10034,     4384,     1309,
        1161,      150,      105,     -388,     -411,     -624,     -622,     -636,     -610,     -474,
        -433,     -190,     -162,      149,      129,      475,      368,      715,      480,      810,
         407,      708,       80,      356,     -618,     -367,    -2150,    -2294,   -11863,

  };
 * */

#endif /* LIBAGCPROCESS_H_ */