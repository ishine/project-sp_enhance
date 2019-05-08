/*

3rd ButterWorth High pass Filter:  Wc: 400hz Fs: 16000hz

Function: remove DC offset and noise.    

Author: yang.yang@nanosic.com
Data:  2017-01-14
*/  
  

void WnHighPassFilter400hz(short* In,  //input data   
					  short *Out, //output data
					  int length, //data length
					  int mutelen, // mute sample number,  set 1600, about 0.1s
					  int VoiceStartFlag); // VoiceStartFlag, set "1" as the start of voice record, 


/*  Example for API using:

while(1)
{

   //  voice start flag get
  WnHighPassFilter400hz(input, output, 128, 1600,1);


   //  otherwise
    WnHighPassFilter400hz(input, output, 128, 1600,0);
}

*/