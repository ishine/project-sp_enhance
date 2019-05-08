



#include <stdio.h>

//#include <conio.h>
//#include "main.h"
#include <math.h>
#include "main.h"

 

#define RAW_FRAME_LENGTH	 (320) * 2		 

__int64 testtype;

unsigned char wavhadtab[]={
0x52,0x49,0x46,0x46,0xac,0xae,0x28,0x00,0x57,0x41,0x56,0x45,0x66,0x6d,0x74,0x20,
0x10,0x00,0x00,0x00,0x01,0x00,0x02,0x00,

0x44,0xac,
0x00,0x00,0x10,0xb1,0x02,0x00,
0x04,0x00,0x10,0x00,0x64,0x61,0x74,0x61,0x88,0xae,0x28,0x00
};





short wav500DC[320] = {
	500,
		13040,23669,30772,33267,30772,23669,13039,501,-12039,-22670,-29773,-32267,-29773,-22670,-12039,499,13039,23669,30773,33267,30773,23670,13040,501,-12040,-22670,-29773,-32267,-29774,-22670,-12039,500,
		13039,23670,30773,33267,30773,23669,13039,499,-12040,-22670,-29773,-32267,-29773,-22670,-12040,501,13039,23670,30773,33267,30773,23669,13039,499,-12040,-22670,-29772,-32266,-29772,-22670,-12039,500,
		13039,23669,30774,33267,30772,23670,13040,500,-12040,-22670,-29772,-32267,-29772,-22670,-12040,500,13040,23670,30772,33267,30773,23670,13039,500,-12039,-22670,-29773,-32267,-29774,-22670,-12039,499,
		13039,23669,30773,33267,30773,23670,13039,501,-12039,-22670,-29773,-32267,-29772,-22670,-12039,500,13040,23670,30772,33267,30772,23670,13040,499,-12040,-22670,-29772,-32267,-29773,-22669,-12040,500,
		13039,23670,30773,33267,30773,23671,13039,500,-12039,-22669,-29773,-32267,-29773,-22670,-12039,500,13039,23671,30773,33267,30773,23670,13039,501,-12039,-22669,-29773,-32267,-29772,-22670,-12040,500,
		13039,23670,30773,33267,30773,23669,13039,500,-12039,-22669,-29773,-32267,-29773,-22670,-12039,500,13039,23670,30772,33267,30773,23669,13039,500,-12039,-22669,-29773,-32267,-29773,-22670,-12040,500,
		13039,23669,30773,33267,30773,23670,13039,500,-12039,-22669,-29773,-32267,-29773,-22669,-12039,500,13040,23669,30772,33267,30772,23670,13039,500,-12039,-22670,-29772,-32267,-29772,-22670,-12040,500,
		13039,23671,30773,33267,30772,23670,13039,500,-12040,-22670,-29772,-32267,-29773,-22670,-12040,500,13039,23670,30772,33267,30773,23670,13039,501,-12040,-22669,-29773,-32267,-29773,-22669,-12040,500,
		13039,23669,30773,33267,30772,23670,13039,500,-12040,-22670,-29773,-32266,-29772,-22670,-12039,500,13039,23671,30772,33267,30773,23670,13039,500,-12039,-22670,-29773,-32267,-29772,-22670,-12039,500,
        13040,23670,30773,33267,30773,23670,13040,500,-12040,-22670,-29773,-32267,-29773,-22670,-12039,501,13040,23670,30773,33267,30773,23670,13040,500,-12039,-22669,-29772,-32267,-29772,-22670,-12039,		
}
;       

void test(void);

extern void CoeffInit();
void  SetTestWav(short* wdata,signed short val);

int EncFramNumIndex = 0 ;int encfram,decfram;
 
	short  RawData[RAW_FRAME_LENGTH/2];
 
	char  wrwav[RAW_FRAME_LENGTH]; 

	extern   void WnHighPassFilter(short* In, short *Out , int length );

	
extern	void addzerooffset(short *input, short *output,int inputlen);
extern void CoefInit();

void main(void)
{	 
	FILE * fInput;
	FILE * fSBC;
	FILE * fOutput;
	
	char oval;
	char  name[50];
	
	int encdaleng,decdaleng,backleng;
	int SBC_Fram_leng ;
 
	void * hEnc ;
	int FrameSize;
	signed short teval = 0 ;
 
	printf("SBC 编解码测试\n");
	printf("Please input the file name:");
//	scanf("%s",name);

	if((fInput = fopen("R2.wav","rb")) == NULL)
	{
		printf("Fail to open file\n");
	}
	else{
		printf("%s is open\n",name);
	}

	if(remove("SBC.bin") == -1)   
		perror("Could not delete SBC.bin");   
	else   
		printf("Deleted SBC.bin\n"   );   

	if(remove("Output.wav") == -1)   
		perror("Could not delete Output.wav ");   
	else   
		printf("Deleted Output.wav \n"   );  

 
	hEnc = (void *)0;
//	hEnc = sbc_init( );
  	fOutput = fopen("Output.wav","wb");
 
 
//	CoefInit();
//------------------------------------------------------
// encode
	encfram = 0;
	if(fInput != NULL)
	{
 
		printf("Encoding\n");
		fread(RawData,0x2c,1,fInput);
	 	encdaleng = (int)RawData[0x15]*0x10000 + RawData[0x14];
	//	encdaleng = 3401*512;
		backleng = encdaleng ;
	    wavhadtab[0x28] = backleng &0xff;
		wavhadtab[0x29] = (backleng>>8) &0xff;
		wavhadtab[0x2a] = (backleng>>16) &0xff;
		wavhadtab[0x2b] = (backleng>>24) &0xff;
		fwrite(RawData,0x2C,1,fOutput);
		
		
		while((encdaleng > RAW_FRAME_LENGTH)&&(feof(fInput)==NULL))
		{
			//------------------------------------
			//       读取数据
			fread(RawData,RAW_FRAME_LENGTH,1,fInput);
			EncFramNumIndex++;
			if(EncFramNumIndex==3055)
			{
				EncFramNumIndex = 0 ;
			}
 
			//	addzerooffset(RawData,wrwav, 128);
	
			//	memcpy(RawData , wav500DC,  640);
if(!encfram)
		WnHighPassFilter400hz(RawData,wrwav, RAW_FRAME_LENGTH/2,  640, 1  );
else
	    WnHighPassFilter400hz(RawData,wrwav, RAW_FRAME_LENGTH/2,  640, 0  );


			fwrite(wrwav,RAW_FRAME_LENGTH,1,fOutput);
 
			encfram ++ ;
			encdaleng -= RAW_FRAME_LENGTH ;
			printf(".");
		}
	}


  	{
		fclose(fOutput);		
	}
	printf("\n");
	if(fInput != NULL)
	{
		fclose(fInput);
		printf("SBC Encode succeed !\n");
	}
 

 

 
 

	
}


//=========================================================================

int  T_sbc_proto_8_80m1[] = {
	0xff7c272c,0xfcb02620,0xda612700, 0xfcb02620,
	0xff7c272c,0xff762170,0xfdbb828c, 0xdac7bb40,
	0xfc1417b8,0xff8b1a31,0xff7d4914, 0xff405e01,
	0xdbf79400,0xfbd8f358,0xff9f3e17, 0xff960e94,
	0x0142291c,0xdde26200,0xfbedadc0, 0xffb54b3b,
	0xffc4e05c,0x03bf7948,0xe071bc00, 0xfc3fbb68,
	0xffca00ed,0x000bb7db,0x06af2308, 0xe3889d20,
	0xfcbc98e8,0xffdba705,0x006c1de4, 0x0a00d410,
	0xe7054ca0,0xfd52986c,0xffe9811d, 0x00e530da,
	0x0d9daee0,0xeac182c0,0xfdf1c8d4, 0xfff5bd1a
};
void test1(void)
{
	unsigned int i,j,temp;
	__int64 templ;
	FILE * fTest;

	fTest = fopen("test.txt","w");
	fprintf(fTest,"SBC temp test:\n");
	j=0;
	for(i=0;i<40;i++)
	{
		j++;
		fprintf(fTest,"0x%8x ",T_sbc_proto_8_80m1[i]);
		if(j==4){
			j=0;
			fprintf(fTest,"\n");
		}	
	}
	fprintf(fTest,"\n");
	for(i=0;i<40;i++)
	{
		j++;
		templ = (__int64)T_sbc_proto_8_80m1[i];
		templ *= -2 ;
		temp = (unsigned int)templ;
		fprintf(fTest,"0x%8x	",temp);
		if(j==4){
			j=0;
			fprintf(fTest,"\n");
		}	
	}
	fclose(fTest);
}

int  T_sbc_proto_8_80m0[] = {
	0x00000000,0xfe8d1970,0xee979f00,0x11686100,
	0x0172e690,0xfff5bd1a,0xfdf1c8d4,0xeac182c0,
	0x0d9daee0,0x00e530da,0xffe9811d,0xfd52986c,
	0xe7054ca0,0x0a00d410,0x006c1de4,0xffdba705,
	0xfcbc98e8,0xe3889d20,0x06af2308,0x000bb7db,
	0xffca00ed,0xfc3fbb68,0xe071bc00,0x03bf7948,
	0xffc4e05c,0xffb54b3b,0xfbedadc0,0xdde26200,
	0x0142291c,0xff960e94,0xff9f3e17,0xfbd8f358,
	0xdbf79400,0xff405e01,0xff7d4914,0xff8b1a31,
	0xfc1417b8,0xdac7bb40,0xfdbb828c,0xff762170
};

void test2(void)
{
	unsigned int i,j,temp;
	__int64 templ;
	FILE * fTest;

	fTest = fopen("test.txt","w");
	fprintf(fTest,"SBC temp test:\n");
	j=0;
	for(i=0;i<40;i++)
	{
		j++;
		fprintf(fTest,"0x%8x ",T_sbc_proto_8_80m0[i]);
		if(j==4){
			j=0;
			fprintf(fTest,"\n");
		}	
	}
	fprintf(fTest,"\n");
	for(i=0;i<40;i++)
	{
		j++;
		templ = (__int64)T_sbc_proto_8_80m0[i];
		templ *= -2 ;
		temp = (unsigned int)templ;
		fprintf(fTest,"0x%8x	",temp);
		if(j==4){
			j=0;
			fprintf(fTest,"\n");
		}	
	}
	fclose(fTest);
}

void test(void)
{
	int i,j,k ;
	unsigned int temp;
	double templ;
	FILE * fTest;
	double val;
	int tempval;

	fTest = fopen("test.txt","w");
	fprintf(fTest,"SBC temp test:\n");
	j=0;
	
	for(k=0;k<16;k++)
	{
		for(i=0;i<8;i++)
		{
			val = ((double)i+0.5)*(k+4)*3.14159265/8 ;
			templ = cos(val);
			tempval = (int)templ*0x8000000 ;
			temp = (unsigned int)tempval;
			fprintf(fTest,"0x%8x	",temp);
			j++;
			if(j==4)
			{
				j=0;
				fprintf(fTest,"\n");
			}	
		}
	//	fprintf(fTest,"\n");
	}
	fclose(fTest);
}


