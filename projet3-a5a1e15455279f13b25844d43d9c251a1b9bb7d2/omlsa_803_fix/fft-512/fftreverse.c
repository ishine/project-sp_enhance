
#if 1
	int ibit;
	unsigned int     ta,  tb,  tc,  td;
	unsigned int  pa,pb,pc, pd  , lenb;

void biterverse(unsigned int *In1,  unsigned short In2,   unsigned short *In3){
	
	  unsigned short *ptr0, *ptr1;
	int i;


	ptr1 = &In3[1]; //r1
	lenb = (In2+1)>>2;
	for(ibit = 0; ibit< lenb; ibit ++){
		
		pa = ptr1[2] >>2;
		pb = ptr1[1] >>2;
		pc = ptr1[0] >>2;	
		pd = ptr1[-1]>>2;
		 
		ta =  In1[pa]; 
        tb =  In1[pb]; 		
		tc =  In1[pc];
		td =  In1[pd];
				   
	  In1[pa]  =   tb  ;
	  In1[pb]  =   ta  ;
	  In1[pc]  =   td  ;
	  In1[pd]  =   tc  ;
		
	  ta =  In1[pa+1]; 
	  tb =  In1[pb+1]; 
	  tc =  In1[pc+1];
	  td =  In1[pd+1];
	 	
	  In1[pa+1]  =    tb ;
	  In1[pb+1]  =    ta ;
	  In1[pc+1]  =    td ;
	  In1[pd+1]  =    tc ;	
		 
		ptr1+=4;
	 
	}
}

#else 

int ibit;
	unsigned int    ta, tb, tc, td;
	unsigned int  pa,pb,pc, pd  , lenb;

void biterverse(unsigned int *In1,  unsigned short In2,   unsigned short *In3){
	
	  unsigned short *ptr0, *ptr1;
	int i;


	ptr1 = &In3[1]; //r1
	lenb = (In2+1)>>2;
	for(ibit = 0; ibit< lenb; ibit ++){
		
		pa = ptr1[2]  ;
		pb = ptr1[1]  ;
		pc = ptr1[0]  ;		
		pd = ptr1[-1] ;	
		
		pa =pa/4;
		pb =pb/4;
		pc =pc/4;		
		pd =pd/4;	



		ta =  In1[pa]; 
        tb =  In1[pb]; 		
		tc =  In1[pc];
		td =  In1[pd];
		
	  In1[pa]  =   tb ;
	  In1[pb]  =   ta ;
	  In1[pc]  =   td ;
	  In1[pd]  =   tc ;
		
	 ta = In1[pa+1]; 
	 tb = In1[pb+1]; 
	 tc = In1[pc+1];
	 td = In1[pd+1];
		
	  In1[pa+1]  =   tb ;
	  In1[pb+1]  =   ta ;
	  In1[pc+1]  =   td ;
	  In1[pd+1]  =   tc ;	
		 
		ptr1+=4;
	 
	}
}
#endif