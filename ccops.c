/*
===============================================================================
                                                                        CCOPS.C
                                            Calculate number of molecule copies
===============================================================================
*/
#include <stdio.h>
#include "cdim.h"

float cmvol(int);

extern float mfra[];
extern float vfrac,csize;
extern int   ncop[];
extern int   cellt,nml,ier2;
extern FILE  *fflog;

void ccops(void)
{
float mvol,cvol,mvo,rati;
int   i;

      mvo=0.;
      for(i=1;i<=nml;i++)
        {         
         mvol=cmvol(i);                                   // Volume of molecule
printf("\nprot %d  vol %6.0f",i,mvol);
         mvo+=(mfra[i]*mvol);
        }
      if (cellt==1) cvol=csize*csize*csize;
      if (cellt==2) cvol=0.52*csize*csize*csize;
      
      rati=(cvol*vfrac)/mvo;
      
      for(i=1;i<=nml;i++)
        {
         ncop[i]=(int)(mfra[i]*rati);
printf("\n[ccops] ncop[%d]=%d",i,ncop[i]);
         if (ncop[i]==0)
           { fprintf(fflog,"\n*** ERROR: cell %2.0f and fraction ",csize);
             fprintf(fflog,"%3.1f too small to fit molecules\n",vfrac);
             ier2=1; return; }
         if (ncop[i]>DIMC)
           { fprintf(fflog,"\n*** ERROR: no. of molecule copies %5d",ncop[i]);
             fprintf(fflog," exceeds array DIMC %5d\n",DIMC);
             ier2=1; return; }
        }
}
