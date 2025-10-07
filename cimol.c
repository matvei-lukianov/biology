/*
===============================================================================
                                                                        CIMOL.C
                                            Input molecules names and fractions
===============================================================================
*/
#include <stdio.h>
#include <string.h>

extern float mfra[];
extern int   nml,ier1;
extern char  mnam[][5];
extern FILE  *fopen(),*finn,*fflog;

void cimol(void)
{
int   i,is,fr;
char  lin[250],dum[250],fnam[10],*tr;

      finn=fopen("cmol.gc","r");  
      if (finn==NULL)
        { fprintf(fflog,"\n*** ERROR IN OPENING cmol.gc\n\n"); 
          ier1=1; return;}

      nml=0;
      while(1)
        {
         tr=fgets(lin,200,finn); if(tr==NULL) break;       // read current line

         if(lin[0]=='#') continue;                                   // comment

         is=sscanf(lin,"%s",dum); if(is==EOF) continue;           // empty line
              
         sscanf(lin,"%s%d",fnam,&fr); nml++;                       // data line

         for (i=0;i<4;i++) mnam[nml][i]=fnam[i];
         mfra[nml]=(float)fr;
        }
      fclose(finn);
}
