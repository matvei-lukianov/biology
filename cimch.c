/*
===============================================================================
                                                                        CIMCH.C
                                                                      RES input
===============================================================================
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cdim.h"

extern float Emch[][DIMN+1][DIMH+1];
extern float dx[][DIMN+1][DIMH+1],dy[][DIMN+1][DIMH+1],dz[][DIMN+1][DIMH+1];
extern float d11[][DIMN+1][DIMH+1],d12[][DIMN+1][DIMH+1],d13[][DIMN+1][DIMH+1];
extern float d21[][DIMN+1][DIMH+1],d22[][DIMN+1][DIMH+1],d23[][DIMN+1][DIMH+1];
extern float d31[][DIMN+1][DIMH+1],d32[][DIMN+1][DIMH+1],d33[][DIMN+1][DIMH+1];
extern int   mtch,repr,ier1,ier2;
extern FILE  *fopen(),*finm,*fflog;

void cimch(char *fmch,int im,int jm)
{
float en,a,b,g,u,v,w,sa,sb,sg,ca,cb,cg;
int   tri1,trii,mi,cali,nm;
char  fnn[1000],edr[910],cal[900],dum[910];
char  *resd;

      finm=fopen(fmch,"r");
      if (finm==NULL)                      // RES file not in current directory
        {                                                // Check RES directory
         resd=getenv("CELLRES");
         if (resd==NULL)
           { fprintf(fflog,"\n*** ERROR IN OPENING %s",fmch);ier1=1; return; }
         strcpy(fnn,resd); strcat(fnn,"/"); strcat(fnn,fmch);
         finm=fopen(fnn,"r");
         if (finm==NULL) 
           { fprintf(fflog,"\n*** ERROR IN OPENING %s",fnn); ier1=1; return; }
        } 
//---------------------------------------------- Search for [molecules] section
      nm=0;
      while(1)
        {
         tri1=fscanf(finm,"%s",edr); trii=strcmp(edr,"[molecules]");
         if (tri1==EOF)
           {
            fprintf(fflog,"\n*** ERROR:[molecules] section not found");
            fprintf(fflog," in %s",fmch); ier2=1; return;
           }
         if (trii==0) { fgets(dum,100,finm); break; }
        }
      fscanf(finm,"%s%s%s%s%s%s%s",dum,dum,dum,dum,dum,dum,cal);
      fgets(dum,100,finm);
      fscanf(finm,"%s%s%s%s%s%s%s",dum,dum,dum,dum,dum,dum,cal);
      fgets(dum,100,finm);
      if ((cali=strcmp(cal,"CA"))==0) repr=2; else repr=0;
      
//-------------------------------------------------- Search for [match] section
      while(1)
        {
         tri1=fscanf(finm,"%s",edr); trii=strcmp(edr,"[match]");
         if (tri1==EOF)
           {fprintf(fflog,"\n*** ERROR: [match] section not found\n\n");
            ier1=1; return;}
         if (trii==0) { fgets(dum,100,finm); break; }
        }
//--------------------------------------------------------------- Input matches
      while(1)
        {
         tri1=fscanf(finm,"%d%f%f%f%f%f%f%f",&mi,&en,&a,&b,&g,&u,&v,&w);
         fgets(dum,100,finm); if(tri1==EOF)break;
         
         nm++; Emch[im][jm][nm]=-en;                       // E with minus sign

         dx[im][jm][nm]=u; dy[im][jm][nm]=v; dz[im][jm][nm]=w;

         a*=0.017453; b*=0.017453; g*=0.017453;
         sa=(float)sin((double)a); ca=(float)cos((double)a);
         sb=(float)sin((double)b); cb=(float)cos((double)b);
         sg=(float)sin((double)g); cg=(float)cos((double)g);
         d11[im][jm][nm]=cg*ca;         d12[im][jm][nm]=cg*sa;         d13[im][jm][nm]=-sg;
         d21[im][jm][nm]=sb*sg*ca-cb*sa;d22[im][jm][nm]=sb*sg*sa+cb*ca;d23[im][jm][nm]=cg*sb;
         d31[im][jm][nm]=cb*sg*ca+sb*sa;d32[im][jm][nm]=cb*sg*sa-sb*ca;d33[im][jm][nm]=cg*cb;
         
         if (nm > DIMH)                                 // Check no. of matches
           {fprintf(fflog,"\n*** ERROR: cannot process > %d matches\n\n",DIMH);
            ier1=1; return;}
        }
      if (nm < mtch)
        { fprintf(fflog,"\n*** ERROR: requested %d matches ",mtch);
          fprintf(fflog,"exceed that in %s\n\n",fmch); ier1=1; return; }
      fclose(finm);
}
