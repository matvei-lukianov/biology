/*
===============================================================================
                                                                        COPER.C
                                                               Operating module
===============================================================================
*/
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "cdim.h"

extern float X [][DIMA+1],Y [][DIMA+1],Z [][DIMA+1];
extern float Xs[][DIMA+1],Ys[][DIMA+1],Zs[][DIMA+1];
extern float Xp0[],Yp0[],Zp0[];
extern float mrad[];
extern int   ncop[],nat[],mnum[];
extern int   atnum[][DIMA+1],rsnum[][DIMA+1];
extern int   cellt,nml,mov,nmol,nstep,nout,stat_flag,wst,scren,ier1,ier2;
extern char  mnam[][5];
extern char  atnam[][DIMA+1][4],rsnam[][DIMA+1][4],atnamo[][4],rsnamo[][4];

void cipar(void);
void cimol(void);
void cipdb(char *,float *,float *,float *,int *,float *,int *,int *);
void ccops(void);
void cimch(char *,int,int);
void cpack(void);
void cpacc(void);
void cmove(void);
void ccout(int);
void cstat(int,int,int);

void coper(void)
{
float xk[DIMA+2],yk[DIMA+2],zk[DIMA+2];
float irad;
int   atnumo[DIMA+2],rsnumo[DIMA+2];
int   i,ii,j,jj,ic,im,nt,k,l,nst,nco,nso,nsw,ns1;
char  fn[10],fmch[15];

      srand(time(0));                       // Initiate random number generator
//------------------------------------------------------------ PARAMETERS INPUT

      cipar(); if(ier1)return;

//------------------------------------------------------------- MOLECULES INPUT
      cimol(); if(ier1)return;                // Read molecule names, fractions

      for (i=1; i<=nml; i++)                                // Read coordinates
        {
         for(ii=0;ii<4;ii++) fn[ii]=mnam[i][ii];
         fn[4]='.';fn[5]='p';fn[6]='d';fn[7]='b';fn[8]='\0';
                                                           
         cipdb(fn,xk,yk,zk,&nt,&irad,atnumo,rsnumo); if(ier1)return;
         nat[i]=nt; mrad[i]=irad;
printf("\n[coper] nat[%d]=%d",i,nat[i]);
         for (k=0; k<=nt; k++)
           {atnum[i][k]=atnumo[k]; rsnum[i][k]=rsnumo[k];
            atnam[i][k][0]=atnamo[k][0];atnam[i][k][1]=atnamo[k][1];
            atnam[i][k][2]=atnamo[k][2];atnam[i][k][3]=atnamo[k][3];
            rsnam[i][k][0]=rsnamo[k][0];rsnam[i][k][1]=rsnamo[k][1];
            rsnam[i][k][2]=rsnamo[k][2];rsnam[i][k][3]=rsnamo[k][3];
            Xs[i][k]=xk[k];Ys[i][k]=yk[k];Zs[i][k]=zk[k];}
         }
      ccops(); if(ier2)return;                    // Calculate number of copies

      im=1;                                     // Assign coordinates to copies
      for (i=1; i<=nml; i++)
        {        
         nt=nat[i];
         for (ic=0; ic<ncop[i]; ic++)
           {
            nmol=im+ic; mnum[nmol]=i;   // Number of all molecules, incl.copies
            for (k=0; k<=nt; k++)
              { X[nmol][k]=Xs[i][k];Y[nmol][k]=Ys[i][k];Z[nmol][k]=Zs[i][k]; }
           }
         im+=ncop[i];
        }
//--------------------------------------------------------------------- PACKING
      if (cellt==1) cpack();
      if (cellt==2) cpacc();
      if (mov==0) {ccout(0); return;}         // If packing only, output & exit

//--------------------------------------------------- LANDSCAPE/RES FILES INPUT
      for (i=1; i<=nml; i++)
        { for(ii=0;ii<4;ii++) fmch[ii]=mnam[i][ii]; fmch[4]='-';
          for (j=1; j<=nml; j++) 
            { for(jj=5;jj<9;jj++) fmch[jj]=mnam[j][jj-5]; fmch[9]='.';
              fmch[10]='r';fmch[11]='e';fmch[12]='s';fmch[13]='\0';
              cimch(fmch,i,j); if(ier1||ier2)return; }}
//------------------------------------------------------------------------- RUN
                                                        // CM corrected for PBC
      for (l=1;l<=nmol;l++) { Xp0[l]=X[l][0];Yp0[l]=Y[l][0];Zp0[l]=Z[l][0]; }
      nst=nco=nso=nsw=ns1=0;
      
      while (nst<nstep) 
        {
         cmove(); nst++;                                         // *** M O V E
         
         ns1++; if (ns1==scren) {printf("\n----------\nSTEP %2d ",nst); ns1=0;}
         nco++; if (nco==nout) { ccout(nst); nco=0; }        // Coordinates out
         nso++;nsw++; cstat(nst,nso,nsw); if(nso==stat_flag) nso=0;    // Stats  out
                                          if(nsw== wst) nsw=0;
        }
}
