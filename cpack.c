/*
===============================================================================
                                                                        CPACK.C
                                              Molecules packing on cubical grid
===============================================================================
*/
#include <stdio.h>  // TEMP
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cdim.h"

extern float X [][DIMA+1],Y [][DIMA+1],Z [][DIMA+1];
extern float Xs[][DIMA+1],Ys[][DIMA+1],Zs[][DIMA+1];
extern float U11[],U12[],U13[],U21[],U22[],U23[],U31[],U32[],U33[];
extern float csize;
extern int   mnum[],nat[],E[];
extern int   nmol;

void cpack(void)
{
float st,st2,X0,Y0,Z0,dx0,dy0,dz0,xl,yl,zl,d;
float an,sa,sb,sg,ca,cb,cg,a11,a12,a13,a21,a22,a23,a31,a32,a33;
int   id[DIMM+2];
int   l,i,k,mn,nt,nmt,all;

// Regular grid, corner 0,0,0, gridstep = cell size / no. of mol. 
// Collis. don't matter since every molecule has to move to a dimer in 1st st

      st=csize*csize*csize/(float)nmol; 
      st=(float)pow((double)st,(double)(1.0/3.0));
      for(l=1;l<=nmol;l++) E[l]=0;                            // All energies 0
      all=0;                                 // Trigger to all molecules packed
      
      while(1)
        {
printf("\n[cpack] st=%3.0f",st);
         st2=st/2.; for (k=1;k<=nmol;k++) id[k]=0;
         for (Z0=st2; Z0<csize; Z0+=st)
           {
            for (Y0=st2; Y0<csize; Y0+=st)
              {
               for (X0=st2; X0<csize; X0+=st)
                 {                                   // Check if l already used
                  while(1) { l=rand()%nmol+1; if(id[l]==0) break; }
//printf("%3d",l);
                  id[l]=1; mn=mnum[l]; nt=nat[mn];
                  d=(float)(rand()%(int)st); d=st2-d;           // Random shift
                  dz0=Z0-Zs[mn][0]; dy0=Y0-Ys[mn][0]; dx0=X0-Xs[mn][0];
                  for (i=0;i<=nt;i++) { Z[l][i]=Zs[mn][i]+dz0+d;
                                        Y[l][i]=Ys[mn][i]+dy0+d;
                                        X[l][i]=Xs[mn][i]+dx0+d; }
//------------------------------------------------------------- Random rotation
                  an=(float)(rand()%360+1); an*=0.017453;       // Random angle

                  sa=(float)sin((double)an); ca=(float)cos((double)an);
                  sb=(float)sin((double)an); cb=(float)cos((double)an);
                  sg=(float)sin((double)an); cg=(float)cos((double)an);
 
                  a11=cg*ca;          a12=cg*sa;           a13=-sg;
                  a21=sb*sg*ca-cb*sa; a22=sb*sg*sa+cb*ca;  a23= cg*sb;
                  a31=cb*sg*ca+sb*sa; a32=cb*sg*sa-sb*ca;  a33= cg*cb;
                  
                  U11[l]=a11; U12[l]=a12; U13[l]=a13;
                  U21[l]=a21; U22[l]=a22; U23[l]=a23;
                  U31[l]=a31; U32[l]=a32; U33[l]=a33;
 
                  for (i=1; i<=nt; i++)
                    {
                     xl=X[l][i]-X[l][0];yl=Y[l][i]-Y[l][0];zl=Z[l][i]-Z[l][0];
 
                     X[l][i]=a11*xl+a12*yl+a13*zl;
                     Y[l][i]=a21*xl+a22*yl+a23*zl;
                     Z[l][i]=a31*xl+a32*yl+a33*zl;
 
                     X[l][i]+=X[l][0]; Y[l][i]+=Y[l][0]; Z[l][i]+=Z[l][0];
                    }
 //-----------------------------------------------------------------------------
                  nmt=0; for(k=1;k<=nmol;k++) nmt+=id[k];        // Control sum
                  if(nmt>=nmol) { all=1; break; }       // All molecules packed
                 } if(all)break;
              } if(all)break;
           } if(all)break;
         st=0.97*st;              // Not all molecules packed - shrink gridstep
        }
}
