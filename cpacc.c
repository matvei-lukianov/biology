/*
===============================================================================
                                                                        CPACC.C
                                                    Molecules packing in sphere
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

void cpacc(void)
{
float R,R2,X0,Y0,Z0,dist,dX,dY,dZ,xl,yl,zl;
float an,sa,sb,sg,ca,cb,cg,a11,a12,a13,a21,a22,a23,a31,a32,a33;
int   l,i,mn,nt;

// Random place in shpere, coord origin 0,0,0, sphere center = cell size / 2 
// Collis. don't matter since every molecule has to move to a dimer in 1st st

      R=csize/2.; R2=R*R;
      
      for(l=1; l<=nmol; l++)
        {
         E[l]=0;                                              // All energies 0
         mn=mnum[l]; nt=nat[mn];
         while(1)
           {
            X0=(float)(rand()%(int)csize);
            Y0=(float)(rand()%(int)csize);
            Z0=(float)(rand()%(int)csize);
         
            dist=(X0-R)*(X0-R)+(Y0-R)*(Y0-R)+(Z0-R)*(Z0-R);
            if (dist<R2) break;                                    // In bounds
           }
         dX=X0-Xs[mn][0]; dY=Y0-Ys[mn][0]; dZ=Z0-Zs[mn][0];
         
         for (i=0; i<=nt; i++) { X[l][i]=Xs[mn][i]+dX;
                                 Y[l][i]=Ys[mn][i]+dY;
                                 Z[l][i]=Zs[mn][i]+dZ; }
//------------------------------------------------------------- Random rotation
         an=(float)(rand()%360+1); an*=0.017453;                // Random angle

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
            xl=X[l][i]-X0; yl=Y[l][i]-Y0; zl=Z[l][i]-Z0;
 
            X[l][i]=a11*xl+a12*yl+a13*zl;
            Y[l][i]=a21*xl+a22*yl+a23*zl;
            Z[l][i]=a31*xl+a32*yl+a33*zl;
 
            X[l][i]+=X0; Y[l][i]+=Y0; Z[l][i]+=Z0;
           }
        }
}
