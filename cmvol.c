/*
===============================================================================
                                                                        CMVOL.C
                                          Molecule volume by projection on grid
===============================================================================
*/
#include "cdim.h"

extern float Xs[][DIMA+1],Ys[][DIMA+1],Zs[][DIMA+1];
extern int   nat[];

float cmvol(int im)
{
float X[DIMA+2],Y[DIMA+2],Z[DIMA+2];
float vol,v,step,rada,dist,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,xi,yi,zi,xii,yii,zii;
int   nt,i;

      rada=2.65;     // Radius of atom (fit to match 3V Volume Assessor values)
      step=1.;                                                     // Grid step

      nt=nat[im];
      Xmin=Ymin=Zmin=Xmax=Ymax=Zmax=0.; vol=0.;
      rada=rada*rada;
      
      for (i=1; i<=nt; i++)
        {
         X[i]=Xs[im][i]; Y[i]=Ys[im][i]; Z[i]=Zs[im][i];      // Box boundaries
         if(X[i]<Xmin)Xmin=X[i];if(Y[i]<Ymin)Ymin=Y[i];if(Z[i]<Zmin)Zmin=Z[i];
         if(X[i]>Xmax)Xmax=X[i];if(Y[i]>Ymax)Ymax=Y[i];if(Z[i]>Zmax)Zmax=Z[i];
        }
      for (xi=Xmin; xi<Xmax; xi+=step)
        { for (yi=Ymin; yi<Ymax; yi+=step)
           { for (zi=Zmin; zi<Zmax; zi+=step)
              { for (i=1; i<=nt; i++)
                 {
                  xii=xi-X[i]; yii=yi-Y[i]; zii=zi-Z[i];
                  dist=xii*xii+yii*yii+zii*zii;
                  if (dist<rada) { vol+=1.; break; }
                 }}}}
      v=vol*step*step*step;
      return(v);
}
