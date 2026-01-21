/*
===============================================================================
                                                                        CSTAT.C
                                                                     Statistics
===============================================================================
*/
#include <stdio.h>
#include <math.h>
#include "cdim.h"

extern float E[],sh[],Xp0[],Yp0[],Zp0[],Xs0[],Ys0[],Zs0[];
extern float X[][DIMA+1],Y[][DIMA+1],Z[][DIMA+1];
extern int   mnum[],ncop[];
extern int   nmol,stat_flag,wst,refst,clus,ndr,pp,Nac;
extern FILE  *fp_stat,*fdrot;

static float Eaw[DIMW+2],Saw[DIMW+2],MSDw[DIMW+2],RMSDw[DIMW+2],Nartw[DIMW+2];
static int   tg,nclu;

// Globals for Interface
float LastEavg=0., LastSavg=0., LastAccAvg=0., LastMSDavg=0.;

void cclus(int);

void cstat(int nst,int nso,int nsw)
{
float Ea,sl,Sa,MSD,RMSD,Nart,Eawa,Sawa,MSDwa,RMSDwa,Nartwa,dx,dy,dz;
int   l,i,ncl;

//-------- Output data for rotational diffusion (supresses stat_flag & clust output)      
      if (ndr)
        {
         if (nst==stat_flag) fprintf(fdrot,"Step     Protein/vector XYZ(o)-XYZ(v)");
         if (nso==stat_flag) fprintf(fdrot,"\n%5d ",nst);
         i=1;
         for (l=1; l<=nmol; l++)
           {              
            if ((mnum[l]==pp)&&(nso==stat_flag))
              {
               fprintf(fdrot,"| %3.0f %3.0f %3.0f - %3.0f %3.0f %3.0f ",
                              X[l][0],Y[l][0],Z[l][0],X[l][1],Y[l][1],Z[l][1]);
               i++; if (i>ndr) break;
              }
           }
         return;
        }
//-----------------------------------------------------------------------------
      Ea=Sa=MSD=0.; Eawa=Sawa=MSDwa=RMSDwa=Nartwa=0.; 

      if(nst==stat_flag){fprintf(fp_stat,  "----------------------------------------");
                    fprintf(fp_stat,"\n Step    E    Shift Nacc     MSD    RMSD");
                    fprintf(fp_stat,"\n----------------------------------------");} 

      if (nst==refst) tg=0;
      for (l=1; l<=nmol; l++) 
        {
         Ea+=E[l];                                                    // Energy
         sl=(float)sqrt((double)sh[l]); Sa+=sl;                       // Shift
         if ((mnum[l]==pp)&&(nst>=refst))                             
           { 
            if (tg==0) { Xs0[l]=Xp0[l]; Ys0[l]=Yp0[l]; Zs0[l]=Zp0[l]; } // MSD
            dx=Xp0[l]-Xs0[l]; dy=Yp0[l]-Ys0[l]; dz=Zp0[l]-Zs0[l];
            MSD+=(dx*dx+dy*dy+dz*dz);
           }
        }
      Ea/=((float)nmol);                                      // Energy average
      Sa/=((float)nmol);                                       // Shift average
      MSD/=((float)ncop[pp]);RMSD=(float)sqrt((double)MSD);tg=1; // MSD average
      Nart=100.*(float)Nac/(float)ncop[pp];                  // Acceptance rate

      if (nst<=wst)                                         // Window averaging
        {Eaw[nsw]=Ea;Saw[nsw]=Sa;MSDw[nsw]=MSD;RMSDw[nsw]=RMSD;Nartw[nsw]=Nart;}
      else 
        {for (i=1; i<wst; i++) {Eaw[i]=  Eaw[i+1];   Saw[i]=  Saw[i+1];
                               MSDw[i]= MSDw[i+1]; RMSDw[i]=RMSDw[i+1];
                              Nartw[i]=Nartw[i+1]; }
         Eaw[wst]=Ea;Saw[wst]=Sa;MSDw[wst]=MSD;RMSDw[wst]=RMSD;Nartw[wst]=Nart;}

      for (i=1; i<=wst; i++)
        {Eawa+=Eaw[i];Sawa+=Saw[i];MSDwa+=MSDw[i];RMSDwa+=RMSDw[i];Nartwa+=Nartw[i];}

      Eawa  /=(float)wst;Sawa/=(float)wst;MSDwa/=(float)wst;RMSDwa/=(float)wst;
      Nartwa/=(float)wst;

      if ((nso==stat_flag)&&(nst>=wst))
        fprintf(fp_stat,"\n%5d %6.2f %5.2f %5.2f %8.0f %6.0f",
                         nst,Eawa, Sawa, Nartwa,MSDwa,RMSDwa);
      
      // Store for Python interface
      LastEavg = Eawa;
      LastSavg = Sawa;
      LastAccAvg = Nartwa;
      LastMSDavg = MSDwa;
      
      if (nso==stat_flag)                                        // Cluster data out
        {
         if (nst==stat_flag) nclu=0; 
         nclu++; ncl=nclu*stat_flag;
         if (ncl==clus) { cclus(nst); nclu=0; }
        }
 }

 // Getter for Python
 void get_latest_stats(float* out_stats) {
     out_stats[0] = LastEavg;
     out_stats[1] = LastSavg;
     out_stats[2] = LastAccAvg;
     out_stats[3] = LastMSDavg;
 }
