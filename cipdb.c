/*
===============================================================================
                                                                        CIPDB.C
                                                     Molecule coordinates input
                                                                     PDB format
===============================================================================
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cdim.h"

int  catom(char *);
int  cresi(char *);

extern int   repr,ier1;
extern char  atnamo[][4],rsnamo[][4];

extern FILE  *fopen(),*fin,*fflog;

void cipdb(char *fn,float *xk,float *yk,float *zk,int *nt,float *irad,
           int *atnumo,int *rsnumo)
{
float xl,yl,zl,rx,ry,rz,ird;
int   i,tri,si,rner,rn,an,nnt,atnb;
char  fnn[1000],linename[110],dum[110],dup[10],sub[10],tatt[10],trest[10];
char  *pdbd;

      fin=fopen(fn,"r");
      if (fin==NULL)                       // PDB file not in current directory
        {                                                // Check PDB directory
         pdbd=getenv("CELLPDB");
         if (pdbd==NULL)
           { fprintf(fflog,"\n*** ERROR IN OPENING %s",fn); ier1=1; return; }
         strcpy(fnn,pdbd); strcat(fnn,"/"); strcat(fnn,fn);
         fin=fopen(fnn,"r");
         if (fin==NULL) 
           { fprintf(fflog,"\n*** ERROR IN OPENING %s",fnn);ier1=1; return; }
        }
//-------------------------------------------------------------------- Read PDB
      nnt=0; xk[0]=yk[0]=zk[0]=0.; ird=0.;
      while(1)                                                  
        {
         tri=fscanf(fin,"%s",linename); if (tri==EOF) break;
         if ((si=strcmp(linename,"ATOM"))==0)
           {
            fscanf(fin,"%d",&atnb); fgets(dum,2,fin);    // atom number
            fgets(tatt, 5,fin);                          // atom name    
            fgets(trest,5,fin);                          // residue name   
            fgets(sub,3,fin);                            // subunit
            fscanf(fin,"%d",&rner); fgets(dup,3,fin);    // resude number
            fscanf(fin,"%f%f%f",&xl,&yl,&zl);            // coordinates
            fgets(dum,80,fin);

            if ((rn=cresi(trest))==0)        continue;   // skip double atom
            if ((dup[0]=='2')||(dup[0]=='B'))continue;   // skip double residue
            if ((an=catom(tatt))==0)         continue;   // skip hydrogen   
            if ((repr==2)&&(an!=11))         continue;   // skip non c-alpha

            nnt++;                                       // ACCEPT 
            if (nnt>(DIMA-2))                            // check max.no.
              { fprintf(fflog,"\n*** ERROR: too many atoms in %s\n",fn);
                ier1=1; return; }

//---------------------------------------------------------- Assign output data
            
            atnumo[nnt]=atnb; rsnumo[nnt]=rner;
            atnamo[nnt][0]=tatt [0];atnamo[nnt][1]=tatt [1];
            atnamo[nnt][2]=tatt [2];atnamo[nnt][3]=tatt [3];
            rsnamo[nnt][0]=trest[0];rsnamo[nnt][1]=trest[1];
            rsnamo[nnt][2]=trest[2];rsnamo[nnt][3]=trest[3];
            
            xk[nnt]=xl; yk[nnt]=yl; zk[nnt]=zl;
            xk[0] +=xl; yk[0] +=yl; zk[0] +=zl;
           }
         else fgets(dum,100,fin);
        }
      xk[0]/=(float)nnt; yk[0]/=(float)nnt; zk[0]/=(float)nnt;// center of mass
      
      for (i=1; i<=nnt; i++)                                  // radius
        {
         rx=xk[i]-xk[0]; ry=yk[i]-yk[0]; rz=zk[i]-zk[0];
         rx=rx*rx; ry=ry*ry; rz=rz*rz;
         ird+=(rx+ry+rz);
        }
      ird/=(float)nnt; *irad=(float)sqrt((double)ird); *nt=nnt;
      
      fclose(fin);
 }
