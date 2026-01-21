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

// Safely read a PDB line
// ATOM      1  N   ASP     1      -6.204  12.029  14.615  1.00  0.00
// 0123456789012345678901234567890123456789012345678901234567890123456789

void cipdb(char *fn,float *xk,float *yk,float *zk,int *nt,float *irad,
           int *atnumo,int *rsnumo)
{
      float xl,yl,zl,rx,ry,rz,ird;
      int   i,si,rner,rn,an,nnt,atnb;
      char  fnn[1000],line[256];
      char  sub_str[16]; // buffer for substrings
      char  *pdbd;
      
      fin=fopen(fn,"r");
      if (fin==NULL) {
         pdbd=getenv("CELLPDB");
         if (pdbd==NULL) { 
             fprintf(stderr, "\n*** ERROR IN OPENING %s (CELLPDB not set)\n",fn); 
             ier1=1; return; 
         }
         strcpy(fnn,pdbd); strcat(fnn,"/"); strcat(fnn,fn);
         fin=fopen(fnn,"r");
         if (fin==NULL) { 
             fprintf(stderr, "\n*** ERROR IN OPENING %s\n",fnn);
             ier1=1; return; 
         }
      }

      nnt=0; xk[0]=yk[0]=zk[0]=0.; ird=0.;
      
      while(fgets(line, sizeof(line), fin) != NULL) {
          if (strncmp(line, "ATOM", 4) == 0) {
              // Parse using fixed width positions standard for PDB
              // Atom Serial: 6-11
              strncpy(sub_str, line+6, 5); sub_str[5] = 0;
              atnb = atoi(sub_str);
              
              // Atom Name: 12-16
              char tatt[5];
              strncpy(tatt, line+12, 4); tatt[4] = 0;
              
              // Residue Name: 17-20
              char trest[5];
              strncpy(trest, line+17, 3); trest[3] = 0; 
              // Note: strict PDB is 17-20 (3 chars). simplified here.
              
              // Chain ID: 21
              char chain = line[21];
              
              // Res Seq: 22-26
              strncpy(sub_str, line+22, 4); sub_str[4] = 0;
              rner = atoi(sub_str);
              
              // X, Y, Z: 30-38, 38-46, 46-54
              strncpy(sub_str, line+30, 8); sub_str[8] = 0; xl = atof(sub_str);
              strncpy(sub_str, line+38, 8); sub_str[8] = 0; yl = atof(sub_str);
              strncpy(sub_str, line+46, 8); sub_str[8] = 0; zl = atof(sub_str);
              
              // Checks
              if (cresi(trest)==0) continue; 
              if (catom(tatt)==0)  continue;
              if ((repr==2)&&(catom(tatt)!=11)) continue; // CA only

              nnt++;
              if (nnt>(DIMA-2)) { 
                  fprintf(stderr, "\n*** ERROR: too many atoms in %s\n",fn);
                  ier1=1; fclose(fin); return; 
              }

              atnumo[nnt]=atnb; rsnumo[nnt]=rner;
              // Clean names
              atnamo[nnt][0]=tatt[0]; atnamo[nnt][1]=tatt[1];
              atnamo[nnt][2]=tatt[2]; atnamo[nnt][3]=tatt[3];
              rsnamo[nnt][0]=trest[0]; rsnamo[nnt][1]=trest[1];
              rsnamo[nnt][2]=trest[2]; rsnamo[nnt][3]=trest[3];
              
              xk[nnt]=xl; yk[nnt]=yl; zk[nnt]=zl;
              xk[0] +=xl; yk[0] +=yl; zk[0] +=zl;
          }
      }
      
      if (nnt > 0) {
          xk[0]/=(float)nnt; yk[0]/=(float)nnt; zk[0]/=(float)nnt;
          for (i=1; i<=nnt; i++) {
             rx=xk[i]-xk[0]; ry=yk[i]-yk[0]; rz=zk[i]-zk[0];
             ird += (rx*rx + ry*ry + rz*rz);
          }
          ird /= (float)nnt; 
          *irad = (float)sqrt((double)ird); 
          *nt = nnt;
      } else {
          *nt = 0; *irad = 0;
      }
      
      fclose(fin);
 }
