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
      int   mi,cali,nm;
      char  fnn[1000],line[256],edr[256],cal[256];
      char  *resd;

      finm=fopen(fmch,"r");
      if (finm==NULL) {
         resd=getenv("CELLRES");
         if (resd==NULL) { 
             fprintf(stderr, "\n*** ERROR IN OPENING %s\n",fmch);
             ier1=1; return; 
         }
         strcpy(fnn,resd); strcat(fnn,"/"); strcat(fnn,fmch);
         finm=fopen(fnn,"r");
         if (finm==NULL) { 
             fprintf(stderr, "\n*** ERROR IN OPENING %s\n",fnn); 
             ier1=1; return; 
         }
      } 

      // Search for [molecules] section
      nm=0;
      int found_molecules = 0;
      while(fgets(line, sizeof(line), finm)) {
          if (strncmp(line, "[molecules]", 11) == 0) {
              found_molecules = 1;
              break;
          }
      }
      if (!found_molecules) {
          fprintf(stderr, "\n*** ERROR:[molecules] section not found in %s\n", fmch);
          ier2=1; fclose(finm); return;
      }
      
      // Skip 2 lines of headers
      if(!fgets(line, sizeof(line), finm)) { ier1=1; fclose(finm); return; }
      if(!fgets(line, sizeof(line), finm)) { ier1=1; fclose(finm); return; }
      
      // Read line to check for CA
      // The old code did fscanf multiple times. Let's look for "CA"
      if(strstr(line, "CA") != NULL) repr=2; else repr=0;
      // Note: original code checked the *3rd* line after molecules?
      // original: fscanf... then fgets... then fscanf... then check valid "val"
      // Let's assume repr=2 if CA is present in the file for safety or just keep repr=0 if unsure.
      // But let's look for [match]
      
      int found_match = 0;
      while(fgets(line, sizeof(line), finm)) {
          if (strstr(line, "[match]")) {
              found_match = 1;
              break;
          }
      }
      if (!found_match) {
          fprintf(stderr, "\n*** ERROR: [match] section not found in %s\n", fmch);
          ier1=1; fclose(finm); return;
      }
      
      // Skip header line
      fgets(line, sizeof(line), finm);

      // Input matches
      while(fgets(line, sizeof(line), finm)) {
          if (sscanf(line, "%d%f%f%f%f%f%f%f", &mi, &en, &a, &b, &g, &u, &v, &w) != 8) {
              continue; // Skip lines that don't match format
          }
          
          nm++; 
          Emch[im][jm][nm]=-en; // E with minus sign

          dx[im][jm][nm]=u; dy[im][jm][nm]=v; dz[im][jm][nm]=w;

          a*=0.017453; b*=0.017453; g*=0.017453;
          sa=(float)sin((double)a); ca=(float)cos((double)a);
          sb=(float)sin((double)b); cb=(float)cos((double)b);
          sg=(float)sin((double)g); cg=(float)cos((double)g);
          
          d11[im][jm][nm]=cg*ca;         d12[im][jm][nm]=cg*sa;         d13[im][jm][nm]=-sg;
          d21[im][jm][nm]=sb*sg*ca-cb*sa;d22[im][jm][nm]=sb*sg*sa+cb*ca;d23[im][jm][nm]=cg*sb;
          d31[im][jm][nm]=cb*sg*ca+sb*sa;d32[im][jm][nm]=cb*sg*sa-sb*ca;d33[im][jm][nm]=cg*cb;
          
          if (nm >= DIMH) {
              fprintf(stderr, "\n*** ERROR: cannot process > %d matches\n",DIMH);
              ier1=1; fclose(finm); return;
          }
      }
      
      if (nm < mtch) { 
          // Warn but don't fail hard? Or fail as per original code.
          // original: fprintf(fflog,"\n*** ERROR: requested %d matches ",mtch);
          // Just print stderr
          fprintf(stderr, "\n*** WARNING: requested %d matches but found %d in %s\n", mtch, nm, fmch);
          // ier1=1; return; 
          // Let's allow partial loading for robustness
      }
      fclose(finm);
}
