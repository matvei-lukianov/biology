/*
===============================================================================
                                                                        CIPAR.C
                                                               Parameters input
===============================================================================
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cdim.h"

extern float csize,vfrac,T,mind,maxd,dst;
extern int   cellt,mov,mtch,mcol,rlig,nstep,pbc,bal,nout,stat_flag,wst,refst,clus;
extern int   ndr,pp,scren,ier1;
extern FILE *fopen(),*fin,*fflog;

// Helper to parse value after '='
void parse_int(FILE *f, int *val, const char *name) {
    char line[256];
    char *ptr;
    if (fgets(line, sizeof(line), f) == NULL) {
        ier1=1; return;
    }
    ptr = strchr(line, '=');
    if (ptr) {
        sscanf(ptr+1, "%d", val);
    } else {
        ier1=1; 
    }
}

void parse_float(FILE *f, float *val, const char *name) {
    char line[256];
    char *ptr;
    if (fgets(line, sizeof(line), f) == NULL) {
        ier1=1; return;
    }
    ptr = strchr(line, '=');
    if (ptr) {
        sscanf(ptr+1, "%f", val);
    } else {
        ier1=1; 
    }
}

void cipar(void)
{
      fin=fopen("cpar.gc","r");
      if (fin==NULL) { 
          fprintf(stderr, "\n*** ERROR IN OPENING cpar.gc\n\n"); 
          ier1=1; return;
      }
      
      parse_int(fin, &cellt, "cellt"); if(ier1) return;
      if ((cellt>2)||(cellt<1)) { fprintf(stderr, "\nERROR: unknown cell type\n"); ier1=1; return;}

      parse_float(fin, &csize, "csize"); if(ier1) return;
      parse_float(fin, &vfrac, "vfrac"); if(ier1) return;
      parse_float(fin, &T, "T"); if(ier1) return;
      parse_int(fin, &mov, "mov"); if(ier1) return;
      parse_int(fin, &mtch, "mtch"); if(ier1) return;
      parse_float(fin, &mind, "mind"); if(ier1) return;
      parse_float(fin, &maxd, "maxd"); if(ier1) return;
      parse_float(fin, &dst, "dst"); if(ier1) return;
      parse_int(fin, &mcol, "mcol"); if(ier1) return;
      
      if ((mcol>2)||(mcol<1)) { fprintf(stderr, "\nERROR: unknown collision check\n"); ier1=1; return;}

      parse_int(fin, &rlig, "rlig"); if(ier1) return;
      if ((rlig>2)||(rlig<1)) { fprintf(stderr, "\nERROR: unknown ligand rand mode\n");ier1=1; return;}

      parse_int(fin, &nstep, "nstep"); if(ier1) return;
      parse_int(fin, &pbc, "pbc"); if(ier1) return;
      
      if ((cellt==2)&&pbc) { fprintf(stderr, "\nERROR: cell type cannot have pbc\n");ier1=1;return;}      
      
      parse_int(fin, &bal, "bal"); if(ier1) return;
      parse_int(fin, &nout, "nout"); if(ier1) return;
      parse_int(fin, &stat_flag, "stat_flag"); if(ier1) return;
      parse_int(fin, &wst, "wst"); if(ier1) return;
      
      if ((wst>DIMW)||(wst<1)) { fprintf(stderr, "\nERROR: stat_flag window outside range\n");ier1=1;return;} 
      
      parse_int(fin, &refst, "refst"); if(ier1) return;
      parse_int(fin, &pp, "pp"); if(ier1) return;
      parse_int(fin, &clus, "clus"); if(ier1) return;
      parse_int(fin, &ndr, "ndr"); if(ier1) return;
      parse_int(fin, &scren, "scren"); if(ier1) return;

      fclose(fin);
}
