/*
===============================================================================
                                                                        CIPAR.C
                                                               Parameters input
===============================================================================
*/
#include <stdio.h>
#include "cdim.h"
 
#define SKIP  fgets(dum,100,fin)
#define SKIP1 while(1){fscanf(fin,"%c",&ch);if(ch=='=')break;}

extern float csize,vfrac,T,mind,maxd,dst;
extern int   cellt,mov,mtch,mcol,rlig,nstep,pbc,bal,nout,stat,wst,refst,clus;
extern int   ndr,pp,scren,ier1;
extern FILE *fopen(),*fin,*fflog;
 
void cipar(void)
{
char  dum[110],ch;
 
      fin=fopen("cpar.gc","r");
      if (fin==NULL)
        { fprintf(fflog,"\n*** ERROR IN OPENING cpar.gc\n\n"); ier1=1; return;}

      SKIP1;fscanf(fin,"%d",&cellt); SKIP;      // cell type
      if ((cellt>2)||(cellt<1)) 
        { fprintf(fflog,"\nERROR: unknown cell type\n"); ier1=1; return;}
      SKIP1;fscanf(fin,"%f",&csize); SKIP;      // cell size
      SKIP1;fscanf(fin,"%f",&vfrac); SKIP;      // volume fraction
      SKIP1;fscanf(fin,"%f",&T);     SKIP;      // temperature
      SKIP1;fscanf(fin,"%d",&mov);   SKIP;      // move or not           
      SKIP1;fscanf(fin,"%d",&mtch);  SKIP;      // no. of matches           
      SKIP1;fscanf(fin,"%f",&mind);  SKIP;      // min move of molecule           
      SKIP1;fscanf(fin,"%f",&maxd);  SKIP;      // max move of molecule           
      SKIP1;fscanf(fin,"%f",&dst);   SKIP;      // add to sum of radii           
      SKIP1;fscanf(fin,"%d",&mcol);  SKIP;      // collision check mode  
      if ((mcol>2)||(mcol<1)) 
        { fprintf(fflog,"\nERROR: unknown collision check\n"); ier1=1; return;}
      SKIP1;fscanf(fin,"%d",&rlig);  SKIP;      // ligand random mode  
      if ((rlig>2)||(rlig<1)) 
        { fprintf(fflog,"\nERROR: unknown ligand rand mode\n");ier1=1; return;}
      SKIP1;fscanf(fin,"%d",&nstep); SKIP;      // no. of steps in trajectory          
      SKIP1;fscanf(fin,"%d",&pbc);   SKIP;      // periodic boundary condition
      if ((cellt==2)&&pbc)
        { fprintf(fflog,"\nERROR: cell type cannot have pbc\n");ier1=1;return;}      
      SKIP1;fscanf(fin,"%d",&bal);   SKIP;      // balance condition          
      SKIP1;fscanf(fin,"%d",&nout);  SKIP;      // coordinates output          
      SKIP1;fscanf(fin,"%d",&stat);  SKIP;      // statistics  output         
      SKIP1;fscanf(fin,"%d",&wst);   SKIP;      // window for statistics         
      if ((wst>DIMW)||(wst<1)) 
        { fprintf(fflog,"\nERROR: stat window outside range\n");ier1=1;return;} 
      SKIP1;fscanf(fin,"%d",&refst); SKIP;      // msd,clust etc reference step         
      SKIP1;fscanf(fin,"%d",&pp);    SKIP;      // protein num ID to follow         
      SKIP1;fscanf(fin,"%d",&clus);  SKIP;      // cluster data output         
      SKIP1;fscanf(fin,"%d",&ndr);   SKIP;      // rotational diff data output         
      SKIP1;fscanf(fin,"%d",&scren); SKIP;      // screen output         

      fclose(fin);
}
