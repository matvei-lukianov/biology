/*
===============================================================================
                                                                        CCLUS.C
                     Hierarchical tree clustering of adjecency (contact) matrix
===============================================================================
*/
#include <stdio.h>
#include "cdim.h"

extern int   Cn[][DIMM+1],Cm[][DIMM+1];
extern int   nmol,clus;
extern FILE  *fclst;

void cclus(int nst)
{
float  Nca;
int    Nod[DIMM+2],Nc[100];
int    i,j,k,noc,nt;
int    i1,i2,i3,i4,i5,i6,i7,i8;

      noc=20;            // Number of different cluster occupancies 1,2,... noc
      Nca=0.;                                      // Average cluster occupancy
      nt=0;                                               // Number of clusters
      for(k=1;k<=nmol;k++) Nod[k]=1;  // Initial clusters (nodes) all single
      for(k=1;k<=noc; k++) Nc[k] =0;  // Initial clusters stats
      for(i=1;i<=nmol;i++) {for(j=1;j<=nmol;j++) Cm[i][j]=Cn[i][j];} // Copy Cn

      if (nst==clus)
        {
//for (i=1;i<=nmol;i++)                                  // Output Contact matrix
//{for (j=1;j<=nmol;j++) {if(Cm[i][j]) fprintf(fclst," 1"); else fprintf(fclst," .");}
//fprintf(fclst,"\n");}         
         fprintf(fclst,  "--------------------------------------------------");
         fprintf(fclst,"\n Step                Assembly size");
         fprintf(fclst,"\n   N");
         for (k=1;k<=noc;k++)fprintf(fclst,"%5d",k);fprintf(fclst,"  Average");
         fprintf(fclst,"\n--------------------------------------------------"); 
        }
//--------------------------------------------------------------- Building tree

      for(i1=1;i1<=nmol;i1++)
        {
      for(i2=1;i2<=nmol;i2++) {if(Cm[i1][i2]) {Nod[i1]++; Nod[i2]=0; Cm[i1][i2]=Cm[i2][i1]=0;
      for(i3=1;i3<=nmol;i3++) {if(Cm[i2][i3]) {Nod[i1]++; Nod[i3]=0; Cm[i2][i3]=Cm[i3][i2]=0;
      for(i4=1;i4<=nmol;i4++) {if(Cm[i3][i4]) {Nod[i1]++; Nod[i4]=0; Cm[i3][i4]=Cm[i4][i3]=0;
      for(i5=1;i5<=nmol;i5++) {if(Cm[i4][i5]) {Nod[i1]++; Nod[i5]=0; Cm[i4][i5]=Cm[i5][i4]=0;
      for(i6=1;i6<=nmol;i6++) {if(Cm[i5][i6]) {Nod[i1]++; Nod[i6]=0; Cm[i5][i6]=Cm[i6][i5]=0;
      for(i7=1;i7<=nmol;i7++) {if(Cm[i6][i7]) {Nod[i1]++; Nod[i7]=0; Cm[i6][i7]=Cm[i7][i6]=0;
      for(i8=1;i8<=nmol;i8++) {if(Cm[i7][i8]) {Nod[i1]++; Nod[i8]=0; Cm[i7][i8]=Cm[i8][i7]=0;
                                                     }}}}}}}}}}}}}}} // Deep 8

      for (i=1; i<=nmol; i++)
        { for (k=1; k<=noc; k++) { if (Nod[i]==k) Nc[k]++; }}

//-------------------------------------------------------- Output cluster stats      
      fprintf(fclst,"\n%5d",nst);
      for (k=1; k<=noc; k++) 
        {
         fprintf(fclst,"%5d",Nc[k]);
         Nca+=(float)(Nc[k]*k); nt+=Nc[k];
        }
      Nca/=(float)nt;
      fprintf(fclst,"%7.1f",Nca);
}
