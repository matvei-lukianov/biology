/*
===============================================================================
G R A M M Cell 2.0                                                      CMAIN.C
Global Range Cell Modeling                                         Main program
===============================================================================
*/
#include <stdio.h>
#include "cdim.h"
#include "cglo.h"

void coper(void);

int main(void)
{

//------------------------------------------- Initialize log, stats, parameters
      fflog=fopen("gcell.log","w"); 
      fstat=fopen("gcell.st", "w"); 
      fclst=fopen("clust.st", "w");
      fdrot=fopen("d_rot.st", "w");
//-----------------------------------------------------------------------------

      coper(); if(ier1||ier2)return 0;
      
      fclose(fflog); fclose(fstat); fclose(fclst); fclose(fdrot);
}
