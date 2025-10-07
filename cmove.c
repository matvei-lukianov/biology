/*
===============================================================================
                                                                        CMOVE.C
                                                                 Move molecules
===============================================================================
*/
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <immintrin.h>
#include "cdim.h"

extern float X   [][DIMA+1],Y   [][DIMA+1],Z   [][DIMA+1];
extern float Xs  [][DIMA+1],Ys  [][DIMA+1],Zs  [][DIMA+1];
extern float X_ca[][DIMA+1],Y_ca[][DIMA+1],Z_ca[][DIMA+1];
extern float Xp0[],Yp0[],Zp0[];
extern float Emch[][DIMN+1][DIMH+1],E[];
extern float dx[][DIMN+1][DIMH+1],dy[][DIMN+1][DIMH+1],dz[][DIMN+1][DIMH+1];
extern float d11[][DIMN+1][DIMH+1],d12[][DIMN+1][DIMH+1],d13[][DIMN+1][DIMH+1];
extern float d21[][DIMN+1][DIMH+1],d22[][DIMN+1][DIMH+1],d23[][DIMN+1][DIMH+1];
extern float d31[][DIMN+1][DIMH+1],d32[][DIMN+1][DIMH+1],d33[][DIMN+1][DIMH+1];
extern float U11[],U12[],U13[],U21[],U22[],U23[],U31[],U32[],U33[];
extern float mrad[],sh[];
extern float mind,maxd,dst,T,csize;
extern int   Cn[][DIMM+1],mnum[],nat[],ncop[],nat_ca[];
extern int   cellt,nmol,rlig,mtch,mcol,pbc,bal,Nac,pp;
extern char  atnam[][DIMA+1][4];

// ================================================== COLLISION CHECKS VERSIONS 

// ---------------------------------------------------- Preload xi for CA atoms
void preload_xi_ca(float *xi,float *yi,float *zi,int mnl,int nt,
        float *out_xi_ca,float *out_yi_ca,float *out_zi_ca,size_t *out_ca_size)
{
  (*out_ca_size) = 0;
  for (int i_atom = 1; i_atom <= nt; ++i_atom)
  {
    if ((atnam[mnl][i_atom][1]!='C') || (atnam[mnl][i_atom][2]!='A')) continue;

    out_xi_ca[*out_ca_size] = xi[i_atom];
    out_yi_ca[*out_ca_size] = yi[i_atom];
    out_zi_ca[*out_ca_size] = zi[i_atom];
    (*out_ca_size)++;
  }
}
// ---------------------------------------------------------------- SSE version

void check_collision_sse (int l,int ll,float *xi,float *yi,float *zi,
     int mnl,int mnlc,int nt,int ntlc,float mrdl,float cad,
     float *xi_ca,float *yi_ca,float *zi_ca,size_t ca_size,int *Nln,int *col)
{
  int lc,i,j;
  float x0,y0,z0,dist,smrad,smr;
  
  for (lc=1; lc<=nmol; lc++)
    {
    if ((lc==l)||(lc==ll)) continue;
    mnlc=mnum[lc]; ntlc=nat_ca[mnlc];
    smr=mrdl+mrad[mnlc]; smrad=smr+dst; 
    smr=smr*smr; smrad=smrad*smrad;
    x0=xi[0]-X[lc][0]; y0=yi[0]-Y[lc][0]; z0=zi[0]-Z[lc][0];
    dist=x0*x0+y0*y0+z0*z0;

    if (dist>smrad) continue;                                   // Not neighbor
          
    if (bal) (*Nln)+=1;                         // No. of possible moves for ln
          
    if (mcol==1) {if (dist<smr) {*col=1;break;}            // Collision sum rad
                  else          {*col=0;continue;}}
          
    __m128 v_cad = _mm_set1_ps(cad);
    for (i=0; i<ca_size/4*4; i+=4)               // Collision check CA distance
    {
      __m128 v_x1 = _mm_load_ps(&xi_ca[i]);
      __m128 v_y1 = _mm_load_ps(&yi_ca[i]);
      __m128 v_z1 = _mm_load_ps(&zi_ca[i]);
      for (j=0; j<ntlc; j++)
      {
        __m128 v_x = _mm_set1_ps(X_ca[lc][j]);
        __m128 v_y = _mm_set1_ps(Y_ca[lc][j]);
        __m128 v_z = _mm_set1_ps(Z_ca[lc][j]);
        v_x = _mm_sub_ps(v_x1, v_x);
        v_y = _mm_sub_ps(v_y1, v_y);
        v_z = _mm_sub_ps(v_z1, v_z);

        v_x = _mm_mul_ps(v_x, v_x);
        v_y = _mm_mul_ps(v_y, v_y);
        v_z = _mm_mul_ps(v_z, v_z);

        __m128 v_dist = _mm_add_ps(_mm_add_ps(v_x, v_y), v_z);
        int mask = _mm_movemask_ps(_mm_cmpgt_ps(v_cad, v_dist));
        if (mask) {*col=1; break;}                         // Collision CA dist
      }
      if (*col) break;                                 // Colliding molecule lc
    }
    if (*col) break;                                         // Colliding match

    // fallback for ca_size % 4 last elements
    for (; i<ca_size; i++)                       // Collision check CA distance
      {
        float x1=xi_ca[i]; float y1=yi_ca[i]; float z1=zi_ca[i];
        for (j=0; j<ntlc; j++)
          {
          float x2=x1-X_ca[lc][j];float y2=y1-Y_ca[lc][j];float z2=z1-Z_ca[lc][j];
          float dist=x2*x2+y2*y2+z2*z2;
          if (dist<cad) {*col=1; break;}                   // Collision CA dist
          }
        if (*col) break;                               // Colliding molecule lc
      }
    if (*col) break;                                         // Colliding match
    }
}
// ---------------------------------------------------------------- AVX version

#ifdef __AVX2__      
void check_collision_avx (int l,int ll,float *xi,float *yi,float *zi,
     int mnl,int mnlc,int nt,int ntlc,float mrdl,float cad,
     float *xi_ca,float *yi_ca,float *zi_ca,size_t ca_size,int *Nln,int *col)
{
  int lc,i,j;
  float x0,y0,z0,dist,smrad,smr;
  
  for (lc=1; lc<=nmol; lc++)
    {
    if ((lc==l)||(lc==ll)) continue;
    mnlc=mnum[lc]; ntlc=nat_ca[mnlc];
    smr=mrdl+mrad[mnlc]; smrad=smr+dst; 
    smr=smr*smr; smrad=smrad*smrad;
    x0=xi[0]-X[lc][0]; y0=yi[0]-Y[lc][0]; z0=zi[0]-Z[lc][0];
    dist=x0*x0+y0*y0+z0*z0;

    if (dist>smrad) continue;                                   // Not neighbor
          
    if (bal) (*Nln)+=1;                         // No. of possible moves for ln
          
    if (mcol==1) {if (dist<smr) {*col=1;break;}            // Collision sum rad
                  else          {*col=0;continue;}}
          
    __m256 v_cad = _mm256_set1_ps(cad);
    for (i=0; i<ca_size/8*8; i+=8)               // Collision check CA distance
    {
      __m256 v_x1 = _mm256_load_ps(&xi_ca[i]);
      __m256 v_y1 = _mm256_load_ps(&yi_ca[i]);
      __m256 v_z1 = _mm256_load_ps(&zi_ca[i]);
      for (j=0; j<ntlc; j++)
      {
        __m256 v_x = _mm256_set1_ps(X_ca[lc][j]);
        __m256 v_y = _mm256_set1_ps(Y_ca[lc][j]);
        __m256 v_z = _mm256_set1_ps(Z_ca[lc][j]);
        v_x = _mm256_sub_ps(v_x1, v_x);
        v_y = _mm256_sub_ps(v_y1, v_y);
        v_z = _mm256_sub_ps(v_z1, v_z);

        v_x = _mm256_mul_ps(v_x, v_x);
        v_y = _mm256_mul_ps(v_y, v_y);
        v_z = _mm256_mul_ps(v_z, v_z);

        __m256 v_dist = _mm256_add_ps(_mm256_add_ps(v_x, v_y), v_z);
        int mask = _mm256_movemask_ps(_mm256_cmp_ps(v_cad, v_dist, _CMP_GT_OQ));
        if (mask) {*col=1; break;}                         // Collision CA dist
      }
      if (*col) break;                                 // Colliding molecule lc
    }
    if (*col) break;                                         // Colliding match

    { // fallback for ca_size % 8 last elements.
      __m128 v_cad = _mm_set1_ps(cad);
      for (; i<ca_size/4*4; i+=4)                // Collision check CA distance
      {
        __m128 v_x1 = _mm_load_ps(&xi_ca[i]);
        __m128 v_y1 = _mm_load_ps(&yi_ca[i]);
        __m128 v_z1 = _mm_load_ps(&zi_ca[i]);
        for (j=0; j<ntlc; j++)
        {
          __m128 v_x = _mm_set1_ps(X_ca[lc][j]);
          __m128 v_y = _mm_set1_ps(Y_ca[lc][j]);
          __m128 v_z = _mm_set1_ps(Z_ca[lc][j]);
          v_x = _mm_sub_ps(v_x1, v_x);
          v_y = _mm_sub_ps(v_y1, v_y);
          v_z = _mm_sub_ps(v_z1, v_z);

          v_x = _mm_mul_ps(v_x, v_x);
          v_y = _mm_mul_ps(v_y, v_y);
          v_z = _mm_mul_ps(v_z, v_z);

          __m128 v_dist = _mm_add_ps(_mm_add_ps(v_x, v_y), v_z);
          int mask = _mm_movemask_ps(_mm_cmpgt_ps(v_cad, v_dist));
          if (mask) {*col=1; break;}                       // Collision CA dist
        }
        if (*col) break;                               // Colliding molecule lc
      }
      if (*col) break;                                       // Colliding match
    }

    // fallback for ca_size % 4 last elements.
    for (; i<ca_size; i++)           // Collision check CA distance
      {
        float x1=xi_ca[i]; float y1=yi_ca[i]; float z1=zi_ca[i];
        for (j=0; j<ntlc; j++)
          {
          float x2=x1-X_ca[lc][j];float y2=y1-Y_ca[lc][j];float z2=z1-Z_ca[lc][j];
          float dist=x2*x2+y2*y2+z2*z2;
          if (dist<cad) {*col=1; break;}                   // Collision CA dist
          }
        if (*col) break;                               // Colliding molecule lc
      }
    if (*col) break;                                        // Colliding match
    }
}
#endif
// ======================================================= END COLLISION CHECKS


void cmove(void)
{
float xi[DIMA+2],yi[DIMA+2],zi[DIMA+2];
float xi_ca[DIMA+2] __attribute__ ((aligned (32)));
float yi_ca[DIMA+2] __attribute__ ((aligned (32)));
float zi_ca[DIMA+2] __attribute__ ((aligned (32)));
size_t ca_size = 0;
float x0,y0,z0,x00,y00,z00,xl0,yl0,zl0,xl1,yl1,zl1,xl,yl,zl,xll0,yll0,zll0;
float x1,y1,z1,x2,y2,z2,u,v,w;
float a11,a12,a13,a21,a22,a23,a31,a32,a33,b11,b12,b13,b21,b22,b23,b31,b32,b33;
float smrad,smr,dist,mrdl,rms,maxd2,mind2,cad,el,eln,p,pm,xii,yii,zii,sx,sy,sz;
float R,R2;
int   id[DIMM+2],neib[DIMB+2];
int   l,ll,lc,la,lx,k,i,j,mnl,mnll,mnlc,nt,n2,ntlc,col,nmt,Nl,Nln;
int   pbcx0,pbcy0,pbcz0,pbcxc,pbcyc,pbczc;

      cad=8.; cad=cad*cad;               // CA distance
      mind2=mind*mind; maxd2=maxd*maxd;  // Min & Max RMSD move of molecule
      nmt=0;                             // Number of molecules attempted move
      Nac=0;                             // Number of accepted moves
      R=csize/2.; R2=R*R;                // Radius of sphere
      for (l=1;l<=nmol;l++) { id[l]=0;   // No molecules attempted move yet
                              sh[l]=0.;} // Shift per step

      // Caching CA coordinates
      for (lc=1; lc<=nmol; lc++)
      {
        int n_ca = 0;
        int m = mnum[lc]; 
        int n = nat[m];
        for (j=1; j<=n; j++)
        {
          if ((atnam[m][j][1]!='C') || (atnam[m][j][2]!='A')) continue;

          X_ca[lc][n_ca]=X[lc][j];Y_ca[lc][n_ca]=Y[lc][j];Z_ca[lc][n_ca]=Z[lc][j];
          n_ca++;
        }
        nat_ca[m] = n_ca;
      }
      while(nmt<nmol)
        {
         while(1)
           {l=rand()%nmol+1;
            if((rlig==1)&&(id[l]==0))        // Random order - check if l tried
              {id[l]=1; nmt=0; for(k=1;k<=nmol;k++)nmt+=id[k]; break;}//Ctr sum
            if (rlig==2)                   // Random choice - may repeat & skip
              { nmt++; break; }
           } 
         mnl=mnum[l]; mrdl=mrad[mnl]; nt=nat[mnl]; n2=nt/2; col=Nl=0;
         el=E[l]; xl0=X[l][0]; yl0=Y[l][0]; zl0=Z[l][0];
//-------------------------------------------------- Check neighbors to move to
         for (ll=1; ll<=nmol; ll++)                // Determine neighbors for l
           { 
            if (ll==l) continue;                        // Can't bind to itself
            mnll=mnum[ll];
            smrad=mrdl+mrad[mnll]+dst; smrad=smrad*smrad;
            xll0=X[ll][0];yll0=Y[ll][0];zll0=Z[ll][0];
            x0=xl0-xll0;  y0=yl0-yll0;  z0=zl0-zll0;
            dist=x0*x0+y0*y0+z0*z0;
            if (dist>smrad) continue;                           // Not neighbor
            Nl++;                                // No. of possible moves for l
            neib[Nl]=ll;                                     // Neighbor num ID
           }
         if (Nl==0) continue;                 // l has no neighbors, can't move
         if (Nl> 1) lx=rand()%Nl+1;             // Random R in the neighborhood
         else       lx=1;
         ll=neib[lx];         
         mnll=mnum[ll];
         xll0=X[ll][0];yll0=Y[ll][0];zll0=Z[ll][0];
         b11=U11[ll];b12=U12[ll];b13=U13[ll];
         b21=U21[ll];b22=U22[ll];b23=U23[ll];
         b31=U31[ll];b32=U32[ll];b33=U33[ll];
         col=0;

         k=rand()%mtch+1;                                       // Random match
            
         if (bal==0)                            // Metropolis (implcit balance)
           {
            eln=Emch[mnll][mnl][k]; 
            pm=(float)exp((double)(-(eln-el)/T)); 
            p=(float)(rand()%100+1)/100.; 
            if (p>pm) continue;
           }
//--------------------------------------------------------------- Putative move
// Translate & rotate L in original coordinates; ends up in R original coord
         a11=d11[mnll][mnl][k];a12=d12[mnll][mnl][k];a13=d13[mnll][mnl][k];
         a21=d21[mnll][mnl][k];a22=d22[mnll][mnl][k];a23=d23[mnll][mnl][k];
         a31=d31[mnll][mnl][k];a32=d32[mnll][mnl][k];a33=d33[mnll][mnl][k];
         u=dx[mnll][mnl][k];   v=dy[mnll][mnl][k];   w=dz[mnll][mnl][k];

         x0 =Xs[mnl] [0]; y0 =Ys[mnl] [0]; z0 =Zs[mnl] [0];
         x00=Xs[mnll][0]; y00=Ys[mnll][0]; z00=Zs[mnll][0];
         rms=0.;
               
         for (i=0; i<=nt; i++)
           {                                                // Relative to CM L
            xl=Xs[mnl][i]-x0; yl=Ys[mnl][i]-y0; zl=Zs[mnl][i]-z0;
                 
            xl1=a11*xl+a12*yl+a13*zl+u+x0;
            yl1=a21*xl+a22*yl+a23*zl+v+y0;
            zl1=a31*xl+a32*yl+a33*zl+w+z0;
// Translate & rotate L with R to R current coordinates
            xl1-=x00; yl1-=y00; zl1-=z00;                   // Relative to CM R

            xi[i]=b11*xl1+b12*yl1+b13*zl1+xll0;
            yi[i]=b21*xl1+b22*yl1+b23*zl1+yll0;
            zi[i]=b31*xl1+b32*yl1+b33*zl1+zll0;
           }
//--------------------------------------------------------- Check move validity

         if((cellt==1)&&(pbc==0)&&                                    // Cube
            ((xi[0]<0.)   ||(yi[0]<0.)   ||(zi[0]<0.)  ||       
             (xi[0]>csize)||(yi[0]>csize)||(zi[0]>csize))) continue;
         if((cellt==2)&&                                              // Sphere    
            (((xi[0]-R)*(xi[0]-R)+(yi[0]-R)*(yi[0]-R)+(zi[0]-R)*(zi[0]-R))>R2))
                                                           continue;
                                            // Check move length by 3-atom RMSD
         xii=xi[1] -X[l][1];  yii=yi[1] -Y[l][1];  zii=zi[1] -Z[l][1];
         xii=xii*xii;yii=yii*yii;zii=zii*zii; rms+=(xii+yii+zii);
         xii=xi[n2]-X[l][n2]; yii=yi[n2]-Y[l][n2]; zii=zi[n2]-Z[l][n2];
         xii=xii*xii;yii=yii*yii;zii=zii*zii; rms+=(xii+yii+zii);
         xii=xi[nt]-X[l][nt]; yii=yi[nt]-Y[l][nt]; zii=zi[nt]-Z[l][nt];
         xii=xii*xii;yii=yii*yii;zii=zii*zii; rms+=(xii+yii+zii);rms/=3.;

         if ((rms<mind2)||(rms>maxd2)) continue;         // Move too short/long

         pbcx0=pbcxc=pbcy0=pbcyc=pbcz0=pbczc=0; // Periodic Boundary Cond
         if(xi[0]<0.)   {for(i=0;i<=nt;i++){xi[i]=xi[i]+csize; pbcx0=1;}}
         if(yi[0]<0.)   {for(i=0;i<=nt;i++){yi[i]=yi[i]+csize; pbcy0=1;}}
         if(zi[0]<0.)   {for(i=0;i<=nt;i++){zi[i]=zi[i]+csize; pbcz0=1;}}
         if(xi[0]>csize){for(i=0;i<=nt;i++){xi[i]=xi[i]-csize; pbcxc=1;}}
         if(yi[0]>csize){for(i=0;i<=nt;i++){yi[i]=yi[i]-csize; pbcyc=1;}}
         if(zi[0]>csize){for(i=0;i<=nt;i++){zi[i]=zi[i]-csize; pbczc=1;}}
 
         col=0; Nln=0;                    // Calculate Nln and Check collisions

         preload_xi_ca(&xi[0],&yi[0],&zi[0],mnl,nt,&xi_ca[0],&yi_ca[0],&zi_ca[0],&ca_size);

#ifdef __AVX2__
//#warning AVX implementation is used.
         check_collision_avx(l,ll,&xi[0],&yi[0],&zi[0],mnl,mnlc,nt,ntlc,mrdl,
                          cad,&xi_ca[0],&yi_ca[0],&zi_ca[0],ca_size,&Nln,&col);
#else
//#warning AVX is not supported. SSE implementation is used.
         check_collision_sse(l,ll,&xi[0],&yi[0],&zi[0],mnl,mnlc,nt,ntlc,mrdl,
                          cad,&xi_ca[0],&yi_ca[0],&zi_ca[0],ca_size,&Nln,&col);
#endif
         if (col) continue;                                      // Reject move
             
         if (bal)       // Metropolos expl balance;for non-coll attempted moves
           {
            eln=Emch[mnll][mnl][k]; 
            pm=(float)exp((double)(-(eln-el)/T));
            pm*=((float)Nl/(float)Nln);                   // Balance correction
            p=(float)(rand()%100+1)/100.; 
            if (p>pm) continue;                                  // Reject move
           }
//--------------------------------------------------------------- Accepted move
         if (mnum[l]==pp) Nac++;                            // Accepted counter
         for (la=1; la<=nmol; la++)                                   // Detach
           { if (Cn[l][la]) { Cn[l][la]=Cn[la][l]=0; E[la]=E[la]-E[l]; }}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - Stats
         sx=xi[0]-X[l][0];sy=yi[0]-Y[l][0];sz=zi[0]-Z[l][0]; //Shift/step
         Xp0[l]=Xp0[l]+sx;Yp0[l]=Yp0[l]+sy;Zp0[l]=Zp0[l]+sz; //CM
                                                              // PBC correction
         if (pbc) { if(pbcx0) { Xp0[l]=Xp0[l]-csize; sx-=csize; }
                    if(pbcy0) { Yp0[l]=Yp0[l]-csize; sy-=csize; }
                    if(pbcz0) { Zp0[l]=Zp0[l]-csize; sz-=csize; }
                    if(pbcxc) { Xp0[l]=Xp0[l]+csize; sx+=csize; }
                    if(pbcyc) { Yp0[l]=Yp0[l]+csize; sy+=csize; }
                    if(pbczc) { Zp0[l]=Zp0[l]+csize; sz+=csize; } }
         sh[l]=sx*sx+sy*sy+sz*sz;
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

         E[l]=eln;                                               // New E for l
         Cn[l][ll]=Cn[ll][l]=1; E[ll]=E[ll]+eln;                      // Attach
 
                                                       // New coordinates for l
         for (i=0;i<=nt;i++)
         {
            X[l][i]=xi[i];Y[l][i]=yi[i];Z[l][i]=zi[i];
         }
                                       // Update cached coordinates accordingly
         for (i=0;i<nat_ca[mnl];i++)
         {
            X_ca[l][i]=xi_ca[i]; Y_ca[l][i]=yi_ca[i]; Z_ca[l][i]=zi_ca[i];
         }
                                    // Update U for L (Current R U * prior L U)
         U11[l]=b11*a11+b12*a21+b13*a31;
         U12[l]=b11*a12+b12*a22+b13*a32;
         U13[l]=b11*a13+b12*a23+b13*a33;
         U21[l]=b21*a11+b22*a21+b23*a31;
         U22[l]=b21*a12+b22*a22+b23*a32;
         U23[l]=b21*a13+b22*a23+b23*a33;
         U31[l]=b31*a11+b32*a21+b33*a31;
         U32[l]=b31*a12+b32*a22+b33*a32;
         U33[l]=b31*a13+b32*a23+b33*a33;
        }
}
