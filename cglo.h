/*
===============================================================================
                                                                        CGLO.H
                                                               Global variables
===============================================================================

X,Y,Z                   atom coordinates (X,Y,Z[][0] - center of mass)
Xs,Ys,Zs                atom coordinates original
X_ca,Y_ca,Z_ca          CA cached coordinates
Xp0,Yp0,Zp0             CM coordinates corrected for PBC
Xs0,Ys0,Zs0             CM reference coordinates for MSD
U11,...,U33             rotation matrices from original coordinates
dx,dy,dz                translations
d11,...,d33             rotation matrices
Emch                    energy of match
E                       energy of molecule
T                       temperature
cellt                   cell type
csize                   cell size
vfrac                   volume fraction
mind                    min single move of molecule (not to stick to same mtch)
maxd                    max single move of molecule
dst                     add to sum or radii to define neighbor
mcol                    collision check mode
rlig                    ligand random mode
mov                     move or just pack
nat_ca                  sizes for cached X/Y/Z_ca
ncop                    number of molecule copies
mrad                    molecule radius
sh                      CM shift per step
Nac                     number of accepted moves per step
mfra                    molecule fraction
mnum                    numerical ID of molecule for entry to arrays
nat                     number of atoms in molecule
Cn[][]=1,0 (Cm copy)    connectivity matrix (copy for clustering)
atnum,rsnum             atom, residue number
atnam(o),rsnam(o)       atom, residue name
nml                     number of different molecules (names)
nmol                    number of all molecules
mtch                    number of matches in complex
nstep                   number of steps in trajectory
nout                    coordinates output increment
stat_flag                    statistics  output increment
wst                     window for statistics averaging
repr                    all/c-alpha
pbc                     periodic boundary condition switch
bal                     balance condition implcit or explicit
refst                   reference step for MSD
clus                    cluster data output
ndr                     rotational diffusion output
scren                   screen output interval
pp                      protein to follow for MSD, accept (serial # in cmol.gr)
ier1,2                  error codes
mnam                    molecule names

---------------------------------------------------------------------------- */

float X   [DIMM+2][DIMA+1],Y   [DIMM+2][DIMA+1],Z   [DIMM+2][DIMA+1];
float Xs  [DIMN+2][DIMA+1],Ys  [DIMN+2][DIMA+1],Zs  [DIMN+2][DIMA+1];
float X_ca[DIMM+2][DIMA+1],Y_ca[DIMM+2][DIMA+1],Z_ca[DIMM+2][DIMA+1];
float Xp0[DIMM+2],Yp0[DIMM+2],Zp0[DIMM+2];
float Xs0[DIMM+2],Ys0[DIMM+2],Zs0[DIMM+2];
float U11[DIMM+2],U12[DIMM+2],U13[DIMM+2];
float U21[DIMM+2],U22[DIMM+2],U23[DIMM+2];
float U31[DIMM+2],U32[DIMM+2],U33[DIMM+2];
float Emch[DIMN+2][DIMN+1][DIMH+1],E[DIMM+2];
float dx[DIMN+2][DIMN+1][DIMH+1];
float dy[DIMN+2][DIMN+1][DIMH+1];
float dz[DIMN+2][DIMN+1][DIMH+1];
float d11[DIMN+2][DIMN+1][DIMH+1],d12[DIMN+2][DIMN+1][DIMH+1],d13[DIMN+2][DIMN+1][DIMH+1];
float d21[DIMN+2][DIMN+1][DIMH+1],d22[DIMN+2][DIMN+1][DIMH+1],d23[DIMN+2][DIMN+1][DIMH+1];
float d31[DIMN+2][DIMN+1][DIMH+1],d32[DIMN+2][DIMN+1][DIMH+1],d33[DIMN+2][DIMN+1][DIMH+1];
float mrad[DIMN+2],mfra[DIMN+2];
float sh[DIMM+2];
float csize,vfrac,mind,maxd,dst,T;
int   ncop[DIMN+2];
int   nat_ca[DIMN+2];
int   mnum[DIMM+2];
int   nat[DIMN+2];
int   Cn[DIMM+2][DIMM+1],Cm[DIMM+2][DIMM+1];
int   atnum[DIMN+2][DIMA+1],rsnum[DIMN+2][DIMA+1];
int   nml,mov,nmol,mtch,nstep,nout,repr,mcol,rlig,stat_flag,wst,pbc,bal,refst,clus;
int   ndr,pp,cellt,Nac,scren;
int   ier1,ier2;
char  mnam[DIMN+2][5];
char  atnam[DIMN+2][DIMA+1][4],rsnam[DIMN+2][DIMA+1][4];
char  atnamo[DIMA+2][4],rsnamo[DIMA+2][4];
FILE  *fopen(),*fin,*finn,*finm,*fflog,*fout,*fp_stat,*fclst,*fdrot;
