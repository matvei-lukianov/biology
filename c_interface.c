/*
===============================================================================
                                                                 C_INTERFACE.C
                                                Bridge for Python ctypes driver
===============================================================================
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cdim.h"
#include "cglo.h"

// Externs from existing files
extern void cipar(void);
extern void cimol(void);
extern void cipdb(char *,float *,float *,float *,int *,float *,int *,int *);
extern void ccops(void);
extern void cimch(char *,int,int);
extern void cpack(void);
extern void cpacc(void);
extern void cmove(void);
extern void ccout(int);
extern void cstat(int,int,int);
extern void cstat(int,int,int);
extern void get_latest_stats(float*);
extern float mrad[]; // Molecule radius from cimol/cipdb

// Globals from cglo.h (already avail via linking, but we might need helpers)
// X, Y, Z, E, etc. are global.

// ---------------------------------------------------------------- Initialization
int init_simulation_interface(void) {
    // This replicates the startup sequence from coper.c
    // Return 0 on success, >0 on error

    float xk[DIMA+2],yk[DIMA+2],zk[DIMA+2];
    float irad;
    int   atnumo[DIMA+2],rsnumo[DIMA+2];
    int   i,ii,j,jj,ic,im,nt,k;
    char  fn[10],fmch[15];

    srand(time(0));
    
    fprintf(stderr, "TRACE: Start init\n");
    // Initialize output pointers
    fflog = stderr;
    
    fflog=fopen("gcell.log","w"); 
    if(!fflog) fflog=stderr;
    
    FILE* fstat_chk = fopen("gcell.st", "w"); if(fstat_chk) fp_stat=fstat_chk;
    FILE* fclst_chk = fopen("clust.st", "w"); if(fclst_chk) fclst=fclst_chk;
    FILE* fdrot_chk = fopen("d_rot.st", "w"); if(fdrot_chk) fdrot=fdrot_chk;
    
    fprintf(stderr, "TRACE: calling cipar\n");
    cipar();
    if(ier1) { fprintf(stderr, "Error in cipar\n"); return 1; }
    fprintf(stderr, "TRACE: cipar done\n");

    // 2. Molecules
    printf("TRACE: calling cimol\n"); fflush(stdout);
    cimol(); 
    if(ier1) { fprintf(stderr, "Error in cimol\n"); return 1; }
    printf("TRACE: cimol done, nml=%d\n", nml); fflush(stdout);

    for (i=1; i<=nml; i++) {
         for(ii=0;ii<4;ii++) fn[ii]=mnam[i][ii];
         fn[4]='.';fn[5]='p';fn[6]='d';fn[7]='b';fn[8]='\0';
         
         printf("TRACE: loading pdb %s\n", fn); fflush(stdout);
         cipdb(fn,xk,yk,zk,&nt,&irad,atnumo,rsnumo); 
         if(ier1) { fprintf(stderr, "Error in cipdb for %s\n", fn); return 1; }
         
         nat[i]=nt; mrad[i]=irad;
         for (k=0; k<=nt; k++) {
            atnum[i][k]=atnumo[k]; rsnum[i][k]=rsnumo[k];
            atnam[i][k][0]=atnamo[k][0];atnam[i][k][1]=atnamo[k][1];
            atnam[i][k][2]=atnamo[k][2];atnam[i][k][3]=atnamo[k][3];
            rsnam[i][k][0]=rsnamo[k][0];rsnam[i][k][1]=rsnamo[k][1];
            rsnam[i][k][2]=rsnamo[k][2];rsnam[i][k][3]=rsnamo[k][3];
            Xs[i][k]=xk[k];Ys[i][k]=yk[k];Zs[i][k]=zk[k];
         }
    }
    printf("TRACE: pdbs loaded\n"); fflush(stdout);
    
    printf("TRACE: calling ccops\n"); fflush(stdout);
    ccops(); 
    if(ier2) { fprintf(stderr, "Error in ccops\n"); return 2; }
    printf("TRACE: ccops done. nmol=%d (global)\n", nmol); fflush(stdout);

    // Assign coordinates to copies
    im=1;                                     
    for (i=1; i<=nml; i++) {        
         nt=nat[i];
         printf("TRACE: mol %d nat=%d ncop=%d\n", i, nt, ncop[i]); fflush(stdout);
         for (ic=0; ic<ncop[i]; ic++) {
            nmol=im+ic; mnum[nmol]=i;   
            for (k=0; k<=nt; k++) { 
                X[nmol][k]=Xs[i][k];Y[nmol][k]=Ys[i][k];Z[nmol][k]=Zs[i][k]; 
            }
         }
         im+=ncop[i];
    }
    
    // 3. Packing
    printf("TRACE: calling cpack/cpacc cellt=%d\n", cellt); fflush(stdout);
    if (cellt==1) cpack();
    if (cellt==2) cpacc();
    printf("TRACE: packing done\n"); fflush(stdout);
    
    // 4. Load Matches (Energy Terms)
    printf("TRACE: loading matches\n"); fflush(stdout);
    for (i=1; i<=nml; i++) { 
        for(ii=0;ii<4;ii++) fmch[ii]=mnam[i][ii]; fmch[4]='-';
        for (j=1; j<=nml; j++) { 
            for(jj=5;jj<9;jj++) fmch[jj]=mnam[j][jj-5]; fmch[9]='.';
            fmch[10]='r';fmch[11]='e';fmch[12]='s';fmch[13]='\0';
            printf("TRACE: cimch %s\n", fmch); fflush(stdout);
            cimch(fmch,i,j); 
            if(ier1||ier2) { fprintf(stderr, "Error in cimch %s\n", fmch); return 3; } 
        }
    }
    printf("TRACE: matches loaded\n"); fflush(stdout);
    
    // Initial PBC correction for Center of Mass
    for (int l=1;l<=nmol;l++) { Xp0[l]=X[l][0];Yp0[l]=Y[l][0];Zp0[l]=Z[l][0]; }

    return 0; // Success
}

// Stat globals to track state between calls
static int nst=0, nco=0, nso=0, nsw=0, ns1=0;

// ---------------------------------------------------------------- Simulation Steps
void run_c_steps(int n_steps) {
    int s;
    for(s=0; s<n_steps; s++) {
        cmove();
        
        // Stats/Output Logic (Replicated from coper.c)
        nst++;
        ns1++; 
        if (ns1==scren) {
             // Optional: Print to stdout/stderr if desired, matching coper.c
             // printf("\n----------\nSTEP %2d ",nst); 
             ns1=0;
        }
        
        nco++; 
        if (nco==nout) { 
            ccout(nst); 
            nco=0; 
        } 
        
        nso++; 
        nsw++; 
        cstat(nst,nso,nsw); 
        
        if(nso==stat_flag) nso=0;
        if(nsw==wst) nsw=0;
    }
}

// ---------------------------------------------------------------- Energy Evaluation
// Logic: For a global jump, we have new coordinates X,Y,Z.
// We must compute Total Energy = Sum over pairs (Energy of best match).
// If no match is close enough, Energy contribution is 0 (or collision penalty).

// Helper: Distance squared
float dist2(float x1, float y1, float z1, float x2, float y2, float z2) {
    return (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
}

// ---------------------------------------------------------------- Energy Evaluation
// Revised: Checks Match Geometry (Dist + Rotation)
float calculate_total_energy(void) {
    float total_E = 0.0;
    int l1, l2, m1, m2, k;
    float x1, y1, z1, x2, y2, z2;
    float dx_rel, dy_rel, dz_rel, d2;
    float tx, ty, tz, td2;
    float u11, u12, u13, u21, u22, u23, u31, u32, u33; // U[l1]
    float v11, v12, v13, v21, v22, v23, v31, v32, v33; // U[l2]
    float w11, w12, w13, w21, w22, w23, w31, w32, w33; // D[k]
    float p11, p12, p13, p21, p22, p23, p31, p32, p33; // Predicted U[l2]
    float diff, err;
    int pair_count = 0;
    
    // Loop over unique pairs
    for (l1=1; l1<=nmol; l1++) {
        // Load Rotation l1
        u11=U11[l1]; u12=U12[l1]; u13=U13[l1];
        u21=U21[l1]; u22=U22[l1]; u23=U23[l1];
        u31=U31[l1]; u32=U32[l1]; u33=U33[l1];
        
        for (l2=l1+1; l2<=nmol; l2++) {
            
            m1 = mnum[l1];
            m2 = mnum[l2];
            
            // Current relative position
            x1 = X[l1][0]; y1 = Y[l1][0]; z1 = Z[l1][0];
            x2 = X[l2][0]; y2 = Y[l2][0]; z2 = Z[l2][0];
            
            dx_rel = x2 - x1; 
            dy_rel = y2 - y1; 
            dz_rel = z2 - z1;

            // PBC Minimal Image
            if (pbc) {
                if (dx_rel >  csize/2.0f) dx_rel -= csize;
                if (dx_rel < -csize/2.0f) dx_rel += csize;
                if (dy_rel >  csize/2.0f) dy_rel -= csize;
                if (dy_rel < -csize/2.0f) dy_rel += csize;
                if (dz_rel >  csize/2.0f) dz_rel -= csize;
                if (dz_rel < -csize/2.0f) dz_rel += csize;
            }
            d2 = dx_rel*dx_rel + dy_rel*dy_rel + dz_rel*dz_rel;
            
            // Spatial Cutoff Optimization
            // If molecules are farther than sum of radii + buffer (50.0f), simple matches are impossible.
            float cutoff = mrad[mnum[l1]] + mrad[mnum[l2]] + 50.0f; 
            // Using 50.0f buffer to match cpar.gc dst parameter
            if (d2 > cutoff*cutoff) continue;
            
            // Load Rotation l2
            v11=U11[l2]; v12=U12[l2]; v13=U13[l2];
            v21=U21[l2]; v22=U22[l2]; v23=U23[l2];
            v31=U31[l2]; v32=U32[l2]; v33=U33[l2];
            
            // Search Matches
            float best_e = 0.0;
            int found_match = 0;
            
            // Limit k? We don't know N matches per pair. Use DIMH? Or check if E is set.
            for(k=1; k<=DIMH; k++) {
                // Heuristic: Empty slot usually has 0 entries
                // check if dx=0, dy=0, dz=0. 
                if(Emch[m1][m2][k]==0.) break; // End of matches list (optimization)
                
                tx = dx[m1][m2][k];
                ty = dy[m1][m2][k];
                tz = dz[m1][m2][k]; 
                
                // 1. Distance Check (Fast)
                // Need to rotate target vector? 
                // In cmove, 'u,v,w' (dx,dy,dz) are ADDED to rotated L.
                // So they are in the frame of L?
                // Let's assume dx,dy,dz is target relative vector in GLOBAL frame IF L was at identity?
                // No, matches are relative. 
                // Vector T_target = U1 * T_rel_match.
                // So actual target diff should be U1 * (tx, ty, tz).
                
                // Rotate match vector by U1
                float rtx = u11*tx + u12*ty + u13*tz;
                float rty = u21*tx + u22*ty + u23*tz;
                float rtz = u31*tx + u32*ty + u33*tz;
                
                // Distance squared to rotated target
                float dxd = (dx_rel - rtx);
                float dyd = (dy_rel - rty);
                float dzd = (dz_rel - rtz);
                float dist_err = dxd*dxd + dyd*dyd + dzd*dzd;
                
                if (dist_err > 25.0f) continue; // 5 Angstrom tolerance (coarse)
                
                // 2. Rotation Check
                // D[k] maps L frame to R frame? 
                // Or are they absolute?
                // cmove: U[l] is updated by multiplying.
                // U_new = D[k] * U_old? (Lines 404-412 in cmove.c)
                // U11[l] = b11*a11 ...
                // This implies U[l] becomes some combo.
                
                // Let's verify relation: U[l2] should relate to U[l1] via D[k].
                // U[l2] ~ D[k] * U[l1].
                // Let's compute P = D[k] * U[l1] and compare to U[l2].
                
                w11=d11[m1][m2][k]; w12=d12[m1][m2][k]; w13=d13[m1][m2][k];
                w21=d21[m1][m2][k]; w22=d22[m1][m2][k]; w23=d23[m1][m2][k];
                w31=d31[m1][m2][k]; w32=d32[m1][m2][k]; w33=d33[m1][m2][k];
                
                p11 = w11*u11 + w12*u21 + w13*u31;
                p12 = w11*u12 + w12*u22 + w13*u32;
                p13 = w11*u13 + w12*u23 + w13*u33;
                
                p21 = w21*u11 + w22*u21 + w23*u31;
                p22 = w21*u12 + w22*u22 + w23*u32;
                p23 = w21*u13 + w22*u23 + w23*u33;
                
                p31 = w31*u11 + w32*u21 + w33*u31;
                p32 = w31*u12 + w32*u22 + w33*u32;
                p33 = w31*u13 + w32*u23 + w33*u33;
                
                // Element-wise difference sum square
                err = (p11-v11)*(p11-v11) + (p12-v12)*(p12-v12) + (p13-v13)*(p13-v13) +
                      (p21-v21)*(p21-v21) + (p22-v22)*(p22-v22) + (p23-v23)*(p23-v23) +
                      (p31-v31)*(p31-v31) + (p32-v32)*(p32-v32) + (p33-v33)*(p33-v33);
                      
                if (err < 1.0f) { // Tolerance on rotation matrix
                    // Match found!
                    // Emch stores NEGATIVE energy (binding E).
                    float e = Emch[m1][m2][k];
                    if (e < best_e) {
                         best_e = e;
                         found_match = 1;
                    }
                }
            }
            if (found_match) total_E += best_e;
            
            pair_count++;
            if (pair_count % 100000 == 0) {
                 fprintf(stderr, "TRACE: Energy calc pair %d\n", pair_count);
            }
        }
    }
    return total_E;
}



// Expose pointers for Python
float* get_ptr_X(void) { return &X[0][0]; } // Flattened access might be tricky if not contiguous 2D
// X is defined as X[DIMM+2][DIMA+1]. This IS contiguous in C row-major.
// So X[0][0] is start. X[1][0] is start of molecule 1.

int* get_ptr_nat(void) { return &nat[0]; }
int* get_ptr_mnum(void) { return &mnum[0]; }
float* get_ptr_E(void) { return &E[0]; }

int get_nmol(void) { return nmol; }



// ---------------------------------------------------------------- Stats Access
void get_stats_interface(float* out_array) {
    get_latest_stats(out_array);
}

// ---------------------------------------------------------------- Rotation Access
// Get pointers to U rotation matrices
float* get_ptr_U11(void) { return &U11[0]; }
float* get_ptr_U12(void) { return &U12[0]; }
float* get_ptr_U13(void) { return &U13[0]; }
float* get_ptr_U21(void) { return &U21[0]; }
float* get_ptr_U22(void) { return &U22[0]; }
float* get_ptr_U23(void) { return &U23[0]; }
float* get_ptr_U31(void) { return &U31[0]; }
float* get_ptr_U32(void) { return &U32[0]; }
float* get_ptr_U33(void) { return &U33[0]; }

// Set rotation from Euler angles (ZYX convention)
// angles: array of [alpha, beta, gamma] for each molecule (1-indexed in C)
// n: number of molecules
void set_rotations_euler(float* alpha, float* beta, float* gamma, int n) {
    for (int l = 1; l <= n; l++) {
        float a = alpha[l-1];
        float b = beta[l-1];
        float g = gamma[l-1];
        
        float sa = sinf(a), ca = cosf(a);
        float sb = sinf(b), cb = cosf(b);
        float sg = sinf(g), cg = cosf(g);
        
        // ZYX Euler angles to rotation matrix
        U11[l] = ca*cb;
        U12[l] = ca*sb*sg - sa*cg;
        U13[l] = ca*sb*cg + sa*sg;
        U21[l] = sa*cb;
        U22[l] = sa*sb*sg + ca*cg;
        U23[l] = sa*sb*cg - ca*sg;
        U31[l] = -sb;
        U32[l] = cb*sg;
        U33[l] = cb*cg;
    }
}

// Get Euler angles from current rotation matrices
void get_rotations_euler(float* alpha, float* beta, float* gamma, int n) {
    for (int l = 1; l <= n; l++) {
        // Extract Euler angles from rotation matrix (ZYX convention)
        float r31 = U31[l];
        float r32 = U32[l];
        float r33 = U33[l];
        float r21 = U21[l];
        float r11 = U11[l];
        
        float b = asinf(-r31); // beta
        float cb = cosf(b);
        
        float a, g;
        if (fabsf(cb) > 1e-6) {
            a = atan2f(r21, r11);      // alpha
            g = atan2f(r32, r33);      // gamma
        } else {
            // Gimbal lock - set alpha=0 and compute gamma
            a = 0.0f;
            g = atan2f(-U12[l], U22[l]);
        }
        
        alpha[l-1] = a;
        beta[l-1] = b;
        gamma[l-1] = g;
    }
}
