#pragma once

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <sys/io.h>
#include "Utility.h"
#include "Calc.h"

#define _USE_TIME_DEFINES

void MCtrial_GG(int x,int y,int z,int**** grid,int dir[],int pbc[],double jc[],double tbc[][2],double pmob,double realE);
void MCtrial_SOL(int i,int j,int k,int**** grid,double**** fdmt,int dir[],int pbc[],int mode[],double melt[],double Tmax[],double dHfus[],double jc[],double tbc[][2],double Vsite,double Asite,double factor[],double Nsite[],double pcps[],double Qgg);
double gbenergy(int px,int py,int pz,int now[],int**** const grid,double jc[],int dir[],int pbc[],double tbc[][2]); // GBE calculator
double fixed_dt_set(int mcs,int* P_mode,int *fdm_loop,double dtsave[],double dtdMCS[],double dHfus[],double Tinfo[],double factor[],double Vsite,double Asite,int out[],double Tmax[],double pcps[],int**** grid,double**** fdmt,int dir[],int pbc[],double rscale[],double melt[]);
double active_dt_set(int mcs,int* P_mode,int *fdm_loop,double dtsave[],double dtdMCS[],double dHfus[],double Tinfo[],double factor[],double Vsite,double Asite,int out[],double Tmax[],double pcps[],int**** grid,double**** fdmt,int dir[],int pbc[],double rscale[],double melt[]);
double P_grain_growth(double Q,double revTm,double revT,double dx2,double dG);
double P_nucl(double revT,double dGfus,double gam,double factor,double Nc,double dG,double pcps[],double ps);
double P_growth(double revT,double pcps[],double ps,double dG);
double solidi_df(double dHfus,double Tnow,double Tliq,double Tfreez,double T0);
void heat_transfer(double**** tgrid,double tbc[][2],int dir[],int pbc[],double mcdx,double const melt[],int fdm_loop,double k_per_cp,double Tmelt,int mode);//,double dHfus){
void MeltPool_movement(double pos[],double direction,double displ);
void MeltPool(int**** grid,double**** tgrid,double am[],double mp[],int dir[],int pbc[],double Tm,double dH,double pos[]);
int teardrop(double b, double t,int geo); // Melt pool shape curve
void check_melting(int**** grid,double**** tgrid,double Tm,int dir[]);
void Particle_Form(int**** grid,int dir[],int pbc[],double fraction);
void LatHeatDiffus(double* T,double *L,double dT);
double check_growth_dir(int**** grid,int l,int m,int n,int vector,int dir[],int pbc[],double factor);
double N_Tip(double IFE,double dGv,double area);


