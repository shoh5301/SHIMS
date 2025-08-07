#pragma once

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <sys/io.h>
#include "Utility.h"
#include "Calc.h"
#include "MC.h"

void GGmodule_MPI(int**** grid,double jc[],int dir[],int out[],char ttem0[],double melt[],int pbc[],double const sang[],double rscale[],double tbc[][2]);// Grain growth
void SOLmodule_MPI(int**** grid,double jc[],int dir[],int out[],double melt[],int pbc[],double rscale[],double tbc[][2]);
void AMmodule_MPI(int**** grid,double jc[],int dir[],int out[],double melt[],int pbc[],double rscale[],double tbc[][2],double am[]);
double**** MPI_Tdistr(double**** grid,double**** mpigrid,int mpix0,int mpixlen[],int dir[],int xpbc);
int**** MPI_distr(int**** grid,int**** mpigrid,int mpix0,int mpixlen[],int dir[],int xpbc);
void MPI_BCsync(int**** mpigrid,int mpixlen,int dir[],int xpbc,int mpimode);
void MPI_gather(int**** grid, int**** mpigrid,int mpix0,int mpixlen[],int dir[]);
void MPI_TBCsync(double*** mpigrid,int mpixlen,int dir[],int xpbc);
void MPI_Tsync(double*** fdmt,double*** mpigrid,int mpix0,int mpixlen[],int dir[]);
void heat_transfer_MPI(double**** tgrid,double tbc[][2],int mpidir[],int dir[],int pbc[],double mcdx,double const melt[],int fdm_loop,double k_per_cp,double Tmelt,int mpixlen);
int MPI_MeltPoolCheck(double mode,double direction,double mp[],double radius,double xpos,int mpixlen,int mpidir);
int find_min_rank(double rtim);

