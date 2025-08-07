#pragma once

#ifndef MPI_FUNCTIONS_H
#define MPI_FUNCTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "Calc.h"

extern int n_size;
extern int n_rank;

#define INCAR "input.in"
#define TNOSOLID -1
#define TNOLIQUID 1E6

#define REVRGAS 0.12028
#define RGAS 8.314
#define kB 1.38065E-23
#define REVkB 7.243E22
#define NAV 6.022E23
#define REVNAV 1.6606E-24
//#define HNANA 2.403E14

#define _USE_MATH_DEFINES
#define DTOR M_PI/180
#define RTOD 180/M_PI
#define NNEI 26.0
#define	MAXDG1 900
#define	MAXDG2 1
#define	MAXDG3 1
#define MAXMIS 450.0
#define NINFO 3

#define LIQORI -1
#define POWORI (((MAXDG1+1)*MAXDG2+1)*MAXDG3+10)
#define PTCLORI (((MAXDG1+1)*MAXDG2+1)*MAXDG3+5)
#define MPISURF -99999

#define FHOMO 16.75516
//#define FHOMO 15.08 13.404
//#define SUPERHEATING 1.231
#define SUPERHEATING 1.99
#define FDMSAFE 0.1
#define PROBSAFE 1.0
#define PROBMIN 1.0E-15
#define KAUZMANN 0.48  // Kauzmann T ~ 0.48 Tm // set 0 to make maximum supercooling Tm
#define GAUSS3D 0.06349363593
#define GAUSS2D 0.1591549431

#define MLAT 0
#define MFIT 1
#define MTEM 2
#define MPRO 3

int prob_check(double prob);
void voxcpy(int paste[],int copy[]);
//void nei6sel(int i,int *x,int *y, int *z,int dir[],int pbc[]);
void nei26sel(int i,int *x,int *y,int *z,int dir[],int pbc[]);// NN random selecter
void nei26sel_nodeter(int i,int *x,int *y,int *z,int dir[],int pbc[]);
int get_nei26_index(int dx,int dy,int dz);
int VN_to_Moore(int VNori);
int deter(int a,int lim,int pbc);// PBC determiner
double misorient(int a,int b);
double GBE(int delta,double cutoff,double GBEmax);
int nucl_cond(int**** grid,int x,int y,int z,int dir[],int pbc[],int situation);
int mpi_surf_check(int i,int x,int y,int z,int**** grid);
int surf_check(int i,int x,int y,int z,int dir[],int pbc[]); // return >=0 if surface
void pow2grain(int ori[]);
int latent_heat_check(double*** tgrid,int dir[]);
double average_temp(double*** temp,int x,int y,int z,int inter,int dir[],int pbc[]);
double recover_temp(double*** temp,int x,int y,int z,int del,int pts[]);
double distance(int p11,int p12,int p13,int p21,int p22,int p23);
//double distance(double p11,double p12,double p13,double p21,double p22,double p23);
int get_file(char fnam[],int mode); // Return existing file name
int**** read_info(int**** grid,int dir[],char file[]); // Read grain info.
double*** read_Tmap(double*** temp,int dir[],char file[],double range[],int mode,double Tmelt); // Load temperature distribution info.
int**** alloc4d(int x,int y,int z,int**** grid); // 3DxN_INFO Memory allocation
int**** free4d(int x,int y,int z,int**** grid); // 3DxN_INFO memory de-allocation
double*** alloc3d(int x,int y,int z,double*** grid); // 3D Memory allocation
double*** free3d(int x,int y,int z,double*** grid); // 3D memory de-allocation
void neworient(int grid[]);
void neworient_powder(int grid[]);
void neworient_melting(int grid[]);
void datout(int dir[],int**** grid,int mcs,int out[],double const realtime[],double rt); // Write result.dat file
double lfrac(int**** grid,int dir[]);
void tproout(int dir[],double*** temp,int mcs); // Make MCS.tdat file
double grainsize(int**** const grid,int dir[],char info[],int plane,int cor); // AGI method
void fileout(int dir[],int**** grid,int mcs); // Make MCS.dat file
//void fileoutPrev(int dir[],int**** grid,int mcs); // Make MCS.dat file, previous version
void QS(double arr[], int left, int right); //quick sort
void tok(char* line,char* word); // tok line by " " & "\t"
double grainsize_ff(int**** const grid,int dir[],int mode); // Grain distribution by flood-fill
void DFS_array_simple(int con[],int x,int y,int z,int gori,int**** const grid,int*** checker,int* cnt,int pbc[]); // Flood-fill DFS
void DFS_array(int con[],int x,int y,int z,int gori,int**** const grid,int*** checker,int lim[],int* cnt,int pbc[]); // Flood-fill DFS, max.&min. position check
void IFposition(int**** const grid,int dir[],int cor,int start); // find S-L interface position
void pdas(int**** tgrid,char fnam[],int dir[]);
double MCS_time_acceleration(double dt[],double factor[],double Pmax[],int mode);
int find_keyword(FILE *file,char word[],char option[]);
void find_range(int**** grid,double*** tgrid,int dir[],int pbc[],double temp[]);
double Tcalc_MCS(char func[],int mcs); // calculate T from t(MCS)
double Tcalc_sec(char func[],double time); // calculate T from t(sec)

#endif
