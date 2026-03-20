#pragma once

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <sys/io.h>
#include "Utility.h"
#include "Calc.h"
#include "Visual.h"
#include "MC.h"

#define _USE_TIME_DEFINES

void GGmodule(int**** grid,double jc[],int dir[],int out[],char ttem0[],double melt[],int pbc[],double const sang[],double rscale[],double tbc[][2]);// Grain growth
void SOLmodule(int**** grid,double jc[],int dir[],int out[],double melt[],int pbc[],double rscale[],double tbc[][2],char dfeq[]);
void AMmodule(int**** grid,double jc[],int dir[],int out[],double melt[],int pbc[],double rscale[],double tbc[][2],double am[],char dfeq[]);
void heat_transfer(int**** grid,double**** tgrid,double tbc[][2],int dir[],int pbc[],double mcdx,double const melt[],int fdm_loop,double k_per_cp,double Tmelt,int mode);//,double dHfus){
