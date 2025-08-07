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
void SOLmodule(int**** grid,double jc[],int dir[],int out[],double melt[],int pbc[],double rscale[],double tbc[][2]);
void AMmodule(int**** grid,double jc[],int dir[],int out[],double melt[],int pbc[],double rscale[],double tbc[][2],double am[]);

