#pragma once

#include <stdio.h>
#include "Utility.h"

void Database(double melt[],double jc[],double rscale[]); // Display information included in material.sdb
int Load_DB(char file[],char sys[],double melt[],double jc[],double rscale[]); // Reading material.sdb for simulation

