#pragma once
#include <stdio.h>
#include <time.h>
#include "Utility.h"

void gmake(int dir[],int**** grid); // Grain maker module
void initial_size_setting(int gnum[],int dir[]);
int lcm(int a,int b); // Return LCM value
int**** multiply(int**** grid,int dir[],char fnam[],int multi[],int mod); // Multiply grain structures

