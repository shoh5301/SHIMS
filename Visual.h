#include <stdio.h>
#include <sys/io.h>
#include "Utility.h"

#define DRAWMAP "\
#!/bin/bash\n\
gnuplot4 < gtemp.in\n\
if ! test -e XY2.dat\n\
then\n\
rm XY1.dat gtemp.in\n\
else\n\
gnuplot4 < gtemp2.in\n\
rm XY[1-3].dat YZ[1-3].dat XZ[1-3].dat gtemp.in gtemp2.in\n\
fi\n\
"
#define DRAWT "\
#!/bin/bash\n\
gnuplot4 < gT1.in \n\
gnuplot4 < gT2.in \n\
rm tXY[1-3].dat tYZ[1-3].dat tXZ[1-3].dat gT1.in gT2.in\n\
"
#define DRAWPOLE "\
#!/bin/bash\n\
gnuplot4 < gtemp3.in \n\
rm gtemp3.in texture.dat\n\
"

void graphout(int mode,int out[],int dir[],int melt,double rscale[],double frozen); // Graph module
int gnu_plot(int**** tgrid,char fnam[],int dir[],int option); // Draw grain morphology
void tgnu_plot(double*** temp,char fnam[],int dir[],int melt,double haz[]); // Draw temperature distribution
int max_size(int a,int b);//{ // Return bigger value
int angle(int a,int b,int c); // Return colour for specific angle
void text_plot(int**** tgrid,char fnam[],int dir[]); // Draw texture pole figure
double avgtcalc(double*** temp,char fnam[],int dir[],int mode);
void convert_tdat(double*** grid,int dir[]);
void convert_dat(int**** grid,int dir[]);

