#include "Visual.h"

void graphout(int mode,int out[],int dir[],int melt,double rscale[],double frozen){// Graphics module
	int i=0,j,k,n,module,tdir[3][2]={0},dir2[3]={0},ori[MAXDG1][MAXDG2][MAXDG3]={0};
	int**** tgrid=NULL;
	double haz[2]={0};
	double*** temp=NULL;
	char fnam[30],info[30];
	FILE *tfile=NULL;

	if(mode==0 || mode==2){	// Draw graph after simulation
		for(n=0;n<=out[0];n=n+out[1]){
			sprintf(fnam,"MCS%05d.dat",n);
			tgrid=read_info(tgrid,dir,fnam);
			sprintf(fnam,"MCS%05d",n);
			i=gnu_plot(tgrid,fnam,dir,0);
			tgrid=free4d(dir[0],dir[1],dir[2],tgrid);
		}if(n-out[1]!=out[0]){	// Check unsaved MCSs (last one)
			sprintf(fnam,"MCS%05d.dat",out[0]);
			tgrid=read_info(tgrid,dir,fnam);
			sprintf(fnam,"MCS%05d",out[0]);
			i=gnu_plot(tgrid,fnam,dir,0);
			tgrid=free4d(dir[0],dir[1],dir[2],tgrid);
		}
		if(mode==2 || rscale[0]!=0 && frozen>=0){  // Draw Tmap
			printf("\nDraw T distribution? (GNUPLOT program is required, y=1 / n=0)\n >> ");
			scanf("%d",&n);

			if(n==1){
			    for(n=0;n<=out[0];n=n+out[1]){
				sprintf(fnam,"MCS%05d.tdat",n);
				temp=read_Tmap(temp,dir,fnam,haz,-1,rscale[2]);
				if(temp==NULL){
					printf("No .tdat file detected...\n");
					if(n==0)
					    return;
				}else{
					sprintf(fnam,"MCS%05d",n);
					tgnu_plot(temp,fnam,dir,melt,haz);
					temp=free3d(dir[0],dir[1],dir[2],temp);
					temp=NULL;
				}

			    }
			    if(n-out[1]!=out[0]){	// Check unsaved MCSs (last one)
				sprintf(fnam,"MCS%05d.tdat",out[0]);
				temp=read_Tmap(temp,dir,fnam,haz,-1,rscale[2]);
				sprintf(fnam,"MCS%05d",out[0]);
				tgnu_plot(temp,fnam,dir,melt,haz);
				temp=free3d(dir[0],dir[1],dir[2],temp);
				temp=NULL;
			    }
			}
		}
	}else{		// Load sample & draw graph
		while(1){
			while(1){
printf(" Which module? (1 to draw grain map / 2 to draw T map / 3 to measure grain size ");
printf("/ 4 to measure interface position / 5 to calculate PDAS / 6 to calculate avg. T ");
printf("/ 7 to measure fraction of orientations / 8 to convert .dat (into .vtk) / 9 to convert .tdat (into .vtk) / 0 to quit \n");
// / 5 to draw texture pole figure / 0 to quit)\n);

				printf(" * NOTICE : the script was based on GNUPLOT 4.6 ver. ...\n >> ");
				scanf("%d",&module);
				if(module==0)
					return;
				else if(module>0 && module<=9)
					break;
				printf("Wrong input...\n\n");
					return;
			}
			n=get_file(fnam,0);
			if(n==0)
				return;

			if(module==2 || module==6 || module==9){
			    temp=read_Tmap(temp,dir,fnam,haz,0,rscale[2]);
			    if(module==2){
				printf("Start to draw T map...\n");
				printf("Melting temperature (K, integer)?\n T = ");
				scanf("%d",&n);
				tgnu_plot(temp,fnam,dir,n,haz);
			    }else if(module==6){
				printf("Avg. T across the sample is: %.1f K: maximum %.1f K & minimun %.1f K\n\n",\
					avgtcalc(temp,fnam,dir,0),avgtcalc(temp,fnam,dir,1),avgtcalc(temp,fnam,dir,2));
			    }else if(module==9)
				convert_tdat(temp,dir);

			    temp=free3d(dir[0],dir[1],dir[2],temp);
			    temp=NULL;
			}else{
				tgrid=read_info(tgrid,dir,fnam);

				dir2[0]=dir[0];
				dir2[1]=dir[1];
				dir2[2]=dir[2];
				if(module==1){
					printf("Start to draw grain map...\n");
					printf("  Draw whole sample? (y = 1 / n = 0)\n >> ");
					scanf("%d",&n);
					if(n==0){
						printf("\n Cutting coordinate:\n");
						printf(" X range [0,%d]",dir[0]-1);
						while(1){
							printf("\nX from = ");
							scanf("%d",&i);
							if(i>=0 && i<dir[0]-1)
								break;
							printf("Wrong initial coordinate: not in the range [0,%d)\n",dir[0]-1);
						}while(1){
							printf("X to = ");
							scanf("%d",&j);
							if(j>i && j<dir[0])
								break;
							printf("Wrong initial coordinate: not in the range (%d,%d]\n",i,dir[0]-1);
						}
						tdir[0][0]=i; //min
						tdir[0][1]=j; //max
						printf("\n Y range [0,%d]",dir[1]-1);
						while(1){
							printf("\nY from = ");
							scanf("%d",&i);
							if(i>=0 && i<dir[1]-1)
								break;
							printf("Wrong initial coordinate: not in the range [0,%d)\n",dir[1]-1);
						}while(1){
							printf("Y to = ");
							scanf("%d",&j);
							if(j>i && j<dir[1])
								break;
							printf("Wrong initial coordinate: not in the range (%d,%d]\n",i,dir[1]-1);
						}
						tdir[1][0]=i; //min
						tdir[1][1]=j; //max
						printf("\n Z range [0,%d]",dir[2]-1);
						while(1){
							printf("\nZ from = ");
							scanf("%d",&i);
							if(i>=0 && i<dir[2]-1)
								break;
							printf("Wrong initial coordinate: not in the range [0,%d)\n",dir[2]-1);
						}while(1){
							printf("Z to = ");
							scanf("%d",&j);
							if(j>i && j<dir[2])
								break;
							printf("Wrong initial coordinate: not in the range (%d,%d]\n",i,dir[2]-1);
						}
						tdir[2][0]=i; //min
						tdir[2][1]=j; //max
						printf("\n X range [%d, %d]\n",tdir[0][0],tdir[0][1]);
						printf(" Y range [%d, %d]\n",tdir[1][0],tdir[1][1]);
						printf(" Z range [%d, %d]\n\n",tdir[2][0],tdir[2][1]);

						dir2[0]=tdir[0][1]-tdir[0][0]+1;
						dir2[1]=tdir[1][1]-tdir[1][0]+1;
						dir2[2]=tdir[2][1]-tdir[2][0]+1;

						for(i=tdir[0][0];i<=tdir[0][1];i++){
							for(j=tdir[1][0];j<=tdir[1][1];j++){
								for(k=tdir[2][0];k<=tdir[2][1];k++){
									tgrid[i-tdir[0][0]][j-tdir[1][0]][k-tdir[2][0]][0]=tgrid[i][j][k][0];
									tgrid[i-tdir[0][0]][j-tdir[1][0]][k-tdir[2][0]][1]=tgrid[i][j][k][1];
									tgrid[i-tdir[0][0]][j-tdir[1][0]][k-tdir[2][0]][2]=tgrid[i][j][k][2];
						}	}	}
					}
					printf("Draw specific plane or multiple planes? (1 to specific, 0 to multiple planes)\n >> ");
					scanf("%d",&n);
					i=gnu_plot(tgrid,fnam,dir2,n);
					if(n==1&&i==0)
						printf("\"%s.*.eps\" file was created...\n\n",fnam,fnam);
					else if(i==0)
						printf("\"%s.eps\" & \"3d_%s.eps\" file was created...\n\n",fnam,fnam);
				}else if(module==4){
					printf("Solid growth direction? (+1 = +x, -1 = -x)\n >> ");
					scanf("%d",&i);
					IFposition(tgrid,dir,i,j);
				}else if(module==5){
					pdas(tgrid,fnam,dir);
//				}else if(module==5){	// texture plot
//					text_plot(tgrid,fnam,dir);
				}else if(module==7){	// orientation reading
					n=0;
					for(i=0;i<dir[0];i++){
						for(j=0;j<dir[1];j++){
							for(k=0;k<dir[2];k++){
								if(tgrid[i][j][k][0]!=LIQORI){
									ori[tgrid[i][j][k][0]][tgrid[i][j][k][1]][tgrid[i][j][k][2]]++;
									n++;
								}
					}	}	}
					tfile=fopen("Orient.dat","w");
					fprintf(tfile,"Angle\tNumber of voxels\tFraction\n");
					for(i=0;i<MAXDG1;i++){
						for(j=0;j<MAXDG2;j++){
							for(k=0;k<MAXDG3;k++){
								fprintf(tfile,"( %3d %3d %3d )\t",i,j,k);
								if(ori[i][j][k]!=0)
									fprintf(tfile,"%d\t\t%.2E\n",ori[i][j][k],ori[i][j][k]/(double)n);
								else
									fprintf(tfile,"%d\t\t0.00\n",ori[i][j][k]);
					}	}	}
					fclose(tfile);
					printf("\n * Orient.dat file is created...\n\n");
				}else if(module==8){
				    convert_dat(tgrid,dir);
				}else if(module==3){
				    printf(" Which method? (1 = AGI method / 2 = Flood-fill method)\n");
				    printf(" ### Note: for computational efficiency, AGI method is utilized during MC simulations\n >> ");
				    scanf("%d",&i);

				    if(i==1){ // AGI
					while(1){
						printf("Grain size measurement by AGI Method\n Z Position? (from 0 to %d // -1 to quit) >> ",dir[2]-1);
						scanf("%d",&out[2]);
						if(out[2]<0)
							break;
						else if(out[2]>=0 && out[2]<dir[2]){
							grainsize(tgrid,dir,info,out[2],2);
							printf("%s\n",info);
						}else
							printf("Out of range...\n");
					}
					info[0]='\0';
				    }else // Flood-Fill
					i=(int)grainsize_ff(tgrid,dir,1);
				}
				tgrid=free4d(dir[0],dir[1],dir[2],tgrid);
			}
			printf("Please input other...\n");
	}	}
	return;
}

int gnu_plot(int**** tgrid,char fnam[],int dir[],int option){// Draw grain morphology
	FILE* gpl=fopen("gtemp.in","w");
	FILE* tfile=NULL;
	int i,j,center,max;
	double maxdg=0;

	maxdg=((MAXDG1+1)*MAXDG2+1)*MAXDG3;

	for(i=0;i<dir[0];i++){
	    for(j=0;j<dir[1];j++){
		for(center=0;center<dir[2];center++){
		    if(tgrid[i][j][center][0]==LIQORI){ // liquid
			tgrid[i][j][center][0]=LIQORI;//MAXDG1+1;
			tgrid[i][j][center][1]=LIQORI;//MAXDG2+1;
			tgrid[i][j][center][2]=LIQORI;//MAXDG3+2;
		    }else if(tgrid[i][j][center][1]==POWORI){
			tgrid[i][j][center][0]=POWORI;//MAXDG1+1;
			tgrid[i][j][center][1]=POWORI;//MAXDG2+1;
			tgrid[i][j][center][2]=POWORI;//MAXDG3+2;
		    }
		}
	    }	
	}

// 2D map
fprintf(gpl,"set term postscript portrait \"Times-New-Roman\" 10 background rgb 'white'\n");
fprintf(gpl,"set pm3d map\nset style data lines\nset nokey\n");
fprintf(gpl,"set xtics border in scale 0,0 mirror norotate  autojustify\n");
fprintf(gpl,"set ytics border in scale 0,0 mirror norotate  autojustify\n");
fprintf(gpl,"set rtics axis in scale 0,0 nomirror norotate  autojustify\n");
fprintf(gpl,"set cblabel \"Grain Orientation\"\n");
fprintf(gpl,"set palette defined (%.1f \"red\", -0.01 \"black\", %.1f \"blue\", %.1f \"green\",%.1f \"orange\",%.1f \"yellow\", %.1f \"grey\")\n",LIQORI-1.0,0.25*maxdg,0.5*maxdg,0.75*maxdg,maxdg,(double)POWORI);
fprintf(gpl,"set cbtics (\"(0,0,0)\" 0,\"(%d,%d,%d)\" %.1f)\n",MAXDG1,MAXDG2,MAXDG3,maxdg);
fprintf(gpl,"set cbrange [ 0 : %.1f ] #noreverse nowriteback\n",maxdg);

	if(option==1){
	    fprintf(gpl,"set size 1.5,1.0\n");

	    printf("Which plane? (XY = 1, YZ = 2, XZ = 3)\n >> ");
	    scanf("%d",&i);
	  switch(i){
	    case 1:
		printf("Which plane? (0 ~ %d)\n >> ",dir[2]-1);
		scanf("%d",&center);
		if(center>=dir[2]){
			printf("Wrong input... Quitting the module\n\n");
			return 1;
		}
		fprintf(gpl,"set output '%s.xy.z%d.eps'\n",fnam,center);
		tfile=fopen("XY1.dat","w");
		for(i=0;i<dir[0];i++){
		  for(j=0;j<dir[1];j++)
		    fprintf(tfile,"%d\t%d\t%d\n",i+1,j+1,angle(tgrid[i][j][center][0],tgrid[i][j][center][1],tgrid[i][j][center][2]));
		}
		max=max_size(dir[1],dir[0]);
		break;
	    case 2:
		printf("which plane? (0 ~ %d)\n >> ",dir[0]-1);
		scanf("%d",&center);
		if(center>=dir[0]){
			printf("wrong input... quitting the module\n\n");
			return 1;
		}
		fprintf(gpl,"set output '%s.yz.x%d.eps'\n",fnam,center);
		tfile=fopen("XY1.dat","w");
		for(i=0;i<dir[1];i++){
		  for(j=0;j<dir[2];j++)
		    fprintf(tfile,"%d\t%d\t%d\n",i+1,j+1,angle(tgrid[center][i][j][0],tgrid[center][i][j][1],tgrid[center][i][j][2]));
		}
		max=max_size(dir[1],dir[2]);
		break;
	    case 3:
		printf("which plane? (0 ~ %d)\n >> ",dir[1]-1);
		scanf("%d",&center);
		if(center>=dir[1]){
			printf("wrong input... quitting the module\n\n");
			return 1;
		}
		fprintf(gpl,"set output '%s.xz.y%d.eps'\n",fnam,center);
		tfile=fopen("XY1.dat","w");
		for(i=0;i<dir[0];i++){
		  for(j=0;j<dir[2];j++)
		    fprintf(tfile,"%d\t%d\t%d\n",i+1,j+1,angle(tgrid[i][center][j][0],tgrid[i][center][j][1],tgrid[i][center][j][2]));
		}
		max=max_size(dir[0],dir[2]);
		break;
	    default:
		printf("Wrong input... Quitting the module\n\n");
		return 1;
		break;
	  }
	  fprintf(gpl,"set xrange[0:%d]\nset yrange[0:%d]\n",max,max);
	  fprintf(gpl,"splot 'XY1.dat' u 1:2:3 w image\n");
	  fclose(tfile);
	  fclose(gpl);
	}else{
fprintf(gpl,"set size 1.8,1.0\n");
fprintf(gpl,"set output '%s.eps'\n",fnam);
//fprintf(gpl,"set multiplot\n");
fprintf(gpl,"set multiplot title \"%s\"\n",fnam);
fprintf(gpl,"unset colorbox\n");

max=max_size(dir[0],dir[1]);
fprintf(gpl,"set xrange[0:%d]\nset yrange[0:%d]\n",max,max);	// X-Y
fprintf(gpl,"set title \"Upper X-Y Surface\"\nset size 0.6,0.4\n");
fprintf(gpl,"set origin 0,0.6\nsplot 'XY2.dat' u 1:2:4 w image\n");	// Upper XY Surface
fprintf(gpl,"set title \"Centre X-Y Surface\"\nset size 0.6,0.4\n");
fprintf(gpl,"set origin 0,0.3\nsplot 'XY3.dat' u 1:2:4 w image\n");	// Centre XY Surface
fprintf(gpl,"set title \"Lower X-Y Surface\"\nset size 0.6,0.4\n");
fprintf(gpl,"set origin 0,0\nsplot 'XY1.dat' u 1:2:4 w image\n");	// Lower XY Surface

max=max_size(dir[2],dir[1]);
fprintf(gpl,"set xrange[0:%d]\nset yrange[0:%d]\n",max,max);	// Y-Z
fprintf(gpl,"set title \"Upper Y-Z Surface\"\nset size 0.6,0.4\n");
fprintf(gpl,"set origin 0.5,0.6\nsplot 'YZ2.dat' u 2:3:4 w image\n");	// Upper YZ Surface
fprintf(gpl,"set title \"Centre Y-Z Surface\"\nset size 0.6,0.4\n");
fprintf(gpl,"set origin 0.5,0.3\nsplot 'YZ3.dat' u 2:3:4 w image\n");	// Centre YZ Surface
fprintf(gpl,"set title \"Lower Y-Z Surface\"\nset size 0.6,0.4\n");
fprintf(gpl,"set origin 0.5,0\nsplot 'YZ1.dat' u 2:3:4 w image\n");	// Lower YZ Surface

max=max_size(dir[2],dir[0]);
fprintf(gpl,"set colorbox\n");
fprintf(gpl,"set xrange[0:%d]\nset yrange[0:%d]\n",max,max);	// Z-X
fprintf(gpl,"set title \"Upper Z-X Surface\"\nset size 0.6,0.4\n");
fprintf(gpl,"set origin 1,0.6\nsplot 'XZ2.dat' u 1:3:4 w image\n");	// Upper XZ Surface
fprintf(gpl,"set title \"Centre Z-X Surface\"\nset size 0.6,0.4\n");
fprintf(gpl,"set origin 1,0.3\nsplot 'XZ3.dat' u 1:3:4 w image\n");	// Centre XZ Surface
fprintf(gpl,"set title \"Lower Z-X Surface\"\nset size 0.6,0.4\n");
fprintf(gpl,"set origin 1,0\nsplot 'XZ1.dat' u 1:3:4 w image\n");	// Lower XZ Surface
fclose(gpl);
/*	Lower	Upper
	XY1	XY2
	YZ1	YZ2
	XZ1	XZ2	*/
	tfile=fopen("XY2.dat","w");	//Upper
	for(i=0;i<dir[0];i++){
		for(j=0;j<dir[1];j++)
			fprintf(tfile,"%d\t%d\t%d\t%d\n",i+1,j+1,dir[2]+1,angle(tgrid[i][j][dir[2]-1][0],tgrid[i][j][dir[2]-1][1],tgrid[i][j][dir[2]-1][2]));
	}fclose(tfile);
	tfile=fopen("XY3.dat","w");	//Center
	center=(int)((dir[2]-1)/2);
	for(i=0;i<dir[0];i++){
		for(j=0;j<dir[1];j++)
			fprintf(tfile,"%d\t%d\t%d\t%d\n",i+1,j+1,center,angle(tgrid[i][j][center][0],tgrid[i][j][center][1],tgrid[i][j][center][2]));
	}fclose(tfile);
	tfile=fopen("XY1.dat","w");	//Lower
	for(i=0;i<dir[0];i++){
		for(j=0;j<dir[1];j++)
			fprintf(tfile,"%d\t%d\t%d\t%d\n",i+1,j+1,0,angle(tgrid[i][j][0][0],tgrid[i][j][0][1],tgrid[i][j][0][2]));
	}fclose(tfile);

	tfile=fopen("YZ2.dat","w");	//Upper
	for(i=0;i<dir[1];i++){
		for(j=0;j<dir[2];j++)
			fprintf(tfile,"%d\t%d\t%d\t%d\n",dir[0]+1,i+1,j+1,angle(tgrid[dir[0]-1][i][j][0],tgrid[dir[0]-1][i][j][1],tgrid[dir[0]-1][i][j][2]));
	}fclose(tfile);
	tfile=fopen("YZ3.dat","w");	//Center
	center=(int)((dir[0]-1)/2);
	for(i=0;i<dir[1];i++){
		for(j=0;j<dir[2];j++){
			fprintf(tfile,"%d\t%d\t%d\t%d\n",center,i+1,j+1,angle(tgrid[center][i][j][0],tgrid[center][i][j][1],tgrid[center][i][j][2]));
		}
	}fclose(tfile);
	tfile=fopen("YZ1.dat","w");	//Lower
	for(i=0;i<dir[1];i++){
		for(j=0;j<dir[2];j++)
			fprintf(tfile,"%d\t%d\t%d\t%d\n",0,i+1,j+1,angle(tgrid[0][i][j][0],tgrid[0][i][j][1],tgrid[0][i][j][2]));
	}fclose(tfile);

	tfile=fopen("XZ2.dat","w");	//Upper
	for(i=0;i<dir[0];i++){
		for(j=0;j<dir[2];j++)
			fprintf(tfile,"%d\t%d\t%d\t%d\n",i+1,dir[1]+1,j+1,angle(tgrid[i][dir[1]-1][j][0],tgrid[i][dir[1]-1][j][1],tgrid[i][dir[1]-1][j][2]));
	}fclose(tfile);
	tfile=fopen("XZ3.dat","w");	//Center
	center=(int)((dir[1]-1)/2);
	for(i=0;i<dir[0];i++){
		for(j=0;j<dir[2];j++)
			fprintf(tfile,"%d\t%d\t%d\t%d\n",i+1,center,j+1,angle(tgrid[i][center][j][0],tgrid[i][center][j][1],tgrid[i][center][j][2]));
	}fclose(tfile);
	tfile=fopen("XZ1.dat","w");	//Lower
	for(i=0;i<dir[0];i++){
		for(j=0;j<dir[2];j++)
			fprintf(tfile,"%d\t%d\t%d\t%d\n",i+1,0,j+1,angle(tgrid[i][0][j][0],tgrid[i][0][j][1],tgrid[i][0][j][2]));
	}fclose(tfile);
//3D map
gpl=fopen("gtemp2.in","w");
fprintf(gpl,"set term postscript portrait \"Times-New-Roman\" 14 background rgb 'white'\n");
fprintf(gpl,"set output '3d_%s.eps'\n",fnam);
fprintf(gpl,"set size 1.5,1.0\n");
fprintf(gpl,"set title \"%s\"\n",fnam);
fprintf(gpl,"set pm3d\nset style data lines\nset nokey\n");
max=max_size(dir[0],max_size(dir[1],dir[2]));
fprintf(gpl,"set xrange[0:%d]\nset yrange[0:%d]\nset zrange[0:%d]\n",max+2,max+2,max+2);
fprintf(gpl,"set xtics border in scale 0,0 mirror norotate# autojustify\n");
fprintf(gpl,"set ytics border in scale 0,0 mirror norotate# autojustify\n");
fprintf(gpl,"set ztics border in scale 0,0 mirror norotate# autojustify\n");
fprintf(gpl,"set rtics axis in scale 0,0 nomirror norotate autojustify\n");
fprintf(gpl,"set cblabel \"Grain Orientation\"\n");
fprintf(gpl,"set cbrange [ 0:%.1f ] #noreverse nowriteback\n",maxdg);
fprintf(gpl,"set palette defined (%.1f \"red\", -0.01 \"black\", %.1f \"blue\", %.1f \"green\",%.1f \"orange\",%.1f \"yellow\", %.1f \"grey\")\n",LIQORI-1.0,0.25*maxdg,0.5*maxdg,0.75*maxdg,maxdg,(double)POWORI);
fprintf(gpl,"set cbtics (\"(0,0,0)\" 0,\"(%d,%d,%d)\" %.1f)\n",MAXDG1,MAXDG2,MAXDG3,maxdg);
fprintf(gpl,"set xyplane at 0\n#set view 65, 145, 1, 1 #set hidden3d\n");
fprintf(gpl,"splot 'XY2.dat' u 1:2:3:4 w image,'XZ1.dat' u 1:2:3:4 w image,'YZ2.dat' u 1:2:3:4 w image");
fclose(gpl);
  }
	system(DRAWMAP);
	return 0;
}

void tgnu_plot(double*** temp,char fnam[],int dir[],int melt,double haz[]){// Draw temperature distribution
	FILE* gpl=fopen("gT1.in","w");
	FILE* tfile=NULL;
	int i,j,center,max;
	double mint=0,maxt=0;
/*	for(i=0;i<dir[0];i++){
		for(j=0;j<dir[1];j++){
			for(center=0;center<dir[2];center++){
				if(temp[i][j][center]<0)
					temp[i][j][center]=melt;
				else if(temp[i][j][center]>0)
					temp[i][j][center]=(int)(melt*(double)temp[i][j][center]/haz);
	}	}	}*/

	if((int)haz[1] < melt-50)
		maxt=(double)melt+100;
	else
		maxt=haz[1]+100;
	if((int)haz[0] < 100)
		mint=0;
	else
		mint=haz[0]-100;

fprintf(gpl,"set term postscript portrait \"Times-New-Roman\" 10 background rgb 'white'\n");
fprintf(gpl,"set output 'Tmap_%s.eps'\n",fnam);
fprintf(gpl,"set size 1.8,1.0\n");
fprintf(gpl,"set pm3d map\nset style data lines\nset nokey\n");
fprintf(gpl,"set xtics border in scale 0,0 mirror norotate  autojustify\n");
fprintf(gpl,"set ytics border in scale 0,0 mirror norotate  autojustify\n");
fprintf(gpl,"set rtics axis in scale 0,0 nomirror norotate  autojustify\n");

fprintf(gpl,"set cblabel \"Temperature\"\n");
fprintf(gpl,"set cbtics (\"%d K\" %f, \"Tm(%d K)\" %f, \"%d K\" %f)\n",(int)(0.5*melt),0.5*melt,(int)melt,melt,(int)maxt,maxt);
fprintf(gpl,"set cbrange [ %.1f : %.1f ] noreverse nowriteback\n",0.5*melt,maxt);
fprintf(gpl,"set palette defined (0 \"blue\", 1 \"green\", 2 \"yellow\", 3 \"orange\", 4\"red\")\n");
fprintf(gpl,"set multiplot title \"%s\"\n",fnam);
//fprintf(gpl,"set multiplot\n");
fprintf(gpl,"unset colorbox\n");

max=max_size(dir[0],dir[1]);
fprintf(gpl,"set xrange[0:%d]\nset yrange[0:%d]\n",max,max);	// X-Y
fprintf(gpl,"set title \"Upper X-Y Surface\"\nset size 0.6,0.4\n");
fprintf(gpl,"set origin 0,0.6\nsplot 'tXY2.dat' u 1:2:4 w image\n");	// Upper XY Surface
fprintf(gpl,"set title \"Centre X-Y Surface\"\nset size 0.6,0.4\n");
fprintf(gpl,"set origin 0,0.3\nsplot 'tXY3.dat' u 1:2:4 w image\n");	// Centre XY Surface
fprintf(gpl,"set title \"Lower X-Y Surface\"\nset size 0.6,0.4\n");
fprintf(gpl,"set origin 0,0\nsplot 'tXY1.dat' u 1:2:4 w image\n");	// Lower XY Surface

max=max_size(dir[2],dir[1]);
fprintf(gpl,"set xrange[0:%d]\nset yrange[0:%d]\n",max,max);	// Y-Z
fprintf(gpl,"set title \"Upper Y-Z Surface\"\nset size 0.6,0.4\n");
fprintf(gpl,"set origin 0.5,0.6\nsplot 'tYZ2.dat' u 2:3:4 w image\n");	// Upper YZ Surface
fprintf(gpl,"set title \"Centre Y-Z Surface\"\nset size 0.6,0.4\n");
fprintf(gpl,"set origin 0.5,0.3\nsplot 'tYZ3.dat' u 2:3:4 w image\n");	// Centre YZ Surface
fprintf(gpl,"set title \"Lower Y-Z Surface\"\nset size 0.6,0.4\n");
fprintf(gpl,"set origin 0.5,0\nsplot 'tYZ1.dat' u 2:3:4 w image\n");	// Lower YZ Surface

max=max_size(dir[2],dir[0]);
fprintf(gpl,"set colorbox\n");
fprintf(gpl,"set xrange[0:%d]\nset yrange[0:%d]\n",max,max);	// Z-X
fprintf(gpl,"set title \"Upper Z-X Surface\"\nset size 0.6,0.4\n");
fprintf(gpl,"set origin 1,0.6\nsplot 'tXZ2.dat' u 1:3:4 w image\n");	// Upper XZ Surface
fprintf(gpl,"set title \"Centre Z-X Surface\"\nset size 0.6,0.4\n");
fprintf(gpl,"set origin 1,0.3\nsplot 'tXZ3.dat' u 1:3:4 w image\n");	// Centre XZ Surface
fprintf(gpl,"set title \"Lower Z-X Surface\"\nset size 0.6,0.4\n");
fprintf(gpl,"set origin 1,0\nsplot 'tXZ1.dat' u 1:3:4 w image\n");	// Lower XZ Surface
fclose(gpl);

/*	Lower	Upper
	XY1	XY2
	YZ1	YZ2
	XZ1	XZ2	*/
	tfile=fopen("tXY2.dat","w");	//Upper
	for(i=0;i<dir[0];i++){
		for(j=0;j<dir[1];j++)
			fprintf(tfile,"%d\t%d\t%d\t%f\n",i+1,j+1,dir[2]+1,temp[i][j][dir[2]-1]);
	}fclose(tfile);
	tfile=fopen("tXY3.dat","w");	//Center
	center=(int)((dir[2]-1)/2);
	for(i=0;i<dir[0];i++){
		for(j=0;j<dir[1];j++)
			fprintf(tfile,"%d\t%d\t%d\t%f\n",i+1,j+1,center,temp[i][j][center]);
	}fclose(tfile);
	tfile=fopen("tXY1.dat","w");	//Lower
	for(i=0;i<dir[0];i++){
		for(j=0;j<dir[1];j++)
			fprintf(tfile,"%d\t%d\t%d\t%f\n",i+1,j+1,0,temp[i][j][0]);
	}fclose(tfile);

	tfile=fopen("tYZ2.dat","w");	//Upper
	for(i=0;i<dir[1];i++){
		for(j=0;j<dir[2];j++)
			fprintf(tfile,"%d\t%d\t%d\t%f\n",dir[0]+1,i+1,j+1,temp[dir[0]-1][i][j]);
	}fclose(tfile);
	center=(int)((dir[0]-1)/2);
	tfile=fopen("tYZ3.dat","w");	//Center
	for(i=0;i<dir[1];i++){
		for(j=0;j<dir[2];j++)
			fprintf(tfile,"%d\t%d\t%d\t%f\n",center,i+1,j+1,temp[center][i][j]);
	}fclose(tfile);
	tfile=fopen("tYZ1.dat","w");	//Lower
	for(i=0;i<dir[1];i++){
		for(j=0;j<dir[2];j++)
			fprintf(tfile,"%d\t%d\t%d\t%f\n",0,i+1,j+1,temp[0][i][j]);
	}fclose(tfile);

	tfile=fopen("tXZ2.dat","w");	//Upper
	for(i=0;i<dir[0];i++){
		for(j=0;j<dir[2];j++)
			fprintf(tfile,"%d\t%d\t%d\t%f\n",i+1,dir[1]+1,j+1,temp[i][dir[1]-1][j]);
	}fclose(tfile);
	center=(int)((dir[1]-1)/2);
	tfile=fopen("tXZ3.dat","w");	//Lower
	for(i=0;i<dir[0];i++){
		for(j=0;j<dir[2];j++)
			fprintf(tfile,"%d\t%d\t%d\t%f\n",i+1,center,j+1,temp[i][center][j]);
	}fclose(tfile);
	tfile=fopen("tXZ1.dat","w");	//Lower
	for(i=0;i<dir[0];i++){
		for(j=0;j<dir[2];j++)
			fprintf(tfile,"%d\t%d\t%d\t%f\n",i+1,0,j+1,temp[i][0][j]);
	}fclose(tfile);

gpl=fopen("gT2.in","w");
fprintf(gpl,"set term postscript portrait \"Times-New-Roman\" 14 background rgb 'white'\n");
fprintf(gpl,"set output '3d_Tmap_%s.eps'\n",fnam);
fprintf(gpl,"set size 1.5,1.0\n");
fprintf(gpl,"set title \"%s\"\n",fnam);
fprintf(gpl,"set pm3d\nset style data lines\nset nokey\n");
max=max_size(dir[0],max_size(dir[1],dir[2]));
fprintf(gpl,"set xrange[0:%d]\nset yrange[0:%d]\nset zrange[0:%d]\n",max+1,max+1,max+1);
fprintf(gpl,"set xtics border in scale 0,0 mirror norotate  autojustify\n");
fprintf(gpl,"set ytics border in scale 0,0 mirror norotate  autojustify\n");
fprintf(gpl,"set ztics border in scale 0,0 mirror norotate  autojustify\n");
fprintf(gpl,"set rtics axis in scale 0,0 nomirror norotate  autojustify\n");
fprintf(gpl,"set cblabel \"Temperature\"\n");
fprintf(gpl,"set cbtics (\"%d K\" %f, \"Tm(%d K)\" %f, \"%d K\" %f)\n",(int)(0.5*melt),0.5*melt,(int)melt,melt,(int)maxt,maxt);
fprintf(gpl,"set cbrange [ %.1f : %.1f ] noreverse nowriteback\n",0.5*melt,maxt);
fprintf(gpl,"set palette defined (0 \"blue\", 1 \"green\", 2 \"yellow\", 3 \"orange\", 4\"red\")\n");
fprintf(gpl,"set xyplane at 0\n");
fprintf(gpl,"splot 'tXY2.dat' u 1:2:3:4 w image, 'tYZ2.dat' u 1:2:3:4 w image, 'tXZ1.dat' u 1:2:3:4 w image");
fclose(gpl);

	system(DRAWT);
	return;
}

int max_size(int a,int b){// Return bigger value
	return (a>b)? a:b;
}

int angle(int a,int b,int c){// Return colour for specific angle
	return (a>=0)? a*MAXDG2*MAXDG3+b*MAXDG3+c : LIQORI;
}

void text_plot(int**** tgrid,char fnam[],int dir[]){ // Draw texture pole figure
	FILE* gpl=fopen("gtemp3.in","w");
	FILE* tfile=fopen("texture.dat","w");
	int i,j,k,l=0,deg[MAXDG1][MAXDG2][MAXDG3]={0},m,n;

	for(i=0;i<dir[0];i++){
		for(j=0;j<dir[1];j++){
			for(k=0;k<dir[2];k++){
				if(tgrid[i][j][k][0]>0)
					deg[tgrid[i][j][k][0]][tgrid[i][j][k][1]][tgrid[i][j][k][2]]++;
	}	}	}

	for(i=0;i<MAXDG1;i++){
		for(j=0;j<MAXDG2;j++){
//			if(deg[i][j]!=0)
				fprintf(tfile,"%f\t%f\t%d\n",sin(DTOR*(double)j)*cos(DTOR*(double)i),sin(DTOR*(double)i)*sin(DTOR*(double)j),deg[i][j]);//,i,j);
		}
	}

	fprintf(gpl,"set term jpeg size 680,640\n");
	fprintf(gpl,"set output 'poleTexture_%s.eps'\n",fnam);
	fprintf(gpl,"set dgrid3d\nset contour base\nset cntrparam bspline\nset param\n");
	fprintf(gpl,"set view map\nunset surf\nset style data lp\nunset key\nset xrange [-1:1]\nset yrange [-1:1]\n");
	fprintf(gpl,"set multiplot\nset lmarg at scr 0.05\nset rmarg at scr 0.95\nset tmarg at scr 0.95\nset bmarg at scr 0.05\n");
	fprintf(gpl,"splot 'texture.dat' u 1:2:3\nunset xtics\nunset ytics\nplot cos(t),sin(t)");

	fclose(gpl);
	fclose(tfile);

	system(DRAWPOLE);

	printf("Texture pole figure 'Pole_%s.eps' was created..\n\n",fnam);

	return;
}

double avgtcalc(double*** temp,char fnam[],int dir[],int mode){
	int i,j,k;
	double avg=0;
	
    if(mode==0){
	for(i=0;i<dir[0];i++){
		for(j=0;j<dir[1];j++){
			for(k=0;k<dir[2];k++)
				avg+=temp[i][j][k];
	}	}
	return avg/(double)(dir[0]*dir[1]*dir[2]);
    }
    else if(mode==1){
	for(i=0;i<dir[0];i++){
		for(j=0;j<dir[1];j++){
			for(k=0;k<dir[2];k++){
				if(avg<temp[i][j][k])
				    avg=temp[i][j][k];
	}	}	}
	return avg;
    }else{ // if mode==2
	avg=temp[0][0][0];
	for(i=0;i<dir[0];i++){
		for(j=0;j<dir[1];j++){
			for(k=0;k<dir[2];k++){
				if(avg>temp[i][j][k])
				    avg=temp[i][j][k];
	}	}	}
	return avg;
    }
}

void convert_tdat(double*** grid,int dir[]){
	int i,j,k,n;
	char name[50]={'\0'};
	FILE *vtkfile=NULL;
	
	printf("  Output file name?\n >> ");
	scanf("%s",name);
	vtkfile=fopen(name,"w");

	fprintf(vtkfile,"# vtk DataFile Version 3.0\nConverted temperature from SHIM\nASCII\nDATASET POLYDATA\nPOINTS %d float\n",dir[0]*dir[1]*dir[2]);
	for(i=0;i<dir[0];i++){
	    for(j=0;j<dir[1];j++){
		for(k=0;k<dir[2];k++)
		    fprintf(vtkfile,"%d %d %d\n",i,j,k);
	    }
	}
	fprintf(vtkfile,"VERTICES %d %d\n",dir[0]*dir[1]*dir[2],2*dir[0]*dir[1]*dir[2]);
	n=0;
	for(i=0;i<dir[0];i++){
	    for(j=0;j<dir[1];j++){
		for(k=0;k<dir[2];k++)
		    fprintf(vtkfile,"1 %d\n",n++);
	    }
	}
	fprintf(vtkfile,"POINT_DATA %d\nSCALARS Temperature float 1\nLOOKUP_TABLE default\n",dir[0]*dir[1]*dir[2]);
	for(i=0;i<dir[0];i++){
	    for(j=0;j<dir[1];j++){
		for(k=0;k<dir[2];k++)
		    fprintf(vtkfile,"%f\n",grid[i][j][k]);
	    }
	}

	return;
}

void convert_dat(int**** grid,int dir[]){
	int i,j,k,n;
	char name[50]={'\0'};
	FILE *vtkfile=NULL;
	
	printf("  Output file name?\n >> ");
	scanf("%s",name);
	vtkfile=fopen(name,"w");

	fprintf(vtkfile,"# vtk DataFile Version 3.0\nConverted grain morphology from SHIM\nASCII\nDATASET POLYDATA\nPOINTS %d float\n",dir[0]*dir[1]*dir[2]);
	for(i=0;i<dir[0];i++){
	    for(j=0;j<dir[1];j++){
		for(k=0;k<dir[2];k++)
		    fprintf(vtkfile,"%d %d %d\n",i,j,k);
	    }
	}
	fprintf(vtkfile,"VERTICES %d %d\n",dir[0]*dir[1]*dir[2],2*dir[0]*dir[1]*dir[2]);
	n=0;
	for(i=0;i<dir[0];i++){
	    for(j=0;j<dir[1];j++){
		for(k=0;k<dir[2];k++)
		    fprintf(vtkfile,"1 %d\n",n++);
	    }
	}
	fprintf(vtkfile,"POINT_DATA %d\nSCALARS Orientation_colorized float 1\nLOOKUP_TABLE default\n",dir[0]*dir[1]*dir[2]);
	for(i=0;i<dir[0];i++){
	    for(j=0;j<dir[1];j++){
		for(k=0;k<dir[2];k++)
		    fprintf(vtkfile,"%d\n",angle(grid[i][j][k][0],grid[i][j][k][1],grid[i][j][k][2]));
	    }
	}
	
	fclose(vtkfile);
	return;
}

