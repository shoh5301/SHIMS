#include "Utility.h"

int prob_check(double prob){
	if(prob<PROBMIN)
	    return 0;
	else
	    return (prob>(double)rand()/(RAND_MAX+1.0))? 1:0;
//	    return (prob>(double)(rand()+1)/(RAND_MAX+1.0))? 1:0;
}

void voxcpy(int paste[],int copy[]){
	paste[0]=copy[0];
	paste[1]=copy[1];
	paste[2]=copy[2];
	return;
}
/*
void nei6sel(int i,int *x,int *y, int *z,int dir[],int pbc[]){
	switch(i){ // +x, -x, +y, -y, +z, -z
		case 0 :
			*x=deter(*x-1,dir[0],pbc[0]);
			break;
		case 1 :
			*x=deter(*x+1,dir[0],pbc[0]);
			break;
		case 2 :
			*y=deter(*y-1,dir[1],pbc[1]);
			break;
		case 3 :
			*y=deter(*y+1,dir[1],pbc[1]);
			break;
		case 4 :
			*z=deter(*z-1,dir[2],pbc[2]);
			break;
		case 5 :
			*z=deter(*z+1,dir[2],pbc[2]);
			break;
		default:
puts("ERROR: wrong case for neigh6sel function\n");
			break;
	}
}*/

void nei26sel(int i,int *x,int *y,int *z,int dir[],int pbc[]){
	switch(i){
		case 0 :  // x-1 y-1 z-1
			*x=deter(*x-1,dir[0],pbc[0]);
			*y=deter(*y-1,dir[1],pbc[1]);
			*z=deter(*z-1,dir[2],pbc[2]);
			break;
		case 1 :  // x-1 y-1 z
			*x=deter(*x-1,dir[0],pbc[0]);
			*y=deter(*y-1,dir[1],pbc[1]);
			break;
		case 2 :  // x-1 y-1 z+1
			*x=deter(*x-1,dir[0],pbc[0]);
			*y=deter(*y-1,dir[1],pbc[1]);
			*z=deter(*z+1,dir[2],pbc[2]);
                        break;
		case 3 :  // x-1 y   z-1
       			*x=deter(*x-1,dir[0],pbc[0]);
			*z=deter(*z-1,dir[2],pbc[2]);
                        break;
		case 4 :  // x-1 y   z
       			*x=deter(*x-1,dir[0],pbc[0]);
                        break;
		case 5 :  // x-1 y   z+1
       			*x=deter(*x-1,dir[0],pbc[0]);
			*z=deter(*z+1,dir[2],pbc[2]);
                        break;
		case 6 :  // x-1 y+1 z-1
       			*x=deter(*x-1,dir[0],pbc[0]);
			*y=deter(*y+1,dir[1],pbc[1]);
			*z=deter(*z-1,dir[2],pbc[2]);
                        break;
		case 7 :  // x-1 y+1 z  
       			*x=deter(*x-1,dir[0],pbc[0]);
			*y=deter(*y+1,dir[1],pbc[1]);
                        break;
		case 8 :  // x-1 y+1 z+1
       			*x=deter(*x-1,dir[0],pbc[0]);
			*y=deter(*y+1,dir[1],pbc[1]);
			*z=deter(*z+1,dir[2],pbc[2]);
                        break;
		case 9 :  // x   y-1 z-1
			*y=deter(*y-1,dir[1],pbc[1]);
			*z=deter(*z-1,dir[2],pbc[2]);
                        break;
		case 10 : // x   y-1 z  
       			*y=deter(*y-1,dir[1],pbc[1]);
                        break;
		case 11 : // x   y-1 z+1
       			*y=deter(*y-1,dir[1],pbc[1]);
       			*z=deter(*z+1,dir[2],pbc[2]);
                        break;
		case 12 : // x   y   z-1
       			*z=deter(*z-1,dir[2],pbc[2]);
                        break;
		case 13 : // x   y   z+1
       			*z=deter(*z+1,dir[2],pbc[2]);
                        break;
		case 14 : // x   y+1 z-1
       			*y=deter(*y+1,dir[1],pbc[1]);
			*z=deter(*z-1,dir[2],pbc[2]);
                        break;
		case 15 : // x   y+1 z  
       			*y=deter(*y+1,dir[1],pbc[1]);
                        break;
		case 16 : // x   y+1 z+1
       			*y=deter(*y+1,dir[1],pbc[1]);
			*z=deter(*z+1,dir[2],pbc[2]);
                        break;
		case 17 : // x+1 y-1 z-1
			*x=deter(*x+1,dir[0],pbc[0]);
			*y=deter(*y-1,dir[1],pbc[1]);
			*z=deter(*z-1,dir[2],pbc[2]);
                        break;
		case 18 : // x+1 y-1 z  
       			*x=deter(*x+1,dir[0],pbc[0]);
			*y=deter(*y-1,dir[1],pbc[1]);
                        break;
		case 19 : // x+1 y-1 z+1
       			*x=deter(*x+1,dir[0],pbc[0]);
			*y=deter(*y-1,dir[1],pbc[1]);
			*z=deter(*z+1,dir[2],pbc[2]);
                        break;
		case 20 : // x+1 y   z-1
       			*x=deter(*x+1,dir[0],pbc[0]);
			*z=deter(*z-1,dir[2],pbc[2]);
                        break;
		case 21 : // x+1 y   z  
       			*x=deter(*x+1,dir[0],pbc[0]);
                        break;
		case 22 : // x+1 y   z+1
       			*x=deter(*x+1,dir[0],pbc[0]);
			*z=deter(*z+1,dir[2],pbc[2]);
                        break;
		case 23 : // x+1 y+1 z-1
       			*x=deter(*x+1,dir[0],pbc[0]);
			*y=deter(*y+1,dir[1],pbc[1]);
			*z=deter(*z-1,dir[2],pbc[2]);
                        break;
		case 24 : // x+1 y+1 z  
       			*x=deter(*x+1,dir[0],pbc[0]);
			*y=deter(*y+1,dir[1],pbc[1]);
                        break;
		case 25 ://25x+1 y+1 z+1
       			*x=deter(*x+1,dir[0],pbc[0]);
			*y=deter(*y+1,dir[1],pbc[1]);
			*z=deter(*z+1,dir[2],pbc[2]);
                        break;
		default:
puts("ERROR: wrong random number over 25 for nei26sel function...");
			break;
	}
	return;
}

void nei26sel_nodeter(int i,int *x,int *y,int *z,int dir[],int pbc[]){
	switch(i){
		case 0 :  // x-1 y-1 z-1
			*x=*x-1;
			*y=*y-1;
			*z=*z-1;
			break;
		case 1 :  // x-1 y-1 z
			*x=*x-1;
			*y=*y-1;
			break;
		case 2 :  // x-1 y-1 z+1
			*x=*x-1;
			*y=*y-1;
			*z=*z+1;
                        break;
		case 3 :  // x-1 y   z-1
       			*x=*x-1;
			*z=*z-1;
                        break;
		case 4 :  // x-1 y   z
       			*x=*x-1;
                        break;
		case 5 :  // x-1 y   z+1
       			*x=*x-1;
			*z=*z+1;
                        break;
		case 6 :  // x-1 y+1 z-1
       			*x=*x-1;
			*y=*y+1;
			*z=*z-1;
                        break;
		case 7 :  // x-1 y+1 z  
       			*x=*x-1;
			*y=*y+1;
                        break;
		case 8 :  // x-1 y+1 z+1
       			*x=*x-1;
			*y=*y+1;
			*z=*z+1;
                        break;
		case 9 :  // x   y-1 z-1
			*y=*y-1;
			*z=*z-1;
                        break;
		case 10 : // x   y-1 z  
       			*y=*y-1;
                        break;
		case 11 : // x   y-1 z+1
       			*y=*y-1;
       			*z=*z+1;
                        break;
		case 12 : // x   y   z-1
       			*z=*z-1;
                        break;
		case 13 : // x   y   z+1
       			*z=*z+1;
                        break;
		case 14 : // x   y+1 z-1
       			*y=*y+1;
			*z=*z-1;
                        break;
		case 15 : // x   y+1 z  
       			*y=*y+1;
                        break;
		case 16 : // x   y+1 z+1
       			*y=*y+1;
			*z=*z+1;
                        break;
		case 17 : // x+1 y-1 z-1
			*x=*x+1;
			*y=*y-1;
			*z=*z-1;
                        break;
		case 18 : // x+1 y-1 z  
       			*x=*x+1;
			*y=*y-1;
                        break;
		case 19 : // x+1 y-1 z+1
       			*x=*x+1;
			*y=*y-1;
			*z=*z+1;
                        break;
		case 20 : // x+1 y   z-1
       			*x=*x+1;
			*z=*z-1;
                        break;
		case 21 : // x+1 y   z  
       			*x=*x+1;
                        break;
		case 22 : // x+1 y   z+1
       			*x=*x+1;
			*z=*z+1;
                        break;
		case 23 : // x+1 y+1 z-1
       			*x=*x+1;
			*y=*y+1;
			*z=*z-1;
                        break;
		case 24 : // x+1 y+1 z  
       			*x=*x+1;
			*y=*y+1;
                        break;
		case 25 ://25x+1 y+1 z+1
       			*x=*x+1;
			*y=*y+1;
			*z=*z+1;
                        break;
		default:
puts("ERROR: wrong random number over 25 for nei26sel_nodeter function...");
printf("factor %d x %d y %d z %d\n",i,*x,*y,*z);
			break;
	}
	return;
}
int get_nei26_index(int dx,int dy,int dz) {
	int index;
//    int dx=a - x0;	    int dy=b - y0;	    int dz=c - z0;
// dx, dy, dz must be -1, 0, or 1
	if (abs(dx)>1)
	    if(dx>0)
		dx=-1;
	    else
		dx=1;
	if (abs(dy)>1)
	    if(dy>0)
		dy=-1;
	    else
		dy=1;
	if (abs(dz)>1)
	    if(dz>0)
		dz=-1;
	    else
		dz=1;
    // 3진법 순서 인덱싱 (0~26), 자기 자신 (13)은 제외해야 하므로 -1 보정
	index=(dx+1)*9+(dy+1)*3+(dz+1);
	if (index>13)
	    return index-1;
	else
	    return index;
}


int VN_to_Moore(int VNori){
	switch(VNori){ // +x, -x, +y, -y, +z, -z
		case 0 : // -x
			return 4;
			break;
		case 1 : // +x
			return 21;
			break;
		case 2 : // -y
			return 10;
			break;
		case 3 : // +y
			return 15;
			break;
		case 4 : // -z
			return 12;
			break;
		case 5 : // +z
			return 13;
			break;
		default:
printf("ERROR: wrong case for VN_to_Moore functional: %d\n",VNori);
			return -1;
			break;
	}
}

int deter(int a,int lim,int pbc){// PBC determiner
	if(a<0){
		return (pbc!=0) ? lim+a : 0;
	}else if(a >= lim){
		return(pbc!=0) ? a-lim : lim-1;
	}else
		return a;
}

double misorient(int a,int b){
	return ((double)abs(a-b)<=MAXMIS) ? (double)abs(a-b) : 2*MAXMIS-(double)abs(a-b);
}

double GBE(int delta,double cutoff,double GBEmax){
	return (delta<cutoff) ? GBEmax*(1-log(delta/cutoff))*delta/cutoff : GBEmax;
}

/*double dendrite_dir(int const now[],int dir,double factor){ // return 0.5 for preferred direction
// now = grid[x][y][z] // at this stage only [0] is used
// dir: refer to nei26sel
	int i=0,vector[2][3]={0};
	double class=(double)MAXDG1/NNEI;

	for(i=25;i>0;i--){
		if(now[0]>=class*i)
		    break;
	}
	// equivalent perpendicular directions
//printf("i %d class %f vs. now %d vs. dir %d\n",i,class,now[0],dir);
	if(i==dir) // preferred direction
		return factor;
	num_to_vector(i,vector[0]);
	num_to_vector(dir,vector[1]);
	if(vector[0][0]+vector[1][0]==0\
	 &&vector[0][1]+vector[1][1]==0\
	 &&vector[0][2]+vector[1][2]==0) // reverse direction
		return factor;
	if(inner_product(vector[0],vector[1])==0){ // perpendicular relationship
		if((int)(rand()/(RAND_MAX/2+1))!=0) // preferred direction
			return factor;
		else // not preferred direction
			return 1.0;
	}
	return 1.0;
}

void num_to_vector(int ori,int v[]){
// Moore neighborhood (#26)
	switch(ori){
		case 0 :  // x-1 y-1 z-1
			v[0]=-1;
			v[1]=-1;
			v[2]=-1;
			break;
		case 1 :  // x-1 y-1 z
			v[0]=-1;
			v[1]=-1;
			break;
		case 2 :  // x-1 y-1 z+1
			v[0]=-1;
			v[1]=-1;
			v[2]=+1;
                        break;
		case 3 :  // x-1 y   z-1
			v[0]=-1;
			v[2]=-1;
                        break;
		case 4 :  // x-1 y   z
			v[0]=-1;
                        break;
		case 5 :  // x-1 y   z+1
			v[0]=-1;
			v[2]=+1;
                        break;
		case 6 :  // x-1 y+1 z-1
			v[0]=-1;
			v[1]=+1;
			v[2]=-1;
                        break;
		case 7 :  // x-1 y+1 z  
			v[0]=-1;
			v[1]=+1;
                        break;
		case 8 :  // x-1 y+1 z+1
			v[0]=-1;
			v[1]=+1;
			v[2]=+1;
                        break;
		case 9 :  // x   y-1 z-1
			v[1]=-1;
			v[2]=-1;
                        break;
		case 10 : // x   y-1 z  
			v[1]=-1;
                        break;
		case 11 : // x   y-1 z+1
			v[1]=-1;
			v[2]=+1;
                        break;
		case 12 : // x   y   z-1
			v[2]=-1;
                        break;
		case 13 : // x   y   z+1
			v[2]=+1;
                        break;
		case 14 : // x   y+1 z-1
			v[1]=+1;
			v[2]=-1;
                        break;
		case 15 : // x   y+1 z  
			v[1]=+1;
                        break;
		case 16 : // x   y+1 z+1
			v[1]=+1;
			v[2]=+1;
                        break;
		case 17 : // x+1 y-1 z-1
			v[0]=+1;
			v[1]=-1;
			v[2]=-1;
                        break;
		case 18 : // x+1 y-1 z  
			v[0]=+1;
			v[1]=-1;
                        break;
		case 19 : // x+1 y-1 z+1
			v[0]=+1;
			v[1]=-1;
			v[2]=+1;
                        break;
		case 20 : // x+1 y   z-1
			v[0]=+1;
			v[2]=-1;
                        break;
		case 21 : // x+1 y   z  
			v[0]=+1;
                        break;
		case 22 : // x+1 y   z+1
			v[0]=+1;
			v[2]=+1;
                        break;
		case 23 : // x+1 y+1 z-1
			v[0]=+1;
			v[1]=+1;
			v[2]=-1;
                        break;
		case 24 : // x+1 y+1 z  
			v[0]=+1;
			v[1]=+1;
                        break;
		case 25 ://25x+1 y+1 z+1
			v[0]=+1;
			v[1]=+1;
			v[2]=+1;
                        break;
		default:
puts("ERROR: wrong random number over 25 for num_to_vector function...");
			break;
	}
	return;
}

double IFE_aniso(double gamma,double factor,int predir){ //maxdir==1 for prefered direction
	return (predir==0)? gamma : gamma*factor;
}*/

int nucl_cond(int**** grid,int x,int y,int z,int dir[],int pbc[],int situation){
	int i,l,m,n;

	for(i=25;i>=0;i--){
	    if(mpi_surf_check(i,x,y,z,grid)>=0)
		return -mpi_surf_check(i,x,y,z,grid)-1;
	    if(surf_check(i,x,y,z,dir,pbc)>=0) // surface condition
		return -surf_check(i,x,y,z,dir,pbc)-1;
	}
	if(situation==0){ // solidification
	    for(i=25;i>=0;i--){
		l=x;
		m=y;
		n=z;
		nei26sel(i,&l,&m,&n,dir,pbc);
		if(grid[l][m][n][0]==PTCLORI) // PARTICLE exist
		    return PTCLORI;
	    }
	    for(i=25;i>=0;i--){
		l=x;
		m=y;
		n=z;
		nei26sel(i,&l,&m,&n,dir,pbc);
		if(grid[l][m][n][0]!=LIQORI) // solid exist
		    return grid[l][m][n][0];
	    }
	    return 0; // no solid around
	}
	puts("ERROR in NUCL_CONDITION function...");
	return -9999999;
}

int mpi_surf_check(int i,int x,int y,int z,int**** grid){
	if(n_size==1)
	    return -1;
	if(i>=0 && i<=8){
	    if(grid[x-1][y][z][0]==MPISURF)
		return 0;
	    else
		return -1;
	}else if(i>=17){
	    if(grid[x+1][y][z][0]==MPISURF)
		return 1;
	    else
		return -1;
	}else
	    return -1;
printf("Surface Pattern Error...: i is %d\n",i);
return -999999; //ERROR CODE
}

int surf_check(int i,int x,int y,int z,int dir[],int pbc[]){ //return N>=0 when surface
	switch(i){
		case 0 :  // x-1 y-1 z-1
			if(x==0&&pbc[0]==0)
				return 0;
			else if(y==0&&pbc[1]==0)
				return 2;
			else if(z==0&&pbc[2]==0)
				return 4;
			break;
		case 1 :  // x-1 y-1 z
			if(x==0&&pbc[0]==0)
				return 0;
			else if(y==0&&pbc[1]==0)
				return 2;
			break;
		case 2 :  // x-1 y-1 z+1
			if(x==0&&pbc[0]==0)
				return 0;
			else if(y==0&&pbc[1]==0)
				return 2;
			else if(z==dir[2]-1&&pbc[2]==0)
				return 5;
                        break;
		case 3 :  // x-1 y   z-1
			if(x==0&&pbc[0]==0)
				return 0;
			else if(z==0&&pbc[2]==0)
				return 4;
			break;
		case 4 :  // x-1 y   z
			if(x==0&&pbc[0]==0)
				return 0;
                        break;
		case 5 :  // x-1 y   z+1
			if(x==0&&pbc[0]==0)
				return 0;
			else if(z==dir[2]-1&&pbc[2]==0)
				return 5;
                        break;
		case 6 :  // x-1 y+1 z-1
			if(x==0&&pbc[0]==0)
				return 0;
			else if(y==dir[1]-1&&pbc[1]==0)
				return 3;
			else if(z==0&&pbc[2]==0)
				return 4;
                        break;
		case 7 :  // x-1 y+1 z  
			if(x==0&&pbc[0]==0)
				return 0;
			else if(y==dir[1]-1&&pbc[1]==0)
				return 3;
                        break;
		case 8 :  // x-1 y+1 z+1
			if(x==0&&pbc[0]==0)
				return 0;
			else if(y==dir[1]-1&&pbc[1]==0)
				return 3;
			else if(z==dir[2]-1&&pbc[2]==0)
				return 5;
                        break;
		case 9 :  // x   y-1 z-1
			if(y==0&&pbc[1]==0)
				return 2;
			else if(z==0&&pbc[2]==0)
				return 4;
                        break;
		case 10 : // x   y-1 z  
			if(y==0&&pbc[1]==0)
				return 2;
                        break;
		case 11 : // x   y-1 z+1
			if(y==0&&pbc[1]==0)
				return 2;
			else if(z==dir[2]-1&&pbc[2]==0)
				return 5;
                        break;
		case 12 : // x   y   z-1
			if(z==0&&pbc[2]==0)
				return 4;
                        break;
		case 13 : // x   y   z+1
			if(z==dir[2]-1&&pbc[2]==0)
				return 5;
                        break;
		case 14 : // x   y+1 z-1
			if(y==dir[1]-1&&pbc[1]==0)
				return 3;
			else if(z==0&&pbc[2]==0)
				return 4;
                        break;
		case 15 : // x   y+1 z  
			if(y==dir[1]-1&&pbc[1]==0)
				return 3;
                        break;
		case 16 : // x   y+1 z+1
			if(y==dir[1]-1&&pbc[1]==0)
				return 3;
			else if(z==dir[2]-1&&pbc[2]==0)
				return 5;
                        break;
		case 17 : // x+1 y-1 z-1
			if(x==dir[0]-1&&pbc[0]==0)
				return 1;
			else if(y==0&&pbc[1]==0)
				return 2;
			else if(z==0&&pbc[2]==0)
				return 4;
                        break;
		case 18 : // x+1 y-1 z  
			if(x==dir[0]-1&&pbc[0]==0)
				return 1;
			else if(y==0&&pbc[1]==0)
				return 2;
                        break;
		case 19 : // x+1 y-1 z+1
			if(x==dir[0]-1&&pbc[0]==0)
				return 1;
			else if(y==0&&pbc[1]==0)
				return 2;
			else if(z==dir[2]-1&&pbc[2]==0)
				return 5;
                        break;
		case 20 : // x+1 y   z-1
			if(x==dir[0]-1&&pbc[0]==0)
				return 1;
			else if(z==0&&pbc[2]==0)
				return 4;
                        break;
		case 21 : // x+1 y   z  
			if(x==dir[0]-1&&pbc[0]==0)
				return 1;
                        break;
		case 22 : // x+1 y   z+1
			if(x==dir[0]-1&&pbc[0]==0)
				return 1;
			else if(z==dir[2]-1&&pbc[2]==0)
				return 5;
                        break;
		case 23 : // x+1 y+1 z-1
			if(x==dir[0]-1&&pbc[0]==0)
				return 1;
			else if(y==dir[1]-1&&pbc[1]==0)
				return 3;
			else if(z==0&&pbc[2]==0)
				return 4;
                        break;
		case 24 : // x+1 y+1 z  
			if(x==dir[0]-1&&pbc[0]==0)
				return 1;
			else if(y==dir[1]-1&&pbc[1]==0)
				return 3;
                        break;
		default ://25 x+1 y+1 z+1
			if(x==dir[0]-1&&pbc[0]==0)
				return 1;
			else if(y==dir[1]-1&&pbc[1]==0)
				return 3;
			else if(z==dir[2]-1&&pbc[2]==0)
				return 5;
                        break;
	}
	return -1;
}

void pow2grain(int ori[]){
	ori[0]-=POWORI;
	ori[1]=1;
	ori[2]=1;
	return;
}

int latent_heat_check(double*** tgrid,int dir[]){
	int i,j,k;
	for(i=dir[0]-1;i>=0;i--){
	    for(j=dir[1]-1;j>=0;j--){
		for(k=dir[2]-2;k>=0;k--){
		    if(tgrid[i][j][k]!=0)
			return 1;
		}
	    }
	}
	return 0;
}

double average_temp(double*** temp,int x,int y,int z,int inter,int dir[],int pbc[]){
	int i,j,k,o;
	double sum=0,n=1;
	
	sum=temp[x][y][z];
/*	for(o=5;o>=0;o--){ // for FDM average, VN neighborhood is only considered
	// Assuming that their influence is half of original positiona
	    i=x;
	    j=y;
	    k=z;
	    nei6sel(o,&i,&j,&k,dir,pbc);
	    if(i!=x || j!=y || k!=z){
		sum=0.5*temp[i][j][k];
		n+=0.5;
	    }
	}
*/
	o=inter/2;
	for(i=deter(x-o,dir[0],pbc[0]);i<=deter(x+o,dir[0],pbc[0]);i++){
	    for(j=deter(y-o,dir[1],pbc[1]);j<=deter(y+o,dir[1],pbc[1]);j++){
		for(k=deter(z-o,dir[2],pbc[2]);k<=deter(z+o,dir[2],pbc[2]);k++){
		    sum+=temp[i][j][k];
		    n++;
		}
	    }
	}
//	*/
	return sum/n;
}

double recover_temp(double*** temp,int x,int y,int z,int del,int pts[]){
	int i;
	double dist[8],sum=0;
//printf("%d %d %d vs. %d %d %d //",x,y,z,pts[0]-del,pts[1]-del,pts[2]-del);
// temp-del < x < temp
	dist[0]=distance(x,y,z,pts[0]-del,pts[1]-del,pts[2]-del); // - - -
	dist[1]=distance(x,y,z,pts[0],pts[1]-del,pts[2]-del); // + - -
	dist[2]=distance(x,y,z,pts[0]-del,pts[1],pts[2]-del); // - + -
	dist[3]=distance(x,y,z,pts[0]-del,pts[1]-del,pts[2]); // - - +
	dist[4]=distance(x,y,z,pts[0],pts[1],pts[2]-del); // + + -
	dist[5]=distance(x,y,z,pts[0],pts[1]-del,pts[2]); // + - +
	dist[6]=distance(x,y,z,pts[0]-del,pts[1],pts[2]); // - + +
	dist[7]=distance(x,y,z,pts[0],pts[1],pts[2]); // + + +
//printf("dist %.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",dist[0],dist[1],dist[2],dist[3],dist[4],dist[5],dist[6],dist[7]);
	for(i=7;i>=0;i--){
	    if(dist[i]>=del)
		dist[i]=0;
	    else
		sum+=dist[i];
	}
	return (temp[pts[0]-del][pts[1]-del][pts[2]-del]*dist[0]\
		+temp[pts[0]][pts[1]-del][pts[2]-del]*dist[1]\
		+temp[pts[0]-del][pts[1]][pts[2]-del]*dist[2]\
		+temp[pts[0]-del][pts[1]-del][pts[2]]*dist[3]\
		+temp[pts[0]][pts[1]][pts[2]-del]*dist[4]\
		+temp[pts[0]][pts[1]-del][pts[2]]*dist[5]\
		+temp[pts[0]-del][pts[1]][pts[2]]*dist[6]\
		+temp[pts[0]][pts[1]][pts[2]]*dist[7])/sum;
}

double distance(int p11,int p12,int p13,int p21,int p22,int p23){
	return sqrt((double)((p11-p21)*(p11-p21)+(p12-p22)*(p12-p22)+(p13-p23)*(p13-p23)));
}

int get_file(char fnam[],int mode){// Return existing file name
	FILE* inp=NULL;
	while(1){
		if(mode==0){
			printf("\nSample file name? (Exit = input 'q')\n >> ");
			scanf("%s",fnam);
		}else if(mode==-1){
			printf("\nNew sample file name?\n >> ");
			scanf("%s",fnam);
			return 1;
		}else if(mode==-2){
			printf("\nNew T map file name? (*.tdat)\n >> ");
			scanf("%s",fnam);
			return 1;
		}else if(mode==-3){
			printf("\nT map file name? (*.tdat)\n >> ");
			scanf("%s",fnam);
		}
/*		else if(mode==-4){
			printf("\nInput file name?\n >> ");
			scanf("%s",fnam);
		}*/
		if(strcmp(fnam,"q")==0){
			printf("\nExit program\n");
			return 0;
		}

		inp=fopen(fnam,"r");
		if(inp==NULL){
			printf("\nNo file named \"%s\": Please input other(exit=input \"q\")\n",fnam);
			fflush(stdin);
			mode=0;
		}else if(mode==0){
			printf("\nStart to read data\n");
			fclose(inp);
			return 1;
		}else if(mode==-3){
			fclose(inp);
			if(strstr(fnam,"tdat")!=NULL){
				printf("\nStart to read T map data\n");
				return 1;
			}
			fflush(stdin);
			printf("File format must be *.tdat ... \n");
		}
		/*else if(mode=-4){
			printf("\nStart to read simulation condition from %s file...\n",fnam);
			fclose(inp);
			return 1;
		}*/
	}
}

int**** read_info(int**** grid,int dir[],char file[]){// Read grain info.
	int i,j,k,x[2],y[2],z[2],now,prev[3]={0};
	FILE* input=NULL;

	input=fopen(file,"r");
	printf("\nLoad input file \"%s\"...\n",file);
	fscanf(input,"%d %d %d\n",&dir[0],&dir[1],&dir[2]);
	printf(" ( Site = %d x %d x %d ) ",dir[0],dir[1],dir[2]);
	grid=alloc4d(dir[0],dir[1],dir[2],grid); //grid allocation

	fscanf(input,"%d %d %d %d\n",&x[0],&y[0],&z[0],&prev[0]); //,&prev[1],&prev[2] N/A yet
	if(x[0]!=0||y[0]!=0||z[0]!=0){
		printf("Error: wrong file format...\n");
		grid=free4d(dir[0],dir[1],dir[2],grid);
		return NULL;
	}
	if(prev[0]==LIQORI){ // 0 0 0 is for liquid
		prev[1]=LIQORI;
		prev[2]=LIQORI;
	}else if(prev[0]>=POWORI){ // Powder
		prev[1]=POWORI;
		prev[2]=POWORI;
	}else{
		prev[1]=1;
		prev[2]=1;
	}
	grid[0][0][0][0]=prev[0];
	grid[0][0][0][1]=prev[1];
	grid[0][0][0][2]=prev[2];

	while(feof(input)==0){
	    fscanf(input,"%d %d %d %d\n",&x[1],&y[1],&z[1],&now);
	    if(x[0]==x[1] && y[0]==y[1]){
		    for(k=z[0];k<z[1];k++){
			grid[x[0]][y[0]][k][0]=prev[0];
			grid[x[0]][y[0]][k][1]=prev[1];
			grid[x[0]][y[0]][k][2]=prev[2];
		    }
	    }else{
		for(k=z[0];k<dir[2];k++){
		    grid[x[0]][y[0]][k][0]=prev[0];
		    grid[x[0]][y[0]][k][1]=prev[1];
		    grid[x[0]][y[0]][k][2]=prev[2];
		}
	    }
	    if(now==0){ // 0 0 0 is for liquid
		prev[0]=-1;
		prev[1]=-1;
		prev[2]=-1;
	    }else{
		prev[0]=now;
		prev[1]=1;
		prev[2]=1;
	    }
	    x[0]=x[1];
	    y[0]=y[1];
	    z[0]=z[1];
	}
	for(k=z[0];k<dir[2];k++){
	    grid[x[0]][y[0]][k][0]=prev[0];
	    grid[x[0]][y[0]][k][1]=prev[1];
	    grid[x[0]][y[0]][k][2]=prev[2];
	}

	for(i=dir[0]-1;i>=0;i--){
	    for(j=dir[1]-1;j>=0;j--){
		for(k=dir[2]-1;k>=0;k--){
		    if(grid[i][j][k][0]>=POWORI){
//			neworient_powder(grid[i][j][k]);
			grid[i][j][k][1]=POWORI;
			grid[i][j][k][2]=POWORI;
		    }
		}
	    }
	}

	printf("successfully loaded...!!!\n",file);
	fclose(input);

	return grid;
}
/*
int**** read_infoPrev(int**** grid,int dir[],char file[]){// Read grain info.
	// old version of MCS.dat
	int i=0,j=0,k=0,maxi=0,maxj=0,maxk=0,deg[3]={0};
	FILE* input=NULL;

	input=fopen(file,"r");
	printf("\nLoad input file \"%s\"...\n",file);
if(dir[0]==0){
	while(feof(input)==0){
		fscanf(input,"%d\t%d\t%d\t%d\t%d\t%d\n",&i,&j,&k,&deg[0],&deg[1],&deg[2]);
		if(i>maxi)
			maxi=i;
		if(j>maxj)
			maxj=j;
		if(k>maxk)
			maxk=k;
	}
	rewind(input);
	dir[0]=maxi+1;
	dir[1]=maxj+1;
	dir[2]=maxk+1;
	printf(" ( Site = %d x %d x %d ) ",dir[0],dir[1],dir[2]);
}
	grid=alloc4d(dir[0],dir[1],dir[2],grid); //grid allocation

	while(feof(input)==0){
		fscanf(input,"%d\t%d\t%d\t%d\t%d\t%d\n",&i,&j,&k,&deg[0],&deg[1],&deg[2]);
		if(deg[0]==0){ // 0 0 0 is for liquid
			grid[i][j][k][0]=-1;
			grid[i][j][k][1]=-1;
			grid[i][j][k][2]=-1;
		}else{
			grid[i][j][k][0]=deg[0];
			grid[i][j][k][1]=deg[1];
			grid[i][j][k][2]=deg[2];
		}
	}
	printf("successfully loaded...!!!\n",file);
	fclose(input);

	return grid;
}
// */

double*** read_Tmap(double*** temp,int dir[],char file[],double range[],int mode,double Tmelt){
// mode 0 for drawing 1 to apply for sim. -1 for drawing after sim.
// Load temperature distribution info.
	int i=0,j=0,k=0,maxi=0,maxj=0,maxk=0,forcematch=0;
	double maxt=0,mint=1E6,t;
	FILE* input=NULL;

	input=fopen(file,"r");
	if(input==NULL){
		printf("No file named %s...\n",file);
		return NULL;
	}
	printf("\nLoad input file \"%s\"...\n",file);
	if(mode>=0){ // load sample for drawing or simulation
	    while(feof(input)==0){
		fscanf(input,"%d\t%d\t%d\t%lf\n",&i,&j,&k,&t);
		if(i>maxi)
			maxi=i;
		if(j>maxj)
			maxj=j;
		if(k>maxk)
			maxk=k;
		if(t>maxt)
			maxt=t;
		if(t<mint)
			mint=t;
	    }
	    rewind(input);
	    if(mode==0){
		dir[0]=maxi+1;
		dir[1]=maxj+1;
		dir[2]=maxk+1;
		range[0]=mint;
		range[1]=maxt;
		printf(" ( Site = %d x %d x %d ) ",dir[0],dir[1],dir[2]);
	    }
	    if(mode==1){
		printf(" ( Site = %d x %d x %d ) ",maxi+1,maxj+1,maxk+1);
		if(dir[0]!=maxi+1||dir[1]!=maxj+1||dir[2]!=maxk+1){
			printf(" @@@ WARNING: number of voxels are different...\n     Forced-matching activated...\n");
			forcematch=1;
		}
	    }
	}else{
// load sample after simulation for plot graph
	    rewind(input);
	    range[0]=Tmelt*KAUZMANN;
	    range[1]=Tmelt;
	}
	temp=alloc3d(dir[0],dir[1],dir[2],temp);

	while(feof(input)==0){
		fscanf(input,"%d\t%d\t%d\t%lf\n",&i,&j,&k,&t);
		if(forcematch==0)
			temp[i][j][k]=t;
		else{
			if(i<dir[0]&&j<dir[1]&&k<dir[2])
				temp[i][j][k]=t;
		}
	}
	if(forcematch!=0){ //forcematching mode
	// larger information amount -> crop
	// smaller information amount -> expand from its last value
	    if(maxi+1<dir[0]){
		for(i=dir[0]-1;i>=maxi;i--){
		    for(j=dir[1]-1;j>=0;j--){
			for(k=dir[2]-1;k>=0;k--){
			    temp[i][j][k]=temp[dir[0]-1][j][k];
		}   }	}
	    }
	    if(maxj+1<dir[1]){
		for(i=dir[0]-1;i>=0;i--){
		    for(j=dir[1]-1;j>=maxj;j--){
			for(k=dir[2]-1;k>=0;k--){
			    temp[i][j][k]=temp[i][dir[1]-1][k];
		}   }	}
	    }
	    if(maxk+1<dir[2]){
		for(i=dir[0]-1;i>=0;i--){
		    for(j=dir[1]-1;j>=0;j--){
			for(k=dir[2]-1;k>=maxk;k--){
			    temp[i][j][k]=temp[i][j][dir[2]-1];
		}   }	}
	    }
	}
	printf("Successfully loaded...!!!\n",file);
	fclose(input);

	return temp;
}

int**** alloc4d(int x,int y,int z,int**** grid){// 3DxN_INFO Memory allocation
	int i,j,k;
	int* data;
// for continuous allocation
	data=(int*)calloc(x*y*z*NINFO,sizeof(int));
	grid=(int****)malloc(x*sizeof(int***));
	for(i=0;i<x;i++){
		grid[i]=(int***)malloc(y*sizeof(int**));
		for(j=0;j<y;j++){
			grid[i][j]=(int**)malloc(z*sizeof(int*));
			for(k=0;k<z;k++)
			    grid[i][j][k]=data+(((i*y*z)+(j*z)+k)*NINFO);
//			    grid[i][j][k]=(int*)calloc(NINFO,sizeof(int));
		}	
	}
	return grid;
}

int**** free4d(int x,int y,int z,int**** grid){// 3Dx3 memory de-allocation
	int i,j,k;

	free(grid[0][0][0]); //pointing contiunous data allocation
	for(i=0;i<x;i++){ //Memory deallocation
		for(j=0;j<y;j++){
/*			for(k=0;k<z;k++)
				free(grid[i][j][k]);*/
			free(grid[i][j]);
		}free(grid[i]);
	}free(grid);
	grid=NULL;
	return grid;
}

double*** alloc3d(int x,int y,int z,double*** grid){// 3DxN Memory allocation
	int i,j,k;
	double* data;
// for continuous allocation
	data=(double*)calloc(x*y*z,sizeof(double));
	grid=(double***)malloc(x*sizeof(double**));
	for(i=0;i<x;i++){
		grid[i]=(double**)malloc(y*sizeof(double*));
		for(j=0;j<y;j++)
			grid[i][j]=data+((i*y*z)+(j*z));
//			grid[i][j][k]=(int*)calloc(NINFO,sizeof(int));
	}
	return grid;
}

double*** free3d(int x,int y,int z,double*** grid){// 3D memory de-allocation
	int i,j,k;

	free(grid[0][0]); //pointing contiunous data allocation
	for(i=0;i<x;i++){ //Memory deallocation
		free(grid[i]);
	}free(grid);
	grid=NULL;
	return grid;
}

void neworient(int grid[]){	// new Euler angle info. : 1~90 / 1 ~ 90 / 1 ~ 90
	grid[0]=(int)(rand()/(RAND_MAX/MAXDG1+1))+1;
//	grid[1]=(int)(rand()/(RAND_MAX/MAXDG2+1))+1;
//	grid[2]=(int)(rand()/(RAND_MAX/MAXDG3+1))+1;
	grid[1]=1;
	grid[2]=1;
	return;
}

void neworient_powder(int grid[]){
	grid[0]=(int)(rand()/(RAND_MAX/MAXDG1+1))+1+POWORI;
	grid[1]=POWORI;
	grid[2]=POWORI;
	return;
}

void neworient_melting(int grid[]){
	grid[0]=LIQORI;
	grid[1]=LIQORI;
	grid[2]=LIQORI;
	return;
}

void datout(int dir[],int**** grid,int mcs,int out[],double const realtime[],double rt){// Write result.dat file
	int i;
	double avg=0;
	char info[50]="\0";
	char* temp;
	FILE* list=NULL;

	avg=grainsize_ff(grid,dir,0);
	sprintf(info,"%-6.2lf",avg);
	if(mcs==0){
		list=fopen("result.dat","w");
		fprintf(list,"$ MCS\t GrainSize(Voxels)");
		if(rt!=0){ // Real-time fitting
			fprintf(list,"     time(s)\t  GrainSize(m)");
			if(realtime[2]!=0) // solid volume fraction
				fprintf(list,"    f_liquid\n");
			else
				fprintf(list,"\n");
			fprintf(list,"%4d\t %s\t\t%12.3E    %14.3f",mcs,info,realtime[0],realtime[3]*avg);
		}else
			fprintf(list,"\n%4d\t %s",mcs,info);
	}else{
		list=fopen("result.dat","a");
		if(rt!=0) // Real-time fitting
			fprintf(list,"%4d\t %s\t\t%12.3E    %14.3f",mcs,info,realtime[0],realtime[3]*avg);
		else
			fprintf(list,"%4d\t %s",mcs,info);
	}
	if(realtime[2]!=0)
		fprintf(list,"       %.6f",lfrac(grid,dir));
	fprintf(list,"\n");
	fclose(list);

	return;
}

double lfrac(int**** grid,int dir[]){// Count liquid fraction
	int i,j,k,l=0;
	for(i=dir[0]-1;i>=0;i--){
		for(j=dir[1]-1;j>=0;j--){
			for(k=dir[2]-1;k>=0;k--){
				if(grid[i][j][k][0]==LIQORI)
					l++;
	}	}	}
	return (double)l/(double)(dir[0]*dir[1]*dir[2]);
}

void tproout(int dir[],double*** temp,int mcs){// Make MCS.tdat file
	char fnam[30];
	int i,j,k;
	FILE* data=NULL;

	if(mcs>=0){
		sprintf(fnam,"MCS%05d.tdat",mcs);
		data=fopen(fnam,"w");
		for(i=0;i<dir[0];i++){
		for(j=0;j<dir[1];j++){
		for(k=0;k<dir[2];k++)
			fprintf(data,"%d\t%d\t%d\t%f\n",i,j,k,temp[i][j][k]);
		}}
		printf("    \"%s\" T file created...\n",fnam);
		fclose(data);
	}else{
		i=get_file(fnam,-2);
		data=fopen(fnam,"w");
		for(i=0;i<dir[0];i++){
		for(j=0;j<dir[1];j++){
		for(k=0;k<dir[2];k++)
			fprintf(data,"%d\t%d\t%d\t%f\n",i,j,k,temp[i][j][k]);
		}}
		printf("    \"%s\" file created...\n",fnam);
		fclose(data);
	}

	return;
}
double grainsize(int**** const grid,int dir[],char info[],int plane,int cor){
//prev: Average Grain Intercept Method
	int i,j,k,now[3]={0},prev[3]={0},bc[3]={0},gb,max;
	double sy,smax=0,smin,dev=0;
	double* tavg=NULL;

	switch(cor){
		case 0: // x : x=N yz plane, line along z direction
	max=dir[1];
	smin=dir[1];
	tavg=(double*)calloc(max,sizeof(double));
	for(j=max-1;j>=0;j--){
		gb=0;
		voxcpy(prev,grid[plane][j][dir[2]-1]);
		voxcpy(bc,prev); // to consider PBC
		for(k=dir[2]-2;k>=0;k--){
			voxcpy(now,grid[plane][j][k]);
			if(now[0]!=prev[0] || now[1]!=prev[1]){
			//if(now[0]!=prev[0] || now[1]!=prev[1] || now[2]!=prev[2]){
				voxcpy(prev,now);
				gb++;
			}
		}
		if(now[0]==bc[0])
			gb--;
		if(gb>0)
			sy=max/(double)gb;
		else
			sy=max;
		tavg[j]=sy;
		if(smax<sy)
			smax=sy;
		else if(smin>sy)
			smin=sy;
	}
			break;
		case 1: // y : y=N xz plane, scan along x direction
	max=dir[2];
	smin=dir[2];
	tavg=(double*)calloc(max,sizeof(double));
	for(j=max-1;j>=0;j--){
		gb=0;
		voxcpy(prev,grid[dir[0]-1][plane][j]);
		voxcpy(bc,prev); // to consider PBC
		for(k=dir[0]-2;k>=0;k--){
			voxcpy(now,grid[k][plane][j]);
			if(now[0]!=prev[0] || now[1]!=prev[1]){
			//if(now[0]!=prev[0] || now[1]!=prev[1] || now[2]!=prev[2]){
				voxcpy(prev,now);
				gb++;
			}
		}
		if(now[0]==bc[0])
			gb--;
		if(gb>0)
			sy=max/(double)gb;
		else
			sy=max;
		tavg[j]=sy;
		if(smax<sy)
			smax=sy;
		else if(smin>sy)
			smin=sy;
	}
		case 2: // z : z=N xy plane, scan along x direction
	max=dir[1];
	smin=dir[1];
	tavg=(double*)calloc(max,sizeof(double));
	for(j=max-1;j>=0;j--){
		gb=0;
		voxcpy(prev,grid[dir[0]-1][j][plane]);
		voxcpy(bc,prev); // to consider PBC
		for(k=dir[0]-2;k>=0;k--){
			voxcpy(now,grid[k][j][plane]);
			if(now[0]!=prev[0] || now[1]!=prev[1]){
			//if(now[0]!=prev[0] || now[1]!=prev[1] || now[2]!=prev[2]){
				voxcpy(prev,now);
				gb++;
			}
		}
		if(now[0]==bc[0])
			gb--;
		if(gb>0)
			sy=max/(double)gb;
		else
			sy=max;
		tavg[j]=sy;
		if(smax<sy)
			smax=sy;
		else if(smin>sy)
			smin=sy;
	}
			break;
		default:
			printf("Error @ Old grain size module...\n\n");
			break;
	}

	sy=0;
	for(i=max-1;i>=0;i--)
		sy=sy+tavg[i];
	sy=sy/(double)max;
	for(i=max-1;i>=0;i--)
		dev=dev+(sy-tavg[i])*(sy-tavg[i]);
	dev=sqrt(dev/((double)(max*max)));
	sprintf(info,"%-6.1lf  %-6.2lf ±%4.2lf  %-6.1lf",smin,sy,dev,smax);
	free(tavg);

	return sy;
}

void fileout(int dir[],int**** grid,int mcs){// Make MCS.dat file
	char fnam[30];
	int i,j,k,now;
	FILE* data=NULL;

	if(mcs>=0)
		sprintf(fnam,"MCS%05d.dat",mcs);
	else
		i=get_file(fnam,-1);
	printf("    \"%s\" file was created...\n",fnam);

	data=fopen(fnam,"w");
// new version of MCS.dat file
	fprintf(data,"%d %d %d\n",dir[0],dir[1],dir[2]); // xmax, ymax, zmax

// POWDER check
/* 	for(i=dir[0]-1;i>=0;i--){
	    for(j=dir[1]-1;j>=0;j--){
 		for(k=dir[2]-1;k>=0;k--){
		    if(grid[i][j][k][1]==POWORI){
			now=grid[i][j][k][0];
			grid[i][j][k][0]=grid[i][j][k][2];
			grid[i][j][k][2]=now;
		    }
		}
	    }
	}*/

	for(i=0;i<dir[0];i++){
	    now=grid[i][0][0][0];
	    fprintf(data,"%d %d %d %d\n",i,0,0,now);//,grid[i][j][k][1],grid[i][j][k][2]); 2nd, 3rd euler angle N/A yet
	    for(j=0;j<dir[1];j++){
		if(j!=0){
		    now=grid[i][j][0][0];
		    fprintf(data,"%d %d %d %d\n",i,j,0,now);//,grid[i][j][k][1],grid[i][j][k][2]); 2nd, 3rd euler angle N/A yet
		}
		for(k=0;k<dir[2];k++){
			if(now!=grid[i][j][k][0]){
				now=grid[i][j][k][0];
				fprintf(data,"%d %d %d %d\n",i,j,k,now);
			}
		}
	}   }
	fclose(data);

/* 	for(i=dir[0]-1;i>=0;i--){
	    for(j=dir[1]-1;j>=0;j--){
 		for(k=dir[2]-1;k>=0;k--){
		    if(grid[i][j][k][1]==POWORI)
			neworient_powder(grid[i][j][k]);
		}
	    }
	}
*/
	return;
}

/*
void fileoutPrev(int dir[],int**** grid,int mcs){// Make MCS.dat file
	char fnam[30];
	int i,j,k;
	FILE* data=NULL;

	if(mcs>=0){
		sprintf(fnam,"MCS%05d.dat",mcs);
		printf("    \"%s\" file was created...\n",fnam);
		data=fopen(fnam,"w");
		for(i=0;i<dir[0];i++){
		for(j=0;j<dir[1];j++){
			for(k=0;k<dir[2];k++){
				if(grid[i][j][k][0]<=0) // liquid is 0 0 0
					fprintf(data,"%d\t%d\t%d\t%d\t%d\t%d\n",i,j,k,0,0,0);
				else
					fprintf(data,"%d\t%d\t%d\t%d\t%d\t%d\n",i,j,k,grid[i][j][k][0],grid[i][j][k][1],grid[i][j][k][2]);
			}
		}}
		fclose(data);
	}else{
		i=get_file(fnam,-1);
		printf("    \"%s\" file was created...\n\n",fnam);
		data=fopen(fnam,"w");
		for(i=0;i<dir[0];i++){
		for(j=0;j<dir[1];j++){
			for(k=0;k<dir[2];k++){
				if(grid[i][j][k][0]<=0) // liquid is 0 0 0
					fprintf(data,"%d\t%d\t%d\t%d\t%d\t%d\n",i,j,k,0,0,0);
				else
					fprintf(data,"%d\t%d\t%d\t%d\t%d\t%d\n",i,j,k,grid[i][j][k][0],grid[i][j][k][1],grid[i][j][k][2]);
			}
		}}
		fclose(data);
	}
	return;
}
// */

void QS(double arr[], int left, int right){//quick sort
	double pivot=arr[left];
	double temp;
	int low=left+1;
	int high=right;

    if (left<=right){
    	while (low <= high){
			while (low<=right && pivot>=arr[low])
				low++;
			while (high>=(left+1) && pivot<=arr[high])
				high--;
			if (low<=high){
				temp=arr[low];
				arr[low]=arr[high];
				arr[high]=temp;
        	}
    	}
    	temp=arr[left];
		arr[left]=arr[high];
		arr[high]=temp;
		QS(arr,left,high-1);
		QS(arr,high+1,right);
    }
    return;
}

void tok(char* line,char* word){
	int j;
	for(j=0;word[j]!='\0';j++) //init.
		word[j]='\0';
	while(line[0]!='\0'){//tok
		if(line[0]!=' ' && line[0]!='\t'){
			sprintf(word,"%s%c",word,line[0]);
			for(j=1;line[j]!='\0';j++)
				line[j-1]=line[j];
			line[j-1]=line[j];
		}else{
			while(line[0]==' '||line[0]=='\t'){ //delete isspaces
				for(j=1;line[j]!='\0';j++)
					line[j-1]=line[j];
				line[j-1]=line[j];
			}
			if(word[0]!='\0')
				break;
	}	}
	return;
}

double grainsize_ff(int**** const grid,int dir[],int mode){ // Grain distribution by flood-fill algorithm
	int i,j,k,l=0,gnum=0,num[2]={0},sum[2]={0},lim[6]={0},cnt,len[3],pbc[3];
	double aratio[2],asum=0,ratio=0,vavg[3],ravg[2],dev[2]={0},size;
	int*** checker=NULL;
	char term[10]={'\0'};
// size 0=number of sites in x dir // 1=in y dir // 2=total MCS // 3=mapout

	checker=(int***)malloc(dir[0]*sizeof(int**));
	for(i=0;i<dir[0];i++){
		checker[i]=(int**)calloc(dir[1],sizeof(int*));
		for(j=0;j<dir[1];j++)
			checker[i][j]=(int*)calloc(dir[2],sizeof(int));
	}

    if(mode==0){
	for(k=0;k<dir[2];k++){
	    for(j=0;j<dir[1];j++){
		for(i=0;i<dir[0];i++){
			if(grid[i][j][k][0]!=LIQORI && checker[i][j][k]==0){ // grains not yet measured, no liquid
				gnum++;
				cnt=0;
				DFS_array_simple(dir,i,j,k,grid[i][j][k][0],grid,checker,&cnt,pbc);
				l+=cnt;
			}
		}
	    }
	}

	if(l!=0) // No full-liquid sample
//	ratio=pow(0.75*(double)l/((double)gnum*M_PI),0.3333333)*2.0; // avg. grain size (diameter of sphere)
		ratio=pow((double)l/(double)gnum,0.3333333); // avg. grain size (edge length of cube)
	// else then ratio=0

    }else{ // from Auxiliary module
	printf(" ### NOTE: PBC is not considered ###\n");
/*	while(j==0){
		j=1;
		printf("\nPBC? (XYZ, on=1 / off=0)\nex) 111: PBC in every direction\n    000: no PBC in any direction\n");
		scanf("%s",&temp);
		for(i=0;i<3;i++){
			pbc[i]=temp[i]-'0';
			if(pbc[i]!=0&&pbc[i]!=1){
				printf("Error: Invalid PBC condition..\n\n");
				sprintf(temp,"");
				j=0;
			}
		}
	}
*/
	printf("Minimum aspect ratio for columnar grains?\n");
	scanf("%lf",&ratio);

	aratio[0]=1/ratio;
	aratio[1]=ratio;

	for(k=0;k<dir[2];k++){
	    for(j=0;j<dir[1];j++){
		for(i=0;i<dir[0];i++){
			if(grid[i][j][k][0]!=LIQORI && checker[i][j][k]==0){ // grains not yet measured, not liquid
				gnum++;
				lim[0]=i; // x min
				lim[1]=i; // x max
				lim[2]=j; // y min
				lim[3]=j; // y max
				lim[4]=k; // z min
				lim[5]=k; // z max
				cnt=0;
				DFS_array(dir,i,j,k,grid[i][j][k][0],grid,checker,lim,&cnt,pbc);
				len[0]=lim[1]-lim[0]+1; // x len
				len[1]=lim[3]-lim[2]+1; // y len
				len[2]=lim[5]-lim[4]+1; // z len
				if(len[0]>=len[1]){
					if(len[0]<=len[2]) // len[2] > len[0] > len[1]
						ratio=(double)len[2]/(double)len[1];
					else{ // len [2] < len [0] > len[1]
						if(len[1]>=len[2])
							ratio=(double)len[0]/(double)len[2];
						else
							ratio=(double)len[0]/(double)len[1];
					}
				}else{ // len[0] < len[1]
					if(len[1]<=len[2]) // len[2] > len[1] > len[0]
						ratio=(double)len[2]/(double)len[0];
					else{ // len [2] < len [1] > len[0]
						if(len[0]>=len[2])
							ratio=(double)len[1]/(double)len[2];
						else
							ratio=(double)len[1]/(double)len[0];
					}
				}
				if(ratio>aratio[0]&&ratio<aratio[1]){ // equiaxed
					num[0]++;
					sum[0]=sum[0]+cnt;
				}else{ //columnar
					num[1]++;
					sum[1]=sum[1]+cnt;
					asum=asum+ratio;
				}
			}
		}
	    }
	}

if(gnum==0)
	printf("No solid grains in the sample...\n");
else{
	vavg[0]=(double)(sum[0]+sum[1])/(double)gnum; //total
	vavg[1]=(double)sum[0]/(double)num[0]; //equiaxed
	vavg[2]=(double)sum[1]/(double)num[1]; //columnar

// STDEV
	printf("Assumption of Grain Shape? (1 = sphere, 2 = cube)\n >> ");
	scanf("%d",&i);

	if(i==1){
	    ravg[0]=pow(0.75*vavg[0]/M_PI,1.0/3.0)*2.0; //total
	    ravg[1]=pow(0.75*vavg[1]/M_PI,1.0/3.0)*2.0; // equiaxed
	    sprintf(term,"diameter");
	}else{
	    ravg[0]=pow(vavg[0],1.0/3.0); //total
	    ravg[1]=pow(vavg[1],1.0/3.0); // equiaxed
	    sprintf(term,"size");
	}

	for(k=dir[2]-1;k>=0;k--){
	    for(j=dir[1]-1;j>=0;j--){
		for(i=dir[0]-1;i>=0;i--){
			checker[i][j][k]=0;
	}   }   }

	for(k=0;k<dir[2];k++){
	    for(j=0;j<dir[1];j++){
		for(i=0;i<dir[0];i++){
		    if(grid[i][j][k][0]!=LIQORI && checker[i][j][k]==0){ // grains not yet measured, not liquid
			lim[0]=i; // x min
			lim[1]=i; // x max
			lim[2]=j; // y min
			lim[3]=j; // y max
			lim[4]=k; // z min
			lim[5]=k; // z max
			cnt=0;
			DFS_array(dir,i,j,k,grid[i][j][k][0],grid,checker,lim,&cnt,pbc);
			len[0]=lim[1]-lim[0]+1; // x len
			len[1]=lim[3]-lim[2]+1; // y len
			len[2]=lim[5]-lim[4]+1; // z len
			if(len[0]>=len[1]){
				if(len[0]<=len[2]) // len[2] > len[0] > len[1]
					ratio=(double)len[2]/(double)len[1];
				else{ // len [2] < len [0] > len[1]
					if(len[1]>=len[2])
						ratio=(double)len[0]/(double)len[2];
					else
						ratio=(double)len[0]/(double)len[1];
				}
			}else{ // len[0] < len[1]
				if(len[1]<=len[2]) // len[2] > len[1] > len[0]
					ratio=(double)len[2]/(double)len[0];
				else{ // len [2] < len [1] > len[0]
					if(len[0]>=len[2])
						ratio=(double)len[1]/(double)len[2];
					else
						ratio=(double)len[1]/(double)len[0];
				}
			}
			// cnt == volume (voxels)
/*			if(term[0]=='d') // sphere assumption
			    size=pow(0.75*cnt/M_PI,1.0/3.0)*2.0;
			else
			    size=pow(cnt,1.0/3.0);
			dev[0]+=((size-ravg[0])*(size-ravg[0]));
			if(ratio>aratio[0]&&ratio<aratio[1]) // equiaxed
			    dev[1]+=((size-ravg[1])*(size-ravg[1]));*/
			dev[0]+=((cnt-vavg[0])*(cnt-vavg[0]));
			if(ratio>aratio[0]&&ratio<aratio[1]) // equiaxed
			    dev[1]+=((cnt-vavg[1])*(cnt-vavg[1]));
		    }
		}
	    }
	}

	if(term[0]=='d'){ // sphere assumption
	    dev[0]=(2.0/3.0)*sqrt(dev[0]/(double)gnum)/(ravg[0]*ravg[0]*0.25); // for sphere diameter
	    dev[1]=(2.0/3.0)*sqrt(dev[1]/(double)gnum)/(ravg[1]*ravg[1]*0.25); // for sphere diameter
	}else{ // cube
	    dev[0]=(1.0/3.0)*sqrt(dev[0]/(double)gnum)/(ravg[0]*ravg[0]); // for cube size
	    dev[1]=(1.0/3.0)*sqrt(dev[1]/(double)gnum)/(ravg[1]*ravg[1]); // for cube size
	}

	printf("\n# Total %d grains, avg. grain %s %.2f ±%.2E  x unit length #\n",\
			gnum,term,ravg[0],dev[0]);
	printf("Equiaxed: total %d grains, total volume fraction %f\n       \
	avg. volume %.2f voxels (avg. %s %.2f ±%.2E  x unit length)\n",\
			num[0],(double)sum[0]/(double)(sum[0]+sum[1]),\
	vavg[1],term,ravg[1],dev[1]);
	printf("Columnar: total %d grains, total volume fraction %f\n      \
	avg. volume %.2f voxels (avg. aspect ratio %.2f)\n\n",num[1],(double)sum[1]/(double)(sum[0]+sum[1]),\
	vavg[1],asum/(double)num[1]);
}
    }
	for(i=0;i<dir[0];i++){
		for(j=0;j<dir[1];j++)
			free(checker[i][j]);
		free(checker[i]);
	}
	free(checker);
	return ratio;
}

void DFS_array_simple(int con[],int x,int y,int z,int gori,int**** const grid,int*** checker,int* cnt,int pbc[]){
// Does not consider diagonal direction connection since such a connect bridge may unstable considering GB curvature effect
	checker[x][y][z]=1;
	*cnt=*cnt+1;
	
    if(z<con[2]-1){
	if(gori==grid[x][y][z+1][0]&&checker[x][y][z+1]==0)
		DFS_array_simple(con,x,y,z+1,gori,grid,checker,cnt,pbc);
    }
    if(y<con[1]-1){
	if(gori==grid[x][y+1][z][0]&&checker[x][y+1][z]==0)
		DFS_array_simple(con,x,y+1,z,gori,grid,checker,cnt,pbc);
    }
    if(x<con[0]-1){
	if(gori==grid[x+1][y][z][0]&&checker[x+1][y][z]==0)
		DFS_array_simple(con,x+1,y,z,gori,grid,checker,cnt,pbc);
    }
    if(z>0){
	if(gori==grid[x][y][z-1][0]&&checker[x][y][z-1]==0)
		DFS_array_simple(con,x,y,z-1,gori,grid,checker,cnt,pbc);
    }
    if(y>0){
	if(gori==grid[x][y-1][z][0]&&checker[x][y-1][z]==0)
		DFS_array_simple(con,x,y-1,z,gori,grid,checker,cnt,pbc);
    }
    if(x>0){
	if(gori==grid[x-1][y][z][0]&&checker[x-1][y][z]==0)
		DFS_array_simple(con,x-1,y,z,gori,grid,checker,cnt,pbc);
    }
    return;
}

void DFS_array(int con[],int x,int y,int z,int gori,int**** const grid,int*** checker,int lim[],int* cnt,int pbc[]){
// Does not consider diagonal direction connection since such a connect bridge may unstable considering GB curvature effect
	checker[x][y][z]=1;
	
    if(z<con[2]-1){
	if(gori==grid[x][y][z+1][0]&&checker[x][y][z+1]==0){
		if(lim[5]<z)
			lim[5]=z;
		DFS_array(con,x,y,z+1,gori,grid,checker,lim,cnt,pbc);
	}
    }
    if(y<con[1]-1){
	if(gori==grid[x][y+1][z][0]&&checker[x][y+1][z]==0){
		if(lim[3]<y)
			lim[3]=y;
		DFS_array(con,x,y+1,z,gori,grid,checker,lim,cnt,pbc);
	}
    }
    if(x<con[0]-1){
	if(gori==grid[x+1][y][z][0]&&checker[x+1][y][z]==0){
		if(lim[1]<x)
			lim[1]=x;
		DFS_array(con,x+1,y,z,gori,grid,checker,lim,cnt,pbc);
	}
    }
    if(z>0){
	if(gori==grid[x][y][z-1][0]&&checker[x][y][z-1]==0){
		if(lim[4]>y)
			lim[4]=y;
		DFS_array(con,x,y,z-1,gori,grid,checker,lim,cnt,pbc);
	}
    }
    if(y>0){
	if(gori==grid[x][y-1][z][0]&&checker[x][y-1][z]==0){
		if(lim[2]>y)
			lim[2]=y;
		DFS_array(con,x,y-1,z,gori,grid,checker,lim,cnt,pbc);
	}
    }
    if(x>0){
	if(gori==grid[x-1][y][z][0]&&checker[x-1][y][z]==0){
		if(lim[0]>x)
			lim[0]=x;
		DFS_array(con,x-1,y,z,gori,grid,checker,lim,cnt,pbc);
	}
    }
	*cnt=*cnt+1;
    return;
}

void IFposition(int**** const grid,int dir[],int cor,int start){ // find S-L interface position
	int i,j,k,tot,save[2];
	double sum[2]={0};

	if(abs(cor)==1){ // x dir
		save[0]=dir[0]; // min pos
		save[1]=0;	// max pos
		sum[0]=0;	// sum to avg
		for(k=0;k<dir[2];k++){
		    for(j=0;j<dir[1];j++){
		      if(cor>0){ // solid grows to + direction
			for(i=dir[0]-1;i>=0;i--){
			    if(grid[i][j][k][0]!=LIQORI){
				if(i>save[1])
				    save[1]=i;
				else if(i<save[0])
				    save[0]=i;
				sum[0]+=i;
				break;
			    }
		    	}
		      }else{ // reverse direction
			for(i=0;i<dir[0];i++){
			    if(grid[i][j][k][0]!=LIQORI){
				if(i>save[1])
				    save[1]=i;
				else if(i<save[0])
				    save[0]=i;
				sum[0]+=i;
				break;
			    }
		    	}
		      }
		    }
		}
		if(save[0]==dir[0]) // no liquid
			save[0]=0;

		sum[0]=sum[0]/(double)(dir[1]*dir[2]); // avg.
		// STDEV
		for(k=0;k<dir[2];k++){
		    for(j=0;j<dir[1];j++){
			for(i=dir[0]-1;i>=0;i--){
			    if(grid[i][j][k][0]!=LIQORI){
				sum[1]+=(sum[0]-i)*(sum[0]-i);
				break;
			    }
		    	}
			if(i<0) // No liquid at that line
				sum[1]+=sum[0]*sum[0];
		    }
		}
		sum[1]=sqrt(sum[1]/(double)(dir[1]*dir[2])); // stdev.
	}
	printf("  # Miminum position %d, Maximum position %d, Avg. position %.1f, STDEV %.2f\n",save[0],save[1],sum[0],sum[1]);
	return;
}

void pdas(int**** tgrid,char fnam[],int dir[]){
	int i,cor,max;
	char c,info[50]="\0";

	printf("Dendrite growth direction? (1 = X, 2 = Y, 3 = Z)\n");
	scanf("%d",&cor);
	cor--;

	switch(cor){
		case 0 :
			c='X';
			max=dir[0];
			break;
		case 1 :
			c='Y';
			max=dir[1];
			break;
		case 2 :
			c='Z';
			max=dir[2];
			break;
		default:
			printf("Wrong input...\n\n");
			return;
			break;
	}

//	printf("From where? (%c 0 ~ %d, unit voxel)\n >> ",c,max);
//	scanf("%d",&start);

	for(i=0;i<max;i++)
		printf("%d    %d\n",i,(int)grainsize(tgrid,dir,info,i,cor));

	return;
}

double MCS_time_acceleration(double dt[],double factor[],double Pmax[],int mode){
	if(mode==0){ // Maximum dt method
	    if(Pmax[1]==0){ // No liquid (to be solidified) in the simulation
	    // Acceleration of pure grain growth simulation
		factor[0]=PROBSAFE/Pmax[0];
		factor[1]=1.0;
		return dt[0]*factor[0];
	    }else if(Pmax[1]==1.0){ // Liquid with T<Tk condition exist
		if(dt[0]>dt[1]){
		    factor[0]=dt[1]/dt[0]; //<1
		    factor[1]=1.0;
		    return dt[1];
		}else{ // dt[0] < dt[1]
		    factor[0]=PROBSAFE/Pmax[0];
		    if(dt[1]<factor[0]*dt[0]){ // 다시 조정; maximum은 무조건 <= dt[1]
			factor[0]=dt[1]/dt[0];
			return dt[1];
		    }else{ // factor[0] * dt[0] == factor[1] * dt[1] 만들기 위해 factor[1] 지정
			factor[1]=factor[0]*dt[0]/dt[1];
			return factor[0]*dt[0];
		    }
		}
	    }

	    if(Pmax[0]==0){ // No solid in the simulation
	    // Acceleration of pure solidification
		factor[0]=1.0;
		factor[1]=PROBSAFE/Pmax[1];
		return dt[1]*factor[1];
	    }else if(Pmax[0]==1.0){ // Solid with T>Tm condition exist
		if(dt[0]<dt[1]){
		    factor[0]=1.0;
		    factor[1]=dt[0]/dt[1]; //<1
		    return dt[0];
		}else{ // dt[0] > dt[1]
		    factor[1]=PROBSAFE/Pmax[1];
		    if(dt[0]<factor[1]*dt[1]){ // 다시 조정; maximum은 무조건 <=dt[0]
			factor[0]=dt[1]/dt[0];
			return dt[0];
		    }else{ // factor[0] * dt[0] == factor[1] * dt[1] 만들기 위해 factor[0] 지정
			factor[0]=factor[1]*dt[1]/dt[0];
			return factor[0]*dt[0];
		    }
		}
	    }
// Normal P values, nor 0 and 1
	    factor[0]=PROBSAFE/Pmax[0];
	    factor[1]=PROBSAFE/Pmax[1];
	    if(factor[0]*dt[0]<factor[1]*dt[1]){ // Smallest among two maximum
		factor[1]=factor[0]*dt[0]/dt[1];
		return factor[0]*dt[0];
	    }else{
		factor[0]=factor[1]*dt[1]/dt[0];
		return factor[1]*dt[1];
	    }
	}else if(mode==1){ // Smallest dt method
	    if(dt[0]>dt[1]){
		factor[0]=dt[1]/dt[0];
		factor[1]=1.0;
		return dt[1];
	    }else{
		factor[0]=1.0;
		factor[1]=dt[0]/dt[1];
		return dt[0];
	    }
	}else if(mode==2){// dt = grain growth
	    factor[0]=1.0;
	    factor[1]=1.0;
	    return dt[0];
	}
}

int find_keyword(FILE *file,char word[],char option[]){
	char line[100],temp[100];
	int i=0,j=0;

	rewind(file);
	while(fgets(line,sizeof(line),file)!=NULL){
	    if(line[0]!='#'){
		j=0;
		for(i=0;line[i]!='\0';i++){ // delete '='
		    if(line[i]!='=')
			line[j++]=line[i];
		}
		line[j]='\0';
		tok(line,temp);
		if(strcmp(temp,word)==0){
			i=0;
			break;
		}
	    }
	}
	if(i==0){
		line[strcspn(line,"\n")]='\0';
		strcpy(option,line);
		return 0;
	}else
		return 1;
}

void find_range(int**** grid,double*** tgrid,int dir[],int pbc[],double temp[]){
	int i,j,k,l,m,n,o,p,liq=0;
	double ptcl=0;

	  temp[0]=TNOLIQUID;
	  temp[1]=TNOSOLID;
	  for(i=dir[0]-1;i>=0;i--){
	    for(j=dir[1]-1;j>=0;j--){
		for(k=dir[2]-1;k>=0;k--){
		    if(grid[i][j][k][0]==LIQORI){ //LIQUID
			if(temp[0]>tgrid[i][j][k])
			    temp[0]=tgrid[i][j][k]; //Lowest T among liquid
		    }else if(grid[i][j][k][0]!=PTCLORI){ // SOLID
			if(temp[1]<tgrid[i][j][k])
			    temp[1]=tgrid[i][j][k]; // Highest T among solid
		    }
		}
	    }
	  }

	return;
}

double Tcalc_MCS(char func[],int mcs){ // calculate T from t(MCS)
	char temp[50];
	sprintf(temp,"%s",func);
	time_ins(temp,(double)mcs);
	return eval(temp);
}

double Tcalc_sec(char func[],double time){ // calculate T from t(sec)
	char temp[50];
	sprintf(temp,"%s",func);
	time_ins(temp,time);
	return eval(temp);
}


