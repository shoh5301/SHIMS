#include "SampleMaker.h"

void gmake(int dir[],int**** grid){// Grain maker module
	int dir2[3]={0},tdir[3]={0},multi[3]={1},gnum[3]={0},orient,i,j,k,l,mode=0,seed[3]={0};
	int gnz[2]={-1},rot[2][2]={0};
	int* nsd=NULL;
	int**** tgrid=NULL;
	int**** tgrid2=NULL;
	char fnam[100],fnam2[100];

	printf("\nEdit mode? (1 = Make new / 2 = Multiply / 3 = Crop / 4 = Merge / 5 = Rotate // 0 = Quit )\n >> ");
	scanf("%d",&mode);

	srand(time(NULL));

if(mode==0)
	return;
else if(mode==5){	// Rotater
	i=get_file(fnam,0);
	if(i==0)
		return;
	tgrid=read_info(tgrid,tdir,fnam);
	printf("\nRotation axis? (1 = X axis / 2 = Y axis / 3 = Z axis)\n");
	orient=0;
	while(orient<1 || orient>4){
		scanf("%d",&orient);
		if(orient>0 && orient<4)
			break;
		printf("Wrong input... \n");
	}
	printf("\n90' Rotation...\n");

	if(orient==1){  // X axis
		dir[0]=tdir[0];
		dir[1]=tdir[2];
		dir[2]=tdir[1];
		grid=alloc4d(dir[0],dir[1],dir[2],grid); // grid allocation
		for(i=0;i<tdir[0];i++){
			for(j=0;j<tdir[1];j++){
				for(k=0;k<tdir[2];k++)
					voxcpy(grid[i][tdir[2]-k-1][j],tgrid[i][j][k]);
		}	}
	}else if(orient==2){ // Y axis
		dir[0]=tdir[2];
		dir[1]=tdir[1];
		dir[2]=tdir[0];
		grid=alloc4d(dir[0],dir[1],dir[2],grid); // grid allocation
		for(i=0;i<tdir[0];i++){
			for(j=0;j<tdir[1];j++){
				for(k=0;k<tdir[2];k++)
					voxcpy(grid[tdir[2]-k-1][j][i],tgrid[i][j][k]);
		}	}
	}else{		// Z axis
		dir[0]=tdir[1];
		dir[1]=tdir[0];
		dir[2]=tdir[2];
		grid=alloc4d(dir[0],dir[1],dir[2],grid); // grid allocation
		for(i=0;i<tdir[0];i++){
			for(j=0;j<tdir[1];j++){
				for(k=0;k<tdir[2];k++)
					voxcpy(grid[tdir[1]-j-1][i][k],tgrid[i][j][k]);
		}	}
	}

	tgrid=free4d(tdir[0],tdir[1],tdir[2],tgrid);
	fileout(dir,grid,-1);
	grid=free4d(dir[0],dir[1],dir[2],grid);
	return;

}else if(mode==4){  // Paster
	printf("\nInput file name of 1st file");
	i=get_file(fnam,0);
	if(i==0)
		return;
	printf("\n1st file confirmed...\nPlease input 2nd file");
	i=get_file(fnam2,0);
	if(i==0)
		return;
	printf("\n Merge '%s' and '%s'...\n",fnam,fnam2);

	while(1){
		printf("What is the direction to attach 2nd file?\n(x = 1, y = 2, z = 3)\n >> ");
		scanf("%d",&orient);
		if(orient<1||orient>3)
			printf("Wrong input...\n\n");
		else
			break;
	}
	tgrid=read_info(tgrid,dir,fnam);
	tgrid2=read_info(tgrid2,dir2,fnam2);

	tdir[0]=lcm(dir[0],dir2[0]);
	tdir[1]=lcm(dir[1],dir2[1]);
	tdir[2]=lcm(dir[2],dir2[2]);

	multi[0]=tdir[0]/dir[0];
	multi[1]=tdir[1]/dir[1];
	multi[2]=tdir[2]/dir[2];

	if(orient==1){
		multi[0]=1;
		tdir[0]=dir[0]+dir2[0];
	}else if(orient==2){
		multi[1]=1;
		tdir[1]=dir[1]+dir2[1];
	}else{
		multi[2]=1;
		tdir[2]=dir[2]+dir2[2];
	}
	if(multi[0]!=1 || multi[1]!=1 || multi[2]!=1){
		tgrid=free4d(dir[0],dir[1],dir[2],tgrid);
		printf("Re-create sample1 (matching PBC)\n");
		tgrid=multiply(tgrid,dir,fnam,multi,1);
	}
	multi[0]=tdir[0]/dir2[0];
	multi[1]=tdir[1]/dir2[1];
	multi[2]=tdir[2]/dir2[2];
	if(orient==1)
		multi[0]=1;
	else if(orient==2)
		multi[1]=1;
	else
		multi[2]=1;
	if(multi[0]!=1 || multi[1]!=1 || multi[2]!=1){
		tgrid2=free4d(dir2[0],dir2[1],dir2[2],tgrid2);
		printf("Re-create sample2 (matching PBC)\n");
		tgrid2=multiply(tgrid2,dir2,fnam2,multi,1);
	}

	grid=alloc4d(tdir[0],tdir[1],tdir[2],grid); // grid allocation
	for(i=0;i<tdir[0];i++){
		for(j=0;j<tdir[1];j++){
			for(k=0;k<tdir[2];k++){
				if(orient==1){
					if(i<dir[0])
						voxcpy(grid[i][j][k],tgrid[i][j][k]);
					else
						voxcpy(grid[i][j][k],tgrid2[i-dir[0]][j][k]);
				}else if(orient==2){
					if(j<dir[1])
						voxcpy(grid[i][j][k],tgrid[i][j][k]);
					else
						voxcpy(grid[i][j][k],tgrid2[i][j-dir[1]][k]);
				}else{
					if(k<dir[2])
						voxcpy(grid[i][j][k],tgrid[i][j][k]);
					else
						voxcpy(grid[i][j][k],tgrid2[i][j][k-dir[2]]);
				}
	}	}	}
	tgrid=free4d(dir[0],dir[1],dir[2],tgrid);
	tgrid2=free4d(dir2[0],dir2[1],dir2[2],tgrid2);
	dir[0]=tdir[0];
	dir[1]=tdir[1];
	dir[2]=tdir[2];

	fileout(dir,grid,-1);
	grid=free4d(dir[0],dir[1],dir[2],grid);
	return;
}else if(mode==3){	// Cutter
	i=get_file(fnam,0);
	if(i==0)
		return;
	tgrid=read_info(tgrid,tdir,fnam);

	printf("\nCutting coordinate:\n");
	printf("\n X range [0,%d]",tdir[0]-1);
	while(1){
		printf("\nX from = ");
		scanf("%d",&i);
		if(i>=0 && i<tdir[0]-1)
			break;
		printf("Wrong initial coordinate: not in the range [0,%d)\n",tdir[0]-1);
	}while(1){
		printf("X to = ");
		scanf("%d",&j);
		if(j>i && j<tdir[0])
			break;
		printf("Wrong initial coordinate: not in the range (%d,%d]\n",i,tdir[0]-1);
	}
	printf("\n Y range [0,%d]",tdir[1]-1);
	dir2[0]=i; //min
	tdir[0]=j; //max
	while(1){
		printf("\nY from = ");
		scanf("%d",&i);
		if(i>=0 && i<tdir[1]-1)
			break;
		printf("Wrong initial coordinate: not in the range [0,%d)\n",tdir[1]-1);
	}while(1){
		printf("Y to = ");
		scanf("%d",&j);
		if(j>i && j<tdir[1])
			break;
		printf("Wrong initial coordinate: not in the range (%d,%d]\n",i,tdir[1]-1);
	}
	dir2[1]=i; //min
	tdir[1]=j; //max
	printf("\n Z range [0, %d]",tdir[2]-1);
	while(1){
		printf("\nZ from = ");
		scanf("%d",&i);
		if(i>=0 && i<tdir[2]-1)
			break;
		printf("Wrong initial coordinate: not in the range [0,%d)\n",tdir[2]-1);
	}while(1){
		printf("Z to = ");
		scanf("%d",&j);
		if(j>i && j<tdir[2])
			break;
		printf("Wrong initial coordinate: not in the range (%d,%d]\n",i,tdir[2]-1);
	}
	dir2[2]=i; //min
	tdir[2]=j; //max
	printf("X range [%d, %d]\n",dir2[0],tdir[0]);
	printf("Y range [%d, %d]\n",dir2[1],tdir[1]);
	printf("Z range [%d, %d]\n\n",dir2[2],tdir[2]);

	dir[0]=tdir[0]-dir2[0]+1;
	dir[1]=tdir[1]-dir2[1]+1;
	dir[2]=tdir[2]-dir2[2]+1;

	grid=alloc4d(dir[0],dir[1],dir[2],grid); // grid allocation

	for(i=dir2[0];i<=tdir[0];i++){
		for(j=dir2[1];j<=tdir[1];j++){
			for(k=dir2[2];k<=tdir[2];k++)
				voxcpy(grid[i-dir2[0]][j-dir2[1]][k-dir2[2]],tgrid[i][j][k]);
	}	}

	tgrid=free4d(tdir[0],tdir[1],tdir[2],tgrid);
	fileout(dir,grid,-1);
	grid=free4d(dir[0],dir[1],dir[2],grid);
	return;
    }else if(mode==2){	// Multiplier
	i=get_file(fnam,0);
	if(i==0)
		return;
	grid=read_info(grid,dir,fnam);
	grid=free4d(dir[0],dir[1],dir[2],grid);
	dir[0]=0;
	printf("\nMultiply factor to X direction?(>=1, integer)\n >> ");
	scanf("%d",&multi[0]);
	printf("Multiply factor to Y direction?(>=1, integer)\n >> ");
	scanf("%d",&multi[1]);
	printf("Multiply factor to Z direction?(>=1, integer)\n >> ");
	scanf("%d",&multi[2]);
	grid=multiply(grid,dir,fnam,multi,0);
	grid=free4d(dir[0],dir[1],dir[2],grid);
	return;
    }else if(mode==1){	// Maker
	printf("\nSample size (x direction, voxels)?: ");
	scanf("%d",&dir[0]);
	printf("Sample size (y direction, voxels)?: ");
	scanf("%d",&dir[1]);
	printf("Sample size (z direction, voxels)?: ");
	scanf("%d",&dir[2]);

	printf("\nSample type? (0 = Solid, 1 = Liquid, 2 = Mould-liquid, 3 = Powder)\n >> ");
//	, 4 = single nucleus among liquid)\n >> ");
	scanf("%d",&gnum[0]);

	if(gnum[0]!=0)
		grid=alloc4d(dir[0],dir[1],dir[2],grid); // grid allocation

	if(gnum[0]==1){
	    for(i=0;i<dir[0];i++){
		for(j=0;j<dir[1];j++){
			for(k=0;k<dir[2];k++)
				neworient_melting(grid[i][j][k]);
	    }	}
	    printf("Liquid box created...\n");
	    fileout(dir,grid,-1);
	    return;
	}else if(gnum[0]==2){
	    for(i=0;i<dir[0];i++){
		for(j=0;j<dir[1];j++){
			for(k=0;k<dir[2];k++)
				neworient_melting(grid[i][j][k]);
	    }	}
	    printf("Small nuclei at X-lower plane? (1 to y)\n >> ");
	    scanf("%d",&gnum[0]);
	    if(gnum[0]==1){
		for(j=0;j<dir[1];j++){
			for(k=0;k<dir[2];k++)
				neworient(grid[0][j][k]);
		}
	    }
	    printf("Small nuclei at X-upper plane? (1 to y)\n >> ");
	    scanf("%d",&gnum[0]);
	    if(gnum[0]==1){
		i=dir[0]-1;
		for(j=0;j<dir[1];j++){
			for(k=0;k<dir[2];k++)
				neworient(grid[i][j][k]);
		}
	    }
	    printf("Small nuclei at Y-lower plane? (1 to y)\n >> ");
	    scanf("%d",&gnum[0]);
	    if(gnum[0]==1){
		for(i=0;i<dir[0];i++){
			for(k=0;k<dir[2];k++)
				neworient(grid[i][0][k]);
		}
	    }
	    printf("Small nuclei at Y-upper plane? (1 to y)\n >> ");
	    scanf("%d",&gnum[0]);
	    if(gnum[0]==1){
		j=dir[1]-1;
		for(i=0;i<dir[0];i++){
			for(k=0;k<dir[2];k++)
				neworient(grid[i][j][k]);
		}
	    }
	    printf("Small nuclei at Z-lower plane? (1 to y)\n >> ");
	    scanf("%d",&gnum[0]);
	    if(gnum[0]==1){
		for(i=0;i<dir[0];i++){
			for(j=0;j<dir[1];j++)
				neworient(grid[i][j][0]);
		}
	    }
	    printf("Small nuclei at Z-upper plane? (1 to y)\n >> ");
	    scanf("%d",&gnum[0]);
	    if(gnum[0]==1){
		k=dir[2]-1;
		for(i=0;i<dir[0];i++){
			for(j=0;j<dir[1];j++)
				neworient(grid[i][j][k]);
		}
	    }
	    fileout(dir,grid,-1);
	    return;
/*	}else if(gnum[0]==3){ // Powder
	    for(i=0;i<dir[0];i++){
		for(j=0;j<dir[1];j++){
		    for(k=0;k<dir[2];k++)
			neworient_powder(grid[i][j][k]);
		}
	    }
	    fileout(dir,grid,-1);
	    return;*/
	}else if(gnum[0]==4){ // single nucleus
	    for(i=0;i<dir[0];i++){
		for(j=0;j<dir[1];j++){
			for(k=0;k<dir[2];k++)
				neworient_melting(grid[i][j][k]);
	    }	}
	    printf("Size of a nucleus? (unit voxel)\n r = ");
	    scanf("%d",&gnum[0]);

	    printf("At x-low surface vs. center? (1 vs 2)\n >> ");
	    scanf("%d",&i);
	    if(i==1){ // x-low surface
		seed[0]=0;
		seed[1]=dir[1]/2;
		seed[2]=dir[2]/2;
	    }else{
		seed[0]=dir[0]/2;
		seed[1]=dir[1]/2;
		seed[2]=dir[2]/2;
	    }
	    neworient(grid[seed[0]][seed[1]][seed[2]]);
	    orient=grid[seed[0]][seed[1]][seed[2]][0];
	    for(i=seed[0];i<seed[0]+gnum[0];i++){
		for(j=seed[1];j<seed[1]+gnum[0];j++){
		    for(k=seed[2];k<seed[2]+gnum[0];k++)
			grid[i][j][k][0]=orient;
	    }	}
	    fileout(dir,grid,-1);
	    return;
	}else{ // gnum==0 && gnum==3
	    mode=gnum[0];
	    initial_size_setting(gnum,dir);

	    grid=alloc4d(dir[0],dir[1],dir[2],grid); // grid allocation

	    for(l=0;l<dir[0];l=l+gnum[0]){
		for(i=0;i<dir[1];i=i+gnum[1]){
			gnz[0]=-1;
			for(k=0;k<dir[2];k++){ //Z direction
				gnz[1]=k/gnum[2];
				if(gnz[0]!=gnz[1])
					neworient(seed);
				voxcpy(grid[l][i][k],seed);
				gnz[0]=gnz[1];
			}
			for(j=i+1;j<i+gnum[1];j++){ // copy to Y direction
				for(k=0;k<dir[2];k++)
					voxcpy(grid[l][j][k],grid[l][i][k]);
			}
		}
		for(i=l+1;i<l+gnum[0];i++){ // copy to X diection
			for(j=0;j<dir[1];j++){
				for(k=0;k<dir[2];k++)
					voxcpy(grid[i][j][k],grid[l][j][k]);
		}	}
	    }
	    if(mode==0)
		printf("\nRandom Grain Created...\n");
	    else{
		for(i=dir[0]-1;i>=0;i--){
		    for(j=dir[1]-1;j>=0;j--){
			for(k=dir[2]-1;k>=0;k--){
			    grid[i][j][k][0]+=POWORI;
			    grid[i][j][k][1]=POWORI;
			    grid[i][j][k][1]=POWORI;
		}   }	}
		printf("\nRandom Powder Pool Created...\n");
	    }
	    fileout(dir,grid,-1);
	    return;
	    }
	}
	return;
}

int lcm(int a,int b){ // Return LCM value, a is bigger
	int tmp,tmpa,tmpb;
	if(a==b){
		return a;
	}else if(a<b){ //Swap
		tmp=a;
		a=b;
		b=tmp;
	}
	tmpa=a;
	tmpb=b;
	while(b!=0){
		tmp=a%b;
		a=b;
		b=tmp;
	} // a==gcd
	return (tmpa*tmpb)/a; //lcm
}

void initial_size_setting(int gnum[],int dir[]){
	int i,sn;

	if(gnum[0]==3) // Powder
	    printf("\nSize of powder particles? (1 powder particle will take Nx x Ny x Nz voxels)\n Nx = ");
	else
	    printf("\nSize of initial grains? (1 grain will take Nx x Ny x Nz voxels)\n Nx = ");

	    scanf("%d",&gnum[0]);
	    if(gnum[0]==0){
		printf(" ### WARNING: initial size should not be 0...\n     Automatically set to 1\n");
		gnum[0]=1;
	    }
	    printf(" Ny = ");
	    scanf("%d",&gnum[1]);
	    if(gnum[1]==0){
		printf(" ### WARNING: initial size should not be 0...\n     Automatically set to 1\n");
		gnum[1]=1;
	    }
	    printf(" Nz = ");
	    scanf("%d",&gnum[2]);
	    if(gnum[2]==0){
		printf(" ### WARNING: initial size should not be 0...\n     Automatically set to 1\n");
		gnum[2]=1;
	    }

	    i=0;
	    sn=dir[0]%gnum[0];
	    if(sn!=0){
		printf(" ### WARNING: sample size cannot be devided by particle size (X direction)...\n\tSample size is automatically re-adjusted: %d →",dir[0]);
		i++;
	    }
	    while(sn!=0){ // correct dir number
		dir[0]++;
		sn=dir[0]%gnum[0];
	    }
	    if(i!=0){
		printf("%d\n",dir[0]);
		i=0;
	    }
	    sn=dir[1]%gnum[1];
	    if(sn!=0){
		printf(" ### WARNING: sample size cannot be devided by particle size (Y direction)...\n\tSample size is automatically re-adjusted: %d →",dir[1]);
		i++;
	    }
	    while(sn!=0){ // correct dir number
		dir[1]++;
		sn=dir[1]%gnum[1];
	    }
	    if(i!=0){
		printf("%d\n",dir[1]);
		i=0;
	    }
	    sn=dir[2]%gnum[2];
	    if(sn!=0){
		printf(" ### WARNING: sample size cannot be devided by particle size (Z direction)...\n\tSample size is automatically re-adjusted: %d →",dir[2]);
		i++;
	    }
	    while(sn!=0){ // correct dir number
		dir[2]++;
		sn=dir[2]%gnum[2];
	    }
	    if(i!=0){
		printf("%d\n",dir[2]);
		i=0;
	    }
	    return;
}

int**** multiply(int**** grid,int dir[],char fnam[],int multi[],int mod){// Multiply grain structures
	int i,j,k,l,seed,m[3]={0};
	int* nsd0=NULL;
	int* nsd1=NULL;
	int* nsd2=NULL;
	int**** tgrid=NULL;
	int**** nori=NULL;

	tgrid=read_info(tgrid,dir,fnam);
	m[0]=multi[0]*dir[0];
	m[1]=multi[1]*dir[1];
	m[2]=multi[2]*dir[2];

	printf("Apply new orientation? (1 to yes)\n >> ");
	scanf("%d",&seed);

if(seed==1){
	seed=0;
	printf("Start to count orientations...\n");
	nsd0=(int*)malloc(dir[0]*dir[1]*dir[2]*sizeof(int));
	nsd1=(int*)malloc(dir[0]*dir[1]*dir[2]*sizeof(int));
	nsd2=(int*)malloc(dir[0]*dir[1]*dir[2]*sizeof(int));
	for(i=0;i<dir[0];i++){
		for(j=0;j<dir[1];j++){
			for(k=0;k<dir[2];k++){ // Count # of orientations
				if(seed==0){
					nsd0[seed]=tgrid[i][j][k][0];
					nsd1[seed]=tgrid[i][j][k][1];
					nsd2[seed]=tgrid[i][j][k][2];
	                                seed++;
				}else{
					for(l=0;l<seed;l++){
						if(nsd0[l]==tgrid[i][j][k][0] && nsd1[l]==tgrid[i][j][k][1] && nsd2[l]==tgrid[i][j][k][2])
							break;
					}
					if(l==seed){
						nsd0[seed]=tgrid[i][j][k][0];
						nsd1[seed]=tgrid[i][j][k][1];
						nsd2[seed]=tgrid[i][j][k][2];
						seed++;
					}
				}
			}
		}	
	}

	nori=alloc4d(multi[0],multi[1],multi[2],nori);
	free(nsd0);
	free(nsd1);
	free(nsd2);

	for(i=0;i<multi[0];i++){
		for(j=0;j<multi[1];j++){
			for(k=0;k<multi[2];k++)
				neworient(nori[i][j][k]);
	}	}
	seed=1;
}else
	seed=0;

	grid=alloc4d(m[0],m[1],m[2],grid);
	for(i=0;i<m[0];i++){
		for(j=0;j<m[1];j++){
			for(k=0;k<m[2];k++){
				if(((i/dir[0]==0)&&(j/dir[1]==0))&&(k/dir[2]==0))
					voxcpy(grid[i][j][k],tgrid[i][j][k]);
				else{
					if(tgrid[i%dir[0]][j%dir[1]][k%dir[2]][1]==LIQORI)
						neworient_melting(grid[i][j][k]);
					else if(seed==0)
						voxcpy(grid[i][j][k],tgrid[i%dir[0]][j%dir[1]][k%dir[2]]);
					else{
						grid[i][j][k][0]=tgrid[i%dir[0]][j%dir[1]][k%dir[2]][0]+nori[i/dir[0]][j/dir[1]][k/dir[2]][0];
						if(grid[i][j][k][0]>=MAXDG1)
							grid[i][j][k][0]=grid[i][j][k][0]-MAXDG1;
						grid[i][j][k][1]=tgrid[i%dir[0]][j%dir[1]][k%dir[2]][1]+nori[i/dir[0]][j/dir[1]][k/dir[2]][1];
						if(grid[i][j][k][1]>=MAXDG2)
							grid[i][j][k][1]=grid[i][j][k][1]-MAXDG2;
						grid[i][j][k][2]=tgrid[i%dir[0]][j%dir[1]][k%dir[2]][2]+nori[i/dir[0]][j/dir[1]][k/dir[2]][2];
						if(grid[i][j][k][2]>=MAXDG2)
							grid[i][j][k][2]=grid[i][j][k][2]-MAXDG2;
					}
				}
	}	}	}

	tgrid=free4d(dir[0],dir[1],dir[2],tgrid);
	if(seed==1)
		nori=free4d(multi[0],multi[1],multi[2],nori);

	dir[0]=m[0];
	dir[1]=m[1];
	dir[2]=m[2];

	printf("\nSample was re-created, as %d x %d x %d\n",dir[0],dir[1],dir[2]);
	if(mod==0)
		fileout(dir,grid,-1);
	return grid;
}

