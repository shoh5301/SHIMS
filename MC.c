#include "MC.h"

void MCtrial_GG(int x,int y,int z,int**** grid,int dir[],int pbc[],double jc[],double tbc[][2],double pmob,double realE){
	int l,m,n,o;
	double dE,prob;

	l=x;
	m=y;
	n=z;
	o=(int)(rand()/(RAND_MAX/26+1));
	nei26sel(o,&l,&m,&n,dir,pbc);
	if(grid[l][m][n][0]!=LIQORI && grid[x][y][z][0]!=grid[l][m][n][0]){
		dE=gbenergy(x,y,z,grid[l][m][n],grid,jc,dir,pbc,tbc)-gbenergy(x,y,z,grid[x][y][z],grid,jc,dir,pbc,tbc);
		if(dE<=0)
			prob=pmob;
		else
			prob=pmob*exp(-dE*realE);
//		if(prob>(double)rand()/(RAND_MAX)) // Change orientation
		if(prob_check(prob)!=0)
			voxcpy(grid[x][y][z],grid[l][m][n]);
	}

	return;
}

void MCtrial_SOL(int i,int j,int k,int**** grid,double**** fdmt,int dir[],int pbc[],int mode[],double melt[],double Tmax[],double dHfus[],double jc[],double tbc[][2],double Vsite,double Asite,double factor[],double Nsite[],double pcps[],double Qgg){
	int l,m,n,o,situ,tempgrid[3]={LIQORI,LIQORI,LIQORI};
	double Tnow,dGfus,dEif,dGfor,prob,fhetero_s;

	l=i;
	m=j;
	n=k;
	if(mode[MTEM]==0)
	    Tnow=fdmt[0][i][j][k];
	else // Frozen temperature approximation
	    Tnow=i*melt[4]-melt[3];

	o=(int)(rand()/(RAND_MAX/26+1));
	nei26sel(o,&l,&m,&n,dir,pbc);

	if(grid[i][j][k][0]==LIQORI){
	// Target voxel is liquid -> solidification
	    if(Tnow<=Tmax[1]){
		if(grid[l][m][n][0]==LIQORI && mode[MFIT]!=0){ 
		// Liquid - liquid -> solid nucleation // fittingmode==1 for dendrite growth fitting mode
		// NOTE: in principle, it should be considered in any cases for liquid-liquid interaction,
		//	  but only for T < Tm, for computational efficiency
		    neworient(tempgrid);
		    dGfus=solidi_df(dHfus[0],Tnow,Tmax[0],Tmax[2],melt[8]); //>0 for prefered reaction

		  if(jc[7]!=0 && jc[7]>(double)rand()/(RAND_MAX)){ // Nucl. by Particle
		    dEif=gbenergy(i,j,k,tempgrid,grid,jc,dir,pbc,tbc)-gbenergy(i,j,k,grid[i][j][k],grid,jc,dir,pbc,tbc);
		    dGfor=Asite*dEif-Vsite*dGfus; //<0 for prefered reaction

		    fhetero_s=cos(jc[8]*DTOR); // jc[8] == PANG
		    fhetero_s=0.25*(2+fhetero_s)*(1-fhetero_s)*(1-fhetero_s)*FHOMO;
		   
//		    prob=factor[1]*P_nucl(1/Tnow,dGfus,jc[4],fhetero_s,Nsite[0],dGfor,pcps,dGfus*Vsite);
		    prob=factor[1]*P_nucl(1/Tnow,dGfus,jc[4],fhetero_s,Nsite[0]*jc[7],dGfor,pcps,dGfus*Vsite);
		  }else{
		    situ=nucl_cond(grid,i,j,k,dir,pbc,0);
			    
		    if(situ==0){ // homog. nucl. among liquid atmosphere
			dEif=gbenergy(i,j,k,tempgrid,grid,jc,dir,pbc,tbc)-gbenergy(i,j,k,grid[i][j][k],grid,jc,dir,pbc,tbc);
			dGfor=Asite*dEif-Vsite*dGfus; //<0 for prefered reaction
			prob=factor[1]*P_nucl(1/Tnow,dGfus,jc[4],FHOMO,Nsite[0],dGfor,pcps,dGfus*Vsite);
		    }else if(situ>0){// heterog. nucl. @ S-L interface / powder surface / particle surface
			while(1){ // select interacting voxel
			    o=(int)(rand()/(RAND_MAX/26+1));
			    nei26sel(o,&l,&m,&n,dir,pbc);
			    if(grid[l][m][n][0]!=LIQORI)
				break;
			}
			if(mode[MTEM]==0){
			    if(fdmt[0][l][m][n]!=fdmt[0][i][j][k])
				Tnow=((double)rand()/(RAND_MAX))*(fdmt[0][l][m][n]-fdmt[0][i][j][k])+fdmt[0][i][j][k];
			}else // Frozen temperature approximation; T btw. T* ~ T*-G
			    Tnow=(i+l)*0.5*melt[4]-melt[3];

			dEif=gbenergy(i,j,k,tempgrid,grid,jc,dir,pbc,tbc)-gbenergy(i,j,k,grid[i][j][k],grid,jc,dir,pbc,tbc);
			dGfor=Asite*dEif-Vsite*dGfus; //<0 for prefered reaction

//			jc[7]=jc[4]*dendrite_dir(tempgrid,o,jc[2]); //IFE anisotropy
			o=misorient(tempgrid[0],grid[l][m][n][0]);
			fhetero_s=1-(GBE(o,jc[3],jc[1])/jc[4]);
			if(fhetero_s<-2) // wetting angle > 180 degree -> almost sphere (~homog.)
			    fhetero_s=FHOMO;
			else
			    fhetero_s=0.25*(2+fhetero_s)*(1-fhetero_s)*(1-fhetero_s)*FHOMO;
			prob=factor[1]*P_nucl(1/Tnow,dGfus,jc[4],fhetero_s,Nsite[1],dGfor,pcps,dGfus*Vsite);
//			prob=factor[1]*P_nucl(1/Tnow,dGfus,jc[7],fhetero_s,Nsite[1],dGfor,pcps,dGfus*Vsite);
		    }else{ // situ<0, heterog. nucl. at surface
			situ=-situ-1; // 0 x-1 1 x+1 2 y-1 3 y+1 4 z-1 5 z+1
//			jc[7]=(double)VN_to_Moore(situ);
			dEif=gbenergy(i,j,k,tempgrid,grid,jc,dir,pbc,tbc)-gbenergy(i,j,k,grid[i][j][k],grid,jc,dir,pbc,tbc);
			dGfor=Asite*dEif-Vsite*dGfus; //<0 for prefered reaction

//			jc[7]=jc[4]*dendrite_dir(tempgrid,VN_to_Moore(situ),jc[2]); //IFE anisotropy
			if(tbc[situ][0]<0 || tbc[situ][1]<0){
			// tbc[][0]<0 -> convection // tbc[][1]<0 -> vapor condition
			// surface condition: solid-vapor 1face
			// jc 4 = S-L 5 = S-V 6 = L-V
//			    if(jc[6]>jc[7]+jc[5]) // actually this does not happens since S-V > S-L + L-V generally
			    if(jc[6]>jc[4]+jc[5]) // actually this does not happens since S-V > S-L + L-V generally
				fhetero_s=0;
			    else{
				fhetero_s=(jc[6]-jc[5])/jc[4]; //cos_theta
				fhetero_s=0.25*(2+fhetero_s)*(1-fhetero_s)*(1-fhetero_s)*FHOMO;
			    }
			}else{
			// surface condition: mould surface
			    if(tbc[situ][1]>=180 || tbc[situ][1]<0) // l-v 1face condition
				fhetero_s=FHOMO;
			    else{
				// fhetero_s=cos(DTOR*tbc[situ][1]); // wetting angle fixed
				fhetero_s=cos(DTOR*(((double)rand()/(RAND_MAX))*(180-tbc[situ][1])+tbc[situ][1]));
				// wetting angle among minimum angle ~ 180
				fhetero_s=0.25*(2+fhetero_s)*(1-fhetero_s)*(1-fhetero_s)*FHOMO;
			    }
//printf("Mould fhetero %f vs fHOMO %f ",fhetero_s,FHOMO);
			}
//			prob=factor[1]*P_nucl(1/Tnow,dGfus,jc[7],fhetero_s,Nsite[1],dGfor,pcps,dGfus*Vsite);
			prob=factor[1]*P_nucl(1/Tnow,dGfus,jc[4],fhetero_s,Nsite[1],dGfor,pcps,dGfus*Vsite);
//printf("prob %E where situ %d\n",prob,situ);
		    }
		  }
//		    if(prob>(double)rand()/(RAND_MAX)){ // Nucl.
		    if(prob_check(prob)!=0){
//printf("Accepted! prob %E\n",prob);
			voxcpy(grid[i][j][k],tempgrid);
			if(mode[MLAT]!=0)
			    fdmt[2][i][j][k]+=dHfus[1]; // dHfus release
		    }
		}else if(grid[l][m][n][0]!=LIQORI && grid[l][m][n][0]!=PTCLORI && grid[l][m][n][0]!=MPISURF){
// interacting voxel is solid
/////////////////////////////////// solid-liquid interaction: solid epitaxial growth ////////////////////////////////////////
		    if(mode[MTEM]==0){
			if(fdmt[0][l][m][n]!=fdmt[0][i][j][k])
			    Tnow=((double)rand()/(RAND_MAX))*(fdmt[0][l][m][n]-fdmt[0][i][j][k])+fdmt[0][i][j][k];
			else
			    Tnow=fdmt[0][i][j][k];
		    }else
			Tnow=((double)rand()/(RAND_MAX))*melt[4]-melt[3];

		    if(Tnow<Tmax[1]){ // only consider when tip T < (Tm-Tcs)

			dEif=gbenergy(i,j,k,grid[l][m][n],grid,jc,dir,pbc,tbc)-gbenergy(i,j,k,grid[i][j][k],grid,jc,dir,pbc,tbc);
			dGfus=solidi_df(dHfus[0],Tnow,Tmax[0],Tmax[2],melt[8]); //>0 for prefered reaction
			dGfor=Asite*dEif-Vsite*dGfus; //<0 for prefered reaction
//			prob=factor[1]*P_growth(1/Tnow,pcps,dGfus*Vsite,dGfor);
			prob=factor[1]*P_growth(1/Tnow,pcps,dGfus*Vsite,dGfor)*check_growth_dir(grid,l,m,n,get_nei26_index(l-i,m-j,n-k),dir,pbc,jc[2]);
//printf("check %f\n",check_growth_dir(grid,l,m,n,get_nei26_index(l-i,m-j,n-k),dir,pbc,jc[2]));

//			if(prob>(double)rand()/(RAND_MAX)){ // growth
			if(prob_check(prob)!=0){
			    if(grid[l][m][n][1]==POWORI)
				pow2grain(grid[l][m][n]);
			    voxcpy(grid[i][j][k],grid[l][m][n]);
			    if(mode[MLAT]!=0)
				fdmt[2][i][j][k]+=dHfus[1]; // dHfus release
			}
		    }
		}
	    }//else T > Tm : no 1est
	}else if(grid[i][j][k][0]!=PTCLORI){ // target voxel is solid
//	    if(Tnow<=Tmax[0]){
		if(grid[l][m][n][0]!=PTCLORI && grid[l][m][n][0]!=LIQORI && grid[l][m][n][0]!=MPISURF){// && grid[l][m][n][1]!=POWORI){
		    // T < Tm: GRAIN GROWTH
		    if(grid[i][j][k][0]!=grid[l][m][n][0]){
			dEif=gbenergy(i,j,k,grid[l][m][n],grid,jc,dir,pbc,tbc)-gbenergy(i,j,k,grid[i][j][k],grid,jc,dir,pbc,tbc);
			prob=factor[0]*P_grain_growth(Qgg,1/Tmax[0],1/Tnow,Asite,dEif);

			if(prob>(double)rand()/(RAND_MAX))
//			if(prob_check(prob)!=0) // case that P < PROBMIN by active dt ctrl
			    voxcpy(grid[i][j][k],grid[l][m][n]);
		    }
		    else if(grid[i][j][k][1]==POWORI || grid[l][m][n][1]==POWORI){
 // powder coarsening -> solid grain
			prob=factor[0]*P_grain_growth(Qgg,1/Tmax[0],1/Tnow,Asite,0); //dEif==0
//			if(prob>(double)rand()/(RAND_MAX)){
			if(prob_check(prob)!=0){
			    if(grid[i][j][k][1]==POWORI)
				pow2grain(grid[i][j][k]);
			    if(grid[l][m][n][1]==POWORI)
				pow2grain(grid[l][m][n]);
			}
		    }
		}
	}

	return;
}

double gbenergy(int px,int py,int pz,int now[],int**** const grid,double jc[],int dir[],int pbc[],double tbc[][2]){ // GBE calculator
	double ten=0,del,temp,f_IFE;
	int i,j,k,l,m,n,ori;

	if(now[1]==POWORI)
	    now[0]-=POWORI;

	l=26;
    if(now[0]!=LIQORI){ // now is solid
//	f_IFE=dendrite_dir(now,(int)jc[7],jc[2]); // return aniso. factor for IFE
	for(l=25;l>=0;l--){
	    i=px;
	    j=py;
	    k=pz;
	    m=surf_check(l,i,j,k,dir,pbc);
	    n=mpi_surf_check(l,i,j,k,grid);
	    if(m>=0){
		if(tbc[m][0]<0 || tbc[m][1]<0) // tbc[][0]<0 -> convection // tbc[][1]<0 -> vapor condition
			ten=ten+jc[5]; // surface condition: solid-vapor 1face
		else
			ten=ten+jc[4]*cos(DTOR*tbc[m][1]); // IFE_s-l * cos theta * aniso
	    }else if(n>=0){ //MPI x direction surface condition
		if(tbc[n][0]<0 || tbc[n][1]<0) // tbc[][0]<0 -> convection // tbc[][1]<0 -> vapor condition
			ten=ten+jc[5]; // surface condition: solid-vapor 1face
		else
			ten=ten+jc[4]*cos(DTOR*tbc[n][1]); // IFE_s-l * cos theta * aniso
	    }else{ // not a surface
		nei26sel(l,&i,&j,&k,dir,pbc);
		if(grid[i][j][k][1]==POWORI)
		    ori=grid[i][j][k][0]-POWORI;
		else
		    ori=grid[i][j][k][0];
		if(ori==LIQORI) // solid-liquid interface
			ten=ten+jc[4]; // IFE_s-l * aniso
//			ten=ten+jc[4]*f_IFE; // IFE_s-l * aniso
//		else if(grid[i][j][k][1]==PTCLORI)
//			ten=ten+jc[4]*(1-cos(jc[8]*DTOR));
		else if(ori!=now[0]){// || grid[i][j][k][1]!=now[1]){ // Kronoker Delta
			del=misorient(now[0],ori);
			ten=ten+GBE(del,jc[3],jc[1]);
		}
	    }
	}
    }else{	// now is liquid
	for(l=25;l>=0;l--){
		i=px;
		j=py;
		k=pz;
		m=surf_check(l,i,j,k,dir,pbc);
		n=mpi_surf_check(l,i,j,k,grid);
	    if(m>=0){ // surface condition: liquid-vapor 1face
		if(tbc[m][0]<0 || tbc[m][1]<0) // tbc[][0]<0 -> convection // tbc[][1]<0 -> vapor condition
			ten=ten+jc[6]; // surface condition: solid-vapor 1face
		else
			ten=ten+jc[4]*cos(DTOR*tbc[m][1]); // IFE_s-l * cos theta
	    }else if(n>=0){ // surface condition: liquid-vapor 1face
		if(tbc[n][0]<0 || tbc[n][1]<0) // tbc[][0]<0 -> convection // tbc[][1]<0 -> vapor condition
			ten=ten+jc[6]; // surface condition: solid-vapor 1face
		else
			ten=ten+jc[4]*cos(DTOR*tbc[m][1]); // IFE_s-l * cos theta
	    }else{ // ASSUMPTION: PTCL also have the same IFE with l (g_s-l == g_p-l)
		nei26sel(l,&i,&j,&k,dir,pbc);
		if(grid[i][j][k][0]!=LIQORI) // liquid - solid
			ten=ten+jc[4];
		// else: liquid - liquid: no IFE
	    }
	}
    }
	if(now[1]==POWORI)
	    now[0]+=POWORI;
	return ten;
}

double fixed_dt_set(int mcs,int *P_mode,int *fdm_loop,double dtsave[],double dtdMCS[],double dHfus[],double Tinfo[],double factor[],double Vsite,double Asite,int out[],double Tmax[],double pcps[],int**** grid,double**** fdmt,int dir[],int pbc[],double rscale[],double melt[]){
// Acceleration scheme // Tcheck 0 for Tmin, 1 for Tmax ////////////////////////////////////
	double Pmax[2],rtim;
	if(out[5]==0 && mcs>0){
	    find_range(grid,fdmt[0],dir,pbc,Tinfo);
	    if(out[3]<0 && Tinfo[0]==TNOLIQUID)
		return -1;
// P_mode ==0 일때 그냥 FDM, +면 single phase mode인데 여기서 온도까지 똑같으면 isothermal mode, -면 완벽한 isothermal
	    if(fabs(Tinfo[0]-Tinfo[1])<1E-3){ // isothermal condition
	      if(*P_mode>0)
		*P_mode=-1;
	      else
		*P_mode--;
	    }else{
	      if(*P_mode<0){
		printf("  @@@ ISOTHERMAL CONDITION IS NO MORE VALID at %d MCS: FDM ON @@@\n",mcs);
		*P_mode=0;
	      }
	      if(Tinfo[1]==TNOSOLID){ //then no solid in the simulation, so we can ignore
//printf(" @@@ No solid in the simulation... MCS %d @@@",mcs);
		Pmax[0]=0;
/*		if(Tinfo[0]==TNOLIQUID){ // SOLID X, LIQUID X -> impossible
		    Pmax[1]=0;
		}else */if(Tinfo[0]>=Tmax[1]) // SOLID X, T of LIQUID > T*
		    Pmax[1]=0;
		else{ // SOLID X, LIQUID < T*
		    if(Tinfo[0]<=Tmax[0]*KAUZMANN)
			Pmax[1]=1.0;
		    else
			Pmax[1]=P_growth(1/Tinfo[0],pcps,solidi_df(dHfus[0],Tinfo[0],Tmax[0],Tmax[2],melt[8])*Vsite,-1.0);
		}
	      }else if(Tinfo[0]>=Tmax[1]){ //then no liquid that can be solidified in the simulation, so we can ignore
		Pmax[1]=0;		// i.e., SOLID O, LIQUID X or > T*
		if(Tinfo[1]>=Tmax[0])
		    Pmax[0]=1.0;
		else
		    Pmax[0]=P_grain_growth(rscale[3],1/Tmax[0],1/Tinfo[1],Asite,-1.0);
	      }else{ // SOLID O, LIQUID O
		if(Tinfo[1]>=Tmax[0])
		    Pmax[0]=1.0;
		else
		    Pmax[0]=P_grain_growth(rscale[3],1/Tmax[0],1/Tinfo[1],Asite,-1.0);
		if(Tinfo[0]<=Tmax[0]*KAUZMANN)
		    Pmax[1]=1.0;
		else
		    Pmax[1]=P_growth(1/Tinfo[0],pcps,solidi_df(dHfus[0],Tinfo[0],Tmax[0],Tmax[2],melt[8])*Vsite,-1.0);
	      }
	      rtim=MCS_time_acceleration(dtdMCS,factor,Pmax,out[5]);
	    }
	}else{ // for fisrt MC step
	    if(out[5]==0)
		rtim=MCS_time_acceleration(dtdMCS,factor,Pmax,1);
	    else
		rtim=MCS_time_acceleration(dtdMCS,factor,Pmax,out[5]);
	}

	if(latent_heat_check(fdmt[2],dir)!=0) // Latent heat is not fully consumed 
	    *P_mode=0;

	if(*P_mode<0){ // isothermal
	  if(*P_mode==-1){ // first step
	    printf("  @@@ PERFECT ISOTHERMAL CONDITION IS ACHIEVED at %d MCS: FDM OFF @@@\n",mcs);
	    *P_mode--;
	    if(Tinfo[0]==TNOLIQUID){
		Pmax[1]=0; // No liquid
		if(Tinfo[0]>=Tmax[0])
		    Pmax[0]=1.0;
		else
		    Pmax[0]=P_grain_growth(rscale[3],1/Tmax[0],1/Tinfo[0],Asite,-1.0);
		rtim=MCS_time_acceleration(dtdMCS,factor,Pmax,out[5]);
	    }else if(Tinfo[0]==TNOSOLID){ // No solid
		Pmax[0]=0; // No solid
		if(Tinfo[0]<=Tmax[0]*KAUZMANN)
		    Pmax[1]=1.0;
		else
		    Pmax[1]=P_growth(1/Tinfo[0],pcps,solidi_df(dHfus[0],Tinfo[0],Tmax[0],Tmax[2],melt[8])*Vsite,-1.0);
		rtim=MCS_time_acceleration(dtdMCS,factor,Pmax,out[5]);
	    }
	    dtsave[1]=rtim;
	    melt[1]=dtsave[0];
	    if(rtim<melt[1]){ // fdm time > MC time
		melt[1]=rtim;
		*fdm_loop=0;
		// in FDM loop, automatically ++
	    }else
		*fdm_loop=(int)(rtim/melt[1]);
		melt[6]=dHfus[1]*melt[1]/(rscale[1]*rscale[1])*melt[3]*rscale[9]/melt[4];
		if(melt[6]>dHfus[1])
		    melt[6]=dHfus[1];
	    }
	}else if(rtim<=rscale[12] && rtim<=rscale[13]*dtsave[1]){
	    if(dtsave[1]!=rtim){
		melt[1]=dtsave[0];
		if(rtim<melt[1]){ // fdm time > MC time
		    melt[1]=rtim;
		    *fdm_loop=0;
		// in FDM loop, automatically ++
		}else
		    *fdm_loop=(int)(rtim/melt[1]);
		melt[6]=dHfus[1]*melt[1]/(rscale[1]*rscale[1])*melt[3]*rscale[9]/melt[4];
		if(melt[6]>dHfus[1])
		    melt[6]=dHfus[1];
		dtsave[1]=rtim;
	    }
	}else{
	    if(rtim>rscale[12]){
		factor[0]=factor[0]*rscale[12]/rtim;
		factor[1]=factor[1]*rscale[12]/rtim;
		rtim=rscale[12];
	    }
	    if(rtim>rscale[13]*dtsave[1]){
		factor[0]=factor[0]*rscale[13]*dtsave[1]/rtim;
		factor[1]=factor[1]*rscale[13]*dtsave[1]/rtim;
		rtim=rscale[13]*dtsave[1];
	    }
	    if(rtim<melt[1]){ // fdm time > MC time
		melt[1]=rtim;
		*fdm_loop=0;
		// in FDM loop, automatically ++
	    }else
		*fdm_loop=(int)(rtim/melt[1]);

	    dtsave[1]=rtim;
	}
//printf("mcs %d rtim %E fdm loop %d\n",mcs,rtim,*fdm_loop);
	return rtim;
}

double active_dt_set(int mcs,int* P_mode,int *fdm_loop,double dtsave[],double dtdMCS[],double dHfus[],double Tinfo[],double factor[],double Vsite,double Asite,int out[],double Tmax[],double pcps[],int**** grid,double**** fdmt,int dir[],int pbc[],double rscale[],double melt[]){
// Acceleration scheme // Tcheck 0 for Tmin, 1 for Tmax ////////////////////////////////////
	double Pmax[2],rtim;
	if(out[5]==0 && mcs>0){
	    find_range(grid,fdmt[0],dir,pbc,Tinfo);
	    if(out[3]<0 && Tinfo[0]==TNOLIQUID)
		return -1;
// P_mode ==0 일때 그냥 FDM, +면 single phase mode인데 여기서 온도까지 똑같으면 isothermal mode, -면 완벽한 isothermal
	    if(fabs(Tinfo[0]-Tinfo[1])<1E-3){ // isothermal condition
	      if(*P_mode>0)
		*P_mode=-1;
	      else
		*P_mode--;
	    }else{
	      if(*P_mode<0){
		printf("  @@@ ISOTHERMAL CONDITION IS NO MORE VALID at %d MCS: FDM ON @@@\n",mcs);
		*P_mode=0;
	      }
	      if(Tinfo[1]==TNOSOLID){ //then no solid in the simulation, so we can ignore
		Pmax[0]=0;
		*P_mode++;
/*			if(Tinfo[0]==TNOLIQUID){ // SOLID X, LIQUID X -> impossible
			    Pmax[1]=0;
			    *P_mode++;
		}else */if(Tinfo[0]>=Tmax[1]) // SOLID X, T of LIQUID > T*
		    Pmax[1]=0;
		else{ // SOLID X, LIQUID < T*
		    if(Tinfo[0]<=Tmax[0]*KAUZMANN)
			Pmax[1]=1.0;
		    else
			Pmax[1]=P_growth(1/Tinfo[0],pcps,solidi_df(dHfus[0],Tinfo[0],Tmax[0],Tmax[2],melt[8])*Vsite,-1.0);
		}
	      }else if(Tinfo[0]>=Tmax[1]){ //then no liquid that can be solidified in the simulation, so we can ignore
		Pmax[1]=0;		// i.e., SOLID O, LIQUID X or > T*
		if(Tinfo[0]==TNOLIQUID) // SOLID O, LIQUID X
		    *P_mode++;
			if(Tinfo[1]>=Tmax[0])
		    Pmax[0]=1.0;
		else
		    Pmax[0]=P_grain_growth(rscale[3],1/Tmax[0],1/Tinfo[1],Asite,-1.0);
	      }else{ // SOLID O, LIQUID O
	        *P_mode=0;
		if(Tinfo[1]>=Tmax[0])
		    Pmax[0]=1.0;
		else
		    Pmax[0]=P_grain_growth(rscale[3],1/Tmax[0],1/Tinfo[1],Asite,-1.0);
			if(Tinfo[0]<=Tmax[0]*KAUZMANN)
		    Pmax[1]=1.0;
		else
		    Pmax[1]=P_growth(1/Tinfo[0],pcps,solidi_df(dHfus[0],Tinfo[0],Tmax[0],Tmax[2],melt[8])*Vsite,-1.0);
	      }
	      rtim=MCS_time_acceleration(dtdMCS,factor,Pmax,out[5]);
	    }
	}else{ // for first MC step
	    if(out[5]==0)
		rtim=MCS_time_acceleration(dtdMCS,factor,Pmax,1);
	    else
		rtim=MCS_time_acceleration(dtdMCS,factor,Pmax,out[5]);
	}
	if(latent_heat_check(fdmt[2],dir)!=0) // Latent heat is not fully consumed 
	    *P_mode=0;

	if(*P_mode<0){ // isothermal
	  if(*P_mode==-1){ // first step
	    printf("  @@@ PERFECT ISOTHERMAL CONDITION IS ACHIEVED at %d MCS: FDM OFF @@@\n",mcs);
	    *P_mode--;
	    if(Tinfo[0]==TNOLIQUID){
		Pmax[1]=0; // No liquid
		if(Tinfo[0]>=Tmax[0])
		    Pmax[0]=1.0;
		else
		    Pmax[0]=P_grain_growth(rscale[3],1/Tmax[0],1/Tinfo[0],Asite,-1.0);
		rtim=MCS_time_acceleration(dtdMCS,factor,Pmax,out[5]);
	    }else if(Tinfo[0]==TNOSOLID){ // No solid
		Pmax[0]=0; // No solid
		if(Tinfo[0]<=Tmax[0]*KAUZMANN)
		    Pmax[1]=1.0;
		else
		    Pmax[1]=P_growth(1/Tinfo[0],pcps,solidi_df(dHfus[0],Tinfo[0],Tmax[0],Tmax[2],melt[8])*Vsite,-1.0);
		rtim=MCS_time_acceleration(dtdMCS,factor,Pmax,out[5]);
	    }
	    dtsave[1]=rtim;
	    melt[1]=dtsave[0];
	    if(rtim<melt[1]){ // fdm time > MC time
		melt[1]=rtim;
		*fdm_loop=0;
		// in FDM loop, automatically ++
	    }else
		*fdm_loop=(int)(rtim/melt[1]);
	    melt[6]=dHfus[1]*melt[1]/(rscale[1]*rscale[1])*melt[3]*rscale[9]/melt[4];
	    if(melt[6]>dHfus[1])
		melt[6]=dHfus[1];
	  }//else{ // isothermal condition remains
//		  }
	}else{
	  if(rtim<=rscale[12] && rtim<=rscale[13]*dtsave[1]){
	    if(dtsave[1]!=rtim){
		if(*P_mode!=0) // Acceleration for FDM
		    melt[1]=dtsave[0]*melt[7]*melt[7];
		else
		    melt[1]=dtsave[0];
			if(rtim<melt[1]){ // fdm time > MC time
		    melt[1]=rtim;
		    *fdm_loop=0;
		// in FDM loop, automatically ++
		}else
		    *fdm_loop=(int)(rtim/melt[1]);
			if(*P_mode!=0) // Acceleration for FDM
		    melt[6]=dHfus[1]*melt[1]/(rscale[1]*rscale[1]*melt[7]*melt[7])*melt[3]*rscale[9]/melt[4];
		else
		    melt[6]=dHfus[1]*melt[1]/(rscale[1]*rscale[1])*melt[3]*rscale[9]/melt[4];
		if(melt[6]>dHfus[1])
		    melt[6]=dHfus[1];
		dtsave[1]=rtim;
	    }
	  }else{
	    if(rtim>rscale[12]){
		factor[0]=factor[0]*rscale[12]/rtim;
		factor[1]=factor[1]*rscale[12]/rtim;
		rtim=rscale[12];
	    }
	    if(rtim>rscale[13]*dtsave[1]){
		factor[0]=factor[0]*rscale[13]*dtsave[1]/rtim;
		factor[1]=factor[1]*rscale[13]*dtsave[1]/rtim;
		rtim=rscale[13]*dtsave[1];
	    }
	    if(rtim<melt[1]){ // fdm time > MC time
		melt[1]=rtim;
		*fdm_loop=0;
		// in FDM loop, automatically ++
	    }else
		*fdm_loop=(int)(rtim/melt[1]);
		dtsave[1]=rtim;
	    }
	}
	return rtim;
}

double P_grain_growth(double Q,double revTm,double revT,double dx2,double dG){
return (revT<revTm) ? ((dG<=0) ? 1.0 : exp(-dx2*dG*REVkB*revT)) :\
		      ((dG<=0) ? exp(Q*REVRGAS*(revTm-revT)) :exp(Q*REVRGAS*(revTm-revT)-dx2*dG*REVkB*revT));
}

double P_nucl(double revT,double dGfus,double gam,double factor,double Nc,double dG,double pcps[],double ps){
	double temp;
	temp=exp(log(Nc)-revT*REVkB*factor*gam*gam*gam/(dGfus*dGfus))-1.0;
	return (temp>=0) ? P_growth(revT,pcps,ps,dG) : (temp+1.0)*P_growth(revT,pcps,ps,dG);
}

double P_growth(double revT,double pcps[],double ps,double dG){
    return (dG<=0) ? (pcps[0]+ps)*pcps[2] : exp(-dG*REVkB*revT)*pcps[0]*pcps[2];
}

double solidi_df(double dHfus,double Tnow,double Tliq,double Tfreez,double T0){
	double Tsol;
	if(Tfreez==0) // for pure element, no freezing range
	    return dHfus*(1-Tnow/Tliq);
	else{
	    Tsol=Tliq-Tfreez;
	    if(Tnow<Tsol)
		return dHfus*(1-Tnow/T0); // assume linear proportional to T, where df=0 at T0
	    else
		return dHfus*(1-Tsol/T0)*(Tliq-Tnow)/Tfreez;
	}
}

void heat_transfer(double**** tgrid,double tbc[][2],int dir[],int pbc[],double mcdx,double const melt[],int fdm_loop,double k_per_cp,double Tmelt,int mode){ //,double dHfus){
	int i,j,k,l,inter=1,temp[3],fdmax[3];//imax,jmax,kmax;
	double hconv,revdx2;
	
	fdm_loop++; // considering round down when (double -> int)
	if(fdm_loop==0){
		printf("FDM ERROR: wrong FDM loop number %d\n  -> Automatically changed to 1\n");
		fdm_loop=1;
	}

	fdmax[0]=dir[0]-1;
	fdmax[1]=dir[1]-1;
	fdmax[2]=dir[2]-1;

	if(mode==0){
	    revdx2=1/(mcdx*mcdx);
	    inter=1;
	}else{
	    mcdx=mcdx*melt[7];
	    revdx2=1/(mcdx*mcdx);
	    inter=(int)melt[7];
	    if(mode==1){
// Take T values as the average of surrounding values
	      for(i=fdmax[0];i>=0;i-=inter){
		for(j=fdmax[1];j>=0;j-=inter){
		    for(k=fdmax[2];k>=0;k-=inter)
	 		tgrid[0][i][j][k]=average_temp(tgrid[0],i,j,k,inter,dir,pbc);
		}
	      }
// Indicate voxels that are not included in the calculation
	      for(i=fdmax[0];i>=0;i--){
		for(j=fdmax[1];j>=0;j--){
		    for(k=fdmax[2];k>=0;k--)
			tgrid[1][i][j][k]=-1;
		}
	      }
	    }
	}

// Latent heat divided by loop number : retroactive accumulation of latent heat for 1 MCS
// Latent heat multiplied with phase transf. fraction : when dx_fdm != dx_mc, delta_fs != 1: delta_fs == V_mc / V_fdm
// FDM loop
	for(l=fdm_loop;l>0;l--){ // FDM time flow to 1 MCS
//printf("fdm_loop %d l %d\n",fdm_loop,l);
	    if(mode==0){ // if P_mode>0, single phase mode w/o accumulated latent heat
// Latent heat consideration
		for(i=fdmax[0];i>=0;i--){
		    for(j=fdmax[1];j>=0;j--){
			for(k=fdmax[2];k>=0;k--){
			    if(tgrid[2][i][j][k]!=0)
				LatHeatDiffus(&(tgrid[0][i][j][k]),&(tgrid[2][i][j][k]),melt[6]);
			}
		    }
		}
	    }
//	0. inner part
		for(i=fdmax[0]-inter;i>0;i-=inter){
			for(j=fdmax[1]-inter;j>0;j-=inter){
				for(k=fdmax[2]-inter;k>0;k-=inter)
				  tgrid[1][i][j][k]=tgrid[0][i][j][k]+melt[1]*k_per_cp*\
						   (revdx2*(tgrid[0][i+inter][j][k]+tgrid[0][i-inter][j][k]-2*tgrid[0][i][j][k])+\
						    revdx2*(tgrid[0][i][j+inter][k]+tgrid[0][i][j-inter][k]-2*tgrid[0][i][j][k])+\
						    revdx2*(tgrid[0][i][j][k+inter]+tgrid[0][i][j][k-inter]-2*tgrid[0][i][j][k]));
			}	
		}

// 1. X direction
		for(i=fdmax[0]-inter;i>0;i-=inter){
			j=fdmax[1];
			for(k=fdmax[2]-inter;k>0;k-=inter){ // bot & top xy planes
				tgrid[1][i][0][k]=tgrid[0][i][0][k]+melt[1]*k_per_cp*\
					   revdx2*(tgrid[0][i+inter][0][k]+tgrid[0][i-inter][0][k]-2*tgrid[0][i][0][k]);
				tgrid[1][i][j][k]=tgrid[0][i][j][k]+melt[1]*k_per_cp*\
					   revdx2*(tgrid[0][i+inter][j][k]+tgrid[0][i-inter][j][k]-2*tgrid[0][i][j][k]);
			}
			k=fdmax[2];
			for(j=fdmax[1];j>=0;j-=inter){ // bot & top xz planes, including edges
				tgrid[1][i][j][0]=tgrid[0][i][j][0]+melt[1]*k_per_cp*\
					   revdx2*(tgrid[0][i+inter][j][0]+tgrid[0][i-inter][j][0]-2*tgrid[0][i][j][0]);
				tgrid[1][i][j][k]=tgrid[0][i][j][k]+melt[1]*k_per_cp*\
					   revdx2*(tgrid[0][i+inter][j][k]+tgrid[0][i-inter][j][k]-2*tgrid[0][i][j][k]);
			}
		}
		i=fdmax[0];
 		if(pbc[0]!=0){// X pbc on
			for(j=fdmax[1];j>=0;j-=inter){
				for(k=fdmax[2];k>=0;k-=inter){ // [0][j][k] = [fdmax[0]][j][k]
					tgrid[1][0][j][k]=tgrid[0][0][j][k]+melt[1]*k_per_cp*revdx2*\
							  (tgrid[0][inter][j][k]+tgrid[0][i][j][k]-2*tgrid[0][0][j][k]);
					tgrid[1][i][j][k]=tgrid[0][i][j][k]+melt[1]*k_per_cp*revdx2*\
							  (tgrid[0][0][j][k]+tgrid[0][i-inter][j][k]-2*tgrid[0][i][j][k]);
				}
			}
		}else{ // X pbc off
// Lower X plane
			if(tbc[0][0]==0){ // adiabatic
			    for(j=fdmax[1];j>=0;j-=inter){
				for(k=fdmax[2];k>=0;k-=inter)
				  tgrid[1][0][j][k]=tgrid[0][0][j][k]+melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][inter][j][k]-tgrid[0][0][j][k]);
			    }
			}else if(tbc[0][0]<0){ // convection
			    hconv=tbc[0][1]*k_per_cp*melt[1]/(mcdx*melt[3]); // h k dt / (cp dx k) = h dt / (cp dx)
			    for(j=fdmax[1];j>=0;j-=inter){
				for(k=fdmax[2];k>=0;k-=inter)
				  tgrid[1][0][j][k]=tgrid[0][0][j][k]+melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][inter][j][k]-tgrid[0][0][j][k])-\
						    hconv*(tbc[0][0]+tgrid[0][0][j][k]); // - (-Ta+Tnow)
			    }
			}else{ // tbc > 0 : heat sink
			    for(j=fdmax[1];j>=0;j-=inter){
				for(k=fdmax[2];k>=0;k-=inter)
				  tgrid[1][0][j][k]=tgrid[0][0][j][k]+melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][inter][j][k]+tbc[0][0]-2*tgrid[0][0][j][k]);
			    }
			}
// Upper X plane
			if(tbc[1][0]==0){ // adiabatic
			    for(j=fdmax[1];j>=0;j-=inter){
				for(k=fdmax[2];k>=0;k-=inter)
				  tgrid[1][i][j][k]=tgrid[0][i][j][k]+melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i-inter][j][k]-tgrid[0][i][j][k]);
			    }
			}else if(tbc[1][0]<0){ // convection
			    hconv=tbc[1][1]*k_per_cp*melt[1]/(mcdx*melt[3]); // h k dt / (cp dx k) = h dt / (cp dx)
			    for(j=fdmax[1];j>=0;j-=inter){
				for(k=fdmax[2];k>=0;k-=inter)
				  tgrid[1][i][j][k]=tgrid[0][i][j][k]+melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i-inter][j][k]-tgrid[0][i][j][k])-\
						    hconv*(tbc[1][0]+tgrid[0][i][j][k]); // - (-Ta+Tnow)
			    }
			}else{ // tbc > 0 : heat sink
			    for(j=fdmax[1];j>=0;j-=inter){
				for(k=fdmax[2];k>=0;k-=inter)
				  tgrid[1][i][j][k]=tgrid[0][i][j][k]+melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i-inter][j][k]+tbc[1][0]-2*tgrid[0][i][j][k]);
			    }
			}
		}

// 2. Y direction
		for(j=fdmax[1]-inter;j>0;j-=inter){
			i=fdmax[0];
			for(k=fdmax[2]-inter;k>0;k-=inter){ // bot & top xz planes
				tgrid[1][0][j][k]+=melt[1]*k_per_cp*\
					   revdx2*(tgrid[0][0][j+inter][k]+tgrid[0][0][j-inter][k]-2*tgrid[0][0][j][k]);
				tgrid[1][i][j][k]+=melt[1]*k_per_cp*\
					   revdx2*(tgrid[0][i][j+inter][k]+tgrid[0][i][j-inter][k]-2*tgrid[0][i][j][k]);
			}
			k=fdmax[2];
			for(i=fdmax[0];i>=0;i-=inter){ // bot & top xz planes, including edges
				tgrid[1][i][j][0]+=melt[1]*k_per_cp*\
					   revdx2*(tgrid[0][i][j+inter][0]+tgrid[0][i][j-inter][0]-2*tgrid[0][i][j][0]);
				tgrid[1][i][j][k]+=melt[1]*k_per_cp*\
					   revdx2*(tgrid[0][i][j+inter][k]+tgrid[0][i][j-inter][k]-2*tgrid[0][i][j][k]);
			}
		}

		j=fdmax[1];
 		if(pbc[1]!=0){// Y pbc on
			for(i=fdmax[0];i>=0;i-=inter){ // [i][fdmax][k] = [i][0][k]
				for(k=fdmax[2];k>=0;k-=inter){
					tgrid[1][i][0][k]+=melt[1]*k_per_cp*revdx2*\
							  (tgrid[0][i][inter][k]+tgrid[0][i][j][k]-2*tgrid[0][i][0][k]);
					tgrid[1][i][j][k]+=melt[1]*k_per_cp*revdx2*\
							  (tgrid[0][i][0][k]+tgrid[0][i][j-inter][k]-2*tgrid[0][i][j][k]);
				}
			}
		}else{
// Lower Y plane
			if(tbc[2][0]==0){ // adiabatic
			    for(i=fdmax[0];i>=0;i-=inter){
				for(k=fdmax[2];k>=0;k-=inter)
				  tgrid[1][i][0][k]+=melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i][inter][k]-tgrid[0][i][0][k]);
			    }
			}else if(tbc[2][0]<0){ // convection
			    hconv=tbc[2][1]*k_per_cp*melt[1]/(mcdx*melt[3]); // h k dt / (cp dx k) = h dt / (cp dx)
			    for(i=fdmax[0];i>=0;i-=inter){
				for(k=fdmax[2];k>=0;k-=inter)
				  tgrid[1][i][0][k]+=melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i][inter][k]-tgrid[0][i][0][k])-\
						    hconv*(tbc[2][0]+tgrid[0][i][0][k]); // - (-Ta+Tnow) = Ta-Tnow
			    }
			}else{ // tbc > 0 : heat sink
			    for(i=fdmax[0];i>=0;i-=inter){
				for(k=fdmax[2];k>=0;k-=inter)
				  tgrid[1][i][0][k]+=melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i][inter][k]+tbc[2][0]-2*tgrid[0][i][0][k]);
			    }
			}
// Upper Y plane
			if(tbc[3][0]==0){ // adiabatic
			    for(i=fdmax[0];i>=0;i-=inter){
				for(k=fdmax[2];k>=0;k-=inter)
				  tgrid[1][i][j][k]+=melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i][j-inter][k]-tgrid[0][i][j][k]);
			    }
			}else if(tbc[3][0]<0){ // convection
			    hconv=tbc[3][1]*k_per_cp*melt[1]/(mcdx*melt[3]); // h k dt / (cp dx k) = h dt / (cp dx)
			    for(i=fdmax[0];i>=0;i-=inter){
				for(k=fdmax[2];k>=0;k-=inter)
				  tgrid[1][i][j][k]+=melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i][j-inter][k]-tgrid[0][i][j][k])-\
						    hconv*(tbc[3][0]+tgrid[0][i][j][k]); // - (-Ta+Tnow)
			    }
			}else{ // tbc > 0 : heat sink
			    for(i=fdmax[0];i>=0;i-=inter){
				for(k=fdmax[2];k>=0;k-=inter)
				  tgrid[1][i][j][k]+=melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i][j-inter][k]+tbc[3][0]-2*tgrid[0][i][j][k]);
			    }
			}
		}

// 3. Z direction
		for(k=fdmax[2]-inter;k>0;k-=inter){
			i=fdmax[0];
			for(j=fdmax[1];j>=0;j-=inter){ // bot & top yz planes
				tgrid[1][0][j][k]+=melt[1]*k_per_cp*\
					   revdx2*(tgrid[0][0][j][k+inter]+tgrid[0][0][j][k-inter]-2*tgrid[0][0][j][k]);
				tgrid[1][i][j][k]+=melt[1]*k_per_cp*\
					   revdx2*(tgrid[0][i][j][k+inter]+tgrid[0][i][j][k-inter]-2*tgrid[0][i][j][k]);
			}
			j=fdmax[1];
			for(i=fdmax[0]-inter;i>0;i-=inter){ // bot & top xz planes
				tgrid[1][i][0][k]+=melt[1]*k_per_cp*\
					   revdx2*(tgrid[0][i][0][k+inter]+tgrid[0][i][0][k-inter]-2*tgrid[0][i][0][k]);
				tgrid[1][i][j][k]+=melt[1]*k_per_cp*\
					   revdx2*(tgrid[0][i][j][k+inter]+tgrid[0][i][j][k-inter]-2*tgrid[0][i][j][k]);
			}
		}

		k=fdmax[2];
 		if(pbc[2]!=0){// Z pbc on
			for(i=fdmax[0];i>=0;i-=inter){
				for(j=fdmax[1];j>=0;j-=inter){ // [i][k][0] = [i][j][fdmax]
					tgrid[1][i][j][0]+=melt[1]*k_per_cp*revdx2*\
							  (tgrid[0][i][j][inter]+tgrid[0][i][j][k]-2*tgrid[0][i][j][0]);
					tgrid[1][i][j][k]+=melt[1]*k_per_cp*revdx2*\
							  (tgrid[0][i][j][0]+tgrid[0][i][j][k-inter]-2*tgrid[0][i][j][k]);
				}
			}
		}else{
// Lower Z plane
			if(tbc[4][0]==0){ // adiabatic
			    for(i=fdmax[0];i>=0;i-=inter){
				for(j=fdmax[1];j>=0;j-=inter)
				  tgrid[1][i][j][0]+=melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i][j][inter]-tgrid[0][i][j][0]);
			    }
			}else if(tbc[4][0]<0){ // convection
			    hconv=tbc[4][1]*k_per_cp*melt[1]/(mcdx*melt[3]); // h k dt / (cp dx k) = h dt / (cp dx)
			    for(i=fdmax[0];i>=0;i-=inter){
				for(j=fdmax[1];j>=0;j-=inter)
				  tgrid[1][i][j][0]+=melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i][j][inter]-tgrid[0][i][j][0])-\
						    hconv*(tbc[4][0]+tgrid[0][i][j][0]); // - (-Ta+Tnow) = Ta-Tnow
			    }
			}else{ // tbc > 0 : heat sink
			    for(i=fdmax[0];i>=0;i-=inter){
				for(j=fdmax[1];j>=0;j-=inter)
				  tgrid[1][i][j][0]+=melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i][j][inter]+tbc[4][0]-2*tgrid[0][i][j][0]);
			    }
			}
// Upper Z plane
			if(tbc[5][0]==0){ // adiabatic
			    for(i=fdmax[0];i>=0;i-=inter){
				for(j=fdmax[1];j>=0;j-=inter)
				  tgrid[1][i][j][k]+=melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i][j][k-inter]-tgrid[0][i][j][k]);
			    }
			}else if(tbc[5][0]<0){ // convection
			    hconv=tbc[5][1]*k_per_cp*melt[1]/(mcdx*melt[3]); // h k dt / (cp dx k) = h dt / (cp dx)
//printf("hconv %E for time %E if melt 1 is 1 then hconv %E\n",hconv,melt[1],hconv/melt[1]);
			    for(i=fdmax[0];i>=0;i-=inter){
				for(j=fdmax[1];j>=0;j-=inter)
//printf("T %.f -> ",tgrid[1][i][j][k]);
				  tgrid[1][i][j][k]+=melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i][j][k-inter]-tgrid[0][i][j][k])-\
						    hconv*(tbc[5][0]+tgrid[0][i][j][k]); // - (-Ta+Tnow)
//printf("T %.f -> ",tgrid[1][i][j][k]);
			    }
			}else{ // tbc > 0 : heat sink
			    for(i=fdmax[0];i>=0;i-=inter){
				for(j=fdmax[1];j>=0;j-=inter)
				  tgrid[1][i][j][k]+=melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i][j][k-inter]+tbc[5][0]-2*tgrid[0][i][j][k]);
			    }
			}
		}
// next -> now
		for(i=fdmax[0];i>=0;i-=inter){
			for(j=fdmax[1];j>=0;j-=inter){
				for(k=fdmax[2];k>=0;k-=inter)
					tgrid[0][i][j][k]=tgrid[1][i][j][k]; // dHfus is dH_fus / Cp
		}	}
	} // FDM 1 iteration
	if(mode>0){ // if P_mode>0, single phase mode w/o accumulated latent heat
	// Take averaged values for other voxels that are not included in the calculation
	    temp[0]=fdmax[0];
	    for(i=fdmax[0];i>=0;i--){
		if(i<temp[0]-inter)
		    temp[0]-=inter;
		temp[1]=fdmax[1];

		for(j=fdmax[1];j>=0;j--){
		    if(j<temp[1]-inter)
			temp[1]-=inter;
		    temp[2]=fdmax[2];

		    for(k=fdmax[2];k>=0;k--){
			if(k<temp[2]-inter)
			    temp[2]-=inter;
			if(tgrid[1][i][j][k]<0)
			    tgrid[0][i][j][k]=recover_temp(tgrid[0],i,j,k,inter,temp);
		    }
		}
	    }
	}

	return;
}

void MeltPool_movement(double pos[],double direction,double displ){
	if(direction<0)
	    pos[0]+=displ;
	else
	    pos[1]+=displ;
	return;
}

void MeltPool(int**** grid,double**** tgrid,double am[],double mp[],int dir[],int pbc[],double Tm,double dH,double pos[]){
	int i,j,k,xmax,ymax,zmax,tear[2],mode=0;
	double temp,x,y;

//printf("position x %f y %f // rank %d\n",pos[0],pos[1],n_rank);
//MPI_Barrier(MPI_COMM_WORLD);
	if(am[0]>0){//GAUSSIAN // only consider 2 sigma range
	    xmax=deter((int)(pos[0]+am[2]),dir[0],0);
	    x=deter((int)(pos[0]-am[2]),dir[0],0);
	    ymax=deter((int)(pos[1]+am[3]),dir[1],0);
	    y=deter((int)(pos[1]-am[3]),dir[1],0);

/*	  if(mp[0]==0){ // 1st mcs // Arbitrary asume a primarily formed melt pool, for comput. efficiency
	    k=dir[2]-1;
	    for(i=xmax;i>=x;i--){
		temp=am[3]*sqrt(1-(double)((i-pos[0])*(i-pos[0]))/(am[2]*am[2]));
		for(j=deter((int)(pos[1]+temp),dir[1],0);j>=deter((int)(pos[1]-temp),dir[1],0);j--){
		    tgrid[0][i][j][k]=Tm;
		    if(grid[i][j][k][0]!=LIQORI){ // LIQ melt at FDM stage?
			if(grid[i][j][k][0]!=MPISURF && grid[i][j][k][0]!=PTCLORI)
			    neworient_melting(grid[i][j][k]);
//			    tgrid[2][i][j][k]-=dH; //dH == dHfus[1]
		    } //
		}
	    }
	  }*/
	    for(i=xmax;i>=x;i--){
		for(j=ymax;j>=y;j--){
		    // 2D gaussian
		    tgrid[0][i][j][dir[2]-1]+=mp[1]*exp(-0.5*((i-am[2])*(i-am[2])*mp[2]+(j-am[3])*(j-am[3])*mp[3]));
	/*		    for(k=dir[2]-1;k>=0;k--){ // 3D gaussian
	tgrid[0][i][j][k]+=mp[1]*exp(-0.5*((i-am[2])*(i-am[2])*mp[2]+(j-am[3])*(j-am[3])*mp[3]+(dir[2]-1-k)*(dir[2]-1-k)*mp[4]));
		    }*/
		}
	    }
/*	    for(i=xmax;i>=x;i--){
		for(j=ymax;j>=y;j--){
		    for(k=dir[2]-1;k>=0;k--){
			if(tgrid[0][i][j][k]>=Tm && grid[i][j][k][0]!=LIQORI){
			    if(grid[i][j][k][0]!=MPISURF && grid[i][j][k][0]!=PTCLORI)
				neworient_melting(grid[i][j][k]);
//				tgrid[2][i][j][k]-=dH; //dH == dHfus[1]
			} //
		    }
		}
	    }*/
	}else if(am[0]==0){ // Teardrop shaped Melt pool
	    if(am[9]<Tm) // MP boundary temperature < Tmax[0]
		mode=1;
//W	am[2]	mp[2]=am[2]/(2 sin amxt (sin maxt/2)^am[5])
//C+T L	am[3]	mp[1]=am[3]/2
//dep	am[4]	mp[3]=am[4]/(sin 4atan pow(sin 4atan/2)^am[5])
	    if(am[1]<0){ // +x direction
		xmax=deter((int)(pos[0]+mp[1]),dir[0],0);
		for(i=deter((int)(pos[0]-mp[1]),dir[0],0);i<xmax;i++){ // Cap+Tail Melt Pool
		    if(pos[0]-(double)i<=mp[1])
			temp=acos((pos[0]-(double)i)/mp[1]); // t
		    else
			temp=0; // t=0

		    tear[0]=teardrop(mp[2],temp,am[5]); // y value on the curve when x=i & o=(0,0)
		    ymax=deter((int)(pos[1]+tear[0]),dir[1],0);
		    tear[1]=teardrop(mp[3],temp,am[5]); // z value on the curve: dept of melt pool
		    for(j=deter(pos[1]-tear[0],dir[1],0);j<=ymax;j++){
			if(tear[0]!=0){
			    zmax=(int)(tear[1]*(1-(double)((j-pos[1])*(j-pos[1]))/(tear[0]*tear[0]))); // Bezier for Melt Pool
			    for(k=deter(dir[2]-zmax,dir[2],0);k<dir[2];k++){ // Melt Pool
				tgrid[0][i][j][k]=am[9];
				if(mode==0 && grid[i][j][k][0]!=LIQORI){ // LIQ melt at FDM stage?
				    if(grid[i][j][k][0]!=MPISURF && grid[i][j][k][0]!=PTCLORI)
					neworient_melting(grid[i][j][k]);
//				    tgrid[2][i][j][k]-=dH; //dH == dHfus[1]
				} //
			    }
			}
		    }
		}
	    }else{ // am[1]>0 +y direction
		xmax=deter((int)(pos[1]+mp[1]),dir[1],0);
		for(i=deter((int)(pos[1]-mp[1]),dir[1],0);i<xmax;i++){ // Cap+Tail Melt Pool
		    temp=acos((pos[1]-(double)i)/mp[1]); // t
		    tear[0]=teardrop(mp[2],temp,am[5]); // y value on the curve when x=i & o=(0,0)
		    ymax=deter((int)(pos[0]+tear[0]),dir[0],0);
		    tear[1]=teardrop(mp[3],temp,am[5]); // z value on the curve: dept of melt pool
		    for(j=deter(pos[0]-tear[0],dir[0],0);j<=ymax;j++){
			if(tear[0]!=0){
			    zmax=(int)(tear[1]*(1-(double)((j-pos[0])*(j-pos[0]))/(tear[0]*tear[0]))); // Bezier for Melt Pool
			    for(k=deter(dir[2]-zmax,dir[2],0);k<dir[2];k++){ // Melt Pool
				tgrid[0][j][i][k]=am[9];
				if(mode==0 && grid[j][i][k][0]!=LIQORI){ // LIQ melt at FDM stage?
				    if(grid[j][i][k][0]!=MPISURF && grid[j][i][k][0]!=PTCLORI)
					neworient_melting(grid[i][j][k]);
//				    tgrid[2][i][j][k]-=dH; //dH == dHfus[1]
				} //
			    }
			}
		    }
		}
	    }
	}else{ //if(am[0]<0){ // Circle-shaped melt pool
	    if(am[9]<Tm) // MP boundary temperature < Tmax[0]
		mode=1;
	    zmax=dir[2]-1-(int)am[4];
	    if(zmax<0)
		zmax=0;

	    if(am[1]<0){ // +x direction am[6] -> mp6
		x=round(pos[0]);
		y=pos[1];
	    }else{ // +y direction am[7] -> mp6
		x=pos[0];
		y=round(pos[1]);
	    }

	    for(k=dir[2]-1;k>zmax;k--){
		temp=sqrt(am[2]*am[2]*(1-(double)(dir[2]-1-k)/am[4])); // radius at each z coordinate
		xmax=deter((int)(x+(int)temp),dir[0],0);
		for(i=deter((int)(x-(int)temp),dir[0],0);i<xmax;i++){ // Cap+Tail Melt Pool
		    ymax=deter((int)(y+sqrt(temp*temp-(i-x)*(i-x))),dir[1],0);
		    for(j=deter((int)(y-sqrt(temp*temp-(i-x)*(i-x))),dir[1],0);j<ymax;j++){
			tgrid[0][i][j][k]=am[9];
			if(mode==0 && grid[i][j][k][0]!=LIQORI){ // LIQ melt at FDM stage?
			    if(grid[i][j][k][0]!=MPISURF && grid[i][j][k][0]!=PTCLORI)
				neworient_melting(grid[i][j][k]);
//			    tgrid[2][i][j][k]-=dH; //dH == dHfus[1]
			} // */
		    }
		}
	    }
	}
	return;
}

int teardrop(double b, double t,int geo){  // Melt pool shape curve
return (int)(b*sin(t)*(pow(sin(t/2),geo)));
}

void check_melting(int**** grid,double**** tgrid,double Tm,int dir[]){
	int i,j,k;
	for(i=dir[0]-1;i>=0;i--){
	    for(j=dir[1];j>=0;j--){
		for(k=dir[2];k>=0;k--){
		    if(tgrid[0][i][j][k]>=Tm && grid[i][j][k][0]!=LIQORI){
			if(grid[i][j][k][0]!=MPISURF && grid[i][j][k][0]!=PTCLORI)
			    neworient_melting(grid[i][j][k]);
		    }
		}
	    }
	}

	return;
}

void Particle_Form(int**** grid,int dir[],int pbc[],double fraction){
	int i,j,k,l,m,n,o,liq=0,p,imax,imin;
	double ptcl=0;

	if(n_size==1){
	    imax=dir[0]-1;
	    imin=0;
	}else{ // MPI
	    imax=dir[0];
	    imin=1;
	}

	for(i=imax;i>=imin;i--){
	    for(j=dir[1]-1;j>=0;j--){
		for(k=dir[2]-1;k>=0;k--){
		    if(grid[i][j][k][0]==LIQORI) //LIQUID
			liq++;
		    else if(grid[i][j][k][0]==PTCLORI){
			if(nucl_cond(grid,i,j,k,dir,pbc,0)==PTCLORI){
			    p=0;
			    for(o=25;o>=0;o--){
				l=i;
				m=j;
				n=k;
				nei26sel(o,&l,&m,&n,dir,pbc);
				if(grid[l][m][n][0]==LIQORI) // particle in liquid
				    p++;
			    }
			    if(p!=0)
				ptcl++;
			}
		    }
		}
	    }
	}
// particle introduction
	if(liq!=0){
	  ptcl=fraction*((double)liq+ptcl)-ptcl;
	  while(ptcl>0){
		i=(int)(rand()/(RAND_MAX/(imax-imin+1)+1))+imin;
		j=(int)(rand()/(RAND_MAX/dir[1]+1));
		k=(int)(rand()/(RAND_MAX/dir[2]+1));
		if(grid[i][j][k][0]==LIQORI){ //LIQUID
		    if(ptcl>=1.0){
			grid[i][j][k][0]=PTCLORI;
			grid[i][j][k][1]=PTCLORI;
			grid[i][j][k][2]=PTCLORI;
			ptcl--;
		    }else{
			if(ptcl>=((double)rand()/(RAND_MAX))){
			    grid[i][j][k][0]=PTCLORI;
			    grid[i][j][k][1]=PTCLORI;
			    grid[i][j][k][2]=PTCLORI;
			}
			ptcl--;
		    }
		}
	  }
	}
	return;
}

void LatHeatDiffus(double* T,double *L,double dT){
	if(*L>0){ // Latent heat release
	    if((*L)*(*L-dT)<0){ // due to changing L value
		*T+=*L;
		*L=0;
	    }else{
		*T+=dT;
		*L-=dT;
		if(fabs(*L)<PROBMIN)
		    *L=0;
	    }
	}else{ // Latent heat absorption
	    if((*L)*(*L-dT)<0){ // due to changing L value
		*T-=*L;
		*L=0;
	    }else{
		*T-=dT;
		*L+=dT;
		if(fabs(*L)<PROBMIN)
		    *L=0;
	    }
	}
	return;
}

double check_growth_dir(int**** grid,int l,int m,int n,int vector,int dir[],int pbc[],double factor){
	int i,j,k;
	i=l;
	j=m;
	k=n;
	nei26sel_nodeter(vector,&i,&j,&k,dir,pbc);
	if(n_size==1){
	    if(i<0||i>=dir[0])
		if(pbc[0]==0) // surface
		    return factor;
		else
		    i=deter(i,dir[0],pbc[0]);
	}else{ // MPI
	    if(i<0||i>=dir[0]){
	    // Boundary 1 voxel이라 그 이상은 알수가 없음 -> 그냥 랜덤하게 처리
	    // preferred direction을 확실히 정하는 방법에서는 에러가 없음
		if((double)rand()/(RAND_MAX)<0.5)
		    return factor;
		else
		    return 1.0;
	    }
	    if(n_rank==0){
		if(i>=dir[0]-1){
		    if(pbc[0]==0) // surface
			return factor;
		}
	    }else if(n_rank==n_size-1){
		if(i<=0){
		    if(pbc[0]==0)// surface
			return factor;
		}
	    }
	}
	if((j<0||j>=dir[1])){
	    if(dir[1]==0)// surface
		return factor;
	    else
		j=deter(j,dir[1],pbc[1]);
	}
	if((k<0||k>=dir[2])){
	    if(dir[2]==0)// surface
		return factor;
	    else
		k=deter(k,dir[2],pbc[2]);
	}
	if(grid[l][m][n][0]==grid[i][j][k][0])
	    return 1.0;
	else
	    return factor;
}

