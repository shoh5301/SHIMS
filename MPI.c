#include "MPI.h"

void GGmodule_MPI(int**** grid,double jc[],int dir[],int out[],char ttem0[],double melt[],int pbc[],double const sang[],double rscale[],double tbc[][2]){// Grain growth
	clock_t start,end;
	char tcheck[2],temp[100],Tfun[50]; //Tfun: postfix equation for T-t rel.
	int i,j,k,l,o,mcs=0,total=0,step=0,stepdiv[2],check=2,dmc[3]={1};
	FILE* log=NULL; //     rtim 0 = total time // 1 = dt // 2 = junk // 3 = Lamda (mm / sites) // 4 = time range for T
	double realt,rtim[5]={0},realE,pmob,Tmelt,factor;

	int**** mpigrid=NULL;
	int mpixlen[3],mpix0[3],mpimode=0,mpixnum,mpidir[3];//count=0;
// mpixlen 0 = avg. length 1 = half avg. length 2 = actual length for each process
// mpix0 0 = original start point for each process 1 = 1st step start point 2 = 2nd step start point

	if(out[3]<=0){ // Time criteria is MCS
		dmc[0]=out[4]; // data print
		dmc[1]=out[1]; // map print
	}
    if(n_rank==0){
	printf("\n####### MODULE : GG #########\n");
// Grian growth only
	factor=melt[1];
	if(rscale[0]>0)
		rtim[3]=rscale[1]*1E3;//(m to mm)
	datout(dir,grid,0,out,rtim,rscale[0]); // Write data file & get L0
	rtim[2]=0;
    }
	Tmelt=rscale[2];
	total=dir[0]*dir[1]*dir[2];
	tok(ttem0,tcheck); //tcheck[0]: 0 for MCS, 1 for time(s)
	rtim[4]=-1; // junk value for MCS=0
	if(rscale[0]!=0)
	    rtim[1]=rscale[1]*rscale[1]*rscale[7]/rscale[5]; // Time is independent to temperature!
		//     lamda^2     *  K0_mc  /  K0_e

// MPI processing
	// dir, pbc, jc, rtim, factor, rscale, melt, out
	// from rank #0 to others
	MPI_Bcast(&factor,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(rtim,5,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(tcheck,2,MPI_CHAR,0,MPI_COMM_WORLD);

	mpixlen[0]=(int)((double)dir[0]/(double)n_size);
	while(n_size*mpixlen[0]<dir[0])
		mpixlen[0]++;

	mpix0[0]=mpixlen[0]*(n_size-n_rank-1);
	// when size 4:
	// rank 3 -> mpix0 = 0, rank 2 -> xlen, 1 -> xlen*2, 0 -> xlen*3 < dir[0] (<= xlen*4)
	// mpixlen 0 = original avg. length 1 = half value 2 = actual length in each process
	if(n_rank!=0)
	    mpixlen[2]=mpixlen[0];
	else // n_rank ==0
	    mpixlen[2]=(dir[0]-(n_size-1)*mpixlen[0]);
	mpixlen[1]=mpixlen[2]/2;

	stepdiv[0]=mpixlen[1]*dir[1]*dir[2];
	stepdiv[1]=(mpixlen[2]-mpixlen[1])*dir[1]*dir[2];
	mpidir[0]=mpixlen[2]+2;
	mpidir[1]=dir[1];
	mpidir[2]=dir[2];

	mpigrid=MPI_distr(grid,mpigrid,mpix0[0],mpixlen,dir,pbc[0]);

    if(n_rank==0){
	printf("\n << Monte Carlo Simulation Started >>\n");
	printf("    Sample size : %d X %d X %d\n",dir[0],dir[1],dir[2]);
	if(out[3]<=0) // MCS criteria
		printf("    Total MCS %d steps, Print Grain Map every %d steps, Print Result every %d steps\n",out[0],out[1],out[4]);
	else // sec criteria
		printf("    Total time %.3f sec, Print Grain Map every %.3f sec, Print Result every %.3f sec\n",(float)out[0]/1000,(float)out[1]/1000,(float)out[4]/1000);
	if(rscale[0]!=0){
	    if(factor==1.0)
		printf("    1 MCS = %E sec\n",rtim[1]);
	   else{
		printf("    1 MCS = %E sec x %E = %E sec\n",rtim[1],factor,rtim[1]*factor);
		rtim[1]=rtim[1]*factor;
	    }
	}
	if(out[3]<0)
		puts("    # Simulation will automatically stop when there is no liquid... #");
	printf("    %d processes: MPI activated\n\n",n_size);
	start=clock();
    }
	srand(time(NULL)+n_rank);

	MPI_Bcast(rtim,5,MPI_DOUBLE,0,MPI_COMM_WORLD);
	while(1){//	srand(time(NULL));
		// Temperature & real time setting
		if(rtim[4]!=0&&rtim[0]>rtim[4]){
			tok(ttem0,temp);
			rtim[4]=atof(temp); //rtim[4] is time range for T func.
			tok(ttem0,temp);
			change_postfix(temp,Tfun); // Tfun is postfix eq.
			if(n_rank==0){
			    printf(" # New T(t) function accepted: T(t) = %s",temp);
			    if(tcheck[0]=='0') // MCS variable
				printf(" after t = %d MCS\n",mcs);
			    else
				printf(" after t = %.2f sec\n",rtim[0]);
			}
			check=-1;
		}
		if(check!=0){ // if T changes
			if(tcheck[0]=='0'){ // MCS variable
				realt=Tcalc_MCS(Tfun,mcs);
				if(realt!=Tcalc_MCS(Tfun,mcs+1)) // Not isothermal
					check=1;
				else // New T, isothermal
					check=0;
			}else{	// sec variable
				realt=Tcalc_sec(Tfun,rtim[0]);
				if(realt!=Tcalc_sec(Tfun,rtim[0]+1)) // Not isothermal
					check=1;
				else // New T, isothermal
					check=0;
			}
			if(realt>Tmelt)
				printf("@@@ WARNING: Higher temperature then melting temperature: T = %.1f K vs. Tm = %.1f K ... @@@\n",realt,Tmelt);
			realE=rscale[1]*rscale[1]/(kB*realt);
//			     dx^2 / (K T)
			pmob=factor*exp(-rscale[3]*(1/realt-1/Tmelt)*REVRGAS);
//					Q  * (1/T-1/Tm)  / R
		}
		mcs++;

//Initially, first mpimode=0 @ function declaration
// MC trial #1 (0.5 MCS) ////////////////////////////////////////////////////////////////////////////////
		if(mpimode==0){ // former half first, later half second
		    mpix0[1]=1;
		    mpix0[2]=mpixlen[1]+1;
		    step=stepdiv[0];
		}else{ //mpimode==1
		    mpix0[1]=mpixlen[1]+1;
		    mpix0[2]=1;
		    step=stepdiv[1];
		}
		mpixnum=step/(dir[1]*dir[2]);
//if(n_rank==0);printf("mpix0[1] %d [2]%d mpixnum %d stepdiv 0 %d 1 %d\n",mpix0[1],mpix0[2],mpixnum,stepdiv[0],stepdiv[1]);
		while(step!=0){	
		    i=(int)(rand()/(RAND_MAX/mpixnum+1))+mpix0[1];
		    j=(int)(rand()/(RAND_MAX/dir[1]+1));
		    k=(int)(rand()/(RAND_MAX/dir[2]+1));

		    MCtrial_GG(i,j,k,mpigrid,mpidir,pbc,jc,tbc,pmob,realE);
		    step--;	// 1 Attempt
		}// 0.5 MCS
/////////////////// boundary 동기화
		MPI_BCsync(mpigrid,mpixlen[2],dir,pbc[0],mpimode);

// MC trial #2 (last 0.5 MCS) ////////////////////////////////////////////////////////////////////////////////
		if(mpimode==0)
		    step=stepdiv[1];
		else
		    step=stepdiv[0];
		mpixnum=step/(dir[1]*dir[2]);

		while(step!=0){	
		    i=(int)(rand()/(RAND_MAX/mpixnum+1))+mpix0[2];
		    j=(int)(rand()/(RAND_MAX/dir[1]+1));
		    k=(int)(rand()/(RAND_MAX/dir[2]+1));

		    MCtrial_GG(i,j,k,mpigrid,mpidir,pbc,jc,tbc,pmob,realE);
		    step--;	// 1 Attempt
		}// 0.5 MCS
///  나머지 boundary condition 동기화 ////////////////////////////////
		MPI_BCsync(mpigrid,mpixlen[2],dir,pbc[0],mpimode-1);
	
		if(n_rank==0) // for next mpimode
		    mpimode=(int)(rand()/(RAND_MAX/2+1));
		MPI_Bcast(&mpimode,1,MPI_INT,0,MPI_COMM_WORLD);

///////////////////////////////////////////////////puts("1 MCS Done");
		if(rscale[0]>0)
		    rtim[0]=rtim[0]+rtim[1];
		else
		    rtim[0]=mcs;

		if(out[3]<=0){ // Time criteria is MCS
		    dmc[0]--;
		    dmc[1]--;
		    if(dmc[0]==0){ // Data out
			dmc[0]=out[4];
			MPI_gather(grid,mpigrid,mpix0[0],mpixlen,dir);
			if(n_rank==0){
			  printf(" %dth MC Step was taken...\n",mcs);
			  datout(dir,grid,mcs,out,rtim,rscale[0]); // Write data file
			}
		    }
		    if(dmc[1]==0){ // Create graphics file
			dmc[1]=out[1];
			if(dmc[0]!=out[4])
			    MPI_gather(grid,mpigrid,mpix0[0],mpixlen,dir);
			if(n_rank==0)
			    fileout(dir,grid,mcs);
		    }
		    if(mcs==out[0])
		  	break;
		}else{ // Time criteria is sec
		    if(dmc[0]*out[4]<=rtim[0]*1000){ // Data out
			dmc[0]++;
			MPI_gather(grid,mpigrid,mpix0[0],mpixlen,dir);
			if(n_rank==0){
			    printf(" %.5f sec is passed... where T = %.1f\n",rtim[0],realt);
			    datout(dir,grid,mcs,out,rtim,rscale[0]); // Write data file
			}
		    }
		    if(dmc[1]*out[1]<=rtim[0]*1000){ // Create graphics file
			dmc[1]++;
			MPI_gather(grid,mpigrid,mpix0[0],mpixlen,dir);
			if(n_rank==0)
			    fileout(dir,grid,mcs);
		    }
		    if(rtim[0]*1000>=out[0])
		  	break;
		}
	}

    if(n_rank==0){
	if(out[3]<=0){ // Time criteria is MCS
		if(dmc[0]!=out[4]) // Check unsaved MCSs
		    datout(dir,grid,mcs,out,rtim,rscale[0]);
		if(dmc[1]!=out[1])
		    fileout(dir,grid,mcs);
	}else{ // Time criteria is sec
		if(dmc[0]*out[4]!=rtim[0]*1000)
			fileout(dir,grid,mcs);
		if(dmc[1]*out[1]<=rtim[0]*1000) // Create graphics file
			datout(dir,grid,mcs,out,rtim,rscale[0]);
	}

	end=clock();
	printf("\nExecution time: %d hour %d min %d sec\n",(int)((end-start)/CLOCKS_PER_SEC)/3600,(int)(((end-start)/CLOCKS_PER_SEC)%3600)/60,(int)(((end-start)/CLOCKS_PER_SEC)%3600)%60);
	printf("== Simulation is done ==\n");
    }
	free4d(mpixlen[2]+2,dir[1],dir[2],mpigrid);

    return;
}

void SOLmodule_MPI(int**** grid,double jc[],int dir[],int out[],double melt[],int pbc[],double rscale[],double tbc[][2]){
	clock_t start,end,ttime;
	char temp[100],Tfun[50]; //Tfun: postfix equation for T-t rel.
	int i,j,k,l,m,n,o,situ;
	int mcs=0,total=0,step=0,dmc[3]={1},fdm_loop=0,mpidir[3];
	int mode[4]={1,1,0,0};
	// [0] == lat_mode,[1] == fit_mode,[2] == T_mode, [3] == P_mode;
	// P_mode > 0 -> homogeneous phases, so accelerate FDM // P_mode < 0 -> fixed (isothermal)
	double k_per_cp,Tmax[3],Tinfo[2]={0},dtsave[2]; //dtsave 0 = dtfdm, 1 = prevdt
	// Tmax 0 == Tliq, 1 == T_init for solidification, 2 == freezing range
	// Tinfo 0 = minimum among liquid, 1 = maximum among solid, Pmax 0 for GG, 1 for solidification
	double Vsite,Asite,fnucl,factor[2]; // factor 0 = g.g. 1 = solid nucl. & growth
	double pcps[3]; // 0 == pcmax, 1 == psmax, 2 == 1/(pcmax+psmax)
	double rtim[5]={0},dEif,prob,dGfor,dGfus,dtdMCS[2]={0}; //dtdMCS 0 = g.g. 1 = solid nucl. & growth
	double Tnow,dHfus[2],fhetero_s,Nsite[2];
	// dHfus 0 == dHfus_m3, dHfus 1 == dH_per_Cp
	//Nsite[2]; //Nsite 0 == volumic density * voxel volume, 1 == plane density *voxel area
	double ****fdmt=NULL; // 0 is T for now, 1 is T for next, 2 is for latent heat info.
	FILE* ftemp=NULL;
//  rtim 0 = total time // 1 = dt_MCS in simul // 2 = solid fraction measurement (1 to yes) // 3 = Lamda (mm / sites) // 4 = time range for T

        int**** mpigrid=NULL;
	double**** mpitgrid=NULL;
        int mpixlen[3],mpix0[3],mpimode=0,mpixnum,stepdiv[2];//count=0;
// mpixlen 0 = avg. length 1 = half avg. length 2 = actual length for each process
// mpix0 0 = original start point for each process 1 = 1st step start point 2 = 2nd step start point

	if(out[5]>1)
	    rtim[2]=1;
	else
	    rtim[2]=0;

	// latent heat is considered only when lat_mode==1
	mode[MLAT]=(int)melt[6];
	// nucleation is considered only when fit_mode==1
	mode[MFIT]=(int)rscale[11];

	Vsite=rscale[1]*rscale[1]*rscale[1];
	Asite=rscale[1]*rscale[1];

	Nsite[0]=NAV/rscale[9]*Vsite;
	Nsite[1]=pow(Nsite[0],2.0/3.0);//*Asite;

	dHfus[0]=rscale[8]/rscale[9]; // J/mol -> J/m3
	dHfus[1]=rscale[8]/melt[4];
//		     dHfus (J/mol) / Cp (J/K mol) = K
	k_per_cp=melt[3]/melt[4]*rscale[9]; // (W / m K ) / (J/ m3 K) = W m2 / J
//	alpha=dt*melt[3]/(cp_m3*rscale[1]*rscale[1]);
	// Latent heat retroactive application scheme
/*	// first, thermal diffusivity alpha = k / (Cp / Vm) = k Vm / Cp
	melt[6]=melt[3]*rscale[9]/melt[4];
	// next, time for thermal diffusion to voxel size
	melt[6]=rscale[1]*rscale[1]/melt[6];
	// finally, 
	melt[6]=dHfus_per_Cp*(melt[1]/melt[6]);*/
	dtsave[0]=melt[1]; // original dt for FDM
	melt[6]=dHfus[1]*melt[1]/(rscale[1]*rscale[1])*melt[3]*rscale[9]/melt[4];
	if(melt[6]>dHfus[1])
	    melt[6]=dHfus[1];

	//      k    *dt / Cp    * dx^2
	Tmax[0]=rscale[2]; //Tl
	Tmax[1]=Tmax[0]-rscale[14]; // Critical supercooling
	Tmax[2]=rscale[2]-melt[0]; // freezing range

	// for solidification
	pcps[0]=Asite*(26*(jc[1]-jc[4])); // heterog. nucleation surrounded by other grains: maximum pc value for "a solid nucleus"
	pcps[1]=Vsite*dHfus[0]*(1-KAUZMANN); // 1-Tk/Tm where Tk is Kauzmann Temperature for glass transition
//	pcps[1]=Vsite*dHfus[0];
	pcps[2]=1/(pcps[0]+pcps[1]);

	if(out[3]<=0){ // Time criteria is MCS
	    dmc[0]=out[4]; // data print
	    dmc[1]=out[1]; // map print
	    dmc[2]=out[6]; // dtdMCs optimization
	}
	total=dir[0]*dir[1]*dir[2];
    if(n_rank==0){
	printf("\n####### MODULE : SOL #########\n");
	rtim[3]=rscale[1];//(m to mm)
	datout(dir,grid,0,out,rtim,rscale[0]); // Write data file
    }

       // dir, pbc, jc, rtim, factor, rscale, melt, out
        // from rank #0 to others
	MPI_Bcast(factor,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(rtim,5,MPI_DOUBLE,0,MPI_COMM_WORLD);

        mpixlen[0]=(int)((double)dir[0]/(double)n_size);
        while(n_size*mpixlen[0]<dir[0])
                mpixlen[0]++;
        mpix0[0]=mpixlen[0]*(n_size-n_rank-1);
        // when size 4:
        // rank 3 -> mpix0 = 0, rank 2 -> xlen, 1 -> xlen*2, 0 -> xlen*3 < dir[0] (<= xlen*4)
        // mpixlen 0 = original avg. length 1 = half value 2 = actual length in each process
        if(n_rank!=0)
            mpixlen[2]=mpixlen[0];
        else // n_rank ==0
            mpixlen[2]=(dir[0]-(n_size-1)*mpixlen[0]);
        mpixlen[1]=mpixlen[2]/2;

        stepdiv[0]=mpixlen[1]*dir[1]*dir[2];
        if(mpixlen[2]%2==0)
            stepdiv[1]=stepdiv[0];
        else
            stepdiv[1]=(mpixlen[1]+1)*dir[1]*dir[2];

	mpidir[0]=mpixlen[2]+2;
	mpidir[1]=dir[1];
	mpidir[2]=dir[2];

	mpigrid=MPI_distr(grid,mpigrid,mpix0[0],mpixlen,dir,pbc[0]);

    if(n_rank==0){
	// T at each simul. voxel 
	if(melt[2]<0){
	    ftemp=fopen(INCAR,"r");
	    if(find_keyword(ftemp,"T0",temp)!=0){ // ==0 when there is no error
		printf(" ### INPUT FILE ERROR in reading temperature profile file name from INCAR... \n");
		printf("     Program aborted...\n");
		return;
	    }
	    fclose(ftemp);
	    fdmt=(double****)calloc(3,sizeof(double***));
	    fdmt[0]=read_Tmap(fdmt[0],dir,temp,Nsite,1,rscale[2]); //Nsite is just a dummy
	    if(fdmt[0]==NULL){
		printf(" ### ERROR in temperature profile..\n       Program aborted...\n");
		return;
	    }
	    printf("\n");
	}
    }

	if(melt[2]<0){
	  if(n_rank==0){
	    // FDM T profile grid allocation & T0 setting
	    for(i=1;i<=2;i++)
		fdmt[i]=alloc3d(dir[0],dir[1],dir[2],fdmt[i]);
	  }else
		fdmt=(double****)calloc(3,sizeof(double***));

	  mpitgrid=MPI_Tdistr(fdmt,mpitgrid,mpix0[0],mpixlen,dir,pbc[0]);
	}else if(melt[2]==0){ //frozen temperature condition
	    mode[MTEM]=1; //frozen temperature approximation mode on
//	T*-G*
	    mode[MLAT]=0;
	}else{
	  fdmt=(double****)calloc(3,sizeof(double***));
	  if(n_rank==0){
	    for(i=2;i>=0;i--)
		fdmt[i]=alloc3d(dir[0],dir[1],dir[2],fdmt[i]);
	    for(i=dir[0]-1;i>=0;i--){
		for(j=dir[1]-1;j>=0;j--){
			for(k=dir[2]-1;k>=0;k--)
			    fdmt[0][i][j][k]=melt[2];
	    	}
	    }
	  }
	  mpitgrid=MPI_Tdistr(fdmt,mpitgrid,mpix0[0],mpixlen,dir,pbc[0]);
	}

	if(n_rank==0 && mode[MTEM]==0)
	    tproout(dir,fdmt[0],0);

///////// dt for solidification & melting?/////////////////////////////////////////
	dtdMCS[0]=rscale[1]*rscale[1]*rscale[7]/rscale[5];	// For grain growth
	dtdMCS[1]=rscale[1]*rscale[10];	// dt for solid nucleation & growth
//		   dx * (K_MC / K_e)
// Dynamic MCS setting
	out[5]=3-out[5];
// Now, out[5] 0 -> = t_max, 1 -> t_small, 2 -> t_gg

	if(n_rank==0)
	    printf(" ### NOTE: EVEN IN THE ACCELERATION MODE, FIRST 1 MCS IS IN MINIMUM MODE... ###\n\n");
	if(out[5]==2) // Only t_gg is considered as dt
	    rtim[1]=MCS_time_acceleration(dtdMCS,factor,factor,2); // here, second factor[2] is junk value
	else // smallest dt
	    rtim[1]=MCS_time_acceleration(dtdMCS,factor,factor,1);

// NOTE: P_Mode==0 at first MCS
	if(rtim[1]<=rscale[12]){
	    melt[1]=dtsave[0];
	    if(rtim[1]<melt[1]){ // fdm time > MC time
		melt[1]=rtim[1];
		fdm_loop=0;
		// in FDM loop, automatically ++
	    }else
		fdm_loop=(int)(rtim[1]/melt[1]);
	    melt[6]=dHfus[1]*melt[1]/(rscale[1]*rscale[1])*melt[3]*rscale[9]/melt[4];
	    if(melt[6]>dHfus[1])
		melt[6]=dHfus[1];
	    dtsave[2]=rtim[1];
	}else{
//	    if(rtim[1]>rscale[12]){
	    factor[0]=factor[0]*rscale[12]/rtim[1];
	    factor[1]=factor[1]*rscale[12]/rtim[1];
	    rtim[1]=rscale[12];
//	    }
	}
//////////////////////////////////////////////////////////////////////////////////

    if(n_rank==0){
	printf(" << Monte Carlo Simulation Started >>\n");
	printf("    Sample size : %d X %d X %d\n",dir[0],dir[1],dir[2]);
	if(out[3]<=0) // MCS criteria
		printf("    Total MCS %d steps, Print Grain Map every %d steps, Print Result every %d steps \n",out[0],out[1],out[4]);
	else // sec criteria
		printf("    Total time %.3f sec, Print Grain Map every %.3f sec, Print Result every %.3f sec \n",(float)out[0]/1000,(float)out[1]/1000,(float)out[4]/1000);

	printf("    1 MCS = %E sec,",rtim[1]);
	if(mode[MLAT]!=0){
	    printf(" Latent heat for S-L reaction is considered\n");
	    printf("    t_gg  = %E sec, t_solidi = %E sec",dtdMCS[0],dtdMCS[1]);
	    if(mode[MTEM]==0)
		printf(", t_FDM = %E sec ",melt[1]);
	    else
		printf(" // No FDM, as fronze temperature approximation is applied...\n");
	    if(rtim[1]<melt[1]){ // fdm time > MC time
		printf("(this will be adjusted because MCS < t_FDM)\n");
		melt[1]=rtim[1];
		melt[6]=dHfus[1]*melt[1]/(rscale[1]*rscale[1])*melt[3]*rscale[9]/melt[4];
		if(melt[6]>dHfus[1])
		    melt[6]=dHfus[1];

		fdm_loop=0;
		// in FDM loop, automatically ++
	    }else{
		fdm_loop=(int)(rtim[1]/melt[1]);
		printf("\n    Execute %d FDM loops for 1 MCS\n",fdm_loop);
	    }
	}else{
	    printf(" Latent heat for S-L reaction is ignored\n");
	    printf("    t_gg  = %E sec, t_solidi = %E sec\n",dtdMCS[0],dtdMCS[1]);
	}
	if(mode[MTEM]!=0)
	    printf("    Frozen temperature approximation is applied: T @IF = %.1f K + (%.2E K / vox) * distance\n",-melt[3],melt[4]); // T*-G

	if(mode[MFIT]==0)
	    printf("    Growth velocity fitting mode: no nucleation in the simulation...\n");

	printf("    %d processes: MPI activated\n\n",n_size);
    }

	MPI_Bcast(&mode,4,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&fdm_loop,1,MPI_INT,0,MPI_COMM_WORLD);
	if(fdm_loop==0){
	    MPI_Bcast(melt,8,MPI_DOUBLE,0,MPI_COMM_WORLD);
	    MPI_Bcast(&(rtim[1]),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	    MPI_Bcast(&(melt[6]),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
	rtim[4]=-1; // junk value for MCS=0
	dtsave[1]=rtim[1];

	MPI_Bcast(&(rtim[1]),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(factor,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&(mode[MPRO]),1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&fdm_loop,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&(melt[1]),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&(melt[6]),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(dtsave,2,MPI_DOUBLE,0,MPI_COMM_WORLD);

	start=clock();
	srand(time(NULL)+n_rank);
/////////////////////////////// MC trial & FDM /////////////////////////////////////////////////
	while(1){
	    mcs++;
//if(n_rank==0)
//printf("mcs %d rtim %E\n",mcs,rtim[1]);

	    if(mpimode==0){ // former half first, later half second
		mpix0[1]=1;
		mpix0[2]=mpixlen[1]+1;
		step=stepdiv[0];
	    }else{ //mpimode==1
		mpix0[1]=mpixlen[1]+1;
		mpix0[2]=1;
		step=stepdiv[1];
	    }
	    mpixnum=step/(dir[1]*dir[2]);

	    while(step!=0){
		i=(int)(rand()/(RAND_MAX/mpixnum+1))+mpix0[1];
		j=(int)(rand()/(RAND_MAX/dir[1]+1));
		k=(int)(rand()/(RAND_MAX/dir[2]+1));
		MCtrial_SOL(i,j,k,mpigrid,mpitgrid,mpidir,pbc,mode,melt,Tmax,dHfus,jc,tbc,Vsite,Asite,factor,Nsite,pcps,rscale[3]);
		step--;     // 1 Attempt
	    }// 0.5 MCS
	    MPI_BCsync(mpigrid,mpixlen[2],dir,pbc[0],mpimode-1);
// MC trial #2 (last 0.5 MCS) ////////////////////////////////////////////////////////////////////////////////
	    if(mpimode==0)
		step=stepdiv[1];
	    else
		step=stepdiv[0];
	    mpixnum=step/(dir[1]*dir[2]);

	    while(step!=0){
		i=(int)(rand()/(RAND_MAX/mpixnum+1))+mpix0[2];
		j=(int)(rand()/(RAND_MAX/dir[1]+1));
		k=(int)(rand()/(RAND_MAX/dir[2]+1));

		MCtrial_SOL(i,j,k,mpigrid,mpitgrid,mpidir,pbc,mode,melt,Tmax,dHfus,jc,tbc,Vsite,Asite,factor,Nsite,pcps,rscale[3]);
		step--;     // 1 Attempt
	    }// 0.5 MCS
///  나머지 boundary condition 동기화 ////////////////////////////////
	    MPI_BCsync(mpigrid,mpixlen[2],dir,pbc[0],mpimode-1);

	    if(n_rank==0) // for next mpimode
		mpimode=(int)(rand()/(RAND_MAX/2+1));
            MPI_Bcast(&mpimode,1,MPI_INT,0,MPI_COMM_WORLD);
// 1 MCS DONE

// mpixlen 0 = avg. length 1 = half avg. length 2 = actual length for each process
// mpix0 0 = original start point for each process 1 = 1st step start point 2 = 2nd step start point
//mpix0[0] <= i < mpix0[0]+mpixlen[2]

	    if(mode[MTEM]==0 && mode[MPRO]>=0)
		heat_transfer_MPI(mpitgrid,tbc,mpidir,dir,pbc,rscale[1],melt,fdm_loop,k_per_cp,Tmax[0],mpixlen[2]);
	    
	    rtim[0]+=rtim[1];

	    if(out[3]<=0){ // Time criteria is MCS
		dmc[0]--;
		dmc[1]--;
		if(mode[MTEM]==0 && out[5]==0){
		    dmc[2]--;
		    if(dmc[2]==0){ // time OPT
			dmc[2]=out[6];
			MPI_gather(grid,mpigrid,mpix0[0],mpixlen,dir);
			MPI_Tsync(fdmt[0],mpitgrid[0],mpix0[0],mpixlen,dir);
			MPI_Tsync(fdmt[2],mpitgrid[2],mpix0[0],mpixlen,dir);
// Active dt determination is not valid in MPI mode
//		if(n_rank==0)
//		    rtim[1]=active_dt_set(mcs,&(mode[MPRO]),&fdm_loop,dtsave,dtdMCS,dHfus,Tinfo,factor,Vsite,Asite,out,Tmax,pcps,grid,fdmt,dir,pbc,rscale,melt);
			if(n_rank==0)
			    rtim[1]=fixed_dt_set(mcs,&(mode[MPRO]),&fdm_loop,dtsave,dtdMCS,dHfus,Tinfo,factor,Vsite,Asite,out,Tmax,pcps,grid,fdmt,dir,pbc,rscale,melt); // */
//			rtim[1]=fixed_dt_set(mcs,&fdm_loop,dtsave,dtdMCS,dHfus,Tinfo,factor,Vsite,Asite,out,Tmax,pcps,mpigrid,mpitgrid,mpidir,pbc,rscale,melt);
//			situ=find_min_rank(rtim[1]);
			
			situ=0;
			MPI_Bcast(&(mode[MPRO]),1,MPI_DOUBLE,situ,MPI_COMM_WORLD);
			MPI_Bcast(&(rtim[1]),1,MPI_DOUBLE,situ,MPI_COMM_WORLD);
			MPI_Bcast(factor,2,MPI_DOUBLE,situ,MPI_COMM_WORLD);
			MPI_Bcast(&fdm_loop,1,MPI_INT,situ,MPI_COMM_WORLD);
			MPI_Bcast(&(melt[1]),1,MPI_DOUBLE,situ,MPI_COMM_WORLD);
			MPI_Bcast(&(melt[6]),1,MPI_DOUBLE,situ,MPI_COMM_WORLD);
			MPI_Bcast(dtsave,2,MPI_DOUBLE,situ,MPI_COMM_WORLD);

			if(rtim[1]<0){
			    if(n_rank==0)
				printf(" @@@ Liquid fraction is 0, solidification is DONE @@@\n");
			    break;
			}
		    }
		}
		if(dmc[0]==0){ // Data out
		    dmc[0]=out[4];
		    MPI_gather(grid,mpigrid,mpix0[0],mpixlen,dir);
		    if(n_rank==0){
			printf(" %dth MC Step was taken...\n",mcs);
			datout(dir,grid,mcs,out,rtim,rscale[0]); // Write data file
		    }
		}
		if(dmc[1]==0){ // Create graphics file
		    dmc[1]=out[1];
		    if(dmc[0]!=out[4])
			MPI_gather(grid,mpigrid,mpix0[0],mpixlen,dir);
		    if(mode[MTEM]==0)
			MPI_Tsync(fdmt[0],mpitgrid[0],mpix0[0],mpixlen,dir);
		    if(n_rank==0){
			fileout(dir,grid,mcs);
			if(mode[MTEM]==0)
			    tproout(dir,fdmt[0],mcs);
		    }
		}
		if(mcs==out[0])
		    break;
	    }else{ // Time criteria is sec
		if(dmc[0]*out[4]<=rtim[0]*1000){ // Data out
		    dmc[0]++;
		    MPI_gather(grid,mpigrid,mpix0[0],mpixlen,dir);
		    if(n_rank==0){
			printf(" %.5f sec is passed...\n",rtim[0]);
			datout(dir,grid,mcs,out,rtim,rscale[0]); // Write data file
		    }
		}
		if(dmc[1]*out[1]<=rtim[0]*1000){ // Create graphics file
		    dmc[1]++;
		    MPI_gather(grid,mpigrid,mpix0[0],mpixlen,dir);
		    if(mode[MTEM]==0)
			MPI_Tsync(fdmt[0],mpitgrid[0],mpix0[0],mpixlen,dir);
		    if(n_rank==0){
			fileout(dir,grid,mcs);
			if(mode[MTEM]==0)
			    tproout(dir,fdmt[0],mcs);
		    }
	        }
		if(rtim[0]*1000>=out[0])
		    break;
	    }
	}

	if(out[3]<=0){ // Time criteria is MCS
	    if(dmc[1]!=out[1]){
		MPI_gather(grid,mpigrid,mpix0[0],mpixlen,dir);
		if(mode[MTEM]==0)
		    MPI_Tsync(fdmt[0],mpitgrid[0],mpix0[0],mpixlen,dir);
	    }
	    if(n_rank==0){
		if(out[0]!=mcs)
		    out[0]=mcs;

		if(dmc[0]!=out[4]) // Check unsaved MCSs
			datout(dir,grid,mcs,out,rtim,rscale[0]);
		if(dmc[1]!=out[1]){
		    fileout(dir,grid,mcs);
		    if(mode[MTEM]==0)
			tproout(dir,fdmt[0],mcs);
		}
	    }
	}else{ // Time criteria is sec
	    if(dmc[0]*out[4]!=rtim[0]*1000){
		MPI_gather(grid,mpigrid,mpix0[0],mpixlen,dir);
		if(mode[MTEM]==0)
		    MPI_Tsync(fdmt[0],mpitgrid[0],mpix0[0],mpixlen,dir);
	    }
	    if(n_rank==0){
		if(out[0]!=rtim[0]*1000)
		    out[0]=rtim[0]*1000;

		if(dmc[0]*out[4]!=rtim[0]*1000)
		    datout(dir,grid,mcs,out,rtim,rscale[0]);
		if(dmc[1]*out[1]!=rtim[0]*1000){ // Create graphics file
		    fileout(dir,grid,mcs);
		    if(mode[MTEM]==0)
			tproout(dir,fdmt[0],mcs);
		}
	    }
	}

	if(n_rank==0){
	    end=clock();
	    printf("\nExecution time: %d hour %d min %d sec\n",(int)((end-start)/CLOCKS_PER_SEC)/3600,(int)(((end-start)/CLOCKS_PER_SEC)%3600)/60,(int)(((end-start)/CLOCKS_PER_SEC)%3600)%60);
	    printf("== Simulation is done ==\n");
	}

	free4d(mpixlen[2]+2,dir[1],dir[2],mpigrid);

	if(mode[MTEM]==0){
	    for(i=2;i>=0;i--)
		mpitgrid[i]=free3d(mpixlen[2]+2,dir[1],dir[2],mpitgrid[i]);
	    free(mpitgrid);
	    if(n_rank==0){
	      for(i=2;i>=0;i--)
		fdmt[i]=free3d(dir[0],dir[1],dir[2],fdmt[i]);
	    }
	    free(fdmt);
	}

	if(mode[MLAT]==0)
	    melt[3]=-1;
	else if(melt[3]<0)
	    melt[3]=0;
	return;
}

void AMmodule_MPI(int**** grid,double jc[],int dir[],int out[],double melt[],int pbc[],double rscale[],double tbc[][2],double am[]){
	clock_t start,end,ttime;
	char temp[100],Tfun[50]; //Tfun: postfix equation for T-t rel.
	int i,j,k,l,m,n,o,situ;
	int mcs=0,total=0,step=0,dmc[3]={1},fdm_loop=0,mpidir[3];
	int mode[4]={1,1,0,0};
	// [0] == lat_mode,[1] == fit_mode,[2] == T_mode, [3] == P_mode;
	// P_mode > 0 -> homogeneous phases, so accelerate FDM // P_mode < 0 -> fixed (isothermal)
	double k_per_cp,Tmax[3],Tinfo[2]={0},dtsave[2],mp[7]; //dtsave 0 = dtfdm, 1 = prevdt
	// Tmax 0 == Tmelt, 1 == T_init for solidification, 2 == Tmax for melting, 
	// Tinfo 0 = minimum among liquid, 1 = maximum among solid, Pmax 0 for GG, 1 for solidification
	double Vsite,Asite,fnucl,factor[2],pcps[3]; // factor 0 = g.g. 1 = solid nucl. & growth
	// 0 == pcmax, 1 == psmax, 2 == 1/(pcmax+psmax)
	double rtim[5]={0},dEif,prob,dGfor,dGfus,dtdMCS[2]={0}; //dtdMCS 0 = g.g. 1 = solid nucl. & growth
	double Tnow,dHfus[2],fhetero_s,Nsite[2],pos[2];
	// dHfus 0 == dHfus_m3, dHfus 1 == dH_per_Cp
	//Nsite[2]; //Nsite 0 == volumic density * voxel volume, 1 == plane density *voxel area
	double ****fdmt=NULL; // 0 is T for now, 1 is T for next, 2 is for latent heat info.
	FILE* ftemp=NULL;
//  rtim 0 = total time // 1 = dt_MCS in simul // 2 = solid fraction measurement (1 to yes) // 3 = Lamda (mm / sites) // 4 = time range for T

        int**** mpigrid=NULL;
	double**** mpitgrid=NULL;
        int mpixlen[3],mpix0[3],mpimode=0,mpixnum,stepdiv[2];//count=0;
// mpixlen 0 = avg. length 1 = half avg. length 2 = actual length for each process
// mpix0 0 = original start point for each process 1 = 1st step start point 2 = 2nd step start point

	if(out[5]>1)
	    rtim[2]=1;
	else
	    rtim[2]=0;

	// latent heat is considered only when lat_mode==1
	mode[MLAT]=(int)melt[6];
	// nucleation is considered only when fit_mode==1
	mode[MFIT]=(int)rscale[11];

	Vsite=rscale[1]*rscale[1]*rscale[1];
	Asite=rscale[1]*rscale[1];

	Nsite[0]=NAV/rscale[9]*Vsite;
	Nsite[1]=pow(Nsite[0],2.0/3.0);//*Asite;

	dHfus[0]=rscale[8]/rscale[9]; // J/mol -> J/m3
	dHfus[1]=rscale[8]/melt[4];
//		     dHfus (J/mol) / Cp (J/K mol) = K
	k_per_cp=melt[3]/melt[4]*rscale[9]; // (W / m K ) / (J/ m3 K) = W m2 / J
//	alpha=dt*melt[3]/(cp_m3*rscale[1]*rscale[1]);
	// Latent heat retroactive application scheme
/*	// first, thermal diffusivity alpha = k / (Cp / Vm) = k Vm / Cp
	melt[6]=melt[3]*rscale[9]/melt[4];
	// next, time for thermal diffusion to voxel size
	melt[6]=rscale[1]*rscale[1]/melt[6];
	// finally, 
	melt[6]=dHfus_per_Cp*(melt[1]/melt[6]);*/
	dtsave[0]=melt[1]; // original dt for FDM
	melt[6]=dHfus[1]*melt[1]/(rscale[1]*rscale[1])*melt[3]*rscale[9]/melt[4];
	if(melt[6]>dHfus[1])
		melt[6]=dHfus[1];

	//      k    *dt / Cp    * dx^2
	Tmax[0]=rscale[2];
	Tmax[1]=Tmax[0]-rscale[14]; // Critical supercooling
	Tmax[2]=rscale[2]-melt[0]; // freezing range

	// for solidification
	pcps[0]=Asite*(26*(jc[1]-jc[4])); // heterog. nucleation surrounded by other grains: maximum pc value for "a solid nucleus"
	pcps[1]=Vsite*dHfus[0]*(1-KAUZMANN); // 1-Tk/Tm where Tk is Kauzmann Temperature for glass transition
	pcps[2]=1/(pcps[0]+pcps[1]);

	if(out[3]<=0){ // Time criteria is MCS
	    dmc[0]=out[4]; // data print
	    dmc[1]=out[1]; // map print
	    dmc[2]=out[6]; // dtdMCs optimization
	}
	total=dir[0]*dir[1]*dir[2];
    if(n_rank==0){
	printf("\n####### MODULE : AM #########\n");
	rtim[3]=rscale[1];//(m to mm)
	datout(dir,grid,0,out,rtim,rscale[0]); // Write data file
    }

       // dir, pbc, jc, rtim, factor, rscale, melt, out
        // from rank #0 to others
	MPI_Bcast(factor,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(rtim,5,MPI_DOUBLE,0,MPI_COMM_WORLD);

        mpixlen[0]=(int)((double)dir[0]/(double)n_size);
        while(n_size*mpixlen[0]<dir[0])
                mpixlen[0]++;
        mpix0[0]=mpixlen[0]*(n_size-n_rank-1);
        // when size 4:
        // rank 3 -> mpix0 = 0, rank 2 -> xlen, 1 -> xlen*2, 0 -> xlen*3 < dir[0] (<= xlen*4)
        // mpixlen 0 = original avg. length 1 = half value 2 = actual length in each process
        if(n_rank!=0)
            mpixlen[2]=mpixlen[0];
        else // n_rank ==0
            mpixlen[2]=(dir[0]-(n_size-1)*mpixlen[0]);
        mpixlen[1]=mpixlen[2]/2;

        stepdiv[0]=mpixlen[1]*dir[1]*dir[2];
        if(mpixlen[2]%2==0)
            stepdiv[1]=stepdiv[0];
        else
            stepdiv[1]=(mpixlen[1]+1)*dir[1]*dir[2];

	mpidir[0]=mpixlen[2]+2;
	mpidir[1]=dir[1];
	mpidir[2]=dir[2];

	mpigrid=MPI_distr(grid,mpigrid,mpix0[0],mpixlen,dir,pbc[0]);

    if(n_rank==0){
	// T at each simul. voxel 
	if(melt[2]<0){
	    ftemp=fopen(INCAR,"r");
	    if(find_keyword(ftemp,"T0",temp)!=0){ // ==0 when there is no error
		printf(" ### INPUT FILE ERROR in reading temperature profile file name from INCAR... \n");
		printf("     Program aborted...\n");
		return;
	    }
	    fclose(ftemp);
	    fdmt=(double****)calloc(3,sizeof(double***));
	    fdmt[0]=read_Tmap(fdmt[0],dir,temp,Nsite,1,rscale[2]); //Nsite is just a dummy
	    if(fdmt[0]==NULL){
		printf(" ### ERROR in temperature profile..\n       Program aborted...\n");
		return;
	    }
	    printf("\n");
	}
    }

	if(melt[2]<0){
	  if(n_rank==0){
	    // FDM T profile grid allocation & T0 setting
	    for(i=1;i<=2;i++)
		fdmt[i]=alloc3d(dir[0],dir[1],dir[2],fdmt[i]);
	  }else
	    fdmt=(double****)calloc(3,sizeof(double***));
	  mpitgrid=MPI_Tdistr(fdmt,mpitgrid,mpix0[0],mpixlen,dir,pbc[0]);
	}else{
	  fdmt=(double****)calloc(3,sizeof(double***));
	  if(n_rank==0){
	    for(i=2;i>=0;i--){
		fdmt[i]=alloc3d(dir[0],dir[1],dir[2],fdmt[i]);
	    }
	    for(i=dir[0]-1;i>=0;i--){
		for(j=dir[1]-1;j>=0;j--){
			for(k=dir[2]-1;k>=0;k--)
			    fdmt[0][i][j][k]=melt[2];
	    	}
	    }
	  }
	  mpitgrid=MPI_Tdistr(fdmt,mpitgrid,mpix0[0],mpixlen,dir,pbc[0]);
	}

	if(n_rank==0)
	    tproout(dir,fdmt[0],0);

///////// dt for solidification & melting?/////////////////////////////////////////
	dtdMCS[0]=rscale[1]*rscale[1]*rscale[7]/rscale[5];	// For grain growth
	dtdMCS[1]=rscale[1]*rscale[10];	// dt for solid nucleation & growth
//		   dx * (K_MC / K_e)
// Dynamic MCS setting
	out[5]=3-out[5];
// Now, out[5] 0 -> = t_max, 1 -> t_small, 2 -> t_gg

	if(n_rank==0)
	    printf(" ### NOTE: EVEN IN THE ACCELERATION MODE, FIRST 1 MCS IS IN MINIMUM MODE... ###\n\n");
	if(out[5]==2) // Only t_gg is considered as dt
	    rtim[1]=MCS_time_acceleration(dtdMCS,factor,factor,2); // here, second factor[2] is junk value
	else // smallest dt
	    rtim[1]=MCS_time_acceleration(dtdMCS,factor,factor,1);

// NOTE: P_Mode==0 at first MCS
	if(rtim[1]<=rscale[12]){
	    melt[1]=dtsave[0];
	    if(rtim[1]<melt[1]){ // fdm time > MC time
		melt[1]=rtim[1];
		fdm_loop=0;
		// in FDM loop, automatically ++
	    }else
		fdm_loop=(int)(rtim[1]/melt[1]);
	    melt[6]=dHfus[1]*melt[1]/(rscale[1]*rscale[1])*melt[3]*rscale[9]/melt[4];
	    dtsave[2]=rtim[1];
	}else{
//	    if(rtim[1]>rscale[12]){
	    factor[0]=factor[0]*rscale[12]/rtim[1];
	    factor[1]=factor[1]*rscale[12]/rtim[1];
	    rtim[1]=rscale[12];
//	    }
	}
//////////////////////////////////////////////////////////////////////////////////

    if(n_rank==0){
	printf(" << Monte Carlo Simulation Started >>\n");
	printf("    Sample size : %d X %d X %d\n",dir[0],dir[1],dir[2]);
	if(out[3]<=0) // MCS criteria
		printf("    Total MCS %d steps, Print Grain Map every %d steps, Print Result every %d steps \n",out[0],out[1],out[4]);
	else // sec criteria
		printf("    Total time %.3f sec, Print Grain Map every %.3f sec, Print Result every %.3f sec \n",(float)out[0]/1000,(float)out[1]/1000,(float)out[4]/1000);

	printf("    1 MCS = %E sec,",rtim[1]);
	if(mode[MLAT]!=0){
	    printf(" Latent heat for S-L reaction is considered\n");
	    printf("    t_gg  = %E sec, t_solidi = %E sec",dtdMCS[0],dtdMCS[1]);
	    if(mode[MTEM]==0)
		printf(", t_FDM = %E sec ",melt[1]);
	    else
		printf(" // No FDM, as fronze temperature approximation is applied...\n");
	    if(rtim[1]<melt[1]){ // fdm time > MC time
		printf("(this will be adjusted because MCS < t_FDM)\n");
		melt[1]=rtim[1];
		melt[6]=dHfus[1]*melt[1]/(rscale[1]*rscale[1])*melt[3]*rscale[9]/melt[4];

		fdm_loop=0;
		// in FDM loop, automatically ++
	    }else{
		fdm_loop=(int)(rtim[1]/melt[1]);
		printf("\n    Execute %d FDM loops for 1 MCS\n",fdm_loop);
	    }
	}else{
	    printf(" Latent heat for S-L reaction is ignored\n");
	    printf("    t_gg  = %E sec, t_solidi = %E sec\n",dtdMCS[0],dtdMCS[1]);
	}
	if(am[0]<0)
		printf("    Circular shaped melt pool\n");
	else if(am[0]>0)
		printf("    Gaussian LASER irradiation\n");
	else
		printf("    Teardropshaped melt pool\n");
	if(out[3]<0)
		puts("    # Simulation will automatically stop when there is no liquid... #");
	printf("    %d processes: MPI activated\n\n",n_size);

    }

// Melt pool info.
// am 0 = Melt pool mode (-1 = Circle, 0 = TEARDROP,  1 = Gaussian) 
// 1 = Melt pool direction (negative -> X, positive -> Y)
// Start position (x = 6, y = 7), 8 = scan speed
// when am[0]==-1: 2 = Radius, 9 = Boundary Temperature
// when am[0]==0:  2 = Width  3 = Cap+Tail length  4 = Depth  5 = Geometry (Teardrop shape), 9 = Boundary Temperature
// when am[0]==1:  2 = SX     3 = SY               4 = SZ    5 = POW  9 = ABS    (Gaussian)
	if(am[0]==0){
	    mp[0]=4*atan(sqrt(2-sqrt(3))); // temporarily, mp[0] = theta
	    mp[1]=am[3]/2;
	    mp[2]=am[2]/(2*sin(mp[0])*pow(sin(mp[0]*0.5),am[5]));
	    mp[3]=am[4]/(sin(mp[0])*pow(sin(mp[0]*0.5),am[5]));
	}else if(am[0]>0){ // GAUSSIAN
// mp 0 -> distance for jumping
//    5 -> aP/(2 pi)^1.5/(sx sy sz)
//    2 -> 1/sx^2    3 -> 1/sy^2   4 -> 1/sz^2 // all in voxel^2 unit
//    1 -> ap/(2 pi)^1.5/(sx sy sz) converted to K / s
    
//    GAUSS3D = (2 pi)^-1.5
//	    mp[1]=am[9]*am[5]*GAUSS3D/(am[2]*am[3]*am[4]); // W / m3 = J / m3 s
//	    mp[1]=mp[1]*rscale[9]*rtim[1]/melt[4]; // J / m3 s * m3 / mol * s / J/K mol = K
/*	    mp[5]=am[9]*am[5]*GAUSS3D/(am[2]*am[3]*am[4])*rscale[9]/melt[4];
	    mp[1]=mp[5]*rtim[1];*/
//    GAUSS2D = (2 pi) ^ -1
//	    mp[1]=am[9]*am[5]*GAUSS2D/(am[2]*am[3]); // W / m2 = J / m2 s
//	    mp[1]=mp[1]*pow(rscale[9],2/3)*rtim[1]/melt[4]; // J / m2 s * m2 / mol * s / J/K mol = K
	    mp[5]=am[9]*am[5]*GAUSS2D/(am[2]*am[3])*pow(rscale[9],2/3)/melt[4];
	    mp[1]=mp[5]*rtim[1];
// sx, sy, sz unit m -> sx, sy, sz unit voxel
	    am[2]=am[2]/rscale[1];
	    am[3]=am[3]/rscale[1];
	    am[4]=am[4]/rscale[1];
	    mp[2]=1/(am[2]*am[2]);
	    mp[3]=1/(am[3]*am[3]);
	    mp[4]=1/(am[4]*am[4]);
	}
	mp[0]=-am[8]/rscale[1]*rtim[1]-1;
	pos[0]=am[6]-(double)mpix0[0]+1;
	pos[1]=am[7];

	MPI_Bcast(&mode,4,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&fdm_loop,1,MPI_INT,0,MPI_COMM_WORLD);
	if(fdm_loop==0){
	    MPI_Bcast(melt,8,MPI_DOUBLE,0,MPI_COMM_WORLD);
	    MPI_Bcast(&(rtim[1]),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	    MPI_Bcast(&(melt[6]),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
	rtim[4]=-1; // junk value for MCS=0
	dtsave[1]=rtim[1];

	MPI_Bcast(&(rtim[1]),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(factor,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&(mode[MPRO]),1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&fdm_loop,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&(melt[1]),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&(melt[6]),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(dtsave,2,MPI_DOUBLE,0,MPI_COMM_WORLD);

	start=clock();
	srand(time(NULL)+n_rank);
/////////////////////////////// MC trial & FDM /////////////////////////////////////////////////

	while(1){
	    mcs++;

//if(n_rank==0)
//printf("mcs %d\n",mcs);
//printf("[rank %d] Start am9 %f\n",n_rank,am[9]);
	    if(am[9]!=0){ // TMP or ABS == 0 -> No melt pool calculation
	      if(am[0]<=0){ // TEARDROP mode or CIRCLE mode
		if(mp[0]<0){ // 1st MCS
		    mp[0]=0;
		    if(MPI_MeltPoolCheck(am[0],am[1],mp,am[7],pos[0],mpixlen[2],mpidir[0])!=0)
			MeltPool(mpigrid,mpitgrid,am,mp,mpidir,pbc,Tmax[0],dHfus[1],pos);
		}else{
		    mp[0]=am[8]/rscale[1]*rtim[1];
		    MeltPool_movement(pos,am[1],mp[0]);
		    if(MPI_MeltPoolCheck(am[0],am[1],mp,am[7],pos[0],mpixlen[2],mpidir[0])!=0)
			MeltPool(mpigrid,mpitgrid,am,mp,mpidir,pbc,Tmax[0],dHfus[1],pos);
		}
		if(am[6]<0 || am[6]>dir[0] || am[7]<0 || am[7]>dir[1])
		    am[9]=0; // 

	      }else if(am[0]>0){ // GAUSSIAN mode // 짧은 timestep 가져가면서 계속 지지기
		if(mp[0]<0) // 1st MCS
		    mp[0]=0;
		else{
		    mp[0]=am[8]/rscale[1]*rtim[1];
		    mp[1]=mp[5]*rtim[1];
		    MeltPool_movement(pos,am[1],mp[0]);
		}
		MeltPool(mpigrid,mpitgrid,am,mp,mpidir,pbc,Tmax[0],dHfus[1],pos);
/*		if(am[6]<0 || am[6]>dir[0] || am[7]<0 || am[7]>dir[1])
		    am[9]=0;*/
		if(pos[0]>mpixlen[2])
		    am[9]=0;
	      }
	    }
//	    if(jc[7]!=0) // PTCL introduction
//		Particle_Form(mpigrid,mpidir,pbc,jc[7]);

	    if(mode[MPRO]>=0){
		heat_transfer_MPI(mpitgrid,tbc,mpidir,dir,pbc,rscale[1],melt,fdm_loop,k_per_cp,Tmax[0],mpixlen[2]);
		if(am[9]>Tmax[0])
		    check_melting(mpigrid,mpitgrid,Tmax[0],mpidir);
	    }

	    if(mpimode==0){ // former half first, later half second
		mpix0[1]=1;
		mpix0[2]=mpixlen[1]+1;
		step=stepdiv[0];
	    }else{ //mpimode==1
		mpix0[1]=mpixlen[1]+1;
		mpix0[2]=1;
		step=stepdiv[1];
	    }
	    mpixnum=step/(dir[1]*dir[2]);

	    while(step!=0){
		i=(int)(rand()/(RAND_MAX/mpixnum+1))+mpix0[1];
		j=(int)(rand()/(RAND_MAX/dir[1]+1));
		k=(int)(rand()/(RAND_MAX/dir[2]+1));
		MCtrial_SOL(i,j,k,mpigrid,mpitgrid,mpidir,pbc,mode,melt,Tmax,dHfus,jc,tbc,Vsite,Asite,factor,Nsite,pcps,rscale[3]);
		step--;     // 1 Attempt
	    }// 0.5 MCS
	    MPI_BCsync(mpigrid,mpixlen[2],dir,pbc[0],mpimode-1);
// MC trial #2 (last 0.5 MCS) ////////////////////////////////////////////////////////////////////////////////
	    if(mpimode==0)
		step=stepdiv[1];
	    else
		step=stepdiv[0];
	    mpixnum=step/(dir[1]*dir[2]);

	    while(step!=0){
		i=(int)(rand()/(RAND_MAX/mpixnum+1))+mpix0[2];
		j=(int)(rand()/(RAND_MAX/dir[1]+1));
		k=(int)(rand()/(RAND_MAX/dir[2]+1));

		MCtrial_SOL(i,j,k,mpigrid,mpitgrid,mpidir,pbc,mode,melt,Tmax,dHfus,jc,tbc,Vsite,Asite,factor,Nsite,pcps,rscale[3]);
		step--;     // 1 Attempt
	    }// 0.5 MCS
///  나머지 boundary condition 동기화 ////////////////////////////////
	    MPI_BCsync(mpigrid,mpixlen[2],dir,pbc[0],mpimode-1);

	    if(n_rank==0) // for next mpimode
		mpimode=(int)(rand()/(RAND_MAX/2+1));
            MPI_Bcast(&mpimode,1,MPI_INT,0,MPI_COMM_WORLD);
// 1 MCS DONE
//printf("[rank %d] MCS DONE\n",n_rank);
// mpixlen 0 = avg. length 1 = half avg. length 2 = actual length for each process
// mpix0 0 = original start point for each process 1 = 1st step start point 2 = 2nd step start point
//mpix0[0] <= i < mpix0[0]+mpixlen[2]

	    rtim[0]+=rtim[1];

	    if(out[3]<=0){ // Time criteria is MCS
		dmc[0]--;
		dmc[1]--;

		if(out[5]==0){
		    dmc[2]--;
		    if(dmc[2]==0){ // time OPT
			dmc[2]=out[6];
			MPI_gather(grid,mpigrid,mpix0[0],mpixlen,dir);
			MPI_Tsync(fdmt[0],mpitgrid[0],mpix0[0],mpixlen,dir);
			MPI_Tsync(fdmt[2],mpitgrid[2],mpix0[0],mpixlen,dir);
// Active dt determination is not valid in MPI mode
			    if(n_rank==0)
				rtim[1]=fixed_dt_set(mcs,&(mode[MPRO]),&fdm_loop,dtsave,dtdMCS,dHfus,Tinfo,factor,Vsite,Asite,out,Tmax,pcps,grid,fdmt,dir,pbc,rscale,melt); // */
//			rtim[1]=fixed_dt_set(mcs,&fdm_loop,dtsave,dtdMCS,dHfus,Tinfo,factor,Vsite,Asite,out,Tmax,pcps,mpigrid,mpitgrid,mpidir,pbc,rscale,melt);
//			situ=find_min_rank(rtim[1]);
			    situ=0;
			    MPI_Bcast(&(mode[MPRO]),1,MPI_DOUBLE,situ,MPI_COMM_WORLD);
			    MPI_Bcast(&(rtim[1]),1,MPI_DOUBLE,situ,MPI_COMM_WORLD);
			    MPI_Bcast(factor,2,MPI_DOUBLE,situ,MPI_COMM_WORLD);
			    MPI_Bcast(&fdm_loop,1,MPI_INT,situ,MPI_COMM_WORLD);
			    MPI_Bcast(&(melt[1]),1,MPI_DOUBLE,situ,MPI_COMM_WORLD);
			    MPI_Bcast(&(melt[6]),1,MPI_DOUBLE,situ,MPI_COMM_WORLD);
			    MPI_Bcast(dtsave,2,MPI_DOUBLE,situ,MPI_COMM_WORLD);
			if(rtim[1]<0){
			    if(n_rank==0)
				printf(" @@@ Liquid fraction is 0, solidification is DONE @@@\n");
			    break;
			}
		    }
		}
		if(dmc[0]==0){ // Data out
		    dmc[0]=out[4];
		    MPI_gather(grid,mpigrid,mpix0[0],mpixlen,dir);
		    if(n_rank==0){
			printf(" %dth MC Step was taken...\n",mcs);
			datout(dir,grid,mcs,out,rtim,rscale[0]); // Write data file
		    }
		}
		if(dmc[1]==0){ // Create graphics file
		    dmc[1]=out[1];
		    if(dmc[0]!=out[4])
			MPI_gather(grid,mpigrid,mpix0[0],mpixlen,dir);
		    MPI_Tsync(fdmt[0],mpitgrid[0],mpix0[0],mpixlen,dir);
		    if(n_rank==0){
			fileout(dir,grid,mcs);
			tproout(dir,fdmt[0],mcs);
		    }
		}
		if(mcs==out[0])
		    break;
	    }else{ // Time criteria is sec
		if(dmc[0]*out[4]<=rtim[0]*1000){ // Data out
		    dmc[0]++;
		    MPI_gather(grid,mpigrid,mpix0[0],mpixlen,dir);
		    if(n_rank==0){
			printf(" %.5f sec is passed...\n",rtim[0]);
			datout(dir,grid,mcs,out,rtim,rscale[0]); // Write data file
		    }
		}
		if(dmc[1]*out[1]<=rtim[0]*1000){ // Create graphics file
		    dmc[1]++;
		    MPI_gather(grid,mpigrid,mpix0[0],mpixlen,dir);
		    MPI_Tsync(fdmt[0],mpitgrid[0],mpix0[0],mpixlen,dir);
		    if(n_rank==0){
			fileout(dir,grid,mcs);
			tproout(dir,fdmt[0],mcs);
		    }
	        }
		if(rtim[0]*1000>=out[0])
		    break;
	    }
	}

	if(out[3]<=0){ // Time criteria is MCS
	    if(dmc[1]!=out[1]){
		MPI_gather(grid,mpigrid,mpix0[0],mpixlen,dir);
		if(mode[MTEM]==0)
		    MPI_Tsync(fdmt[0],mpitgrid[0],mpix0[0],mpixlen,dir);
	    }
	    if(n_rank==0){
		if(out[0]!=mcs)
		    out[0]=mcs;

		if(dmc[0]!=out[4]) // Check unsaved MCSs
			datout(dir,grid,mcs,out,rtim,rscale[0]);
		if(dmc[1]!=out[1]){
		    fileout(dir,grid,mcs);
		    if(mode[MTEM]==0)
			tproout(dir,fdmt[0],mcs);
		}
	    }
	}else{ // Time criteria is sec
	    if(dmc[0]*out[4]!=rtim[0]*1000){
		MPI_gather(grid,mpigrid,mpix0[0],mpixlen,dir);
		if(mode[MTEM]==0)
		    MPI_Tsync(fdmt[0],mpitgrid[0],mpix0[0],mpixlen,dir);
	    }
	    if(n_rank==0){
		if(out[0]!=rtim[0]*1000)
		    out[0]=rtim[0]*1000;

		if(dmc[0]*out[4]!=rtim[0]*1000)
		    datout(dir,grid,mcs,out,rtim,rscale[0]);
		if(dmc[1]*out[1]!=rtim[0]*1000){ // Create graphics file
		    fileout(dir,grid,mcs);
		    if(mode[MTEM]==0)
			tproout(dir,fdmt[0],mcs);
		}
	    }
	}

	if(n_rank==0){
	    end=clock();
	    printf("\nExecution time: %d hour %d min %d sec\n",(int)((end-start)/CLOCKS_PER_SEC)/3600,(int)(((end-start)/CLOCKS_PER_SEC)%3600)/60,(int)(((end-start)/CLOCKS_PER_SEC)%3600)%60);
	    printf("== Simulation is done ==\n");
	}

	free4d(mpixlen[2]+2,dir[1],dir[2],mpigrid);
	for(i=2;i>=0;i--)
	    mpitgrid[i]=free3d(mpixlen[2]+2,dir[1],dir[2],mpitgrid[i]);
	free(mpitgrid);
	if(n_rank==0){
	    for(i=2;i>=0;i--)
		fdmt[i]=free3d(dir[0],dir[1],dir[2],fdmt[i]);
	}
	free(fdmt);

	if(mode[MLAT]==0)
	    melt[3]=-1;
	else if(melt[3]<0)
	    melt[3]=0;

	return;
}

int**** MPI_distr(int**** grid,int**** mpigrid,int mpix0,int mpixlen[],int dir[],int xpbc){
	int i,j,k,l,rank;
	int *data;

	// actual MC trial in each cores // Only the size of each part + boundaries
	// Actual MC: 1 <= i <= mpixlen[2] // i==0 & i=mpixlen[0] for boundaries
	if(n_rank==0){
	    data=(int*)calloc((mpixlen[2]+2)*dir[1]*dir[2]*NINFO,sizeof(int));
	    mpigrid=(int****)malloc((mpixlen[2]+2)*sizeof(int***));
	    for(i=mpixlen[2];i>0;i--){ // mpixl0 <= i <= dir[0]-1 
		mpigrid[i]=(int***)malloc(dir[1]*sizeof(int**));
		for(j=dir[1]-1;j>=0;j--){
		    mpigrid[i][j]=(int**)malloc(dir[2]*sizeof(int*));
		    for(k=dir[2]-1;k>=0;k--){
			mpigrid[i][j][k]=data+(((i*dir[1]*dir[2])+(j*dir[2])+k)*NINFO);
//			mpigrid[i][j][k]=(int*)calloc(3,sizeof(int));
			voxcpy(mpigrid[i][j][k],grid[mpix0+i-1][j][k]);
		    }
		}
	    }
// first & last arrays
	    i=mpixlen[2]+1; 
	    mpigrid[0]=(int***)malloc(dir[1]*sizeof(int**));
	    mpigrid[i]=(int***)malloc(dir[1]*sizeof(int**));
	    for(j=dir[1]-1;j>=0;j--){
		mpigrid[0][j]=(int**)malloc(dir[2]*sizeof(int*));
		mpigrid[i][j]=(int**)malloc(dir[2]*sizeof(int*));
		for(k=dir[2]-1;k>=0;k--){
		    mpigrid[0][j][k]=data+(((j*dir[2])+k)*NINFO);
		    mpigrid[i][j][k]=data+(((i*dir[1]*dir[2])+(j*dir[2])+k)*NINFO);
//		    mpigrid[i][j][k]=(int*)calloc(3,sizeof(int));
		}
	    }
// sending grid info. to other processes
	    for(rank=n_size-1;rank>0;rank--){ // rank
		l=mpixlen[0]*(n_size-rank-1); // mpix0 at each process
		MPI_Send(&(grid[l][0][0][0]),NINFO*mpixlen[0]*dir[1]*dir[2],MPI_INT,rank,rank,MPI_COMM_WORLD);
            }
	}else{
	    mpigrid=alloc4d(mpixlen[2]+2,dir[1],dir[2],mpigrid);
	    MPI_Recv(&(mpigrid[1][0][0][0]),NINFO*mpixlen[2]*dir[1]*dir[2],MPI_INT,0,n_rank,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	    // mpixlen[2] == mpixlen[0] except for rank==0
	}

// x dir PBC
	if(n_rank==0){ 
	    if(xpbc==0){ // X pbc off
		i=mpixlen[2]+1;
		for(j=dir[1]-1;j>=0;j--){
		    for(k=dir[2]-1;k>=0;k--)
			mpigrid[i][j][k][0]=MPISURF;
		}
	    }
	}else if(n_rank==n_size-1){
	    if(xpbc==0){ // X pbc off
		for(j=dir[1]-1;j>=0;j--){
		    for(k=dir[2]-1;k>=0;k--)
			mpigrid[0][j][k][0]=MPISURF;
		}
	    }	
	}

	MPI_BCsync(mpigrid,mpixlen[2],dir,xpbc,0);
	MPI_BCsync(mpigrid,mpixlen[2],dir,xpbc,1);
	return mpigrid;
}

void MPI_BCsync(int**** mpigrid,int mpixlen,int dir[],int xpbc,int mpimode){
	int i,j,k,o;

	MPI_Barrier(MPI_COMM_WORLD);
	if(mpimode==0){ // -x 방향 동기화
	    i=1; // min x
	    o=mpixlen+1; // max x +1
	    if(n_rank==0)
		MPI_Send(&(mpigrid[i][0][0][0]),NINFO*dir[1]*dir[2],MPI_INT,n_rank+1,n_rank,MPI_COMM_WORLD);
	    else if(n_rank!=n_size-1){
		MPI_Recv(&(mpigrid[o][0][0][0]),NINFO*dir[1]*dir[2],MPI_INT,n_rank-1,n_rank-1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Send(&(mpigrid[i][0][0][0]),NINFO*dir[1]*dir[2],MPI_INT,n_rank+1,n_rank,MPI_COMM_WORLD);
	    }else // n_rank==n_size-1
		MPI_Recv(&(mpigrid[o][0][0][0]),NINFO*dir[1]*dir[2],MPI_INT,n_rank-1,n_rank-1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

	    if(xpbc!=0){ // PBC ON on x direction
		 if(n_rank==0)
		    MPI_Recv(&(mpigrid[o][0][0][0]),NINFO*dir[1]*dir[2],MPI_INT,n_size-1,999,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		else if(n_rank==n_size-1)
		    MPI_Send(&(mpigrid[i][0][0][0]),NINFO*dir[1]*dir[2],MPI_INT,0,999,MPI_COMM_WORLD);
	    }
	}else{ // +x 방향 동기화
		i=mpixlen; // max x
		if(n_rank==0)
		    MPI_Recv(&(mpigrid[0][0][0][0]),NINFO*dir[1]*dir[2],MPI_INT,n_rank+1,n_rank+1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		else if(n_rank!=n_size-1){
		    MPI_Recv(&(mpigrid[0][0][0][0]),NINFO*dir[1]*dir[2],MPI_INT,n_rank+1,n_rank+1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		    MPI_Send(&(mpigrid[i][0][0][0]),NINFO*dir[1]*dir[2],MPI_INT,n_rank-1,n_rank,MPI_COMM_WORLD);
		}else // n_rank==n_size-1
		    MPI_Send(&(mpigrid[i][0][0][0]),NINFO*dir[1]*dir[2],MPI_INT,n_rank-1,n_rank,MPI_COMM_WORLD);

		if(xpbc!=0){ // PBC ON on x direction
		    if(n_rank==0)
			MPI_Send(&(mpigrid[i][0][0][0]),NINFO*dir[1]*dir[2],MPI_INT,n_size-1,999,MPI_COMM_WORLD);
		    else if(n_rank==n_size-1)
			MPI_Recv(&(mpigrid[0][0][0][0]),NINFO*dir[1]*dir[2],MPI_INT,0,999,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
	}
	return;
}

void MPI_gather(int**** grid,int**** mpigrid,int mpix0,int mpixlen[],int dir[]){
	int i,j,k,l,rank;

	MPI_Barrier(MPI_COMM_WORLD);
	if(n_rank==0){
	    for(i=mpixlen[2];i>0;i--){ // mpix0 <= i <= dir[0]-1, w/o boundaries
		for(j=dir[1]-1;j>=0;j--){
		    for(k=dir[2]-1;k>=0;k--)
			voxcpy(grid[mpix0+i-1][j][k],mpigrid[i][j][k]);
		}
	    }
	    for(rank=1;rank<n_size;rank++){ // rank
		l=mpixlen[0]*(n_size-rank-1); // mpix0 at each process
		MPI_Recv(&(grid[l][0][0][0]),NINFO*mpixlen[0]*dir[1]*dir[2],MPI_INT,rank,rank,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
	}else
	    MPI_Send(&(mpigrid[1][0][0][0]),NINFO*mpixlen[0]*dir[1]*dir[2],MPI_INT,0,n_rank,MPI_COMM_WORLD);

	return;
}

double**** MPI_Tdistr(double**** grid,double**** mpigrid,int mpix0,int mpixlen[],int dir[],int xpbc){
	int i,j,k,l,rank;
	double *data;

	// actual MC trial in each cores // Only the size of each part + boundaries
	// Actual MC: 1 <= i <= mpixlen[2] // i==0 & i=mpixlen[0] for boundaries
	if(n_rank==0){
	    mpigrid=(double****)calloc(3,sizeof(double***));
	    for(i=2;i>=0;i--)
		mpigrid[i]=alloc3d(mpixlen[2]+2,dir[1],dir[2],mpigrid[i]);

	    for(j=mpixlen[2];j>0;j--){ // mpixl0 <= i <= dir[0]-1 
		for(k=dir[1]-1;k>=0;k--){
		    for(i=dir[2]-1;i>=0;i--)
			mpigrid[0][j][k][i]=grid[0][mpix0+j-1][k][i];
		}
	    }
// sending grid info. to other processes
	    for(rank=n_size-1;rank>0;rank--){ // rank
		l=mpixlen[0]*(n_size-rank-1); // mpix0 at each process
		MPI_Send(&(grid[0][l][0][0]),mpixlen[0]*dir[1]*dir[2],MPI_DOUBLE,rank,rank,MPI_COMM_WORLD);
            }
	}else{
	    mpigrid=(double****)calloc(3,sizeof(double***));
	    for(i=2;i>=0;i--)
		mpigrid[i]=alloc3d(mpixlen[2]+2,dir[1],dir[2],mpigrid[i]);
	    MPI_Recv(&(mpigrid[0][1][0][0]),mpixlen[0]*dir[1]*dir[2],MPI_DOUBLE,0,n_rank,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}
// x dir PBC
	if(n_rank==0){ 
	    if(xpbc==0){ // X pbc off
		i=mpixlen[2]+1;
		for(j=dir[1]-1;j>=0;j--){
		    for(k=dir[2]-1;k>=0;k--)
			mpigrid[0][i][j][k]=-1;
		}
	    }
	}else if(n_rank==n_size-1){
	    if(xpbc==0){ // X pbc off
		for(j=dir[1]-1;j>=0;j--){
		    for(k=dir[2]-1;k>=0;k--)
			mpigrid[0][0][j][k]=-1;
		}
	    }	
	}
	MPI_TBCsync(mpigrid[0],mpixlen[2],dir,xpbc);

	return mpigrid;
}

void MPI_TBCsync(double*** mpigrid,int mpixlen,int dir[],int xpbc){
	int i,j,o;
	// 참고: 
	// rank N	| rank N-1		| ... | rank 0 
	// 0 ~ mpixlen-1, mpixlen ~ 2*mpixlen-1, ..., (~~~) ~ dir[0]-1
//	0=0
	i=1; // min x
	j=mpixlen; // max x
	o=mpixlen+1; // max x +1
	if(n_rank==0){
	    MPI_Send(&(mpigrid[i][0][0]),dir[1]*dir[2],MPI_DOUBLE,n_rank+1,n_rank,MPI_COMM_WORLD);

	    MPI_Recv(&(mpigrid[0][0][0]),dir[1]*dir[2],MPI_DOUBLE,n_rank+1,n_rank+1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}else if(n_rank!=n_size-1){
	    MPI_Recv(&(mpigrid[o][0][0]),dir[1]*dir[2],MPI_DOUBLE,n_rank-1,n_rank-1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	    MPI_Send(&(mpigrid[i][0][0]),dir[1]*dir[2],MPI_DOUBLE,n_rank+1,n_rank,MPI_COMM_WORLD);

	    MPI_Recv(&(mpigrid[0][0][0]),dir[1]*dir[2],MPI_DOUBLE,n_rank+1,n_rank+1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	    MPI_Send(&(mpigrid[j][0][0]),dir[1]*dir[2],MPI_DOUBLE,n_rank-1,n_rank,MPI_COMM_WORLD);
	}else{ // n_rank==n_size-1
	    MPI_Recv(&(mpigrid[o][0][0]),dir[1]*dir[2],MPI_DOUBLE,n_rank-1,n_rank-1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

	    MPI_Send(&(mpigrid[j][0][0]),dir[1]*dir[2],MPI_DOUBLE,n_rank-1,n_rank,MPI_COMM_WORLD);
	}
	if(xpbc!=0){ // PBC ON on x direction
	     if(n_rank==0){
		MPI_Recv(&(mpigrid[o][0][0]),dir[1]*dir[2],MPI_DOUBLE,n_rank-1,999,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

		MPI_Send(&(mpigrid[j][0][0]),dir[1]*dir[2],MPI_DOUBLE,n_size-1,999,MPI_COMM_WORLD);
	     }else if(n_rank==n_size-1){
		MPI_Send(&(mpigrid[i][0][0]),dir[1]*dir[2],MPI_DOUBLE,0,999,MPI_COMM_WORLD);

		MPI_Recv(&(mpigrid[0][0][0]),dir[1]*dir[2],MPI_DOUBLE,0,999,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	    }
	}
	return;
}


void MPI_Tsync(double*** fdmt,double*** mpigrid,int mpix0,int mpixlen[],int dir[]){
	int i,j,k,l,rank;

	if(n_rank==0){
	    for(i=mpixlen[2];i>0;i--){ // mpix0 <= i <= dir[0]-1, w/o boundaries
		for(j=dir[1]-1;j>=0;j--){
		    for(k=dir[2]-1;k>=0;k--)
			fdmt[mpix0+i-1][j][k]=mpigrid[i][j][k];
		}
	    }
	    for(rank=n_size-1;rank>0;rank--){ // rank
		l=mpixlen[0]*(n_size-rank-1); // mpix0 at each process
		MPI_Recv(&(fdmt[l][0][0]),mpixlen[0]*dir[1]*dir[2],MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
	}else
	    MPI_Send(&(mpigrid[1][0][0]),mpixlen[0]*dir[1]*dir[2],MPI_DOUBLE,0,n_rank,MPI_COMM_WORLD);
	return;
}

void heat_transfer_MPI(double**** tgrid,double tbc[][2],int mpidir[],int dir[],int pbc[],double mcdx,double const melt[],int fdm_loop,double k_per_cp,double Tmelt,int mpixlen){
	int i,j,k,l,temp[3],fdmax[3];//imax,jmax,kmax;
	double hconv,revdx2;
	
	fdm_loop++; // considering round down when (double -> int)
	if(fdm_loop==0 && n_rank==0){
	    printf("FDM ERROR: wrong FDM loop number %d\n  -> Automatically changed to 1\n",fdm_loop);
	    fdm_loop=1;
	}
	fdmax[0]=mpidir[0]-2; // 0, mpidir-1 voxels are only for BC sync
	fdmax[1]=mpidir[1]-1;
	fdmax[2]=mpidir[2]-1;

	revdx2=1/(mcdx*mcdx);

// Latent heat divided by loop number : retroactive accumulation of latent heat for 1 MCS
// Latent heat multiplied with phase transf. fraction : when dx_fdm != dx_mc, delta_fs != 1: delta_fs == V_mc / V_fdm
// FDM loop
	for(l=fdm_loop;l>0;l--){ // FDM time flow to 1 MCS
// Latent heat consideration
		for(i=fdmax[0];i>0;i--){ // Except for hypothetical boundary voxels
		    for(j=fdmax[1];j>=0;j--){
			for(k=fdmax[2];k>=0;k--){
			    if(tgrid[2][i][j][k]!=0)
				LatHeatDiffus(&(tgrid[0][i][j][k]),&(tgrid[2][i][j][k]),melt[6]);
			}
		    }
		}
//	0. inner part
		for(i=fdmax[0]-1;i>1;i--){
		    for(j=fdmax[1]-1;j>0;j--){
			for(k=fdmax[2]-1;k>0;k--)
			    tgrid[1][i][j][k]=tgrid[0][i][j][k]+melt[1]*k_per_cp*\
					     (revdx2*(tgrid[0][i+1][j][k]+tgrid[0][i-1][j][k]-2*tgrid[0][i][j][k])+\
					      revdx2*(tgrid[0][i][j+1][k]+tgrid[0][i][j-1][k]-2*tgrid[0][i][j][k])+\
					      revdx2*(tgrid[0][i][j][k+1]+tgrid[0][i][j][k-1]-2*tgrid[0][i][j][k]));
		    }
		}

// PBC treatment
// 1. X direction; for PBC, only for n_rank==0 & n_rank=n_size-1
		for(i=fdmax[0]-1;i>1;i--){
			j=fdmax[1];
			for(k=fdmax[2]-1;k>0;k--){ // bot & top xy planes
				tgrid[1][i][0][k]=tgrid[0][i][0][k]+melt[1]*k_per_cp*\
					   revdx2*(tgrid[0][i+1][0][k]+tgrid[0][i-1][0][k]-2*tgrid[0][i][0][k]);
				tgrid[1][i][j][k]=tgrid[0][i][j][k]+melt[1]*k_per_cp*\
					   revdx2*(tgrid[0][i+1][j][k]+tgrid[0][i-1][j][k]-2*tgrid[0][i][j][k]);
			}
			k=fdmax[2];
			for(j=fdmax[1];j>=0;j--){ // bot & top xz planes, including edges
				tgrid[1][i][j][0]=tgrid[0][i][j][0]+melt[1]*k_per_cp*\
					   revdx2*(tgrid[0][i+1][j][0]+tgrid[0][i-1][j][0]-2*tgrid[0][i][j][0]);
				tgrid[1][i][j][k]=tgrid[0][i][j][k]+melt[1]*k_per_cp*\
					   revdx2*(tgrid[0][i+1][j][k]+tgrid[0][i-1][j][k]-2*tgrid[0][i][j][k]);
			}
		}

		i=fdmax[0];
		if(n_rank==n_size-1){ // Lower X plane
		    for(j=fdmax[1];j>=0;j--){ // upper: mpi BC
			for(k=fdmax[2];k>=0;k--)
			    tgrid[1][i][j][k]=tgrid[0][i][j][k]+melt[1]*k_per_cp*revdx2*\
					     (tgrid[0][i+1][j][k]+tgrid[0][i-1][j][k]-2*tgrid[0][i][j][k]);
		    }

		    if(pbc[0]!=0){// X pbc on
			for(j=fdmax[1];j>=0;j--){
			    for(k=fdmax[2];k>=0;k--)
				tgrid[1][1][j][k]=tgrid[0][1][j][k]+melt[1]*k_per_cp*revdx2*\
						  (tgrid[0][0][j][k]+tgrid[0][2][j][k]-2*tgrid[0][1][j][k]);
			}

		    }else{ // X pbc off
// Lower X plane
			if(tbc[0][0]==0){ // adiabatic
			    for(j=fdmax[1];j>=0;j--){
				for(k=fdmax[2];k>=0;k--)
				  tgrid[1][1][j][k]=tgrid[0][1][j][k]+melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][2][j][k]-tgrid[0][1][j][k]);
			    }
			}else if(tbc[0][0]<0){ // convection
			    hconv=tbc[0][1]*k_per_cp*melt[1]/(mcdx*melt[3]); // h k dt / (cp dx k) = h dt / (cp dx)
			    for(j=fdmax[1];j>=0;j--){
				for(k=fdmax[2];k>=0;k--)
				  tgrid[1][1][j][k]=tgrid[0][1][j][k]+melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][2][j][k]-tgrid[0][1][j][k])-\
						    hconv*(tbc[0][0]+tgrid[0][1][j][k]); // - (-Ta+Tnow)
			    }
			}else{ // tbc > 0 : heat sink
			    for(j=fdmax[1];j>=0;j--){
				for(k=fdmax[2];k>=0;k--)
				  tgrid[1][1][j][k]=tgrid[0][1][j][k]+melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][2][j][k]+tbc[0][0]-2*tgrid[0][1][j][k]);
			    }
			}
		    }
		}else if(n_rank==0){ // upper X plane
		    for(j=fdmax[1];j>=0;j--){ // lower: mpi BC
			for(k=fdmax[2];k>=0;k--)
			    tgrid[1][1][j][k]=tgrid[0][1][j][k]+melt[1]*k_per_cp*revdx2*\
					     (tgrid[0][0][j][k]+tgrid[0][2][j][k]-2*tgrid[0][1][j][k]);
		    }
		    if(pbc[0]!=0){// X pbc on
			for(j=fdmax[1];j>=0;j--){
			    for(k=fdmax[2];k>=0;k--) // [0][j][k] = [fdmax[0]][j][k]
				tgrid[1][i][j][k]=tgrid[0][i][j][k]+melt[1]*k_per_cp*revdx2*\
						  (tgrid[0][i+1][j][k]+tgrid[0][i-1][j][k]-2*tgrid[0][i][j][k]);
			}
		    }else{//X pbc off
			if(tbc[1][0]==0){ // adiabatic
			    for(j=fdmax[1];j>=0;j--){
				for(k=fdmax[2];k>=0;k--)
				  tgrid[1][i][j][k]=tgrid[0][i][j][k]+melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i-1][j][k]-tgrid[0][i][j][k]);
			    }
			}else if(tbc[1][0]<0){ // convection
			    hconv=tbc[1][1]*k_per_cp*melt[1]/(mcdx*melt[3]); // h k dt / (cp dx k) = h dt / (cp dx)
			    for(j=fdmax[1];j>=0;j--){
				for(k=fdmax[2];k>=0;k--)
				  tgrid[1][i][j][k]=tgrid[0][i][j][k]+melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i-1][j][k]-tgrid[0][i][j][k])-\
						    hconv*(tbc[1][0]+tgrid[0][i][j][k]); // - (-Ta+Tnow)
			    }
			}else{ // tbc > 0 : heat sink
			    for(j=fdmax[1];j>=0;j--){
				for(k=fdmax[2];k>=0;k--)
				  tgrid[1][i][j][k]=tgrid[0][i][j][k]+melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i-1][j][k]+tbc[1][0]-2*tgrid[0][i][j][k]);
			    }
			}
		    }
		}else{ // inner MPI grids
		    for(j=fdmax[1];j>=0;j--){
			for(k=fdmax[2];k>=0;k--){
			    tgrid[1][1][j][k]=tgrid[0][1][j][k]+melt[1]*k_per_cp*revdx2*\
					     (tgrid[0][0][j][k]+tgrid[0][2][j][k]-2*tgrid[0][1][j][k]);
			    tgrid[1][i][j][k]=tgrid[0][i][j][k]+melt[1]*k_per_cp*revdx2*\
					     (tgrid[0][i+1][j][k]+tgrid[0][i-1][j][k]-2*tgrid[0][i][j][k]);
			}
		    }
		}

// 2. Y direction
		for(j=fdmax[1]-1;j>0;j--){
		    i=fdmax[0];
		    for(k=fdmax[2]-1;k>0;k--){ // bot & top xz planes
			tgrid[1][1][j][k]+=melt[1]*k_per_cp*\
					   revdx2*(tgrid[0][1][j+1][k]+tgrid[0][1][j-1][k]-2*tgrid[0][1][j][k]);
			tgrid[1][i][j][k]+=melt[1]*k_per_cp*\
					   revdx2*(tgrid[0][i][j+1][k]+tgrid[0][i][j-1][k]-2*tgrid[0][i][j][k]);
		    }
		    k=fdmax[2];
		    for(i=fdmax[0];i>0;i--){ // bot & top xz planes, including edges
			tgrid[1][i][j][0]+=melt[1]*k_per_cp*\
					   revdx2*(tgrid[0][i][j+1][0]+tgrid[0][i][j-1][0]-2*tgrid[0][i][j][0]);
			tgrid[1][i][j][k]+=melt[1]*k_per_cp*\
					   revdx2*(tgrid[0][i][j+1][k]+tgrid[0][i][j-1][k]-2*tgrid[0][i][j][k]);
		    }
		}

		j=fdmax[1];
 		if(pbc[1]!=0){// Y pbc on
			for(i=fdmax[0];i>0;i--){ // [i][fdmax][k] = [i][0][k]
				for(k=fdmax[2];k>=0;k--){
					tgrid[1][i][0][k]+=melt[1]*k_per_cp*revdx2*\
							  (tgrid[0][i][1][k]+tgrid[0][i][j][k]-2*tgrid[0][i][0][k]);
					tgrid[1][i][j][k]+=melt[1]*k_per_cp*revdx2*\
							  (tgrid[0][i][0][k]+tgrid[0][i][j-1][k]-2*tgrid[0][i][j][k]);
				}
			}
		}else{
// Lower Y plane
			if(tbc[2][0]==0){ // adiabatic
			    for(i=fdmax[0];i=0;i--){
				for(k=fdmax[2];k>=0;k--)
				  tgrid[1][i][0][k]+=melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i][1][k]-tgrid[0][i][0][k]);
			    }
			}else if(tbc[2][0]<0){ // convection
			    hconv=tbc[2][1]*k_per_cp*melt[1]/(mcdx*melt[3]); // h k dt / (cp dx k) = h dt / (cp dx)
			    for(i=fdmax[0];i>0;i--){
				for(k=fdmax[2];k>=0;k--)
				  tgrid[1][i][0][k]+=melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i][1][k]-tgrid[0][i][0][k])-\
						    hconv*(tbc[2][0]+tgrid[0][i][0][k]); // - (-Ta+Tnow) = Ta-Tnow
			    }
			}else{ // tbc > 0 : heat sink
			    for(i=fdmax[0];i>0;i--){
				for(k=fdmax[2];k>=0;k--)
				  tgrid[1][i][0][k]+=melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i][1][k]+tbc[2][0]-2*tgrid[0][i][0][k]);
			    }
			}
// Upper Y plane
			if(tbc[3][0]==0){ // adiabatic
			    for(i=fdmax[0];i>0;i--){
				for(k=fdmax[2];k>=0;k--)
				  tgrid[1][i][j][k]+=melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i][j-1][k]-tgrid[0][i][j][k]);
			    }
			}else if(tbc[3][0]<0){ // convection
			    hconv=tbc[3][1]*k_per_cp*melt[1]/(mcdx*melt[3]); // h k dt / (cp dx k) = h dt / (cp dx)
			    for(i=fdmax[0];i>0;i--){
				for(k=fdmax[2];k>=0;k--)
				  tgrid[1][i][j][k]+=melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i][j-1][k]-tgrid[0][i][j][k])-\
						    hconv*(tbc[3][0]+tgrid[0][i][j][k]); // - (-Ta+Tnow)
			    }
			}else{ // tbc > 0 : heat sink
			    for(i=fdmax[0];i>0;i--){
				for(k=fdmax[2];k>=0;k--)
				  tgrid[1][i][j][k]+=melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i][j-1][k]+tbc[3][0]-2*tgrid[0][i][j][k]);
			    }
			}
		}

// 3. Z direction
		for(k=fdmax[2]-1;k>0;k--){
			i=fdmax[0];
			for(j=fdmax[1];j>=0;j--){ // bot & top yz planes
				tgrid[1][0][j][k]+=melt[1]*k_per_cp*\
					   revdx2*(tgrid[0][0][j][k+1]+tgrid[0][0][j][k-1]-2*tgrid[0][0][j][k]);
				tgrid[1][i][j][k]+=melt[1]*k_per_cp*\
					   revdx2*(tgrid[0][i][j][k+1]+tgrid[0][i][j][k-1]-2*tgrid[0][i][j][k]);
			}
			j=fdmax[1];
			for(i=fdmax[0]-1;i>1;i--){ // bot & top xz planes
				tgrid[1][i][0][k]+=melt[1]*k_per_cp*\
					   revdx2*(tgrid[0][i][0][k+1]+tgrid[0][i][0][k-1]-2*tgrid[0][i][0][k]);
				tgrid[1][i][j][k]+=melt[1]*k_per_cp*\
					   revdx2*(tgrid[0][i][j][k+1]+tgrid[0][i][j][k-1]-2*tgrid[0][i][j][k]);
			}
		}

		k=fdmax[2];
 		if(pbc[2]!=0){// Z pbc on
		    for(i=fdmax[0];i>0;i--){
			for(j=fdmax[1];j>=0;j--){ // [i][k][0] = [i][j][fdmax]
				tgrid[1][i][j][0]+=melt[1]*k_per_cp*revdx2*\
						  (tgrid[0][i][j][1]+tgrid[0][i][j][k]-2*tgrid[0][i][j][0]);
				tgrid[1][i][j][k]+=melt[1]*k_per_cp*revdx2*\
						  (tgrid[0][i][j][0]+tgrid[0][i][j][k-1]-2*tgrid[0][i][j][k]);
			}
		    }
		}else{
// Lower Z plane
		    if(tbc[4][0]==0){ // adiabatic
			for(i=fdmax[0];i>0;i--){
			    for(j=fdmax[1];j>=0;j--)
				tgrid[1][i][j][0]+=melt[1]*k_per_cp*revdx2*\
						  (tgrid[0][i][j][1]-tgrid[0][i][j][0]);
			}
		    }else if(tbc[4][0]<0){ // convection
			hconv=tbc[4][1]*k_per_cp*melt[1]/(mcdx*melt[3]); // h k dt / (cp dx k) = h dt / (cp dx)
			for(i=fdmax[0];i>0;i--){
			    for(j=fdmax[1];j>=0;j--)
				tgrid[1][i][j][0]+=melt[1]*k_per_cp*revdx2*\
						   (tgrid[0][i][j][1]-tgrid[0][i][j][0])-\
						    hconv*(tbc[4][0]+tgrid[0][i][j][0]); // - (-Ta+Tnow) = Ta-Tnow
			}
		    }else{ // tbc > 0 : heat sink
			for(i=fdmax[0];i>0;i--){
			    for(j=fdmax[1];j>=0;j--)
				  tgrid[1][i][j][0]+=melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i][j][1]+tbc[4][0]-2*tgrid[0][i][j][0]);
			}
		    }
// Upper Z plane
		    if(tbc[5][0]==0){ // adiabatic
			for(i=fdmax[0];i>0;i--){
			    for(j=fdmax[1];j>=0;j--)
				  tgrid[1][i][j][k]+=melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i][j][k-1]-tgrid[0][i][j][k]);
			}
		    }else if(tbc[5][0]<0){ // convection
			hconv=tbc[5][1]*k_per_cp*melt[1]/(mcdx*melt[3]); // h k dt / (cp dx k) = h dt / (cp dx)
			for(i=fdmax[0];i>0;i--){
			    for(j=fdmax[1];j>=0;j--)
				  tgrid[1][i][j][k]+=melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i][j][k-1]-tgrid[0][i][j][k])-\
						    hconv*(tbc[5][0]+tgrid[0][i][j][k]); // - (-Ta+Tnow)
			}
		    }else{ // tbc > 0 : heat sink
			for(i=fdmax[0];i>0;i--){
			    for(j=fdmax[1];j>=0;j--)
				  tgrid[1][i][j][k]+=melt[1]*k_per_cp*revdx2*\
						    (tgrid[0][i][j][k-1]+tbc[5][0]-2*tgrid[0][i][j][k]);
			}
		    }
		}
// next -> now
		for(i=fdmax[0];i>0;i--){
		    for(j=fdmax[1];j>=0;j--){
			for(k=fdmax[2];k>=0;k--)
				tgrid[0][i][j][k]=tgrid[1][i][j][k];
		    }
		}
		MPI_TBCsync(tgrid[0],mpixlen,dir,pbc[0]);
//		MPI_TBCsync(tgrid[2],mpixlen,dir,pbc[0]);
	} // FDM 1 iteration

	return;
}

int MPI_MeltPoolCheck(double mode,double direction,double mp[],double radius,double xpos,int mpixlen,int mpidir){
	int xmin,xmax;
	if(mode==0){ // TD shape
	    if(direction<0){ // x direction movement
		xmin=(int)(xpos-mp[1]);
		xmax=(int)(xpos+mp[1]);
	    }else{ // y direction movement
		xmin=(int)(xpos-0.7698*mp[2]);; //maximum of teardrop(mp[2],t,1);
		xmax=(int)(xpos+0.7698*mp[2]);
	    }
	}else{ // circular shape
	    xmin=(int)(xpos-radius);
	    xmax=(int)(xpos+radius);
	}
//	return !(xmax<mpix0||xmin>mpix0+mpixlen);
	return !(xmax<=0 || xmin>mpixlen);
}

int find_min_rank(double rtim){
	double local_pair[2],global_pair[2];
	local_pair[0]=rtim;
	local_pair[1]=(double)n_rank;

	MPI_Allreduce(local_pair,global_pair,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);

	return (int)global_pair[1];
}


