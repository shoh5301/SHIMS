#include "Single.h"

void GGmodule(int**** grid,double jc[],int dir[],int out[],char ttem0[],double melt[],int pbc[],double const sang[],double rscale[],double tbc[][2]){// Grain growth
	clock_t start,end;
	char tcheck[2],temp[100],Tfun[50]; //Tfun: postfix equation for T-t rel.
	int i,j,k,l,o,mcs=0,total=0,step=0,stepdiv[2],check=2,dmc[3]={1};
	FILE* log=NULL;
	double realt,rtim[5]={0},realE,pmob,Tmelt,factor;

	if(out[3]<=0){ // Time criteria is MCS
		dmc[0]=out[4]; // data print
		dmc[1]=out[1]; // map print
	}
	printf("\n####### MODULE : GG #########\n");
// Grian growth only
	factor=melt[1];
	if(rscale[0]>0)
		rtim[3]=rscale[1]*1E3;//(m to mm)
	datout(dir,grid,0,out,rtim,rscale[0]); // Write data file & get L0
	rtim[2]=0;

	Tmelt=rscale[2];
	total=dir[0]*dir[1]*dir[2];
	tok(ttem0,tcheck); //tcheck[0]: 0 for MCS, 1 for time(s)
	rtim[4]=-1; // junk value for MCS=0
	if(rscale[0]!=0)
	    rtim[1]=rscale[1]*rscale[1]*rscale[7]/rscale[5]; // Time is independent to temperature!
		//     lamda^2     *  K0_mc  /  K0_e

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
	printf("    Single core mode\n\n");

	start=clock();
	srand(time(NULL));

	while(1){//	srand(time(NULL));
		// Temperature & real time setting
		if(rtim[4]!=0&&rtim[0]>rtim[4]){
			tok(ttem0,temp);
			rtim[4]=atof(temp); //rtim[4] is time range for T func.
			tok(ttem0,temp);
			change_postfix(temp,Tfun); // Tfun is postfix eq.
			printf(" # New T(t) function accepted: T(t) = %s",temp);
			if(tcheck[0]=='0') // MCS variable
				printf(" after t = %d MCS\n",mcs);
			else
				printf(" after t = %.2f sec\n",rtim[0]);
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
				printf("@@@ WARNING: Higher temperature then melting temperature: T = %.1f K vs. Tm = %.1f K ... @@@\n",realt,melt[0]);
			realE=rscale[1]*rscale[1]/(kB*realt);
//			     dx^2 / (K T)
			pmob=factor*exp(-rscale[3]*(1/realt-1/Tmelt)*REVRGAS);
//					Q  * (1/T-1/Tm)  / R
		}
		mcs++;
		step=total;
		while(step!=0){	// MC trial ////////////////////////////////////////////////////////////////////////////////
			i=(int)(rand()/(RAND_MAX/dir[0]+1));
			j=(int)(rand()/(RAND_MAX/dir[1]+1));
			k=(int)(rand()/(RAND_MAX/dir[2]+1));

			MCtrial_GG(i,j,k,grid,dir,pbc,jc,tbc,pmob,realE);
			step--;	// 1 Attempt
		} // 1 MCS
		if(rscale[0]>0)
			rtim[0]=rtim[0]+rtim[1];
		else
			rtim[0]=mcs;
		if(out[3]<=0){ // Time criteria is MCS
		    dmc[0]--;
		    dmc[1]--;
		    if(dmc[0]==0){ // Data out
			printf(" %dth MC Step was taken...\n",mcs);
			datout(dir,grid,mcs,out,rtim,rscale[0]); // Write data file
			dmc[0]=out[4];
		    }
		    if(dmc[1]==0){ // Create graphics file
			fileout(dir,grid,mcs);
			dmc[1]=out[1];
		    }
		    if(mcs==out[0])
			break;
		}else{ // Time criteria is sec
		    if(dmc[0]*out[4]<=rtim[0]*1000){ // Data out
			printf(" %.5f sec is passed... where T = %.1f\n",rtim[0],realt);
			datout(dir,grid,mcs,out,rtim,rscale[0]); // Write data file
			dmc[0]++;
//			count=0;
		    }
		    if(dmc[1]*out[1]<=rtim[0]*1000){ // Create graphics file
			fileout(dir,grid,mcs);
			dmc[1]++;
		    }
		    if(rtim[0]*1000>=out[0])
			break;
		}
	} // End a MC trial

	if(out[3]<=0){ // Time criteria is MCS
		if(dmc[0]!=out[4]) // Check unsaved MCSs
		    datout(dir,grid,mcs,out,rtim,rscale[0]);
		if(dmc[1]!=out[1])
		    fileout(dir,grid,mcs);
	}else{ // Time criteria is sec
		if(dmc[0]*out[4]!=rtim[0]*1000){
			fileout(dir,grid,mcs);
			datout(dir,grid,mcs,out,rtim,rscale[0]);
		}
	}

	end=clock();
	printf("\nExecution time: %d hour %d min %d sec\n",(int)((end-start)/CLOCKS_PER_SEC)/3600,(int)(((end-start)/CLOCKS_PER_SEC)%3600)/60,(int)(((end-start)/CLOCKS_PER_SEC)%3600)%60);
	printf("== Simulation is done ==\n");

	return;
}

void SOLmodule(int**** grid,double jc[],int dir[],int out[],double melt[],int pbc[],double rscale[],double tbc[][2],char dfeq[]){
	clock_t start,end,ttime;
	char temp[100],Tfun[50];
	int i,j,k,l,m,n,o,situ;
	int mcs=0,total=0,step=0,dmc[3]={1},fdm_loop=0;
	int mode[4]={1,1,0,0};
	double k_per_cp,Tmax[3],Tinfo[2]={0},dtsave[2]; 
	double Vsite,Asite,fnucl,factor[2]; 
	double pcps[3]; 
	double rtim[5]={0},dEif,prob,dGfor,dGfus,dtdMCS[2]={0};
	double Tnow,dHfus[2],fhetero_s,Nsite[2];
	double ****fdmt=NULL;
	FILE* ftemp=NULL;

        int**** mpigrid=NULL;
        int** mpicount=NULL; // yz plane: duplicated cell checker
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
	melt[6]=LatHeatPerStep(dHfus[1],melt[1],rscale[1],melt[3],rscale[9],melt[4]);

	//      k    *dt / Cp    * dx^2
	Tmax[0]=rscale[2]; // Tliq
	Tmax[1]=Tmax[0]-rscale[14]; // Critical supercooling
	Tmax[2]=rscale[2]-melt[0]; // Tsol

	// for solidification
	pcps[0]=Asite*(26*(jc[1]-jc[4])); // heterog. nucleation surrounded by other grains: maximum pc value for "a solid nucleus"
	if(dfeq[0]=='d')
		pcps[1]=Vsite*dHfus[0]*(1-KAUZMANN); // 1-Tk/Tm where Tk is Kauzmann Temperature for glass transition
	else
		pcps[1]=Vsite*Xl_df(dfeq,KAUZMANN*melt[8],dHfus[0],Tmax,melt[8]); //dHfus[0]*(1-KAUZMANN); // 1-Tk/Tm where Tk is Kauzmann Temperature for glass transition
	pcps[2]=1/(pcps[0]+pcps[1]);

	if(out[3]<=0){ // Time criteria is MCS
	    dmc[0]=out[4]; // data print
	    dmc[1]=out[1]; // map print
	}
	total=dir[0]*dir[1]*dir[2];

	printf("\n####### MODULE : SOL #########\n");
	rtim[3]=rscale[1];//(m to mm)
	datout(dir,grid,0,out,rtim,rscale[0]); // Write data file

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

	if(melt[2]<0){
	    // FDM T profile grid allocation & T0 setting
	    for(i=1;i<=2;i++)
		fdmt[i]=alloc3d(dir[0],dir[1],dir[2],fdmt[i]);

	}else if(melt[2]==0){ //frozen temperature condition
	    mode[MTEM]=1; //frozen temperature approximation mode on
//	T*-G*
	    mode[MLAT]=0;
	}else{
	    fdmt=(double****)calloc(3,sizeof(double***));
	    for(i=2;i>=0;i--)
		fdmt[i]=alloc3d(dir[0],dir[1],dir[2],fdmt[i]);
	    for(i=dir[0]-1;i>=0;i--){
		for(j=dir[1]-1;j>=0;j--){
			for(k=dir[2]-1;k>=0;k--)
			    fdmt[0][i][j][k]=melt[2];
		}
	    }
	}

	if(mode[MTEM]==0)
	    tproout(dir,fdmt[0],0);

///////// dt for solidification & melting?/////////////////////////////////////////
	dtdMCS[0]=rscale[1]*rscale[1]*rscale[7]/rscale[5];	// For grain growth
	dtdMCS[1]=rscale[1]*rscale[10];	// dt for solid nucleation & growth
//		   dx * (K_MC / K_e)
// Dynamic MCS setting
	out[5]=3-out[5];
// Now, out[5] 0 -> = t_max, 1 -> t_small, 2 -> t_gg

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
	    melt[6]=LatHeatPerStep(dHfus[1],melt[1],rscale[1],melt[3],rscale[9],melt[4]);
	    dtsave[2]=rtim[1];
	}else{
//	    if(rtim[1]>rscale[12]){
	    factor[0]=factor[0]*rscale[12]/rtim[1];
	    factor[1]=factor[1]*rscale[12]/rtim[1];
	    rtim[1]=rscale[12];
//	    }
	}

//////////////////////////////////////////////////////////////////////////////////

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
		melt[6]=LatHeatPerStep(dHfus[1],melt[1],rscale[1],melt[3],rscale[9],melt[4]);

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
	if(out[3]<0)
		puts("    # Simulation will automatically stop when there is no liquid... #");
	printf("    Single core mode\n\n");

	rtim[4]=-1; // junk value for MCS=0
	dtsave[1]=rtim[1];

	start=clock();
	srand(time(NULL));
/////////////////////////////// MC trial & FDM /////////////////////////////////////////////////

	while(1){
	    mcs++;
	    step=total;
	    while(step!=0){
		i=(int)(rand()/(RAND_MAX/dir[0]+1));
		j=(int)(rand()/(RAND_MAX/dir[1]+1));
		k=(int)(rand()/(RAND_MAX/dir[2]+1));
		MCtrial_SOL(i,j,k,grid,fdmt,dir,pbc,mode,melt,Tmax,dHfus,jc,tbc,Vsite,Asite,factor,Nsite,pcps,rscale[3],dfeq);
		step--;     // 1 Attempt
	    }

	    if(mode[MLAT]!=0 && mode[MTEM]==0)
		heat_transfer(grid,fdmt,tbc,dir,pbc,rscale[1],melt,fdm_loop,k_per_cp,Tmax[0],mode[MPRO]); // Loop inside the function
	    else{
		if(mode[MPRO]>=0)
		    heat_transfer(grid,fdmt,tbc,dir,pbc,rscale[1],melt,fdm_loop,k_per_cp,Tmax[0],mode[MPRO]);
		else{
		    if(latent_heat_check(fdmt[2],dir)!=0) // Latent heat is not fully consumed 
		        heat_transfer(grid,fdmt,tbc,dir,pbc,rscale[1],melt,fdm_loop,k_per_cp,Tmax[0],mode[MPRO]);
		    // else, isothermal condition
		}
	    }

	    if(out[5]==0 && mode[MLAT]!=0 && mode[MTEM]==0){ // Active dt determination during simulation
		dmc[2]--;
		if(dmc[2]==0){ // time OPT
		    dmc[2]=out[6];
		    rtim[1]=fixed_dt_set(mcs,&(mode[MPRO]),&fdm_loop,dtsave,dtdMCS,dHfus,Tinfo,factor,Vsite,Asite,out,Tmax,pcps,grid,fdmt,dir,pbc,rscale,melt);
		}
		if(rtim[1]<0){
		    printf(" @@@ Liquid fraction is 0, solidification is DONE @@@\n");
		    break;
		}
	    }

	    if(rscale[0]!=0)
		rtim[0]=rtim[0]+rtim[1];
	    else
		rtim[0]=mcs;
	    if(out[3]<=0){ // Time criteria is MCS
		dmc[0]--;
		dmc[1]--;
		if(dmc[0]==0){ // Data out
		    printf(" %dth MC Step was taken...\n",mcs);
		    datout(dir,grid,mcs,out,rtim,rscale[0]); // Write data file
		    dmc[0]=out[4];
		}
		if(dmc[1]==0){ // Create graphics file
		    fileout(dir,grid,mcs);
		    if(mode[MTEM]==0)
			tproout(dir,fdmt[0],mcs);
		    dmc[1]=out[1];
		}
		if(mcs==out[0])
		    break;
	    }else{ // Time criteria is sec
		if(dmc[0]*out[4]<=rtim[0]*1000){ // Data out
		    printf(" %.5f sec is passed...\n",rtim[0]);
		    datout(dir,grid,mcs,out,rtim,rscale[0]); // Write data file
		    dmc[0]++;
		}
		if(dmc[1]*out[1]<=rtim[0]*1000){ // Create graphics file
		    fileout(dir,grid,mcs);
		    if(mode[MTEM]==0)
			tproout(dir,fdmt[0],mcs);
		    dmc[1]++;
		}
		if(rtim[0]*1000>=out[0])
		    break;
	    }
	} // End a MCS

	if(out[3]<=0){ // Time criteria is MCS
	    if(out[0]!=mcs)
		out[0]=mcs;

		if(dmc[0]!=out[4]) // Check unsaved MCSs
			datout(dir,grid,mcs,out,rtim,rscale[0]);
		if(dmc[1]!=out[1]){
			fileout(dir,grid,mcs);
			if(mode[MTEM]==0)
			    tproout(dir,fdmt[0],mcs);
		}
	}else{ // Time criteria is sec
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

	end=clock();
	printf("\nExecution time: %d hour %d min %d sec\n",(int)((end-start)/CLOCKS_PER_SEC)/3600,(int)(((end-start)/CLOCKS_PER_SEC)%3600)/60,(int)(((end-start)/CLOCKS_PER_SEC)%3600)%60);
	printf("== Simulation is done ==\n");

	if(mode[MTEM]==0){
	    for(i=2;i>=0;i--)
		fdmt[i]=free3d(dir[0],dir[1],dir[2],fdmt[i]);
	    free(fdmt);
	}

	if(mode[MLAT]==0)
	    melt[3]=-1;
	else if(melt[3]<0)
	    melt[3]=0;

	return;
}

void AMmodule(int**** grid,double jc[],int dir[],int out[],double melt[],int pbc[],double rscale[],double tbc[][2],double am[],char dfeq[]){
	clock_t start,end,ttime;
	char temp[100],Tfun[50]; //Tfun: postfix equation for T-t rel.
	int i,j,k,l,m,n,o,situ;
	int mcs=0,total=0,step=0,dmc[3]={1},fdm_loop=0;
	int mode[4]={1,1,0,0};
	double k_per_cp,Tmax[3],Tinfo[2]={0},dtsave[2],Pmax[2];
	double Vsite,Asite,fnucl,factor[2];
	double pcps[3];
	double rtim[5]={0},dEif,prob,dGfor,dGfus,dtdMCS[2]={0};
	double Tnow,dHfus[2],fhetero_s,Nsite[2],mp[5],pos[2];
	double ****fdmt=NULL;
	FILE* ftemp=NULL;

	printf("\n####### MODULE : AM #########\n");
	rtim[3]=rscale[1];//(m to mm)
	datout(dir,grid,0,out,rtim,rscale[0]); // Write data file
// Initial setting is similar to SOL Module

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
	melt[6]=dHfus[1]*(melt[1]/melt[6]);*/
	dtsave[0]=melt[1]; // original dt for FDM
	melt[6]=LatHeatPerStep(dHfus[1],melt[1],rscale[1],melt[3],rscale[9],melt[4]);

	//      k    *dt / Cp    * dx^2
	Tmax[0]=rscale[2];
	Tmax[1]=Tmax[0]-rscale[14]; // Critical supercooling
	Tmax[2]=rscale[2]-melt[0]; // Tsol

	// for solidification
	pcps[0]=Asite*(26*(jc[1]-jc[4])); // heterog. nucleation surrounded by other grains: maximum pc value for "a solid nucleus"
	if(dfeq[0]=='d')
		pcps[1]=Vsite*dHfus[0]*(1-KAUZMANN); // 1-Tk/Tm where Tk is Kauzmann Temperature for glass transition
	else
		pcps[1]=Vsite*Xl_df(dfeq,KAUZMANN*melt[8],dHfus[0],Tmax,melt[8]); //dHfus[0]*(1-KAUZMANN); // 1-Tk/Tm where Tk is Kauzmann Temperature for glass transition
	pcps[2]=1/(pcps[0]+pcps[1]);

	if(out[3]<=0){ // Time criteria is MCS
	    dmc[0]=out[4]; // data print
	    dmc[1]=out[1]; // map print
	}
	total=dir[0]*dir[1]*dir[2];

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

	if(melt[2]<0){
	    // FDM T profile grid allocation & T0 setting
	    for(i=2;i>0;i--)
		fdmt[i]=alloc3d(dir[0],dir[1],dir[2],fdmt[i]);
	}else{
	    fdmt=(double****)calloc(3,sizeof(double***));
	    for(i=2;i>=0;i--)
		fdmt[i]=alloc3d(dir[0],dir[1],dir[2],fdmt[i]);
	    for(i=dir[0]-1;i>=0;i--){
		for(j=dir[1]-1;j>=0;j--){
		    for(k=dir[2]-1;k>=0;k--)
			fdmt[0][i][j][k]=melt[2];
	    	}
	    }
	}

	if(mode[MLAT]!=0)
	    tproout(dir,fdmt[0],0);

///////// dt for solidification & melting?/////////////////////////////////////////
	dtdMCS[0]=rscale[1]*rscale[1]*rscale[7]/rscale[5];	// For grain growth
	dtdMCS[1]=rscale[1]*rscale[10];	// dt for solid nucleation & growth
//		   dx * (K_MC / K_e)
// Dynamic MCS setting
	out[5]=3-out[5];
// Now, out[5] 0 -> = t_max, 1 -> t_small, 2 -> t_gg

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
	    melt[6]=LatHeatPerStep(dHfus[1],melt[1],rscale[1],melt[3],rscale[9],melt[4]);
	    dtsave[2]=rtim[1];
	}else{
//	    if(rtim[1]>rscale[12]){
	    factor[0]=factor[0]*rscale[12]/rtim[1];
	    factor[1]=factor[1]*rscale[12]/rtim[1];
	    rtim[1]=rscale[12];
//	    }
	}

	if(out[3]<=0){ // Time criteria is MCS
		dmc[0]=out[4]; // data print
		dmc[1]=out[1]; // map print
	}
	rtim[3]=rscale[1]*1E3;//(m to mm)

	datout(dir,grid,0,out,rtim,rscale[0]); // Write data file & get L0
	total=dir[0]*dir[1]*dir[2];

	fdmt=(double****)calloc(3,sizeof(double***));
	for(i=2;i>=0;i--)
	    fdmt[i]=alloc3d(dir[0],dir[1],dir[2],fdmt[i]);
	for(i=dir[0]-1;i>=0;i--){
		for(j=dir[1]-1;j>=0;j--){
			for(k=dir[2]-1;k>=0;k--)
			    fdmt[0][i][j][k]=melt[2];
	}	}
	tproout(dir,fdmt[0],0);

	printf(" << Monte Carlo Simulation Started >>\n");
	printf("    Sample size : %d X %d X %d\n",dir[0],dir[1],dir[2]);
	if(out[3]<=0) // MCS criteria
		printf("    Total MCS %d steps, Print Grain Map every %d steps, Print Result every %d steps \n",out[0],out[1],out[4]);
	else // sec criteria
		printf("    Total time %.3f sec, Print Grain Map every %.3f sec, Print Result every %.3f sec \n",(float)out[0]/1000,(float)out[1]/1000,(float)out[4]/1000);
	start=clock();
	srand(time(NULL));
	rtim[4]=-1; // junk value for MCS=0

	if(rtim[1]<=rscale[12]){
	    melt[1]=dtsave[0];
	    if(rtim[1]<melt[1]){ // fdm time > MC time
		melt[1]=rtim[1];
		fdm_loop=0;
		// in FDM loop, automatically ++
	    }else
		fdm_loop=(int)(rtim[1]/melt[1]);
	    melt[6]=LatHeatPerStep(dHfus[1],melt[1],rscale[1],melt[3],rscale[9],melt[4]);
	    dtsave[1]=rtim[1];
	}else{// if(rtim[1]>rscale[12])
	    factor[0]=factor[0]*rscale[12]/rtim[1];
	    factor[1]=factor[1]*rscale[12]/rtim[1];
	    rtim[1]=rscale[12];
	}

	printf("    1 MCS = %E sec,",rtim[1]);
	if(mode[MLAT]!=0){
		printf(" Latent heat for S-L reaction is considered\n");
		printf("    t_gg  = %E sec, t_solidi = %E sec",dtdMCS[0],dtdMCS[1]);
		printf(", t_FDM = %E sec ",melt[1]);
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
		printf("    t_gg  = %E sec, t_solidi = %E sec",dtdMCS[0],dtdMCS[1]);
	}
	if(am[0]<0)
		printf("    Circular shaped melt pool\n");
	else if(am[0]>0)
		printf("    Gaussian LASER irradiation\n");
	else
		printf("    Teardropshaped melt pool\n");
	if(out[3]<0)
		printf("    # Simulation will automatically stop when there is no liquid... #\n\n");
		
	dtsave[1]=rtim[1];

	if(am[0]==0){
	    mp[0]=4*atan(sqrt(2-sqrt(3))); // temporarily, mp[0] = theta
	    mp[1]=am[3]/2;
	    mp[2]=am[2]/(2*sin(mp[0])*pow(sin(mp[0]*0.5),am[5]));
	    mp[3]=am[4]/(sin(mp[0])*pow(sin(mp[0]*0.5),am[5]));
	}else if(am[0]>0){ // GAUSSIAN
	    mp[5]=am[9]*am[5]*GAUSS2D/(am[2]*am[3])*pow(rscale[9],2/3)/melt[4];
	    mp[1]=mp[5]*rtim[1];
	    am[2]=am[2]/rscale[1];
	    am[3]=am[3]/rscale[1];
	    am[4]=am[4]/rscale[1];
	    mp[2]=1/(am[2]*am[2]);
	    mp[3]=1/(am[3]*am[3]);
	    mp[4]=1/(am[4]*am[4]);
	}
	mp[0]=-am[8]/rscale[1]*rtim[1]-1;
	pos[0]=am[6];
	pos[1]=am[7];

/////////////////////////////////////////////////////////////   MC Loop start   /////
	while(1){//	srand(time(NULL));
	    if(am[9]!=0){ // TMP or ABS == 0 -> No melt pool calculation
	      if(am[0]<=0){ // TEARDROP mode or CIRCLE mode
		if(mp[0]<0){ // 1st MCS
		    mp[0]=0;
		    MeltPool(grid,fdmt,am,mp,dir,pbc,Tmax[0],dHfus[1],pos);
		}else{
		    mp[0]=am[8]/rscale[1]*rtim[1];
		    MeltPool_movement(pos,am[1],mp[0]);
		    MeltPool(grid,fdmt,am,mp,dir,pbc,Tmax[0],dHfus[1],pos);
		}
// Meltpool이 SAMPLE 밖으로 나갔을때 예외처리
		if(pos[0]<0 || pos[0]>dir[0] || pos[1]<0 || pos[1]>dir[1])
		    am[9]=0;

	      }else if(am[0]>0){ // GAUSSIAN mode
// 짧은 timestep 가져가면서 계속 지지기
		if(mp[0]<0) // 1st MCS
		    mp[0]=0;
		else{
		    mp[0]=am[8]/rscale[1]*rtim[1];
		    mp[1]=mp[5]*rtim[1];
		    MeltPool_movement(pos,am[1],mp[0]);
		}
		MeltPool(grid,fdmt,am,mp,dir,pbc,Tmax[0],dHfus[1],pos);
// Meltpool이 SAMPLE 밖으로 나갔을때 예외처리
		if(pos[0]<0 || pos[0]>dir[0] || pos[1]<0 || pos[1]>dir[1])
		    am[9]=0; // 
	      }
	    }

//	    if(jc[7]!=0) // PTCL introduction
//		Particle_Form(grid,dir,pbc,jc[7]);
	    if(mode[MPRO]>=0){
		heat_transfer(grid,fdmt,tbc,dir,pbc,rscale[1],melt,fdm_loop,k_per_cp,Tmax[0],mode[MPRO]); // Loop inside the function
		if(am[9]>=Tmax[0])
		    check_melting(grid,fdmt,Tmax[0],dir);
	    }
	    mcs++;
	    step=total;
//printf("MCS %d\n",mcs);
	    while(step!=0){	// 1 MCS ////////////////////////////////////////////////////////////////////////////////
		i=(int)(rand()/(RAND_MAX/dir[0]+1));
		j=(int)(rand()/(RAND_MAX/dir[1]+1));
		k=(int)(rand()/(RAND_MAX/dir[2]+1));

		MCtrial_SOL(i,j,k,grid,fdmt,dir,pbc,mode,melt,Tmax,dHfus,jc,tbc,Vsite,Asite,factor,Nsite,pcps,rscale[3],dfeq);
		step--;
	    }
//puts("MCS DONE");
	    rtim[0]+=rtim[1];

	    if(out[5]==0){
		dmc[2]--;
		if(dmc[2]==0){ // time OPT
		    dmc[2]=out[6];
		    rtim[1]=fixed_dt_set(mcs,&(mode[MPRO]),&fdm_loop,dtsave,dtdMCS,dHfus,Tinfo,factor,Vsite,Asite,out,Tmax,pcps,grid,fdmt,dir,pbc,rscale,melt);
		}
		if(rtim[1]<0){
		    printf(" @@@ Liquid fraction is 0, solidification is DONE @@@\n");
		    break;
		}
	    }

	    if(out[3]<=0){ // Time criteria is MCS
		dmc[0]--;
		dmc[1]--;
		if(dmc[0]==0){ // Data out
		    printf(" %dth MC Step was taken...\n",mcs);
		    datout(dir,grid,mcs,out,rtim,rscale[0]); // Write data file
		    dmc[0]=out[4];
		}
		if(dmc[1]==0){ // Create graphics file
		    fileout(dir,grid,mcs);
		    tproout(dir,fdmt[0],mcs);
		    dmc[1]=out[1];
		}
		if(mcs==out[0])
		    break;
	    }else{ // Time criteria is sec
		if(dmc[0]*out[4]<=rtim[0]*1000){ // Data out
		    printf(" %.5f sec is passed...\n",rtim[0]);
		    datout(dir,grid,mcs,out,rtim,rscale[0]); // Write data file
		    dmc[0]++;
		}
		if(dmc[1]*out[1]<=rtim[0]*1000){ // Create graphics file
		    fileout(dir,grid,mcs);
		    dmc[1]++;
		}
		if(rtim[0]*1000>=out[0])
		    break;
	    }
	} // End a MCS

	if(out[3]<=0){ // Time criteria is MCS
	    if(out[0]!=mcs)
		out[0]=mcs;

		if(dmc[0]!=out[4]) // Check unsaved MCSs
			datout(dir,grid,mcs,out,rtim,rscale[0]);
		if(dmc[1]!=out[1]){
			fileout(dir,grid,mcs);
			tproout(dir,fdmt[0],mcs);
		}
	}else{ // Time criteria is sec
	    if(out[0]!=rtim[0]*1000)
		out[0]=rtim[0]*1000;

		if(dmc[0]*out[4]!=rtim[0]*1000)
		    datout(dir,grid,mcs,out,rtim,rscale[0]);
		if(dmc[1]*out[1]!=rtim[0]*1000) // Create graphics file
		    fileout(dir,grid,mcs);
	}
	end=clock();

	for(i=2;i>=0;i--)
	    fdmt[i]=free3d(dir[0],dir[1],dir[2],fdmt[i]);
	free(fdmt);

	printf("\nExecution time: %d hour %d min %d sec\n",(int)((end-start)/CLOCKS_PER_SEC)/3600,(int)(((end-start)/CLOCKS_PER_SEC)%3600)/60,(int)(((end-start)/CLOCKS_PER_SEC)%3600)%60);
	printf("== Simulation is done ==\n");

	if(mode[MLAT]==0)
	    melt[3]=-1;
	else if(melt[3]<0)
	    melt[3]=0;


	return;
}

void heat_transfer(int**** grid,double**** tgrid,double tbc[][2],int dir[],int pbc[],double mcdx,double const melt[],int fdm_loop,double k_per_cp,double Tmelt,int mode){ //,double dHfus){
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
				for(k=fdmax[2]-inter;k>0;k-=inter){
				  if(grid[i][j][k][1]==POWORI)
				  tgrid[1][i][j][k]=tgrid[0][i][j][k]+melt[1]*k_per_cp*melt[9]*\
						   (revdx2*(tgrid[0][i+inter][j][k]+tgrid[0][i-inter][j][k]-2*tgrid[0][i][j][k])+\
						    revdx2*(tgrid[0][i][j+inter][k]+tgrid[0][i][j-inter][k]-2*tgrid[0][i][j][k])+\
						    revdx2*(tgrid[0][i][j][k+inter]+tgrid[0][i][j][k-inter]-2*tgrid[0][i][j][k]));
				  else
				  tgrid[1][i][j][k]=tgrid[0][i][j][k]+melt[1]*k_per_cp*\
						   (revdx2*(tgrid[0][i+inter][j][k]+tgrid[0][i-inter][j][k]-2*tgrid[0][i][j][k])+\
						    revdx2*(tgrid[0][i][j+inter][k]+tgrid[0][i][j-inter][k]-2*tgrid[0][i][j][k])+\
						    revdx2*(tgrid[0][i][j][k+inter]+tgrid[0][i][j][k-inter]-2*tgrid[0][i][j][k]));
				}
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


