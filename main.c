#include "SampleMaker.h"
#include "Single.h"
#include "MPI.h"
#include "MC.h"
#include "Database.h"
#include "Visual.h"
#include "Utility.h"
#include "Calc.h"

/*	Update log
 *	May   15, 2020	// MC Potts + Metropolis method			// version 0.1
 *	May   18, 2020	// Change simulation name: "SHGGSim"		// version 0.2
 *			// Include heatmap graph plotting module
 *			//		by utilizing gnuplot program
 *	May   20, 2020	// Add MC trial mode:				// version 0.3
 *			//		scan whole feature or Random pick
 *			// Calculate Volume Fraction of Grains
 *			// Execution time display
 *			// 3D mode + PBC on/off
 *			// Graph module for 3D sample included
 *	May   26, 2020	// Seperatre Modules
 *			// Improve grain make module			// version 0.34
 *	June   1, 2020	// 1) Improve grain make module:		// version 0.4
 *						multiply, cut, paste
 *			   2) Update search method
 *	June   2, 2020	// 1) Log file control
 *			   2) Grain size measurement by AGI method
 *			   3) Speed Optimization
 *	June   3, 2020	// Take real orientation infomation		// version 0.5
 *					&& delet mctrial_whole module
 *			// Change graphics scheme
 *	June  11, 2020	// Speed Optimization, Error fix(theta=0 case)	// version 0.52
 *	June  12, 2020	// Take GB Mobility Scheme (~T dependence)	// version 0.6
 *	June  18, 2020	// 1) Change GBE format
 *			// 2) Make Input-Log File (input.txt)
 *			// 3) Hetero-Material (Different Tm)
 *	June  29, 2020	// Input T as a function
 *	July   3, 2020	// 1) High interface energy at liq-sol interface
 *			// 2) Nucleation at interface
 *			// 3) Draw temperature profile
 *	July   4, 2020	// Take nucleation possibility scheme
 *	July  22, 2020	// T input with MCS range			// version 0.66
 *	July  30, 2020  // Rotate sample module included		// version 0.70
 *	August 1, 2020	// Add melt pool mode (AM simulation)
 *	August 4, 2020	// Add HAZ zone around melt pool
 *			// Speed optimization
 *			// Modify melt pool temperature
 *			// Edit plot scheme
 *	Sept.  3, 2020	// Change melt pool geometry
 *	Sept. 15, 2020	// Optimization of AGI method
 *			// Error fix
 *			// Real T profile applied for melt pool
 *	Oct.  22, 2020	// Melt pool movement				// version 0.77
 *	Nov.   4, 2020	// Code seperation				// version 0.80
 *			// Grain growth factor modification
 *	Nov.  25, 2020	// Applying real-scale time model
 *	Apr.   1, 2021	// Change schemes for M and G (probability)
 *	Nov.   9, 2021	// Solidification with epitaxial growth
 *	Nov.  19, 2021	// Real T profile as a function of real time
 *	Mar.   8, 2022	// Change DB file format
 *			// Add MCS-Time plot module
 *			// Error fixed (in proflie function)		// version 0.83
 *	Apr.   2, 2022	// Error fixed (for changing temperature)	// version 0.90
 *	May    5, 2022	// Change orientation description: Euler angle	// version 0.91
 *	Jun.   1, 2022	// Change scheme for real time conversion	// version 1.00
 *	Aug.  23, 2022	// Error fix (Grain growth at middle of str.)
 *	Sep.  29, 2022	// Error fix					// version 1.02
 *			// Add solidification model & FDM for heat eq.	// version 1.1
 *	Nov.  28, 2022	// Error fix (FDM)
 *			// Add melting model
 *	Feb.  20, 2022	// Metropolis function convergence (step ftn)	// version 1.2
 *			// Error fix (Melting model)
 *	May   26, 2023	// Error fix (Metropolis function convergence)	// version 1.21
 *			// Real time conversion for phase transf.
 *			// Combination of real time conversion(GG+TRF.)
 *	Oct.   5, 2023	// Nucleation model				// version 1.30
 *	Nov.   7, 2023	// Different dx value at FDM and MC
 *			// Speed optimization				// version 1.31
 *	Mar.  14, 2024	// Final nucleation model
 *			// Final realistic time assignment	
 *			// Mould interface condition (wetting angle)	// version 1.40
 *	Mar.  22, 2024	// Materials Information module updated
 *			// Flood-fill code for measuring grain size
 *			// Solid-liquid IF position finding module
 *			// Real time assignment for complex situation
 *			// Primary Dendrite Arm Spacing module
 *			// Avg. T calculation module
 *	May    5, 2024	// Fix compile error due to inline functions
 *	Aug.  27, 2024	// Change save format of data file (MCS*.dat)
 *			// Change read & write method of data file
 *	Sep.  12, 2024	// Support MPI computing (GG module)		// version 1.45
 *	Oct.   8, 2024	// MPI optimization				// version 1.5
 *	Nov.  22, 2024	// Solid-liquid IF energy anisotropy?
 *	Nov.  26, 2024	// FDM roll-back: voxel-grid 1:1 matching,
 *			//	         to avoid remaining-L artifacts
 *			// && latent heat retroactive application	// version 1.52
 *	Feb.   2, 2025	// INCAR script format is established		// version 1.60
 *			// SOL module optimization: remove tgrid[][][]
 *	Feb.  19, 2025	// AM module update
 *	Apr.   1, 2025	// Speed optimization scheme revised
 *	May    8, 2025	// Gaussian LASER temperature scheme updated
 *	May   20, 2025	// Support .vtk file format
 *	Jun.  19, 2025	// AM module update
 *	Jun.  27, 2025	// Module files (codes & headers) rearrangement
 *	Jul.   1, 2025	// MPI optimization
 *			// Supprot MPI computing (SOL module)
 *			// Supprot MPI computing (AM module)
 *	Jul.   9, 2025	// MPI debugging (SOL, AM module)		// version 1.68
 *	Jul.  22, 2025	// New solidification model			// version 1.7
 *	Aug.   1, 2025	// Directional term in solidification model	// version 1.71
 */
int ver=171;
int n_size;
int n_rank;

int**** init_setting(int**** grid,double jc[],int dir[],int out[],char ttem[],double melt[],int pbc[],double rscale[],double tbc[][2],double am[]); // initial setting
int read_incar(double jc[],int dir[],int out[],char ttem[],double melt[],int pbc[],double rscale[],double tbc[][2],char str[],double am[]); // input file reading
int set_bc(FILE *incar,int pbc[],double tbc[][2],double mode);
int set_wetting(FILE *incar,double* angle,char word[],double mode);
int is_valid_digits(const char *str);
void example_incar();
void input_log(int dir[],int pbc[],int out[],char ttem[],double melt[],double jc[],double rscale[],double tbc[][2],double am[]); // Write input information

int main(int argc,char** argv){
	MPI_Init(&argc,&argv);

	MPI_Comm_rank(MPI_COMM_WORLD,&n_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&n_size);

	int pbc[3]={0},out[8]={0},dir[3]={0},i,j=1;
// dir 0 1 2 for MC (x y z length) // 3 4 5 for FDM (x y z number of grids) in previous version
// out 0 = total -> real total, 1 = step to print map 2 = AGI lines 3 = criteria for ending simul. (sec = 1 vs MCS = 0), 4 = step to print data
//     5 = Timestep Acceleartion scheme (0, 1, 2), 6 = step to do active dt setting, 7 = Liquid fraction criteria? (1 == Yes // 0 == No)
	char func,ttem[100]={'\0'};
	double sang[11]={0},tbc[6][2]={0},jc[9]={0},melt[8]={0},rscale[15]={0},am[10]={0};
	// am 0 = Melt pool mode (-1 = Circle, 0 = TEARDROP,  1 = Gaussian) 
	// 1 = Melt pool direction (negative -> X, positive -> Y)
	// Start position (x = 6, y = 7), 8 = scan speed
	// when am[0]==-1: 2 = Radius                      4 = Depth                                9 = Boundary Temperature
	// when am[0]==0:  2 = Width  3 = Cap+Tail length  4 = Depth  5 = Geometry (Teardrop shape) 9 = Boundary Temperature
	// when am[0]==1:  2 = SX     3 = SY               4 = SZ    5 = POW  9 = ABS    (Gaussian)
	// jc 0=Qcutoff(degree) // 1 = GBEmax // 2 = IFE_aniso factor // 3=C(low angle boundary,degree) // 4=solid/liquid interface energy(J/m2) // 5=solid-vapor IFE (J/m2) // 6=liquid-vapor IFE 
	// 			// 7 = particle volume fraction in liquid, 8 = wetting angle // 9 =temp (empty) for dendr. dir.
	// ttem: (time range 1) (t expression) (time range 2) (t expression) ...
	// melt 0 = Max. T  1 = dt for FDM in SOL module // accel. factor in GG module
	//      2 = heat transfer mode (0 = excluded // T = included (T=initial temperature, T>0) // <0 = load T profile)
	// (for FDM mode)     3 = thermal conductivity (k)  4 = Cp  5 = Tboil, 6 = latent heat retroactive portion; dt_fdm / (dx2/alpha) * dHfus
	//		      7 = Acceleration factor for FDM (integer)
	// (for FTC approx)   3 = T* at IF 4 = -(G K/voxel)
	// tbc [xyz][] 0 = x, lower  1 = x, upper  2 = y, lower  3 = y, upper  4 = z, lower  5 = z, upper
	// tbc [][0]  0 is for adiabatic, negative number for convection (= -T_surr), positive number for heat sink (=T_surr)
	// tbc [][1]  (for convection) h (for adiabatic or heat sink) wetting angle (0 < theta < 180, degree)
	// rscale 0 = does it do real scale conversion? (1 to yes, 2 to AM module)
	//	L = Lamda L_MC
	// 	  1 = Lamda
	//	P = exp (- dx2 gamma/kT) exp (- Qmob/R (1/T-1/Tm))
	//	  2 = Tmelt  3 = Qmob
	//	t= Lamda^2 KMC/Ke tMCS
	// 	  4 = Qsol_diffu  5 = K0_e  6 = Qliq_diffu  7 = K0_mc
	// 	Solidification & melting
	// 	  8 = dH_fus  9 = Vm 
	// 	  10 = KMC_solidi/Ke_solidi; time conversion factor for solidification / melting
	// 	  11 = Nucleation mode (1 == On)
	// 	  12 = Maximum timestep length for acceleration
	// 	  13 = Maximum timestep scaling factor for acceleration
	// 	  14 = Critical supercooling to consider solidification
	int**** grid=NULL;

    if(n_rank==0){
	printf("=============================================================\n");
	printf("!  Stochastic HIgh-speed Microstructure Simulation          !\n");
	printf("!                      based on the Monte Carlo Potts Model !\n");
	printf("!                                         SHIMS  ver. %5.2f !\n",(float)ver/100);
	printf("!                                                           !\n");
	printf("!      Coded by  S.-H. Oh                                   !\n");
	printf("!                CMSE, POSTECH                              !\n");
	printf("!                shoh97@postech.ac.kr                       !\n");
	printf("=============================================================\n\n");
    }
	while(1){
	    if(n_rank==0){
		printf("\nSelect mode:\n");
		printf("  1 / m = Sample Maker\n");
		printf("  2 / s = Monte Carlo Potts simulation\n");
		printf("  4 / t = Auxiliary Tools: Graphics & Measurement\n");
		printf("  5 / d = Database Checker\n");
		printf("  0 / q = Quit\n  >> ");
		scanf("%c",&func);

		     if(func=='5' || func=='d')
			Database(melt,jc,rscale);
		else if(func=='4' || func=='t')  // Draw graph with GNUPLOT
			graphout(1,out,dir,0,rscale,0);
		else if(func=='1' || func=='m')  // Make Grain Sample
			gmake(dir,grid);
		else if(func=='0' || func=='q'){
			if(n_size==1)
				return 0;
		}else if(func=='2' || func=='s'){
			grid=init_setting(grid,jc,dir,out,ttem,melt,pbc,rscale,tbc,am);
			if(grid==NULL){
				printf(" Simulation aborted...\n");
				return 0;
			}
		}else
			printf("Wrong input...\n");
	    }
	    if(n_size>1){
		MPI_Barrier(MPI_COMM_WORLD);

		MPI_Bcast(dir,3,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(pbc,3,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(out,7,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(jc,9,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(rscale,15,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(melt,8,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(am,10,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(ttem,100,MPI_CHAR,0,MPI_COMM_WORLD);
		MPI_Bcast(&func,1,MPI_CHAR,0,MPI_COMM_WORLD);
		for(i=0;i<6;i++){
		    MPI_Bcast(&(tbc[i][0]),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		    MPI_Bcast(&(tbc[i][1]),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		}
		if(func=='0' || func=='q')
			return n_rank;
	    }
	    if((func=='2' || func=='s')){ // to avoid MPI issue
	      if(n_size==1){
		if(melt[2]==0 && melt[3]>=0)
		    GGmodule(grid,jc,dir,out,ttem,melt,pbc,sang,rscale,tbc);	//MC simulation without FDM
		else{
		    if(rscale[0]==1)
			SOLmodule(grid,jc,dir,out,melt,pbc,rscale,tbc);	//MC simulation with FDM
		    else if(rscale[0]==2)
			AMmodule(grid,jc,dir,out,melt,pbc,rscale,tbc,am);	//AM simulation
		}
	      }else{
		if(melt[2]==0 && melt[3]>=0)
		    GGmodule_MPI(grid,jc,dir,out,ttem,melt,pbc,sang,rscale,tbc);	//MC simulation without FDM
		else{
		    if(rscale[0]==1)
			SOLmodule_MPI(grid,jc,dir,out,melt,pbc,rscale,tbc);	//MC simulation with FDM
		    else if(rscale[0]==2)
			AMmodule_MPI(grid,jc,dir,out,melt,pbc,rscale,tbc,am);	//AM simulation
		}	      
	      }
		if(n_rank==0){
		    grid=free4d(dir[0],dir[1],dir[2],grid);
		    if(out[3]<=0)
			graphout(0,out,dir,(int)rscale[2],rscale,melt[3]);  // Draw graph with GNUPLOT
		}
	    }
	if(n_rank==0)
		while ((getchar()) != '\n');

		// Initialization
		for(i=0;i<3;i++)
			pbc[i]=0;
		for(i=0;i<4;i++){
			dir[i]=0;
			jc[i]=0;
			out[i]=0;
		}
		out[4]=0;
		for(i=0;i<100;i++){
			if(ttem[i]=='\0')
				break;
			ttem[i]='\0';
		}
		for(i=0;i<11;i++){
			if(sang[i]==0)
				break;
			sang[i]=0;
		}
		for(i=0;i<11;i++){
			rscale[i]=0;
		}
		melt[0]=0;
		melt[1]=0;
		rscale[0]=0;
		if(j>100){
			printf(" ### INF. LOOP... SYSOUT ###\n");
			break;
		}
		j++;
	}
	MPI_Finalize();

	return 0;
}

int**** init_setting(int**** grid,double jc[],int dir[],int out[],char ttem[],double melt[],int pbc[],double rscale[],double tbc[][2],double am[]){ // initial setting
	char temp[100],ri[100];
	int i;

	printf("\n");
	if(read_incar(jc,dir,out,ttem,melt,pbc,rscale,tbc,ri,am)!=0)
		return NULL;

	puts(" # Initial setting input is completed! #");
	grid=read_info(grid,dir,ri); // read sample information
	if(grid==NULL){
		printf("Error: Couldn't load sample file...\n\n");
		return NULL;
	}

	for(i=0;i<3;i++){
	    if((double)(dir[i]-1)/melt[7]!=floor((dir[i]-1)/melt[7])){
		printf(" INPUT ERROR for 'FACL' option... (Total sample size is not divisible...)\n");
		grid=free4d(dir[0],dir[1],dir[2],grid);
		return NULL;
	    }
	}

	fileout(dir,grid,0); // make MCS0.dat file
	input_log(dir,pbc,out,ttem,melt,jc,rscale,tbc,am);
	return grid;
}

int read_incar(double jc[],int dir[],int out[],char ttem[],double melt[],int pbc[],double rscale[],double tbc[][2],char str[],double am[]){ // input file reading
	int i=0,j,k;
	char temp[100],temp2[500],option[100];
	FILE* incar=NULL;

	incar=fopen(INCAR,"r");  // Input log file
	if(incar==NULL){
	    printf(" ### ERROR: CANNOT FIND INPUT FILE '%s'... ###\n",INCAR);
	    printf("     Example input file 'example.in' is created...\n\n");
	    example_incar();
	    return 1;
	}

	if(find_keyword(incar,"STR",option)!=0){ // ==0 when there is no error
		printf(" INPUT FILE ERROR in reading 'STR' option... \n");
		fclose(incar);
		return 1;
	}else
	    tok(option,str);

	if(find_keyword(incar,"REAL",option)!=0){ // ==0 when there is no error
		printf(" INPUT FILE ERROR in reading 'REAL' option... \n");
		fclose(incar);
		return 1;
	}
	if(option[0]=='Y' || option[0]=='y')
		rscale[0]=1;  // full-scale MC with time conversion
	else
		rscale[0]=0;  // No realistic time assignment
	
	if(find_keyword(incar,"MDB",option)!=0){ // ==0 when there is no error
		printf(" INPUT FILE ERROR in reading 'MDB' option... \n");
		fclose(incar);
		return 1;
	}

	if(strcmp(option,"NONE")==0){ // input materials parameter manually
	    printf("Manual input for simulation parameters\n");
	    printf("Maximum temperature for grain growth(K) = ");
	    scanf("%lf",&melt[0]);
	    printf("Melting point of Material (liquidus temperature) (K) = ");
	    scanf("%lf",&rscale[2]);
	    printf("Max Grain Boundary Energy for Misorientation?(J/m2, def=1.0)\n >> ");
	    scanf("%lf",&jc[1]);
	    printf("\nGB Energy Anisotropy Level(Degree) = "); // GBE information
	    scanf("%lf",&jc[3]);
	    printf("Grain boundary mobility energy barrier? (J/mol)\n Q = ");
	    scanf("%lf",&rscale[3]);
	    printf("Solid-Liquid interface energy? (J/m2)\n γ(sl) = ");
	    scanf("%lf",&jc[4]);
	    printf("Solid-Vapor interface energy? (J/m2)?\nγ(sv) = ");
	    scanf("%lf",&jc[5]);
	    printf("Liquid-Vapor interface energy? (J/m2)?\nγ(lv) = ");
	    scanf("%lf",&jc[6]);
	    printf("Enthalpy of fusion? (J/g-atom)\n dH_fus = ");
	    scanf("%lf",&rscale[8]);
	    printf("Molar volume? (m3/g-atom)\n Vm = ");
	    scanf("%lf",&rscale[9]);

	    if(rscale[0]!=0){
		printf(" 2) L^2 - L0^2 = K0_e * EXP( - Qe / R T ) * t\n");
		printf("                 K0_e (m^2/s) = ");
		scanf("%lf",&rscale[5]);
		printf("\n Lmc^2 - L0^2 = Kmc_0 * EXP( - Qmc_0 / R T ) * MCS\n");
		printf("                  Kmc_0 (site^2/s)= ");
		scanf("%lf",&rscale[7]);
//For solidification model, manual input mode is not supported...
	    }
	}else{
	    if(find_keyword(incar,"SYS",temp)!=0){ // ==0 when there is no error
		printf(" INPUT FILE ERROR in reading 'SYS' option... \n");
		fclose(incar);
		return 1;
	    }

	    if(Load_DB(option,temp,melt,jc,rscale)!=0){
		fclose(incar);
		return 1;
	    }
	}

	if(find_keyword(incar,"DX",option)!=0){ // ==0 when there is no error
	    printf(" INPUT FILE ERROR in reading 'DX' option...\n");
	    fclose(incar);
	    return 1;
	}
	tok(option,temp);
	if(is_valid_digits(temp))
	    rscale[1]=atof(temp);
	else{
	    printf(" INPUT ERROR for 'DX' option...\n");
	    fclose(incar);
	    return 1;
	}

	if(rscale[0]!=0){
	    if(find_keyword(incar,"MODE",option)!=0){ // ==0 when there is no error
		printf(" INPUT FILE ERROR in reading 'MODE' option... \n");
		fclose(incar);
		return 1;
	    }
	    if(option[0]=='C' || option[0]=='c') // SOL module
		rscale[0]=1;
	    else if(option[0]=='A' || option[0]=='a') // AM module
		rscale[0]=2;
	    else{
		printf(" INPUT ERROR for 'MODE' option... \n");
		fclose(incar);
		return 1;
	    }
	}

	if(rscale[0]==0){//melt[1]!=0){ // No real time conversion case
		printf("### NOTE: heat transfer model requires real time conversion...\n          Automatically disabled...\n ");
		melt[2]=0;
	}else{
		if(find_keyword(incar,"FDM",option)!=0){ // ==0 when there is no error
			printf(" INPUT FILE ERROR in reading 'FDM' option... \n");
			fclose(incar);
			return 1;
		}
		if(option[0]=='Y' || option[0]=='y')
			melt[2]=1;
		else
			melt[2]=0;
		if(melt[2]==0&&rscale[0]==2){
			printf(" ### ERROR: AM module requires FDM option ON... ###\n");
			fclose(incar);
			return 1;
		}

		if(find_keyword(incar,"FACL",option)!=0) // ==0 when there is no error
		    melt[7]=1;
		else{
		    melt[7]=atof(option);
		    if(melt[7]!=floor(melt[7])){
			printf(" INPUT ERROR for 'FACL' option... \n");
			fclose(incar);
			return 1;
		    }
		}
	}
	if(find_keyword(incar,"PBC",option)!=0){ // ==0 when there is no error
		printf(" INPUT FILE ERROR in reading 'PBC' option... \n");
		fclose(incar);
		return 1;
	}
	for(i=0;i<3;i++){
		if(option[i]=='1')
		    pbc[i]=1;
		else if(option[i]=='0')
		    pbc[i]=0;
	        else{
		    printf("### ERROR: Invalid PBC condition..\n\n");
		    fclose(incar);
		    return 1;
		}
	}

// NUCLEATION option ON or OFF
	if(find_keyword(incar,"NUCL",option)!=0) // ==0 when there is no error
		rscale[11]=1;
	else{
		if(option[0]=='Y' || option[0]=='y')
		    rscale[11]=1;
		else
		    rscale[11]=0;
	}

	if(melt[2]!=0){ // solving heat transfer equation
	    if(strcmp(option,"111")==0){
		printf("### ERROR: at least ONE SURFACE must exist for solving heat transfer problem...!\n\n");
		fclose(incar);
		return 1;
	    }
	    melt[1]=FDMSAFE*0.5*rscale[1]*rscale[1]*melt[4]/(rscale[9]*melt[3]); // dt for FDM
//		printf(" # NOTE: FDM grid point locates at the center of each MC voxels...\n");
//		printf("\n # dt for FDM is %E sec (safe factor %.3f)\n",melt[1],FDMSAFE);
	    if(find_keyword(incar,"T0",option)!=0){ // ==0 when there is no error
		printf(" INPUT FILE ERROR in reading 'T0' option... \n");
		fclose(incar);
		return 1;
	    }
	    tok(option,temp);
	    if(is_valid_digits(temp))
		melt[2]=atof(temp);
	    else
		melt[2]=-1; // loading will be done in SOL module
	    if(set_bc(incar,pbc,tbc,rscale[11])!=0){
		fclose(incar);
		return 1;
	    }
	}else{ // no heat transfer
	    if(find_keyword(incar,"FROZ",option)!=0){ // ==0 when there is no error
		printf(" INPUT FILE ERROR in reading 'FROZ' option... \n");
		fclose(incar);
		return 1;
	    }
	    if(option[0]=='Y' || option[0]=='y'){
		if(find_keyword(incar,"TINT",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in reading 'TINT' option... \n");
		    fclose(incar);
		    return 1;
		}
		tok(option,temp);
		melt[3]=-1*atof(temp);
		if(find_keyword(incar,"TGRA",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in reading 'TGRA' option... \n");
		    fclose(incar);
		    return 1;
		}
		tok(option,temp);
		melt[4]=atof(temp);
	    }else{ // NO FDM, NO FROZEN T APPROX. -> ISOTHERMAL ASSUMPTION
		i=0;
		while(++i){ // RANGE
		    sprintf(temp,"t(%d)",i);
		    if(find_keyword(incar,temp,option)!=0){ // ==0 when there is no error
			if(i==1){
			    printf(" INPUT FILE ERROR in reading '%s' option... \n",temp);
			    fclose(incar);
			    return 1;
			}else
			    break;
		    }
		    tok(option,temp2);
		    if(is_valid_digits(temp2)){
			if(i==1)
			    strcpy(ttem,temp2);
			else
			    sprintf(ttem,"%s %s",ttem,temp2);
		    }else{
			printf(" INPUT ERROR for '%s' option... \n",temp);
			fclose(incar);
			return 1;
		    }
		    // EQUATION
		    sprintf(temp,"T(%d)",i);
		    if(find_keyword(incar,temp,option)!=0){ // ==0 when there is no error
			printf(" INPUT FILE ERROR in reading '%s' option... \n",temp);
			fclose(incar);
			return 1;
		    }
		    tok(option,temp2);
		    Etodec(temp2);
		    sprintf(ttem,"%s %s",ttem,temp2);
		}

		if(find_keyword(incar,"tUNIT",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in reading 'tUNIT' option... \n");
		    fclose(incar);
		    return 1;
		}
		strcpy(temp,ttem);
		if(option[0]=='M' || option[0]=='m') // in MCS
		    sprintf(ttem,"0 %s",temp);
		else if(option[0]=='S' || option[0]=='s')
		    sprintf(ttem,"1 %s",temp);
		else{ // unit is neither MCS nor sec
		    printf(" INPUT ERROR for 'tUNIT' option unit...\n");
		    fclose(incar);
		    return 1;
		}
	    }
	}

// WETTING ANGLE
	if(tbc[0][0]>=0){
	    if(set_wetting(incar,&tbc[0][1],"WXL",rscale[11])!=0){
		fclose(incar);
		return 1;
	    }
	}
	if(tbc[1][0]>=0){
	    if(set_wetting(incar,&tbc[1][1],"WXU",rscale[11])!=0){
		fclose(incar);
		return 1;
	    }
	}
	if(tbc[2][0]>=0){
	    if(set_wetting(incar,&tbc[2][1],"WYL",rscale[11])!=0){
		fclose(incar);
		return 1;
	    }
	}
	if(tbc[3][0]>=0){
	    if(set_wetting(incar,&tbc[3][1],"WYU",rscale[11])!=0){
		fclose(incar);
		return 1;
	    }
	}
	if(tbc[4][0]>=0){
	    if(set_wetting(incar,&tbc[4][1],"WZL",rscale[11])!=0){
		fclose(incar);
		return 1;
	    }
	}
	if(tbc[5][0]>=0){
	    if(set_wetting(incar,&tbc[5][1],"WZU",rscale[11])!=0){
		fclose(incar);
		return 1;
	    }
	}
	if(rscale[0]==0)
	    out[3]=0;
	else{
	    if(find_keyword(incar,"CRI",option)!=0){ // ==0 when there is no error
		printf(" INPUT FILE ERROR in reading 'CRI' option... \n");
		fclose(incar);
		return 1;
	    }
	    if(option[0]=='M' || option[0]=='m')
		out[3]=0;
	    else if(option[0]=='S' || option[0]=='s')
		out[3]=1;
	    else if(option[0]=='L' || option[0]=='l')
		out[3]=-1;
	    else{
		printf(" INPUT ERROR for 'CRI' option... \n");
		fclose(incar);
		return 1;
	    }
	}

	if(find_keyword(incar,"NSTEP",option)!=0){ // ==0 when there is no error
	    printf(" INPUT FILE ERROR in reading 'NSTEP' option... \n");
	    fclose(incar);
	    return 1;
	}
	tok(option,temp);
	if(is_valid_digits(temp)){
	    if(out[3]<=0)
		out[0]=atoi(temp);
	    else // Use 1E-3 sec unit in GG module for computational efficiency... Maybe corrected later?
		out[0]=(int)(1000*atof(temp));
	}else{
	    printf(" INPUT ERROR for 'NSTEP' option... \n");
	    fclose(incar);
	    return 1;
	}
	if(find_keyword(incar,"NPRNT",option)!=0){ // ==0 when there is no error
	    printf(" INPUT FILE ERROR in reading 'NPRNT' option... \n");
	    fclose(incar);
	    return 1;
	}
	tok(option,temp);
	if(is_valid_digits(temp)){
	    if(out[3]<=0)
		out[1]=atoi(temp);
	    else // Use 1E-3 sec unit in GG module for computational efficiency... Maybe corrected later?
		out[1]=(int)(1000*atof(temp));
	}else{
	    printf(" INPUT ERROR for 'NPRNT' option... \n");
	    fclose(incar);
	    return 1;
	}
	
	if(find_keyword(incar,"NDATA",option)!=0){ // ==0 when there is no error
	    printf(" INPUT FILE ERROR in reading 'NDATA' option... \n");
	    fclose(incar);
	    return 1;
	}
	tok(option,temp);
	if(is_valid_digits(temp)){
	    if(out[3]<=0)
		out[4]=atoi(temp);
	    else // Use 1E-3 sec unit in GG module for computational efficiency... Maybe corrected later?
		out[4]=(int)(1000*atof(temp));
	}else{
	    printf(" INPUT ERROR for 'NDATA' option... \n");
	    fclose(incar);
	    return 1;
	}
// Additional information for each module
	if(melt[2]==0&&melt[3]>=0){ // GG module
	    if(find_keyword(incar,"ACCL",option)!=0) // ==0 when there is no error
		melt[1]=1.0;
	    else{
		tok(option,temp);
		melt[1]=atof(temp);
	    }
	}else{ // SOL module & AM module
	    if(find_keyword(incar,"TIME",option)!=0) // ==0 when there is no error
		out[5]=2;
	    else{
		tok(option,temp);
		out[5]=atoi(temp);
	    }
    	    if(out[5]>3 || out[5]<=0){
		printf(" INPUT ERROR for 'TIME' option...\n");
		fclose(incar);
		return 1;
	    }else{
		if(find_keyword(incar,"OPT",option)!=0) // ==0 when there is no error
		    out[6]=10;
		else{
		    tok(option,temp);
		    out[6]=atoi(temp);
		}
	    }

	    if(find_keyword(incar,"tMAX",option)!=0) // ==0 when there is no error
		rscale[12]=0;
	    else{
		tok(option,temp);
		rscale[12]=atof(temp);
	    }
    	    if(rscale[12]<0){
		printf(" INPUT ERROR for 'tMAX' option...\n");
		fclose(incar);
		return 1;
	    }else if(rscale[12]==0)
		rscale[12]=86400; //MAXIMUM TIME LIMIT <24 hr

	    if(find_keyword(incar,"tSCA",option)!=0) // ==0 when there is no error
		rscale[13]=0;
	    else{
		tok(option,temp);
		rscale[13]=atof(temp);
	    }
    	    if(rscale[13]<0||(rscale[13]!=0&&rscale[13]<1)){
		printf(" INPUT ERROR for 'tSCA' option...\n");
		fclose(incar);
		return 1;
	    }else if(rscale[13]==0)
		rscale[13]=1E10; //MAXIMUM TIME SCALING LIMIT <10^10 hr

	    if(find_keyword(incar,"TCS",option)!=0) // ==0 when there is no error
		rscale[14]=0;
	    else
		rscale[14]=atof(option);

	    if(find_keyword(incar,"LAT",option)!=0) // ==0 when there is no error
		melt[6]=1;
	    else{
		if(option[0]=='Y' || option[0]=='y')
		    melt[6]=1;
		else
		    melt[6]=0;
	    }
	}
	if(rscale[0]==2){ // AM module
	    if(find_keyword(incar,"DIR",option)!=0){ // ==0 when there is no error
		printf(" INPUT FILE ERROR in READING 'DIR' option...\n");
		fclose(incar);
		return 1;
	    }
	    if(option[0]=='X' || option[0]=='x')
		am[1]=-1;
	    else if(option[0]=='Y' || option[0]=='y')
		am[1]=1;
	    else{
		printf(" INPUT ERROR for 'DIR' option...\n");
		fclose(incar);
		return 1;
	    }
	    if(find_keyword(incar,"POS",option)!=0){ // ==0 when there is no error
		printf(" INPUT FILE ERROR in READING 'POS' option...\n");
		fclose(incar);
		return 1;
	    }
	    tok(option,temp);
	    if(is_valid_digits(temp))
		am[6]=atof(temp);
	    else{
		printf(" INPUT ERROR for 'POS' x option...\n");
		fclose(incar);
		return 1;
	    }
	    tok(option,temp);
	    if(is_valid_digits(temp))
		am[7]=atof(temp);
	    else{
		printf(" INPUT ERROR for 'POS' y option...\n");
		fclose(incar);
		return 1;
	    }

	    if(find_keyword(incar,"VEL",option)!=0){ // ==0 when there is no error
		printf(" INPUT FILE ERROR in READING 'VEL' option...\n");
		fclose(incar);
		return 1;
	    }
	    tok(option,temp);
	    if(is_valid_digits(temp))
		am[8]=atof(temp);
	    else{
		printf(" INPUT ERROR for 'VEL' option...\n");
		fclose(incar);
		return 1;
	    }

	    if(find_keyword(incar,"MMOD",option)!=0){ // ==0 when there is no error
		printf(" INPUT ERROR for 'MMOD' option...\n");
		fclose(incar);
		return 1;
	    }
	    if(option[0]=='T' || option[0]=='t')
		am[0]=0;
	    else if(option[0]=='G' || option[0]=='g')
		am[0]=1;
	    else if(option[0]=='C' || option[0]=='c')
		am[0]=-1;
	    else{
		printf(" INPUT ERROR for 'MMOD' option...\n");
		fclose(incar);
		return 1;
	    }
	    if(am[0]==0){ // Teardrop shaped meltpool
		if(find_keyword(incar,"WID",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in READING 'WID' option...\n");
		    fclose(incar);
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    am[2]=atof(temp);
		else{
		    printf(" INPUT ERROR for 'WID' option...\n");
		    fclose(incar);
		    return 1;
		}
		if(find_keyword(incar,"CTL",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in READING 'CTL' option...\n");
		    fclose(incar);
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    am[3]=atof(temp);
		else{
		    printf(" INPUT ERROR for 'CTL' option...\n");
		    fclose(incar);
		    return 1;
		}
		if(find_keyword(incar,"DEP",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in READING 'DEP' option...\n");
		    fclose(incar);
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    am[4]=atof(temp);
		else{
		    printf(" INPUT ERROR for 'DEP' option...\n");
		    fclose(incar);
		    return 1;
		}
		if(find_keyword(incar,"GEO",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in READING 'GEO' option...\n");
		    fclose(incar);
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    am[5]=atof(temp);
		else{
		    printf(" INPUT ERROR for 'GEO' option...\n");
		    fclose(incar);
		    return 1;
		}

		if(find_keyword(incar,"TMP",option)!=0) // ==0 when there is no error
		    am[9]=rscale[2]; // default value == Tm
		else{
		    tok(option,temp);
		    if(is_valid_digits(temp))
			am[9]=atof(temp);
		    else{
			printf(" INPUT ERROR for 'TMP' option...\n");
			fclose(incar);
			return 1;
		    }
		}
	    }else if(am[0]<0){ // Circular shaped Meltpool
		if(find_keyword(incar,"RAD",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in READING 'RAD' option...\n");
		    fclose(incar);
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    am[2]=atof(temp);
		else{
		    printf(" INPUT ERROR for 'RAD' option...\n");
		    fclose(incar);
		    return 1;
		}

		if(find_keyword(incar,"DEP",option)!=0) // ==0 when there is no error
		    am[4]=am[2]; // Default == radius value
		else{
		    tok(option,temp);
		    if(is_valid_digits(temp))
			am[4]=atof(temp);
		    else{
			printf(" INPUT ERROR for 'DEP' option...\n");
			fclose(incar);
			return 1;
		    }
		}
		if(find_keyword(incar,"TMP",option)!=0) // ==0 when there is no error
		    am[9]=rscale[2]; // default value == Tm
		else{
		    tok(option,temp);
		    if(is_valid_digits(temp))
			am[9]=atof(temp);
		    else{
			printf(" INPUT ERROR for 'TMP' option...\n");
			fclose(incar);
			return 1;
		    }
		}
	    }else if(am[0]>0){ // Gaussian melting
		if(find_keyword(incar,"POW",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in READING 'POW' option...\n");
		    fclose(incar);
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    am[5]=atof(temp);
		else{
		    printf(" INPUT ERROR for 'POW' option...\n");
		    fclose(incar);
		    return 1;
		}
		if(find_keyword(incar,"SX",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in READING 'SX' option...\n");
		    fclose(incar);
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    am[2]=atof(temp);
		else{
		    printf(" INPUT ERROR for 'SX' option...\n");
		    fclose(incar);
		    return 1;
		}
		if(find_keyword(incar,"SY",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in READING 'SY' option...\n");
		    fclose(incar);
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    am[3]=atof(temp);
		else{
		    printf(" INPUT ERROR for 'SY' option...\n");
		    fclose(incar);
		    return 1;
		}
/*		if(find_keyword(incar,"SZ",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in READING 'SZ' option...\n");
		    fclose(incar);
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    am[4]=atof(temp);
		else{
		    printf(" INPUT ERROR for 'SZ' option...\n");
		    fclose(incar);
		    return 1;
		}*/
		if(find_keyword(incar,"ABS",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in READING 'ABS' option...\n");
		    fclose(incar);
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    am[9]=atof(temp);
		else{
		    printf(" INPUT ERROR for 'ABS' option...\n");
		    fclose(incar);
		    return 1;
		}
		fclose(incar);
		return 0;
	    }
	    if(find_keyword(incar,"PTCL",option)!=0) // ==0 when there is no error
		jc[7]=0;
	    else{
		tok(option,temp);
		if(is_valid_digits(temp))
			jc[7]=atof(temp);
		else{
			printf(" INPUT ERROR for 'PTCL' option...\n");
			fclose(incar);
			return 1;
		}
	    }
	    if(find_keyword(incar,"PANG",option)!=0) // ==0 when there is no error
		    jc[8]=0;
	    else{
		    tok(option,temp);
		    if(is_valid_digits(temp))
			jc[8]=atof(temp);
		    else{
			printf(" INPUT ERROR for 'PANG' option...\n");
			fclose(incar);
			return 1;
		    }
	    }
	}
	fclose(incar);
	return 0;
}

int set_bc(FILE *incar,int pbc[],double tbc[][2],double mode){
	int i,j;
	char option[100],temp[100];

	if(pbc[0]==0){ // X surface
	    if(find_keyword(incar,"BXL",option)!=0){ // ==0 when there is no error
		printf(" INPUT FILE ERROR in reading 'BXL' option... \n");
		return 1;
	    }
	    if(option[0]=='A' || option[0]=='a')
		i=0;
	    else if(option[0]=='C' || option[0]=='c')
		i=1;
	    else if(option[0]=='F' || option[0]=='f')
		i=2;
	    else{
		printf(" INPUT ERROR for 'BXL' option... \n");
		return 1;
	    }
	    // Parameters for boundary conditoin
	    if(i==0){ // Adiabatic
		if(set_wetting(incar,&tbc[0][1],"WXL",mode)!=0) // Wetting angle
		    return 1;
		tbc[0][0]=0;
	    }else if(i==1){ // Convection
		if(find_keyword(incar,"TXL",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in reading 'TXL' option... \n");
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    tbc[0][0]=-1*atof(temp);
		else{
		    printf(" INPUT ERROR for 'TXL' option... \n");
		    return 1;
		}
		if(find_keyword(incar,"hXL",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in reading 'hXL' option... \n");
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    tbc[0][1]=-1*atof(temp);
		else{
		    printf(" INPUT ERROR for 'hXL' option... \n");
		    return 1;
		}
	    }else if(i==2){ // Fixed T
		if(set_wetting(incar,&tbc[0][1],"WXL",mode)!=0) // Wetting angle
		    return 1;
		if(find_keyword(incar,"TXL",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in reading 'TXL' option... \n");
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    tbc[0][0]=atof(temp);
		else{
		    printf(" INPUT ERROR for 'TXL' option... \n");
		    return 1;
		}
	    }
	    if(find_keyword(incar,"BXU",option)!=0){ // ==0 when there is no error
		printf(" INPUT FILE ERROR in reading 'BXU' option... \n");
		return 1;
	    }
	    if(option[0]=='A' || option[0]=='a')
		i=0;
	    else if(option[0]=='C' || option[0]=='c')
		i=1;
	    else if(option[0]=='F' || option[0]=='f')
		i=2;
	    else{
		printf(" INPUT ERROR for 'BXU' option... \n");
		return 1;
	    }
	    // Parameters for boundary conditoin
	    if(i==0){ // Adiabatic
		if(set_wetting(incar,&tbc[1][1],"WXU",mode)!=0) // Wetting angle
		    return 1;
		tbc[1][0]=0;
	    }else if(i==1){ // Convection
		if(find_keyword(incar,"TXU",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in reading 'TXU' option... \n");
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    tbc[1][0]=-1*atof(temp);
		else{
		    printf(" INPUT ERROR for 'TXU' option... \n");
		    return 1;
		}
		if(find_keyword(incar,"hXU",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in reading 'hXU' option... \n");
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    tbc[1][1]=-1*atof(temp);
		else{
		    printf(" INPUT ERROR for 'hXU' option... \n");
		    return 1;
		}
	    }else if(i==2){ // Fixed T
		if(set_wetting(incar,&tbc[1][1],"WXU",mode)!=0) // Wetting angle
		    return 1;
		if(find_keyword(incar,"TXU",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in reading 'TXU' option... \n");
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    tbc[1][0]=atof(temp);
		else{
		    printf(" INPUT ERROR for 'TXU' option... \n");
		    return 1;
		}
	    }
	}

	if(pbc[1]==0){ // Y surface
	    if(find_keyword(incar,"BYL",option)!=0){ // ==0 when there is no error
		printf(" INPUT FILE ERROR in reading 'BYL' option... \n");
		return 1;
	    }
	    if(option[0]=='A' || option[0]=='a')
		i=0;
	    else if(option[0]=='C' || option[0]=='c')
		i=1;
	    else if(option[0]=='F' || option[0]=='f')
		i=2;
	    else{
		printf(" INPUT ERROR for 'BYL' option... \n");
		return 1;
	    }
	    // Parameters for boundary conditoin
	    if(i==0){ // Adiabatic
		if(set_wetting(incar,&tbc[2][1],"WYL",mode)!=0) // Wetting angle
		    return 1;
		tbc[2][0]=0;
	    }else if(i==1){ // Convection
		if(find_keyword(incar,"TYL",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in reading 'TYL' option... \n");
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    tbc[2][0]=-1*atof(temp);
		else{
		    printf(" INPUT ERROR for 'TYL' option... \n");
		    return 1;
		}
		if(find_keyword(incar,"hYL",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in reading 'hYL' option... \n");
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    tbc[2][1]=-1*atof(temp);
		else{
		    printf(" INPUT ERROR for 'hYL' option... \n");
		    return 1;
		}
	    }else if(i==2){ // Fixed T
		if(set_wetting(incar,&tbc[2][1],"WYL",mode)!=0) // Wetting angle
		    return 1;
		if(find_keyword(incar,"TYL",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in reading 'TYL' option... \n");
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    tbc[2][0]=atof(temp);
		else{
		    printf(" INPUT ERROR for 'TYL' option... \n");
		    return 1;
		}
	    }

	    if(find_keyword(incar,"BYU",option)!=0){ // ==0 when there is no error
		printf(" INPUT FILE ERROR in reading 'BYU' option... \n");
		return 1;
	    }
	    if(option[0]=='A' || option[0]=='a')
		i=0;
	    else if(option[0]=='C' || option[0]=='c')
		i=1;
	    else if(option[0]=='F' || option[0]=='f')
		i=2;
	    else{
		printf(" INPUT ERROR for 'BYU' option... \n");
		return 1;
	    }
	    // Parameters for boundary conditoin
	    if(i==0){ // Adiabatic
		if(set_wetting(incar,&tbc[3][1],"WYU",mode)!=0) // Wetting angle
		    return 1;
		tbc[3][0]=0;
	    }else if(i==1){ // Convection
		if(find_keyword(incar,"TYU",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in reading 'TYU' option... \n");
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    tbc[3][0]=-1*atof(temp);
		else{
		    printf(" INPUT ERROR for 'TYU' option... \n");
		    return 1;
		}
		if(find_keyword(incar,"hYU",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in reading 'hYU' option... \n");
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    tbc[3][1]=-1*atof(temp);
		else{
		    printf(" INPUT ERROR for 'hYU' option... \n");
		    return 1;
		}
	    }else if(i==2){ // Fixed T
		if(set_wetting(incar,&tbc[3][1],"WYU",mode)!=0) // Wetting angle
		    return 1;
		if(find_keyword(incar,"TYU",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in reading 'TYU' option... \n");
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    tbc[3][0]=atof(temp);
		else{
		    printf(" INPUT ERROR for 'TYU' option... \n");
		    return 1;
		}
	    }
	}

	if(pbc[2]==0){ // Z surface
	    if(find_keyword(incar,"BZL",option)!=0){ // ==0 when there is no error
		printf(" INPUT FILE ERROR in reading 'BZL' option... \n");
		return 1;
	    }
	    if(option[0]=='A' || option[0]=='a')
		i=0;
	    else if(option[0]=='C' || option[0]=='c')
		i=1;
	    else if(option[0]=='F' || option[0]=='f')
		i=2;
	    else{
		printf(" INPUT ERROR for 'BZL' option... \n");
		return 1;
	    }
	    // Parameters for boundary conditoin
	    if(i==0){ // Adiabatic
		if(set_wetting(incar,&tbc[4][1],"WZL",mode)!=0) // Wetting angle
		    return 1;
		tbc[4][0]=0;
	    }else if(i==1){ // Convection
		if(find_keyword(incar,"TZL",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in reading 'TZL' option... \n");
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    tbc[4][0]=-1*atof(temp);
		else{
		    printf(" INPUT ERROR for 'TZL' option... \n");
		    return 1;
		}
		if(find_keyword(incar,"hZL",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in reading 'hZL' option... \n");
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    tbc[4][1]=-1*atof(temp);
		else{
		    printf(" INPUT ERROR for 'hZL' option... \n");
		    return 1;
		}
	    }else if(i==2){ // Fixed T
		if(set_wetting(incar,&tbc[4][1],"WZL",mode)!=0) // Wetting angle
		    return 1;
		if(find_keyword(incar,"TZL",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in reading 'TZL' option... \n");
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    tbc[4][0]=atof(temp);
		else{
		    printf(" INPUT ERROR for 'TZL' option... \n");
		    return 1;
		}
	    }

	    if(find_keyword(incar,"BZU",option)!=0){ // ==0 when there is no error
		printf(" INPUT FILE ERROR in reading 'BZU' option... \n");
		return 1;
	    }
	    if(option[0]=='A' || option[0]=='a')
		i=0;
	    else if(option[0]=='C' || option[0]=='c')
		i=1;
	    else if(option[0]=='F' || option[0]=='f')
		i=2;
	    else{
		printf(" INPUT ERROR for 'BZU' option... \n");
		return 1;
	    }
	    // Parameters for boundary conditoin
	    if(i==0){ // Adiabatic
		if(set_wetting(incar,&tbc[5][1],"WZU",mode)!=0) // Wetting angle
		    return 1;
		tbc[5][0]=0;
	    }else if(i==1){ // Convection
		if(find_keyword(incar,"TZU",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in reading 'TZU' option... \n");
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    tbc[5][0]=-1*atof(temp);
		else{
		    printf(" INPUT ERROR for 'TZU' option... \n");
		    return 1;
		}
		if(find_keyword(incar,"hZU",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in reading 'hZU' option... \n");
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    tbc[5][1]=-1*atof(temp);
		else{
		    printf(" INPUT ERROR for 'hZU' option... \n");
		    return 1;
		}
	    }else if(i==2){ // Fixed T
		if(set_wetting(incar,&tbc[5][1],"WZU",mode)!=0) // Wetting angle
		    return 1;
		if(find_keyword(incar,"TZU",option)!=0){ // ==0 when there is no error
		    printf(" INPUT FILE ERROR in reading 'TZU' option... \n");
		    return 1;
		}
		tok(option,temp);
		if(is_valid_digits(temp))
		    tbc[5][0]=atof(temp);
		else{
		    printf(" INPUT ERROR for 'TZU' option... \n");
		    return 1;
		}
	    }
	}

	return 0;
}

int set_wetting(FILE *incar,double *angle,char word[],double mode){
	char temp[100],line[100];

	if(mode==0){ // No nucl. during the simulation...
	    *angle=180.0;
	    return 0;
	}

	if(find_keyword(incar,word,line)!=0) // ==0 when there is no error
	    *angle=180.0;
	else{
	   tok(line,temp);
	    if(is_valid_digits(temp))
		*angle=atof(temp);
	    else{
		printf(" INPUT ERROR IN READING WETTING ANGLE OPTION... \n");
		return 1;
	    }
	}
	return 0;
}

int is_valid_digits(const char *str) {
    int has_digits = 0;    // 숫자가 있는지 확인
    int has_dot = 0;       // 소수점이 있는지 확인
    int has_e = 0;         // 지수 표현이 있는지 확인

    // 문자열의 첫 문자가 부호일 경우 처리
    if (*str == '+' || *str == '-') {
        str++;
    }

    while (*str) {
        if (isdigit(*str)) {
            has_digits = 1; // 숫자가 발견되면 true
        } else if (*str == '.') {
            if (has_dot || has_e) {
                return 0; // 소수점이 이미 있거나 지수 표현 뒤에 소수점이 오면 false
            }
            has_dot = 1; // 소수점이 발견됨
        } else if (*str == 'e' || *str == 'E') {
            if (has_e || !has_digits) {
                return 0; // 지수 표현이 이미 있거나 숫자가 없으면 false
            }
            has_e = 1; // 지수 표현이 발견됨
            str++; // 지수 기호 뒤에 부호가 올 수 있음
            if (*str == '+' || *str == '-') {
                str++; // 부호 처리
            }
        } else {
            return 0; // 숫자, 소수점, 지수 표현이 아닌 경우 false
        }
        str++; // 다음 문자로 이동
    }
    return has_digits; // 숫자가 있었는지 여부 반환
}

void example_incar(){ // Create example INCAR file
	FILE *incar;

	incar=fopen("example.in","w");
	fprintf(incar,"\
##################################################################################################\n\
# SHIMES INCAR FORMAT\n\
# YES = Y, NO = N\n\
# USE \"#\" FOR FOOTNOTES\n\
# NO MORE THAN 100 CHARs\n\
##################################################################################################\n\
# SIMCARD\n\
#\n\
# MDB :  NAME OF MATERIAL INFORMATION DB FILE\n\
#        (SET 'MDB = NONE' TO MANUALLY INPUT MATERIAL INFORMATION)\n\
# SYS :  NAME OF SYSTEM IN DB FILE\n\
# STR :  NAME OF INITIAL MICROSTRUCTURE FILE\n\
# DX  :  UNIT LENGTH (LENGTH OF ONE SIMULATION VOXEL, unit m)\n\
# MODE:  SIMULATION TYPE (CONV: Conventional processes vs. AM: additive manufacturing)\n\
# REAL:  REALISTIC TIME ASSIGNMENT\n\
# PBC :  PERIODIC BOUNDARY CONDITION (XYZ, 1 for on, 0 for off)\n\
#        ex) 111 for setting PBC in all direciton, 100 for setting PBC only in X direction\n\
# CRI :  CRITERIA TO STOP SIMULATION (Available options are: MCS, SEC, LIQ)\n\
#        (MCS = # of MC steps, SEC = realistic time, LIQ = liquid fraction to be 0)\n\
#        (IF REAL = N, AUTOMATICALLY CRI = MCS)\n\
# NSTEP: NUMBER OF TOTAL MCS or SEC TO STOP SIMULATION\n\
# NDATA: RECORD SIMULATION DATA TO 'result.dat' IN EVERY NDATA MCS or SEC\n\
# NPRNT: PRINT MICROSTRUCTURE DATA IN EVERY NPRNT MCS or SEC\n\
#        (NOTE: IF CRI = SEC, MINIMUM ACCURACY OF N*** OPTIONS IS 1E-3 sec...)\n\
#\n\
MDB  = material.sdb\n\
SYS  = fccNi\n\
STR  = str.ini\n\
DX   = 1E-5 m\n\
MODE = CONV\n\
REAL = Y\n\
PBC  = 111\n\
CRI  = MCS\n\
#CRI = LIQ\n\
NSTEP = 2000 MCS\n\
NDATA =  100 MCS\n\
NPRNT =  100 MCS\n\
##################################################################################################\n\
# FDMCARD\n\
#\n\
# T0  :  INITIAL TEMPERATURE OF WHOLE SAMPLE\n\
#        (NOTE: JUST TYPE THE NAME OF *.tdat TEMPERATURE PROFILE FILE TO LOAD THE T PROFILE)\n\
#        ex) T0 = initial.tdat\n\
# FDM :  CONSIDER HEAT TRANSFER? (FDM = Y TO SOLVE HEAT TRANSFER EQUATION USING FDM)\n\
# FACL:  FACTOR FOR ACCELERATION (BY REDUCING GRID LENGTH, INTEGER)\n\
#        ex) FACL = 5, then Lx = 100 -> Lx' = 20\n\
#        (NOTE: (Lx-1), (Ly-1), (Lz-1) SHOULD BE DIVISIBLE BY FACL AND FACL SHOULD NOT EXCEED THE NUMBER OF VOXELS IN A SPECIFIC DIRECTION)\n\
# BAN :  BOUNDARY COUNDTION FOR SOLVING HEAT TRANSFER EQUATION (A = X,Y,Z, N = U,L)\n\
#        ex) BXU = BC AT UPPER-X SURFACE, BZL = BC AT LOWER-Z SURFACE\n\
#        AVAILABLE OPTIONS ARE -\n\
#                    A (ADIABATIC), C (CONVECTION; Robin cond'), F (FIXED T; Dirichlet cond')\n\
#        ex) BXU = A -> ADIABATIC FOR UPPER-X SURFACE\n\
#        (NOTE: IT IS OKAY TO IGNORE IF PBC IS SET IN THAT DIRECTION.)\n\
# TAN :  TEMPERATURE OF OUTER MEDIUM (A = X,Y,Z, N = U,L // VALID ONLY WHEN B** = C or F, unit K)\n\
# hAN :  HEAT TRANSFER COFFICIENT FOR OUTER MEDIUM\n\
#        (A = X,Y,Z, N = U,L // VALID ONLY WHEN B** = C, unit K)\n\
#        ex) BXL = A; TXL = 1000.0; hXL = 25.0\n\
#\n\
T0   = 1000.0 K\n\
FDM  = Y\n\
FACL = 5\n\
BXL  = C\n\
TXL  = 300.0 K\n\
hXL  = 25.0  W/m2K\n\
BXU  = F\n\
TXU  = 400.0\n\
BYL  = F\n\
TYL  = 500.0\n\
BYU  = F\n\
TYU  = 600.0\n\
BZL  = F\n\
TZL  = 700.0\n\
BZU  = F\n\
TZU  = 800.0\n\
h##################################################################################################\n\
# FROZCARD (ONLY AVAILABLE WHEN FDM = N)\n\
#\n\
# FROZ:  FROZEN TEMPERATURE APPROXIMATION?\n\
# TINT:  INTERFACE TEMPERATURE? (ONLY AVAILABLE WHEN FROZ = Y OPTION, in K unit)\n\
# TGRA:  TEMPERATURE GRADIENT?  (ONLY AVAILABLE WHEN FROZ = Y OPTION, in K/voxel unit)\n\
#        TGRA > 0 MEANS dT/dx > 0 IN X DIRECTION\n\
#\n\
FROZ = N\n\
TINT = 000.0 K\n\
TGRA = 000.0 K\n\
##################################################################################################\n\
# TEMPCARD (ONLY AVAILABLE WHEN FDM = N & FROZ = N)\n\
#\n\
# T(N):  FUNCTION FOR TEMPERATURE AT N-th TIME RANGE\n\
#        (UNIT: K, t FOR TIME // !!! NO SPACE DURING WRITING THE EQUATION !!!)\n\
#        ex) 1.53E-3*t+1000\n\
# t(N):  VALID TIME RANGE OF T(N) (FROM t(N-1) ~ t(N) // t(0) IS AUTOMATICALLY 0)\n\
#        (NOTE: t(N) = 0 FOR INPUT INFINITY)\n\
# tUNIT: UNIT FOR TIME (MCS or sec // IF REAL = N, AUTOMATICALLY UNIT = MCS)\n\
# ACCL:  ACCELERATION FACTOR FOR BOTH TIME AND ACCEPTANCE PROBABILITY (DEFUALT = 1.0)\n\
#        (WARNING: ACCL X MAXIMUM PROBABILITY SHOULD BE <= 1.0)\n\
#\n\
T(1) = 0.00E-0*t+000.0\n\
t(1) = 0\n\
tUNIT = MCS\n\
ACCL = 1.0\n\
##################################################################################################\n\
# SOLCARD (ONLY AVAILABLE WHEN REAL = Y)\n\
#\n\
# TCS:  CRITICAL SUPERCOOLING WHERE SOLIDIFICATION STARTS TO BE CONSIDERED (UNIT: K)\n\
#        ex) SOLIDIFICATION ONLY OCCURS WHEN T < (T_MELT - TCS)\n\
# LAT:  CONSIDER LATANT HEAT RELEASE? (DEFAULT OPTION = Y)\n\
# NUCL: CONSIDER NUCLEATION? (DEFAULT OPTION = Y)\n\
# PTCL: Volume fraction of IMC particle? (UNIT: Number of particle voxels / Liquid voxels // DEFAULT = 0)\n\
# PANG: Wetting angle of nucleus at IMC particle surface? (0 ~ 180 DEGREE, DEFAULT = 0)\n\
#       (NOTE: IFE b/w liquid vs. particle ~ IFE b/w liquid vs. solid)\n\
# WAN:  WETTING ANGLE OF HETEROGENEOUS NUCLEI AT SURFACE (A = X,Y,Z // N = U,L, UNIT: DEGREE)\n\
#        ex) WXU = W AT UPPER-X SURFACE, WZL = W AT LOWER-Z SURFACE\n\
#        (NOTE: 180.0 FOR HOMOGENEOUS NUCLEATION CONDITION //\n\
#               IF THE OPTION IS NOT AVAILABLE IN THE SCRIPT, 180.0 WILL BE AUTOMATICALLY SET)\n\
# TIME: REALISTIC TIME MODE? (DEFUALT OPTION = 2)\n\
#                           1 = ONLY CONSIDER t_MCS OF GRAIN GROWTH MODEL (NO ACCELERATION)\n\
#                           2 = COMPLEX t_MCS DETERMINATION SCHEME (SMALLEST METHOD)\n\
#                           3 = PRATICAL t_MCS DETERMINATION SCHEME (OPTIMIZATION AT EACH MCS)\n\
# OPT:  (ONLY VALID WHEN TIME == 3) # of MCSs to modifiy length of MCS (DEFAULT = 10 MCS)\n\
# tMAX: MAXIMUM TIMESTEP LENGTH FOR 1 MCS (UNIT: sec, input 0 to make it unlimited)\n\
# tSCA: MAXIMUM SCALING FACTOR TO INCREASE TIMESTEP LENGTH FOR 1 MCS\n\
#       (NOTE: input 0 to make it unlimited // elsewise, tSCA must be >=1...)\n\
#        ex) t_now/t_prev <= tSCA\n\
#\n\
TCS  = 0 K\n\
LAT  = Y\n\
NUCL = Y\n\
PTCL = 0\n\
PANG = 0.0 deg\n\
WXL  = 0.0 deg\n\
WXU  = 180.0 deg\n\
TIME = 3\n\
OPT  = 10\n\
tMAX = 1E-3 sec\n\
tSCA = 1E2\n\
##################################################################################################\n\
# AMCARD (ONLY AVAILABLE WHEN REAL = Y AND MODE = AM)\n\
#\n\
# DIR:  Direction of LASER scan (X,Y)\n\
# POS:  First position to LASER irradiation (X Y, unit voxel from origin)\n\
#       ex) POS = 50 50\n\
# VEL: LASER SCAN VELOCITY (unit m/s)\n\
# MMOD: Melt pool mode\n\
#      (CIRCLE = Circular shape melt pool assumption, TEAR = teardrop shape melt pool assumption, GAUS = Gaussian LASER simulation)\n\
#       NOTE: NO PBC FOR MELT POOL...\n\
### Below options are for MMOD = CIRCLE ###\n\
# RAD:  Initial melt pool radius (unit voxel)\n\
# DEP:  Maximum depth of melt pool (unit voxel, default = RAD)\n\
# TMP:  Melt pool boundary temperature (unit K, default = Tm)\n\
### Below options are for MMOD = TEAR ###\n\
# WID:  Initial melt pool width (unit voxel)\n\
# CTL:  Initial melt pool cap+tail length (unit voxel)\n\
# DEP:  Initial melt pool maximum depth (unit voxel)\n\
# GEO:  Geometry factor por teardrop shape (0 < GEO < 5)\n\
# TMP:  Melt pool boundary temperature (unit K, defulat = Tm)\n\
### Below options are for MMOD = GAUS ###\n\
# SX :  LASER standard deviation in X direction (unit m)\n\
# SY :  LASER standard deviation in Y direction (unit m)\n\
# SZ :  LASER standard deviation in Z direction (unit m)\n\
# POW:  LASER Power, in W unit\n\
# ABS:  Relative LASER absorptivity (0 < ABS < 1)\n\
#\n\
#\n\
DIR  = X\n\
POS  = 20 20\n\
VEL = 500E-3 m/s\n\
MMOD = TEAR\n\
GEO  = 1.0\n\
WID  = 10 vox\n\
CTL  = 20 vox\n\
DEP  = 5  vox\n\
#RAD  = 20 vox\n\
#TMP  = 2000 K\n\
#MMOD = GAUS\n\
SX   = 100E-6 m\n\
SY   = 100E-6 m\n\
POW  = 400 W\n\
ABS  = 0.25\n\
#EOF\n");
	fclose(incar);
	return;
}

void input_log(int dir[],int pbc[],int out[],char ttem[],double melt[],double jc[],double rscale[],double tbc[][2],double am[]){
	FILE* log=NULL;
	char temp[50]={'\0'},ttem2[100]={'\0'},check[3]={'\0'};
	int i;

	log=fopen("log.input","w");  // Input log file
	fprintf(log,"Sampe size = %d X %d X %d\n",dir[0],dir[1],dir[2]);
	fprintf(log,"PBC[XYZ] = [%d%d%d]\n",pbc[0],pbc[1],pbc[2]);
	fprintf(log,"dx = %.2E \n",rscale[1]);
	if(out[3]<=0)
	    fprintf(log,"Total %d MCS\n",out[0]);
	else
	    fprintf(log,"Total %.3f sec\n",(float)out[0]/1000);
	fprintf(log,"\n=========== 1. Grain Growth =========================================");
	fprintf(log,"\n      γGB (Max) = %.2f J/m2  where θcoh ≤%.1f'",jc[1],jc[3]);
	fprintf(log,"\n      Q GB = %.1f J/mol, Tmax = %.f K\n",rscale[3],rscale[2]);
	fprintf(log,"\n                  2     K0MC = %.2E unit2/MCS",rscale[7]);
	fprintf(log,"\n      t (sec) = dx  * -------------------------- * tMC (MCS)");
	fprintf(log,"\n                        K0_e = %.2E m2/sec",rscale[5]);
	fprintf(log,"\n                =  %E sec / MCS * tMC(MCS)\n",rscale[7]/rscale[5]*rscale[1]*rscale[1]);
	fprintf(log,"\n=========== 2. Solid-Liquid =========================================");
	fprintf(log,"\n                              /");
	fprintf(log,"\n                     L       /");
	fprintf(log,"\n                            / γlv = %.2f J/m2",jc[6]);
	fprintf(log,"\n        γsl = %.2f J/m2   /",jc[4]);
	fprintf(log,"\n       --------------------     V");
	fprintf(log,"\n                           \\");
	fprintf(log,"\n                     S      \\ γsv = %.2f J/m2",jc[5]);
	fprintf(log,"\n                             \\\n");
	fprintf(log,"\n      ΔH_fus = %.1f J/mol = %E J/m3 where Vm = %.2E m3/mol",rscale[8],rscale[8]/rscale[9],rscale[9]);
	fprintf(log,"\n      Tm = %.f K, T* (= Tm - Tcs) = %.2f K,Tsuper = %.2f * Tm = %.f K",rscale[2],rscale[2]-rscale[14],SUPERHEATING,SUPERHEATING*rscale[2]);
	fprintf(log,"\n                        Ks_MC");
	fprintf(log,"\n      t (sec) = dx  * ------- ( = %.2E s unit/m MCS) * tMC (MCS)",rscale[10]);
	fprintf(log,"\n                        Ks_e\n");
	fprintf(log,"\n      Maximum time acceleration: t < %E sec && f_scale < %E\n",rscale[12],rscale[13]);
	fprintf(log,"\n=========== 3. Temperature Condition ================================");
	fprintf(log,"\n   Thermal conductivity k = %.1f W/m K",melt[3]);
	fprintf(log,"\n   Heat capacity Cp = %.3f J/K-mol\n",melt[4]);
	sprintf(ttem2,"%s",ttem);
	tok(ttem2,check);

	if(melt[2]==0){ // no FDM
	    while(1){
	    	if(ttem2[0]=='\0')
		    break;
		tok(ttem2,temp);
		fprintf(log,"T(K) = %-15s",temp);
		tok(ttem2,temp);
		if(temp[0]=='0'){
		    if(check[0]=='0')
			fprintf(log," (unit: MCS) ");
		    else
			fprintf(log," (unit: sec) ");
		    fprintf(log,"...until the end of simulation.\n");
		    break;
		}else{
		    if(check[0]=='0')
			fprintf(log,"until %s MCS...\n",temp);
		    else
			fprintf(log,"until %s sec...\n",temp);
		}
	    }
	}else{ // with FDM
	    fprintf(log,"\n  dt_FDM = %E s (×%d = %E s if applicable)\n",melt[1],(int)(melt[7]*melt[7]),melt[1]*melt[7]*melt[7]);
	    fprintf(log,"  Boundary   Condition    Details\n");
	    fprintf(log,"  X  (up)  | ");
	    if(pbc[0]!=0){
		fprintf(log,"           |    -\n");
		fprintf(log,"  X (down) | ");
		fprintf(log,"           |    -\n");
	    }else{
		if(tbc[1][0]==0)
		    fprintf(log,"Adiabatic  |    -\n");
		else if(tbc[1][0]<0)
		    fprintf(log,"Convection | Ts = %d K, h = %.2f W/m2-K\n",-(int)tbc[1][0],tbc[1][1]);
		else
		    fprintf(log,"Heat sink  | Ts = %d K\n",(int)tbc[1][0]);
		fprintf(log,"  X (down) | ");
		if(tbc[0][0]==0)
		    fprintf(log,"Adiabatic  |    -\n");
		else if(tbc[0][0]<0)
		    fprintf(log,"Convection | Ts = %d K, h = %.2f W/m2-K\n",-(int)tbc[0][0],tbc[0][1]);
		else
		    fprintf(log,"Heat sink  | Ts = %d K\n",(int)tbc[0][0]);
	    }
	    fprintf(log,"  Y  (up)  | ");
	    if(pbc[1]!=0){
		fprintf(log,"           |    -\n");
		fprintf(log,"  Y (down) | ");
		fprintf(log,"           |    -\n");
	    }else{
		if(tbc[3][0]==0)
		    fprintf(log,"Adiabatic  |    -\n");
		else if(tbc[3][0]<0)
		    fprintf(log,"Convection | Ts = %d K, h = %.2f W/m2-K\n",-(int)tbc[3][0],tbc[3][1]);
		else
		    fprintf(log,"Heat sink  | Ts = %d K\n",(int)tbc[3][0]);
		fprintf(log,"  Y (down) | ");
		if(tbc[2][0]==0)
		    fprintf(log,"Adiabatic  |    -\n");
		else if(tbc[2][0]<0)
		    fprintf(log,"Convection | Ts = %d K, h = %.2f W/m2-K\n",-(int)tbc[2][0],tbc[2][1]);
		else
		    fprintf(log,"Heat sink  | Ts = %d K\n",(int)tbc[2][0]);
	    }
	    fprintf(log,"  Z  (up)  | ");
	    if(pbc[2]!=0){
		fprintf(log,"           |    -\n");
		fprintf(log,"  Z (down) | ");
		fprintf(log,"           |    -\n");
	    }else{
		if(tbc[5][0]==0)
		    fprintf(log,"Adiabatic  |    -\n");
		else if(tbc[5][0]<0)
		    fprintf(log,"Convection | Ts = %d K, h = %.2f W/m2-K\n",-(int)tbc[5][0],tbc[5][1]);
		else
		    fprintf(log,"Heat sink  | Ts = %d K\n",(int)tbc[5][0]);
		fprintf(log,"  Z (down) | ");
		if(tbc[4][0]==0)
		    fprintf(log,"Adiabatic  |    -\n");
		else if(tbc[4][0]<0)
		    fprintf(log,"Convection | Ts = %d K, h = %.2f W/m2-K\n",-(int)tbc[4][0],tbc[4][1]);
		else
		    fprintf(log,"Heat sink  | Ts = %d K\n",(int)tbc[4][0]);
	    }
	}
	fprintf(log,"\n======================================================================\n\n");

	if(rscale[0]==2){
		fprintf(log,"\n## Melt pool information ##\n"); // Write input information
		if(am[1]<0)
		    fprintf(log," X direction, Speed(v) = ");
		else
		    fprintf(log," Y direction, Speed(v) = ");
		fprintf(log," %E m/s, start from: ( %d , %d )\n",am[8],(int)am[6],(int)am[7]);
		if(am[0]>0){
		    fprintf(log," Melt pool temperautre is calculated by Gaussian,");
		}else if(am[0]<0){
		    fprintf(log," A circular shaped melt pool is assumed,");
		}else{ //if(am[0]==0){
		    fprintf(log," A teardrop shaped melt pool is assumed,");
		    fprintf(log,"  Width(W) = %d sites = %E m\n",(int)am[2],am[2]*rscale[1]);
		    fprintf(log,"  Cap+tail length(c) = %d sites = %E m\n",(int)am[3],am[3]*rscale[1]);
		    fprintf(log,"  Depth(D) = %d sites = %E m\n",(int)am[4],am[4]*rscale[1]);
		    fprintf(log,"  Melt pool geometry factor(n) = %.3f\n",am[5]);
		}
	}
	fclose(log);
	return;
}

