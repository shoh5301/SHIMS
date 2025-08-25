#include "Database.h"

void DB_main(double melt[],double jc[],double rscale[]){
	int i;

	printf(" Which module? (1 = read database // 2 = create example database file)\n >> ");
	scanf("%d",&i);

	if(i==1)
	    Database(melt,jc,rscale);
	else
	    DB_example();

	return;
}

void DB_example(){
	FILE* db=NULL;

	db=fopen("example.sdb","w");
	fprintf(db,"\
#NAME  Tliq  Tsol  GBEmax Aniso   Qmob    K0_e    K0_MC  beta_SL  dHfus  Vm_solid E_IFsl E_IFsv E_IFlv k_heat Cp(at_mp) IFE_aniso T0\n\
FDMON      1    1   1.00     15      1       1     1E12        1      1     1E-6      1      1      1   500   30        1.0\n\
fccNi   1728 1728  1.732     15 295140 0.05319  3925335   3.0E-3  17160  6.99E-6  0.356  2.104  1.750  88.5   38.14     0.2\n\
Al4Cu    923  845   0.65     15  98220 1.63E-5  11627.6   5.0E-2  10770  9.77E-6  0.154  1.032  0.865   238   37.112    0.2       903\n\
EOF\n\
T(K) GBEmax(J/m2) Aniso(degree) Q(J/mol) D0(m2/s) K0_e(m2/s) K0_MC (site2/s) dHfus(J/mol) Vm (m3/mol-atom) E_IF(J/m2) k_heat (W/m K) Cp (J/K-mol)\n");
	fclose(db);
	printf(" \"example.sdb\" file created...\n");
	
	return;
}

void Database(double melt[],double jc[],double rscale[]){ // Reading material.sdb
	char buffer[200]="",ele[10]="",temp[10]="";
	int i=1,check=0,expo=0;
	double dt=0,tmp;
	FILE* dbfile=NULL;

	printf("DB file name?\n >> ");
	scanf("%s",buffer);

	dbfile=fopen(buffer,"r");

	if(dbfile==NULL){
		printf(" ## No database file (%s)... ##\n\n",buffer);
		return;
	}
	printf("\n * Available materials are...\n");
	fgets(buffer,sizeof(buffer),dbfile);

	while(fgets(buffer,sizeof(buffer),dbfile)!=NULL){ // Reading list of system
		if(strcmp(buffer,"EOF\n")==0)
			break;
		tok(buffer,ele);
		printf(" #%d. %s \n",i++,ele);
	}

	while(check==0){
	    rewind(dbfile);
	    fgets(buffer,sizeof(buffer),dbfile);

	    printf("\nWhich element? (quit = q)\n >> ");
	    scanf("%s",ele);
	    printf("\n");
	    if(strcmp(ele,"q")==0){
		fclose(dbfile);
		return;
	    }
	    while(fgets(buffer,sizeof(buffer),dbfile)!=NULL){ // Reading list of system
		if(strcmp(buffer,"EOF\n")==0)
			break;
		tok(buffer,temp);
		if(strcmp(ele,temp)==0){
			check=1;
			break;
		}
	    }
	if(check==0)
		printf("No material system named as \"%s\"...\nPlease input other...\n\n",ele);
	else{ // Reading information
		printf("\n ### Parameters for the < %s > system ###\n",ele);
		arrange_var(buffer,melt,jc,rscale);
	    }
	}

	printf("\n   1. Grain Growth\n");
	printf("\n      GBE (Max) = %.2f J/m2  where theta_coh <= %.1f'",jc[1],jc[3]);
	printf("\n      Q_GB = %.1f J/mol, Tmax = %.f K\n",rscale[3],rscale[2]);
	printf("\n                  2     K0_e = %.2E m2/s\n",rscale[5]);
	printf("\n      t (sec) = dx  * ----------------------- * t_MC (MCS)");
	printf("\n                        K0MC = %.2E unit2/MCS\n",rscale[7]);
	printf("\n   2. Solid-Liquid Transformation\n");
	printf("\n                              /");
	printf("\n                     LIQ     /");
	printf("\n                            / IFE_lv = %.2f J/m2",jc[6]);
	printf("\n        IFE_sl = %.2f J/m2 /",jc[4]);
	printf("\n       --------------------     VAP(GAS)");
	printf("\n                           \\");
	printf("\n                     SOL    \\ IFE_sv = %.2f J/m2",jc[5]);
	printf("\n                             \\\n");
	printf("\n      dH_fus = %.1f J/mol = %E J/m3  where Vm = %.2E m3/mol",rscale[8],rscale[8]/rscale[9],rscale[9]);
	printf("\n      Tm = %.f K, Tsuper = %.2f * Tm = %.f K",rscale[2],SUPERHEATING,SUPERHEATING*rscale[2]);
	printf("\n                       Ks_MC");
	printf("\n      t (sec) = dx  * ------- ( = %.2E s unit/m MCS) * t_MC (MCS)",rscale[10]);
	printf("\n                       Ks_e\n");
	printf("\n   3. Heat Conduction\n");
	printf("\n   Thermal conductivity k = %.1f W/m K",melt[3]);
	printf("\n   Heat capacity Cp = %.3f J/K-mol",melt[4]);
	printf("\n");
	printf("   ======================================================================\n\n");
	printf(" Press anything to continue..\n >> ");
	scanf("%d",&i);

	fclose(dbfile);
	return;
}


int Load_DB(char file[],char sys[],double melt[],double jc[],double rscale[]){ // Reading material.sdb
	FILE* dbfile=fopen(file,"r");
	char buffer[200],temp[20];
	int check=0;

	if(dbfile==NULL){
		printf(" ## No database file (%s)... ##\n",file);
		return 1;
	}

	while(fgets(buffer,sizeof(buffer),dbfile)!=NULL){
	    if(buffer[0]!='#'){
		if(strcmp(buffer,"EOF\n")==0)
			break;
		tok(buffer,temp);
		if(strcmp(sys,temp)==0){
			check=1;
			break;
		}
	    }
	}
	if(check==0){
	   printf("No material system named as '%s' in DB file '%s'...\n",sys,file);
	   fclose(dbfile);
	   return 1;
	}else
	    arrange_var(buffer,melt,jc,rscale);

	fclose(dbfile);
	return 0;
}

void arrange_var(char buffer[],double melt[],double jc[],double rscale[]){
	char temp[20];

	jc[0]=0.5;
//	tok(buffer,temp);
//	melt[5]=atof(temp); //Tboil
//	tok(buffer,temp);
//	melt[0]=atof(temp); //Tmax 
	tok(buffer,temp);
	rscale[2]=atof(temp); //T_liq
	tok(buffer,temp);
	melt[0]=atof(temp); //T_sol
	tok(buffer,temp);
	jc[1]=atof(temp); // GBEmax
	tok(buffer,temp);
	jc[3]=atof(temp); // Anisotropy angle
	tok(buffer,temp);
	rscale[3]=atof(temp); //Qmob
	tok(buffer,temp);
	rscale[5]=atof(temp); //K0_e for GG
	tok(buffer,temp);
	rscale[7]=atof(temp); //K0_mc for GG
	tok(buffer,temp);
	rscale[10]=atof(temp); //beta_SL = KMC_solidi / Ke_solidi  1/(m/s K)
	tok(buffer,temp);
	rscale[8]=atof(temp); //dH_fus
	tok(buffer,temp);
	rscale[9]=atof(temp); //Vm
	tok(buffer,temp);
	jc[4]=atof(temp); //Sol-Liq. Interf. En.
	tok(buffer,temp);
	jc[5]=atof(temp); //Sol-Vap. Interf. En.
	tok(buffer,temp);
	jc[6]=atof(temp); //Liq-Vap. Interf. En.
	tok(buffer,temp);
	melt[3]=atof(temp); //k (thermal conductivity), W/m K = J/s M K
	tok(buffer,temp);
	melt[4]=atof(temp); //Cp (heat capacity) J / mol
	tok(buffer,temp);
	jc[2]=atof(temp); //Growth direction anisotropy
	if(rscale[2]!=melt[0]){
	    tok(buffer,temp);
	    melt[8]=atof(temp);
	}else
	    melt[8]=rscale[2];

	if(melt[8]==0)
	    printf(" @@@ WARNING: Wrong T0 value! // NOTE: T0 would be T_sol < T0 < T_liq @@@\n");

	return;
}
