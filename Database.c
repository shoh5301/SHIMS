#include "Database.h"

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
		jc[0]=0.5;
//		tok(buffer,temp);
//		melt[5]=atof(temp); //Tboil
//		tok(buffer,temp);
//		melt[0]=atof(temp); //Tmax 
		tok(buffer,temp);
		rscale[2]=atof(temp); //Tmelt
		melt[0]=SUPERHEATING*rscale[2]; // superheating limit T
//		melt[0]=rscale[2]; // just Tmelt
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
		rscale[10]=atof(temp); //KMC_solidi / Ke_solidi  1/(m/s K)
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
	    }
	}

	printf("\n   1. Grain Growth\n");
	printf("\n      γGB (Max) = %.2f J/m2  where θcoh ≤%.1f'",jc[1],jc[3]);
	printf("\n      Q GB = %.1f J/mol, Tmax = %.f K\n",rscale[3],rscale[2]);
	printf("\n                  2     K0_e = %.2E m2/s\n",rscale[5]);
	printf("\n      t (sec) = dx  * ----------------------- * tMC (MCS)");
	printf("\n                        K0MC = %.2E unit2/MCS\n",rscale[7]);
	printf("\n   2. Solid-Liquid Transformation\n");
	printf("\n                              /");
	printf("\n                     L       /");
	printf("\n                            / γlv = %.2f J/m2",jc[6]);
	printf("\n        γsl = %.2f J/m2   /",jc[4]);
	printf("\n       --------------------     V");
	printf("\n                           \\");
	printf("\n                     S      \\ γsv = %.2f J/m2",jc[5]);
	printf("\n                             \\\n");
	printf("\n      ΔH_fus = %.1f J/mol = %E J/m3  where Vm = %.2E m3/mol",rscale[8],rscale[8]/rscale[9],rscale[9]);
	printf("\n      Tm = %.f K, Tsuper = %.2f * Tm = %.f K",rscale[2],SUPERHEATING,SUPERHEATING*rscale[2]);
	printf("\n                        Ks_MC\n");
	printf("\n      t (sec) = dx  * ------- ( = %.2E s unit/m MCS) * tMC (MCS)",rscale[10]);
	printf("\n                        Ks_e\n");
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
	}else{
		jc[0]=0.5;
//		tok(buffer,temp);
//		melt[5]=atof(temp); //Tboil
		tok(buffer,temp);
		rscale[2]=atof(temp); //Tmelt
//		melt[0]=rscale[2]; // just Tmelt
		melt[0]=SUPERHEATING*rscale[2]; // Ideal melting superheating limit
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
		rscale[10]=atof(temp); //beta = KMC_solidi / Ke_solidi  1/(m/s K)
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
		jc[2]=atof(temp); //S-L IFE_anisotropy
	}
	fclose(dbfile);
	return 0;
}

