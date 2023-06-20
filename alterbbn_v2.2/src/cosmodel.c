#include "include.h"

/*--------------------------------------------------------------*/

double heff(double Temp, struct relicparam* paramrelic)
{
	int ie;

	if(Temp>= paramrelic->table_eff[0][0]) return paramrelic->table_eff[0][2];

	if(Temp<2.35e-13) return paramrelic->table_eff[275][2];

	ie=1;
	while(Temp<paramrelic->table_eff[ie][0]) ie++;

	double heff1,heff2,logT1,logT2,logT;
	logT=log(Temp);
	heff1=paramrelic->table_eff[ie][2];
	heff2=paramrelic->table_eff[ie-1][2];
	logT1=log(paramrelic->table_eff[ie][0]);
	logT2=log(paramrelic->table_eff[ie-1][0]);

	return (heff2-heff1)/(logT2-logT1)*(logT-logT1)+heff1;
}

/*--------------------------------------------------------------*/

double sgStar(double Temp, struct relicparam* paramrelic)
{
	int ie;

	if(Temp>=paramrelic->table_eff[0][0]) return paramrelic->table_eff[0][1];

	if(Temp<2.35e-13) return paramrelic->table_eff[275][1];

	ie=1;
	while(Temp<paramrelic->table_eff[ie][0]) ie++;

	double sgStar1,sgStar2,logT1,logT2,logT;
	logT=log(Temp);
	sgStar1=paramrelic->table_eff[ie][1];
	sgStar2=paramrelic->table_eff[ie-1][1];
	logT1=log(paramrelic->table_eff[ie][0]);
	logT2=log(paramrelic->table_eff[ie-1][0]);

	return (sgStar2-sgStar1)/(logT2-logT1)*(logT-logT1)+sgStar1;
}

/*--------------------------------------------------------------*/

double geff(double Temp, struct relicparam* paramrelic)
{
	double heff0=heff(Temp,paramrelic);

	return pow(heff0/sgStar(Temp,paramrelic)*(1.+(heff(Temp*1.001,paramrelic)-heff(Temp*0.999,paramrelic))/0.006/heff0),2.);
}

/*--------------------------------------------------------------*/

void Init_cosmomodel(struct relicparam* paramrelic)
{
	paramrelic->full_comput=0;

	paramrelic->solver=1; /* 1=logarithmic, 2=linear */

	paramrelic->failsafe=1;         // 0=fast, 1=precise (default), ... See stand_cosmo.c
	paramrelic->err=0;
    paramrelic->Tinit=27.;          // Starting at T = 27 x 10^9 K as default
    paramrelic->Tnudec=27.;           // Neutrino decoupling T = 27 x 10^9 K as default
    paramrelic->eta0=6.10e-10;      // Baryon-to-photon ratio (Planck 2015 results XIII)
    paramrelic->Nnu=3.046;          // Number of SM neutrinos, e+- reheating included
    paramrelic->dNnu=0.;            // Number of extra neutrino species (e.g. sterile neutrinos)
    paramrelic->life_neutron=880.2; // Neutron lifetime (PDG2018)
    paramrelic->life_neutron_error=1.0; // Neutron lifetime uncertainty (PDG2017)
	paramrelic->xinu1=0.;
	paramrelic->xinu2=0.;
	paramrelic->xinu3=0.;
	paramrelic->beta_samples=100;	// how accurately to model n<-->p beta reactions

	paramrelic->b_cdm_ratio=0.02242/0.11933;   // current baryon to cold dark matter density ratio (Planck 2018 results VI)

    paramrelic->wimp=0;
	//paramrelic->m_chi=paramrelic->g_chi=paramrelic->SMC_wimp=paramrelic->selfConjugate=paramrelic->fermion=paramrelic->EM_coupled=paramrelic->neut_coupled=paramrelic->neuteq_coupled=0;
	paramrelic->m_chi=0;
	paramrelic->g_chi=0;
	paramrelic->SMC_wimp=0;
	paramrelic->selfConjugate=0;
	paramrelic->fermion=0;
	paramrelic->EM_coupled=0;
	paramrelic->neut_coupled=0;
	paramrelic->neuteq_coupled=0;

	paramrelic->fierz=0.;			// no fierz interface in the standard model
	paramrelic->B_chi=0.;			// default is no branching to dark matter m_p < m_chi < m_n

	paramrelic->dd0=paramrelic->ndd=paramrelic->Tdend=paramrelic->Tddeq=0.;
	paramrelic->sd0=paramrelic->nsd=paramrelic->Tsend=0.;
	paramrelic->nt0=paramrelic->nnt=paramrelic->Tnend=0.;
	paramrelic->Sigmad0=paramrelic->nSigmad=paramrelic->TSigmadend=0.;
	paramrelic->Sigmarad0=paramrelic->nSigmarad=paramrelic->TSigmaradend=0.;
	paramrelic->coupd=0;

	paramrelic->mgravitino=0.;

	paramrelic->phi_model=0;
	paramrelic->eta_phi=paramrelic->Gamma_phi=paramrelic->n_phi=paramrelic->rhot_phi_Tmax=paramrelic->Tphi0=paramrelic->rhot_phi0=0.;

	paramrelic->vs_model=0;
	paramrelic->Gamma_vs=paramrelic->ns0=0.;
	//paramrelic->eta_vs=paramrelic->Gamma_vs=paramrelic->rhot_vs_Tmax=paramrelic->Tvs0=paramrelic->rhot_vs0=0.;

  paramrelic->entropy_model=1;
	paramrelic->energy_model=1;

 	paramrelic->mgravitino=paramrelic->relicmass=0.;

	paramrelic->full_comput=paramrelic->scalar=0;
	paramrelic->Tfo=0.;

	paramrelic->Tmax=50.;

    Init_modeleff(2,paramrelic);

    paramrelic->use_table_rhoPD=paramrelic->size_table_rhoPD=0;

    paramrelic->constraints=2;

	return;
}

/*--------------------------------------------------------------*/

void Init_cosmomodel_param(double eta, double Nnu, double dNnu, double life_neutron, double life_neutron_error, double xinu1,
                           double xinu2, double xinu3, struct relicparam* paramrelic)
/* modifies the values of the baryon-to-photon ratio eta, the number of SM neutrinos Nnu, extra neutrino species dNnu
 *  and the neutron lifetime life_neutron */
{
    paramrelic->eta0=eta;
    paramrelic->Nnu=Nnu;
    paramrelic->dNnu=dNnu;
    paramrelic->life_neutron=life_neutron;
    paramrelic->life_neutron_error=life_neutron_error;
    paramrelic->xinu1=xinu1;
    paramrelic->xinu2=xinu2;
    paramrelic->xinu3=xinu3;
	paramrelic->beta_samples=100;	// how accurately to model n<-->p beta reactions
    return;
}

/*--------------------------------------------------------------*/

void Init_modeleff(int model_eff, struct relicparam* paramrelic)
{
	int ie,je;

	if(model_eff==1)
	{
		const double tableA[276][3]=
		{
#include "sgStar_heff/sgStar_heff_A.tab"
  		};
		for(ie=0;ie<=275;ie++) for(je=0;je<=2;je++) paramrelic->table_eff[ie][je]=tableA[ie][je];
	}
	else if(model_eff==2)
	{
		const double tableB[276][3]=
		{
#include "sgStar_heff/sgStar_heff_B.tab"
  		};
		for(ie=0;ie<=275;ie++) for(je=0;je<=2;je++) paramrelic->table_eff[ie][je]=tableB[ie][je];
	}
	else if(model_eff==3)
	{
	const double tableB2[276][3]=
		{
#include "sgStar_heff/sgStar_heff_B2.tab"
  		};
		for(ie=0;ie<=275;ie++) for(je=0;je<=2;je++) paramrelic->table_eff[ie][je]=tableB2[ie][je];
		}
	else if(model_eff==4)
	{
		const double tableB3[276][3]=
		{
#include "sgStar_heff/sgStar_heff_B3.tab"
  		};
		for(ie=0;ie<=275;ie++) for(je=0;je<=2;je++) paramrelic->table_eff[ie][je]=tableB3[ie][je];
	}
	else if(model_eff==5)
	{
		const double tableC[276][3]=
		{
#include "sgStar_heff/sgStar_heff_C.tab"
  		};
		for(ie=0;ie<=275;ie++) for(je=0;je<=2;je++) paramrelic->table_eff[ie][je]=tableC[ie][je];
	}
	else if(model_eff==6)
	{
		const double tableBonn[276][3]=
		{
#include "sgStar_heff/sgStar_heff_Bonn.tab"
  		};
		for(ie=0;ie<=275;ie++) for(je=0;je<=2;je++) paramrelic->table_eff[ie][je]=tableBonn[ie][je];
	}
	else
	{
		const double tableold[276][3]=
		{
#include "sgStar_heff/sgStar_heff_old.tab"
  		};
		for(ie=0;ie<=275;ie++) for(je=0;je<=2;je++) paramrelic->table_eff[ie][je]=tableold[ie][je];
	}
	return;
}

/*--------------------------------------------------------------*/

void Init_wimp(double mass_wimp, int EM_coupled, int neut_coupled, int neuteq_coupled, int fermion, int selfConjugate, double g_chi, struct relicparam* paramrelic)
/* modifies the parameters of an included light WIMP */
{
    paramrelic->m_chi=mass_wimp;
    paramrelic->g_chi=g_chi;
    paramrelic->fermion=fermion;
    paramrelic->EM_coupled=EM_coupled;
    paramrelic->neut_coupled=neut_coupled;
    paramrelic->neuteq_coupled=neuteq_coupled;
    paramrelic->wimp=1;
    paramrelic->selfConjugate=selfConjugate;
    return;
}

/*--------------------------------------------------------------*/

void Init_dark_density(double dd0, double ndd, double T_end, struct relicparam* paramrelic)
{
	if(paramrelic->phi_model)
	{
		paramrelic->energy_model=0;
		return;
	}

	paramrelic->energy_model=1;
	paramrelic->dd0=dd0;
	paramrelic->ndd=ndd;
	paramrelic->Tdend=T_end;

	paramrelic->Tddeq=0.;

    paramrelic->use_table_rhoPD=paramrelic->size_table_rhoPD=0;

	return;
}

/*--------------------------------------------------------------*/

void Init_dark_density2(double ndd, double Tddeq, double T_end, struct relicparam* paramrelic)
{
	if(paramrelic->phi_model)
	{
		paramrelic->energy_model=0;
		return;
	}

	paramrelic->energy_model=2;
	paramrelic->ndd=ndd;
	paramrelic->Tddeq=Tddeq;
	paramrelic->Tdend=T_end;

	paramrelic->dd0=0.;

    paramrelic->use_table_rhoPD=paramrelic->size_table_rhoPD=0;

	return;
}

/*--------------------------------------------------------------*/

void Init_dark_coupling(int coupD, struct relicparam* paramrelic)
{
	paramrelic->coupd=coupD;

	return;
}

/*--------------------------------------------------------------*/

void Init_quintessence(double T12, double n2, double T23, double n3, double T34, double n4, struct relicparam* paramrelic)
{
	if(paramrelic->phi_model)
	{
		paramrelic->energy_model=0;
		return;
	}

	paramrelic->energy_model=3;
	paramrelic->dd0=0.;
    paramrelic->use_table_rhoPD=paramrelic->size_table_rhoPD=0;

	paramrelic->quintn2=n2;
	paramrelic->quintn3=n3;
	paramrelic->quintn4=n4;
	paramrelic->quintT12=T12;
	paramrelic->quintT23=T23;
	paramrelic->quintT34=T34;

	return;
}

/*--------------------------------------------------------------*/

void Init_dark_entropy(double sd0, double nsd, double T_end, struct relicparam* paramrelic)
{
	if(paramrelic->phi_model) return;

	paramrelic->sd0=sd0;
	paramrelic->nsd=nsd;
	paramrelic->Tsend=T_end;

    paramrelic->use_table_rhoPD=paramrelic->size_table_rhoPD=0;

	return;
}

/*--------------------------------------------------------------*/

void Init_dark_entropySigmaD(double Sigmad0, double nSigmad, double T_end, struct relicparam* paramrelic)
{
	if(paramrelic->phi_model) return;

	paramrelic->Sigmad0=Sigmad0;
	paramrelic->nSigmad=nSigmad;
	paramrelic->TSigmadend=T_end;

    paramrelic->use_table_rhoPD=paramrelic->size_table_rhoPD=0;

	return;
}

/*--------------------------------------------------------------*/

void Init_entropySigmarad(double Sigmarad0, double nSigmarad, double T_end, struct relicparam* paramrelic)
{
	if(paramrelic->phi_model) return;

	paramrelic->Sigmarad0=Sigmarad0;
	paramrelic->nSigmarad=nSigmarad;
	paramrelic->TSigmaradend=T_end;

    paramrelic->use_table_rhoPD=paramrelic->size_table_rhoPD=0;

	return;
}

/*--------------------------------------------------------------*/

void Init_nonthermal(double nt0, double nnt, double T_end, struct relicparam* paramrelic)
{
	if(paramrelic->phi_model) return;

	paramrelic->nt0=nt0;
	paramrelic->nnt=nnt;
	paramrelic->Tnend=T_end;

    paramrelic->use_table_rhoPD=paramrelic->size_table_rhoPD=0;

	return;
}
/*--------------------------------------------------------------*/

void Init_gravitino(double mgravitino, struct relicparam* paramrelic)
{
	paramrelic->mgravitino=mgravitino;

	return;
}

/*--------------------------------------------------------------*/

void Init_scalarfield(double rhotilde_phi_Tmax, double Tmax, double T_RH, double eta_phi, double n_phi, struct relicparam* paramrelic)
{
	paramrelic->phi_model=1;
	paramrelic->full_comput=1;

	paramrelic->T_RH=T_RH;
	paramrelic->eta_phi=eta_phi;  /* = b/m_phi */
	paramrelic->Gamma_phi=sqrt(4.*pow(pi,3.)*geff(T_RH,paramrelic)/45.)*T_RH*T_RH/Mplanck;
	paramrelic->rhot_phi_Tmax=rhotilde_phi_Tmax;
	paramrelic->Tmax=Tmax;
	paramrelic->n_phi=n_phi;

	paramrelic->rhot_phi0=paramrelic->rhot_phi_Tmax;

	Init_dark_density(0.,0.,0.,paramrelic);
	Init_dark_density2(0.,0.,0.,paramrelic);
	Init_dark_entropy(0.,0.,0.,paramrelic);
	Init_dark_entropySigmaD(0.,0.,0.,paramrelic);
	Init_entropySigmarad(0.,0.,0.,paramrelic);
	Init_nonthermal(0.,0.,0.,paramrelic);
	Init_dark_coupling(0,paramrelic);

    paramrelic->use_table_rhoPD=paramrelic->size_table_rhoPD=0;

	return;
}

/*--------------------------------------------------------------*/

void Init_dark_density_table(double table[2][NTABMAX], int nlines, relicparam* paramrelic)
{
	Init_dark_density(0.,0.,0.,paramrelic);
	Init_dark_density2(0.,0.,0.,paramrelic);
	Init_dark_entropy(0.,0.,0.,paramrelic);
	Init_dark_entropySigmaD(0.,0.,0.,paramrelic);
	Init_entropySigmarad(0.,0.,0.,paramrelic);
	Init_nonthermal(0.,0.,0.,paramrelic);

	paramrelic->use_table_rhoPD=1;

	int ie,je;

	paramrelic->size_table_rhoPD=nlines;
	for(ie=0;ie<=1;ie++) for(je=0;je<nlines;je++) paramrelic->table_rhoPD[ie][je]=table[ie][je];

	return;
}


/*--------------------------------------------------------------*/

void Init_neutron_decay(double tau, double tau_err, double fierz, double m_chi, double B_chi, relicparam* paramrelic)
{
    paramrelic->life_neutron=tau;
    paramrelic->life_neutron_error=tau_err; // Neutron lifetime uncertainty (PDG2017)
    paramrelic->wimp=0;
	paramrelic->m_chi=m_chi;
	paramrelic->g_chi=0;
	paramrelic->SMC_wimp=0;
	paramrelic->selfConjugate=0;
	paramrelic->fermion=0;
	paramrelic->EM_coupled=0;
	paramrelic->neut_coupled=0;
	paramrelic->neuteq_coupled=0;
	paramrelic->fierz=fierz;			// no fierz interface in the standard model
	paramrelic->B_chi=B_chi;			// default is no branching to dark matter m_p < m_chi < m_n
	paramrelic->beta_samples=1000;		// if we are gonna simulate beta decay, we're gonna simulate it.
}


/*--------------------------------------------------------------*/

void read_csv(int row, int col, char *filename, double **data)
{
	FILE *file;
	file = fopen(filename, "r"); //what's the r?
    
    if (file == NULL) {printf("missing %s\n", filename);}

	int i = 0;
	char line[4098];
	while (fgets(line, 4098, file) && (i<row))
	{
		//double row[ssParams->nreal+1]; //what??
		char *tmp = strdup(line);

		int j = 0;
		const char* tok;
		for (tok = strtok(line, ","); tok && *tok; j++, tok=strtok(NULL, ","))
		{
			data[i][j] = atof(tok);
		}
		free(tmp);
		i++;
	}
}

void Init_vs(char ms_ch[256], char mix_ch[256], double ms_d, double mix_d, double ns0, int row, struct relicparam* paramrelic)
{
    //printf("got here!\n");
	int col = 3;
	char folder[256]= "../";
	paramrelic->vs_model = 1;
	paramrelic->ms = ms_d;
	paramrelic->mix = mix_d;
	paramrelic->ns0 = ns0/pow(1000,3.); //expects it to be passed in terms of MeV^3, but we need to turn it into GeV^3 for use in the code
	paramrelic->rhot_vs0 = ms_d*ns0/pow(1000,4.); //I think it might be easiest and most accurate to just pass in the initial ns value
	ts(ms_d, mix_d, paramrelic);

	char fname_Tcm_rhonu[256];
	char fname_a_rhonu[256];
	char fname_b_rhonu[256];
	char fname_c_rhonu[256];
	char fname_d_rhonu[256];
	char fname_T_dQdt[256];
	char fname_a_dQdt[256];
	char fname_b_dQdt[256];
	char fname_c_dQdt[256];
	char fname_d_dQdt[256];
	char fname_T_n2p[256];
	char fname_a_n2p[256];
	char fname_b_n2p[256];
	char fname_c_n2p[256];
	char fname_d_n2p[256];
	char fname_T_p2n[256];
	char fname_a_p2n[256];
	char fname_b_p2n[256];
	char fname_c_p2n[256];
	char fname_d_p2n[256];

	double **data_Tcm_rhonu;
	double **data_a_rhonu;
	double **data_b_rhonu;
	double **data_c_rhonu;
	double **data_d_rhonu;
	double **data_T_dQdt;
	double **data_a_dQdt;
	double **data_b_dQdt;
	double **data_c_dQdt;
	double **data_d_dQdt;
	double **data_T_n2p;
	double **data_a_n2p;
	double **data_b_n2p;
	double **data_c_n2p;
	double **data_d_n2p;
	double **data_T_p2n;
	double **data_a_p2n;
	double **data_b_p2n;
	double **data_c_p2n;
	double **data_d_p2n;

	//strcat(folder,ms_ch); strcat(folder,"-"); strcat(folder,mix_ch); strcat(folder,"-FullTestNew/"); 
    
    strcat(folder,"alterbbn_v2.2/"); strcat(folder,"CSVFiles/"); /* //CSV_Files is Hannah's files
    

    strcat(folder,ms_ch); strcat(folder,"-"); strcat(folder,mix_ch); strcat(folder,"-FullTestNew/"); strcat(folder, "mass_");       strcat(folder,ms_ch); strcat(folder, "_mix_"); strcat(folder,mix_ch); */
    

    // for Kathryn's code 
	strcpy(fname_Tcm_rhonu, folder); strcat(fname_Tcm_rhonu, "T_rhonu.csv");
	strcpy(fname_a_rhonu, folder); strcat(fname_a_rhonu, "a_rhonu.csv");
	strcpy(fname_b_rhonu, folder); strcat(fname_b_rhonu, "b_rhonu.csv");
	strcpy(fname_c_rhonu, folder); strcat(fname_c_rhonu, "c_rhonu.csv");
	strcpy(fname_d_rhonu, folder); strcat(fname_d_rhonu, "d_rhonu.csv");
	strcpy(fname_T_dQdt, folder); strcat(fname_T_dQdt, "T_dQdt.csv");
	strcpy(fname_a_dQdt, folder); strcat(fname_a_dQdt, "a_dQdt.csv");
	strcpy(fname_b_dQdt, folder); strcat(fname_b_dQdt, "b_dQdt.csv");
	strcpy(fname_c_dQdt, folder); strcat(fname_c_dQdt, "c_dQdt.csv");
	strcpy(fname_d_dQdt, folder); strcat(fname_d_dQdt, "d_dQdt.csv");
	strcpy(fname_T_n2p, folder); strcat(fname_T_n2p, "T_np.csv");
	strcpy(fname_a_n2p, folder); strcat(fname_a_n2p, "a_np.csv");
	strcpy(fname_b_n2p, folder); strcat(fname_b_n2p, "b_np.csv");
	strcpy(fname_c_n2p, folder); strcat(fname_c_n2p, "c_np.csv");
	strcpy(fname_d_n2p, folder); strcat(fname_d_n2p, "d_np.csv");
	strcpy(fname_T_p2n, folder); strcat(fname_T_p2n, "T_pn.csv");
	strcpy(fname_a_p2n, folder); strcat(fname_a_p2n, "a_pn.csv");
	strcpy(fname_b_p2n, folder); strcat(fname_b_p2n, "b_pn.csv");
	strcpy(fname_c_p2n, folder); strcat(fname_c_p2n, "c_pn.csv");
	strcpy(fname_d_p2n, folder); strcat(fname_d_p2n, "d_pn.csv"); /*
    

    //for Hannah's code
   	strcpy(fname_Tcm_rhonu, folder); strcat(fname_Tcm_rhonu, "_Tcm_rhonu.csv");
	strcpy(fname_a_rhonu, folder); strcat(fname_a_rhonu, "_a_rhonu.csv");
	strcpy(fname_b_rhonu, folder); strcat(fname_b_rhonu, "_b_rhonu.csv");
	strcpy(fname_c_rhonu, folder); strcat(fname_c_rhonu, "_c_rhonu.csv");
	strcpy(fname_d_rhonu, folder); strcat(fname_d_rhonu, "_d_rhonu.csv");
	strcpy(fname_T_dQdt, folder); strcat(fname_T_dQdt, "_T_dQdt.csv");
	strcpy(fname_a_dQdt, folder); strcat(fname_a_dQdt, "_a_dQdt.csv");
	strcpy(fname_b_dQdt, folder); strcat(fname_b_dQdt, "_b_dQdt.csv");
	strcpy(fname_c_dQdt, folder); strcat(fname_c_dQdt, "_c_dQdt.csv");
	strcpy(fname_d_dQdt, folder); strcat(fname_d_dQdt, "_d_dQdt.csv");
	strcpy(fname_T_n2p, folder); strcat(fname_T_n2p, "_T_np.csv");
	strcpy(fname_a_n2p, folder); strcat(fname_a_n2p, "_a_np.csv");
	strcpy(fname_b_n2p, folder); strcat(fname_b_n2p, "_b_np.csv");
	strcpy(fname_c_n2p, folder); strcat(fname_c_n2p, "_c_np.csv");
	strcpy(fname_d_n2p, folder); strcat(fname_d_n2p, "_d_np.csv");
	strcpy(fname_T_p2n, folder); strcat(fname_T_p2n, "_T_pn.csv");
	strcpy(fname_a_p2n, folder); strcat(fname_a_p2n, "_a_pn.csv");
	strcpy(fname_b_p2n, folder); strcat(fname_b_p2n, "_b_pn.csv");
	strcpy(fname_c_p2n, folder); strcat(fname_c_p2n, "_c_pn.csv");
	strcpy(fname_d_p2n, folder); strcat(fname_d_p2n, "_d_pn.csv");*/
    
    
    printf("print names now please: \n");
    printf("%s \n", fname_Tcm_rhonu);
    printf("%s \n",fname_a_rhonu);
    printf("%s \n",fname_b_rhonu);
    printf("%s \n",fname_c_rhonu);
    printf("%s \n",fname_d_rhonu);
    printf("%s \n",fname_T_dQdt);
    printf("%s \n",fname_a_dQdt);
    printf("%s \n",fname_b_dQdt);
    printf("%s \n",fname_c_dQdt);
    printf("%s \n",fname_d_dQdt);
    printf("%s \n",fname_T_n2p);
    printf("%s \n",fname_a_p2n);
    printf("%s \n",fname_b_p2n);
    printf("%s \n",fname_c_p2n);
    printf("%s \n",fname_d_p2n); 
    

	//the number of rows for the three data types might not be the same, but it should always be the same relative to each other (for example, row in dQdt will always be one less than in rhonu I think)
	data_Tcm_rhonu = (double **)malloc(row * sizeof(double *));
	data_a_rhonu = (double **)malloc(row * sizeof(double *));
	data_b_rhonu = (double **)malloc(row * sizeof(double *));
	data_c_rhonu = (double **)malloc(row * sizeof(double *));
	data_d_rhonu = (double **)malloc(row * sizeof(double *));
	data_T_dQdt = (double **)malloc((row-1) * sizeof(double *));
	data_a_dQdt = (double **)malloc((row-1) * sizeof(double *));
	data_b_dQdt = (double **)malloc((row-1) * sizeof(double *));
	data_c_dQdt = (double **)malloc((row-1) * sizeof(double *));
	data_d_dQdt = (double **)malloc((row-1) * sizeof(double *));
	data_T_n2p = (double **)malloc(row * sizeof(double *));
	data_a_n2p = (double **)malloc(row * sizeof(double *));
	data_b_n2p = (double **)malloc(row * sizeof(double *));
	data_c_n2p = (double **)malloc(row * sizeof(double *));
	data_d_n2p = (double **)malloc(row * sizeof(double *));
	data_T_p2n = (double **)malloc(row * sizeof(double *));
	data_a_p2n = (double **)malloc(row * sizeof(double *));
	data_b_p2n = (double **)malloc(row * sizeof(double *));
	data_c_p2n = (double **)malloc(row * sizeof(double *));
	data_d_p2n = (double **)malloc(row * sizeof(double *));

	for (int i = 0; i<row; i++) //hmm, row might matter here...
	{
		data_Tcm_rhonu[i] = (double *)malloc(col * sizeof(double));
		data_a_rhonu[i] = (double *)malloc(col * sizeof(double));
		data_b_rhonu[i] = (double *)malloc(col * sizeof(double));
		data_c_rhonu[i] = (double *)malloc(col * sizeof(double));
		data_d_rhonu[i] = (double *)malloc(col * sizeof(double));
		data_T_n2p[i] = (double *)malloc(col * sizeof(double));
		data_a_n2p[i] = (double *)malloc(col * sizeof(double));
		data_b_n2p[i] = (double *)malloc(col * sizeof(double));
		data_c_n2p[i] = (double *)malloc(col * sizeof(double));
		data_d_n2p[i] = (double *)malloc(col * sizeof(double));
		data_T_p2n[i] = (double *)malloc(col * sizeof(double));
		data_a_p2n[i] = (double *)malloc(col * sizeof(double));
		data_b_p2n[i] = (double *)malloc(col * sizeof(double));
		data_c_p2n[i] = (double *)malloc(col * sizeof(double));
		data_d_p2n[i] = (double *)malloc(col * sizeof(double));
	}
	for (int i = 0; i<row-1; i++) //hmm, row might matter here...
	{
		data_T_dQdt[i] = (double *)malloc(col * sizeof(double));
		data_a_dQdt[i] = (double *)malloc(col * sizeof(double));
		data_b_dQdt[i] = (double *)malloc(col * sizeof(double));
		data_c_dQdt[i] = (double *)malloc(col * sizeof(double));
		data_d_dQdt[i] = (double *)malloc(col * sizeof(double));
	}

	read_csv(row, col, fname_Tcm_rhonu, data_Tcm_rhonu);
	read_csv(row, col, fname_a_rhonu, data_a_rhonu);
	read_csv(row, col, fname_b_rhonu, data_b_rhonu);
	read_csv(row, col, fname_c_rhonu, data_c_rhonu);
	read_csv(row, col, fname_d_rhonu, data_d_rhonu);
	read_csv(row-1, col, fname_T_dQdt, data_T_dQdt);
	read_csv(row-1, col, fname_a_dQdt, data_a_dQdt);
	read_csv(row-1, col, fname_b_dQdt, data_b_dQdt);
	read_csv(row-1, col, fname_c_dQdt, data_c_dQdt);
	read_csv(row-1, col, fname_d_dQdt, data_d_dQdt);
	read_csv(row, col, fname_T_n2p, data_T_n2p);
	read_csv(row, col, fname_a_n2p, data_a_n2p);
	read_csv(row, col, fname_b_n2p, data_b_n2p);
	read_csv(row, col, fname_c_n2p, data_c_n2p);
	read_csv(row, col, fname_d_n2p, data_d_n2p);
	read_csv(row, col, fname_T_p2n, data_T_p2n);
	read_csv(row, col, fname_a_p2n, data_a_p2n);
	read_csv(row, col, fname_b_p2n, data_b_p2n);
	read_csv(row, col, fname_c_p2n, data_c_p2n);
	read_csv(row, col, fname_d_p2n, data_d_p2n);

	for (int i=0; i<row-1; i++)
	{
		paramrelic->Tcm_rho[i] = data_Tcm_rhonu[i+1][1];
		paramrelic->arho[i] = data_a_rhonu[i+1][1];
		paramrelic->brho[i] = data_b_rhonu[i+1][1];
		paramrelic->crho[i] = data_c_rhonu[i+1][1];
		paramrelic->drho[i] = data_d_rhonu[i+1][1];
		paramrelic->Tnp[i] = data_T_n2p[i+1][1];
		paramrelic->anp[i] = data_a_n2p[i+1][1];
		paramrelic->bnp[i] = data_b_n2p[i+1][1];
		paramrelic->cnp[i] = data_c_n2p[i+1][1];
		paramrelic->dnp[i] = data_d_n2p[i+1][1];
		paramrelic->Tpn[i] = data_T_p2n[i+1][1];
		paramrelic->apn[i] = data_a_p2n[i+1][1];
		paramrelic->bpn[i] = data_b_p2n[i+1][1];
		paramrelic->cpn[i] = data_c_p2n[i+1][1];
		paramrelic->dpn[i] = data_d_p2n[i+1][1];
		//printf("%.15e, %.15e, %.15e, %.15e, %.15e \n", paramrelic.Tcm_rho[i], paramrelic.arho[i], paramrelic.brho[i], paramrelic.crho[i], paramrelic.drho[i]);
	}
	for (int i=0; i<row-2; i++)
	{
		paramrelic->TdQdt[i] = data_T_dQdt[i+1][1];
		paramrelic->adQdt[i] = data_a_dQdt[i+1][1];
		paramrelic->bdQdt[i] = data_b_dQdt[i+1][1];
		paramrelic->cdQdt[i] = data_c_dQdt[i+1][1];
		paramrelic->ddQdt[i] = data_d_dQdt[i+1][1];
		//printf("%.15e, %.15e, %.15e, %.15e, %.15e \n", paramrelic->TdQdt[i], paramrelic->adQdt[i], paramrelic->bdQdt[i], paramrelic->cdQdt[i], paramrelic->ddQdt[i]);
	}
    
    for (int i=40; i < 50; i++)
    {
        printf("%.15e, %.15e, %.15e \n", paramrelic->Tnp[i], paramrelic->dnp[i], n2p_vs(paramrelic->Tnp[i]/K_to_eV/1000., paramrelic));
    }
	return;
}
// END OF INIT_VS() !!!!


double rate1_vs(double ms, double mix)
{
	double num, den, Gamma;
	num = 9*pow(Gf,2.)*alphaem*pow(ms,5.)*pow(sin(mix),2.);
  den = 512*pow(pi,4.);
  Gamma = num/den;
  return Gamma;
}

double rate2_vs(double ms, double mix)
{
	double part1, part2, Gamma;
	part1 = pow(Gf,2.)*pow(f_pi,2.)/(16*pi);
	part2 = ms*(pow(ms,2.)-pow(mpi_neutral,2.))*pow(sin(mix),2.);
	Gamma = part1*part2;
	return Gamma;
}

double rate3_vs(double ms, double mix)
{
	double part1, part2, term, Gamma;
	part1 = pow(Gf,2.)*pow(f_pi,2.)/(16*pi);
	term = (pow(ms,2.) - pow(mpi_charged+me,2.))*(pow(ms,2.) - pow(mpi_charged-me,2.));
	part2 = ms * sqrt(term) * pow(sin(mix),2.);
	Gamma = part1*part2;
	return 2*Gamma;
}

double rate4_vs(double ms, double mix)
{
	double part1, part2, term, Gamma;
	part1 = pow(Gf,2.)*pow(f_pi,2.)/(16*pi);
	term = (pow(ms,2.) - pow(mpi_charged+mu,2.))*(pow(ms,2.) - pow(mpi_charged-mu,2.));
	part2 = ms * sqrt(term) * pow(sin(mix),2.);
	Gamma = part1*part2;
	return 2*Gamma;
}

void ts(double ms, double mix, struct relicparam* paramrelic)
{
	double Gam1=rate1_vs(ms, mix);
	double Gam2=0;
	double Gam3=0;
	double Gam4=0;
	if (ms>mpi_neutral) {Gam2 = rate2_vs(ms,mix);} //mass of a neutral pion
	if (ms>(mpi_charged+me)) {Gam3 = rate3_vs(ms,mix);} //mass of a charged pion plus the mass of an electron
  if (ms>(mpi_charged+mu)) {Gam4 = rate4_vs(ms,mix);} //mass of a charged pion plus the mass of a muon
	paramrelic->tau_vs = 1/(Gam1+Gam2+Gam3+Gam4)*1000; //??multiply by 1000 just once since tau is in units of MeV^{-1}
	paramrelic->Gamma_vs = (Gam1+Gam2+Gam3+Gam4)/1000; //??divide by 1000 just once since gamma is in units of MeV
}

double dQdt_vs(double T, struct relicparam* paramrelic)
{
	double Thold = T*1000; //T is passed in units of GeV, but our T array is in terms of MeV, so we need to use MeV to compare the two
	double dQdt = 0.;
	int index = 0;
	double a,b,c,d,T_cs;

	if (Thold<0.6635) return 0; //cubic spline fits poorly after this

  for (int i=0; i<paramrelic->row-1; i++) //want to find the index we'll be using to find dQdt
	{
		if (Thold<paramrelic->TdQdt[i]) index=i; //if x is going up, not down, as i increases, this will need to be changed to if (x>output[i,0]):
	}

  T_cs = paramrelic->TdQdt[index];
	a = paramrelic->adQdt[index]; b = paramrelic->bdQdt[index]; c = paramrelic->cdQdt[index]; d = paramrelic->ddQdt[index];
  dQdt = a*pow((Thold-T_cs),3.) + b*pow((Thold-T_cs),2.) + c*(Thold-T_cs) + d;
	dQdt = dQdt/pow(1000,5.); //turn into GeV^5
	return dQdt;
}

double n2p_vs(double T9, struct relicparam* paramrelic)
{
	double Thold = T9*K_to_eV*1000; // T9 is in Kelvin, so K_to_eV turns it into GeV (not eV as one may think), then multiplying by 1000 gives MeV
	double n2p;
	int index = 0;
	double a,b,c,d,T_cs;

	for (int i=0; i<paramrelic->row; i++) //want to find the index we'll be using to find the neutron to proton rate
	{
		if (Thold<paramrelic->Tnp[i]) index=i; //if x is going up, not down, as i increases, this will need to be changed to if (x>output[i,0]):
	}

	T_cs = paramrelic->Tnp[index];
	a = paramrelic->anp[index]; b = paramrelic->bnp[index]; c = paramrelic->cnp[index]; d = paramrelic->dnp[index];
    
    double log_Trat = log(Thold/T_cs);
    n2p = a * pow(log_Trat, 3.) + b*pow(log_Trat, 2.) + c*log_Trat + d;
//	n2p = a*pow((Thold-T_cs),3.) + b*pow((Thold-T_cs),2.) + c*(Thold-T_cs) + d;
	n2p = exp(n2p); //turn into some other unit??

	return n2p;
}

double p2n_vs(double T9, struct relicparam* paramrelic)
{
	double Thold = T9*K_to_eV*1000; //// T9 is in Kelvin, so K_to_eV turns it into GeV (not eV as one may think), then multiplying by 1000 gives MeV
	double p2n;
	int index = 0;
	double a,b,c,d,T_cs;

	for (int i=0; i<paramrelic->row; i++) //want to find the index we'll be using to find the proton to neutron rate
	{
		if (Thold<paramrelic->Tpn[i]) index=i; //if x is going up, not down, as i increases, this will need to be changed to if (x>output[i,0]):
	}

	T_cs = paramrelic->Tpn[index];
	a = paramrelic->apn[index]; b = paramrelic->bpn[index]; c = paramrelic->cpn[index]; d = paramrelic->dpn[index];
    
    double log_Trat = log(Thold/T_cs);
    p2n = a * pow(log_Trat, 3.) + b*pow(log_Trat, 2.) + c*log_Trat + d;

//	p2n = a*pow((Thold-T_cs),3.) + b*pow((Thold-T_cs),2.) + c*(Thold-T_cs) + d;
	p2n = exp(p2n); //turn into some other unit??

	return p2n;
}


/*--------------------------------------------------------------*/

double dark_density(double T, struct relicparam* paramrelic)
{
	if(paramrelic->phi_model) return 0.;
	if(paramrelic->vs_model) return 0.;

	if(paramrelic->size_table_rhoPD>1&&paramrelic->use_table_rhoPD)
	{
		int ie=1;

		if(T<paramrelic->table_rhoPD[0][paramrelic->size_table_rhoPD-1])
		{
			ie=paramrelic->size_table_rhoPD-1;
		}
		else if(T<paramrelic->table_rhoPD[0][0]&&T>paramrelic->table_rhoPD[0][paramrelic->size_table_rhoPD-1])
		{
			while(T<paramrelic->table_rhoPD[0][ie]&&ie<paramrelic->size_table_rhoPD) ie++;
		}

		double logrhoD1,logrhoD2,rhoD,logT1,logT2,logT;
		logT=log(T);
		logrhoD1=log(paramrelic->table_rhoPD[1][ie]);
		logrhoD2=log(paramrelic->table_rhoPD[1][ie-1]);
		logT1=log(paramrelic->table_rhoPD[0][ie]);
		logT2=log(paramrelic->table_rhoPD[0][ie-1]);

		rhoD=exp((logrhoD2-logrhoD1)/(logT2-logT1)*(logT-logT1)+logrhoD1);

		return rhoD;
	}

	if(T<paramrelic->Tdend) return 0.;


	if(paramrelic->energy_model==3)
	{
		double H0=67.8/3.0856e19; /* Hubble constant in second */
		double rho_Lambda=0.7*H0*H0/(8.*pi*Gn)/2.322e17;

		if(T<=paramrelic->quintT12) return rho_Lambda;

		double rho02=rho_Lambda;
		if(T<=paramrelic->quintT23) return rho02*pow(T/paramrelic->quintT12,paramrelic->quintn2);

		double rho03=rho02*pow(paramrelic->quintT23/paramrelic->quintT12,paramrelic->quintn2);
		if(T<=paramrelic->quintT34) return rho03*pow(T/paramrelic->quintT23,paramrelic->quintn3);

		double rho04=rho03*pow(paramrelic->quintT34/paramrelic->quintT23,paramrelic->quintn3);
		return rho04*pow(T/paramrelic->quintT34,paramrelic->quintn4);
	}

	if(paramrelic->energy_model==2)
	{
		if(paramrelic->Tddeq==0.) return 0.;

		double geffT=geff(T,paramrelic);
		double rhorad=pi*pi/30.*geffT*pow(T,4.);

		return rhorad*(geff(paramrelic->Tddeq,paramrelic)/geffT)*pow(heff(T,paramrelic)/heff(paramrelic->Tddeq,paramrelic),paramrelic->ndd/3.)*pow(T/paramrelic->Tddeq,paramrelic->ndd);
	}

	if(paramrelic->energy_model==1)
	{
		if(paramrelic->dd0==0.) return 0.;

		double rho_photon_1MeV=pi*pi/15.*1.e-12;

		return paramrelic->dd0*rho_photon_1MeV*pow(T/1.e-3,paramrelic->ndd);
	}

	return 0.;
}


double dark_density_pressure(double T, struct relicparam* paramrelic)
{
	if(paramrelic->phi_model) return 0.;
	if(paramrelic->vs_model) return 0.;

	if(T<paramrelic->Tdend) return 0.;

	if(paramrelic->energy_model==1) if(paramrelic->dd0==0.) return 0.;

	if(paramrelic->energy_model==2) if(paramrelic->Tddeq==0.) return 0.;

	double ddark_density_dT=(dark_density(T*1.001,paramrelic)-dark_density(T*0.999,paramrelic))/0.002/T;

	double dentropy_dT=2.*pi*pi/45.*(heff(T*1.001,paramrelic)*pow(T*1.001,3.)-heff(T*0.999,paramrelic)*pow(T*0.999,3.))/0.002/T+dark_entropy_derivative(T,paramrelic);

	double entropy=2.*pi*pi/45.*heff(T,paramrelic)*pow(T,3.)+dark_entropy(T,paramrelic);

	//return (paramrelic->ndd/3.-1.)*dark_density(T,paramrelic); /* outdated */
	return (ddark_density_dT-dentropy_dT/entropy*dark_density(T,paramrelic))*entropy/dentropy_dT;
}

/*--------------------------------------------------------------*/

double sigma_entropy(double T, struct relicparam* paramrelic)
{
	if(paramrelic->phi_model) return 1.;
	if(paramrelic->vs_model) return 1.;

	if(paramrelic->Sigmarad0==0.) return 1.;

	double lnT,dlnT,Ttmp;
	int ie,nmax;
	double integ=0.;

	double heffT,geffT,darkdensitytilde,heffTdT,geffTdT,darkdensitytildeTdT,Htilde;
	double Sigmatildestar,dSigmatildestar_dT;

	nmax=10;

	lnT=log(1.e-15);

	dlnT=(log(T)-lnT)/nmax;

	for(ie=1;ie<nmax;ie++)
	{
		lnT+=dlnT;
		Ttmp=exp(lnT);

		heffT=heff(Ttmp,paramrelic);
		geffT=geff(Ttmp,paramrelic);
		darkdensitytilde=dark_density(Ttmp,paramrelic)/(pi*pi/30.*geffT*pow(Ttmp,4.));

		Htilde=sqrt(1.+darkdensitytilde); /*Htilde = H / sqrt(8 pi / 3 M_P^2) / rho_rad */

		Sigmatildestar=45.*sqrt(5.)/4./pow(pi,3.5)*Mplanck/heffT/sqrt(geffT)/pow(Ttmp,5.)/Htilde*entropy_Sigmarad(Ttmp,paramrelic);

		heffTdT=heff(Ttmp*1.001,paramrelic);
		geffTdT=geff(Ttmp*1.001,paramrelic);
		darkdensitytildeTdT=dark_density(Ttmp*1.001,paramrelic)/(pi*pi/30.*geffTdT*pow(Ttmp*1.001,4.));

		dSigmatildestar_dT=((45.*sqrt(5.)/4./pow(pi,3.5)*Mplanck/heffTdT/sqrt(geffTdT)/pow(Ttmp*1.001,5.)/sqrt(1.+darkdensitytildeTdT)*entropy_Sigmarad(Ttmp*1.001,paramrelic))-Sigmatildestar)/0.001/Ttmp;

		integ+=-dSigmatildestar_dT/pow(1.-Sigmatildestar,2.)*log(heffT*pow(Ttmp,3.));
	}

	heffT=heff(T,paramrelic);
	geffT=geff(T,paramrelic);
	darkdensitytilde=dark_density(T,paramrelic)/(pi*pi/30.*geffT*pow(T,4.));

	Htilde=sqrt(1.+darkdensitytilde); /*Htilde = H / sqrt(8 pi / 3 M_P^2) / rho_rad */

	Sigmatildestar=45.*sqrt(5.)/4./pow(pi,3.5)*Mplanck/heffT/sqrt(geffT)/pow(T,5.)/Htilde*entropy_Sigmarad(T,paramrelic);

	heffTdT=heff(T*1.001,paramrelic);
	geffTdT=geff(T*1.001,paramrelic);
	darkdensitytildeTdT=dark_density(T*1.001,paramrelic)/(pi*pi/30.*geffTdT*pow(T*1.001,4.));

	dSigmatildestar_dT=((45.*sqrt(5.)/4./pow(pi,3.5)*Mplanck/heffTdT/sqrt(geffTdT)/pow(T*1.001,5.)/sqrt(1.+darkdensitytildeTdT)*entropy_Sigmarad(T*1.001,paramrelic))-Sigmatildestar)/0.001/T;

	integ+=-dSigmatildestar_dT/pow(1.-Sigmatildestar,2.)*log(heffT*pow(T,3.))/2.;

	integ*=dlnT;

	return exp(integ);
}


double dark_entropy(double T, struct relicparam* paramrelic)
{
	if(paramrelic->phi_model) return 0.;
	if(paramrelic->vs_model) return 0.;

	if((paramrelic->sd0==0.)&&(paramrelic->Sigmad0==0.)) return 0.;

	if((paramrelic->Sigmad0==0.)&&(T<paramrelic->Tsend)) return 0.;

	if(paramrelic->Sigmad0==0.)
	{
		double s_photon_1MeV=4.*pi*pi/45.*1.e-9;

		return paramrelic->sd0*s_photon_1MeV*pow(T/1.e-3,paramrelic->nsd);
	}
	else
	{
		double lnT,dlnT,Ttmp;
		int ie,nmax;
		double integ=0.;

		double heffT,geffT,darkdensitytilde,Htilde;
		double Sigmatildestar=0.;

		nmax=50;

		lnT=log(1.e-15);

		dlnT=(log(T)-lnT)/nmax;

		for(ie=1;ie<nmax;ie++)
		{
			lnT+=dlnT;
			Ttmp=exp(lnT);

			heffT=heff(Ttmp,paramrelic);
			geffT=geff(Ttmp,paramrelic);
			darkdensitytilde=dark_density(Ttmp,paramrelic)/(pi*pi/30.*geffT*pow(Ttmp,4.));

			Htilde=sqrt(1.+darkdensitytilde); /*Htilde = H / sqrt(8 pi / 3 M_P^2) / rho_rad */

			if(paramrelic->Sigmarad0!=0.) Sigmatildestar=45.*sqrt(5.)/4./pow(pi,3.5)*Mplanck/heffT/sqrt(geffT)/pow(Ttmp,5.)/Htilde*entropy_Sigmarad(Ttmp,paramrelic);

			integ+=sgStar(Ttmp,paramrelic)*dark_entropy_Sigmad(Ttmp,paramrelic)/Htilde/(1.-Sigmatildestar)/sigma_entropy(Ttmp,paramrelic)/pow(heffT*pow(Ttmp,3.),(2.-Sigmatildestar)/(1.-Sigmatildestar));
		}

		heffT=heff(T,paramrelic);
		geffT=geff(T,paramrelic);
		darkdensitytilde=dark_density(T,paramrelic)/(pi*pi/30.*geffT*pow(T,4.));

		Htilde=sqrt(1.+darkdensitytilde); /*Htilde = H / sqrt(8 pi / 3 M_P^2) / rho_rad */

		if(paramrelic->Sigmarad0!=0.) Sigmatildestar=45.*sqrt(5.)/4./pow(pi,3.5)*Mplanck/heffT/sqrt(geffT)/pow(T,5.)/Htilde*entropy_Sigmarad(T,paramrelic);

		integ+=sgStar(T,paramrelic)*dark_entropy_Sigmad(T,paramrelic)/Htilde/(1.-Sigmatildestar)/sigma_entropy(T,paramrelic)/pow(heffT*pow(T,3.),(2.-Sigmatildestar)/(1.-Sigmatildestar))/2.;

		integ*=dlnT;

		return Mplanck*sqrt(45./4./pow(pi,3.))*pow(heffT*pow(T,3.),1./(1.-Sigmatildestar))*sigma_entropy(T,paramrelic)*integ;
	}
}


double dark_entropy_derivative(double T, struct relicparam* paramrelic)
{
	if(paramrelic->phi_model) return 0.;
	if(paramrelic->vs_model) return 0.;

	if((paramrelic->sd0==0.)&&(paramrelic->Sigmad0==0.)) return 0.;

	if((paramrelic->Sigmad0==0.)&&(T<paramrelic->Tsend)) return 0.;

	if(paramrelic->Sigmad0==0.)
	{
		return paramrelic->nsd*dark_entropy(T,paramrelic)/T;
	}
	else
	{
		double heffT=heff(T,paramrelic);
		double geffT=geff(T,paramrelic);
		double sradT=2.*pi*pi/45.*heffT*pow(T,3.);
		double rhoradT=pi*pi/30.*geffT*pow(T,4.);
		double darkdensitytilde=dark_density(T,paramrelic)/rhoradT;
		double Sigmarad=entropy_Sigmarad(T,paramrelic);

		double Htilde=sqrt(1.+darkdensitytilde); /*Htilde = H / sqrt(8 pi / 3 M_P^2) / rho_rad */

		double Sigmatildestar=45.*sqrt(5.)/4./pow(pi,3.5)*Mplanck/heffT/sqrt(geffT)/pow(T,5.)/Htilde*Sigmarad;

		return 3.*sgStar(T,paramrelic)/T/(1.-Sigmatildestar)/heffT
		*(sqrt(geffT)*dark_entropy(T,paramrelic)
		-sqrt(5.*Mplanck/4./pow(pi,3.))/T/T*dark_entropy_Sigmad(T,paramrelic)/Htilde);
	}
}


double dark_entropy_Sigmad(double T, struct relicparam* paramrelic)
{
	if(paramrelic->phi_model) return 0.;
	if(paramrelic->vs_model) return 0.;

	if((paramrelic->sd0==0.)&&(paramrelic->Sigmad0==0.)) return 0.;

	if((paramrelic->Sigmad0==0.)&&(T<paramrelic->Tsend)) return 0.;

	if(paramrelic->Sigmad0==0.)
	{
		double heffT=heff(T,paramrelic);
		double geffT=geff(T,paramrelic);
		double sradT=2.*pi*pi/45.*heffT*pow(T,3.);
		double rhoradT=pi*pi/30.*geffT*pow(T,4.);
		double darkdensitytilde=dark_density(T,paramrelic)/rhoradT;
		double Sigmarad=entropy_Sigmarad(T,paramrelic);

		double Htilde=sqrt(1.+darkdensitytilde); /*Htilde = H / sqrt(8 pi / 3 M_P^2) / rho_rad */

		double Sigmatildestar=45.*sqrt(5.)/4./pow(pi,3.5)*Mplanck/heffT/sqrt(geffT)/pow(T,5.)/Htilde*Sigmarad;

		return sqrt(4.*pow(pi,3.)/5.)/Mplanck*Htilde*T*T*(sqrt(geffT)*dark_entropy(T,paramrelic)-heffT/3./sgStar(T,paramrelic)*T*(1.-Sigmatildestar)*dark_entropy_derivative(T,paramrelic));
	}
	else
	{
		if(T<paramrelic->TSigmadend) return 0.;

		double s_photon_1MeV=4.*pi*pi/45.*1.e-9;

		double Sigma_photon_1MeV= 1./Mplanck*sqrt(8.*pi*pi*pi/5.)*(1.e-6)*s_photon_1MeV;

		return paramrelic->Sigmad0*Sigma_photon_1MeV*pow(T/1.e-3,paramrelic->nSigmad);
	}
}

/*--------------------------------------------------------------*/

double entropy_Sigmarad(double T, struct relicparam* paramrelic)
{
	if(paramrelic->phi_model) return 0.; //paramrelic->Gamma_phi*paramrelic->rho_phi/T;
	if(paramrelic->vs_model) return 0.; //dQdt_vs(T,paramrelic)/T;

	if((paramrelic->Sigmarad0==0.)||(T<paramrelic->TSigmaradend)) return 0.;

	double s_photon_1MeV=4.*pi*pi/45.*1.e-9;

	double Sigma_photon_1MeV=1./Mplanck*sqrt(8.*pi*pi*pi/5.)*(1.e-6)*s_photon_1MeV;

	return paramrelic->Sigmarad0*Sigma_photon_1MeV*pow(T/1.e-3,paramrelic->nSigmarad);
}

/*--------------------------------------------------------------*/

double nonthermal(double T, struct relicparam* paramrelic)
{
	if(paramrelic->phi_model) return 0.; //paramrelic->eta_phi*paramrelic->Gamma_phi*paramrelic->rho_phi;
	if(paramrelic->vs_model) return 0.;

	if((paramrelic->nt0==0.)||(T<paramrelic->Tnend)) return 0.;

	return paramrelic->nt0*1.e-50*pow(T/1.e-3,paramrelic->nnt);
}

/*--------------------------------------------------------------*/


double neutdens_vs(double Tnu, struct relicparam* paramrelic)
{
	double Thold = Tnu*1000; //Tnu is passed in units of GeV, but our Tcm array is in terms of MeV, so we need to use MeV to compare the two
	double rho = 0.;
	int index = 0;
	double a,b,c,d,Tcm;

  for (int i=0; i<paramrelic->row; i++) //want to find the index we'll be using to find rho
	{
		if (Thold<paramrelic->Tcm_rho[i]) {index=i;} //if x is going up, not down, as i increases, this will need to be changed to if (x>output[i,0]):
	}

  Tcm = paramrelic->Tcm_rho[index];
	a = paramrelic->arho[index]; b = paramrelic->brho[index]; c = paramrelic->crho[index]; d = paramrelic->drho[index];

  rho = a*pow((Thold-Tcm),3.) + b*pow((Thold-Tcm),2.) + c*(Thold-Tcm) + d;
	rho = rho/pow(1000,4.); //To turn it back to Gev? hmm, but rho goes like Tcm^4 so honestly maybe need to divide by 1000^4? oof
	return rho;
}


double neutdens(double Tnu, struct relicparam* paramrelic)
/* Computes the neutrino density, including any effects from a neutrino degeneracy */
{
     if((paramrelic->xinu1==0.)&&(paramrelic->xinu2==0.)&&(paramrelic->xinu3==0.))
    {
        // No degeneracy, relativistic approximation
        return 2.*pi*pi/30.*7./8.*paramrelic->Nnu*pow(Tnu,4.);
    }

    int ie,je,n;
    double rho=0.;
    double xinu[4];
    double max1,max2,int1,int2;
    double x;

    xinu[1]=paramrelic->xinu1;
    xinu[2]=paramrelic->xinu2;
    xinu[3]=paramrelic->xinu3;

    // SM neutrinos
    for(ie=1;ie<=3;ie++)
    {
        // The factor (paramrelic->Nnu/3.) includes extra DOF from non-rel. e+- and non-inst. nu decoupling
        if(fabs(xinu[ie])<=0.03)
        {
            rho+=(paramrelic->Nnu/3.)*2.*pi*pi/30.*pow(Tnu,4.)*(7./8.+(15./(4*pi*pi))*xinu[ie]*xinu[ie]+
                                                               (15./(8.*pow(pi,4.)))*pow(xinu[ie],4.));
        }
        else if(fabs(xinu[ie])>=30.)
        {
            rho+=(paramrelic->Nnu/3.)*pow(Tnu,4.)/(8.*pi*pi)*pow(xinu[ie],4.)*(1.+12.*1.645/xinu[ie]/xinu[ie]);
        }
        else
        {
            // Neutrinos
            max1=(88.029+xinu[ie])*Tnu;
            int1=0.;
            n=50;
            for(je=1;je<=n-1;je++)
            {
                x=(double)je/(double)n*max1;
                int1+=1./(2.*pi*pi)*pow(x,3.)/(1.+exp(x/Tnu-xinu[ie]));
            }
            int1+=0.5*1./(2.*pi*pi)*pow(max1,3.)/(1.+exp(max1/Tnu-xinu[ie]));
            int1*=(paramrelic->Nnu/3.)*max1/(double)n;
            rho+=int1;

            // Anti-neutrinos
            max2=(88.029-xinu[ie])*Tnu;
            if(max2>0.)
            {
                int2=0.;
                n=50;
                for(je=1;je<=n-1;je++)
                {
                    x=(double)je/(double)n*max2;
                    int2+=1./(2.*pi*pi)*pow(x,3.)/(1.+exp(x/Tnu+xinu[ie]));
                }
                int2+=0.5/(2.*pi*pi)*pow(max2,3.)/(1.+exp(max2/Tnu+xinu[ie]));
                int2*=(paramrelic->Nnu/3.)*max2/(double)n;
                rho+=int2;
            }
        }
    }

    return rho;
}

/*--------------------------------------------------------------*/


double neutdens_deriv_vs(double Tnu, struct relicparam* paramrelic) //is this in the right units? I think it needs to be in GeV?
{
	double Thold = Tnu*1000; //Tnu is passed in units of GeV
	double drho = 0.;
	int index = 0;
	double a,b,c,Tcm;

  for (int i=0; i<paramrelic->row; i++) //want to find the index we'll be using to find rho
	{
      if (Thold<paramrelic->Tcm_rho[i]) {index=i;} //if x is going up, not down, as i increases, this will need to be changed to if (x>output[i,0]):
	}

  Tcm = paramrelic->Tcm_rho[index];
	a = paramrelic->arho[index]; b = paramrelic->brho[index]; c = paramrelic->crho[index];

  drho = 3*a*pow((Thold-Tcm),2.) + 2*b*(Thold-Tcm) + c;
	drho = drho/pow(1000,3.); //To turn it back to Gev? whyyyy does it calculate everything in GeV
	return drho;
}


double neutdens_deriv(double Tnu, struct relicparam* paramrelic)
/* Computes the temperature (Tnu) derivative of the neutrino energy density */
{
    if((paramrelic->xinu1==0.)&&(paramrelic->xinu2==0.)&&(paramrelic->xinu3==0.))
    {
        return 7.*pi*pi/30.*paramrelic->Nnu*pow(Tnu,3.);
    }

    int ie,je,n;
    double drho=0.;
    double xinu[4];
    double max1,max2,int1,int2;
    double x;

    xinu[1]=paramrelic->xinu1;
    xinu[2]=paramrelic->xinu2;
    xinu[3]=paramrelic->xinu3;

    /* SM neutrinos */
    for(ie=1;ie<=3;ie++)
    {
        if(fabs(xinu[ie])<=0.03)
        {
            drho+=(paramrelic->Nnu/3.)*4.*pi*pi/15.*pow(Tnu,3.)*(7./8.+(15./(4*pi*pi))*xinu[ie]*xinu[ie]
                                                                +(15./(8.*pow(pi,4.)))*pow(xinu[ie],4.));
        }
        else if(fabs(xinu[ie])>=30.)
        {
            drho+=(paramrelic->Nnu/3.)*pow(Tnu,3.)/(2.*pi*pi)*pow(xinu[ie],4.)*(1.+12.*1.645/xinu[ie]/xinu[ie]);
        }
        else
        {
            max1=(88.029+xinu[ie])*Tnu;
            int1=0.;
            n=50;
            for(je=1;je<=n-1;je++)
            {
                x=(double)je/(double)n*max1;
                int1+=1./(2.*pi*pi)*pow(x,3.)/(1.+exp(x/Tnu-xinu[ie]));
            }
            int1+=0.5*1./(2.*pi*pi)*pow(max1,3.)/(1.+exp(max1/Tnu-xinu[ie]));
            int1*=(paramrelic->Nnu/3.)*4.*max1/Tnu/(double)n;
            drho+=int1;

            max2=(88.029-xinu[ie])*Tnu;
            if(max2>0.)
            {
                int2=0.;
                n=50;
                for(je=1;je<=n-1;je++)
                {
                    x=(double)je/(double)n*max2;
                    int2+=1./(2.*pi*pi)*pow(x,3.)/(1.+exp(x/Tnu+xinu[ie]));
                }
                int2+=0.5/(2.*pi*pi)*pow(max2,3.)/(1.+exp(max2/Tnu+xinu[ie]));
                int2*=(paramrelic->Nnu/3.)*4.*max2/Tnu/(double)n;
                drho+=int2;
            }
        }
    }
    return drho;
}

/*--------------------------------------------------------------*/

double neutN(double T)
/* Computes the round N(z) function of the neutrinos - Pisanti et al., 0705.0290, eq. (A24) */
{
	double z=m_e/T;

	if(z>=4.) return 0.;

	double logNz=0.;

	double n[14]={-10.21703221236002,61.24438067531452,-340.3323864212157,1057.2707914654834,-2045.577491331372,2605.9087171012848,-2266.1521815470196,1374.2623075963388,-586.0618273295763,174.87532902234145,-35.715878215468045,4.7538967685808755,-0.3713438862054167,0.012908416591272199};

	int ie;
	for(ie=0;ie<14;ie++) logNz+=n[ie]*pow(z,ie);

	return exp(logNz);
}
