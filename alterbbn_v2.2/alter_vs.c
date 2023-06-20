#include "src/include.h"

/* ----------------------------------------------------- */
/* Calculation of the abundance of the elements from BBN */
/* ----------------------------------------------------- */

// For std cosmological model: ./alter_vs.x 14.87365053664758 6.1e-10 1 180 0 0 0
// For mass = 300 MeV mixing angle = 7e-05 model: ./alter_vs.x 14.87365053664758 6.1e-10 1 255 104.56255360706342 300 7e-05

int main(int argc, char** argv)
{
   /*printf("here\n");
    
   if (argc==8)
	{
        printf("here 1\n");
    }  */ 
    
	struct relicparam paramrelic;
	double ratioH[NNUC+1],cov_ratioH[NNUC+1][NNUC+1];
	double H2_H, He3_H, Yp, Li7_H, Li6_H, Be7_H;
	double sigma_H2_H, sigma_He3_H, sigma_Yp,sigma_Li7_H, sigma_Li6_H, sigma_Be7_H;
	double eta, Tinit;
	int row;
	char ms_ch[256];
	char mix_ch[256];
	double ms_d;
	double mix_d;
	double ns0;

	int failsafe;

	if (argc<8)
	{
		printf(" This program needs 7 parameters:\n"
		"   Tinit       initial temperature in MeV (2.3267 MeV by default)\n"
		"   eta0        initial value of the baryon-to-photon ratio (default: 6.1e-10)\n"
		"   failsafe    0=fast, 1=precise, 6=robust but slow. See stand_cosmo.c for more options.\n"
		"   row         number of rows in the rho_nu file, including the header row.\n"
		"   ns0         initial number density of sterile neutrinos in units of MeV^3.\n"
		"   mass        The mass of the sterile neutrino model of interest.\n"
	  "   mix         The mixing angle of the sterile neutrino model of interset.\n");
        //printf("hi");
		exit(1);
	}
	else
	{
		sscanf(argv[1],"%lf",&Tinit); //default is Tinit=27.*K_to_eV*1.e3
		sscanf(argv[2],"%lf",&eta); //default is eta=6.1e-10;
		sscanf(argv[3],"%d",&failsafe); //default is failsafe=1;
		sscanf(argv[4],"%d",&row); //number of rows in rho_nu and drho_nu
		sscanf(argv[5],"%lf",&ns0);
		sscanf(argv[6],"%s",ms_ch);
		sscanf(argv[7],"%s",mix_ch);
		sscanf(argv[6],"%lf",&ms_d);
		sscanf(argv[7],"%lf",&mix_d);
	}
    
    //printf("here 2\n");
    
	Init_cosmomodel(&paramrelic); //this sets the values of paramrelic that *every* model will use, other init functions set specific values that only that model uses

	paramrelic.failsafe = failsafe;
	paramrelic.row = row;

	Init_cosmomodel_param(eta,paramrelic.Nnu,0.,paramrelic.life_neutron,paramrelic.life_neutron_error,0.,0.,0.,&paramrelic);
	// why are Nnu, life_neutron, and life_neutron_error defined above by Init_cosmomodel_param when alter_eta supposedly just... alters eta?
	paramrelic.Tinit=Tinit*1.e-3/K_to_eV; //setting Tinit, why isn't this done in Init_cosmomodel_param?
    
    //printf("here 4\n");
    
	Init_vs(ms_ch, mix_ch, ms_d, mix_d, ns0, row, &paramrelic);
    
    //printf("here 3 !!!!! \n");
    
	printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
	//paramrelic.err=2;
	//nucl(&paramrelic,ratioH);
	//H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
	//printf("low:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);

	paramrelic.err=0;
	nucl(&paramrelic,ratioH);
	H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
	printf("cent:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);

	//paramrelic.err=1;
	//nucl(&paramrelic,ratioH);
	//H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
	//printf("high:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);


	/*paramrelic.err=3;
	if(nucl_err(&paramrelic,ratioH,cov_ratioH)==1)
	{
		printf("--------------------\n\n");
		printf("With uncertainties:\n");
		H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
		sigma_H2_H=sqrt(cov_ratioH[3][3]);sigma_Yp=sqrt(cov_ratioH[6][6]);sigma_Li7_H=sqrt(cov_ratioH[8][8]);sigma_Be7_H=sqrt(cov_ratioH[9][9]);sigma_He3_H=sqrt(cov_ratioH[5][5]);sigma_Li6_H=sqrt(cov_ratioH[7][7]);
		printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");

		printf("value:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);
		printf(" +/- :\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n\n",sigma_Yp,sigma_H2_H,sigma_He3_H,sigma_Li7_H,sigma_Li6_H,sigma_Be7_H);

		double corr_ratioH[NNUC+1][NNUC+1];
		for(int ie=1;ie<=NNUC;ie++) for(int je=1;je<=NNUC;je++) corr_ratioH[ie][je]=cov_ratioH[ie][je]/sqrt(cov_ratioH[ie][ie]*cov_ratioH[je][je]);
		printf("Correlation matrix:\n");
		printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
		printf("Yp\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[6][6],corr_ratioH[6][3],corr_ratioH[6][5],corr_ratioH[6][8],corr_ratioH[6][7],corr_ratioH[6][9]);
		printf("H2/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[3][6],corr_ratioH[3][3],corr_ratioH[3][5],corr_ratioH[3][8],corr_ratioH[3][7],corr_ratioH[3][9]);
		printf("He3/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[5][6],corr_ratioH[5][3],corr_ratioH[5][5],corr_ratioH[5][8],corr_ratioH[5][7],corr_ratioH[5][9]);
		printf("Li7/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[8][6],corr_ratioH[8][3],corr_ratioH[8][5],corr_ratioH[8][8],corr_ratioH[8][7],corr_ratioH[8][9]);
		printf("Li6/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[7][6],corr_ratioH[7][3],corr_ratioH[7][5],corr_ratioH[7][8],corr_ratioH[7][7],corr_ratioH[7][9]);
		printf("Be7/H\t %f\t %f\t %f\t %f\t %f\t %f\n\n",corr_ratioH[9][6],corr_ratioH[9][3],corr_ratioH[9][5],corr_ratioH[9][8],corr_ratioH[9][7],corr_ratioH[9][9]);
	}
	else printf("Uncertainty calculation failed\n\n");
	*/

	/*paramrelic.err=4;
	if(nucl_err(&paramrelic,ratioH,cov_ratioH))
	{
		printf("--------------------\n\n");
		printf("With MC uncertainties:\n");
				H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
		sigma_H2_H=sqrt(cov_ratioH[3][3]);sigma_Yp=sqrt(cov_ratioH[6][6]);sigma_Li7_H=sqrt(cov_ratioH[8][8]);sigma_Be7_H=sqrt(cov_ratioH[9][9]);sigma_He3_H=sqrt(cov_ratioH[5][5]);sigma_Li6_H=sqrt(cov_ratioH[7][7]);
		printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");

		printf("mean:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);
		printf(" +/- :\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n\n",sigma_Yp,sigma_H2_H,sigma_He3_H,sigma_Li7_H,sigma_Li6_H,sigma_Be7_H);

		double corr_ratioH[NNUC+1][NNUC+1];
		for(int ie=1;ie<=NNUC;ie++) for(int je=1;je<=NNUC;je++) corr_ratioH[ie][je]=cov_ratioH[ie][je]/sqrt(cov_ratioH[ie][ie]*cov_ratioH[je][je]);
		printf("Correlation matrix:\n");
		printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
		printf("Yp\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[6][6],corr_ratioH[6][3],corr_ratioH[6][5],corr_ratioH[6][8],corr_ratioH[6][7],corr_ratioH[6][9]);
		printf("H2/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[3][6],corr_ratioH[3][3],corr_ratioH[3][5],corr_ratioH[3][8],corr_ratioH[3][7],corr_ratioH[3][9]);
		printf("He3/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[5][6],corr_ratioH[5][3],corr_ratioH[5][5],corr_ratioH[5][8],corr_ratioH[5][7],corr_ratioH[5][9]);
		printf("Li7/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[8][6],corr_ratioH[8][3],corr_ratioH[8][5],corr_ratioH[8][8],corr_ratioH[8][7],corr_ratioH[8][9]);
		printf("Li6/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[7][6],corr_ratioH[7][3],corr_ratioH[7][5],corr_ratioH[7][8],corr_ratioH[7][7],corr_ratioH[7][9]);
		printf("Be7/H\t %f\t %f\t %f\t %f\t %f\t %f\n\n",corr_ratioH[9][6],corr_ratioH[9][3],corr_ratioH[9][5],corr_ratioH[9][8],corr_ratioH[9][7],corr_ratioH[9][9]);
	}
	else printf("Uncertainty calculation failed\n\n");*/

	/* paramrelic.err=0;
	int compat=bbn_excluded(&paramrelic);

	if(compat==1) printf("Excluded by BBN constraints (chi2 without correlations)\n");
	else if(compat==0) printf("Compatible with BBN constraints (chi2 without correlations)\n");
	else printf("Computation failed (chi2 without correlations)\n");

	paramrelic.err=3;
	compat=bbn_excluded(&paramrelic);

	if(compat==1) printf("Excluded by BBN constraints (chi2 including correlations)\n");
	else if(compat==0) printf("Compatible with BBN constraints (chi2 including correlations)\n");
	else printf("Computation failed (chi2 including correlations)\n"); */

	return 1;
}
