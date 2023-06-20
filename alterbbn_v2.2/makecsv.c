#include<stdio.h>
#include<string.h>


void enter_line_csv(FILE *fp,double a[6],int m)
{
	//fopen(fp);
        for(int i=0;i<m;i++)
	{
		fprintf(fp,",%e ",a[i]);
	}
	//fclose(fp);
}

void create_csv(char *filename,double a[][6],int m,int n)
{
	printf("Creating %s.csv file\n",filename);
	FILE *fp;
	int i,j;
	filename=strcat(filename,".csv");
	fp=fopen(filename,"w+");
	fprintf(fp,"bound, Yp, H2/H, He3/H, Li7/H, Li6/H, He7/H");

	for(i=0;i<m;i++)
	{
		if (i==0) fprintf(fp,"\n%s","low");
		else if (i==1) fprintf(fp, "\n%s", "cent");
		else fprintf(fp, "\n%s", "high");
    		for(j=0;j<n;j++)
		{
        		fprintf(fp,",%e ",a[i][j]);

		}
	}

	fclose(fp);
}

int main()
{
    	//double a[3][6]={{2.474e-01, 2.526e-05, 1.025e-05, 5.028e-10, 1.689e-15, 4.745e-10},
	//	{2.473e-01, 2.463e-05, 1.034e-05, 5.376e-10, 1.085e-14, 5.087e-10},
	//	{2.473e-01, 2.404e-05, 1.044e-05, 5.746e-10, 3.522e-14, 5.454e-10}};
	//char str[100] = "trythis";
	//create_csv(str,a,3,6);

	char str2[100] = "trythis2";
	printf("Creating %s.csv file\n",str2);
	FILE *fp;
	char *filename=strcat(str2,".csv");
	fp=fopen(filename,"w+");
	fprintf(fp,"Temperature, Yp, H2/H, He3/H, Li7/H, Li6/H, He7/H");
	double abundance[7] = {3,1,2,3,4,5,6};

	for (int i=0; i<100; i++)
	{
		//enter_line_csv(*fp,abundance,6);
		fprintf(fp, "\n");
		for (int j=0; j<7; j++)
		{
			fprintf(fp,"%e,", abundance[j]);
		}
		abundance[0]-=0.01;
		for (int j=1; j<7; j++)
		{
			abundance[j]++;
		}
	}
	fclose(fp);
	return 0;
}
