#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define X 0
#define Y 1
#define Z 2

double a(int i, int j, int k, int d, double r[][d], double k_constant);

int main()
{
    //char path[200] = "/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS/metcompc/Aula9_0310/Results/ManySprings.txt";    
	//n -> particles
	//d -> dimensions
    int t, i, j, k;
	const int n = 2, d = 2;
	double dt = 0.1, tmax = 200;
	double k_constant = 1., b = 1.;
	double m[n];
	double vh[n][d];
	double v[n][d];
	double r[n][d];
	for(i = 0; i < n; i++)
		for(k = 0; k < d; k++)
		{
			vh[i][k] = 0.;
			v[i][k] = 0.;
			r[i][k] = 0.;
		}
	//particle 0 starts at (0,b)
	//particle 1 starts at (10,0)
	r[0][Y] = b;
	r[1][X] = 10.;
	v[0][X] = 1.;
    
    FILE *file = fopen("data.txt", "w");
    for(t = 0; t < tmax; t++)
    {
		//Report
        for(i = 0; i < n; i++) 						//for each particle
			for(k = 0; k < d; k++)					//for each dimension
	        	fprintf(file, "%.2f ", r[i][k]);
		fprintf(file, "\n");
        
        //Velocity-verlet
		//TODO: generalize this loop...
		for(i = 0; i < n; i++) 						//for each particle
			for(j = 0; j < n; j++) 					//for each other particle
				if(j != i)
					for(k = 0; k < d; k++)			//for each dimension
						vh[i][k] = v[i][k] + a(i, j, k, d, r, k_constant)*dt/2.;

        for(i = 0; i < n; i++) 						//for each particle
			for(k = 0; k < d; k++)					//for each dimension
	            r[i][k] = r[i][k] + vh[i][k]*dt;

		//TODO: ...and this loop into one?
		for(i = 0; i < n; i++) 						//for each particle
			for(j = 0; j < n; j++) 					//for each other particle
				if(j != i)
					for(k = 0; k < d; k++)			//for each dimension
						v[i][k] = vh[i][k] + a(i, j, k, d, r, k_constant)*dt/2.;
    }
    fclose(file);
    return 0;
}

double a(int i, int j, int k, int d, double r[][d], double k_constant)
{
	int l;
	double m = 0.;
	for(l = 0; l < d; l++)
		m += pow(r[j][l] - r[i][l],2.);
	return -k_constant*(r[j][k]-r[i][k])/pow(m,3./2.);
}



