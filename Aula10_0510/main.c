#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define X 0
#define Y 1
#define Z 2

double a_x(double *x, double *y, int i, int j, double kk);
double a_y(double *x, double *y, int i, int j, double kk);

int main()
{
    //char path[200] = "/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS/metcompc/Aula9_0310/Results/ManySprings.txt";    
	//n -> particles
	//d -> dimensions
    int t, i, j, k, n = 2, d = 2;
	double dt = 0.1, tmax = 200;
	double kk = 1., b = 1.;
	double m[n];
	double vh[n][d];
	double v[n][d];
	double r[n][d];
	for(i = 0; i < n; i++)
		for(k = 0; k < d; i++)
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
		for(i = 0; i < n; i++)
        	fprintf(file, "%.2f %.2f ", x[i], y[i]);
		fprintf(file, "\n");
        
        //Velocity-verlet, TODO: generalize such loops
		for(i = 0; i < n; i++) //for each particle
			for(j = 0; j < n; j++) //for each other particle
				if(j != i)
					for(k = 0; k < d; k++)	//for each dimension
						vh[i][k] = v[i][k] + a_x(x, y, i, j, kk)*dt/2.;

        for(i = 0; i < n; i++)
			for(k = 0; k < d; k++)	//for each dimension
	            r[i][k] = r[i][k] + vh[i][k]*dt;

		for(i = 0; i < n; i++) //for each particle
			for(j = 0; j < n; j++) //for each other particle
				if(j != i)
					for(k = 0; k < d; k++)	//for each dimension
						v[i][k] = vh[i][k] + a_x(x, y, i, j, kk)*dt/2.;
    }
    fclose(file);
    return 0;
}

double a_x(double *x, double *y, int i, int j, double kk)
{
	double dx = x[j]-x[i];
	double dy = y[j]-y[i];
	double m = pow(pow(dx,2.)+pow(dy,2.),3./2.);
	return -k*dx/m;
}

double a_y(double *x, double *y, int i, int j, double kk)
{
	double dx = x[j]-x[i];
	double dy = y[j]-y[i];
	double m = pow(pow(dx,2.)+pow(dy,2.),3./2.);
	return -k*dy/m;
}




