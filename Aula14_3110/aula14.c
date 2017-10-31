#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#define X 0
#define Y 1
#define Z 2

double a(int i, int j, int kk, int d, double r[][d], double sides[d], double k_constant, double m);
double distance(int i, int j, int k, int d, double r[][d], double sides[d]);

int main()
{
    srand((unsigned int)time(NULL));
    int t, i, j, k;
    const int n = 2, d = 2;
    double dt = 0.01, tmax = 2000, T = 0.1, m = 1., k_constant = 1.;
    double v[n][d], vh[n][d], r[n][d], sides[d], spacing[d];
    sides[0] = 100.;
    sides[1] = 100.;
	for(k = 0; k < d; k++)        
		spacing[k] = sides[k]/sqrt(n);
    
    //default values
    for(i = 0; i < n; i++)
        for(k = 0; k < d; k++)
        {
            v[i][k] = 0.;
            r[i][k] = 0.;
        }

	//initial positions
	r[0][X] = 40;
	r[0][Y] = 50;
	r[1][X] = 60;
	r[1][Y] = 50;

	//initial velocities
	v[0][X] = 10;
	v[0][Y] = 0;
	v[1][X] = -10;
	v[1][Y] = 0;
 
    for(t = 0; t < tmax; t++)
    {
		printf("set size square\n");
		printf("set xrange [0:%d]\n", (int)sides[0]);
		printf("set yrange [0:%d]\n", (int)sides[1]);
		printf("plot \"-\" w p pt 7 ps 2\n");

        for(i = 0; i < n; i++){
            for(k = 0; k < d; k++)
				printf("%.2f ", r[i][k]);
			printf("\n");
		}
        printf("e\n");

         //Velocity-verlet
        for(i = 0; i < n; i++)                                                      //for each particle
            for(j = 0; j < n; j++)                                                  //for each other particle
                if(j != i)                                                          //for all particles other than particle i
                    for(k = 0; k < d; k++)                                          //for each dimension
                        vh[i][k] = v[i][k] + a(i, j, k, d, r, sides, k_constant, m)*dt/2.;
        
        for(i = 0; i < n; i++)                                                      //for each particle
            for(k = 0; k < d; k++) {                                                //for each dimension
                r[i][k] = r[i][k] + vh[i][k]*dt;
                if(r[i][k] >= sides[k])                                             //PCB
                    r[i][k] = r[i][k]-sides[k];
                if(r[i][k] <= 0)                                                    //PCB
                    r[i][k] = r[i][k]+sides[k];
            }
        
        for(i = 0; i < n; i++)                                                      //for each particle
            for(j = 0; j < n; j++)                                                  //for each other particle
                if(j != i)                                                          //for all particles other than particle i
                    for(k = 0; k < d; k++)                                          //for each dimension
                        v[i][k] = vh[i][k] + a(i, j, k, d, r, sides, k_constant, m)*dt/2.;
   }
	
	return 0;
}

double a(int i, int j, int kk, int d, double r[][d], double sides[d], double k_constant, double m)
{
    int k;
    double sum_of_squares = 0.;
    for(k = 0; k < d; k++)
        sum_of_squares += pow(distance(i,j,k,d,r,sides),2.);
    return -k_constant/m*distance(i,j,kk,d,r,sides)/pow(sum_of_squares,3./2.);
}

double distance(int i, int j, int k, int d, double r[][d], double sides[d]){
    double u = r[j][k] - r[i][k];
    //don't know if this handles negative coordinates well, it seems to
    return fmin(u, sides[k]-u);
}
