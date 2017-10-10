#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#define X 0
#define Y 1
#define Z 2

double a(int i, int j, int kk, int d, double r[][d], double sides[d], double k_constant);
double distance(int i, int j, int k, int d, double r[][d], double sides[d]);
double randomDoubleInclusive(double lowerLimit, double upperLimit);

int main()
{
    srand((unsigned int)time(NULL));
    //char path[200] = "/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS/metcompc/Aula10_0510/Results/Box.txt";
    //n = particles, d = dimensions
    //i = particle index, j = other particle index, k = dimension index
    //vh = velocity at half step (part of the velocity verlet algorithm)
    int t, i, j, k;
    const int n = 20, d = 2;
    double dt = 0.1, tmax = 400, k_constant = 1., b = 1., m = 1., T = 1.;
    double vh[n][d], v[n][d], r[n][d], sides[d];
    sides[0] = 100.;
    sides[1] = 100.;    
    
    //default values
    for(i = 0; i < n; i++)
        for(k = 0; k < d; k++)
        {
            vh[i][k] = 0.;
            v[i][k] = 0.;
            r[i][k] = 0.;
        }
    
    //sort random positions
    //for(i = 0; i < n; i++)                                                      	//for each particle
        //for(k = 0; k < d; k++) {                                                	//for each dimension
            //r[i][k] = randomDoubleInclusive(0, sides[k]);
	//square lattice
	//double volume = 1.;	
	//double parameters[d];	
	//double parameter = sides[0];
	//double z = 0.;
    //for(k = 0; k < d; k++) 		                                                	//for each dimension
	//	volume *= sides[k];
    //for(k = 0; k < d; k++) 		                                                	//for each dimension
    //    parameters[k] = ceil(sides[k]/sqrt(volume/n));
    //for(k = 0; k < d; k++) 		                                                	//for each dimension
    //    parameter = fmin(parameter, parameters[k]);	

	//start in 2D then generalize	
	//for(i = 0; i < n; i++){															//for each particle
	//	z = i/parameters[k];
	//	r[i][X] = parameter*(parameters[0]*(z-floor(z))+.5);
	//	r[i][Y] = parameter*(ceil(n/parameters[0])-.5);
	//}

	//sort random positions
	for(i = 0; i < n; i++)															//for each particle
		for(k = 0; k < d; k++)														//for each dimension
			r[i][k] = randomDoubleInclusive(0,sides[k]);

	//for(k = 0; k < d; k++)															//for each particle
	//	printf("%.2f ", parameters[k]);
	//printf("%.2f \n", parameter);

	//for(i = 0; i < n; i++)															//for each particle
	//	printf("(%.2f, %.2f)\n", r[i][X], r[i][Y]);

	//sort gaussian velocities
    for(i = 0; i < n; i++)                                                      //for each particle
        for(k = 0; k < d; k++)                                                 	//for each dimension
            v[i][k] = sqrt(T)*randomDoubleInclusive(0, 1);

	//printf("set size square\n");
	//printf("set xrange [0:%d]\n", sides[0]);
	//printf("set yrange [0:%d]\n", sides[1]);
	//printf("plot \"-\" w p pt 7 ps 2\n");
    //for(i = 0; i < n; i++){                                                     //for each particle
    //    for(k = 0; k < d; k++)                                                  //for each dimension
	//		printf("%.2f ", r[i][k]);
	//	printf("\n");
	//}
	//printf("&\n");
    //
	    
    //FILE *file = fopen("mru.txt", "w");
    for(t = 0; t <= tmax; t++)
    {
		printf("set size square\n");
		printf("set xrange [0:%d]\n", sides[0]);
		printf("set yrange [0:%d]\n", sides[1]);
		printf("plot \"-\" w p pt 7 ps 2\n");
        //Report
        for(i = 0; i < n; i++){                                                      //for each particle
            for(k = 0; k < d; k++)                                                  //for each dimension
				printf("%.2f ", r[i][k]);
			printf("\n");
//                fprintf(file, "%.2f ", r[i][k]);
  //      fprintf(file, "\n");
		}
        printf("e\n");
        
        //Velocity-verlet
        for(i = 0; i < n; i++)                                                      //for each particle
            for(j = 0; j < n; j++)                                                  //for each other particle
                if(j != i)                                                          //for all particles other than particle i
                    for(k = 0; k < d; k++)                                          //for each dimension
                        vh[i][k] = v[i][k] + a(i, j, k, d, r, sides, k_constant)*dt/2.;
        
        for(i = 0; i < n; i++)                                                      //for each particle
            for(k = 0; k < d; k++) {                                                //for each dimension
                r[i][k] = r[i][k] + vh[i][k]*dt;
                //if(r[i][k] >= sides[k])                                             //PCB
                //    r[i][k] = r[i][k]-sides[k];
                //if(r[i][k] <= 0)                                                    //PCB
                //    r[i][k] = r[i][k]+sides[k];
            }
        
        for(i = 0; i < n; i++)                                                      //for each particle
            for(j = 0; j < n; j++)                                                  //for each other particle
                if(j != i)                                                          //for all particles other than particle i
                    for(k = 0; k < d; k++){                                         //for each dimension
                        v[i][k] = vh[i][k] + a(i, j, k, d, r, sides, k_constant)*dt/2.;
						if(v[i][k] >= sides[k] || v[i][k] <= 0)
							v[i][k] = -v[i][k];
					}
    }
    //fclose(file);
	for(i = 0; i < 10; i++)
		printf("%.2f ", randomDoubleInclusive(0,10));
    return 0;
}

double a(int i, int j, int kk, int d, double r[][d], double sides[d], double k_constant)
{
	return 0;
    //int k;
    //double sum_of_squares = 0.;
    //for(k = 0; k < d; k++)
    //    sum_of_squares += pow(distance(i,j,k,d,r,sides),2.);
    //return -k_constant*distance(i,j,kk,d,r,sides)/pow(sum_of_squares,3./2.);
}

double distance(int i, int j, int k, int d, double r[][d], double sides[d]){
    double u = r[j][k] - r[i][k];
    //don't know if this handles negative coordinates well, it seems to
    //return fmin(u, sides[k]-u);//PCB
	return u;
}

double randomDoubleInclusive(double lowerLimit, double upperLimit){
    return lowerLimit + (double)rand() / (double)RAND_MAX * (upperLimit - lowerLimit);
}

int randomIntegerInclusive(int lowerLimit, int upperLimit){
    return lowerLimit + rand() % (upperLimit+1 - lowerLimit);
}

int randomIntegerLeftInclusive(int lowerLimit, int upperLimit){
    if(lowerLimit == upperLimit) return 0;
    return lowerLimit + rand() % (upperLimit - lowerLimit);
}
