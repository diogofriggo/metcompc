#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#define X 0
#define Y 1
#define Z 2

double randomDoubleInclusive(double lowerLimit, double upperLimit);

int main()
{
    srand((unsigned int)time(NULL));
    int t, i, k;
    const int n = 20, d = 2;
    double dt = 0.1, tmax = 2000, T = 0.1;
    double v[n][d], r[n][d], sides[d], spacing[d];
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

	printf("%.2f %.2f\n", spacing[0], spacing[1]);
	for(i = 0; i < n; i++)
		printf("%d \% %d = \n", i, );
	//rectangular positions
	for(i = 0; i < n; i++)														//for each particle
		for(k = 0; k < d; k++)													//for each dimension
			r[i][k] = (i%100)*spacing[k];
	
	//sort random velocities
    for(i = 0; i < n; i++)                                                      //for each particle
        for(k = 0; k < d; k++)                                                 	//for each dimension
            v[i][k] = sqrt(T)*randomDoubleInclusive(-1, 1);

	for(i = 0; i < n; i++)
		printf("(%.2f, %.2f) (%.2f, %.2f)\n", r[i][X], r[i][Y], v[i][X], v[i][Y]);

	double excess = 0.;
    for(t = 0; t <= -1; t++)
    {
		printf("set size square\n");
		printf("set xrange [0:%d]\n", (int)sides[0]);
		printf("set yrange [0:%d]\n", (int)sides[1]);
		printf("plot \"-\" w p pt 7 ps 2\n");

        for(i = 0; i < n; i++){                                                 //for each particle
            for(k = 0; k < d; k++)                                              //for each dimension
				printf("%.2f ", r[i][k]);
			printf("\n");
		}
        printf("e\n");

/*        for(i = 0; i < n; i++)                                                  //for each particle*/
/*            for(k = 0; k < d; k++) {                                            //for each dimension*/
/*                r[i][k] += v[i][k]*dt;*/
/*				*/
/*                if(r[i][k] >= sides[k]){*/
/*                    v[i][k] = -v[i][k];*/
/*					r[i][k] -= r[i][k]-sides[k];*/
/*				}*/
/*				if(r[i][k] <= 0.){*/
/*                    v[i][k] = -v[i][k];*/
/*					r[i][k] = -r[i][k];*/
/*				}*/
/*            }*/
    }
	
	return 0;
}

double randomDoubleInclusive(double lowerLimit, double upperLimit){
    return lowerLimit + (double)rand() / (double)RAND_MAX * (upperLimit - lowerLimit);
}
