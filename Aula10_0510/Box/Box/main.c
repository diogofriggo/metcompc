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
double randomDoubleInclusive(double lowerLimit, double upperLimit);

int main()
{
    srand((unsigned int)time(NULL));
    char path[200] = "/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS/metcompc/Aula10_0510/Results/Box.txt";
    //n = particles, d = dimensions
    //i = particle index, j = other particle index, k = dimension index
    //vh = velocity at half step (part of the velocity verlet algorithm)
    int t, i, j, k;
    const int n = 2, d = 2;
    double dt = 0.1, tmax = 290, k_constant = 1., b = 1., m = 1.;
    double vh[n][d], v[n][d], r[n][d], sides[d];
    
    //default values
    for(i = 0; i < n; i++)
        for(k = 0; k < d; k++)
        {
            vh[i][k] = 0.;
            v[i][k] = 0.;
            r[i][k] = 0.;
        }
    
    sides[0] = 10.;
    sides[1] = 10.;
    
    //initial conditions: incoming x particle hits still particle
    r[0][Y] = 7*b;
    v[0][X] = 1.;
    r[1][X] = 4.;
    r[1][Y] = 5*b;
    
    //flux of incoming particles separated by b
    //TODO: randomly place particles?
//    r[0][Y] = b;
//    r[1][Y] = 2*b;
//    r[2][Y] = 3*b;
//    r[3][Y] = 4*b;
//    r[4][Y] = 0*b;
//    r[5][Y] = -b;
//    r[6][Y] = -2*b;
//    r[7][Y] = -3*b;
//    r[8][Y] = -4*b;
//    r[9][X] = 10.;
//    v[0][X] = 1.;
//    v[1][X] = 1.;
//    v[2][X] = 1.;
//    v[3][X] = 1.;
//    v[4][X] = 1.;
//    v[5][X] = 1.;
//    v[6][X] = 1.;
//    v[7][X] = 1.;
//    v[8][X] = 1.;
    
    FILE *file = fopen(path, "w");
    for(t = 0; t <= tmax; t++)
    {
        //Report
        for(i = 0; i < n; i++)                                                      //for each particle
            for(k = 0; k < d; k++)                                                  //for each dimension
                fprintf(file, "%.2f ", r[i][k]);
        fprintf(file, "\n");
        
        //Velocity-verlet
        for(i = 0; i < n; i++)                                                      //for each particle
            for(j = 0; j < n; j++)                                                  //for each other particle
                if(j != i)                                                          //for all particles other than particle i
                    for(k = 0; k < d; k++)                                          //for each dimension
                        vh[i][k] = v[i][k] + a(i, j, k, d, r, sides, k_constant, m)*dt/2.;
        
        //TODO: implement PCB
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
    fclose(file);
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
