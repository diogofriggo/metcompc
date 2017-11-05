#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#define X 0
#define Y 1
#define Z 2

double a(int i, int j, int kk, int n, int d, double r[][d], double sides[d], double epsilon, double sigma);
double distance(int i, int j, int k, int d, double r[][d], double sides[d]);
double randomDoubleInclusive(double lowerLimit, double upperLimit);
void timeit(void (*fun)(void));
void program();

//double intermediates[100][49];
//int intermediate_curried[100][49];
double accelerations[100][49][2];
int main()
{
//    int i, j, n = 100;
//    for(i = 0; i < n; i++)
//        for(j = 0; j < 49; j++){
//            intermediates[i][j] = 0.;
//            intermediate_curried[i][j] = 0;
//        }
    //n=100, t=2000, dt = 0.001, d = 2
    //with partial memoization: 17.56s
    //without: 16.3562s
    timeit(program);
    return 0;
}

void program()
{
    srand((unsigned int)time(NULL));
    const char path[] = "/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS/metcompc/Aula14_3110/Results/LennardJones.txt";
    int t, i, j, k, ix, iy;
    const int n = 100, d = 2;
    //dt = 0.01 doesn't conserve energy! why? beats me!
    double dt = 0.001, tmax = 4000, epsilon = 10000., sigma = 1.;
    double v[n][d], vh[n][d], r[n][d];
    double sides[d], spacing = 10.;
    double total_acceleration = 0.;
    sides[0] = 100.;
    sides[1] = 100.;
    
//    int (*ptr)[100][49];
//    ptr = &intermediate_curried;
    
    //random velocities
    for(i = 0; i < n; i++)
        for(k = 0; k < d; k++)
            v[i][k] = randomDoubleInclusive(-10, 10);
    
    //position particles equally spaced from each other
    i = 0;
    for(iy = 0; iy < 10; iy++)
        for(ix = 0; ix < 10; ix++){
            r[i][X] = (ix+0.5)*spacing;
            r[i][Y] = (iy+0.5)*spacing;
            i++;
        }
    
    FILE *file = fopen(path,"w");
    
    for(t = 0; t < tmax; t++)
    {
        //Report
        for(i = 0; i < n; i++)
            for(k = 0; k < d; k++)
                fprintf(file, "%.2f ", r[i][k]);
        fprintf(file, "\n");
        
//        for(i = 0; i < n; i++)
//            for(j = 0; j < 49; j++)
//                (*ptr)[i][j] = 0;
        
        //Velocity-verlet
        for(i = 0; i < n; i++)
            for(k = 0; k < d; k++) {
                total_acceleration = 0.;
                for(j = 0; j < n; j++)
                    if(j != i)
                        total_acceleration += a(i, j, k, n, d, r, sides, epsilon, sigma)*dt/2.;
                vh[i][k] = v[i][k] + total_acceleration;
            }
        
        for(i = 0; i < n; i++)
            for(k = 0; k < d; k++){
                r[i][k] = r[i][k] + vh[i][k]*dt;
                if(r[i][k] >= sides[k])
                    r[i][k] = r[i][k]-sides[k];
                if(r[i][k] <= 0)
                    r[i][k] = r[i][k]+sides[k];
            }
        
//        for(i = 0; i < n; i++)
//            for(j = 0; j < 49; j++)
//                (*ptr)[i][j] = 0;

        for(i = 0; i < n; i++)
            for(k = 0; k < d; k++) {
                total_acceleration = 0.;
                for(j = 0; j < n; j++)
                    if(j != i)
                        total_acceleration += a(i, j, k, n, d, r, sides, epsilon, sigma)*dt/2.;
                v[i][k] = vh[i][k] + total_acceleration;
            }

        //PCB (can it be done together with the loop before the above?)
        for(i = 0; i < n; i++)
            for(k = 0; k < d; k++) {
                if(r[i][k] >= sides[k])
                    r[i][k] = r[i][k]-sides[k];
                if(r[i][k] <= 0)
                    r[i][k] = r[i][k]+sides[k];
            }
    }
}

double a(int i, int j, int kk, int n, int d, double r[][d], double sides[d], double epsilon, double sigma)
{
    if(j > i)
        return -accelerations[j][i][kk];
    
    double A = 0.;
//    double (*ptr)[100][49];
//    ptr = &intermediates;
//    if(intermediate_curried[i][j])
//        A = (*ptr)[i][j];
//    else{
        double sum_of_squares = 0.;
        int k;
        for(k = 0; k < d; k++)
            sum_of_squares += pow(distance(i, j, k, d, r, sides), 2.);
        A = 48*epsilon*pow(sum_of_squares,-7.)*pow(sigma,14)-0.5*pow(sum_of_squares,-4.)*pow(sigma,8);
//        intermediates[i][j] = A;
//        intermediate_curried[i][j] = 1;
//    }
    double acceleration = A*distance(i,j,kk,d,r,sides);
    //curry acceleration: they are antisymmetrical for i->j and j->i
    accelerations[i][j][kk] = acceleration;
    return acceleration;
}

double distance(int i, int j, int k, int d, double r[][d], double sides[d]){
    double u = r[i][k] - r[j][k];
    //return u > 0 ? fmin(u, sides[k]/2.) : -fmin(-u, sides[k]/2.)
    if(u > 0)
        return u < sides[k]/2. ? u : sides[k]/2.;
    return u > -sides[k]/2. ? u : -sides[k]/2.;
}

double randomDoubleInclusive(double lowerLimit, double upperLimit){
    return lowerLimit + (double)rand() / (double)RAND_MAX * (upperLimit - lowerLimit);
}

void timeit(void (*fun)(void)){
    clock_t begin = clock();
    fun();
    clock_t end = clock();
    double timeElapsed = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("%.6fs\n", timeElapsed);
}
