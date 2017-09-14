#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void stationaryDiffusion();
void describeStationaryDiffusion(int n, double a[n-1], double b[n], double c[n-1], double d[n]);
void describeStationaryDiffusionBoundaryConditions1(int n, double a[n-1], double b[n], double c[n-1], double d[n]);
void describeStationaryDiffusionBoundaryConditions2(int n, double a[n-1], double b[n], double c[n-1], double d[n]);
void reportSolution(int n, double x[n]);

void nonStationaryDiffusion();
void describeNonStationaryDiffusionBoundaryConditions(double k, int nXs, double a[nXs-1], double b[nXs], double c[nXs-1], double d[nXs]);
void runNonStationaryDiffusion(double *a, double *b, double *c, double *d, double *x, int nXs, int nTs, double k);
void reportNonStationaryDiffusion(double *x, int nXs, double dx, double k, double maxT);
void thomasSolveTridiagonalMatrix(int n, double a[n-2], double b[n], double c[n-1], double d[n], double x[n]);

int main(int argc, const char * argv[]) {
    stationaryDiffusion();
    nonStationaryDiffusion();
    return 0;
}

/****************************
 * NON-STATIONARY DIFFUSION *
 ****************************/

void nonStationaryDiffusion(){
    int k, i, L = 100;
    double dt = 0.1, dx = 1;
    //if time is divided in intervals of 0.1 between 0 and 100 then nTs = [0, 0.1, ... , 999.9] -> 1000 steps
    //if space is divided in intervals of 0.1 between 0 and 100 then nXs = [0, 0.1, ... , 999.9, 1000] -> 1001 positions
    int nTs, nXs = (int)round(L/dx) + 1;
    const int maxTsLength = 6, ksLength = 10;
    double maxTs[maxTsLength], ks[ksLength];
    
    double *a, *b, *c, *d, *x;
    a = (double*)calloc(nXs-1, sizeof(double));
    b = (double*)calloc(nXs, sizeof(double));
    c = (double*)calloc(nXs-1, sizeof(double));
    d = (double*)calloc(nXs, sizeof(double));
    x = (double*)calloc(nXs, sizeof(double));

    for(k = 0; k < ksLength; k++)
        ks[k] = 0.1*(k+1);
    
    maxTs[0] = 2;
    maxTs[1] = 10;
    maxTs[2] = 50;
    maxTs[3] = 200;
    maxTs[4] = 1000;
    maxTs[5] = 2000;
    
    for(k = 0; k < ksLength; k++){
        for(i = 0; i < maxTsLength; i++){
            nTs = (int)round(maxTs[i]/dt);
            runNonStationaryDiffusion(a, b, c, d, x, nXs, nTs, ks[k]);
            reportNonStationaryDiffusion(x, nXs, dx, ks[k], maxTs[i]);
        }
    }
    
    free(a);
    free(b);
    free(c);
    free(d);
    free(x);
}

void reportNonStationaryDiffusion(double *x, int nXs, double dx, double k, double maxT){
    char path[200] = "/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS COMPUTACIONAIS C/metcompc/Aula3_0509/Results/";
    char name[50];
    sprintf(name, "CThomasDiffusionK%.1ft%.0f.txt", k, maxT);
    strcat(path, name);
    FILE *file = fopen(path, "w");
    int j;
    for(j = 0; j < nXs; j++)
        fprintf(file, "%.4f %.4f\n", j*dx, x[j]);
    fclose(file);
}

void runNonStationaryDiffusion(double *a, double *b, double *c, double *d, double *x, int nXs, int nTs, double k){
    describeNonStationaryDiffusionBoundaryConditions(k, nXs, a, b, c, d);
    int t;
    for(t = 0; t < nTs; t++){
        thomasSolveTridiagonalMatrix(nXs, a, b, c, d, x);
        d = x;
    }
}

void describeNonStationaryDiffusionBoundaryConditions(double k, int nXs, double a[nXs-1], double b[nXs], double c[nXs-1], double d[nXs]){
    int i;
    for(i = 0; i < nXs-2; i++)
        a[i] = -k;
    for(i = 1; i < nXs-1; i++)
        b[i] = 1+2*k;
    for(i = 1; i < nXs-1; i++)
        c[i] = -k;
    for(i = 1; i < nXs-1; i++)
        d[i] = 0;
    
    //specifying left boundary condition
    b[0] = 1;
    c[0] = 0;
    d[0] = 1;
    
    //specifying right boundary condition
    a[nXs-2] = 0;
    b[nXs-1] = 1;
    d[nXs-1] = 0;
}

/************************
 * STATIONARY DIFFUSION *
 ************************/

void stationaryDiffusion(){
    int n = 6;
    double *a;
    double *b;
    double *c;
    double *d;
    double *x;
    a = (double*)calloc(n-1, sizeof(double));
    b = (double*)calloc(n, sizeof(double));
    c = (double*)calloc(n-1, sizeof(double));
    d = (double*)calloc(n, sizeof(double));
    x = (double*)calloc(n, sizeof(double));
    
    describeStationaryDiffusion(n, a, b, c, d);
    
    describeStationaryDiffusionBoundaryConditions1(n, a, b, c, d);
    thomasSolveTridiagonalMatrix(n, a, b, c, d, x);
    reportSolution(n, x);
    
    describeStationaryDiffusionBoundaryConditions2(n, a, b, c, d);
    thomasSolveTridiagonalMatrix(n, a, b, c, d, x);
    reportSolution(n, x);

    free(a);
    free(b);
    free(c);
    free(d);
    free(x);
}

void describeStationaryDiffusion(int n, double a[n-1], double b[n], double c[n-1], double d[n]){
    int i;
    for(i = 1; i < n-1; i++){
        a[i-1] = -1;
        b[i] = 2;
        c[i] = -1;
        d[i] = 0;
    }
}

void describeStationaryDiffusionBoundaryConditions1(int n, double a[n-1], double b[n], double c[n-1], double d[n]){
    //specifying left boundary condition
    b[0] = 1;
    c[0] = 0;
    d[0] = 1;
    
    //specifying right boundary condition
    a[n-2] = 0;
    b[n-1] = 1;
    d[n-1] = 0;
}

void describeStationaryDiffusionBoundaryConditions2(int n, double a[n-1], double b[n], double c[n-1], double d[n]){
    //specifying left boundary condition
    b[0] = 1;
    c[0] = 0;
    d[0] = 0.5;
    
    //specifying right boundary condition
    a[n-2] = 0;
    b[n-1] = 1;
    d[n-1] = 1;
}

void reportSolution(int n, double x[n]){
    int i;
    for(i = 0; i < n; i++)
        printf("%d: %.4f ", i, x[i]);
    printf("\n");
}

//it writes the result to the x[n] vector
void thomasSolveTridiagonalMatrix(int n, double a[n-1], double b[n], double c[n-1], double d[n], double x[n]){
    int i;
    double denominator;
    double c_[n];
    double d_[n];
    
    //1 - compute coefficients
    c_[0] = c[0]/b[0];
    d_[0] = d[0]/b[0];
    for(i = 1; i < (n-1); i++){
        //there's a very subtle issue when manipulating a because it is one index behind
        //so to instead of a[i] as prescribed by thomas' formula we must use a[i-1]
        denominator = b[i]-c_[i-1]*a[i-1];
        c_[i] = c[i]/denominator;
        d_[i] = (d[i]-d_[i-1]*a[i-1])/denominator;
    }
    //d[] has one further step:
    denominator = b[n-1]-c_[n-2]*a[n-2];
    d_[n-1] = (d[n-1]-d_[n-2]*a[n-2])/denominator;
    
    //2 - substitute back to find solution vector
    x[n-1] = d_[n-1];
    for(i = n-2; i >= 0; i--)
        x[i] = d_[i]-c_[i]*x[i+1];
}

