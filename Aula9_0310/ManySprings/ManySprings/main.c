#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double a(double *x, int i, double d, int n, double k[][n], double M);
double positionOf(double *x, int index, double d);
double forceConstantOf(int i, int j, int n, double k[][n]);

int main()
{
    char path[200] = "/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS/metcompc/Aula9_0310/Results/ManySprings.txt";
    
    double dt = 0.001, tmax = 100000, M = 1., d = 1.;
    int t, i, j;
    const int n = 10;
    
    //default force constants = 1.
    double k[n][n];
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            k[i][j] = 1.;
    
    //default displacements and velocities = 0.
    double x[n], v[n], vh[n];
    for(i = 0; i < n; i++)
        x[i] = v[i] = vh[i] = 0.;

    x[1] = -0.4;
    x[2] = -0.1;
    x[3] = 0.4;
    
    FILE *file = fopen(path, "w");
    for(t = 0; t < tmax; t++)
    {
        //Report
        fprintf(file, "%.2f ", t*dt);
        for(i = 1; i < n-1; i++)
            fprintf(file, "%.2f %.2f ", positionOf(x,i,d), v[i]);
        fprintf(file, "\n");
        
        //Velocity-verlet
        for(i = 1; i < n-1; i++)
            vh[i] = v[i] + a(x, i, d, n, k, M)*dt/2.;
        for(i = 1; i < n-1; i++)
            x[i] = x[i] + vh[i]*dt;
        for(i = 1; i < n-1; i++)
            v[i] = vh[i] + a(x, i, d, n, k, M)*dt/2.;
    }
    fclose(file);
    return 0;
}

double a(double *x, int i, double d, int n, double k[][n], double M)
{
    return forceConstantOf(i,i-1,n,k)*positionOf(x,i-1,d)
         + forceConstantOf(i,i+1,n,k)*positionOf(x,i+1,d)
        - (forceConstantOf(i,i-1,n,k)+forceConstantOf(i,i+1,n,k))*positionOf(x,i,d);
}

double positionOf(double *x, int i, double d)
{
    return x[i]+i*d;
}

double forceConstantOf(int i, int j, int n, double k[][n])
{
    return k[i][j];
}
