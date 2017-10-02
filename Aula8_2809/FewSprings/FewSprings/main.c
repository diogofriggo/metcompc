#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double a(double M, double k, double k01, double *x, int this);

int main()
{
    char path[200] = "/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS/metcompc/Aula8_2809/Results/FewSprings.txt";
    double dt = 0.01;
    double tmax = 10000;
    double M = 1.;
    double L = 1.;
    //double aprev[2];
    
    double x[2];
    double v[2];
    double vh[2];
    
    x[0] = 0.;
    v[0] = 10.;
    vh[0] = 0.0;
    
    x[1] = 0.;
    v[1] = 10.;
    vh[1] = 0.0;
    
    double k = 1.;
    double k01 = 1.;
    int t = 0;
    
    FILE *file = fopen(path, "w");
    fprintf(file, "%.2f %.2f %.2f %.2f %.2f\n", t*dt, x[0]-L/.2, v[0], x[1]+L/.2, v[1]);
    for(t += dt; t < tmax; t++)
    {
        //IMPLEMENTATION 1: GENERAL
        vh[0] = v[0] + a(M, k, k01, x, 0)*dt/2.;
        vh[1] = v[1] + a(M, k, k01, x, 1)*dt/2.;
        
        x[0] = x[0] + vh[0]*dt;
        x[1] = x[1] + vh[1]*dt;
        
        v[0] = vh[0] + a(M, k, k01, x, 0) * dt/2.;
        v[1] = vh[1] + a(M, k, k01, x, 1) * dt/2.;
        
        //UNSTABLE: DAMPED OSCILLATION
        //IMPLEMENTATION 2: FORCE INDEPENDENT OF VELOCITY
//        aprev[0] = a(M, k, k01, x, 0) * dt/2.;
//        aprev[1] = a(M, k, k01, x, 1) * dt/2.;
//        
//        x[0] = x[0] + v[0]*dt + a(M, k, k01, x, 0) * pow(dt,2)/2.;
//        x[1] = x[1] + v[1]*dt + a(M, k, k01, x, 1) * pow(dt,2)/2.;
//        v[0] = v[0] + (aprev[0] + a(M, k, k01, x, 0)) * dt/2.;
//        v[1] = v[1] + (aprev[1] + a(M, k, k01, x, 1)) * dt/2.;
        
        fprintf(file, "%.2f %.2f %.2f %.2f %.2f\n", t*dt, x[0]-L/.2, v[0], x[1]+L/.2, v[1]);
    }
    fclose(file);
    return 0;
}

double a(double M, double k, double k01, double *x, int this)
{
    return k01*x[abs(1-this)] - (k+k01)*x[this];
}