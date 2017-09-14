#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void diffusion();
void run(double *x, double *d, int nXs, int nTs, double k);
void describe(int nXs, double *x);
void report(char *fileName, double *x, int nXs, double dx, double k, double maxT);

int main() {
    diffusion();
    return 0;
}

void diffusion(void (*runExample)(double*, double*, int, int, double, double, double)){
    int k, i, L = 100;
    double dt = 0.1, dx = 1;
    //if time is divided in intervals of 0.1 between 0 and 100 then nTs = [0, 0.1, ... , 999.9] -> 1000 steps
    //if space is divided in intervals of 0.1 between 0 and 100 then nXs = [0, 0.1, ... , 999.9, 1000] -> 1001 positions
    int nTs, nXs = (int)round(L/dx) + 1;
    const int maxTsLength = 6, ksLength = 10;
    double maxTs[maxTsLength], ks[ksLength];
    
    double *d, *x;
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
    
    for(k = 0; k < ksLength; k++)
        for(i = 0; i < maxTsLength; i++){
            nTs = (int)round(maxTs[i]/dt);
            describe(nXs, x);
            run(x, d, nXs, nTs, ks[k]);
            report("CDiffusionPCBK%.1ft%.0f.txt", x, nXs, dx, ks[k], maxTs[i]);

        }
    free(x);
    free(d);
}

void describe(int nXs, double *x){
    int i;
    
    //initial condition
    //f(L/4 < x < 3L/4, t = 0) = 1
    //Obs: performing integer division
    int start = (nXs-1)/4 + 1;
    int end = 3*(nXs-1)/4 + 1;
    for(i = 0; i < nXs; i++)
        x[i] = i >= start && i <= end;
}

void run(double *x, double *d, int nXs, int nTs, double k)
{
    register int t, j;
    
    for(t = 0; t < nTs; t++){
        //for some mysterious reason if I use d = x (d a pointer), it hides the instability of the method
        //we have to use instead an explicit for to copy from x to d:
        for(j = 0; j < nXs; j++)
            d[j] = x[j];
        //avoid the two boundaries (x=0, x=L) so as to not mess up boundary conditions
//        for(j = 1; j < nXs-1; j++){
//            x[j] = d[j] + k*(d[j+1] - 2.*d[j] + d[j-1]);
//        }
        //implements periodic boundary condition:
        for(j = 0; j < nXs; j++){
            if(j == 0)
                x[j] = d[j] + k*(d[j+1] - 2.*d[j] + d[nXs-1]);
            else if(j == nXs-1)
                x[j] = d[j] + k*(d[0] - 2.*d[j] + d[j-1]);
            else
                x[j] = d[j] + k*(d[j+1] - 2.*d[j] + d[j-1]);
        }
    }
}

void report(char *fileName, double *x, int nXs, double dx, double k, double maxT){
    char path[200] = "/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS COMPUTACIONAIS C/metcompc/Aula2_3108/Results/";
    char name[50];
    sprintf(name, fileName, k, maxT);
    strcat(path, name);
    FILE *file = fopen(path, "w");
    int j;
    for(j = 0; j < nXs; j++)
        fprintf(file, "%.4f %.4f\n", j*dx, x[j]);
    fclose(file);
}
