#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#define TOP 0
#define BOTTOM 1
#define LEFT 2
#define RIGHT 3
#define DIAGNOSTICS 1

double randomDoubleInclusive(double lowerLimit, double upperLimit);
void run();
void run_snapshots();
void run_magnetization();
void setup();
void timeit(void (*fun)(void));

const int nx = 100, ny = 100, tmax = 20000;
double const kB = 1.;
double T = 1.;
int S[nx*ny+1];
int board[nx][ny];
const int offset = 8;
double energies[2*offset+1];
int v[nx*ny][4];

int main(int argc, const char * argv[]) {
    //timeit(run_snapshots);
    //timeit(run_magnetization);
    T = 1.;
    timeit(run);
    return 0;
}

void run_magnetization()
{
    //for each temperature run as long as necessary to reach equilibrium
    //how to determine if equilibrium has been reached?
    //once there measure magnetization and write to file
    const char path[] = "/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS/metcompc/Ising/Ising/IsingMagnetization.txt";
    FILE *file = fopen(path, "w");
    int magnetization = 0;
    int i, j;
    setup();
    for(i = 0; i < 1000; i++){
        T = i/100.;
        run();
        magnetization = 0;
        for(j = 0; j < nx*ny; j++)
            magnetization += S[j];
        fprintf(file, "%.4f %.4f\n", T, ((double)magnetization)/(nx*ny));
    }
    fclose(file);
}

void run_snapshots()
{
    int i;
    double ts[3];
    ts[0] = 1.;
    ts[1] = 2./log(1.+sqrt(2.));
    //printf("%.4f\n", ts[1]);
    ts[2] = 10.;
    setup();
    for(i = 0; i < 3; i++){
        T = ts[i];
        run();
    }
}

void setup()
{
    int i, j, n;
    
    //board setup
    for(i = 0; i < nx; i++)
        for(j = 0; j < ny; j++)
            board[i][j] = j+i*ny;
    
    //<print board>
    //    for(i = 0; i < nx; i++)
    //        for(j = 0; j < ny; j++)
    //            printf("%2d%s", board[i][j], (j+1)%ny == 0 ? "\n" : " ");
    //</print board>
    
    //neighbours setup
    for(i = 0; i < nx; i++)
        for(j = 0; j < ny; j++){
            n = j+i*ny;
            v[n][LEFT] = j == 0 ? board[i][ny-1] : board[i][j-1];
            v[n][RIGHT] = j == (ny-1) ? board[i][0] : board[i][j+1];
            v[n][TOP] = i == 0 ? board[nx-1][j] : board[i-1][j];
            v[n][BOTTOM] = i == (nx-1) ? board[0][j] : board[i+1][j];
        }
    
    //    for(i = 0; i < nx*ny; i++)
    //        printf("%5d\n%2d %2d %2d\n%5d\n", v[i][TOP], v[i][LEFT], i, v[i][RIGHT], v[i][BOTTOM]);
    
}

void run()
{
    srand((unsigned int)time(NULL));
    int i, j, n, t;
    
    //spins initialization
    for(i = 0; i < nx*ny; i++)
        S[i] = randomDoubleInclusive(0,1) > 0.5 ? 1 : -1;
    
    //<print S>
    //    for(i = 0; i < nx; i++)
    //        for(j = 0; j < ny; j++)
    //            printf("%2d%s", S[j+i*ny], (j+1)%ny == 0 ? "\n" : " ");
    //</print S>
    
    //energies setup
    double beta = 1./(kB*T);
    for(i = 0; i < 2*offset+1; i++)
        energies[i] = 0.;
    energies[offset-8] = exp(-beta*-8); //0
    energies[offset-4] = exp(-beta*-4);
    energies[offset] = 1;
    energies[offset+4] = exp(-beta*+4);
    energies[offset+8] = exp(-beta*+8); //16
    //    printf("%d:%.4f\n%d:%.4f\n%d:%.4f\n%d:%.4f\n%d:%.4f\n", -8, energies[offset-8], -4, energies[offset-4], 0, energies[offset], +4, energies[offset+4], +8, energies[offset+8]);

    //<report>
//    const char path[200] = "/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS/metcompc/Ising/Ising/Ising%.1f.txt";
//    char filePath[300];
//    sprintf(filePath, path, T);
//    FILE *file = fopen(filePath, "w");
//    for(i = 0; i < nx; i++)
//        for(j = 0; j < ny; j++)
//            fprintf(file, "%d %d %d\n", i, j, S[j+i*ny]);
    //</report>

    //<report magnetization>
        const char path[200] = "/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS/metcompc/Ising/Ising/IsingMagnetizationEquilibrium.txt";
        FILE *file = fopen(path, "w");
    //</report magnetization>

    //dynamics
    int energy, magnetization = 0;

    for(t = 0; t < tmax; t++){
        for(i = 0; i < nx; i++)
            for(j = 0; j < ny; j++){
                n = j+i*ny;
                energy = 2*S[n]*(S[v[n][LEFT]]+S[v[n][RIGHT]]+S[v[n][TOP]]+S[v[n][BOTTOM]]);
                //printf("%d = 2*%d*(%d + %d + %d + %d)", energy, S[n], S[v[n][LEFT]], S[v[n][RIGHT]], S[v[n][TOP]], S[v[n][BOTTOM]]);
                if(energy <= 0){
                    S[n] = -S[n];
                    //printf(": flipped\n");
                }
                else if(energies[offset+energy] < randomDoubleInclusive(0, 1)){
                    S[n] = -S[n];
                    //printf(": metropolis flip\n");
                }
                //            else
                //                printf(": didn't flip (%.4f)\n", energies[offset+energy]);
                //report spin
                //fprintf(file, "%d %d %d\n", i, j, S[n]);
            }
        //<report magnetization>
        magnetization = 0;
        for(j = 0; j < nx*ny; j++)
            magnetization += S[j];
        fprintf(file, "%d %.4f\n", t, ((double)magnetization)/(nx*ny));
        //</report magnetization>
    }
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