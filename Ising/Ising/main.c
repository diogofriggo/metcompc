#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#define TOP 0
#define BOTTOM 1
#define LEFT 2
#define RIGHT 3
#define DIAGNOSTICS 1

double randomDoubleInclusive(double lowerLimit, double upperLimit);
    
const int nx = 10, ny = 10, tmax = 100;

int main(int argc, const char * argv[]) {
    srand((unsigned int)time(NULL));
    const char path[] = "/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS/metcompc/Ising/Ising/Ising.txt";
    int i, j, n, t;
    
    //board setup
    int board[nx][ny];
    for(i = 0; i < nx; i++)
        for(j = 0; j < ny; j++)
            board[i][j] = j+i*ny;

//    for(i = 0; i < nx; i++)
//        for(j = 0; j < ny; j++)
//            printf("%2d%s", board[i][j], (j+1)%ny == 0 ? "\n" : " ");

    //neighbours setup
    int v[nx*ny][4];
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
    
    //energies setup
    int beta = 1;
    int offset = 8;
    double energies[2*offset+1];
    for(i = 0; i < 2*offset+1; i++)
        energies[i] = 0.;
    energies[offset-8] = exp(-beta*-8); //0
    energies[offset-4] = exp(-beta*-4);
    energies[offset] = 1;
    energies[offset+4] = exp(-beta*+4);
    energies[offset+8] = exp(-beta*+8); //16
//    printf("%d:%.4f\n%d:%.4f\n%d:%.4f\n%d:%.4f\n%d:%.4f\n", -8, energies[offset-8], -4, energies[offset-4], 0, energies[offset], +4, energies[offset+4], +8, energies[offset+8]);
    
    //spins initialization
    int S[nx*ny+1];
    for(i = 0; i < nx*ny; i++)
        S[i] = randomDoubleInclusive(0,1) > 0.5 ? 1 : -1;

    for(i = 0; i < nx; i++)
        for(j = 0; j < ny; j++)
            printf("%2d%s", S[j+i*ny], (j+1)%ny == 0 ? "\n" : " ");

    //dynamics
    FILE *file = fopen(path, "w");
    for(i = 0; i < nx; i++)
        for(j = 0; j < ny; j++)
            fprintf(file, "%d %d %d\n", i, j, S[j+i*ny]);

    int energy;
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
                else if(energies[offset+energy] > 0.01){
                    S[n] = -S[n];
                    //printf(": metropolis flip\n");
                }
    //            else
    //                printf(": didn't flip (%.4f)\n", energies[offset+energy]);
                fprintf(file, "%d %d %d\n", i, j, S[n]);
        }
    }

    return 0;
}
                             
double randomDoubleInclusive(double lowerLimit, double upperLimit){
    return lowerLimit + (double)rand() / (double)RAND_MAX * (upperLimit - lowerLimit);
}