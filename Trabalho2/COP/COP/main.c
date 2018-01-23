#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#define ATTEMPT_FLIP 0
#define PRINT_LATTICE 0

void run();
void stats();
void setup();
void iterate();
void attemptFlip(int i, int j, int in, int jn);
int computePairEnergy(int i, int j, int in, int jn);
int computeEnergy(int i, int j, int in, int jn);
void report(FILE *file);
void setup();
void timeit(void (*fun)(void));
double randomDoubleInclusive(double lowerLimit, double upperLimit);
void printLattice();

const double rho = 0.5;
const int l = 50;
const int l2 = l/2;
const int iterations = 40;
const int J = 1;
const int kB = 1;
const int T = 1;
double beta = 0.7 ;//1/(kB*T);
const int offset = 12;
double exps[25];
int lattice[l][l];
int getPBCIndex(int i);
    
int main() {
    srand((unsigned int)time(NULL));
    timeit(run);
    int i;
    for(i = 1; i < 10; i++)
    {
        beta = i/10.;
        run();
    }
}

void run(){
    int i;
    setup();
    char path[300] = "/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS/metcompc/Trabalho2/COP/COP/";
    char name[50];
    sprintf(name, "cop%.1f.txt", beta);
    strcat(path, name);
    FILE *file = fopen(path, "w");
    report(file);
    for(i = 0; i < iterations; i++){
        iterate();
        report(file);
        #if PRINT_LATTICE
        printLattice();
        #endif
    }
    fclose(file);
}

void setup()
{
    int i, j;
    //precompute possible exponentials to save CPU time
    for(i = 0; i < 25; i++)
        exps[i] = 0.;
    int deltaEs[7] = {12, 8, 4, 0, -4, -8, -12};
    for(i = 0; i < 7; i++)
        exps[offset + deltaEs[i]] = exp(-beta*deltaEs[i]);
    
    for(i = 0; i < l; i++)
        for(j = 0; j < l; j++)
            lattice[i][j] = -1;
    
    for(i = 0; i < l; i++)
        for(j = 0; j < l2; j++)
            lattice[i][j] = 1;
}

//your neighbours are (top,right,bottom,left) = (i,j+1) (i+1,j) (i,j-1) (i-1,j)
void iterate(){
    int i, j;
    for(i = 0; i < l; i++)
        for(j = 1; j < (l-1); j++){ //avoid first and last row due to boundary conditions
            //for vertical movement we prevent a flip with a boundary neighbour since
            //that neighbour's spin is fixed to either +1 (bottom) or -1 (top)
            if((j+1) < (l-1))
                attemptFlip(i, j, i, j+1); //top
            if((j-1) > 0)
                attemptFlip(i, j, i, j-1); //bottom
            
            //for horizontal movement apply periodic boundary conditions
            if((i+1) < l) //if neighbour is not on the right boundary
                attemptFlip(i, j, i+1, j);
            else //else neighbour is at the right boundary, apply boundary condition
                attemptFlip(i, j, 0, j);
            
            if((i-1) >= 0) //if neighbour is not on the left boundary
                attemptFlip(i, j, i-1, j);
            else //else neighbour is at the right boundary, apply boundary condition
                attemptFlip(i, j, l-1, j);
        }
}


//TIP: find out the possible values of energy difference and estimate for what values of beta they would yield a non-zero, non-infinite exponential

//in = neighbour's i, jn = neighbour's j
//all indices are expected to be valid, that is, corrected for boundary conditions, before coming to this function
void attemptFlip(int i, int j, int in, int jn){
    if(lattice[i][j] == lattice[in][jn])
        return;
    int deltaE = computePairEnergy(i, j, in, jn);
//    printf("%d\n", deltaE);
//    printf("exps[12 + %d] = %.4f\n", deltaE, exps[offset + deltaE]);
    double r = randomDoubleInclusive(0.,1.);
    if(deltaE < 0 || exps[offset + deltaE] > r){
        int temp = lattice[i][j];
        lattice[i][j] = lattice[in][jn];
        lattice[in][jn] = temp;
        #if ATTEMPT_FLIP
        printf("success flip. ");
        #endif
    }
    #if ATTEMPT_FLIP
    else
        printf("failed flip. ");
    printf("%d < 0 || %.4f >= %.4f) (%d,%d,%d) <-> (%d,%d,%d)\n", deltaE, exp(-beta*deltaE), r, i, j, lattice[i][j], in, jn, lattice[in][jn]);
    if(deltaE > 1000 || deltaE < -1000)
        printf("ERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\n");
    #endif
}

int computePairEnergy(int i, int j, int in, int jn){
    return 2*J*(lattice[i][j]*computeEnergy(i, j, in, jn) + lattice[in][jn]*computeEnergy(in, jn, i, j));
}

//shit is incomplete
//in = neighbour's i, jn = neighbour's j
int computeEnergy(int i, int j, int in, int jn){
    int energy = 0;
    //compute the energy of the neighbours of (i,j) excluding (in,jn)
    if(in != i || jn != (j+1)) energy += lattice[i][j+1];
    if(in != i || jn != (j-1)) energy += lattice[i][j-1];
    if(in != getPBCIndex(i+1) || jn != j)
    {
        if((i+1) < l)
            energy += lattice[i+1][j];
        else
            energy += lattice[0][j];
    }
    if(in != getPBCIndex(i-1) || jn != j)
    {
        if((i-1) >= 0)
            energy += lattice[i-1][j];
        else
            energy += lattice[l-1][j];
    }
    return energy;
}

int getPBCIndex(int i){
    if(i < 0)
        return 49;
    if(i > (l-1))
        return 0;
    return i;
}

void report(FILE *file){
    int i, j;
    for(i = 0; i < l; i++)
        for(j = 0; j < l; j++)
            fprintf(file, "%d %d %d\n", i, j, lattice[i][j]);
}

void timeit(void (*fun)(void)){
    clock_t begin = clock();
    fun();
    clock_t end = clock();
    double timeElapsed = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("%.6f\n", timeElapsed);
}

double randomDoubleInclusive(double lowerLimit, double upperLimit){
    return lowerLimit + (double)rand() / (double)RAND_MAX * (upperLimit - lowerLimit);
}

void printLattice(){
    printf("Printing lattice...\n");
    int i, j;
    for(j = l-1; j >= 0; j--){
        printf("%2d ", j);
        for(i = 0; i < l; i++)
            printf("%s", lattice[i][j] == 1 ? "u" : "-");
        printf("\n");
    }
}
