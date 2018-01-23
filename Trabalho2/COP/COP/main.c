#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>

void run();
void setup();
void iterate();
void attemptFlip(int i, int j, int in, int jn);
int computePairEnergy(int i, int j, int in, int jn);
int computeEnergy(int i, int j, int in, int jn);
void report(FILE *file);
void setup();
void timeit(void (*fun)(void));
double randomDoubleInclusive(double lowerLimit, double upperLimit);

const double rho = 0.5;
const int l = 50;
const int l2 = l/2;
const int iterations = 100;
const int J = 1;
const int kB = 1;
const int T = 1;
//beta >= 187 forbids any transition
//changing boundary conditions beta >= 94 forbids any transition
const int beta = 1000000000;//1/(kB*T);

int lattice[l][l];

int main() {
    srand((unsigned int)time(NULL));
    timeit(run);
    int i;
    for(i = 0; i < 100; i++){
        printf("%.2f ", randomDoubleInclusive(0.,1.));
    }
    return 0;
}

void run(){
    int i;
    
    setup();
    const char path[] = "/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS/metcompc/Trabalho2/COP/COP/cop.txt";
    FILE *file = fopen(path, "w");
    for(i = 0; i < iterations; i++){
        iterate();
//        if(i%10==0)
            report(file);
    }
    fclose(file);
}

void setup()
{
    int i, j;
    
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
        for(j = 0; j < l; j++){
            if(j+1 < l)
                attemptFlip(i, j, i, j+1); //top
            if(i+1 < l)
                attemptFlip(i, j, i+1, j); //right
            if(j-1 > 0)
                attemptFlip(i, j, i, j-1); //bottom
            if(i-1 > 0)
                attemptFlip(i, j, i-1, j); //left
        }
}

//in = neighbour's i, jn = neighbour's j
void attemptFlip(int i, int j, int in, int jn){
    int deltaE = computePairEnergy(i, j, in, jn);
    if(deltaE < 0 || exp(-beta*deltaE) > randomDoubleInclusive(0.,1.)){
        int temp = lattice[i][j];
        lattice[i][j] = lattice[in][jn];
        lattice[in][jn] = temp;
    }
}

int computePairEnergy(int i, int j, int in, int jn){
    return 2*J*(lattice[i][j]*computeEnergy(i, j, in, jn) + lattice[in][jn]*computeEnergy(in, jn, i, j));
}

//in = neighbour's i, jn = neighbour's j
int computeEnergy(int i, int j, int in, int jn){
    int energy = 0;
    if(!(in == i && jn == j+1)) //top
        energy += j+1 < l ? lattice[i][j+1] : -1; //keeps interface from wandering off
    if(!(in == i+1 && jn == j)) //right
        energy += i+1 < l ? lattice[i+1][j] : lattice[0][j]; //horizontal periodic boundary condition
    if(!(in == i && jn == j-1)) //bottom
        energy +=  j-1 > 0 ? lattice[i][j-1] : 1; //keeps interface from wandering off
    if(!(in == i-1 && jn == j)) //left
        energy += i-1 > 0 ? lattice[i-1][j] : lattice[l-1][j]; //horizontal periodic boundary condition
    return energy;
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