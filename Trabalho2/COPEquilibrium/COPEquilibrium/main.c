#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#define ATTEMPT_FLIP 0
#define PRINT_LATTICE 0

typedef unsigned long uint32;

#define N              (624)                 // length of state vector
#define M              (397)                 // a period parameter
#define K              (0x9908B0DFU)         // a magic constant
#define hiBit(u)       ((u) & 0x80000000U)   // mask all but highest   bit of u
#define loBit(u)       ((u) & 0x00000001U)   // mask all but lowest    bit of u
#define loBits(u)      ((u) & 0x7FFFFFFFU)   // mask     the highest   bit of u
#define mixBits(u, v)  (hiBit(u)|loBits(v))  // move hi bit of u to hi bit of v

static uint32   state[N+1];     // state vector + 1 extra to not violate ANSI C
static uint32   *next;          // next random value is computed from here
static int      left = -1;      // can *next++ this many times before reloading
void seedMT(uint32 seed);
uint32 reloadMT(void);
uint32 randomMT(void);

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
int getPBCIndex(int i);
void measureEquilibrium(int step, FILE *file);

const double rho = 0.6;
const int l = 100;
//const int l2 = l/2;
const int iterations = 2000;
const int J = 1;
const int kB = 1;
const int T = 1;

double beta = 0.2 ;//1/(kB*T);
const int offset = 12;
double exps[25];
int lattice[l][l];
double p_new[l];
double p_old[l];

int main() {
    seedMT(4357U);
//    run();
        timeit(run);
//    int i;
//    double betas[4] = {0.01, 0.35, 0.5};
//    for(i = 0; i < 3; i++)
//    {
//        beta = betas[i];
//        run();
//    }
    return 0;
}

void run(){
    int i;
    setup();
//    char path[300] = "/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS/metcompc/Trabalho2/COPSquare/COPSquare/";
//    char name[50];
//    sprintf(name, "copSquare%.2f.txt", beta);
//    strcat(path, name);
//    FILE *file = fopen(path, "w");
//    report(file);
    char path[] = "/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS/metcompc/Trabalho2/COPEquilibrium/COPEquilibrium/copEquilibrium.txt";
    FILE *file = fopen(path, "w");
    for(i = 0; i < iterations; i++)
    {
        iterate();
        if(i == iterations - 1)
            printLattice();
//            report(file);
        measureEquilibrium(i, file);
#if PRINT_LATTICE
        //        if(i == iterations-1)
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
    int deltaEs[7] = {-12, -8, -4, 0, 4, 8, 12};
    for(i = 0; i < 7; i++)
        exps[offset + deltaEs[i]] = exp(-beta*deltaEs[i]);
    
    for(i = 0; i < l; i++)
        for(j = 0; j < l; j++)
            lattice[i][j] = -1;
    
    double off = 25;
    for(i = off; i < l-off; i++)
        for(j = off; j < l-off; j++)
            lattice[i][j] = 1;
    
    //
    int sum;
    for(j = 0; j < l; j++)
    {
        sum = 0;
        for(i = 0; i < l; i++)
        {
            if(lattice[i][j] == 1)
                sum++;
        }
        p_old[j] = (sum*1.)/(1.*l);
    }

}

//your neighbours are (top,right,bottom,left) = (i,j+1) (i+1,j) (i,j-1) (i-1,j)
void iterate(){
    int i, j;
    for(i = 0; i < l; i++)
        for(j = 0; j < l; j++){
            double r = randomDoubleInclusive(0., 1.);
            
            if(r < 0.25) //top
                attemptFlip(i, j, i, getPBCIndex(j+1));
            else if(r < 0.5) //bottom
                attemptFlip(i, j, i, getPBCIndex(j-1));
            else if(r < 0.75) //left
                attemptFlip(i, j, getPBCIndex(i-1), j);
            else //right
                attemptFlip(i, j, getPBCIndex(i+1), j);
        }
}


//in = neighbour's i, jn = neighbour's j
//all indices are expected to be valid, that is, corrected for boundary conditions, before coming to this function
void attemptFlip(int i, int j, int in, int jn){
    if(lattice[i][j] == lattice[in][jn])
        return;
    int deltaE = computePairEnergy(i, j, in, jn);
    //    printf("%d\n", deltaE);
    //    printf("exps[12 + %d] = %.4f\n", deltaE, exps[offset + deltaE]);
    double r = randomDoubleInclusive(0.,1.);
    if(deltaE < 0 || exps[offset + deltaE] >= r){
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

//in = neighbour's i, jn = neighbour's j
int computeEnergy(int i, int j, int in, int jn){
    int energy = 0;
    //compute the energy of the neighbours of (i,j) excluding (in,jn)
    if(in != i || jn != getPBCIndex(j+1))
    {
        if(j+1 < l)
            energy += lattice[i][j+1];
        else
            energy += lattice[i][0];
    }
    if(in != i || jn != getPBCIndex(j-1))
    {
        if(j-1 >= 0)
            energy += lattice[i][j-1];
        else
            energy += lattice[i][0];
    }
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
        return l-1;
    if(i > (l-1))
        return 0;
    return i;
}

void measureEquilibrium(int step, FILE *file){
    int i, j, sum;
    
    for(j = 0; j < l; j++)
        p_old[j] = p_new[j];
    
    for(j = 0; j < l; j++)
    {
        sum = 0;
        for(i = 0; i < l; i++)
            if(lattice[i][j] == 1)
                sum++;
        p_new[j] = (1.*sum)/(1.*l);
    }
    double diff = 0.;
    for(j = 0; j < l; j++)
        diff += fabs(p_new[j]-p_old[j]);
    fprintf(file, "%d %.6f\n", step, diff/l);
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


void seedMT(uint32 seed)
{
    register uint32 x = (seed | 1U) & 0xFFFFFFFFU, *s = state;
    register int    j;
    
    for(left=0, *s++=x, j=N; --j;
        *s++ = (x*=69069U) & 0xFFFFFFFFU);
}


uint32 reloadMT(void)
{
    register uint32 *p0=state, *p2=state+2, *pM=state+M, s0, s1;
    register int    j;
    
    if(left < -1)
        seedMT(4357U);
    
    left=N-1, next=state+1;
    
    for(s0=state[0], s1=state[1], j=N-M+1; --j; s0=s1, s1=*p2++)
        *p0++ = *pM++ ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);
    
    for(pM=state, j=M; --j; s0=s1, s1=*p2++)
        *p0++ = *pM++ ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);
    
    s1=state[0], *p0 = *pM ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);
    s1 ^= (s1 >> 11);
    s1 ^= (s1 <<  7) & 0x9D2C5680U;
    s1 ^= (s1 << 15) & 0xEFC60000U;
    return(s1 ^ (s1 >> 18));
}


uint32 randomMT(void)
{
    uint32 y;
    
    if(--left < 0)
        return(reloadMT());
    
    y  = *next++;
    y ^= (y >> 11);
    y ^= (y <<  7) & 0x9D2C5680U;
    y ^= (y << 15) & 0xEFC60000U;
    return(y ^ (y >> 18));
}

double randomDoubleInclusive(double lowerLimit, double upperLimit){
    unsigned long randomNumber = (unsigned long) randomMT();
    unsigned long maxRandomNumber = pow(2,32) - 1;
    return lowerLimit + (double)randomNumber / (double)maxRandomNumber * (upperLimit - lowerLimit);
}
