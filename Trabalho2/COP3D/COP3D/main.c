#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#define ATTEMPT_FLIP 0

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
void attemptFlip(int i, int j, int k, int in, int jn, int kn);
int computePairEnergy(int i, int j, int k, int in, int jn, int kn);
int computeEnergy(int i, int j, int k, int in, int jn, int kn);
void report(FILE *file);
void setup();
void timeit(void (*fun)(void));
double randomDoubleInclusive(double lowerLimit, double upperLimit);

const double rho = 0.5;
const int l = 50;
//const int l2 = l/2;
const int iterations = 250;
const int J = 1;
const int kB = 1;
const int T = 1;

double beta = 0.2 ;//1/(kB*T);
const int offset = 20;
double exps[41];
int lattice[l][l][l];
int getPBCIndex(int i);

int main() {
    seedMT(4357U);
//    timeit(run);
    int i;
    double betas[4] = {0.01, 0.2, 0.5};
    for(i = 0; i < 3; i++)
    {
        beta = betas[i];
        run();
    }
    return 0;
}

void run(){
    int i;
    setup();
    char path[300] = "/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS/metcompc/Trabalho2/COP3D/COP3D/";
    char name[50];
    sprintf(name, "cop3D%.2f.txt", beta);
    strcat(path, name);
    FILE *file = fopen(path, "w");
    report(file);
    for(i = 0; i < iterations; i++){
        iterate();
        if(i % 5 == 0)
            report(file);
    }
    fclose(file);
}

void setup()
{
    int i, j, k;
    //precompute possible exponentials to save CPU time
    for(i = 0; i < 41; i++)
        exps[i] = 0.;
    int deltaEs[11] = {-20, -16, -12, -8, -4, 0, 4, 8, 12, 16, 20};
    for(i = 0; i < 11; i++)
        exps[offset + deltaEs[i]] = exp(-beta*deltaEs[i]);
    
    for(i = 0; i < l; i++)
        for(j = 0; j < l; j++)
            for(k = 0; k < l; k++)
            lattice[i][j][k] = -1;
    
    double p = 20;
    for(i = p; i < l-p; i++)
        for(j = p; j < l-p; j++)
            for(k = p; k < l-p; k++)
                lattice[i][j][k] = 1;
}

//your neighbours are (top,right,bottom,left) = (i,j+1) (i+1,j) (i,j-1) (i-1,j)
void iterate(){
    int i, j, k;
    for(i = 0; i < l; i++)
        for(j = 0; j < l; j++)
            for(k = 0; k < l; k++){
                double r = randomDoubleInclusive(0., 1.);
                
                if(r < 0.1666)
                    attemptFlip(i, j, k, getPBCIndex(i+1), j, k);
                else if(r < 0.3333)
                    attemptFlip(i, j, k, getPBCIndex(i-1), j, k);
                else if(r < 0.5)
                    attemptFlip(i, j, k, i, getPBCIndex(j+1), k);
                else if(r < 0.6666)
                    attemptFlip(i, j, k, i, getPBCIndex(j-1), k);
                else if(r < 0.8333)
                    attemptFlip(i, j, k, i, j, getPBCIndex(k+1));
                else
                    attemptFlip(i, j, k, i, j, getPBCIndex(k-1));
            }
}


//in = neighbour's i, jn = neighbour's j
//all indices are expected to be valid, that is, corrected for boundary conditions, before coming to this function
void attemptFlip(int i, int j, int k, int in, int jn, int kn){
    if(lattice[i][j][k] == lattice[in][jn][kn])
        return;
    int deltaE = computePairEnergy(i, j, k, in, jn, kn);
    //    printf("%d\n", deltaE);
    //    printf("exps[12 + %d] = %.4f\n", deltaE, exps[offset + deltaE]);
    double r = randomDoubleInclusive(0.,1.);
    if(deltaE < 0 || exps[offset + deltaE] >= r){
        int temp = lattice[i][j][k];
        lattice[i][j][k] = lattice[in][jn][kn];
        lattice[in][jn][kn] = temp;
#if ATTEMPT_FLIP
        printf("success flip. ");
#endif
    }
#if ATTEMPT_FLIP
    else
        printf("failed flip. ");
    printf("%d < 0 || %.4f >= %.4f) (%d,%d,%d,%d) <-> (%d,%d,%d,%d)\n", deltaE, exp(-beta*deltaE), r, i, j, k, lattice[i][j][k], in, jn, kn, lattice[in][jn][kn]);
    if(deltaE > 1000 || deltaE < -1000)
        printf("ERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\n");
#endif
}

int computePairEnergy(int i, int j, int k, int in, int jn, int kn){
    return 2*J*(lattice[i][j][k]*computeEnergy(i, j, k, in, jn, kn) + lattice[in][jn][kn]*computeEnergy(in, jn, kn, i, j, k));
}

//in = neighbour's i, jn = neighbour's j, kn = neighbour's kn
int computeEnergy(int i, int j, int k, int in, int jn, int kn){
    int energy = 0;
    //compute the energy of the neighbours of (i,j,k) excluding (in,jn,kn)
    if(in != getPBCIndex(i+1) || jn != j || kn != k)
    {
        if((i+1) < l)
            energy += lattice[i+1][j][k];
        else
            energy += lattice[0][j][k];
    }
    if(in != getPBCIndex(i-1) || jn != j || kn != k)
    {
        if((i-1) >= 0)
            energy += lattice[i-1][j][k];
        else
            energy += lattice[l-1][j][k];
    }
    if(in != i || jn != getPBCIndex(j+1) || kn != k)
    {
        if(j+1 < l)
            energy += lattice[i][j+1][k];
        else
            energy += lattice[i][0][k];
    }
    if(in != i || jn != getPBCIndex(j-1) || kn != k)
    {
        if(j-1 >= 0)
            energy += lattice[i][j-1][k];
        else
            energy += lattice[i][0][k];
    }
    if(in != i || jn != j || kn != getPBCIndex(k+1))
    {
        if((k+1) < l)
            energy += lattice[i][j][k+1];
        else
            energy += lattice[i][j][0];
    }
    if(in != i || jn != j || kn != getPBCIndex(k-1))
    {
        if((k-1) >= 0)
            energy += lattice[i][j][k-1];
        else
            energy += lattice[i][j][l-1];
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

void report(FILE *file){
    int i, j, k;
    for(i = 0; i < l; i++)
        for(j = 0; j < l; j++)
            for(k = 0; k < l; k++)
                fprintf(file, "%d %d %d %d\n", i, j, k, lattice[i][j][k]);
}

void timeit(void (*fun)(void)){
    clock_t begin = clock();
    fun();
    clock_t end = clock();
    double timeElapsed = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("%.6f\n", timeElapsed);
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
