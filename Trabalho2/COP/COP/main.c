#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>

void run();
void setup();
void iterate();
void report(FILE *file);
void setup();
void timeit(void (*fun)(void));
int randomIntegerInclusive(int lowerLimit, int upperLimit);

const int l = 50;
const int l2 = l/2;
const int iterations = 1000;
int lattice[l][l];

int main() {
    srand((unsigned int)time(NULL));
    timeit(run);
    return 0;
}

void run(){
    int i;
    
    setup();
    const char path[] = "/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS/metcompc/Trabalho2/COP/COP/cop.txt";
    FILE *file = fopen(path, "w");
    for(i = 0; i < iterations; i++){
        iterate();
        if(i%10==0)
            report(file);
    }
    fclose(file);
}

void setup()
{
    int i, j;
    
    for(i = 0; i < l; i++)
        for(j = 0; j < l; j++)
            lattice[i][j] = 0;
    
    for(i = 0; i < l; i++)
        for(j = 0; j < l2; j++)
            lattice[i][j] = 1;
}

void iterate(){
    for(i = 0; i < l; i++)
        for(j = 0; j < l; j++)
            
}

void calculateEnergy(int x, int y){
    
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

int randomIntegerInclusive(int lowerLimit, int upperLimit){
    return lowerLimit + rand() % (upperLimit+1 - lowerLimit);
}