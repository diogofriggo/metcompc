#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double randomPosition();
int placeDisks();
void numberOfDisksPlaced();
void positionsOfDisksPlaced();
void run();
void timeit(void (*fun)(void));

const int tolerance = 1000;
const int nplacements = 10000;
const int maxDisks = 600; //guess
const double L = 50;
const double R = 1;
const double D = 2*R;
double (*disks)[maxDisks];

//By removing structs and using calloc I got not significant improvement in performance at all!
//Running on 1000 tolerance, 10000 placements I got from 241s to 246s !!!
int main() {
    srand((unsigned int)time(NULL));
    timeit(run);
    return 0;
}

void timeit(void (*fun)(void)){
    clock_t begin = clock();
    fun();
    clock_t end = clock();
    double timeSpent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("%2.f\n", timeSpent);
}

void run(){
    disks = calloc(maxDisks*2, sizeof(double));
    numberOfDisksPlaced();
    positionsOfDisksPlaced();
    free(disks);
}

void numberOfDisksPlaced(){
    int i;
    FILE *file = fopen("/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS COMPUTACIONAIS C/metcompc/Aula1_2908/Results/CNumberOfDisksPlacedFaster.txt", "w");
    for(i = 0; i < nplacements; i++)
        fprintf(file, "%d\n", placeDisks());
    fclose(file);
}

void positionsOfDisksPlaced(){
    char dir[200] = "/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS COMPUTACIONAIS C/metcompc/Aula1_2908/Results/";
    char path[250];
    char name[50];
    FILE *file;
    int k, j, disksPlaced;
    for(k = 0; k < 4; k++){
        disksPlaced = placeDisks();
        strcpy(path, dir);
        sprintf(name, "CPositionsOfDisksPlaced%dFaster.txt", k+1);
        strcat(path, name);
        file = fopen(path, "w");
        for(j = 0; j < disksPlaced; j++)
            fprintf(file, "%.4f %.4f\n", disks[j][0], disks[j][1]);
        fclose(file);
    }
}

int placeDisks() {
    double x, y, xDiff, yDiff;
    register int attempts = 0, diskCounter = 0;
    int i, overlapping = 0;
    
    while(attempts < tolerance) {
        x = randomPosition();
        y = randomPosition();
        overlapping = 0;
        for(i = 0; i < diskCounter; i++){
            xDiff = fabs(disks[i][0] - x);
            yDiff = fabs(disks[i][1] - y);
            if(sqrt(pow(xDiff,2)+pow(yDiff,2)) <= D){
                overlapping = 1;
                break;
            }
        }
        if(overlapping){
            attempts++;
        }
        else {
            disks[diskCounter][0] = x;
            disks[diskCounter][1] = y;
            diskCounter++;
        }
    }
    return diskCounter;
}

double randomPosition(){
    return 1 + rand() % 48;
}