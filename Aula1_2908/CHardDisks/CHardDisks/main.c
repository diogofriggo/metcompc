#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

struct Disk{
    double x, y;
};

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
struct Disk disks[maxDisks];

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
    numberOfDisksPlaced();
    positionsOfDisksPlaced();
}

void numberOfDisksPlaced(){
    int i;
    int disksPlaced[nplacements];
    
    for(i = 0; i < nplacements; i++)
        disksPlaced[i] = placeDisks();
    
    FILE *file = fopen("/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS COMPUTACIONAIS C/metcompc/Aula1_2908/Results/CNumberOfDisksPlaced.txt", "w");
    for(i = 0; i < nplacements; i++)
        //printf("%d%s", disksPlaced[i], i % 50 == 0 ? "\n" : " ");
        fprintf(file, "%d\n", disksPlaced[i]);
    fclose(file);
}

void positionsOfDisksPlaced(){
    char dir[200] = "/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS COMPUTACIONAIS C/metcompc/Aula1_2908/Results/";
    char path[250];
    char name[50];
    FILE *file;
    int i, j, disksPlaced;
    for(i = 0; i < 4; i++){
        disksPlaced = placeDisks();
        strcpy(path, dir);
        sprintf(name, "CPositionsOfDisksPlaced%d.txt", i+1);
        strcat(path, name);
        file = fopen(path, "w");
        for(j = 0; j < disksPlaced; j++)
            fprintf(file, "%.4f %.4f\n", disks[j].x, disks[j].y);
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
            xDiff = fabs(disks[i].x - x);
            yDiff = fabs(disks[i].y - y);
            if(pow(xDiff,2)+pow(yDiff,2) <= D){
                overlapping = 1;
                break;
            }
        }
        if(overlapping){
            attempts++;
        }
        else {
            disks[diskCounter].x = x;
            disks[diskCounter].y = y;
            diskCounter++;
        }
    }
    return diskCounter;
}

double randomPosition(){
    return 1 + rand() % 48;
}