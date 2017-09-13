#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

struct Disk{
    float x;
    float y;
};

float randomPosition();
int placeDisks();
float meanOfNPlacements();
const int tolerance = 1000;
const int nruns = 1000;
const int nplacements = 1000;
const int maxDisks = 30; //guess
const float L = 50;
const float R = 1;
const float D = 2*R;

float randomPosition(){
    return 1 + rand() % 48;
}

int placeDisks() {
    struct Disk disks[maxDisks];
    int attempts = 0;
    int disk_counter = 0;
    while(attempts < tolerance) {
        float x = randomPosition();
        float y = randomPosition();
        int overlapping = 0;
        int i = 0;
        for(; i < disk_counter; i++){
            float x_sup = fabs(disks[i].x - x);
            float y_sup = fabs(disks[i].y - y);
            if(x_sup < D || y_sup < D){
                overlapping = 1;
                break;
            }
        }
        if(overlapping){
            attempts++;
        }
        else {
            disks[disk_counter].x = x;
            disks[disk_counter].y = y;
            disk_counter++;
        }
    }
    return disk_counter;
}

float meanOfNPlacements(){
    int i = 0;
    float sum = 0;
    for(; i < nplacements; i++){
        int disksPlaced = placeDisks();
        sum += disksPlaced;
    }
    return sum/nplacements;
}

int main() {
    srand((unsigned int)time(NULL));
    //int i = 0;
    //float means[nmeans] = {};
    //for(; i < nmeans; i++){
    //    means[i] = run();
    //}
    int i = 0;
    FILE *file = fopen("/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/MÉTODOS COMPUTACIONAIS C/MetCompC/placingcircles.dat", "w");
    for(; i < nruns; i++){
        fprintf(file, "%.4f\n", meanOfNPlacements());
    }
    fclose(file);
    return 0;
}










