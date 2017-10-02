//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include <string.h>
//
//void nonStationaryDiffusion();
//void describeLeftMatrix(double k, int nXs, double a[nXs-1], double b[nXs], double c[nXs-1], double d[nXs]);
//void describeRightMatrix(double k, int nXs, double a[nXs-1], double b[nXs], double c[nXs-1], double d[nXs]);
//void runNonStationaryDiffusion(double *a, double *b, double *c, double *d, double *x, int nXs, int nTs, double k, double x0, double xl);
//void reportNonStationaryDiffusion(double *x, int nXs, double dx, double k, double maxT);
//void thomasSolveTridiagonalMatrix(int n, double a[n-2], double b[n], double c[n-1], double d[n], double x[n]);
//void once();
//void copy(int n, double *d, double *x);
//
//int main(int argc, const char * argv[]) {
//    once();
//    //nonStationaryDiffusion();
//    return 0;
//}
//
//void once(){
//    int i, L = 100, t;
//    double k = 0.5;
//    double dt = 0.1, dx = 1;
//    //if time is divided in intervals of 0.1 between 0 and 100 then nTs = [0, 0.1, ... , 999.9] -> 1000 steps
//    //if space is divided in intervals of 0.1 between 0 and 100 then nXs = [0, 0.1, ... , 999.9, 1000] -> 1001 positions
//    int nTs, nXs = (int)round(L/dx) + 1;
//    
//    double x0 = 1, xl = 0;
//    
//    double *leftA, *leftB, *leftC, *leftD, *leftX;
//    double *rightA, *rightB, *rightC, *rightD, *rightX;
//    
//    leftA = (double*)calloc(nXs-1, sizeof(double));
//    leftB = (double*)calloc(nXs, sizeof(double));
//    leftC = (double*)calloc(nXs-1, sizeof(double));
//    leftD = (double*)calloc(nXs, sizeof(double));
//    leftX = (double*)calloc(nXs, sizeof(double));
//    
//    rightA = (double*)calloc(nXs-1, sizeof(double));
//    rightB = (double*)calloc(nXs, sizeof(double));
//    rightC = (double*)calloc(nXs-1, sizeof(double));
//    rightD = (double*)calloc(nXs, sizeof(double));
//    rightX = (double*)calloc(nXs, sizeof(double));
//    
//    describeLeftMatrix(k, nXs, leftA, leftB, leftC, leftD);
//    describeRightMatrix(k, nXs, rightA, rightB, rightC, rightD);
//    
//    nTs = (int)round(1000/dt);
//    for(t = 0; t < nTs; t++){
//        thomasSolveTridiagonalMatrix(nXs, rightA, rightB, rightC, rightD, rightX);
//        
//        //should I do it here?
//        copy(nXs, rightD, leftX);
//        d[0] = x[0] + k*x0;
//        d[nXs-1] = x[nXs-1] + k*xl;
//
//        
//        copy(nXs, leftD, rightX);
//        //leftD = rightX;
//
//        thomasSolveTridiagonalMatrix(nXs, leftA, leftB, leftC, leftD, leftX);
//
//        copy(nXs, rightD, leftX);
//        //rightD = leftX;
//        if(t == 0)
//            for(i = 0; i < nXs; i++)
//                printf("%d: %.4f%s", i, leftD[i], i % 10 == 0 ? "\n" : " ");
//    }
//    
//    free(leftA);
//    free(leftB);
//    free(leftC);
//    free(leftD);
//    free(leftX);
//    free(rightA);
//    free(rightB);
//    free(rightC);
//    //free(rightD);
//    //free(rightX);
//}
//
//void copy(int n, double *d, double *x){
//    int i;
//    for(i = 0; i < n; i++)
//        d[i] = x[i];
//}
//
//void nonStationaryDiffusion(){
//    int k, i, L = 100;
//    double dt = 0.1, dx = 1;
//    //if time is divided in intervals of 0.1 between 0 and 100 then nTs = [0, 0.1, ... , 999.9] -> 1000 steps
//    //if space is divided in intervals of 0.1 between 0 and 100 then nXs = [0, 0.1, ... , 999.9, 1000] -> 1001 positions
//    int nTs, nXs = (int)round(L/dx) + 1;
//    const int maxTsLength = 6, ksLength = 10;
//    double maxTs[maxTsLength], ks[ksLength];
//    
//    double *a, *b, *c, *d, *x;
//    a = (double*)calloc(nXs-1, sizeof(double));
//    b = (double*)calloc(nXs, sizeof(double));
//    c = (double*)calloc(nXs-1, sizeof(double));
//    d = (double*)calloc(nXs, sizeof(double));
//    x = (double*)calloc(nXs, sizeof(double));
//    
//    for(k = 0; k < ksLength; k++)
//        ks[k] = 0.1*(k+1);
//    
//    maxTs[0] = 2;
//    maxTs[1] = 10;
//    maxTs[2] = 50;
//    maxTs[3] = 200;
//    maxTs[4] = 1000;
//    maxTs[5] = 2000;
//    
//    for(k = 0; k < ksLength; k++){
//        for(i = 0; i < maxTsLength; i++){
//            nTs = (int)round(maxTs[i]/dt);
//            runNonStationaryDiffusion(a, b, c, d, x, nXs, nTs, ks[k]);
//            //reportNonStationaryDiffusion(x, nXs, dx, ks[k], maxTs[i]);
//        }
//    }
//    
//    free(a);
//    free(b);
//    free(c);
//    free(d);
//    free(x);
//}
//
//void reportNonStationaryDiffusion(double *x, int nXs, double dx, double k, double maxT){
//    char path[200] = "/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS/metcompc/Aula4_1209/Results/";
//    char name[50];
//    sprintf(name, "CCrankDiffusionK%.1ft%.0f.txt", k, maxT);
//    strcat(path, name);
//    FILE *file = fopen(path, "w");
//    int j;
//    for(j = 0; j < nXs; j++)
//        fprintf(file, "%.4f %.4f\n", j*dx, x[j]);
//    fclose(file);
//}
//
//void runNonStationaryDiffusion(double *a, double *b, double *c, double *d, double *x, int nXs, int nTs, double k, double x0, double xl){
//    //describeNonStationaryDiffusionBoundaryConditions(k, nXs, a, b, c, d);
//    int t, i;
//    for(t = 0; t < nTs; t++){
//        thomasSolveTridiagonalMatrix(nXs, a, b, c, d, x);
//        d = x;
//    }
//}
//
//void describeLeftMatrix(double k, int nXs, double a[nXs-1], double b[nXs], double c[nXs-1], double d[nXs]){
//    int i;
//    for(i = 0; i < nXs-2; i++)
//        a[i] = -k;
//    for(i = 1; i < nXs-1; i++)
//        b[i] = 2+2*k;
//    for(i = 1; i < nXs-1; i++)
//        c[i] = -k;
//    for(i = 1; i < nXs-1; i++)
//        d[i] = 0;
//    
//    //specifying left boundary condition
//    b[0] = 1;
//    c[0] = 0;
//    d[0] = 1;
//    
//    //specifying right boundary condition
//    a[nXs-2] = 0;
//    b[nXs-1] = 1;
//    d[nXs-1] = 0;
//}
//
//void describeRightMatrix(double k, int nXs, double a[nXs-1], double b[nXs], double c[nXs-1], double d[nXs]){
//    int i;
//    for(i = 0; i < nXs-2; i++)
//        a[i] = k;
//    for(i = 1; i < nXs-1; i++)
//        b[i] = 2-2*k;
//    for(i = 1; i < nXs-1; i++)
//        c[i] = k;
//    for(i = 1; i < nXs-1; i++)
//        d[i] = 0;
//    
//    //specifying left boundary condition
//    b[0] = 1;
//    c[0] = 0;
//    d[0] = 0.5; //WARNING: CHANGED FROM 1 TO 0.5 1809
//    
//    //specifying right boundary condition
//    a[nXs-2] = 0;
//    b[nXs-1] = 1;
//    d[nXs-1] = 0;
//}
//
////it writes the result to the x[n] vector
//void thomasSolveTridiagonalMatrix(int n, double a[n-1], double b[n], double c[n-1], double d[n], double x[n]){
//    int i;
//    double denominator;
//    double c_[n];
//    double d_[n];
//    
//    //1 - compute coefficients
//    c_[0] = c[0]/b[0];
//    d_[0] = d[0]/b[0];
//    for(i = 1; i < (n-1); i++){
//        //there's a very subtle issue when manipulating a because it is one index behind
//        //so to instead of a[i] as prescribed by thomas' formula we must use a[i-1]
//        denominator = b[i]-c_[i-1]*a[i-1];
//        c_[i] = c[i]/denominator;
//        d_[i] = (d[i]-d_[i-1]*a[i-1])/denominator;
//    }
//    
//    //d[] has one further step:
//    denominator = b[n-1]-c_[n-2]*a[n-2];
//    d_[n-1] = (d[n-1]-d_[n-2]*a[n-2])/denominator;
//    
//    //2 - substitute back to find solution vector
//    x[n-1] = d_[n-1];
//    for(i = n-2; i >= 0; i--)
//        x[i] = d_[i]-c_[i]*x[i+1];
//}
//

#include <time.h>

#include<stdio.h>

void swap(int* a, int* b)
{
    int t = *a;
    *a = *b;
    *b = t;
}

int partition (int arr[], int low, int high)
{
    int pivot = arr[high];
    int i = (low - 1);
    
    for (int j = low; j <= high- 1; j++)
    {
        if (arr[j] <= pivot)
        {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}

void quickSort(int arr[], int low, int high)
{
    if (low < high)
    {
        int pi = partition(arr, low, high);
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}

void printArray(int arr[], int size)
{
    int i;
    for (i=0; i < size; i++)
        printf("%d ", arr[i]);
    printf("n");
}

int main()
{
    clock_t begin = clock();
    const int sizee = 100000;
    int i;
    int arr[sizee];// = {10, 7, 8, 9, 1, 5};
    for(i = 0; i < sizee; i++)
        arr[i] = sizee-i;
    int n = sizeof(arr)/sizeof(arr[0]);
    quickSort(arr, 0, n-1);
    //printf("Sorted array: n");
    //printArray(arr, n);
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("%.3f\n", time_spent);
    return 0;
}
