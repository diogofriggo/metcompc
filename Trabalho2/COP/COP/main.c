#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>

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
void printBoard();

const double rho = 0.5;
const int l = 50;
const int l2 = l/2;
const int iterations = 100;
const int J = 1;
const int kB = 1;
const int T = 1;
//beta >= 187 forbids any transition
//changing boundary conditions beta >= 94 forbids any transition
int beta = 530000000;//1/(kB*T);

int lattice[l][l];

int main() {
    srand((unsigned int)time(NULL));
//    stats();
    int i;
    for(i = 1; i < 550000000; i += 1000000)
    {
        beta = i;
        run();
        //timeit(run);
        int beta_stable = 1;
        int j;
        for(j = 0; j < l; j++)
            if(lattice[1][j] == -1){
                beta_stable = 0;
                break;
            }
        if(beta_stable)
            printf("%.2f is stable\n", beta);
    }
    
//    int i;
//    for(i = 0; i < 100; i++){
//        printf("%.2f ", randomDoubleInclusive(0.,1.));
//    }
    return 0;
}

void run(){
    int i;
    setup();
    const char path[] = "/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS/metcompc/Trabalho2/COP/COP/cop.txt";
    FILE *file = fopen(path, "w");
    for(i = 0; i < iterations; i++){
        iterate();
        report(file);
//        if(i%100==0)
            //printBoard();
    }
    fclose(file);
}

//void stats()
//{
//    double m = pow(1. - pow(1./sinh(2*beta*J), 2.), 1./8.);
//    double pplus = (1.+m)/2.;
//    double pminus = (1.-m)/2.;
//    printf("m = %.2f, %.2f <= p <= %.2f\n", m, pminus, pplus);
//}

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

//in = neighbour's i, jn = neighbour's j
void attemptFlip(int i, int j, int in, int jn){
    int deltaE = computePairEnergy(i, j, in, jn);
//    printf("%d\n", deltaE);
    if(deltaE < 0 || exp(-beta*deltaE) < randomDoubleInclusive(0.,1.)){
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
    //compute the energy of the neighbours of (i,j) excluding (in,jn)
    //I think in != i && jn != (j+1) would also work
    if(!(in == i && jn == (j+1))) //top
    {
        if((j+1) < (l-1)) //if not neighbour on the top layer
            energy += lattice[i][j+1];
        else
            energy += -1; //boundary condition which keeps interface from wandering off
    }
    if(!(in == (i+1) && jn == j)) //right
    {
        if((i+1) < l)
            energy += lattice[i+1][j];
        else
            energy += lattice[0][j]; //horizontal periodic boundary condition
    }
    if(!(in == i && jn == (j-1))) //bottom
    {
        if((j-1) > 0) //if not neighbour on the bottom layer
            energy += lattice[i][j-1];
        else
            energy += 1; //boundary condition which keeps interface from wandering off
    }
    if(!(in == (i-1) && jn == j)) //left
    {
        if((i-1) > 0)
            energy += lattice[i-1][j];
        else
            energy += lattice[l-1][j]; //horizontal periodic boundary condition
    }
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

void printBoard(){
    printf("Printing board...\n");
    int i, j;
    for(j = l-1; j >= 0; j--){
        printf("%2d ", j);
        for(i = 0; i < l; i++)
            printf("%s", lattice[i][j] == 1 ? "u" : "-");
        printf("\n");
    }
    printf("\n");
    //I avoided this clever trick because it hinders readability:
    //printf("%3d%s", board[i][j], j%ncols == (ncols-1) ? "\n" : "");
}
