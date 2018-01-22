/**********************
 * important reminders *
 ***********************/

//players[id].offer instead of assigning to struct
//avoid players[id-1] by ignoring the first element of players

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>

/***************
 * DEFINITIONS *
 ***************/

#define DIAGNOSTICS 0

struct Player{
    int id;
    float offer;
    float accept;
    float current;
    float previous;
    float earned;
};

int sub2ind(int i, int j);
int sub2ind_generic(int i, int j, int jmax);
void ind2sub(int ind, int *i, int *j);
void ind2sub_generic(int ind, int *i, int *j, int jmax);
void test_sub2ind_generic();
void test_ind2sub();
int randomIntegerInclusive();
int randomIntegerLeftInclusive();
void test_randomIntegerInclusive();
void test_randomIntegerLeftInclusive();
void test_readingMatrixWithPointer();
void test_writingMatrixWithPointer();
int getBoardValueByIndex(int index);
void setBoardValueByIndex(int index, int value);
int randomBoardIndex();
void setUpBoard();
void printBoard();
void setUpPlayers();
void printPlayers();
struct Player getPlayer(int playerId);
void setPlayer(int playerId, struct Player player);
float randomFloatInclusive(float lowerLimit, float upperLimit);
void test_randomFloatInclusive(float lowerLimit, float upperLimit);
void boundsTest_randomFloatInclusive(float lowerLimit, float upperLimit);
void printIndexedBoard();
void buildFirstNeighbours();
void test_buildFirstNeighbours();
void test_WhetherChangingAStructChangesTheCollectionItCameFrom();
void buildFirstNeighbours();
int fixForPeriodicBoundary(int index, int minValue, int maxValue);
void moveToEmptyPosition(int playerId, int i, int j);
void iterateThroughBoard(void (*apply)(int, int, int));
void playWithNeighboursOf(int playerId, int i, int j);
void play(int playerId, int neighbourId);
void tryMovingToEmptyPosition(int playerId, int i, int j);
void runDiagnosticsMode();
void replicate(int playerId, int i, int j);
void setup();
void run();
void timeit(void (*fun)(void));
void test_flooringFloat();
float roundFloat(float number, int decimalPlaces);
void report();
void testDynamicMatrixAllocation();
void testDynamicMatrixAllocationWithDoublePointer();

const int topLeft = 0, top = 1, topRight = 2;
const int left = 3, right = 4;
const int bottomLeft = 5, bottom = 6, bottomRight = 7;
const int nrows = 1000;
const int ncols = 1000;
const int boardSize = nrows*ncols;
const int maximumBoardSize = 1000000;
int board[nrows][ncols];
const int nplayers = 1000;
const int maximumNumberOfPlayers = 100000;
struct Player players[nplayers];
int neighbours[boardSize][8][2];
int maximumAmount = 1000;
//int iterations = 1000;
int maximumNumberOfIterations = 100;

/******************
 * IMPLEMENTATION *
 *****************/

//10.000 cells, 1.000 players
//100 iterations 0.06606s
//1.000 iterations 0.6214s
//10.000 iterations 6.0889s
//100.000 iterations 61.03257s
//1.000.000 iterations 606.070643s ~ 10min

int main() {
    srand((unsigned int)time(NULL));
    //runDiagnosticsMode();
    //report();
    //timeit(run);
    //test_flooringFloat();
    //testDynamicMatrixAllocation();
    testDynamicMatrixAllocationWithDoublePointer();
    return 0;
}

int *boardTest1;

void testDynamicMatrixAllocation(){
    const size_t nrows = 10;
    const size_t ncols = 10;
    boardTest1 = (int *)malloc(sizeof(int) * nrows*ncols);
    int i, j;
    for(i = 0; i < nrows; i++)
        for(j = 0; j < ncols; j++)
            *(boardTest1 + i*ncols + j) = sub2ind_generic(i, j, ncols);
    for(i = 0; i < nrows; i++)
        for(j = 0; j < ncols; j++){
            printf("%4d ", *(boardTest1 + i*ncols + j));
            if(j == ncols-1)
                printf("\n");
        }
    
    free(boardTest1);
}

//int **boardTest2;

void testDynamicMatrixAllocationWithDoublePointer(){
    int nrows = 10, ncolumns = 10;
    int i, j;
    
    int **array2 = malloc(nrows * sizeof(int *));
    for(i = 0; i < nrows; i++)
        array2[i] = malloc(ncolumns * sizeof(int));
    
    for(i = 0; i < nrows; i++){
        for(j = 0; j < ncolumns; j++)
            array2[i][j] = 1;
        
        for(i = 0; i < nrows; i++)
            for(j = 0; j < ncolumns; j++){
                printf("%4d ", array2[i][j]);
                if(j == ncolumns-1)
                    printf("\n");
            }
    }
    //I think I should iterate and free also array2[i]
    free(array2);
}

//histograma das estrategias
void report(){
    const char dir[200] = "/Users/diogofriggo/Google Drive/Bolsa/UltimatumReport/";
    int i, b, np, it;
    for(b = 100; b <= maximumBoardSize; b *= 10)
        for(np = 10; np <= maximumNumberOfPlayers; np *= 10){
            //prevent playing with more players than there are cells
            if(np >= b) break;
            for(it = 10; it <= maximumNumberOfIterations; it *= 10){
                run(it);
                char path[300];
                strcpy(path, dir);
                char name[50];
                sprintf(name, "ultimatumStrategiesb=%dnp=%dit=%d.txt", b, np, it);
                strcat(path, name);
                FILE *file = fopen(path, "w");
                for(i = 0; i < nplayers; i++)
                    fprintf(file, "%.2f %.2f\n", players[i].offer, players[i].accept);
                fclose(file);
                printf("Written %s\n", path);
            }
        }
    
}


void timeit(void (*fun)(void)){
    clock_t begin = clock();
    fun();
    clock_t end = clock();
    double timeElapsed = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("%.6f\n", timeElapsed);
}

void run(int iterations){
    setup();
    int i;
    for(i = 0; i < iterations; i++){
        iterateThroughBoard(playWithNeighboursOf);
        iterateThroughBoard(replicate);
        iterateThroughBoard(tryMovingToEmptyPosition);
    }
}

void setup()
{
    setUpBoard();
    setUpPlayers();
    buildFirstNeighbours();
}

void runDiagnosticsMode(){
    test_sub2ind_generic();
    test_randomFloatInclusive(0,1);
    test_randomFloatInclusive(2,4);
    //boundsTest_randomFloatInclusive(0,1);
    test_randomIntegerInclusive();
    test_randomIntegerLeftInclusive();
    test_readingMatrixWithPointer();
    test_writingMatrixWithPointer();
    setUpBoard();
    printBoard();
    setUpPlayers();
    printPlayers();
    printIndexedBoard();
    buildFirstNeighbours();
    test_buildFirstNeighbours();
    test_WhetherChangingAStructChangesTheCollectionItCameFrom();
    printf("\nBEFORE DIFFUSION:\n");
    printBoard();
    //iterateThroughBoard(moveToEmptyPosition);
    iterateThroughBoard(tryMovingToEmptyPosition);
    printf("\nAFTER DIFFUSION:\n");
    printBoard();
    //test_ind2sub();
    iterateThroughBoard(playWithNeighboursOf);
    iterateThroughBoard(replicate);
    test_flooringFloat();
}

//way simpler than the approach before it but may not be desirable
void tryMovingToEmptyPosition(int playerId, int i, int j){
    int neighbourIndex = randomIntegerInclusive(0, 7);
    int playerPosition = sub2ind(i, j);
    int row = neighbours[playerPosition][neighbourIndex][0];
    int col = neighbours[playerPosition][neighbourIndex][1];
    int neighbourId = board[row][col];
    if(!neighbourId){
        board[row][col] = playerId;
        board[i][j] = 0;
    }
#if DIAGNOSTICS
    else{
        printf("Couldn't move %d because %d was on the way\n", playerId, neighbourId);
    }
#endif
}

//should I sort a position and if it's occupied not move?
//that is way easier to implement, but I'd need to implement ind2sub
//this implementation is a bit messy, too many conversions
//implements diffusion
//TODO: improve to move to a radius
void moveToEmptyPosition(int playerId, int i, int j){
    int neighbourIds[8][2];
    int emptyPositions = 0;
    int n, row, col;
    int playerPosition, neighbourId;
    //gather empty positions
    for(n = 0; n < 8; n++){
        playerPosition = sub2ind(i, j);
        row = neighbours[playerPosition][n][0];
        col = neighbours[playerPosition][n][1];
        neighbourId = board[row][col];
        if(!neighbourId){
            neighbourIds[emptyPositions][0] = row;
            neighbourIds[emptyPositions][1] = col;
            emptyPositions++;
        }
    }
    //sort which one to diffuse to
    neighbourId = randomIntegerLeftInclusive(0, emptyPositions);
    //printf("randomIntegerSorted: %d, emptyPositions: %d\n", neighbourId, emptyPositions);
    row = neighbourIds[neighbourId][0];
    col = neighbourIds[neighbourId][1];
#if DIAGNOSTICS
    //printf("%d moving from (%d,%d) to (%d,%d)\n", playerId, i, j, row, col);
#endif
    //perform the actual moving
    board[row][col] = playerId;
    board[i][j] = 0;
}

//applies function to all players of the board
void iterateThroughBoard(void (*apply)(int, int, int)){
    int i, j, playerId;
    //iterate over all positions in the board
    for(i = 0; i < nrows; i++){
        for(j = 0; j < ncols; j++){
            playerId = board[i][j];
            //if there's a player at this position look for neighbours
            if(playerId)
                apply(playerId, i, j);
        }
    }
}

//makes transactions with neighbours. i and j are passed to satisfy generic function definition
void playWithNeighboursOf(int playerId, int i, int j){
    //iterate over all neighbour positions
    int n, row, col, neighbourId;
    int playerPosition = sub2ind(i, j);
    //there are only 8 first neighbours
    for(n = 0; n < 8; n++){
        row = neighbours[playerPosition][n][0];
        col = neighbours[playerPosition][n][1];
        neighbourId = board[row][col];
        //interact with neighbour if present
        if(neighbourId)
            play(playerId, neighbourId);
    }
}

//copies best strategy among its neighbours
//do it before diffusion (naturally)?
//again this struct deal might be slow
void replicate(int playerId, int i, int j){
    int n, row, col, neighbourId;
    int playerPosition = sub2ind(i, j);
    int winnerId = playerId;
    struct Player player, neighbour;
    player = getPlayer(playerId);
    //there are only 8 first neighbours
#if DIAGNOSTICS
    printf("%d: (%.2f,%.2f) -> ", player.id, player.offer, player.accept);
#endif
    for(n = 0; n < 8; n++){
        row = neighbours[playerPosition][n][0];
        col = neighbours[playerPosition][n][1];
        neighbourId = board[row][col];
        //interact with neighbour if present
        if(neighbourId){
            neighbour = getPlayer(neighbourId);
            //copy strategy
            if(neighbour.earned > player.earned){
                player.offer = neighbour.offer;
                player.accept = neighbour.accept;
                player.earned = neighbour.earned;
                winnerId = neighbourId;
            }
        }
    }
#if DIAGNOSTICS
    printf("%d: (%.2f, %.2f)\n", winnerId, player.offer, player.accept);
#endif
    setPlayer(playerId, player);
    
}

//I think assigning to structs may have a performance toll
//on the one hand we avoid accessing the collection some 14 times
//but on the other we have to store two structs in memory
//DOUBT: player offers to give x to neighbour or to take x to himself
void play(int playerId, int neighbourId){
    //does this reflect changes in the players collection?
    struct Player player = getPlayer(playerId);
    struct Player neighbour = getPlayer(neighbourId);
    if(player.offer >= neighbour.accept){
        //sort an amount randomly
        float amount = randomFloatInclusive(0, maximumAmount);
        
        //store what each one has before the transaction
        player.previous      = player.current;
        neighbour.previous   = neighbour.current;
        
        //make the transaction
        player.current      += (1-player.offer)*amount;
        neighbour.current   += player.offer*amount;
        
        //compute earnings for both
        player.earned = player.current-player.previous;
        neighbour.earned = neighbour.current-neighbour.previous;
        
        //this may be really inefficient
        setPlayer(playerId, player);
        setPlayer(neighbourId, neighbour);
        
#if DIAGNOSTICS
        printf("\nDiagnostics of business deals...\n\n");
        printf("The business went down!\n");
        printf("Player %d offered %.2f and Neighbour %d accepted anything >= %.2f\n", player.id, player.offer, neighbour.id, neighbour.accept);
        printf("Player %d had %.2f now has %.2f, earned: %.2f\n", player.id, player.previous, player.current, player.earned);
        printf("Neighbour %d had %.2f now has %.2f, earned: %.2f\n", neighbour.id, neighbour.previous, neighbour.current, neighbour.earned);
    }
    else{
        printf("The business didn't go down because %.2f < %.2f\n", player.offer, neighbour.accept);
#endif
        
    }
}

//this can be generalized to r-th nearest neighbours
void buildFirstNeighbours(){
    int ii, jj, r = 1; //radius
    int i, j, playerPosition, counter;
    //iterate over all positions of the board
    for(i = 0; i < nrows; i++)
        for(j = 0; j < ncols; j++){
            playerPosition = sub2ind(i, j);
            counter = 0;
            //for each position iterate over the 3x3 quadrant centered at the given position
            for(ii = i-r; ii <= i+r; ii++)
                for(jj = j-r; jj <= j+r; jj++){
                    //skip the center position
                    if(!(ii == i && jj == j)){
                        neighbours[playerPosition][counter][0] = fixForPeriodicBoundary(ii, 0, nrows-1);
                        neighbours[playerPosition][counter][1] = fixForPeriodicBoundary(jj, 0, ncols-1);
                        counter++;
                    }
                }
        }
}

int fixForPeriodicBoundary(int value, int minValue, int maxValue)
{
    if(value < minValue) return maxValue;
    if(value > maxValue) return minValue;
    return value;
}

void setUpPlayers(){
    //printf("Setting up players (randomly assigning strategies)...\n");
    int playerId;
    //I'm avoiding player with id 0
    for(playerId = 1; playerId <= nplayers; playerId++){
        struct Player player;
        player.id = playerId;
        player.offer = roundFloat(randomFloatInclusive(0,1),2);
        player.accept = roundFloat(randomFloatInclusive(0,1),2);
        setPlayer(playerId, player);
    }
    //printPlayers();
}

void printPlayers(){
    printf("Printing players...\n");
    int playerId;
    //I'm avoiding player with id 0
    for(playerId = 1; playerId <= nplayers; playerId++){
        struct Player player = getPlayer(playerId);
        printf("%d offers %.2f and accepts %.2f\n", player.id, player.offer, player.accept);
    }
}

//Encapsulating access this way, we avoid forgetting that player_id differs from player index
struct Player getPlayer(int playerId){
    return players[playerId-1];
}

void setPlayer(int playerId, struct Player player){
    players[playerId-1] = player;
}

//this assumes nplayers <= nrows*ncols, otherwise we get an infinite loop
void setUpBoard(){
    //printf("Setting up board (randomly distributing players)...\n");
    int playerId;
    //I'm avoiding player with id 0
    for(playerId = 1; playerId <= nplayers; playerId++){
        int index, occupied;
        do{
            index  = randomBoardIndex();
            occupied = getBoardValueByIndex(index);
        }while(occupied);
        setBoardValueByIndex(index, playerId);
    }
}

void printBoard(){
    printf("Printing board...\n");
    int i, j;
    for(i = 0; i < nrows; i++){
        for(j = 0; j < ncols; j++)
            printf("%3d", board[i][j]);
        printf("\n");
    }
    //I avoided this clever trick because it hinders readability:
    //printf("%3d%s", board[i][j], j%ncols == (ncols-1) ? "\n" : "");
    
}

int randomBoardIndex(){
    return randomIntegerLeftInclusive(0,boardSize);
}

int getBoardValueByIndex(int index){
    return *(board[0] + index);
}

void setBoardValueByIndex(int index, int value){
    *(board[0] + index) = value;
}

void ind2sub(int ind, int *i, int *j){
    *i = floor(ind/(ncols));
    *j = ind%ncols;
}

void ind2sub_generic(int ind, int *i, int *j, int jmax){
    *i = floor(ind/(jmax+1));
    *j = ind%(jmax+1);
}

int sub2ind(int i, int j){
    return i*ncols + j;
}

//overloading would be useful here
//this method is used only in testing
int sub2ind_generic(int i, int j, int jmax){
    return i*jmax + j;
}

float roundFloat(float number, int decimalPlaces){
    int power = pow(10,decimalPlaces);
    return round(number*power)/power;
}

float randomFloatInclusive(float lowerLimit, float upperLimit){
    return lowerLimit + (float)rand() / (float)RAND_MAX * (upperLimit - lowerLimit);
}

int randomIntegerInclusive(int lowerLimit, int upperLimit){
    return lowerLimit + rand() % (upperLimit+1 - lowerLimit);
}

int randomIntegerLeftInclusive(int lowerLimit, int upperLimit){
    if(lowerLimit == upperLimit) return 0;
    return lowerLimit + rand() % (upperLimit - lowerLimit);
}

/*********
 * TESTS *
 *********/

//assuming board greater than or equal to 5x5
void test_buildFirstNeighbours(){
    printf("Testing buildFirstNeighbours (compare with indexed board whether the number at the center has the neighbours printed in each 3x3 matrix)...\n");
    buildFirstNeighbours();
    int samplePlayers[9];
    samplePlayers[0] = sub2ind(0,0);
    samplePlayers[1] = sub2ind(nrows-1,ncols-1);
    samplePlayers[2] = sub2ind(0,ncols-1);
    samplePlayers[3] = sub2ind(nrows-1,0);
    samplePlayers[4] = sub2ind(0,4);
    samplePlayers[5] = sub2ind(nrows-1,4);
    samplePlayers[6] = sub2ind(4,0);
    samplePlayers[7] = sub2ind(4,ncols-1);
    samplePlayers[8] = sub2ind(4,4);
    int i, playerPosition;
    for(i = 0; i < 9; i++){
        playerPosition = samplePlayers[i];
        printf("%3d", sub2ind(neighbours[playerPosition][topLeft][0], neighbours[playerPosition][topLeft][1]));
        printf("%3d", sub2ind(neighbours[playerPosition][top][0], neighbours[playerPosition][top][1]));
        printf("%3d\n", sub2ind(neighbours[playerPosition][topRight][0], neighbours[playerPosition][topRight][1]));
        printf("%3d", sub2ind(neighbours[playerPosition][left][0], neighbours[playerPosition][left][1]));
        printf("%3d", playerPosition);
        printf("%3d\n", sub2ind(neighbours[playerPosition][right][0], neighbours[playerPosition][right][1]));
        printf("%3d", sub2ind(neighbours[playerPosition][bottomLeft][0], neighbours[playerPosition][bottomLeft][1]));
        printf("%3d", sub2ind(neighbours[playerPosition][bottom][0], neighbours[playerPosition][bottom][1]));
        printf("%3d\n\n", sub2ind(neighbours[playerPosition][bottomRight][0], neighbours[playerPosition][bottomRight][1]));
    }
}

//assumes players have been set up
//i.e. test assignment by value vs by reference
void test_WhetherChangingAStructChangesTheCollectionItCameFrom(){
    printf("Testing WhetherChangingAStructChangesTheCollectionItCameFrom...\n\n");
    printf("Before: id: %d, earned: %.2f\n", players[10].id, players[10].earned);
    struct Player player = players[10];
    player.id += 3;
    player.earned += 30;
    printf("After (out of collection): id: %d, earned: %.2f\n", players[10].id, players[10].earned);
    players[10].id += 2;
    players[10].earned += 10;
    printf("After (in collection): id: %d, earned: %.2f\n", players[10].id, players[10].earned);
    
    players[10] = player;
    printf("After (out of collection assigning struct back): id: %d, earned: %.2f\n", players[10].id, players[10].earned);
}

void test_randomFloatInclusive(float lowerLimit, float upperLimit){
    printf("Testing randomFloatInclusive [%.2f,%.2f]...\n", lowerLimit, upperLimit);
    int i;
    for(i = 0; i < 100; i++)
        printf("%.2f ", randomFloatInclusive(lowerLimit,upperLimit));
    printf("\n");
}

//not passing for the lower limit!
void boundsTest_randomFloatInclusive(float lowerLimit, float upperLimit){
    printf("Testing for bounds of randomFloatInclusive [%.2f,%.2f]...\n", lowerLimit, upperLimit);
    int lowerLimitFound = 0;
    int upperLimitFound = 0;
    int maximumAttempts = INT_MAX;
    int attempt;
    for(attempt = 0; attempt < maximumAttempts; attempt++){
        float value = randomFloatInclusive(lowerLimit,upperLimit);
        if(value == lowerLimit) lowerLimitFound = 1;
        if(value == upperLimit) upperLimitFound = 1;
        //if(lowerLimitFound) printf("Lower limit found!\n");
        //if(upperLimitFound) printf("Upper limit found!\n");
        if(lowerLimitFound && upperLimitFound){
            printf("Limits found!\n");
            break;
        }
    }
    printf("Finished testing randomFloatInclusive [0,1]...\n");
}

void test_randomIntegerInclusive(){
    printf("Testing randomIntegerInclusive [0,10]...\n");
    int i;
    for(i = 0; i < 100; i++)
        printf("%d ", randomIntegerInclusive(0,10));
    printf("\n");
}

void test_randomIntegerLeftInclusive(){
    printf("Testing randomIntegerLeftInclusive [0,10)...\n");
    int i;
    for(i = 0; i < 100; i++)
        printf("%d ", randomIntegerLeftInclusive(0,10));
    printf("\n");
}

void test_ind2sub(){
    printf("Testing ind2sub (converting index values to row,column values)...\n");
    int n, i, j;
    for(n = 0; n < 100; n++){
        ind2sub(n, &i, &j);
        printf("%d -> (%d,%d)\n", n, i, j);
    }
}

void test_sub2ind_generic(){
    printf("Testing sub2ind_generic (converting row,column values to index values)...\n");
    int i, j, imax = 3, jmax = 3;
    for(i = 0; i < imax; i++)
        for(j = 0; j < jmax; j++)
            printf("(%d,%d) -> %d\n", i, j, sub2ind_generic(i, j, jmax));
}

void test_readingMatrixWithPointer(){
    printf("Testing readingMatrixWithPointer...\n");
    int i, j;
    const int imax = 3, jmax = 3;
    int matrix[imax][jmax];
    for(i = 0; i < imax; i++)
        for(j = 0; j < jmax; j++)
            matrix[i][j] = sub2ind_generic(i, j, jmax);
    for(i = 0; i < imax; i++)
        for(j = 0; j < jmax; j++)
            printf("%d ", *(matrix[0] + sub2ind_generic(i, j, jmax)));
    printf("\n");
}

void test_writingMatrixWithPointer(){
    printf("Testing writingMatrixWithPointer...\n");
    int i, j;
    const int imax = 3, jmax = 3;
    int matrix[imax][jmax];
    for(i = 0; i < imax; i++)
        for(j = 0; j < jmax; j++)
            *(matrix[0] + sub2ind_generic(i, j, jmax)) = sub2ind_generic(i, j, jmax);
    for(i = 0; i < imax; i++)
        for(j = 0; j < jmax; j++)
            printf("%d ", matrix[i][j]);
    printf("\n");
}

void test_flooringFloat(){
    int i= 0;
    for(i = 0; i < 100; i++){
        float number = randomFloatInclusive(0, 100000);
        int power = pow(10,2);
        printf("Normal: %f, Floored: %f\n", number, round(number*power)/power);
    }
}

/***************
 * DIAGNOSTICS *
 ***************/

void printIndexedBoard(){
    printf("Printing indexed board...\n");
    int i, j;
    for(i = 0; i < nrows; i++){
        for(j = 0; j < ncols; j++)
            printf("%3d", sub2ind(i, j));
        printf("\n");
    }
}

////I'm sure there's a neater way to do this,
////I won't explore it further because it's called only once
////suggestion: enlarged the matrix by repeating rows/columns
//void buildFirstNeighboursDeprecated(){
//    int i, j, index;
//    int firstColumn = 0, lastColumn = ncols-1;
//    int firstRow = 0, lastRow = nrows-1;
//    for(i = 0; i < nrows; i++)
//        for(j = 0; j < ncols; j++){
//            index = sub2ind(i, j);
//            neighbours[index][topLeft]      = sub2ind(i-1,j-1);
//            neighbours[index][top]          = sub2ind(i-1,j);
//            neighbours[index][topRight]     = sub2ind(i-1,j+1);
//            neighbours[index][left]         = sub2ind(i  ,j-1);
//            neighbours[index][right]        = sub2ind(i  ,j+1);
//            neighbours[index][bottomLeft]   = sub2ind(i+1,j-1);
//            neighbours[index][bottom]       = sub2ind(i+1,j);
//            neighbours[index][bottomRight]  = sub2ind(i+1,j+1);
//            //Fixing negative indices:
//            if(j == firstColumn){
//                neighbours[index][topLeft] += ncols;
//                neighbours[index][left] += ncols;
//                neighbours[index][bottomLeft] += ncols;
//            }
//            else if(j == lastColumn){
//                neighbours[index][topRight] -= ncols;
//                neighbours[index][right] -= ncols;
//                neighbours[index][bottomRight] -= ncols;
//            }
//            if(i == firstRow){
//                neighbours[index][topLeft] += boardSize;
//                neighbours[index][top] += boardSize;
//                neighbours[index][topRight] += boardSize;
//            }
//            else if(i == lastRow){
//                neighbours[index][bottomLeft] -= boardSize;
//                neighbours[index][bottom] -= boardSize;
//                neighbours[index][bottomRight] -= boardSize;
//            }
//        }
//}

////assuming board greater than or equal to 5x5
//void test_buildFirstNeighboursDeprecated(){
//    printf("Testing buildFirstNeighboursDeprecated (compare with indexed board whether the number at the center has the neighbours printed in each 3x3 matrix)...\n");
//    buildFirstNeighbours();
//    int samplePlayers[9];
//    samplePlayers[0] = sub2ind(0,0);
//    samplePlayers[1] = sub2ind(nrows-1,ncols-1);
//    samplePlayers[2] = sub2ind(0,ncols-1);
//    samplePlayers[3] = sub2ind(nrows-1,0);
//    samplePlayers[4] = sub2ind(0,4);
//    samplePlayers[5] = sub2ind(nrows-1,4);
//    samplePlayers[6] = sub2ind(4,0);
//    samplePlayers[7] = sub2ind(4,ncols-1);
//    samplePlayers[8] = sub2ind(4,4);
//    int i, index;
//    for(i = 0; i < 9; i++){
//        index = samplePlayers[i];
//        printf("%3d", neighbours[index][topLeft]);
//        printf("%3d", neighbours[index][top]);
//        printf("%3d\n", neighbours[index][topRight]);
//        printf("%3d", neighbours[index][left]);
//        printf("%3d", index);
//        printf("%3d\n", neighbours[index][right]);
//        printf("%3d", neighbours[index][bottomLeft]);
//        printf("%3d", neighbours[index][bottom]);
//        printf("%3d\n\n", neighbours[index][bottomRight]);
//    }
//}
//



