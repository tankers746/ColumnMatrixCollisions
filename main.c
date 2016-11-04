#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>
//#include <libiomp/omp.h>
/*
 *
 */
#define COLS 500
#define ROWS 4400
#define BLOCK_SIZE 4
#define KEY_SIZE 14
#define MASTER 0/* root process's rank          */


struct element {
    int index;
    //long long key;
    float value;
};

struct block {
    int index;
    long long signature;
};

//a group of elements from the same column, either a block or neighbourhood
struct colElementGroup {
    int count;
    int col;
    struct element *elements;
    //optional (used in block)
    long long signature;
    int blockCount;
};

//a collection of element-groups
struct elementGroups {
    int count;
    struct colElementGroup *groups;
    //optional (used in neighbourhood)
    int blockCount;
};

struct collisions {
    int count;
    struct elementGroups *collisions;
};

long long *loadKeys() {
    long long *keys = malloc((ROWS+1) * sizeof(long long));
    char buffer[KEY_SIZE+2];
    char *record, *line;
    int i = 0;
    FILE *fstream = fopen("keys.txt", "r");
    if (fstream == NULL)
    {
        printf("\n Opening key vector failed ");
        return keys;
    }
    for(int i = 0; i < ROWS; i++) {
        line = fgets(buffer, sizeof(buffer), fstream);
        if(line != NULL) {
            keys[i] = atoll(line);
        }
    }
    return keys;
}

float **loadMatrix() {
    float **data = malloc(COLS * sizeof(float*));
    for(int i = 0; i < COLS; i++) {
        data[i] = malloc(ROWS * sizeof(float));
    }
    
    int bufsize = COLS * sizeof(char) * 10;
    char buffer[bufsize];
    char *record, *line;
    FILE *fstream = fopen("data.txt", "r");
    if (fstream == NULL)
    {
        printf("\n Opening matrix failed ");
        return data;
    }
    for(int i = 0; i < ROWS; i++) {
        line = fgets(buffer, sizeof(buffer), fstream);
        if(line != NULL) {
            int j = 0;
            record = strtok(line, ",");
            while (record != NULL) {
                data[j++][i] = atof(record);
                record = strtok(NULL, ",");
            }
        }
    }
    return data;
}

struct element **getElementMatrix(float **data, long long *keys) {
    struct element **elementMatrix = malloc(COLS * sizeof(struct element*));
    for(int i = 0; i < COLS; i++) {
        elementMatrix[i] = malloc(ROWS * sizeof(struct element)); 
        for(int j = 0; j < ROWS; j++) {
            elementMatrix[i][j].index = j;
            elementMatrix[i][j].value = data[i][j];
        }
    }  
    return elementMatrix;
}

void printBlocks(struct elementGroups blocks, float **data, long long *keys) {
    for(int i = 0; i < blocks.count; i++) {
        printf("col %i [", blocks.groups[i].col);
        for(int j = 0; j < BLOCK_SIZE; j++) {
            struct colElementGroup b = blocks.groups[i];
            printf("[%i] %f, key %lld ", b.elements[j].index, data[b.col][b.elements[j].index], keys[b.elements[j].index]);
        }
        printf("] - sig %lld\n", blocks.groups[i].signature);
    }
    
}

void printCollisions(struct collisions c, float **data) {
    for(int i = 0; i < c.count; i++) {
        struct elementGroups blocks = c.collisions[i];
        printf("\nsignature %llu\n", blocks.groups[0].signature);
        for(int j = 0; j < blocks.count; j++) {
            struct colElementGroup b = blocks.groups[j];
            printf("col %i: [", b.col);
            for(int k = 0; k  < BLOCK_SIZE; k++) {
                printf("%f, ", data[b.col][b.elements[k].index]);
            }
            printf("]\n");
        }
    }
}

long long getSignature(long long *keys, struct element *elements) {
    long long signature = 0;
    for (int i = 0; i < BLOCK_SIZE; i++) {
        signature += keys[elements[i].index];
    }
    return signature;
}

int elementComp(const void* p1, const void* p2) {
    const struct element *elem1 = p1;
    const struct element *elem2 = p2;
    
    if(elem1->value < elem2->value) {
        return -1;
    }
    return (elem1->value > elem2->value);
}

int blockComp(const void* p1, const void* p2) {
    const struct block *elem1 = p1;
    const struct block *elem2 = p2;
    
    if(elem1->signature < elem2->signature) {
        return -1;
    }
    return (elem1->signature > elem2->signature);
}

int groupComp_sig(const void* p1, const void* p2) {
    const struct colElementGroup *elem1 = p1;
    const struct colElementGroup *elem2 = p2;
    
    if(elem1->signature < elem2->signature) {
        return -1;
    }
    return (elem1->signature > elem2->signature);
}

int groupComp_size(const void* p1, const void* p2) {
    const struct colElementGroup *elem1 = p1;
    const struct colElementGroup *elem2 = p2;
    
    if(elem1->count < elem2->count) {
        return -1;
    }
    return (elem1->count > elem2->count);
}


struct elementGroups getNeighbourhoods(int col, float dia, struct element **elementMatrix) {
    //sort the column by size of the value
    struct element *column = elementMatrix[col];
    qsort(column, ROWS, sizeof(struct element), elementComp);
    
    struct elementGroups neighbourhoods;
    neighbourhoods.groups = malloc(COLS * sizeof(struct colElementGroup));
    neighbourhoods.count = 0;
    neighbourhoods.blockCount = 0;
    
    struct element temp[ROWS];
    memset(&temp, -1, sizeof(temp));
    
    float min = 0, max = 0;
    int neighbourhoodSize = 0;
    int lastNeighbourhoodSize = 0;
    
    for(int i = 0; i < ROWS; i++) {
        //fprintf(stderr,"[%d] %f\n",col[i].index, col[i].value);
        if(temp[0].index == -1) {
            min = max = column[i].value;
            temp[neighbourhoodSize++] = column[i];
        } else {
            if (fabs(column[i].value - min) < dia && fabs(column[i].value - max) < dia) {
                temp[neighbourhoodSize++] = column[i];
                if(column[i].value < min) {
                    min = column[i].value;
                } else if(column[i].value > max) {
                    max = column[i].value;
                }
                
            } else {
                 /*
                 -If the block is not empty (-1)
                 -If the block is larger than 3
                 -If the block is larger or equal in size to the previous block
                 This ensures that the current block is not a sub-block of the previous block
                 */
                if((temp[0].index != -1) && (neighbourhoodSize >= BLOCK_SIZE) && (neighbourhoodSize >= lastNeighbourhoodSize)) {
                    //formula for working out combinations of size k(blocksize) for n(tempcount) values
                    neighbourhoods.groups[neighbourhoods.count].blockCount = round(exp(lgamma(neighbourhoodSize+1)-lgamma(neighbourhoodSize-BLOCK_SIZE+1))/tgamma(BLOCK_SIZE+1));
                    
                    //allocate the memory to store the neighbourhood's elements
                    neighbourhoods.groups[neighbourhoods.count].elements = malloc(neighbourhoodSize * sizeof(struct element));
                    neighbourhoods.groups[neighbourhoods.count].count = neighbourhoodSize;                    
                    neighbourhoods.groups[neighbourhoods.count].col = col;
                    //loops over all the elements in temp and adds them to the neighbourhoods struct
                    for(int j = 0; j < neighbourhoodSize; j++) {
                        neighbourhoods.groups[neighbourhoods.count].elements[j] = temp[j];
                    }
                    neighbourhoods.blockCount += neighbourhoods.groups[neighbourhoods.count++].blockCount;
                }
                min = max = column[i].value;
                memset(&temp, -1, sizeof(temp));
                i = i - neighbourhoodSize;
                lastNeighbourhoodSize = neighbourhoodSize;
                neighbourhoodSize = 0;
            }
        }
    }
    return neighbourhoods;
}

//takes an array of many groups and combines them all into one large struct
struct elementGroups groupArrayToStruct(struct elementGroups temp[], int arrayCount, int totalCount) {
    struct elementGroups combined;
    combined.groups = malloc(totalCount * sizeof(struct colElementGroup));
    combined.count = 0;
    combined.blockCount = 0;
    for(int i = 0; i < arrayCount; i++) {
        struct elementGroups column = temp[i];
        combined.blockCount += column.blockCount;
        for(int j = 0; j < column.count; j++) {
            combined.groups[combined.count++] = column.groups[j];
        }
        free(column.groups);
    }
    return combined;
}

struct elementGroups getAllNeighbourhoods(float dia, struct element **elementMatrix) {
    struct elementGroups temp[COLS];
    int totalNeighbourhoodCount = 0;
    #pragma omp parallel for
    for(int i = 0; i < COLS; i++) {
        temp[i] = getNeighbourhoods(i, dia, elementMatrix);
        totalNeighbourhoodCount += temp[i].count;
        free(elementMatrix[i]);
    }
    return groupArrayToStruct(temp, COLS, totalNeighbourhoodCount);
}

void findCombinations(struct elementGroups *blocks, struct colElementGroup neighbourhood, int start, int currLen, bool used[], long long *keys) {
    if (currLen == BLOCK_SIZE) {
        int blockCount;
        #pragma omp atomic capture
        {
            blockCount = blocks->count;
            blocks->count++;
        }   
        int elementCount = 0;
        blocks->groups[blockCount].elements = malloc(BLOCK_SIZE * sizeof(struct element));
        for (int i = 0; i < neighbourhood.count; i++) {
            if (used[i] == true) {
                blocks->groups[blockCount].elements[elementCount] = neighbourhood.elements[i];
                elementCount++;
            }
        }
        blocks->groups[blockCount].count = elementCount;
        blocks->groups[blockCount].signature = getSignature(keys, blocks->groups[blockCount].elements);
        blocks->groups[blockCount].col = neighbourhood.col;        
        return;
    }
    if (start == neighbourhood.count) {
        return;
    }
    used[start] = true;
    findCombinations(blocks, neighbourhood, start + 1, currLen + 1, used, keys); 
    used[start] = false;
    findCombinations(blocks, neighbourhood, start + 1, currLen, used, keys);
}

struct elementGroups getBlocks(struct elementGroups neighbourhoods, long long *keys) {
    struct elementGroups blocks;
    blocks.groups = malloc(neighbourhoods.blockCount * sizeof(struct colElementGroup));
    blocks.count = 0;
    #pragma omp parallel for
    for(int i = 0; i < neighbourhoods.count; i++) {
        int length = neighbourhoods.groups[i].count;
        bool used[length];
        memset(used, false, sizeof(used));
        findCombinations(&blocks, neighbourhoods.groups[i], 0, 0, used, keys);
        free(neighbourhoods.groups[i].elements);
    }
    free(neighbourhoods.groups);
    return blocks;
}

int triangularNumber(int n) {  
    if (n == 1) return 1;  
    return n + triangularNumber(n-1);  
} 

struct elementGroups getBlocksParallel(struct elementGroups neighbourhoods, long long *keys) {
    qsort(neighbourhoods.groups, neighbourhoods.count, sizeof(struct colElementGroup), groupComp_size);
    int n_threads = omp_get_max_threads();
    int bins[n_threads];
    int bigNumbers[n_threads];
    struct elementGroups nodeBlocks[n_threads];
    
    int allocated = 0;
    for (int i = 0; i < n_threads; i++) {
        int remainder = neighbourhoods.count - allocated;
        int buckets = (n_threads - i);
        bins[i] = remainder / buckets;
        allocated += bins[i];
        //formula for best dividing the amount of big numbers based on the number of threads
        bigNumbers[i] = 1 + (i+1)*50*pow(n_threads, -1.85); 
    }

  
    int totalBlockCount = 0;
    #pragma omp parallel for
    for (int i = 0; i < n_threads; i++) {        
        struct elementGroups node_neighbourhoods;
        node_neighbourhoods.groups = malloc(bins[i] * sizeof(struct colElementGroup));
        node_neighbourhoods.blockCount = 0;
        node_neighbourhoods.count = 0;
        //how many big numbers to take off the bottom  
        //39 nodes int bigNumbers = 1 + (i/(N_NODES*4/9));       
        int topIndex = 0;
        int bottomIndex = 0;
        for(int j = 0; j < i; j++) {
            bottomIndex += bigNumbers[j];
            topIndex += bins[j] - bigNumbers[j];
        }
        //each node takes some from the bottom (big) and the rest from the top (small)  
        for(int j = 0; j < bins[i]; j++) {
            //takes from the bottom (big)
            if(j < bigNumbers[i]) {
                node_neighbourhoods.groups[j] = neighbourhoods.groups[neighbourhoods.count - 1 - (bottomIndex + j)]; 
            //takes from the top (small)
            } else {              
                node_neighbourhoods.groups[j] = neighbourhoods.groups[topIndex + j - bigNumbers[i]];
            }
            node_neighbourhoods.blockCount += node_neighbourhoods.groups[j].blockCount;
            node_neighbourhoods.count++;
        }
        totalBlockCount += node_neighbourhoods.blockCount;
        nodeBlocks[i] = getBlocks(node_neighbourhoods, keys);
        //printf("Thread %d will produce %d blocks\n", i, node_neighbourhoods.blockCount);
        //qsort(nodeBlocks[i].groups, nodeBlocks[i].count, sizeof(struct colElementGroup), groupComp_sig);
    }
    return groupArrayToStruct(nodeBlocks, n_threads, totalBlockCount);
}

struct collisions getCollisions(struct elementGroups blocks, struct block *sortedBlocks) {
    //takes aaaaaaaaaages
    //clock_t start = clock();
    //qsort(blocks.groups, blocks.count, sizeof(struct colElementGroup), groupComp_sig);
    //int msec = (clock() - start) * 1000 / CLOCKS_PER_SEC;
    //printf("Time taken to qqsort %d blocks: %d seconds %d milliseconds\n", blocks.count, msec/1000, msec%1000);       

    struct collisions c;
    c.collisions = malloc(blocks.count * sizeof(struct elementGroups));  
    struct block currentBlock;
    int collisionCount = 0;
    
    int i = 0;
    
    while(i < blocks.count) {
        int blockCount = 0;
        int trueBlockCount = 0;
        //this is used to ensure collisions do not occur inside the same column
        bool columns[COLS];
        memset(columns, false, sizeof(columns));
        do {
            //currentBlock = &blocks.groups[i++];
            currentBlock = sortedBlocks[i++];

            //if no collision in that column
            if(columns[blocks.groups[currentBlock.index].col] == false) {
                trueBlockCount++;
                columns[blocks.groups[currentBlock.index].col] = true;
            }
            blockCount++;
        } while(blocks.groups[currentBlock.index].signature == blocks.groups[i].signature);
        //collision found
        if(trueBlockCount > 1) {
            c.collisions[collisionCount].groups = malloc(trueBlockCount * sizeof(struct colElementGroup));
            c.collisions[collisionCount].count = trueBlockCount;
            for(blockCount; blockCount > 0; blockCount--) {
                //if first instance of collision in that column
                if(columns[blocks.groups[sortedBlocks[i-blockCount].index].col] == true) {
                    c.collisions[collisionCount].groups[--trueBlockCount] = blocks.groups[sortedBlocks[i-blockCount].index];
                    columns[blocks.groups[sortedBlocks[i-blockCount].index].col] = false;
                }                
            }
            collisionCount++;
        }
    }
    c.count = collisionCount;
    return c;
}

void swap(struct block s[], int m, int n) {
    struct block tmp;
    tmp = s[m];
    s[m] = s[n];
    s[n] = tmp;
}

int partition(struct block x[], int first, int last) {
    int pivot;
    int j, i;
    pivot = first;
    i = first;
    j = last;
    while (i < j) {
        //move to the right
        while (x[i].signature <= x[pivot].signature && i < last) {
            i++;
        }
        //move to the left
        while (x[j].signature > x[pivot].signature) {
            j--;
        }
        if (i < j) {
            swap(x, i, j);
        }
    }
    swap(x, pivot, j);
    return j;
}

void quickSortRecursive(struct block x[], int first, int last) {
    int pivot;
    if (first < last) {
        pivot = partition(x, first, last);
        quickSortRecursive(x, first, pivot - 1);
        quickSortRecursive(x, pivot + 1, last);

    }
}

void sort(struct block *globalData, int size) {
    double t_start, t_end;          

    int rank; //process rank   
    int npes; //number of processors   
    MPI_Comm_size(MPI_COMM_WORLD, & npes);
    MPI_Comm_rank(MPI_COMM_WORLD, & rank);

    //start timer and check whether sort can be distributed
    if (rank == MASTER) {
        t_start = MPI_Wtime();
        if (size < npes) {
            printf("Size is less than the number of process, reverting to qsort\n");
            qsort(globalData, size, sizeof(struct block), blockComp);
            return;
        }
        //send the block count from the master node to every other node
        for (int i = 0; i < npes; i++) {
		MPI_Send(&size, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	}       
    }
    //receiving the block count
    MPI_Recv(&size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //getting the size that is used by each process    
    int localsize = size/npes;

    //allocate memory for each node
    struct block *localdata = malloc(localsize * sizeof(struct block));
    if (localdata == NULL) {
        printf("localdata Memory Allocation Failed\n");
        return;
    } 

    //create the block struct as an MPI type
    MPI_Datatype MPI_BLOCK_TYPE;
    MPI_Datatype type[2] = {MPI_INT, MPI_LONG_LONG};    
    int blocklen[2] = {1,1};
    MPI_Aint offsets[2];    
    offsets[0] = offsetof(struct block, index);
    offsets[1] = offsetof(struct block, signature);
    MPI_Type_create_struct(2, blocklen, offsets, type, &MPI_BLOCK_TYPE);
    MPI_Type_commit(&MPI_BLOCK_TYPE);

    //Scatter the blocks to each process
    //Each process receives the array of sructs
    MPI_Scatter(globalData, localsize, MPI_BLOCK_TYPE, localdata, localsize, MPI_BLOCK_TYPE, MASTER, MPI_COMM_WORLD);
    quickSortRecursive(localdata, 0, localsize - 1);

    //merge the locally sorted data
    MPI_Gather(localdata, localsize, MPI_BLOCK_TYPE, globalData, localsize, MPI_BLOCK_TYPE, MASTER, MPI_COMM_WORLD);
    free(localdata);

    if (rank == MASTER) {
        //final sort on master
        quickSortRecursive(globalData, 0, size - 1);
        //stop timer
        t_end = MPI_Wtime();
        printf("Sorting time: %7.4f\n", t_end - t_start);
    }
}

int main(int argc, char* argv[]) {
    clock_t startTotal = clock();    
    long long *keys;
    float **data;    
    struct block * sortedBlocks;
    struct elementGroups b;      

    //initialise MPI
    MPI_Init(NULL, NULL);    
    int rank;
    //get the rank
    MPI_Comm_rank(MPI_COMM_WORLD, & rank);

    if (rank == MASTER) { 
        if(argc > 2) printf("Nodes: %s, PPN: %s\n", argv[1], argv[2]); 

        clock_t start = clock();
        keys = loadKeys();
        data = loadMatrix();
        struct element **elementMatrix = getElementMatrix(data, keys);
        int msec = (clock() - start) * 1000 / CLOCKS_PER_SEC;
        printf("Time taken to load data: %d seconds %d milliseconds\n", msec/1000, msec%1000);

        start = clock();
        struct elementGroups n = getAllNeighbourhoods(0.000001, elementMatrix);
        free(elementMatrix);
        msec = (clock() - start) * 1000 / CLOCKS_PER_SEC;
        printf("Time taken to find %d neighbourhoods: %d seconds %d milliseconds\n", n.count, msec/1000, msec%1000);
    
        start = clock();
        b = getBlocksParallel(n, keys);
        msec = (clock() - start) * 1000 / CLOCKS_PER_SEC;
        printf("Time taken to find %d blocks: %d seconds %d milliseconds\n", b.count, msec/1000, msec%1000);

        sortedBlocks = malloc(b.count * sizeof(struct block));
        for(int i = 0; i < b.count; i++) {
            sortedBlocks[i].index = i;
            sortedBlocks[i].signature = b.groups[i].signature;
        }      
    }
    sort(sortedBlocks, b.count);

    if (rank == MASTER) { 
        struct collisions c = getCollisions(b, sortedBlocks);   
        //printCollisions(c, data);

        int msec = (clock() - startTotal) * 1000 / CLOCKS_PER_SEC;
        printf("Total time taken: %d seconds %d milliseconds\n", msec/1000, msec%1000);
    }
    MPI_Finalize();

    return (EXIT_SUCCESS);
}
