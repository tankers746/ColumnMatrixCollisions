/*
 * This program implements the recursive quicksort algorithm in parallel
 * using the OpenMPI library. The implementation is done in the 
 * following steps:
 *
 * 0. Recursive version of quicksort
 * 1. Read input from file
 * 2. Shows them on the screen
 * 3. Add MPI
 * 4. Start the wall timer
 * 5. Master distribute globaldata to all processes localdata
 * 6. Sort each localdata
 * 7. Gather them to the globaldata
 * 8. Sort semi-sorted globaldata
 * 9. Stop the wall timer
 * 10. Write the duration and final sorted data to the output file
 * 11. Add sorting checker
 * 12. Input and output files entered as command line arguments
 * 13. Get input size from the input filename
 * 14. Change datatype to long long
 *
 * Input:
 * The input file is generated using input_generator.c
 *
 * Compiling:
 *   mpicc   recursive.c -o recursive
 *
 * Running:
 * mpirun -np [number process] <program name> <input file> <output file>
 * e.g: mpirun -np 10 recursive input_100.txt out_recursive.txt
 *
 *
 * File: recursive.c		Author: M. Soulemane
 * Date: 18.01.2016         version: v0.1
 *
 * History: none
 *
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MASTER 0/* root process's rank          */

struct block {
    int index;
    long long signature;
};


void quickSortRecursive(struct block[], long long, long long);
long long partition(struct block[], long long, long long);
void swap(struct block[], long long, long long);
void sortCheckers(long long, struct block[]);

int main(int argc, char * * argv) {
    long long SIZE = 1; /* input size  read from file   */
    /*  rank is the rank of the calling process in the communicator       */
    int rank;
    long long i; /* loop variable                */
    long long retscan; /* return value of scanf        */
    long long tmp; /* temporary variable           */
    double t_start, t_end; /* variable used for clock time */
    long long test_size = 5; /* test loop size's variable    */
    FILE * out, * inp; /* declare output/input stream  */
    int npes; /* number of processes          */
    struct block * globaldata = NULL; /* for global array pointer     */
    struct block * localdata = NULL; /* for local  array pointer     */
    long long localsize; /* for local array size         */

    /*Checking that the run format is properly entered                    */
    if (argc != 3) {
        printf("\n Properly specify the executable, input , output files");
        printf("\nmpirun -np <process nber> %s <input file> <output file>\n", argv[0]);
        exit(1);
    }
    SIZE = 1000000;
    /* Initialize the MPI execution environment                           */
    MPI_Init( & argc, & argv);

    MPI_Comm_size(MPI_COMM_WORLD, & npes);
    MPI_Comm_rank(MPI_COMM_WORLD, & rank);

    if (rank == MASTER) {
        globaldata = malloc(SIZE * sizeof(struct block));
        if (globaldata == NULL) {
            printf("\n\n globaldata Memory Allocation Failed !   \n\n ");
            exit(EXIT_FAILURE);
        }
        inp = fopen(argv[1], "r"); /* Open file for reading        */
        if (inp == NULL) {
            printf("\n\n inp Memory Allocation Failed !  \n\n ");
            exit(EXIT_FAILURE);
        }
        printf("\n\nInput Data \n\n ");
        for (i = 0; i < SIZE; i++) {
            retscan = fscanf(inp, "%lld \t", & tmp);
            globaldata[i].signature = tmp;
        }

        printf("\n\n End Input Data");

        fclose(inp);
        printf("\n\nProcessor %d has data: ", rank);
        for (i = 0; i < test_size; i++) {
            printf("%lld \t", globaldata[i].signature);
        }
        printf("\n");
    }

    /*Start wall  time                                                    */
    if (rank == MASTER) {
        t_start = MPI_Wtime();
    }

    /*Getting the size to be used by each process                         */
    if (SIZE < npes) {
        printf("\n\n SIZE is less than the number of process!  \n\n ");
        exit(EXIT_FAILURE);
    }
    localsize = SIZE / npes;

    /* Allocate memory to localdata of size localsize                     */
    localdata = (struct block * ) malloc(localsize * sizeof(struct block));
    if (localdata == NULL) {
        printf("\n\n localdata Memory Allocation Failed !  \n\n ");
        exit(EXIT_FAILURE);
    }

    /* create a type for struct block */
    MPI_Datatype MPI_BLOCK_TYPE;
    MPI_Datatype type[2] = {MPI_INT, MPI_LONG_LONG};    
    int blocklen[2] = {1,1};
    int offsets[2];    
    offsets[0] = offsetof(struct block, index);
    offsets[1] = offsetof(struct block, signature);
    MPI_Type_create_struct(2, blocklen, offsets, type, &MPI_BLOCK_TYPE);
    MPI_Type_commit(&MPI_BLOCK_TYPE);

    /*Scatter the integers to each number of processes (npes)             */
    MPI_Scatter(globaldata, localsize, MPI_BLOCK_TYPE, localdata,
        localsize, MPI_BLOCK_TYPE, MASTER, MPI_COMM_WORLD);

    /* Perform local sort on each sub data by each process                */
    quickSortRecursive(localdata, 0, localsize - 1);

    /*
     * printf ("\n\nProcessor %d has sorted data \n\n", rank);
     *   for ( i = 0; i<test_size; i++) {
     *     printf ("%lld \t", localdata[i]);
     *  }
     */

    /* Merge locally sorted data of each process by MASTER to globaldata  */
    MPI_Gather(localdata, localsize, MPI_BLOCK_TYPE, globaldata,
        localsize, MPI_BLOCK_TYPE, MASTER, MPI_COMM_WORLD);
    free(localdata);

    if (rank == MASTER) {
        /* Final sorting                                                    */
        quickSortRecursive(globaldata, 0, SIZE - 1);

        /* End wall  time                                                   */
        t_end = MPI_Wtime();

        /* Opening output file to write sorted data                         */
        out = fopen(argv[2], "w");
        if (out == NULL) {
            printf("\n\n out Memory Allocation Failed !  \n\n ");
            exit(EXIT_FAILURE);
        }
        /* Write information to output file                                 */
        fprintf(out, "Recursively Sorted  Data : ");
        fprintf(out, "\n\nInput size : %lld\t", SIZE);
        fprintf(out, "\n\nNber processes : %d\t", npes);
        fprintf(out, "\n\nWall time      : %7.4f\t", t_end - t_start);
        printf("\n\nWall time      : %7.4f\t", t_end - t_start);
        fprintf(out, "\n\n");
        for (i = 0; i < SIZE; i++) {
            fprintf(out, " %lld \t", globaldata[i].signature);
        }
        fclose(out); /* closing the file             */

        /*
        for ( i = 0; i<test_size; i++) {
          printf ("%lld \t", globaldata[i]);
        }
        */
        printf("\n\n");

        /* checking if the final globaldata content is properly sorted      */
        sortCheckers(SIZE, globaldata);

        printf("\n\n");

    }

    if (rank == MASTER) {
        free(globaldata); /* free the allocated memory    */
    }

    /* MPI_Finalize Terminates MPI execution environment                  */
    MPI_Finalize();

    return EXIT_SUCCESS;
}

/* This function divides elements of an array around a pivot element. All
 * elements less than  or equal to the pivot go on the left side and
 * those greater than the pivot go on the right side.
 *
 * Input:		x	input array
 *              first   leftmost element
 *              last    rightmost element
 * Output		none
 * Return value:	j is returned as the pivot element
 * Sideeffects:		none
 *
 */
long long partition(struct block x[], long long first, long long last) {
    long long pivot; /* pivot variable               */
    long long j, i; /* loop variable                */
    pivot = first;
    i = first;
    j = last;

    while (i < j) {
        /* move to the right                                                */
        while (x[i].signature <= x[pivot].signature && i < last) {
            i++;
        }

        /* move to the left                                                 */
        while (x[j].signature > x[pivot].signature) {
            j--;
        }
        if (i < j) {
            swap(x, i, j); /* swap i and j                 */
        }
    }

    swap(x, pivot, j); /* swap pivot and j             */

    return j;
}

/* Swap elements at index m and n of the array s.
 *
 * Input:		s	array
 *              m   left index
 *              n   right index
 * Output		none
 * Return value:	none
 * Sideeffects:		none
 *
 */
void swap(struct block s[], long long m, long long n) {
    struct block tmp; /* temporary variable           */
    tmp = s[m];
    s[m] = s[n];
    s[n] = tmp;
}

void quickSortRecursive(struct block x[], long long first, long long last) {
    long long pivot; /* pivot variable               */

    if (first < last) {
        /* partition the input array x                                        */
        pivot = partition(x, first, last);

        /* recursively sort left side of the pivot                          */
        quickSortRecursive(x, first, pivot - 1);

        /* recursively sort right side of the pivot                         */
        quickSortRecursive(x, pivot + 1, last);

    }
}

/* Checking a list of sorted numbers. Make sure that each number is
 * less or equal to its immediate right neighbours / greater or equal to
 * its immediate left value.
 *
 * input parameters:	SIZE  total number of sorted items
 *                      input array containing sorted items
 *
 * output parameters:	input[index-1], input[index] shown on failure
 * return value:	none
 * side effects:	none
 *
 */
void sortCheckers(long long SIZE, struct block input[]) {
    long long i; /* Variable declaration         */

    for (i = 1; i < SIZE; i++) {
        if (input[i - 1].signature > input[i].signature) {
            printf("\n\n%lld -- %lld \t", input[i - 1].signature, input[i].signature);
            printf("\n\nCheck failed. Array not sorted");
            break;
        }

    }

    printf("\n\nCheck successfully completed. Array Sorted");
}
