#include <mpi.h>
#include <stdio.h>
#include <omp.h>

int main(int argc, char** argv) {

    int required=MPI_THREAD_SERIALIZED;
    int provided;
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
   
    int np;
    int pn; 
    #pragma omp parallel default(shared) private(np, pn)
    {
    // Print off a hello world message    
    np = omp_get_max_threads();
    pn = omp_get_thread_num();   
 
    printf("Currently in thread %d of %d available threads in this proces. Rank is %d, world size is %d.\n", pn, np, world_rank, world_size );
    }

    /*
    printf("Hello world from processor %s, rank %d"
           " out of %d processors\n",
           processor_name, world_rank, world_size);
    */
    
    /*  
    // We are assuming at least 2 processes for this task
    if (world_size < 2) {
        fprintf(stderr, "World size must be greater than 1 for %s\n", argv[0]);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int number;
    if (world_rank == 0) {
    // If we are rank 0, set the number to -1 and send it to process 1
        number = -1;
        MPI_Send(&number, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    } else if (world_rank == 1) {
        MPI_Recv(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Process 1 received number %d from process 0\n", number);
    }
    */
    // Finalize the MPI environment.
    MPI_Finalize();
}


