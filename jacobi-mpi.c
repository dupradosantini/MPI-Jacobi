// to compile: mpicc jacobi-mpi.c -o jacobi-mpi -fopenmp
// to run: mpirun -np 1 ./jacobi-mpi N P T

#include <math.h>
#include <omp.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    int P = atoi(argv[2]); //Usado somente pelo processo pai 'criador' para um check inicial.
    int N = atoi(argv[1]); //Usado somente pelo processo pai 'criador' para um check inicial.
    int linhaEscolhida;
    int my_rank;
    int errcodes[P];

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;

    MPI_Comm parentcomm, intercomm;

    MPI_Init( &argc, &argv );
    MPI_Comm_get_parent( &parentcomm );
    MPI_Get_processor_name(processor_name, &name_len);

    if(parentcomm == MPI_COMM_NULL) //Processo PAI Bootstrap do programa.
    {
        printf("\nEscolha a linha que será usada no teste do resultado: ");
        scanf("%d", &linhaEscolhida);

        int root = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        if(P<=N){
            MPI_Comm_spawn( "./jacobi-mpi", argv, P, MPI_INFO_NULL, root, MPI_COMM_WORLD, &intercomm, errcodes );
        }else{
            MPI_Comm_spawn( "./jacobi-mpi", argv, N, MPI_INFO_NULL, root, MPI_COMM_WORLD, &intercomm, errcodes );
        }

        printf("I'm the parent number %d on processor %s.\n", my_rank, processor_name);
        MPI_Send(&linhaEscolhida, 1, MPI_INT, 0, 1, intercomm);
        fflush(0);
    }else
    {   //Aqui ocorre a execucao dos processos filhos!
        MPI_Comm_rank(parentcomm, &my_rank);

        int orderOfMatrix = atoi(argv[2]);      //Parametro de entrada, ordem da matriz gerada.
        int NumberOfProcecess = atoi(argv[3]); //Parametro de entrada, numero de processos gerados.
        int numberOfThreads = atoi(argv[4]);  //Parametro de entrada, numero de threads por processo.

        if(my_rank==0)
        {
            MPI_Status status;
            MPI_Recv(&linhaEscolhida, 1, MPI_INT, 0, 1, parentcomm, &status); //Recebimento da linha de testagem.
            printf("I'm the spawned process number %d on processor %s.\n A linha escolhida foi: %d \n", my_rank, processor_name,linhaEscolhida);
        }else
        {   //Todos filhos com exceção do 0 executam aqui.
            printf("I'm the spawned process number %d on processor %s.\n", my_rank, processor_name);
        }
    }// Fim do else que separa os spawns do root process
    fflush(stdout);
    MPI_Finalize();
    return 0;
}