// to compile: mpicc jacobi-mpi.c -o jacobi-mpi -fopenmp
// to run: mpirun -np 1 ./jacobi-mpi N P T

#include <math.h>
#include <omp.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define RAND_UPPER_BOUND 100
#define RAND_UPPER_BOUND2 30
#define RAND_LOWER_BOUND 10

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
        //printf("\nEscolha a linha que será usada no teste do resultado: ");
        //scanf("%d", &linhaEscolhida);
        linhaEscolhida=1;

        if(P == 1){
            printf("Use pelo menos 2 processos \n");
            MPI_Finalize();
            return 0;
        }

        if(N%P != 0){
            printf("A ordem da matriz tem q ser divisível pelo numero de processos.\n");
            MPI_Finalize();
            return 0;
        }

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
        int i,j,k;

        //Calculo quantidade de linhas por processo
        int linhasPorProcesso = orderOfMatrix/NumberOfProcecess;
        printf("\n Quantas linhas por processo: %d  e resto: %d \n",linhasPorProcesso, orderOfMatrix%NumberOfProcecess);
        //Fim calculo linhasPorProcesso
        int *vetorRecebLinhas, *vetorRecebimentoB;
        int **matrix;
        int *vetEnvioMatrix,*vetorB;

        if(my_rank==0)
        {
            MPI_Status status;
            MPI_Recv(&linhaEscolhida, 1, MPI_INT, 0, 1, parentcomm, &status); //Recebimento da linha de testagem.
            printf("I'm the spawned process number %d on processor %s.\n A linha escolhida foi: %d \n", my_rank, processor_name,linhaEscolhida);
            //*****************************GERANDO MATRIZ**************************************
            srand(64591);
            vetorB = (int*)malloc(orderOfMatrix*sizeof(int));
            vetEnvioMatrix = (int*)malloc(orderOfMatrix*orderOfMatrix*sizeof(int));
            int* somaProvisoria = (int*)malloc(orderOfMatrix*sizeof(int));
            matrix = (int**)malloc(orderOfMatrix*sizeof(int*));
            for(i=0;i<orderOfMatrix;i++){
                matrix[i]=(int*)malloc(orderOfMatrix*sizeof(int));
                somaProvisoria[i]=0;
            }
            int* linhaTeste = (int*)malloc(orderOfMatrix*sizeof(int));
            //Preenchimento
            for(i=0;i<orderOfMatrix;i++){
                for(j=0;j<orderOfMatrix;j++){
                    if (i != j){
                        matrix[i][j] = (rand() % (RAND_UPPER_BOUND2 - RAND_LOWER_BOUND + 1)) + RAND_LOWER_BOUND;
                       somaProvisoria[i] =somaProvisoria[i] + matrix[i][j];
                    }
                    matrix[i][i] = (rand() % ((RAND_UPPER_BOUND2 + somaProvisoria[i]) - RAND_LOWER_BOUND + 1)) + (somaProvisoria[i] + 1); // posicao diagonal principal > soma dos outros elementos.
                    if(i==linhaEscolhida)
                        linhaTeste[j]=matrix[i][j];
                }
                vetorB[i] = (rand() % (10*RAND_UPPER_BOUND - RAND_LOWER_BOUND + 1)) + 2*somaProvisoria[i];
            }

            j=0;
            k=0;
            for(i=0;i<orderOfMatrix*orderOfMatrix;i++)
            {
                if(j==orderOfMatrix){
                    j=0;
                    k++;
                    printf("\n");
                }
                vetEnvioMatrix[i] = matrix[k][j];
                printf("%d ",vetEnvioMatrix[i]);
                j++;
            }
            //*****************************GERANDO MATRIZ**************************************  
        }
        vetorRecebLinhas = (int*)malloc(linhasPorProcesso*orderOfMatrix*sizeof(int));
        vetorRecebimentoB = (int*)malloc(linhasPorProcesso*sizeof(int));
        MPI_Scatter(vetEnvioMatrix, linhasPorProcesso*orderOfMatrix, MPI_INT, vetorRecebLinhas, linhasPorProcesso*orderOfMatrix, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Scatter(vetorB, linhasPorProcesso, MPI_INT, vetorRecebimentoB, linhasPorProcesso, MPI_INT, 0, MPI_COMM_WORLD);
        //Print debug
        printf("\nVetor da matriz\n");
        j=0;
        k=0;
        for(i=0;i<linhasPorProcesso*orderOfMatrix;i++){
                if(j==orderOfMatrix){
                    j=0;
                    k++;
                    printf("\n");
                }
                printf("%d ",vetorRecebLinhas[i]);
                j++;
        }
        //fim print dos recebimentos
        printf("\n vetor B recebido \n");
        for(i=0;i<linhasPorProcesso;i++){
                printf("%d ",vetorRecebimentoB[i]);
        }
        
        //Matriz auxiliar
        int **matrizAuxiliar;
        matrizAuxiliar= (int**)malloc(linhasPorProcesso*sizeof(int*));
        for(i=0;i<orderOfMatrix;i++){
            matrizAuxiliar[i]=(int*)malloc(orderOfMatrix*sizeof(int));
        }
        k=0;
        printf("\nPrint da matriz auxiliar\n");
        for(i=0;i<linhasPorProcesso;i++){
            for(j=0;i<orderOfMatrix;j++){
                matrizAuxiliar[i][j] = vetorRecebLinhas[k];
                k++;
                printf("%d ",matrizAuxiliar[i][j]);
            }
            printf("\n");
        }

    }// Fim do else que separa os spawns do root process
    fflush(stdout);
    MPI_Finalize();
    return 0;
}