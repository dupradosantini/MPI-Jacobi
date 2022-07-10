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
        int i,j;

        //Calculo proporcao
        int proporcao = orderOfMatrix/NumberOfProcecess;
        if(proporcao == 0){ // Se o numero do processos for maior que a ordem da matriz passaremos 1 linha pra cada e haverao processos ociosos.
            proporcao=1;
            NumberOfProcecess=orderOfMatrix; // Caso haja mais processos que linhas na matriz, define o maximo de processos assim.
        }
        int quantidadePiorCaso = proporcao*orderOfMatrix + (orderOfMatrix%NumberOfProcecess)*orderOfMatrix; //Utilizado para determinar quantos elementos da matriz serão enviados no ultimo processo.
        int quantidadePiorCasoB = proporcao + (orderOfMatrix%NumberOfProcecess); // Utilizado para determinar quantos elementos do vetor B serão enviados ao ultimo processo.
        //Fim calculo proporcao

        if(my_rank==0)
        {
            MPI_Status status;
            MPI_Recv(&linhaEscolhida, 1, MPI_INT, 0, 1, parentcomm, &status); //Recebimento da linha de testagem.
            printf("I'm the spawned process number %d on processor %s.\n A linha escolhida foi: %d \n", my_rank, processor_name,linhaEscolhida);
            //*****************************GERANDO MATRIZ**************************************
            int** matrix = (int**)malloc(orderOfMatrix*sizeof(int*));
            for(i=0;i<orderOfMatrix;i++){
                matrix[i]=(int*)malloc(orderOfMatrix*sizeof(int));
            }
            int* linhaTeste = (int*)malloc(orderOfMatrix*sizeof(int));
            //Preenchimento provisorio
            int count=0;
            for(i=0;i<orderOfMatrix;i++)
                for(j=0;j<orderOfMatrix;j++){
                    matrix[i][j] = ++count;
                    if(i==linhaEscolhida)
                        linhaTeste[j]=matrix[i][j];
                }
            //Print matriz
            printf("\n Print matriz: \n");
            for(i=0;i<orderOfMatrix;i++){
                for(j=0;j<orderOfMatrix;j++)
                    printf("%d ",matrix[i][j]);
                printf("\n");
            }
            //Print linha escolhida
            printf("\n Print linha teste: \n");
            for(i=0;i<orderOfMatrix;i++)
                printf("%d ", linhaTeste[i]);
            //*****************************GERANDO MATRIZ**************************************

            //Tentativa inicial envio variavel
            int* vetor_rec = (int*)malloc(proporcao*sizeof(int));
            int* vetorQuantidades = (int*)malloc(NumberOfProcecess*sizeof(int));
            int* vetorDisplacements = (int*)malloc(NumberOfProcecess*sizeof(int));
            //Preenchimento quantidades e displacements
            for(i=0;i<NumberOfProcecess;i++){
                if(i=NumberOfProcecess-1){
                    vetorQuantidades[i] = quantidadePiorCaso;
                    vetorDisplacements[i] = proporcao*orderOfMatrix*i;
                }else{
                    vetorQuantidades[i] = proporcao*orderOfMatrix;
                    vetorDisplacements[i] = proporcao*orderOfMatrix*i;
                }
            }

            MPI_Scatterv(matrix, vetorQuantidades, vetorDisplacements, MPI_INT, vetor_rec, orderOfMatrix*proporcao, MPI_INT, 0, MPI_COMM_WORLD);
            //DEPOIS DO SCATTER PODE DAR FREE NA MATRIX, vetQuantidades e vetDisplacements.
        }else
        {   //Todos filhos com exceção do 0 executam aqui.
            printf("\nI'm the spawned process number %d on processor %s.\n", my_rank, processor_name);
            //Tentativa inicial envio variavel
            int **matrix;
            int *vetorQuantidades,*vetorDisplacements;
            int quantoReceber = (my_rank == (NumberOfProcecess-1)) ? quantidadePiorCaso : proporcao*orderOfMatrix;
            int* vetor_rec = (int*)malloc(quantoReceber*sizeof(int));
            MPI_Scatterv(matrix, vetorQuantidades, vetorDisplacements, MPI_INT, vetor_rec, quantoReceber, MPI_INT, 0, MPI_COMM_WORLD);
        }
    }// Fim do else que separa os spawns do root process
    fflush(stdout);
    MPI_Finalize();
    return 0;
}