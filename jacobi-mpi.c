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
            srand(64591);
            int* vetorB = (int*)malloc(orderOfMatrix*sizeof(int));
            int* vetEnvioMatrix = (int*)malloc(orderOfMatrix*orderOfMatrix*sizeof(int));
            int* somaProvisoria = (int*)malloc(orderOfMatrix*sizeof(int));
            int** matrix = (int**)malloc(orderOfMatrix*sizeof(int*));
            for(i=0;i<orderOfMatrix;i++){
                matrix[i]=(int*)malloc(orderOfMatrix*sizeof(int));
                somaProvisoria[i]=0;
            }
            int* linhaTeste = (int*)malloc(orderOfMatrix*sizeof(int));
            //Preenchimento provisorio
            int count=0;
            for(i=0;i<orderOfMatrix;i++)
                for(j=0;j<orderOfMatrix;j++){
                    if (i != j){
                        matrix[i][j] = (rand() % (RAND_UPPER_BOUND2 - RAND_LOWER_BOUND + 1)) + RAND_LOWER_BOUND;
                       somaProvisoria[i] =somaProvisoria[i] + matrix[i][j];
                    }
                    matrix[i][i] = (rand() % ((RAND_UPPER_BOUND2 + somaProvisoria[i]) - RAND_LOWER_BOUND + 1)) + (somaProvisoria[i] + 1); // posicao diagonal principal > soma dos outros elementos.
                    vetorB[i] = (rand() % (10*RAND_UPPER_BOUND - RAND_LOWER_BOUND + 1)) + 5*RAND_LOWER_BOUND;
                    if(i==linhaEscolhida)
                        linhaTeste[j]=matrix[i][j];
                }
            j=0;
            int k=0;
            for(i=0;i<orderOfMatrix*orderOfMatrix;i++)
            {
                if(j==orderOfMatrix){
                    j=0;
                    k++;
                }
                vetEnvioMatrix[i] = matrix[k][j];
                printf("%d ",vetEnvioMatrix[i]);
                j++;
            }
            //*****************************GERANDO MATRIZ**************************************
            //printf vetor b
            printf("\n vetorb \n");
            for(i=0;i<orderOfMatrix;i++){
                printf("%d ", vetorB[i]);
            }
            //Tentativa inicial envio variavel
            int* vetorRecebLinhas = (int*)malloc(proporcao*sizeof(int));
            int* vetorRecebimentoB = (int*)malloc(proporcao*sizeof(int));
            int* vetorQuantidades = (int*)malloc(NumberOfProcecess*sizeof(int));
            int* vetorDisplacements = (int*)malloc(NumberOfProcecess*sizeof(int));
            int* vetorQuantidadesB=(int*)malloc(NumberOfProcecess*sizeof(int));
            int* vetorDisplacementsB=(int*)malloc(NumberOfProcecess*sizeof(int));
            //Preenchimento quantidades e displacements
            for(i=0;i<NumberOfProcecess;i++){
                if(i==NumberOfProcecess-1){
                    vetorQuantidades[i] = quantidadePiorCaso;
                    vetorDisplacements[i] = proporcao*orderOfMatrix*i;
                    vetorQuantidadesB[i] = quantidadePiorCasoB;
                    vetorDisplacementsB[i] = proporcao*i;
                }else{
                    vetorQuantidades[i] = proporcao*orderOfMatrix;
                    vetorDisplacements[i] = proporcao*orderOfMatrix*i;
                    vetorQuantidadesB[i] = proporcao;
                    vetorDisplacementsB[i] = proporcao*i;
                }
            }
            MPI_Scatterv(vetEnvioMatrix, vetorQuantidades, vetorDisplacements, MPI_INT, vetorRecebLinhas, orderOfMatrix*proporcao, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Scatterv(vetorB, vetorQuantidadesB, vetorDisplacementsB, MPI_INT, vetorRecebimentoB, proporcao, MPI_INT, 0, MPI_COMM_WORLD);
            //DEPOIS DO SCATTER PODE DAR FREE NA MATRIX, vetQuantidades e vetDisplacements.
            //PRINT VETORB RECEBIDO
            printf("\n vetor B recebido \n");
            for(i=0;i<proporcao;i++){
                printf("%d ",vetorRecebimentoB[i]);
            }
            //PRINT VETORB RECEBIDO
            for(i=0;i<orderOfMatrix;i++)
                free(matrix[i]);
            free(matrix);
            free(vetorDisplacements);
            free(vetorQuantidades);
            free(vetorB);
            free(somaProvisoria);
        }else
        {   //Todos filhos com exceção do 0 executam aqui.
            int **matrix;
            int *vetorQuantidades,*vetorDisplacements,*vetEnvioMatrix,*vetorB,*vetorDisplacementsB,*vetorQuantidadesB;
            int quantoReceber = (my_rank == (NumberOfProcecess-1)) ? quantidadePiorCaso : proporcao*orderOfMatrix;
            int quantoReceberB = (my_rank == (NumberOfProcecess-1)) ? quantidadePiorCasoB : proporcao;
            int* vetorRecebLinhas = (int*)malloc(quantoReceber*sizeof(int));
            int* vetorRecebimentoB = (int*)malloc(quantoReceberB*sizeof(int));
            MPI_Scatterv(vetEnvioMatrix, vetorQuantidades, vetorDisplacements, MPI_INT, vetorRecebLinhas, quantoReceber, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Scatterv(vetorB, vetorQuantidadesB, vetorDisplacementsB, MPI_INT, vetorRecebimentoB, quantoReceberB, MPI_INT, 0, MPI_COMM_WORLD);
            printf("\n vetor B recebido \n");
            for(i=0;i<quantoReceberB;i++){
                printf("%d ",vetorRecebimentoB[i]);
            }
        }

        printf("\n");
    }// Fim do else que separa os spawns do root process
    fflush(stdout);
    MPI_Finalize();
    return 0;
}