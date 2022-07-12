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
        printf("\nEscolha a linha que será usada no teste do resultado: ");
        scanf("%d", &linhaEscolhida);
        //linhaEscolhida=1;

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
        //Fim calculo linhasPorProcesso
        int *vetorRecebLinhas, *vetorRecebimentoB, *vetorRecebimentoDiag;
        int **matrix;
        int *vetEnvioMatrix,*vetorB, *diagPrincipal;

        if(my_rank==0)
        {
            MPI_Status status;
            MPI_Recv(&linhaEscolhida, 1, MPI_INT, 0, 1, parentcomm, &status); //Recebimento da linha de testagem.
            //*****************************GERANDO MATRIZ**************************************
            srand(64591);
            vetorB = (int*)malloc(orderOfMatrix*sizeof(int));
            diagPrincipal = (int*)malloc(orderOfMatrix*sizeof(int));
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
                    //matrix[i][i]=0;
                    diagPrincipal[i] = matrix[i][i];
                    if(i==linhaEscolhida)
                        linhaTeste[j]=matrix[i][j];
                }
                vetorB[i] = (rand() % ( (int)0.1*somaProvisoria[i] - RAND_LOWER_BOUND + 1)) + RAND_LOWER_BOUND;
            }

           /*  matrix[0][0] =4; //Diagonal Princ (geracao hardcoded pra testes.)
		    matrix[0][1] =2;
		    matrix[0][2] =1;
            // matrix[0][3] =0;
            matrix[1][0] =1;
            matrix[1][1] =3; //Diagonal Princ
            matrix[1][2] =1;
        // matrix[1][3] =0;
            matrix[2][0] =2;
            matrix[2][1] =3;
            matrix[2][2] =6; //Diagonal Princ
        // matrix[2][3] =0;
        // matrix[3][0] =5;
        // matrix[3][1] =5;
            //matrix[3][2] =5;
            //matrix[3][3] =16; //Diagonal Princ

            //int vetorB[4] = {7, -8, 6, 5};
            vetorB[0]=7;
            vetorB[1]=-8;
            vetorB[2]=6;
            //vetorB[3]=5;

            diagPrincipal[0]=4;
            diagPrincipal[1]=3;
            diagPrincipal[2]=6; */

            j=0;
            k=0;
            for(i=0;i<orderOfMatrix*orderOfMatrix;i++)
            {
                if(j==orderOfMatrix){
                    j=0;
                    k++;
                }
                vetEnvioMatrix[i] = matrix[k][j];
                j++;
            }
            //*****************************GERANDO MATRIZ**************************************  
        }
        vetorRecebLinhas = (int*)malloc(linhasPorProcesso*orderOfMatrix*sizeof(int));
        vetorRecebimentoB = (int*)malloc(linhasPorProcesso*sizeof(int));
        vetorRecebimentoDiag =(int*)malloc(linhasPorProcesso*sizeof(int));
        MPI_Scatter(vetEnvioMatrix, linhasPorProcesso*orderOfMatrix, MPI_INT, vetorRecebLinhas, linhasPorProcesso*orderOfMatrix, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Scatter(vetorB, linhasPorProcesso, MPI_INT, vetorRecebimentoB, linhasPorProcesso, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Scatter(diagPrincipal, linhasPorProcesso, MPI_INT, vetorRecebimentoDiag, linhasPorProcesso, MPI_INT, 0, MPI_COMM_WORLD);
        j=0;
        k=0;

        //Critério de convergência
        int criterioLinha=1;
        int *somaLinha;
        somaLinha = malloc(linhasPorProcesso*sizeof(int));
        for(i=0;i<linhasPorProcesso;i++){
            somaLinha[i]=0;
        }
        int w = my_rank*linhasPorProcesso;
        int elementoDiagonalAtual;
        j=0;k=0;
        #pragma omp parallel for private(j,k) num_threads(numberOfThreads)
        for(i=0;i<linhasPorProcesso*orderOfMatrix;i++)
        {
            if(criterioLinha) //verificar se da certo msm
                continue;
            if(j==orderOfMatrix){
                    j=0;
                    k++;
            }
            if(j != w+k){
                somaLinha[k] = somaLinha[k] + vetorRecebLinhas[i];
            }
            if(somaLinha[k] > vetorRecebimentoDiag[k]){
                criterioLinha=0;
                printf("\n Falhou no critério, na linha %d",k+w);
            }
            j++;
        }
        int reducaoCriterio;
        MPI_Reduce(&criterioLinha, &reducaoCriterio, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Bcast(&reducaoCriterio, 1, MPI_INT, 0, MPI_COMM_WORLD); //Envia a todos os processos se o criterio foi bem sucedido, se sim continua,
        if(reducaoCriterio){
            if(my_rank==0){
                printf("Passou pelo criterio das linhas.\n");
            }
        }else{
            fflush(stdout);
            MPI_Finalize();
            return 0;
        }
        float *vetorRespostaInicial;
        vetorRespostaInicial = (float*)malloc(orderOfMatrix*sizeof(float));
        if(my_rank==0){
            for(i=0;i<orderOfMatrix;i++){
                vetorRespostaInicial[i] = (float) vetorB[i]/(float) diagPrincipal[i];
            }
        }
        int numIter=0;
        MPI_Bcast(vetorRespostaInicial, orderOfMatrix, MPI_FLOAT, 0, MPI_COMM_WORLD);
        float *resultadosAtuais;
        resultadosAtuais = (float*)malloc(linhasPorProcesso*sizeof(float));
        float maximoDiff, maximoValor;
        float maximoDiffReduzido, maximoValorReduzido;
        w = my_rank*linhasPorProcesso;
        if(my_rank==0){
            printf("\n Inicio das Iterações \n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
        do {
            maximoDiff = 0;
            maximoValor = 0;
            for(i=0;i<linhasPorProcesso;i++){ //Não era pra ser nao, era pra cad aprocesso ter o seu valor independente.
                resultadosAtuais[i]=0;
                for(j=0;j<orderOfMatrix;j++){
                    if(vetorRecebLinhas[i*orderOfMatrix+j] == vetorRecebimentoDiag[i]){
                        resultadosAtuais[i] += ((float)vetorRecebimentoB[i]/(float)vetorRecebLinhas[i*orderOfMatrix + j]);
                    }else{
                        resultadosAtuais[i] -=  ((float) vetorRecebLinhas[i*orderOfMatrix + j] * vetorRespostaInicial[j] / (float) vetorRecebimentoDiag[i]);
                    }
                }
                if(fabs(resultadosAtuais[i]) > maximoValor){
                    maximoValor = fabs(resultadosAtuais[i]);
                }
                if(fabs(resultadosAtuais[i] - vetorRespostaInicial[w+i]) > maximoDiff){
                    maximoDiff = fabs(resultadosAtuais[i] - vetorRespostaInicial[w+i]);
                }
                if(my_rank==0){
                    numIter++;
                }
            }
            MPI_Allreduce(&maximoValor, &maximoValorReduzido, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(&maximoDiff, &maximoDiffReduzido, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allgather(resultadosAtuais, linhasPorProcesso, MPI_FLOAT, vetorRespostaInicial, linhasPorProcesso, MPI_FLOAT, MPI_COMM_WORLD);
        }while(maximoDiffReduzido/maximoValorReduzido >= 0.0015);
        //saindo do do while significa que convergiu,podemos testar a resposta contra a linha escolhida la no inicio
        MPI_Barrier(MPI_COMM_WORLD);
        if(my_rank ==0){ //tendeu
            float resultadoFinal=0;
            if(linhaEscolhida > orderOfMatrix){
                linhaEscolhida = 0;
                printf("\nLinha escolhida fora dos limites, usando a linha 0\n");
            }
            for(i=0;i<orderOfMatrix;i++){
                resultadoFinal += matrix[linhaEscolhida][i] * vetorRespostaInicial[i];
            }
            printf("\nTemos: %.4f = %d",resultadoFinal,vetorB[linhaEscolhida]);
            printf("\nConvergiu em %d iteracoes", numIter);
        }
    }// Fim do else que separa os spawns do root process
    fflush(stdout);
    MPI_Finalize();
    return 0;
}