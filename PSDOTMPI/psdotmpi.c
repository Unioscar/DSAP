#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>

int MAXNPROCS = 8;
int MAXN = 100000000;

int main(int argc,char** argv){
    int myrank,numprocs,n,indice;
    int tam;
    double res = 0, sol= 0;
    double *x,*y;
    MPI_Status estado;


        MPI_Init(&argc,&argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        if(numprocs <= MAXNPROCS){
            if(myrank == 0){
                do{
                    printf("Introduce el tamaño del vector: \n");
                    scanf("%d",&tam);
                    if(tam > MAXN || tam < 0){
                        printf("El tamaño no esta dentro de los rangos permitidos. \n");
                    }
                }while(tam > MAXN || tam < 0);
                x = malloc(tam * sizeof(double));
                y = malloc(tam * sizeof(double));
                for(int i = 0; i < tam; i++){
                    x[i] = 1 / (double)(i + 1);
                    y[i] = i + 1;
                }
                n = tam / numprocs;
                indice = n + tam % numprocs;
            }
                MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
            if(myrank == 0){
                for(int i = 0; i < indice; i++){
                    res += x[i] * y[i];
                }
                for(int i = 1; i < numprocs;i++){
                    MPI_Send(&x[indice],n,MPI_DOUBLE,i,i,MPI_COMM_WORLD);
                    MPI_Send(&y[indice],n,MPI_DOUBLE,i,i,MPI_COMM_WORLD);
                    indice += n;
                }
            }

            if(myrank != 0){
                x = malloc(n * sizeof(double));
                y = malloc(n * sizeof(double));

                MPI_Recv(x,n,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD,&estado);
                MPI_Recv(y,n,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD,&estado);
                for(int i = 0; i < n; i++){
                    res += x[i] * y[i];
                }
            }
            MPI_Reduce(&res,&sol,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
            if(myrank == 0){
                printf("El resultado es: %f \n",sol);
            }
        }
        MPI_Finalize();

}