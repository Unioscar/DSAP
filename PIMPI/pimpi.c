#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(int argc, char** argv){

        int myrank,numprocs,n,indice;
        int intervalos;
        long double res = 3.141592653589793238462643;
        double start_time, end_time, paralelo, secuencial;
        int desde,hasta;
        long double x,h,pi = 0,sol = 0;
        MPI_Status estado;

        MPI_Init(&argc,&argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

        if(myrank == 0){
                printf("Introduce el numero de intervalos: \n");
                scanf("%d",&intervalos);
                start_time = MPI_Wtime();
                n = intervalos / numprocs;
                indice = n + intervalos % numprocs;
                h = 1 / (double)intervalos;

        }
        MPI_Bcast(&h,1,MPI_LONG_DOUBLE,0,MPI_COMM_WORLD);
        if(myrank == 0){
                desde = indice;
                for(int i = 1; i < numprocs; i++){
                        hasta = desde + n;
                        MPI_Send(&desde,1,MPI_INT,i,i,MPI_COMM_WORLD);
                        MPI_Send(&hasta,1,MPI_INT,i,i,MPI_COMM_WORLD);
                        desde += n;
                }
                for(int i = 1; i <= indice; i++){
                        x = (i - 1)* h;
                        pi = pi + h*(4/(1 + pow(x,2)) + 4/(1 + pow(x+h,2)))/2;
                }
        }
        else{
                MPI_Recv(&desde,1,MPI_INT,0,myrank,MPI_COMM_WORLD,&estado);
                MPI_Recv(&hasta,1,MPI_INT,0,myrank,MPI_COMM_WORLD,&estado);
                for(int i = desde+1; i <= hasta;i++){
                        x = (i - 1)*h;
                        pi = pi + h*(4/(1 + pow(x,2)) + 4/(1 + pow(x+h,2)))/2;
                }
        }
        MPI_Reduce(&pi,&sol,1,MPI_LONG_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

        if(myrank == 0){
                printf("El resultado en paralelo es: %0.24Lf \n",sol);
                end_time = MPI_Wtime();
                paralelo = end_time - start_time;
                printf("El error en paralelo es de: %0.24Lf \n", sol - res);
                printf("Tiempo de ejecucion en paralelo= %f seg \n",paralelo);
        }
        if(myrank == 0){
                start_time = MPI_Wtime();
                h = 1 / (double)intervalos;
                pi = 0;
                for(int i = 1; i <= intervalos; i++){
                        x = (i-1)*h;
                        pi = pi + h* (4/(1 + pow(x,2)) + 4/(1 + pow(x+h,2)))/2;
                }
                printf("El resultado en secuencial es: %0.24Lf \n",pi);
                end_time = MPI_Wtime();
                secuencial = end_time - start_time;
                printf("El error en secuencial es de: %0.24Lf \n",pi - res);
                printf("Tiempo de ejecución en secuencial= %f seg \n",secuencial);
                double sp = secuencial / paralelo;
                printf("Speed-up= %f \n",sp);
                printf("Eficiencia= %f \n",sp/(double)numprocs);
        }
        MPI_Finalize();

}