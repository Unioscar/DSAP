#include <mpi.h>
#include <stdio.h>

int main(int argc,char** argv){
    int myrank, numprocs, i;
    int dato;

    MPI_Status estado;
    MPI_Init(&argc,&argv); //Iniciamos MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); //Devuelve el rango del proceso en el que nos encontramos 0...N
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs); //Devuelve el numero de procesos totales
    int desde, hacia;

    if(myrank == 0){
        //Enviamos el dato al primer proceso
        printf("Introduce un entero: \n");
        scanf("%d",&dato);
        hacia = myrank + 1;
        desde = numprocs-1;
        printf("Voy a transmitir el entero: %d \n",dato);
        MPI_Send(&dato,1,MPI_INT,hacia,8,MPI_COMM_WORLD); 
        MPI_Recv(&dato,1,MPI_INT,desde,8,MPI_COMM_WORLD,&estado);
        printf("Soy el proceso %d, el entero que he recibido es: %d \n", myrank,dato);

    }
    else if(myrank == numprocs-1){
        desde = myrank - 1;
        hacia = 0;
        MPI_Recv(&dato,1,MPI_INT,desde,8,MPI_COMM_WORLD,&estado);
        printf("Soy el proceso %d, el entero que he recibido es: %d \n", myrank,dato);
        dato += 1;
        MPI_Send(&dato,1,MPI_INT,hacia,8,MPI_COMM_WORLD);
    }
    else
    {
        desde = myrank - 1;
        hacia = myrank + 1;
        MPI_Recv(&dato,1,MPI_INT,desde,8,MPI_COMM_WORLD,&estado);
        printf("Soy el proceso %d, el entero que he recibido es: %d \n", myrank,dato);
        dato += 1;
        MPI_Send(&dato,1,MPI_INT,hacia,8,MPI_COMM_WORLD);

    }

    MPI_Finalize(); //Finalizamos el proceso

}
