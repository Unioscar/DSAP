#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>

const int RMAX=4;
const int MAXBLOQTAM=100;

void mult(double a[],double b[],double *c,int m)
{
    int i,j,k;
    for(i = 0; i < m ;i++)
        for(j = 0; j < m ;j++)
            for(k = 0; k < m ;k++)
                c[i*m+j]=c[i*m+j]+a[i*m+k]*b[k*m+j];
    return;
};


int main(int argc,char** argv){
    int myrank,numprocs,n,indice;
    int tam,r,bloqtam,fila,columna;
    double *a, *b, *c, *atmp;
    MPI_Status estado;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

    if(myrank == 0)
    {
        do
        {
            printf("Introduce el tamaño de la matriz: \n");
            scanf("%d",&r);
            if(r > RMAX)
            {
                printf("El valor de r sobrepasa los limites \n");
            }
        } while(r > RMAX);

        do
        {
            printf("Introduce el tamaño del bloque: \n");
            scanf("%d",&bloqtam);
            if(bloqtam > MAXBLOQTAM)
            {
                printf("El tamaño del bloque sobrepasa los limites \n");
            }
        } while(bloqtam > MAXBLOQTAM);


        if(r * r != numprocs)
        {
          printf("El numero de procesadores no es el correcto. \n");
        }
    }

    MPI_Bcast(&r,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&bloqtam,1,MPI_INT,0,MPI_COMM_WORLD);

    if( r*r != numprocs )
    {
        MPI_Finalize();
        return 0;
    }

    a = malloc(bloqtam*bloqtam * sizeof(double));
    b = malloc(bloqtam*bloqtam * sizeof(double));
    atmp = malloc(bloqtam*bloqtam * sizeof(double));
    c = malloc(bloqtam*bloqtam * sizeof(double));

    fila = myrank/r;
    columna = myrank%r;


    for(int i = 0; i < bloqtam*bloqtam ; i++){
        a[i] = i*(double)(fila*columna+1)/bloqtam;
    }

    for(int i = 0; i < bloqtam ; i++)
    {
        for(int j = 0; j < bloqtam ; j++)
        {
            if(i == j && columna == fila)
            {
                b[i*bloqtam+j] = 1.0;
            }
            else
            {
                b[i*bloqtam+j] = 0;
            }
        }
    }

    int *vectorRankFilas;
    vectorRankFilas= malloc((r-1) * sizeof(int));
    int pos = 0;

    for(int i = 0; i < numprocs; i++)
    {
        if(fila == i/r && columna != i%r)
        {
            vectorRankFilas[pos] = i;
            pos++;
        }
    }


    int top;
    int bottom;

    if(fila == 0)
    {
        top = (r-1)*r + columna;
        bottom = (fila+1) * r + columna;
    }
    else if(fila == r-1)
    {
        top = (fila-1) * r + columna;
        bottom = columna;
    }
    else
    {
        top = (fila-1)*r + columna;
        bottom = (fila+1)*r + columna;
    }


    for(int k = 0;k < r;k++)
    {
        if(columna == (fila+k)%r)
        {
                for(int i = 0 ; i < r-1 ; i++)
                {
                    MPI_Send(a,bloqtam*bloqtam,MPI_DOUBLE,vectorRankFilas[i],k+4,MPI_COMM_WORLD);
                }
           mult(a,b,c,bloqtam);
        }
        else
        {
           MPI_Recv(atmp,bloqtam*bloqtam,MPI_DOUBLE,r*fila+((fila+k)%r),k+4,MPI_COMM_WORLD,&estado);
           mult(atmp,b,c,bloqtam);
        }

        MPI_Send(b,bloqtam*bloqtam,MPI_DOUBLE,top,0,MPI_COMM_WORLD);
        MPI_Recv(b,bloqtam*bloqtam,MPI_DOUBLE,bottom,0,MPI_COMM_WORLD,&estado);

    }

    bool error = false;

    for(int i = 0 ; i < bloqtam*bloqtam ; i++)
    {
        if(a[i] != c[i])
        {
            if(a[i] != 0 && c[i]!= 0)
               error = true;
        }
    }
    if(error)
    {
        printf("lo valores del proceso %d no son los correctos \n",myrank);
    }
    else
    {
        printf("OK. Desde el proceso %d \n",myrank);
    }

    free(a);
    free(b);
    free(atmp);
    free(c);
    free(vectorRankFilas);
    MPI_Finalize();
}