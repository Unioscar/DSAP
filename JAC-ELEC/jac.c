#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

const int SMAX = 200;
const int MAXPROC = 8;
const double cota = 0.01;

double **CrearArray2D_double_consecutivo(int , int );
void inicializar(double **,int ,int ,double ,double ,double ,double ,double ,int ,int );
void printMatriz(double **, int , int );

int main(int argc,char** argv){
    int myrank,numprocs,s,iterMax,inner,count=0;
    double **x, **xold, **xant, **xinit;
    double vi,fi,fd,fa,fb,sum,norm2=1.0;
    int resto,slocal,slice;
    int *a_chunk_sizes;
    int *a_despla;
    MPI_Status estado;
    FILE *RES;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

    a_chunk_sizes=(int *)malloc(numprocs*sizeof(int));
    a_despla=(int *)malloc(numprocs*sizeof(int));

    if(numprocs < MAXPROC){

        if(myrank == 0){
            do{
                printf("Introduce el valor de s: \n");
                scanf("%d",&s);

                if(s > SMAX){
                    printf("El valor de s no esta dentro de los limites permitidos \n");
                }
            } while(s > SMAX);
            
            printf("Introduce el valor de iterMax: \n");
            scanf("%d",&iterMax);

            printf("Introduce el valor de inner: \n");
            scanf("%d",&inner);

            resto = s %numprocs;

            vi=1.0;
            fi=100.0;
            fd=100.0;
            fa=100.0;
            fb=0.0;

            x = CrearArray2D_double_consecutivo(s+2,s+2);
            xold = CrearArray2D_double_consecutivo(s+2,s+2);
            xant = CrearArray2D_double_consecutivo(s+2,s+2);
            xinit = CrearArray2D_double_consecutivo(s+2,s+2);
        
            inicializar(xinit,s,s,vi,fi,fd,fa,fb,0,1);
            inicializar(xant,s,s,vi,fi,fd,fa,fb,0,1);
            inicializar(xold,s,s,vi,fi,fd,fa,fb,0,1);
        }
        else{

            resto = 0;

            x = CrearArray2D_double_consecutivo(s+2,s+2);
            xold = CrearArray2D_double_consecutivo(s+2,s+2);
            xant = CrearArray2D_double_consecutivo(s+2,s+2);
            xinit = CrearArray2D_double_consecutivo(s+2,s+2);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Bcast(&s,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&iterMax,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&inner,1,MPI_INT,0,MPI_COMM_WORLD);

        slice = s / numprocs;
        slocal = slice + resto;

        for(int i = 0 ; i < numprocs ;i++) {
            a_chunk_sizes[i] = (slocal + 2) * (s + 2);
            a_despla[i] = i * (slice + 1) * (s + 2);
        }

        a_despla[0] = 0;


        MPI_Scatterv(&xinit[0][0],a_chunk_sizes,&a_despla[0],MPI_DOUBLE,&x[0][0],a_chunk_sizes[myrank],MPI_DOUBLE,0,MPI_COMM_WORLD);

        while ((norm2 > cota) && (count < iterMax)){
            count++;
            for(int i = 1 ; i <= slocal; i++){
                for(int j = 1; j <= slocal; j++){
                    xold[i][j] = x[i][j];
                }
            }

            for(int iter = 0; iter < inner; iter++){

                for(int i = 1; i <= slocal; i++){
                    for(int j = 1; j <= slocal; j++){
                        xant[i][j] = x[i][j];
                    }
                }
                for(int i = 1; i <= slocal; i++){
                    for(int j = 1; j <= slocal; j++){
                        x[i][j] = 0.25 * (xant[i+1][j] + xant[i-1][j] + xant[i][j+1] + xant[i][j-1]);
                    }
                }             
            
            }

            for(int i = 1; i <= slocal; i++){
                for(int j = 1; j <= slocal; i++){
                    sum = sum + (x[i][j]-xant[i][j])*(x[i][j]-xant[i][j]);
                }
            }
            norm2 = sqrt(sum);

            MPI_Allreduce(MPI_IN_PLACE,&norm2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

            if(myrank == 0)
            {
                MPI_Send(&x[slocal+1][0],s+2,MPI_DOUBLE,1,0,MPI_COMM_WORLD);
                MPI_Recv(&x[slocal+1][0],s+2,MPI_DOUBLE,1,0,MPI_COMM_WORLD,&estado);
            }
            else if(myrank == numprocs - 1)
            {
                MPI_Send(&x[slocal-1][0],s+2,MPI_DOUBLE,myrank-1,0,MPI_COMM_WORLD);
                MPI_Recv(&x[slocal-1][0],s+2,MPI_DOUBLE,myrank-1,0,MPI_COMM_WORLD,&estado);
            }
            else
            {
                MPI_Send(&x[slocal+1][0],s+2,MPI_DOUBLE,myrank+1,0,MPI_COMM_WORLD);
                MPI_Recv(&x[slocal+1][0],s+2,MPI_DOUBLE,myrank+1,0,MPI_COMM_WORLD,&estado);

                MPI_Send(&x[slocal-1][0],s+2,MPI_DOUBLE,myrank-1,0,MPI_COMM_WORLD);
                MPI_Recv(&x[slocal-1][0],s+2,MPI_DOUBLE,myrank-1,0,MPI_COMM_WORLD,&estado);
            }
        }

        RES = fopen("res.m","w");
        fprintf(RES,"mx=%d;\n",s);
        fprintf(RES,"sol=[\n");
        MPI_Send(&x[1][1],slocal*s,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
        if(myrank == numprocs - 1)
        {
            MPI_Send(&x[1][1],(slocal+1)*s,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
        }
        if(myrank == 0)
        {
            for(int iter = 0; iter < numprocs; iter++)
            {
                for (int i=0; i<=slocal+1; i++) {
                    for (int j=0; j<=s+1; j++)
                        fprintf(RES,"%f ",x[i][j]);
                    fprintf(RES,"\n\n");
                }
                if(iter == numprocs -1)
                {
                    MPI_Recv(&x[1][1],(slice+1)*s,MPI_DOUBLE,iter,0,MPI_COMM_WORLD,&estado);
                }
                else
                {
                    MPI_Recv(&x[1][1],slice*s,MPI_DOUBLE,iter,0,MPI_COMM_WORLD,&estado);
                }
                
            }

            fprintf(RES, "];");
            fclose(RES);
                
        }

    
        free(x);
        free(xold);
        free(xinit);
        free(xant);

        MPI_Finalize();
        return 0;
    }
}

double **CrearArray2D_double_consecutivo(int Filas, int Columnas)
{
// crea un array de 2 dimensiones en posiciones contiguas de memoria
 double *mem_matriz;
 double **matriz;
 int fila, col;
 if (Filas <=0)
    {
        printf("El numero de filas debe ser mayor que cero\n");
//        return;
    }
 if (Columnas <= 0)
    {
        printf("El numero de filas debe ser mayor que cero\n");
//        return;
    }
 mem_matriz = malloc(Filas * Columnas * sizeof(double));
 if (mem_matriz == NULL)
        {
                printf("Insuficiente espacio de memoria\n");
  //              return;
        }
 matriz = malloc(Filas * sizeof(double *));
 if (matriz == NULL)
        {
                printf ("Insuficiente espacio de memoria\n");
//                return;
        }
 for (fila=0; fila<Filas; fila++)
    matriz[fila] = mem_matriz + (fila*Columnas);
 return matriz;
}

void inicializar(double **u,int nf,int nc,double vi,double fi,double fd,double fa,double fb,int inum,int nproc)
{
/*         U=fa
        __________                _________
       |          |              |   0     |
       |Inicial   |              |_________|
  U=fi |U=vi      |U=fd     inum=|   1     | 
       |          |              |_________|
       |          |              |   2     |
       |__________|              |_________|
           U=fb
        nf : filas en U; numero de filas en el array U sin la frontera
        nc : columnas en U; numero de columnas en el array U sin la frontera
        vi : valor del iterado inicial
        fi : valor de U en la frontera izquierda
        fd : valor de U en la frontera derecha
        fa : valor de U en la frontera superior
        fb : valor de U en la frontera inferior
       inum: numero del proceso
      nproc: numero total de procesos
 Output: U : array de dimension (0:nf+1,0:nc+1) que contine:
             -- puntos frontera seleccionados
             -- y vector inicial: U(1:nf,1:nc) */

 int i,j;
// Valor inicial
 for (i=0; i<=nf+1; i++)
    for (j=1; j<=nc; j++)
           u[i][j]=vi;
//  Condiciones de contorno superiores
      if (inum == 0) {
         for (j=1; j<=nc; j++)
             u[0][j]=fa;
      }
// Condiciones de contorno inferiores
      if (inum == nproc-1) {
         for (j=1; j<=nc; j++) 
             u[nf+1][j] = fb;
     }
// Condiciones de contorno izquierda y derecha
   for (i=0; i<=nf+1; i++) {
         u[i][0] = fi;
         u[i][nc+1] = fd;
   }
      return;
}

void printMatriz(double **a, int fila, int col) {
        int i, j;
        char buffer[10];
        printf("      ");
        for (i = 0; i < col; ++i){
                j=sprintf(buffer, "%d",i );
                printf("%6s", buffer);
       }
        printf("\n");
        for (i = 0; i < fila; ++i) {
                j=sprintf(buffer, "%d",i );
                printf("%6s", buffer);
                for (j = 0; j < col; ++j)
                        printf("%6.1f", a[i][j]);
                printf("\n");
        }
        printf("\n");
}

