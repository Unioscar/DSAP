#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<mpi.h>
#include <string.h>

#define nmax 1000

int main(int argc, char *argv[]) {
  float **Crear_matriz_pesos_consecutivo(int, int);
  int **Crear_matriz_caminos_consecutivo(int, int);
  double ctimer(void);
  void printMatrizCaminos(int **, int, int);
  void printMatrizPesos(float **, int, int);
  void calcula_camino(float **, int **, int);
  void Definir_Grafo(int, float **, int **);
  int i, j, k, n=5;
  double t1,t2;
  int sender;
  float **dist;
  float *aux;
  int **caminos;
  int *auxC;
  MPI_Status estado;

  int myrank,numprocs,resto,restoPadre,slice,lm;

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  if(myrank == 0){

    printf("Introduce el número de vertices en el grafo: \n");
    scanf("%d",&n);
    if(n > nmax){
        printf("El número de vertices sobrepasa el maximo \n");
        return 0;
    }

    dist = Crear_matriz_pesos_consecutivo(n,n);
    caminos = Crear_matriz_caminos_consecutivo(n,n);
    Definir_Grafo(n,dist,caminos);
    if (n <= 10) {
        printMatrizPesos(dist,n,n);
        printMatrizCaminos(caminos,n,n);
    }

    resto = n % numprocs;
    restoPadre = resto;
    for(int pos = 1; pos < numprocs;pos++){
        MPI_Send(&restoPadre,1,MPI_INT,pos,8,MPI_COMM_WORLD);
        MPI_Send(&n,1,MPI_INT,pos,8,MPI_COMM_WORLD);
    }
  }
  else{
    resto = 0;
    MPI_Recv(&restoPadre,1,MPI_INT,0,8,MPI_COMM_WORLD,&estado);
    MPI_Recv(&n,1,MPI_INT,0,8,MPI_COMM_WORLD,&estado);
  }
  slice = n / numprocs;
  lm = slice + resto;
  if(myrank == 0){
    MPI_Scatter(&dist[resto][0],slice*n,MPI_FLOAT,MPI_IN_PLACE,slice*n,MPI_FLOAT,0,MPI_COMM_WORLD);
    MPI_Scatter(&caminos[resto][0],slice*n,MPI_INT,MPI_IN_PLACE,slice*n,MPI_INT,0,MPI_COMM_WORLD);
 }
  else{
   // printf("Valor de n: %d \n",n);
   // printf("Valor de lm: %d \n", lm);
    dist= Crear_matriz_pesos_consecutivo(lm,n);
    caminos = Crear_matriz_caminos_consecutivo(lm,n);
    MPI_Scatter(&dist[resto][0],slice*n,MPI_FLOAT,&dist[resto][0],slice*n,MPI_FLOAT,0,MPI_COMM_WORLD);
    MPI_Scatter(&caminos[resto][0],slice*n,MPI_INT,&caminos[resto][0],slice*n,MPI_INT,0,MPI_COMM_WORLD);
  }

        //printf("Vamos a entrar al bucle \n");
  for (k=0; k<n; k++) {
      //printf("voy a mirar las condición de entrada \n");
      for(int m = 0; m < numprocs; m++){
        sender = 0;
        //printf("restoPadre: %d \n",restoPadre);
        if(restoPadre+ slice*(m) < k && k < restoPadre + slice*(m)){
          //  printf("He entrado para encontrar el sender soy k = %d \n",k);
         //   printf("El sender toma el valor de: %d \n",m);
            sender = m;
        }
     }

            aux = malloc(n * sizeof(float));
            auxC = malloc(n * sizeof(int));
        if(myrank == sender){
            for(int pos = 0; pos < n;pos++){
                aux[pos] = dist[k][pos];
        //      printf("valores de aux: %f \n",aux[pos]);
                auxC[pos] = caminos[k][pos];
        //      printf("valores de auxC: %d \n",auxC[pos]);
            }
        }

      //printf("Soy %d voy a hacer el broadcast. \n", sender);
      MPI_Bcast(&aux[0],n,MPI_FLOAT,sender,MPI_COMM_WORLD);
      MPI_Bcast(&auxC[0],n,MPI_INT,sender,MPI_COMM_WORLD);

      for (i = 0; i < resto + slice; i++){
          for (j = 0; j < n; j++){
                //printf("caminos[%d][%d] = %d \n",i,j,caminos[i][j]);
              if ((dist[i][k] * aux[j] != 0) ) {
                 if ((dist[i][k] + aux[j] < dist[i][j]) || (dist[i][j] == 0)){
                    dist[i][j] = dist[i][k] + aux[j];
                    caminos[i][j] = auxC[j];
                //      printf("caminos[%d][%d] = %d \n",i,j,caminos[i][j]);
                 }
              }
          }
      }
  }
  if(myrank == 0){
      MPI_Gather(MPI_IN_PLACE,slice*n,MPI_FLOAT,&dist[resto][0],slice*n,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Gather(MPI_IN_PLACE,slice*n,MPI_INT,&caminos[resto][0],slice*n,MPI_INT,0,MPI_COMM_WORLD);
  }
  else{
      MPI_Gather(&dist[resto][0],slice*n,MPI_FLOAT,&dist[resto][0],slice*n,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Gather(&caminos[resto][0],slice*n,MPI_INT,&caminos[resto][0],slice*n,MPI_INT,0,MPI_COMM_WORLD);
  }

  if(myrank == 0){
    if (n <= 10) {
        printMatrizPesos(dist,n,n);
        printMatrizCaminos(caminos,n,n);
    }
    calcula_camino(dist, caminos, n);
  }
  free(aux);
  free(auxC);
  free(caminos);
  free(dist);
  MPI_Finalize();

}

void Definir_Grafo(int n,float **dist,int **caminos)
{
// Inicializamos la matriz de pesos y la de caminos para el algoritmos de Floyd-Warshall.
// Un 0 en la matriz de pesos indica que no hay arco.
// Para la matriz de caminos se supone que los vertices se numeran de 1 a n.
  int i,j;
  for (i = 0; i < n; ++i) {
      for (j = 0; j < n; ++j) {
          if (i==j)
             dist[i][j]=0;
          else {
             dist[i][j]= (11.0 * rand() / ( RAND_MAX + 1.0 )); // aleatorios 0 <= dist < 11
             dist[i][j] = ((double)((int)(dist[i][j]*10)))/10; // truncamos a 1 decimal
             if (dist[i][j] < 2) dist[i][j]=0; // establecemos algunos a 0
          }
          if (dist[i][j] != 0)
             caminos[i][j] = i+1;
          else
             caminos[i][j] = 0;
      }
  }
}

void calcula_camino(float **a, int **b, int n)
{
 int i,count=2, count2;
 int anterior;
 int *camino;
 int inicio=-1, fin=-1;

 while ((inicio < 0) || (inicio >n) || (fin < 0) || (fin > n)) {
    printf("Vertices inicio y final: (0 0 para salir) \n");
    scanf("%d %d",&inicio, &fin);
 }
 while ((inicio != 0) && (fin != 0)) {
    anterior = fin;
    while (b[inicio-1][anterior-1] != inicio) {
       anterior = b[inicio-1][anterior-1];
       count = count + 1;
    }
    count2 = count;
    camino = malloc(count * sizeof(int));
    anterior = fin;
    camino[count-1]=fin;
    while (b[inicio-1][anterior-1] != inicio) {
       anterior = b[inicio-1][anterior-1];
       count = count - 1;
       camino[count-1]=anterior;
    }
    camino[0] = inicio;
    printf("\nCamino mas corto de %d a %d:\n", inicio, fin);
    printf("          Peso: %5.1f\n", a[inicio-1][fin-1]);
    printf("        Camino: ");
    for (i=0; i<count2; i++) printf("%d ",camino[i]);
    printf("\n");
    free(camino);
    inicio = -1;
    while ((inicio < 0) || (inicio >n) || (fin < 0) || (fin > n)) {
       printf("Vertices inicio y final: (0 0 para salir) \n");
       scanf("%d %d",&inicio, &fin);
    }

 }
}


float **Crear_matriz_pesos_consecutivo(int Filas, int Columnas)
{
// crea un array de 2 dimensiones en posiciones contiguas de memoria
 float *mem_matriz;
 float **matriz;
 int fila, col;
//printf("Numero de filas: %d \n",Filas);
//printf("Numero de columnas: %d \n",Columnas);
 if (Filas <=0)
    {
        printf("El numero de filas debe ser mayor que cero\n");
        return 0;
    }
 if (Columnas <= 0)
    {
        printf("El numero de filas debe ser mayor que cero\n");
        return 0;
    }
 mem_matriz = malloc(Filas * Columnas * sizeof(float));
 if (mem_matriz == NULL)
        {
                printf("Insuficiente espacio de memoria\n");
                return 0;
        }
 matriz = malloc(Filas * sizeof(float *));
 if (matriz == NULL)
        {
                printf ("Insuficiente espacio de memoria\n");
                return 0;
        }
 for (fila=0; fila<Filas; fila++)
    matriz[fila] = mem_matriz + (fila*Columnas);
 return matriz;
}

int **Crear_matriz_caminos_consecutivo(int Filas, int Columnas)
{
// crea un array de 2 dimensiones en posiciones contiguas de memoria
 int *mem_matriz;
 int **matriz;
 int fila, col;
 if (Filas <=0)
    {
        printf("El numero de filas debe ser mayor que cero\n");
        return 0;
    }
 if (Columnas <= 0)
    {
        printf("El numero de filas debe ser mayor que cero\n");
        return 0;
    }
 mem_matriz = malloc(Filas * Columnas * sizeof(int));
 if (mem_matriz == NULL)
        {
                printf("Insuficiente espacio de memoria\n");
                return 0;
        }
 matriz = malloc(Filas * sizeof(int *));
 if (matriz == NULL)
        {
                printf ("Insuficiente espacio de memoria\n");
                return 0;
        }
 for (fila=0; fila<Filas; fila++)
    matriz[fila] = mem_matriz + (fila*Columnas);
 return matriz;
}

void printMatrizCaminos(int **a, int fila, int col) {
        int i, j;
        char buffer[10];
        printf("     ");
        for (i = 0; i < col; ++i){
                j=sprintf(buffer, "%c%d",'V',i+1 );
                printf("%5s", buffer);
       }
        printf("\n");
        for (i = 0; i < fila; ++i) {
                j=sprintf(buffer, "%c%d",'V',i+1 );
                printf("%5s", buffer);
                for (j = 0; j < col; ++j)
                        printf("%5d", a[i][j]);
                printf("\n");
        }
        printf("\n");
}

void printMatrizPesos(float **a, int fila, int col) {
        int i, j;
        char buffer[10];
        printf("     ");
        for (i = 0; i < col; ++i){
                j=sprintf(buffer, "%c%d",'V',i+1 );
                printf("%5s", buffer);
       }
        printf("\n");
        for (i = 0; i < fila; ++i) {
                j=sprintf(buffer, "%c%d",'V',i+1 );
                printf("%5s", buffer);
                for (j = 0; j < col; ++j)
                        printf("%5.1f", a[i][j]);
                printf("\n");
        }
        printf("\n");
}