#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "operaciones_basicas.h"
#include "leer_matriz.h"
#include "metodos_lineales.h"
#include "errores.h"
#include "tools.h"
#include "io_tools.h"

int main() {
    MatrixCOO* A = NULL;
    double tiempoLectura = 0.0;
    char archivoMatriz[256];

    printf("Introduce el nombre del archivo de la matriz (.txt): ");
    if (scanf("%s", archivoMatriz) != 1) return 1;

    A = leerMatrizArchivo(archivoMatriz, &tiempoLectura);
    int n = A->n;
    
    double densidad = (double)A->nnz / (double)(A->n * A->m);
    // Memoria ocupada estimada en KB
    double memoriaKB = (sizeof(MatrixCOO) + A->nnz * (sizeof(int)*2 + sizeof(double))) / 1024.0;

    double* b = (double*)malloc(n * sizeof(double));
    double* x = (double*)calloc(n, sizeof(double));

    int opVector;
    printf("\nSelecciona el vector:\n");
    printf("1. Usar vector de 1s (1, 1, 1...)\n");
    printf("2. Usar vector personalizado\n");
    scanf("%d", &opVector);

    if (opVector == 1) {
        for (int i = 0; i < n; i++) {
            b[i] = 1.0;
        }
    } else {
        int opVecTipo;
        printf("\nSelecciona la entrada del vector:\n");
        printf("1. Leer desde archivo .txt\n");
        printf("2. Introducir elementos por consola separados por comas\n");
        scanf("%d", &opVecTipo);

        if (opVecTipo == 1) {
            char archivoVector[256];
            printf("Introduce el nombre del archivo del vector: ");
            scanf("%s", archivoVector);
            FILE* fv = fopen(archivoVector, "r");
            if (fv) {
                for (int i = 0; i < n; i++) {
                    if (fscanf(fv, "%lf", &b[i]) != 1) {
                        b[i] = 1.0; // Fallback
                    }
                }
                fclose(fv);
            } else {
                printf("No se pudo abrir el archivo, utilizando vector de 1s por defecto.\n");
                for (int i = 0; i < n; i++) b[i] = 1.0;
            }
        } else {
            printf("Introduce los %d elementos separados por comas: ", n);
            char buffer[4096];
            scanf("%s", buffer);
            char* token = strtok(buffer, ",");
            int i = 0;
            while (token != NULL && i < n) {
                b[i] = atof(token);
                token = strtok(NULL, ",");
                i++;
            }
        }
    }

    int opMetodo;
    printf("\nSelecciona el metodo a calcular:\n");
    printf("1. Factorización PLU\n");
    printf("2. Método de Cholesky\n");
    printf("3. Gradiente Conjugado\n");
    printf("4. Todos\n");
    scanf("%d", &opMetodo);

    int opSalida;
    printf("\nSelecciona la opcion de salida:\n");
    printf("1. Imprimir archivo txt o csv con resultados\n");
    printf("2. Imprimir estadisticos en CSV\n");
    printf("3. Imprimir en pantalla solo el tiempo y memoria total usados\n");
    printf("4. Imprimir resultados txt y estadisticos csv\n");
    scanf("%d", &opSalida);

    // Abrir o crear el archivo estadístico (en modo append para sumar filas)
    FILE* csv = fopen("estadisticos_libreriaCOO.csv", "a");
    if (csv) {
        fseek(csv, 0, SEEK_END);
        if (ftell(csv) == 0) {
            fprintf(csv, "Metodo,N,NNZ,MemoriaKB,Densidad,TiempoLectura,TiempoMetodo,TiempoTotal,Iteraciones,Residual\n");
        }
    }

    if (opMetodo == 1 || opMetodo == 4) {
        double tiempoPLU = 0.0;
        int iterPLU = 0;
        double* x_plu = (double*)calloc(n, sizeof(double));

        metodoPLU(A, b, x_plu, &iterPLU, &tiempoPLU);

        double res = calcularResidualFinal(A, b, x_plu);
        double tiempoTotal = tiempoLectura + tiempoPLU;

        if (csv) {
            fprintf(csv, "PLU,%d,%d,%.2f,%.6f,%.6f,%.6f,%.6f,%d,%.2e\n",
                    n, A->nnz, memoriaKB, densidad, tiempoLectura, tiempoPLU, tiempoTotal, iterPLU, res);
        }

        if (opSalida == 1 || opSalida == 4) {
            int formato;
            printf("Selecciona el formato de salida para PLU (1: TXT, 2: CSV): ");
            scanf("%d", &formato);
            if (formato == 1) {
                guardarTxtDetallado("PLU", archivoMatriz, n, x_plu);
            } else {
                guardarCsvResultados("PLU", archivoMatriz, n, x_plu);
            }
        }
        if (opSalida == 3) {
            printf("[PLU] Tiempo Total: %.6f s | Memoria: %.2f KB\n", tiempoTotal, memoriaKB);
        }
        free(x_plu);
    }

    if (opMetodo == 2 || opMetodo == 4) {
        double tiempoCholesky = 0.0;
        int iterCholesky = 0;
        double* x_chol = (double*)calloc(n, sizeof(double));

        metodoCholesky(A, b, x_chol, &iterCholesky, &tiempoCholesky);

        double res = calcularResidualFinal(A, b, x_chol);
        double tiempoTotal = tiempoLectura + tiempoCholesky;

        if (csv) {
            fprintf(csv, "Cholesky,%d,%d,%.2f,%.6f,%.6f,%.6f,%.6f,%d,%.2e\n",
                    n, A->nnz, memoriaKB, densidad, tiempoLectura, tiempoCholesky, tiempoTotal, iterCholesky, res);
        }

        if (opSalida == 1 || opSalida == 4) {
            int formato;
            printf("Selecciona el formato de salida para Cholesky (1: TXT, 2: CSV): ");
            scanf("%d", &formato);
            if (formato == 1) {
                guardarTxtDetallado("Cholesky", archivoMatriz, n, x_chol);
            } else {
                guardarCsvResultados("Cholesky", archivoMatriz, n, x_chol);
            }
        }
        if (opSalida == 3) {
            printf("[Cholesky] Tiempo Total: %.6f s | Memoria: %.2f KB\n", tiempoTotal, memoriaKB);
        }
        free(x_chol);
    }

    if (opMetodo == 3 || opMetodo == 4) {
        double tiempoGC = 0.0;
        int iterGC = 0;
        double* x_gc = (double*)calloc(n, sizeof(double));

        gradienteConjugado(A, b, x_gc, &iterGC, &tiempoGC);

        double res = calcularResidualFinal(A, b, x_gc);
        double tiempoTotal = tiempoLectura + tiempoGC;

        if (csv) {
            fprintf(csv, "GradienteConjugado,%d,%d,%.2f,%.6f,%.6f,%.6f,%.6f,%d,%.2e\n",
                    n, A->nnz, memoriaKB, densidad, tiempoLectura, tiempoGC, tiempoTotal, iterGC, res);
        }

        if (opSalida == 1 || opSalida == 4) {
            int formato;
            printf("Selecciona el formato de salida para Gradiente Conjugado (1: TXT, 2: CSV): ");
            scanf("%d", &formato);
            if (formato == 1) {
                guardarTxtDetallado("GradienteConjugado", archivoMatriz, n, x_gc);
            } else {
                guardarCsvResultados("GradienteConjugado", archivoMatriz, n, x_gc);
            }
        }
        if (opSalida == 3) {
            printf("[Gradiente Conjugado] Tiempo Total: %.6f s | Memoria: %.2f KB\n", tiempoTotal, memoriaKB);
        }
        free(x_gc);
    }

    // Cerrar el archivo estadístico
    if (csv) fclose(csv);

    free(b);
    free(x);
    liberarMatrizCOO(A);

    printf("\n[OK] Ejecucion finalizada con exito.\n");
    return 0;
}