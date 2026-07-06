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

    // ================= EXTRACCIÓN AUTOMÁTICA DE RALIDAD =================
    double ralidadArchivo = 0.0;
    // Buscamos el primer guion bajo '_' en el nombre del archivo
    char* primerGuion = strchr(archivoMatriz, '_');
    if (primerGuion != NULL) {
        // Leemos el numero flotante justo despues del primer guion bajo
        sscanf(primerGuion + 1, "%lf", &ralidadArchivo);
    }
    // ===================================================================

    A = leerMatrizArchivo(archivoMatriz, &tiempoLectura);
    int n = A->n;
    
    double densidad = (double)A->nnz / (double)(A->n * A->m);
    double memMatKB = (sizeof(MatrixCOO) + A->nnz * (sizeof(int)*2 + sizeof(double))) / 1024.0;

    double* b = (double*)malloc(n * sizeof(double));
    double* x = (double*)calloc(n, sizeof(double));

    int opVector;
    printf("\nSelecciona el vector:\n");
    printf("1. Usar vector de 1s (1, 1, 1...)\n");
    printf("2. Usar vector personalizado\n");
    scanf("%d", &opVector);

    if (opVector == 1) {
        for (int i = 0; i < n; i++) b[i] = 1.0;
    } else {
        printf("Introduce los %d elementos del vector b:\n", n);
        for (int i = 0; i < n; i++) {
            scanf("%lf", &b[i]);
        }
    }

    int opMetodo;
    printf("\nSelecciona el metodo a ejecutar:\n");
    printf("1. PLU\n");
    printf("2. Cholesky\n");
    printf("3. Gradiente Conjugado\n");
    printf("4. Ejecutar Todos\n");
    scanf("%d", &opMetodo);

    int opSalida;
    printf("\nSelecciona que deseas hacer con los resultados:\n");
    printf("1. Guardar resultados en archivo (TXT/CSV)\n");
    printf("2. Imprimir estadisticos en pantalla\n");
    printf("3. Guardar estadisticos en CSV general (estadisticos_libreriaCOO.csv)\n");
    printf("4. Todas las anteriores\n");
    scanf("%d", &opSalida);

    FILE* csv = NULL;
    if (opSalida == 3 || opSalida == 4) {
        csv = fopen("estadisticos_libreriaCOO.csv", "a");
        if (csv) {
            fseek(csv, 0, SEEK_END);
            if (ftell(csv) == 0) {
                // NUEVO NUEVO: Cabecera con la columna 'Ralidad' añadida despues de Densidad
                fprintf(csv, "Metodo,N,NNZ,MemMat(KB),MemMet(KB),MemTotal(KB),Densidad,Ralidad,TLectura(s),TMetodo(s),TTotal(s),Iter,Residual\n");
            }
            rewind(csv);
            fseek(csv, 0, SEEK_END);
        }
    }

    // ================= MÉTODO PLU =================
    if (opMetodo == 1 || opMetodo == 4) {
        double tiempoPLU = 0.0;
        int iterPLU = 0;
        double memoriaMetodoKB = 0.0;
        double* x_plu = (double*)calloc(n, sizeof(double));

        metodoPLU(A, b, x_plu, &iterPLU, &tiempoPLU, &memoriaMetodoKB);

        double res = calcularResidualFinal(A, b, x_plu);
        double tiempoTotal = tiempoLectura + tiempoPLU;
        double memTotalKB = memMatKB + memoriaMetodoKB;

        if (csv) {
            // Imprime la ralidadArchivo extraida de forma automatica
            fprintf(csv, "PLU,%d,%d,%.2f,%.2f,%.2f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%.2e\n",
                    n, A->nnz, memMatKB, memoriaMetodoKB, memTotalKB, densidad, ralidadArchivo, tiempoLectura, tiempoPLU, tiempoTotal, iterPLU, res);
        }

        if (opSalida == 1 || opSalida == 4) {
            int formato;
            printf("Selecciona formato para PLU (1: TXT, 2: CSV): ");
            scanf("%d", &formato);
            if (formato == 1) guardarTxtDetallado("PLU", archivoMatriz, n, x_plu);
            else guardarCsvResultados("PLU", archivoMatriz, n, x_plu);
        }
        if (opSalida == 2 || opSalida == 4) {
            printf("[PLU] Tiempo Total: %.6f s | Ralidad: %.4f | Memoria Total: %.2f KB\n", tiempoTotal, ralidadArchivo, memTotalKB);
        }
        free(x_plu);
    }

    // ================= MÉTODO CHOLESKY =================
    if (opMetodo == 2 || opMetodo == 4) {
        double tiempoChol = 0.0;
        int iterChol = 0;
        double memoriaMetodoKB = 0.0;
        double* x_chol = (double*)calloc(n, sizeof(double));

        metodoCholesky(A, b, x_chol, &iterChol, &tiempoChol, &memoriaMetodoKB);

        double res = calcularResidualFinal(A, b, x_chol);
        double tiempoTotal = tiempoLectura + tiempoChol;
        double memTotalKB = memMatKB + memoriaMetodoKB;

        if (csv) {
            fprintf(csv, "Cholesky,%d,%d,%.2f,%.2f,%.2f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%.2e\n",
                    n, A->nnz, memMatKB, memoriaMetodoKB, memTotalKB, densidad, ralidadArchivo, tiempoLectura, tiempoChol, tiempoTotal, iterChol, res);
        }

        if (opSalida == 1 || opSalida == 4) {
            int formato;
            printf("Selecciona formato para Cholesky (1: TXT, 2: CSV): ");
            scanf("%d", &formato);
            if (formato == 1) guardarTxtDetallado("Cholesky", archivoMatriz, n, x_chol);
            else guardarCsvResultados("Cholesky", archivoMatriz, n, x_chol);
        }
        if (opSalida == 2 || opSalida == 4) {
             printf("[Cholesky] Tiempo Total: %.6f s | Ralidad: %.4f | Memoria Total: %.2f KB\n", tiempoTotal, ralidadArchivo, memTotalKB);
        }
        free(x_chol);
    }

    // ================= GRADIENTE CONJUGADO =================
    if (opMetodo == 3 || opMetodo == 4) {
        double tiempoGC = 0.0;
        int iterGC = 0;
        double memoriaMetodoKB = 0.0;
        double* x_gc = (double*)calloc(n, sizeof(double));

        gradienteConjugado(A, b, x_gc, &iterGC, &tiempoGC, &memoriaMetodoKB);

        double res = calcularResidualFinal(A, b, x_gc);
        double tiempoTotal = tiempoLectura + tiempoGC;
        double memTotalKB = memMatKB + memoriaMetodoKB;

        if (csv) {
            fprintf(csv, "GradienteConjugado,%d,%d,%.2f,%.2f,%.2f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%.2e\n",
                    n, A->nnz, memMatKB, memoriaMetodoKB, memTotalKB, densidad, ralidadArchivo, tiempoLectura, tiempoGC, tiempoTotal, iterGC, res);
        }

        if (opSalida == 1 || opSalida == 4) {
            int formato;
            printf("Selecciona formato para Gradiente Conjugado (1: TXT, 2: CSV): ");
            scanf("%d", &formato);
            if (formato == 1) guardarTxtDetallado("GradienteConjugado", archivoMatriz, n, x_gc);
            else guardarCsvResultados("GradienteConjugado", archivoMatriz, n, x_gc);
        }
        if (opSalida == 2 || opSalida == 4) {
             printf("[Gradiente Conjugado] Tiempo Total: %.6f s | Ralidad: %.4f | Memoria Total: %.2f KB\n", tiempoTotal, ralidadArchivo, memTotalKB);
        }
        free(x_gc);
    }

    if (csv) fclose(csv);
    free(b);
    free(x);
    liberarMatrizCOO(A);

    return 0;
}
