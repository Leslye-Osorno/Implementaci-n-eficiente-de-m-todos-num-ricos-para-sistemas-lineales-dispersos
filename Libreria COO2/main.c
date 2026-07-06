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
    char* primerGuion = strchr(archivoMatriz, '_');
    if (primerGuion != NULL) {
        sscanf(primerGuion + 1, "%lf", &ralidadArchivo);
    }

    // 1. Cargar la matriz en formato COO
    A = leerMatrizArchivo(archivoMatriz, &tiempoLectura);
    if (!A) {
        printf("Error: No se pudo cargar la matriz.\n");
        return 1;
    }
    int n = A->n;
    
    double densidad = (double)A->nnz / (double)(A->n * A->m);
    // Cálculo preciso de memoria para la matriz en formato COO (3 arreglos planos + estructura)
    double memMatKB = (sizeof(MatrixCOO) + A->nnz * (sizeof(int) * 2 + sizeof(double))) / 1024.0;

    // 2. GENERACIÓN AUTOMÁTICA DEL VECTOR B (Asegura Solución x = [1,1,...,1])
    // Usamos calloc para inicializar todo el vector b en 0.0 obligatoriamente
    double* b = (double*)calloc(n, sizeof(double));
    if (!b) {
        printf("Error: No hay memoria suficiente para el vector b.\n");
        liberarMatrizCOO(A);
        return 1;
    }
    
    // Al acumular los valores de val en la posición row[k], calculamos la suma por fila de la matriz
    for (int k = 0; k < A->nnz; k++) {
        b[A->row[k]] += A->val[k];
    }

    // 3. Selección de Métodos y Salidas
    int opMetodo, opSalida;
    printf("\nMetodo:\n1. PLU\n2. Cholesky\n3. GC\n4. Todos\nOpcion: ");
    if (scanf("%d", &opMetodo) != 1) opMetodo = 4;
    
    printf("\nSalida:\n1. TXT/Resultados\n2. Solo CSV estadisticos\n3. Solo pantalla\n4. Resultados + Estadisticos\nOpcion: ");
    if (scanf("%d", &opSalida) != 1) opSalida = 4;

    // Preparación o apertura del archivo estadístico CSV
    FILE* csv = fopen("estadisticos_libreriaCOO.csv", "a");
    if (csv) {
        fseek(csv, 0, SEEK_END);
        if (ftell(csv) == 0) {
            fprintf(csv, "Metodo,N,NNZ,MemMat(KB),MemMet(KB),MemTotal(KB),Densidad,Ralidad,TLectura(s),TMetodo(s),TTotal(s),Iter,Residual\n");
        }
    }

    // ================= MÉTODO 1: PLU =================
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
            fprintf(csv, "PLU,%d,%d,%.2f,%.2f,%.2f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%.2e\n",
                    n, A->nnz, memMatKB, memoriaMetodoKB, memTotalKB, densidad, ralidadArchivo, tiempoLectura, tiempoPLU, tiempoTotal, iterPLU, res);
        }
        if (opSalida == 1 || opSalida == 4) guardarTxtDetallado("PLU", archivoMatriz, n, x_plu);
        free(x_plu);
    }

    // ================= MÉTODO 2: CHOLESKY =================
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
        if (opSalida == 1 || opSalida == 4) guardarTxtDetallado("Cholesky", archivoMatriz, n, x_chol);
        free(x_chol);
    }

    // ================= MÉTODO 3: GRADIENTE CONJUGADO =================
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
        if (opSalida == 1 || opSalida == 4) guardarTxtDetallado("GradienteConjugado", archivoMatriz, n, x_gc);
        free(x_gc);
    }

    // Limpieza de memoria
    if (csv) fclose(csv);
    free(b);
    liberarMatrizCOO(A);

    printf("\n[PROCESO COMPLETADO] Resultados del sistema COO calculados y guardados.\n");
    return 0;
}
