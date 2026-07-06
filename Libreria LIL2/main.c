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

// Función auxiliar para ejecutar y reportar resultados de cada método
void ejecutarMetodo(int id, MatrixLIL* A, double* b, const char* nombre, FILE* csv, int opSalida, char* archNom, double tL, double mMat) {
    double tM = 0.0;
    int iter = 0;
    int n = A->n;
    double* x_sol = (double*)calloc(n, sizeof(double));
    
    if (!x_sol) {
        printf("Error: No hay memoria suficiente para el vector solucion.\n");
        return;
    }

    // Selección del método según ID
    if (id == 1) {
        metodoPLU(A, b, x_sol, &iter, &tM);
    } 
    else if (id == 2) {
        if (errorMatrizNoCuadrada(A) || verificarMatrizSimetrica(A) != 0) {
            printf("[SKIP] %s requiere matriz simetrica y cuadrada.\n", nombre);
            free(x_sol);
            return;
        }
        metodoCholesky(A, b, x_sol, &iter, &tM);
    } 
    else if (id == 3) {
        metodoGradienteConjugado(A, b, x_sol, &iter, &tM);
    }
    
    // g_memoria_metodo_kb se actualiza dentro de cada metodo (HPC dinámico)
    double mMet = g_memoria_metodo_kb; 
    
    // El residual usa la matriz original LIL (Cálculo post-factorización)
    double res = calcularResidualFinal(A, b, x_sol);
    double mFinal = mMat + mMet;

    // Cálculo de densidad y ralidad para la matriz
    double densidad = (double)A->nnz / ((double)n * A->m);
    double ralidad = 1.0 - densidad;

    // Guardado en CSV
    if (csv) {
        fprintf(csv, "%s,%d,%d,%.2f,%.2f,%.2f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%.2e\n",
                nombre, n, A->nnz, mMat, mMet, mFinal, densidad, ralidad, tL, tM, (tL + tM), iter, res);
        fflush(csv); // Asegura que se escriba en disco
    }

    // Salida a pantalla
    if (opSalida == 3 || opSalida == 4) {
        printf("[%s] Tiempo: %.6fs | Mem. Total: %.2fKB | Residual: %.2e\n", nombre, tM, mFinal, res);
    }

    // Guardado en archivo detallado
    if (opSalida == 1 || opSalida == 4) {
        guardarTxtDetallado(nombre, archNom, n, x_sol);
    }

    free(x_sol);
}

int main() {
    MatrixLIL* A = NULL;
    double tiempoLectura = 0.0;
    char archivoMatriz[256];

    printf("Introduce el nombre del archivo de la matriz (.txt): ");
    if (scanf("%255s", archivoMatriz) != 1) return 1;

    // 1. Carga de matriz (Se mantiene en LIL para respetar el formato original)
    A = leerMatrizArchivo(archivoMatriz, &tiempoLectura);
    if (!A) {
        printf("Error: No se pudo cargar la matriz.\n");
        return 1;
    }
    
    int n = A->n;
    double memoriaMatrizKB = obtenerMemoriaMatrizKB(A);

    double* b = (double*)malloc(n * sizeof(double));
    if (!b) {
        printf("Error de memoria para el vector b.\n");
        liberarMatrizLIL(A);
        return 1;
    }
    
    // 2. Configuración del vector b
    int opVector;
    printf("\nSelecciona el vector b:\n1. Vector de 1s\n2. Vector personalizado\nOpcion: ");
    if (scanf("%d", &opVector) != 1) opVector = 1;

    if (opVector == 1) {
        for (int i = 0; i < n; i++) b[i] = 1.0;
    } else {
        int opVecTipo;
        printf("\nTipo de entrada:\n1. Archivo .txt\n2. Consola (separados por comas)\nOpcion: ");
        if (scanf("%d", &opVecTipo) != 1) opVecTipo = 1;

        if (opVecTipo == 1) {
            char archivoVector[256];
            printf("Nombre del archivo del vector: ");
            scanf("%255s", archivoVector);
            FILE* fv = fopen(archivoVector, "r");
            if (fv) {
                for (int i = 0; i < n; i++) {
                    if (fscanf(fv, "%lf", &b[i]) != 1) b[i] = 1.0;
                }
                fclose(fv);
            } else {
                printf("Error abriendo archivo, usando vector de 1s por defecto.\n");
                for (int i = 0; i < n; i++) b[i] = 1.0;
            }
        } else {
            printf("Introduce los %d elementos (ej. 1.2,3.4...): ", n);
            // Lectura segura por consola
            for (int i = 0; i < n; i++) {
                if (i < n - 1) scanf("%lf,", &b[i]);
                else scanf("%lf", &b[i]);
            }
        }
    }

    // 3. Selección de Métodos y Salidas
    int opMetodo, opSalida;
    printf("\nMetodo:\n1. PLU\n2. Cholesky\n3. GC\n4. Todos\nOpcion: ");
    if (scanf("%d", &opMetodo) != 1) opMetodo = 4;
    
    printf("\nSalida:\n1. TXT/Resultados\n2. Solo CSV estadisticos\n3. Solo pantalla\n4. Resultados + Estadisticos\nOpcion: ");
    if (scanf("%d", &opSalida) != 1) opSalida = 4;

    // Preparación del CSV 
    FILE* csv = fopen("estadisticos_libreriaLIL.csv", "a");
    if (csv) {
        fseek(csv, 0, SEEK_END);
        if (ftell(csv) == 0) {
            fprintf(csv, "Metodo,N,NNZ,MemMat(KB),MemMet(KB),MemTotal(KB),Densidad,Ralidad,TLectura(s),TMetodo(s),TTotal(s),Iter,Residual\n");
        }
    }

    // 4. Ejecución
    if (opMetodo == 1 || opMetodo == 4) 
        ejecutarMetodo(1, A, b, "PLU", csv, opSalida, archivoMatriz, tiempoLectura, memoriaMatrizKB);
    
    if (opMetodo == 2 || opMetodo == 4) 
        ejecutarMetodo(2, A, b, "Cholesky", csv, opSalida, archivoMatriz, tiempoLectura, memoriaMatrizKB);
    
    if (opMetodo == 3 || opMetodo == 4) 
        ejecutarMetodo(3, A, b, "GradienteConjugado", csv, opSalida, archivoMatriz, tiempoLectura, memoriaMatrizKB);

    // Limpieza final
    if (csv) fclose(csv);
    free(b);
    liberarMatrizLIL(A);
    printf("\n[OK] El proceso ha concluido correctamente.\n");
    
    return 0;
}
