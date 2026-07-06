#include "operaciones_basicas.h"
#include <string.h>


// MÓDULO LIL (Optimizado O(1) para inserción masiva)


MatrixLIL* asignarMatrizLIL(int n, int m) {
    MatrixLIL* A = (MatrixLIL*)malloc(sizeof(MatrixLIL));
    A->n = n; A->m = m; A->nnz = 0;
    A->filas = (FilaLIL*)calloc(n, sizeof(FilaLIL));
    return A;
}

void liberarMatrizLIL(MatrixLIL* A) {
    if (!A) return;
    for (int i = 0; i < A->n; i++) if (A->filas[i].datos) free(A->filas[i].datos);
    free(A->filas);
    free(A);
}

void insertarValorLIL(MatrixLIL* A, int f, int c, double valor) {
    if (fabs(valor) < 1e-18) return; 

    FilaLIL* fila = &A->filas[f];

    // Crecimiento dinámico
    if (fila->tamano >= fila->capacidad) {
        fila->capacidad = (fila->capacidad == 0) ? 4 : fila->capacidad * 2;
        fila->datos = (ElementoLIL*)realloc(fila->datos, fila->capacidad * sizeof(ElementoLIL));
    }
    
    // Inserción ciega O(1) al final
    fila->datos[fila->tamano].col = c;
    fila->datos[fila->tamano].val = valor;
    fila->tamano++;
}

// Función auxiliar para qsort
int cmpElementoLIL(const void* a, const void* b) {
    return ((ElementoLIL*)a)->col - ((ElementoLIL*)b)->col;
}

void ordenarLIL(MatrixLIL* A) {
    int nnz_reales = 0;
    for (int i = 0; i < A->n; i++) {
        if (A->filas[i].tamano > 0) {
            if (A->filas[i].tamano > 1) {
                // Ordenar ultra-rápido por columnas
                qsort(A->filas[i].datos, A->filas[i].tamano, sizeof(ElementoLIL), cmpElementoLIL);
                
                // Compactar duplicados (asegura integridad de los datos)
                int unico = 0;
                for (int j = 1; j < A->filas[i].tamano; j++) {
                    if (A->filas[i].datos[j].col == A->filas[i].datos[unico].col) {
                        A->filas[i].datos[unico].val = A->filas[i].datos[j].val; // Sobrescribir
                    } else {
                        unico++;
                        A->filas[i].datos[unico] = A->filas[i].datos[j];
                    }
                }
                A->filas[i].tamano = unico + 1;
            }
            nnz_reales += A->filas[i].tamano;
        }
    }
    A->nnz = nnz_reales; // Registrar el conteo oficial
}

void ajustarMemoriaLIL(MatrixLIL* A) {
    for (int i = 0; i < A->n; i++) {
        if (A->filas[i].tamano > 0 && A->filas[i].tamano < A->filas[i].capacidad) {
            A->filas[i].capacidad = A->filas[i].tamano;
            A->filas[i].datos = (ElementoLIL*)realloc(A->filas[i].datos, A->filas[i].capacidad * sizeof(ElementoLIL));
        }
    }
}

double obtenerValorLIL(MatrixLIL* A, int f, int c) {
    FilaLIL* fila = &A->filas[f];
    int bajo = 0, alto = fila->tamano - 1;
    while (bajo <= alto) {
        int medio = bajo + (alto - bajo) / 2;
        if (fila->datos[medio].col == c) return fila->datos[medio].val;
        if (fila->datos[medio].col < c) bajo = medio + 1;
        else alto = medio - 1;
    }
    return 0.0;
}

void multiplicarMatrizVectorLIL(MatrixLIL* A, double* x, double* res) {
    for (int i = 0; i < A->n; i++) {
        double suma = 0.0;
        for (int k = 0; k < A->filas[i].tamano; k++) suma += A->filas[i].datos[k].val * x[A->filas[i].datos[k].col];
        res[i] = suma;
    }
}

void sustitucionHaciaAdelante(int n, MatrixLIL* L, double* b, double* y) {
    for (int i = 0; i < n; i++) {
        double suma = 0.0, diag = 1.0;
        for (int k = 0; k < L->filas[i].tamano; k++) {
            if (L->filas[i].datos[k].col < i) suma += L->filas[i].datos[k].val * y[L->filas[i].datos[k].col];
            else if (L->filas[i].datos[k].col == i) diag = L->filas[i].datos[k].val;
        }
        y[i] = (b[i] - suma) / diag;
    }
}

void sustitucionHaciaAtras(int n, MatrixLIL* U, double* y, double* x) {
    for (int i = n - 1; i >= 0; i--) {
        double suma = 0.0, diag = 1.0;
        for (int k = 0; k < U->filas[i].tamano; k++) {
            if (U->filas[i].datos[k].col > i) suma += U->filas[i].datos[k].val * x[U->filas[i].datos[k].col];
            else if (U->filas[i].datos[k].col == i) diag = U->filas[i].datos[k].val;
        }
        x[i] = (y[i] - suma) / diag;
    }
}


// MÓDULO CSR DE ALTO RENDIMIENTO (HPC)


MatrixCSR* convertirLILaCSR(MatrixLIL* A) {
    if (!A) return NULL;
    MatrixCSR* csr = (MatrixCSR*)malloc(sizeof(MatrixCSR));
    csr->n = A->n; csr->m = A->m; csr->nnz = A->nnz;

    csr->valores = (double*)malloc(A->nnz * sizeof(double));
    csr->columnas = (int*)malloc(A->nnz * sizeof(int));
    csr->apuntador_filas = (int*)malloc((A->n + 1) * sizeof(int));
    csr->inv_diag = (double*)malloc(A->n * sizeof(double));

    int idx = 0;
    for (int i = 0; i < A->n; i++) {
        csr->apuntador_filas[i] = idx;
        csr->inv_diag[i] = 1.0; 
        
        for (int k = 0; k < A->filas[i].tamano; k++) {
            int col = A->filas[i].datos[k].col;
            double val = A->filas[i].datos[k].val;
            
            csr->valores[idx] = val;
            csr->columnas[idx] = col;
            
            if (col == i && fabs(val) > 1e-15) csr->inv_diag[i] = 1.0 / val; 
            idx++;
        }
    }
    csr->apuntador_filas[A->n] = A->nnz;
    return csr;
}

void liberarMatrizCSR(MatrixCSR* csr) {
    if (!csr) return;
    free(csr->valores); free(csr->columnas);
    free(csr->apuntador_filas); free(csr->inv_diag); free(csr);
}

void multiplicarMatrizVectorCSR(MatrixCSR* A, double* x, double* res) {
    for (int i = 0; i < A->n; i++) {
        double suma = 0.0;
        int j = A->apuntador_filas[i];
        int fin = A->apuntador_filas[i+1];
        
        for (; j <= fin - 4; j += 4) {
            suma += A->valores[j] * x[A->columnas[j]] +
                    A->valores[j+1] * x[A->columnas[j+1]] +
                    A->valores[j+2] * x[A->columnas[j+2]] +
                    A->valores[j+3] * x[A->columnas[j+3]];
        }
        for (; j < fin; j++) suma += A->valores[j] * x[A->columnas[j]];
        res[i] = suma;
    }
}

void sustitucionHaciaAdelanteCSR(MatrixCSR* L, double* b, double* y) {
    for (int i = 0; i < L->n; i++) {
        double suma = 0.0;
        for (int j = L->apuntador_filas[i]; j < L->apuntador_filas[i+1]; j++) {
            suma += L->valores[j] * y[L->columnas[j]];
        }
        y[i] = (b[i] - suma) * L->inv_diag[i]; 
    }
}

void sustitucionHaciaAtrasCSR(MatrixCSR* U, double* y, double* x) {
    for (int i = U->n - 1; i >= 0; i--) {
        double suma = 0.0;
        for (int j = U->apuntador_filas[i]; j < U->apuntador_filas[i+1]; j++) {
            suma += U->valores[j] * x[U->columnas[j]];
        }
        x[i] = (y[i] - suma) * U->inv_diag[i];
    }
}


// MÓDULO ÁLGEBRA VECTORIAL (Acelerado con restrict)


void sumarVectores(int n, double* restrict v1, double* restrict v2, double* restrict res) {
    for (int i = 0; i < n; i++) res[i] = v1[i] + v2[i];
}

void restarVectores(int n, double* restrict v1, double* restrict v2, double* restrict res) {
    for (int i = 0; i < n; i++) res[i] = v1[i] - v2[i];
}

double productoPunto(int n, double* restrict v1, double* restrict v2) {
    double r = 0.0;
    for (int i = 0; i < n; i++) r += v1[i] * v2[i];
    return r;
}

double norma2(int n, double* restrict v) { 
    return sqrt(productoPunto(n, v, v)); 
}

void saxpy(int n, double a, double* restrict x, double* restrict y, double* restrict res) {
    for (int i = 0; i < n; i++) res[i] = a * x[i] + y[i];
}
