#ifndef OPERACIONES_BASICAS_H
#define OPERACIONES_BASICAS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Definición de la estructura de la matriz COO
typedef struct {
    int n, m, nnz;
    int* row;
    int* col;
    double* val;
    int last_idx;
} MatrixCOO;

// Definición de estructuras para las listas enlazadas (utilizadas en métodos directos)
typedef struct TNode {
    int row, col;
    double val;
    struct TNode *next;
} Node;

typedef struct TList {
    Node *head;
    Node *last;
} List;

MatrixCOO* asignarMatrizCOO(int n, int m, int nnz);
void liberarMatrizCOO(MatrixCOO* A);
void inicializarMatrizDispersa(MatrixCOO* matrix);
void imprimirMatrizCOO(MatrixCOO* A, FILE* archivo);
double obtenerValorCOO(MatrixCOO* A, int row, int col);

double productoPunto(const double* a, const double* b, int n);
double productoPuntoDisperso(int n, double* v1, double* v2);
double normaVectorDisperso(int n, double* v);
void multiplicarVectorEscalar(int n, double escalar, double* v, double* resultado);
void multiplicarMatrizVectorDisperso(MatrixCOO* A, double* x, double* resultado);

MatrixCOO* sumarMatricesDispersas(MatrixCOO* A, MatrixCOO* B);
MatrixCOO* restarMatricesDispersas(MatrixCOO* A, MatrixCOO* B);

void sustitucionHaciaAdelante(int n, MatrixCOO* L, double* b, double* y);
void sustitucionHaciaAtras(int n, MatrixCOO* U, double* y, double* x);

// Prototipos de la lista (utilizados por metodos_lineales.c)
void insert_entry(List *list, int row, int col, double val);
void free_list(List *list);

#endif
