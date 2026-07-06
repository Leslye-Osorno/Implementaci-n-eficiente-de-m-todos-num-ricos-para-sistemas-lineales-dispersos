#ifndef METODOS_LINEALES_H
#define METODOS_LINEALES_H

#include "operaciones_basicas.h"

// Se añade el parámetro double* memoriaUsada a las firmas de las funciones
void metodoPLU(MatrixCOO* A, double* b, double* x, int* iter, double* tiempo, double* memoriaUsada);
void metodoCholesky(MatrixCOO* A, double* b, double* x, int* iter, double* tiempo, double* memoriaUsada);
void gradienteConjugado(MatrixCOO* A, double* b, double* x, int* iter, double* tiempo, double* memoriaUsada);

#endif
