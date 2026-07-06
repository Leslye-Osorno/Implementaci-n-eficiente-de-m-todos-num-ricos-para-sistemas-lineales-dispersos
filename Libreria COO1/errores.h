#ifndef ERRORES_H
#define ERRORES_H

#include "operaciones_basicas.h"

int errorMatrizNoCuadrada(MatrixCOO* A);

int errorDiagonalCeros(MatrixCOO* A);

int verificarMatrizSimetrica(MatrixCOO* A);

int verificarDefinidaPositiva(MatrixCOO* A);

int errorDimensionVector(int n,int m);

#endif

