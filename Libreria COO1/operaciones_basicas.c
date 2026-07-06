#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "operaciones_basicas.h"
#include "errores.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))


// Funciones para la manipulación de la lista enlazada
static void push_back(List *list, Node *item) {
    if (list->head == NULL) {
        list->head = item;
        list->last = item;
        item->next = NULL;
        return;
    }
    
    list->last->next = item;
    list->last = item;
    item->next = NULL;
}

void insert_entry(List *list, int row, int col, double val) {
    Node *temp = (Node*) malloc(sizeof(Node));
    temp->row = row;
    temp->col = col;
    temp->val = val;
    temp->next = NULL;
    
    push_back(list, temp);
}

void free_list(List *list) {
    Node *temp;
    while (list->head != NULL) {
        temp = list->head;
        list->head = list->head->next;
        free(temp);
    }
    free(list);
}

// Estructura para el índice rápido por filas
typedef struct IndexNode {
    int col;
    double val;
    struct IndexNode* next;
} IndexNode;

typedef struct {
    IndexNode* head;
    IndexNode* last;
} RowList;

typedef struct {
    RowList* rows;
    int n;
} MatrixIndex;

static MatrixIndex* construirIndiceFila(MatrixCOO* A) {
    MatrixIndex* idx = (MatrixIndex*)malloc(sizeof(MatrixIndex));
    idx->n = A->n;
    idx->rows = (RowList*)calloc(A->n, sizeof(RowList));

    for (int k = 0; k < A->nnz; k++) {
        int r = A->row[k];
        int c = A->col[k];
        double val = A->val[k];

        IndexNode* node = (IndexNode*)malloc(sizeof(IndexNode));
        node->col = c;
        node->val = val;
        node->next = NULL;

        if (idx->rows[r].head == NULL) {
            idx->rows[r].head = node;
            idx->rows[r].last = node;
        } else {
            idx->rows[r].last->next = node;
            idx->rows[r].last = node;
        }
    }
    return idx;
}

static void liberarIndiceFila(MatrixIndex* idx) {
    for (int i = 0; i < idx->n; i++) {
        IndexNode* curr = idx->rows[i].head;
        while (curr != NULL) {
            IndexNode* temp = curr;
            curr = curr->next;
            free(temp);
        }
    }
    free(idx->rows);
    free(idx);
}

MatrixCOO* asignarMatrizCOO(int n, int m, int nnz) {
    MatrixCOO* A = (MatrixCOO*)malloc(sizeof(MatrixCOO));

    A->n = n;
    A->m = m;
    A->nnz = nnz;

    A->row = (int*)malloc(nnz * sizeof(int));
    A->col = (int*)malloc(nnz * sizeof(int));
    A->val = (double*)malloc(nnz * sizeof(double));
    
    A->last_idx = 0;

    return A;
}

void liberarMatrizCOO(MatrixCOO* A) {
    if (A == NULL) return;
    free(A->row);
    free(A->col);
    free(A->val);
    free(A);
}

void inicializarMatrizDispersa(MatrixCOO* matrix) {
    for (int i = 0; i < matrix->nnz; i++) {
        matrix->row[i] = 0;
        matrix->col[i] = 0;
        matrix->val[i] = 0.0;
    }
}

void imprimirMatrizCOO(MatrixCOO* A, FILE* archivo) {
    fprintf(archivo, "Matriz COO (nnz=%d)\n", A->nnz);
    for (int i = 0; i < A->nnz; i++) {
        fprintf(archivo, "%d %d %.10f\n", A->row[i], A->col[i], A->val[i]);
    }
}

double obtenerValorCOO(MatrixCOO* A, int row, int col) {
    if (A->nnz == 0) return 0.0;

    // Verificar si está en la posición del cursor actual
    int start = A->last_idx;
    if (start < A->nnz && A->row[start] == row && A->col[start] == col) {
        return A->val[start];
    }

    // Si no, buscar desde el inicio de la matriz (o realizar búsqueda secuencial)
    for (int k = 0; k < A->nnz; k++) {
        int idx = (start + k) % A->nnz; // Recorre el arreglo de forma circular
        if (A->row[idx] == row && A->col[idx] == col) {
            A->last_idx = idx; // Actualizamos el cursor al nuevo índice encontrado
            return A->val[idx];
        }
    }

    return 0.0; // Retorna 0.0 si el elemento no se encuentra (es un cero en la matriz)
}
double productoPunto(const double* a, const double* b, int n) {
    double suma = 0;
    for (int i = 0; i < n; i++)
        suma += a[i] * b[i];
    return suma;
}

double productoPuntoDisperso(int n, double* v1, double* v2) {
    double suma = 0;
    for (int i = 0; i < n; i++) {
        if (fabs(v1[i]) > 1e-12 && fabs(v2[i]) > 1e-12) {
            suma += v1[i] * v2[i];
        }
    }
    return suma;
}

double normaVectorDisperso(int n, double* v) {
    double suma = 0;
    for (int i = 0; i < n; i++)
        suma += v[i] * v[i];
    return sqrt(suma);
}

void multiplicarVectorEscalar(int n, double escalar, double* v, double* resultado) {
    for (int i = 0; i < n; i++)
        resultado[i] = escalar * v[i];
}

void multiplicarMatrizVectorDisperso(MatrixCOO* A, double* x, double* resultado) {
    for (int i = 0; i < A->n; i++)
        resultado[i] = 0.0;
    for (int k = 0; k < A->nnz; k++) {
        int r = A->row[k];
        int c = A->col[k];
        resultado[r] += A->val[k] * x[c];
    }
}

static void swapCOO(MatrixCOO *A, int i, int j) {
    int temp_row = A->row[i];
    int temp_col = A->col[i];
    double temp_val = A->val[i];
    
    A->row[i] = A->row[j];
    A->col[i] = A->col[j];
    A->val[i] = A->val[j];
    
    A->row[j] = temp_row;
    A->col[j] = temp_col;
    A->val[j] = temp_val;
}

static void ordenaEntradasMatrixCOO(MatrixCOO *A) {
    int sorted = 0;
    while (!sorted) {
        sorted = 1;
        for (int i = 0; i < A->nnz - 1; i++) {
            if (A->row[i+1] < A->row[i]) {
                swapCOO(A, i, i+1);
                sorted = 0;
            } else if (A->row[i+1] == A->row[i]) {
                if (A->col[i+1] < A->col[i]) {
                    swapCOO(A, i, i+1);
                    sorted = 0;
                }
            } 
        }
    }
}

MatrixCOO* sumarMatricesDispersas(MatrixCOO* A, MatrixCOO* B) {
    if (A->n != B->n || A->m != B->m) {
        printf("Error: Matrices de diferentes dimensiones (A: %dx%d, B: %dx%d)\n", A->n, A->m, B->n, B->m);
        return NULL;
    }

    int n = A->n;
    int m = A->m;
    int counter_A, counter_B;
    List *temp_Matrix = (List*) malloc(sizeof(List));
    temp_Matrix->head = NULL;
    temp_Matrix->last = NULL;
    int k = 0;

    MatrixCOO* tempA = asignarMatrizCOO(A->n, A->m, A->nnz);
    MatrixCOO* tempB = asignarMatrizCOO(B->n, B->m, B->nnz);
    for (int i = 0; i < A->nnz; i++) {
        tempA->row[i] = A->row[i];
        tempA->col[i] = A->col[i];
        tempA->val[i] = A->val[i];
    }
    for (int i = 0; i < B->nnz; i++) {
        tempB->row[i] = B->row[i];
        tempB->col[i] = B->col[i];
        tempB->val[i] = B->val[i];
    }

    ordenaEntradasMatrixCOO(tempA);
    ordenaEntradasMatrixCOO(tempB);

    counter_A = 0;
    counter_B = 0;
    while (counter_A < tempA->nnz && counter_B < tempB->nnz) {
        if (tempA->row[counter_A] < tempB->row[counter_B]) {
            insert_entry(temp_Matrix, tempA->row[counter_A], tempA->col[counter_A], tempA->val[counter_A]);
            counter_A++;
            k++;
            continue;
        }
        if (tempA->row[counter_A] > tempB->row[counter_B]) {
            insert_entry(temp_Matrix, tempB->row[counter_B], tempB->col[counter_B], tempB->val[counter_B]);
            counter_B++;
            k++;
            continue;
        }
        if (tempA->row[counter_A] == tempB->row[counter_B]) {
            if (tempA->col[counter_A] < tempB->col[counter_B]) {
                insert_entry(temp_Matrix, tempA->row[counter_A], tempA->col[counter_A], tempA->val[counter_A]);
                counter_A++;
                k++;
                continue;
            }
            if (tempA->col[counter_A] > tempB->col[counter_B]) {
                insert_entry(temp_Matrix, tempB->row[counter_B], tempB->col[counter_B], tempB->val[counter_B]);
                counter_B++;
                k++;
                continue;
            }
            if (tempA->col[counter_A] == tempB->col[counter_B]) {
                double val = tempA->val[counter_A] + tempB->val[counter_B];
                if (fabs(val) > 1e-12) {
                    insert_entry(temp_Matrix, tempA->row[counter_A], tempA->col[counter_A], val);
                    k++;
                }
                counter_A++;
                counter_B++;
                continue;
            }
        }
    }

    while (counter_A < tempA->nnz) {
        insert_entry(temp_Matrix, tempA->row[counter_A], tempA->col[counter_A], tempA->val[counter_A]);
        counter_A++;
        k++;
    }

    while (counter_B < tempB->nnz) {
        insert_entry(temp_Matrix, tempB->row[counter_B], tempB->col[counter_B], tempB->val[counter_B]);
        counter_B++;
        k++;
    }

    MatrixCOO* C = asignarMatrizCOO(n, m, k);
    Node *it = temp_Matrix->head;
    int i = 0;
    while (it != NULL) {
        C->row[i] = it->row;
        C->col[i] = it->col;
        C->val[i] = it->val;
        i++;
        it = it->next;
    }
    C->nnz = i;

    liberarMatrizCOO(tempA);
    liberarMatrizCOO(tempB);
    free_list(temp_Matrix);
    return C;
}

MatrixCOO* restarMatricesDispersas(MatrixCOO* A, MatrixCOO* B) {
    if (A->n != B->n || A->m != B->m) {
        printf("Error: Matrices de diferentes dimensiones (A: %dx%d, B: %dx%d)\n", A->n, A->m, B->n, B->m);
        return NULL;
    }

    int n = A->n;
    int m = A->m;
    int counter_A, counter_B;
    List *temp_Matrix = (List*) malloc(sizeof(List));
    temp_Matrix->head = NULL;
    temp_Matrix->last = NULL;
    int k = 0;

    MatrixCOO* tempA = asignarMatrizCOO(A->n, A->m, A->nnz);
    MatrixCOO* tempB = asignarMatrizCOO(B->n, B->m, B->nnz);
    for (int i = 0; i < A->nnz; i++) {
        tempA->row[i] = A->row[i];
        tempA->col[i] = A->col[i];
        tempA->val[i] = A->val[i];
    }
    for (int i = 0; i < B->nnz; i++) {
        tempB->row[i] = B->row[i];
        tempB->col[i] = B->col[i];
        tempB->val[i] = B->val[i];
    }

    ordenaEntradasMatrixCOO(tempA);
    ordenaEntradasMatrixCOO(tempB);

    counter_A = 0;
    counter_B = 0;
    while (counter_A < tempA->nnz && counter_B < tempB->nnz) {
        if (tempA->row[counter_A] < tempB->row[counter_B]) {
            insert_entry(temp_Matrix, tempA->row[counter_A], tempA->col[counter_A], tempA->val[counter_A]);
            counter_A++;
            k++;
            continue;
        }
        if (tempA->row[counter_A] > tempB->row[counter_B]) {
            insert_entry(temp_Matrix, tempB->row[counter_B], tempB->col[counter_B], -tempB->val[counter_B]);
            counter_B++;
            k++;
            continue;
        }
        if (tempA->row[counter_A] == tempB->row[counter_B]) {
            if (tempA->col[counter_A] < tempB->col[counter_B]) {
                insert_entry(temp_Matrix, tempA->row[counter_A], tempA->col[counter_A], tempA->val[counter_A]);
                counter_A++;
                k++;
                continue;
            }
            if (tempA->col[counter_A] > tempB->col[counter_B]) {
                insert_entry(temp_Matrix, tempB->row[counter_B], tempB->col[counter_B], -tempB->val[counter_B]);
                counter_B++;
                k++;
                continue;
            }
            if (tempA->col[counter_A] == tempB->col[counter_B]) {
                double val = tempA->val[counter_A] - tempB->val[counter_B];
                if (fabs(val) > 1e-12) {
                    insert_entry(temp_Matrix, tempA->row[counter_A], tempA->col[counter_A], val);
                    k++;
                }
                counter_A++;
                counter_B++;
                continue;
            }
        }
    }

    while (counter_A < tempA->nnz) {
        insert_entry(temp_Matrix, tempA->row[counter_A], tempA->col[counter_A], tempA->val[counter_A]);
        counter_A++;
        k++;
    }

    while (counter_B < tempB->nnz) {
        insert_entry(temp_Matrix, tempB->row[counter_B], tempB->col[counter_B], -tempB->val[counter_B]);
        counter_B++;
        k++;
    }

    MatrixCOO* C = asignarMatrizCOO(n, m, k);
    Node *it = temp_Matrix->head;
    int i = 0;
    while (it != NULL) {
        C->row[i] = it->row;
        C->col[i] = it->col;
        C->val[i] = it->val;
        i++;
        it = it->next;
    }
    C->nnz = i;

    liberarMatrizCOO(tempA);
    liberarMatrizCOO(tempB);
    free_list(temp_Matrix);
    return C;
}

void sustitucionHaciaAdelante(int n, MatrixCOO* L, double* b, double* y) {
    MatrixIndex* idx = construirIndiceFila(L);

    for (int i = 0; i < n; i++) {
        double suma = 0.0;
        double diag = 0.0;

        IndexNode* curr = idx->rows[i].head;
        while (curr != NULL) {
            int col = curr->col;
            double val = curr->val;

            if (col < i) {
                suma += val * y[col];
            } else if (col == i) {
                diag = val;
            }
            curr = curr->next;
        }

        if (fabs(diag) < 1e-15) {
            printf("ERROR: Elemento diagonal cero o muy pequeño en fila %d\n", i);
            liberarIndiceFila(idx);
            exit(1);
        }

        y[i] = (b[i] - suma) / diag;
    }

    liberarIndiceFila(idx);
}

void sustitucionHaciaAtras(int n, MatrixCOO* U, double* y, double* x) {
    MatrixIndex* idx = construirIndiceFila(U);

    for (int i = n - 1; i >= 0; i--) {
        double suma = 0.0;
        double diag = 0.0;

        IndexNode* curr = idx->rows[i].head;
        while (curr != NULL) {
            int col = curr->col;
            double val = curr->val;

            if (col > i) {
                suma += val * x[col];
            } else if (col == i) {
                diag = val;
            }
            curr = curr->next;
        }

        if (fabs(diag) < 1e-15) {
            printf("ERROR: Elemento diagonal cero o muy pequeño en fila %d\n", i);
            liberarIndiceFila(idx);
            exit(1);
        }

        x[i] = (y[i] - suma) / diag;
    }

    liberarIndiceFila(idx);
}
