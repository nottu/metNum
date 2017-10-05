//
// Created by Javier Peralta on 9/16/17.
//

#include "matrix.h"
//#include <omp.h>

void vectorScalar (double *v, double d, int size){
  for(int i = 0; i < size; ++i ){
    v[i] *= d;
  }
}
void restaVector(double *v1, double *v2, double* out, int size){
//#pragma omp parallel for
  for(int i = 0; i < size; ++i ){
    out[i] = v1[i] - v2[i];
  }
}
void sumaVector(double *v1, double *v2, double* out, int size){
//#pragma omp parallel for
  for(int i = 0; i < size; ++i ){
    out[i] = v1[i] + v2[i];
  }
}
void multVector(double *vec1, double *vec2, double* out, int size){
//#pragma omp parallel for
  for (int i = 0; i < size; i++) {
    out[i] = vec1[i] * vec2[i];
  }
}
double productoPunto(double *vec1, double *vec2, int size){
  double c = 0;
  for (int i = 0; i < size; i++) {
    c += vec1[i] * vec2[i];
  }
  return c;
}
void productoPuntoA(double *vec1, double *vec2, double* vec3, int size){
//#pragma omp parallel for
  for (int i = 0; i < size; i++) {
    vec3[i] = vec1[i] * vec2[i];
  }
}

void multMatriz(double **mat1, double **mat2, int n, int m, int p, int q, double **res){
  //fila * columna
  if (m != p) {
    perror("Numero de filas de la primera matriz debe ser igual numero de columnas de la segunda\n");
    return;
  }
//#pragma omp parallel for
  for (int i = 0; i < n; ++i) {
    double *fila = res[i];
//#pragma omp parallel for
    for (int j = 0; j < q; ++j) {
      double c = 0;
//#pragma omp parallel for reduction(+:c)
      for (int k = 0; k < m; ++k) {
        c += mat1[i][k] * mat2[k][j];
      }
      fila[j] = c;
    }
  }
}
void  multMatrizVect(double **mat, double *vec, int n, int m, double* res){
  for (int i = 0; i < n; i++) {
    res[i] = productoPunto(mat[i], vec, m);
  }
}

//other
void printVect(double * a, int n){
  for (int i = 0; i < n; ++i) {
    printf("%g ", a[i]);
  }
  printf("\n");
}
void printMtx(double**a, int nr, int nc){
  for (int i = 0; i < nr; ++i) {
    printVect(a[i], nc);
  }
}

double *readVector(char* name, int* sz){
  FILE *f = fopen(name, "rb");
  if (!f) return NULL;
  fread(sz, sizeof(int), 1, f);
  double *vect = malloc(sizeof(double) * *sz);
  for (int i = 0; i < *sz; ++i) {
    fread(vect, sizeof(double), *sz, f);
  }
  fclose(f);
  return vect;
}
double **readMtx(char* name, int* nr, int* nc){
  FILE *f = fopen(name, "rb");
  if (!f) return NULL;
  fread(nr, sizeof(int), 1, f);
  fread(nc, sizeof(int), 1, f);
  double **mtx = allocMtx(*nr, *nc);
  for (int i = 0; i < *nr; ++i) {
    fread(mtx[i], sizeof(double), (unsigned int)*nc, f);
  }
  fclose(f);
  return mtx;
}
double **allocMtx(int nr, int nc){
  double **mtx = malloc((sizeof(double*)*nr) + sizeof(int));
  int *indi = (int*)mtx;
  mtx = (void*)indi+ sizeof(int);
  if(nr * nc * sizeof(double) < MTXMAXSIZE) {
    indi[0] = 0; //indicate 1 block
    mtx[0] = malloc(sizeof(double) * nr*nc);
    for (int i = 1; i < nr; ++i) {
      mtx[i] = mtx[i-1] + nc;
    }
  } else {
    indi[0] = nr; //indicate nr block
    for (int i = 0; i < nr; ++i) {
      mtx[i] = malloc(sizeof(double) * nc);
    }
  }
  return mtx;
}
void freeMtx(double**a){
  void *indi = (void*)a - sizeof(int);
  int nr = ((int*)indi)[0];
  if(nr){
    for (int i = 0; i < nr; ++i) free(a[i]);
  }
  else free(a[0]);
  free(indi);
}
//
double norma2Vect(double* v, int size){
  return sqrt(norma2VectSq(v, size));
}
double norma2VectSq(double* v, int size){
  double c = 0;
//#pragma omp parallel for reduction(+:c)
  for (int i = 0; i < size; i++) {
    double val =  v[i];
    c += val * val;
  }
  return c;
}
void normalizaVect(double *v, int size){
  double norm = sqrt(norma2VectSq(v, size));
  for (int i = 0; i < size; i++) v[i] /= norm;
}
double diffVectSq(double* v1, double* v2, int size){
  double c = 0;
//#pragma omp parallel for reduction(+:c)
  for (int i = 0; i < size; i++) {
    double val =  v1[i] - v2[i];
    c += val * val;
  }
  return c;
}

double diffMatrizSq(double** m1, double** m2, int nr, int nc){
//#pragma omp parallel for reduction(+:c)
  double c = 0;
  int sz = nr*nc;
  for (int i = 0; i < sz; ++i) {
    double dif = m1[0][i] - m2[0][i];
    c += dif* dif;
  }
  return c;
}

double* diagSol(double*a , double*b, int n){
  double *vect = malloc(sizeof(double) * n);
//#pragma omp parallel
  for (int i = 0; i < n; ++i) {
    if (a[i] == 0){
      if (b[i] != 0){
        printf("Sin solución, X%d no tiene valor\n", i);
        return NULL;
      }
      printf("Multiples Soluciones, X%d puede tener cualquier valor\n", i);
      vect[i] = 0;
      continue;
    }
    vect[i] = b[i]/a[i];
  }
  return vect;
}
double* upperSol(double**a , double*b, int nr, int nc){
  double *vect = malloc(sizeof(double) * nr);

  for (int i = nr -1; i >= 0; i--) {
    double tmp = b[i];
    for (int j = i+1; j < nc; ++j) {
      tmp -= vect[j] * a[i][j];
    }
    vect[i] = tmp / a[i][i];
  }
  return vect;
}
double* lowerSol(double**a , double*b, int nr, int nc){
  double *vect = malloc(sizeof(double) * nr);

  for (int i = 0; i < nr; ++i) {
    double tmp = b[i];
    for (int j = 0; j < i && j < nc; ++j) {
      tmp -= vect[j] * a[i][j];
    }
    tmp /= a[i][i];
    vect[i] = tmp;
  }

  return vect;
}
int luFactor(double** a, double **l, double **u, int nr, int nc){
  for (int i = 0; i < nr; ++i) {
    u[i][i] = 1;
    for (int j = 0; j <= i && j <nc; ++j) {
      double lij = a[i][j];
      for (int k = 0; k < j; ++k) {
        lij -= l[i][k]*u[k][j];
      }
      l[i][j] = lij;
    }
    for (int j = i+1; j < nc; ++j) {
      double lij = a[i][j];
      if(fabs(l[i][i]) < ZERO)
        return 0;
      for (int k = 0; k < i; ++k) {
        lij -= l[i][k]*u[k][j];
      }
      lij /= l[i][i];
      u[i][j] = lij;
    }
  }
  return 1;
}
void luSolver(double **l, double **u, double *b, int nr, int nc){}
//same as lu factor, but in 1 matrix
int luFactor2(double **a, int nr, int nc){
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j <= i && j <nc; ++j) {
      double lij = a[i][j];
      for (int k = 0; k < j; ++k) {
        lij -= a[i][k]*a[k][j];
      }
      a[i][j] = lij;
    }
    for (int j = i+1; j < nc; ++j) {
      double lij = a[i][j];
      if(fabs(a[i][i]) < ZERO)
        return 0;
      for (int k = 0; k < i; ++k) {
        lij -= a[i][k]*a[k][j];
      }
      lij /= a[i][i];
      a[i][j] = lij;
    }
  }
  return 1;
}
void luSolver2(double **a, double *b, int nr, int nc){}

double* triDiagSol(double **a, double *d, int nr, int nc){
  double *xi = malloc(sizeof(double) * nr);
  double *ax = a[0], *bx = a[1], *cx = a[2];
  cx[0] /= bx[0];
  d[0] /= bx[0];
  for (int i = 1; i < nc; ++i) {
    double ptemp = bx[i] - (ax[i] * cx[i-1]);
    cx[i] /= ptemp;
    d[i] = (d[i] - ax[i] * d[i-1])/ptemp;
  }
  xi[nr-1] = d[nr-1];
  for (int i = nr-2; i >= 0; --i) {
    xi[i] = d[i] - cx[i] * xi[i+1];
  }
  return xi;
}
double potencia(double **mat, double *eigvec, int nr, int nc, int maxIter, double toler){
  double error;
  for (int i = 0; i < nr; ++i) eigvec[i] = 1;
  double   *y = malloc(sizeof(double) * nr);
  double  *vt = malloc(sizeof(double) * nr);
  double eigV = 0;
  int i = 0;
  do {
    multMatrizVect(mat, eigvec, nr, nc, y);
    memcpy(eigvec, y, nr * sizeof(double));
    normalizaVect(eigvec, nr);
    multMatrizVect(mat, eigvec, nr, nc, vt);
    eigV = productoPunto(eigvec, vt, nr);
    memcpy(y, eigvec, nr * sizeof(double));
    vectorScalar(vt, eigV, nr);
    restaVector(y, vt, vt, nr);
    error = sqrt(norma2VectSq(vt, nr));
  }
  while(++i < maxIter && error > toler);
  free(y); free(vt);
//  printf("Matriz tam %d x %d\n", nr, nc);
//  printf("Valor lambda %lf\n", eigV);
//  printf("Iteraciones realizadas %d\n", i);
//  printf("Error %g\n", error);
  return eigV;
}
double potenciaInv(double **mat, double *eigvec, int nr, int nc, int maxIter, double toler){
  //factorizar matriz en LU
  double **l = allocMtx(nr, nc);
  double **u = allocMtx(nr, nc);
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nc; ++j) {
      l[i][j] = 0;
      u[i][j] = 0;
    }
  }
  luFactor(mat, l, u, nr, nc);
  printf("----------L---------\n");printMtx(l, nr, nc);
  printf("----------L---------\n");printMtx(u, nr, nc);
  luFactor2(mat, nr, nc);
  printf("----------L---------\n");printMtx(mat, nr, nc);
}