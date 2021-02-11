#ifndef LAPACK_H
#define LAPACK_H

#include <stdint.h>

#ifdef HAVE_OPENBLAS64
#  define F77_CALL(x) x ## _64_
#  define HAVE_BLAS64
#else
#  ifdef HAVE_F77_UNDERSCORE
#   define F77_CALL(x)    x ## _
#  else
#   define F77_CALL(x)    x
#  endif
#endif

#ifdef HAVE_BLAS64
typedef int64_t blas_int;
#else
typedef int32_t blas_int;
#endif //BLAS64

#define F77_NAME(x)    F77_CALL(x)

#ifndef BLAS_extern
#define BLAS_extern extern
#endif

#ifndef La_extern
#define La_extern extern
#endif

BLAS_extern double /* DNRM2 - 2-norm of a vector */
F77_NAME(dnrm2)(const blas_int *n, const double *dx, const blas_int *incx);
BLAS_extern void   /* DAXPY - replace y by da*x + y */
F77_NAME(daxpy)(const blas_int *n, const double *da,
        const double *dx, const blas_int *incx,
        double *dy, const blas_int *incy);
BLAS_extern void   /* DCOPY - copy x to y */
F77_NAME(dcopy)(const blas_int *n, const double *dx, const blas_int *incx,
        double *dy, const blas_int *incy);
BLAS_extern double /* DDOT - inner product of x and y */
F77_NAME(ddot)(const blas_int *n, const double *dx, const blas_int *incx,
           const double *dy, const blas_int *incy);
BLAS_extern void   /* DSCAL - scale a one-dimensional array */
F77_NAME(dscal)(const blas_int *n, const double *alpha, double *dx, const blas_int *incx);

/* DGEMV - perform one of the matrix-vector operations */
/* y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y,  */
BLAS_extern void
F77_NAME(dgemv)(const char *trans, const blas_int *m, const blas_int *n,
        const double *alpha, const double *a, const blas_int *lda,
        const double *x, const blas_int *incx, const double *beta,
        double *y, const blas_int *incy);

/* DGEMM - perform one of the matrix-matrix operations    */
/* C := alpha*op( A )*op( B ) + beta*C */
BLAS_extern void
F77_NAME(dgemm)(const char *transa, const char *transb, const blas_int *m,
        const blas_int *n, const blas_int *k, const double *alpha,
        const double *a, const blas_int *lda,
        const double *b, const blas_int *ldb,
        const double *beta, double *c, const blas_int *ldc);

La_extern void
F77_NAME(dbdsdc)(const char* uplo, const char* compq, blas_int *n,
    double * d, double *e, double *u, blas_int *ldu, double *vt,
    blas_int *ldvt, double *q, blas_int *iq, double *work, blas_int * iwork, blas_int *info);

/* DLACPY - copy all or part of a two-dimensional matrix A to */
/* another matrix B */
La_extern void
F77_NAME(dlacpy)(const char* uplo, const blas_int* m, const blas_int* n,
         const double* a, const blas_int* lda,
         double* b, const blas_int* ldb);

/* DGESDD - compute the singular value decomposition (SVD); of a   */
/* real M-by-N matrix A, optionally computing the left and/or      */
/* right singular vectors.  If singular vectors are desired, it uses a */
/* divide-and-conquer algorithm.                   */
La_extern void
F77_NAME(dgesdd)(const char* jobz,
         const blas_int *m, const blas_int *n,
         double *a, const blas_int *lda, double *s,
         double *u, const blas_int *ldu,
         double *vt, const blas_int *ldvt,
         double *work, const blas_int *lwork, blas_int *iwork, blas_int *info);

/* DGEQRF - compute a QR factorization of a real M-by-N matrix A */
La_extern void
F77_NAME(dgeqrf)(const blas_int* m, const blas_int* n, double* a, const blas_int* lda,
         double* tau, double* work, const blas_int* lwork, blas_int* info);
/* DORGQR - generate an M-by-N real matrix Q with orthonormal */
/* columns, */
La_extern void
F77_NAME(dorgqr)(const blas_int* m, const blas_int* n, const blas_int* k,
         double* a, const blas_int* lda, const double* tau,
         double* work, const blas_int* lwork, blas_int* info);
/* DORMQR - overwrite the general real M-by-N matrix C with   SIDE = */
/* 'L' SIDE = 'R' TRANS = 'N' */
La_extern void
F77_NAME(dormqr)(const char* side, const char* trans,
         const blas_int* m, const blas_int* n, const blas_int* k,
         const double* a, const blas_int* lda,
         const double* tau, double* c, const blas_int* ldc,
         double* work, const blas_int* lwork, blas_int* info);

/* DGEQR - Computes a QR factorization of a real M-by-N matrix, with best performance for tall and skinny matrices */
La_extern void
F77_NAME(dgeqr)(const blas_int* m, const blas_int* n, double* a, const blas_int* lda,
         double*t, const blas_int *tsize, double* work, const blas_int* lwork, blas_int* info);

La_extern void
F77_NAME(dgemqr)(const char *side, const char *trans,
         const blas_int *m, const blas_int *n, const blas_int *k,
         const double *a, const blas_int *lda,
         double *t, blas_int *tsize,
         double *c, const blas_int *ldc,
         double *work, blas_int lwork, blas_int *info);

#endif
