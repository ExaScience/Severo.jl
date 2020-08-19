#ifndef LAPACK_H
#define LAPACK_H

#ifdef HAVE_F77_UNDERSCORE
# define F77_CALL(x)    x ## _
#else
# define F77_CALL(x)    x
#endif
#define F77_NAME(x)    F77_CALL(x)

#ifndef BLAS_extern
#define BLAS_extern extern
#endif

#ifndef La_extern
#define La_extern extern
#endif

BLAS_extern double /* DNRM2 - 2-norm of a vector */
F77_NAME(dnrm2)(const int *n, const double *dx, const int *incx);
BLAS_extern void   /* DAXPY - replace y by da*x + y */
F77_NAME(daxpy)(const int *n, const double *da,
        const double *dx, const int *incx,
        double *dy, const int *incy);
BLAS_extern void   /* DCOPY - copy x to y */
F77_NAME(dcopy)(const int *n, const double *dx, const int *incx,
        double *dy, const int *incy);
BLAS_extern double /* DDOT - inner product of x and y */
F77_NAME(ddot)(const int *n, const double *dx, const int *incx,
           const double *dy, const int *incy);
BLAS_extern void   /* DSCAL - scale a one-dimensional array */
F77_NAME(dscal)(const int *n, const double *alpha, double *dx, const int *incx);

/* DGEMV - perform one of the matrix-vector operations */
/* y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y,  */
BLAS_extern void
F77_NAME(dgemv)(const char *trans, const int *m, const int *n,
        const double *alpha, const double *a, const int *lda,
        const double *x, const int *incx, const double *beta,
        double *y, const int *incy);

/* DGEMM - perform one of the matrix-matrix operations    */
/* C := alpha*op( A )*op( B ) + beta*C */
BLAS_extern void
F77_NAME(dgemm)(const char *transa, const char *transb, const int *m,
        const int *n, const int *k, const double *alpha,
        const double *a, const int *lda,
        const double *b, const int *ldb,
        const double *beta, double *c, const int *ldc);

La_extern void
F77_NAME(dbdsdc)(const char* uplo, const char* compq, int *n,
    double * d, double *e, double *u, int *ldu, double *vt,
    int *ldvt, double *q, int *iq, double *work, int * iwork, int *info);

/* DLACPY - copy all or part of a two-dimensional matrix A to */
/* another matrix B */
La_extern void
F77_NAME(dlacpy)(const char* uplo, const int* m, const int* n,
         const double* a, const int* lda,
         double* b, const int* ldb);

/* DGESDD - compute the singular value decomposition (SVD); of a   */
/* real M-by-N matrix A, optionally computing the left and/or      */
/* right singular vectors.  If singular vectors are desired, it uses a */
/* divide-and-conquer algorithm.                   */
La_extern void
F77_NAME(dgesdd)(const char* jobz,
         const int *m, const int *n,
         double *a, const int *lda, double *s,
         double *u, const int *ldu,
         double *vt, const int *ldvt,
         double *work, const int *lwork, int *iwork, int *info);

/* DGEQRF - compute a QR factorization of a real M-by-N matrix A */
La_extern void
F77_NAME(dgeqrf)(const int* m, const int* n, double* a, const int* lda,
         double* tau, double* work, const int* lwork, int* info);
/* DORGQR - generate an M-by-N real matrix Q with orthonormal */
/* columns, */
La_extern void
F77_NAME(dorgqr)(const int* m, const int* n, const int* k,
         double* a, const int* lda, const double* tau,
         double* work, const int* lwork, int* info);
/* DORMQR - overwrite the general real M-by-N matrix C with   SIDE = */
/* 'L' SIDE = 'R' TRANS = 'N' */
La_extern void
F77_NAME(dormqr)(const char* side, const char* trans,
         const int* m, const int* n, const int* k,
         const double* a, const int* lda,
         const double* tau, double* c, const int* ldc,
         double* work, const int* lwork, int* info);

/* DGEQR - Computes a QR factorization of a real M-by-N matrix, with best performance for tall and skinny matrices */
La_extern void
F77_NAME(dgeqrf)(const int* m, const int* n, double* a, const int* lda,
         double*t, const int *tsize, double* work, const int* lwork, int* info);

La_extern void
F77_NAME(dgemqr)(const char *side, const char *trans,
         const int *m, const int *n, const int *k,
         const double *a, const int *lda,
         double *t, int *tsize,
         double *c, const int *ldc,
         double *work, int lwork, int *info);

#endif
