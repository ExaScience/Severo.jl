#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "lapack.h"

// y = A*x if trans = 'N' or y = A'x if trans = 'T', temp is an n vector
typedef void (*matmul_t)(double *y, char trans, double * x, double *temp, void *data);
typedef void (*randn_t)(double *y, blas_int n, void *data);

/* orthog(X,Y,...)
 * compute Y = Y - X * t(X) * Y
 * xm,xn: nrow, ncol X
 * yn: ASSUMED TO BE 1
 * On entry, number of rows of Y must be xm to compute t(X) * Y and
 * T must be allocated of at least size xn.
 * Modifies contents of Y.
 */
void orthog(double *X, double *Y, double *T, blas_int xm, blas_int xn) {
	double a, b;
	blas_int inc = 1;

	// T = t(X) * Y
	a = 1.0; b = 0.0;
	F77_NAME(dgemv)("t", &xm, &xn, &a, X, &xm, Y, &inc, &b, T, &inc);

	// Y = Y - X * T
	a = -1.0; b = 1.0;
	F77_NAME(dgemv)("n", &xm, &xn, &a, X, &xm, T, &inc, &b, Y, &inc);
}

void ablanzbd(blas_int j, blas_int m, blas_int n, blas_int m_b, double SVTol,
		double *B, // output (m_b x m_b)
		double *V, // output (n x m_b)
		double *W, // output (m x m_b)
		double *F, // output (n)
		double *T, // work (n)
		randn_t randn, matmul_t matmul, void *data) {

	double S, R, SS, RR;
	blas_int inc = 1;

	matmul(W + j * m, 'N', V + j * n, T, data);

	if(j != 0)
		orthog(W, W + j * m, T, m, j);

	S = F77_NAME(dnrm2)(&m, W + j * m, &inc);
	if(S < SVTol) {
		randn(W + (j+1) * m, m, data);
		orthog(W, W + j * m, T, m, j);
		S = F77_NAME(dnrm2)(&m, W + (j + 1) * m, &inc);
	}

	SS = 1.0 / S;
	F77_NAME(dscal)(&m, &SS, W + j * m, &inc);

	while(j < m_b) {
		matmul(F, 'T', W + j * m, T, data);
		SS = -S;
		F77_NAME(daxpy)(&n, &SS, V + j * n, &inc, F, &inc);
		orthog(V, F, T, n, j + 1);

		if(j + 1 < m_b) {
			R = F77_NAME(dnrm2)(&n, F, &inc);
			if(R < SVTol) {
				randn(F, n, data);
				orthog(V, F, T, n, j + 1);
				R = F77_NAME(dnrm2)(&n, F, &inc);
			}

			RR = 1.0 / R;
			F77_NAME(dcopy)(&n, F, &inc, V + (j + 1) * n, &inc);
			F77_NAME(dscal)(&n, &RR, V + (j + 1) * n, &inc);
			B[j * m_b + j] = S;
			B[(j+1) * m_b + j] = R;

			matmul(W + (j + 1) * m, 'N', V + (j + 1) * n, T, data);

			/* One step of block classical Gram-Schmidt process */
			RR = -R;
			F77_NAME(daxpy)(&m, &RR, W + j * m, &inc, W + (j + 1) * m, &inc);

			/* full re-orthogonalization step. "long vectors" */
			orthog(W, W + (j + 1) * m, T, m, j + 1);
			S = F77_NAME(dnrm2)(&m, W + (j + 1) * m, &inc);

			if(S < SVTol) {
				randn(W + (j+1) * m, m, data);
				orthog(W, W + (j + 1) * m, T, m, j + 1);
				S = F77_NAME(dnrm2)(&m, W + (j + 1) * m, &inc);
			}

			SS = 1.0 / S;
			F77_NAME(dscal)(&m, &SS, W + (j + 1) * m, &inc);
		} else {
			B[j * m_b + j] = S;
		}

		j++;
	}
}

int convtest(blas_int m_b, double *res, double *svratio, blas_int n, double tol, double SVtol, blas_int *k) {
	blas_int j;
	int Len_res = 0, Len_res_sv = 0;
	int augment = 3;

	for(j = 0; j < n; j++) {
		if(fabs(res[j]) < tol) {
			if(svratio[j] < SVtol)
				Len_res_sv++;
			Len_res++;
		}
	}

	if(Len_res_sv >= n)
		return 0;

	if(Len_res < augment)
		augment = Len_res;

	if(*k < n + augment)
		*k = n + augment;
	if(*k > m_b - 1)
		*k = m_b - 1;
	return -2;
}

static double eps() {
	return nextafter(0.0, 1.0);
}

void print_matrix(double *A, blas_int m, blas_int n) {
	blas_int i, j;
	for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++)
			printf("%f ", A[j*m + i]);
		puts("");
	}
}

blas_int prep_svd(blas_int m, blas_int n, double *A, double *S, double *U, double *VT, blas_int *iwork) {
	blas_int info;
	blas_int lwork = -1;
	double work;
	blas_int mn = m < n ? m : n;

	F77_NAME(dgesdd)("O", &m, &n, A, &m, S, U, &m, VT, &n, &work, &lwork, iwork, &info);

	if(info == 0) {
		return (int)nextafter(work, work + 1.0);
	} else {
		return 7 * mn * (1 + mn);
	}
}

int do_svd(blas_int m, blas_int n, double *A, double *S, double *U, double *VT,
						 double *work, blas_int lwork, blas_int *iwork) {
	blas_int info;
	F77_NAME(dlacpy)("N", &m, &n, A, &m, U, &m);
	F77_NAME(dgesdd)("O", &m, &n, U, &m, S, U, &m, VT, &n, work, &lwork, iwork, &info);
	return info;
}

blas_int irlba_(blas_int m, blas_int n, blas_int nu, blas_int m_b, blas_int maxit, blas_int restart,
		double tol,
		double *so, // output singular values
		double *Uo, // output left singular vectors
		double *Vo, // output right singular vectors
		double *V, // work (n x m_b)
		double *W, // work (m x m_b)
		double *F, // work (n)
		double *T, // work (n)
		double *B, // work (m_b x m_b)
		double *BS, // work (m_b),
		double *BU, // work (m_b x m_b)
		double *BV, // work (m_b x m_b)
		double *U1, // work (m x m_b)
		double *V1, // work (n x m_b)
		double *work, // work (lwork)
		blas_int lwork,
		blas_int *iwork, // work (8 * m_b)
		randn_t randn, matmul_t matmul, void *data) {
	int iter, i;
	blas_int info = 0;
	blas_int inc = 1;
	blas_int k = restart;
	double R, RR;
	double alpha = 1.0, beta = 0.0;
	double Smax = -INFINITY, Smin = INFINITY;

	double SVtol = sqrt(eps());
	if(tol < SVtol) SVtol = tol;

	memset(so, 0, nu);

	iter = 0;
	while(iter < maxit) {
		/* Compute Lanczos bidiagonalization decomposition B */
		ablanzbd(k, m, n, m_b, SVtol, B, V, W, F, T, randn, matmul, data);

		/* Compute the norm of the vector F, and normalize F. */
		R = F77_NAME(dnrm2)(&n, F, &inc);
		RR = 1.0 / R;
		F77_NAME(dscal)(&n, &RR, F, &inc);

		/* Compute singular triplets of B */
		info = do_svd(m_b, m_b, B, BS, BU, BV, work, lwork, iwork);
		if( info != 0 )
			break;

		if( BS[0] > Smax ) Smax = BS[0];
		if( BS[m_b-1] < Smin ) Smin = BS[m_b-1];

		for(i = 0; i < m_b; ++i)
			T[i] = R * BU[i * m_b + (m_b - 1)];

		for(i = 0; i < nu; ++i)
			so[i] = fabs(so[i] - BS[i]) / BS[i];

		info = convtest(m_b, T, so, nu, tol * Smax, tol, &k);
/*		printf("iter %d: sv[%d] = %f, k = %d -> converged %d\n", iter, nu, BS[nu], k, info == 0);*/

		F77_NAME(dcopy)(&nu, BS, &inc, so, &inc);
		if( info == 0 )
			break;

		/* Update the right approximate singular vectors. */
		F77_NAME(dgemm)("N", "T", &n, &k, &m_b, &alpha, V, &n, BV, &m_b, &beta, V1, &n);
		F77_NAME(dlacpy)("N", &n, &k, V1, &n, V, &n);
		F77_NAME(dcopy)(&n, F, &inc, V + k * n, &inc);

		/* Update the left approximate singular vectors. */
		F77_NAME(dgemm)("N", "N", &m, &k, &m_b, &alpha, W, &m, BU, &m_b, &beta, U1, &m);
		F77_NAME(dlacpy)("N", &m, &k, U1, &m, W, &m);

		/* Update B matrix */
		memset(B, 0, m_b * m_b * sizeof(double));
		for(i = 0; i < k; ++i) {
			B[i * m_b + i] = BS[i];
			B[k * m_b + i] = T[i];
		}

		iter++;
	}

	if(info == 0) {
		F77_NAME(dgemm)("N", "N", &m, &nu, &m_b, &alpha, W, &m, BU, &m_b, &beta, Uo, &m);
		F77_NAME(dgemm)("N", "T", &n, &nu, &m_b, &alpha, V, &n, BV, &m_b, &beta, Vo, &n);
	}

	return info;
}

blas_int irlba(blas_int m, blas_int n, blas_int nu, blas_int m_b, blas_int maxit, blas_int restart, double tol,
		double *init, // input starting vector (n) (can be NULL)
		double *so, // output singular values
		double *Uo, // output left singular vectors
		double *Vo, // output right singular vectors
		randn_t randn, matmul_t matmul, void *data) {

	double *V = malloc(n * m_b * sizeof(double));
	double *W = malloc(m * m_b * sizeof(double));
	double *F = malloc(n * sizeof(double));
	double *T = malloc(n * sizeof(double));
	double *B = malloc(m_b * m_b * sizeof(double));
	double *BS = malloc(m_b * sizeof(double));
	double *BU = malloc(m_b * m_b * sizeof(double));
	double *BV = malloc(m_b * m_b * sizeof(double));
	double *U1 = malloc(m * m_b * sizeof(double));
	double *V1 = malloc(n * m_b * sizeof(double));
	blas_int *iwork = malloc(8 * m_b * sizeof(blas_int));

	blas_int lwork = prep_svd(m_b, m_b, B, BS, BU, BV, iwork);
	double *work = malloc(lwork * sizeof(double));

	blas_int inc = 1;
	blas_int info;
	double S;
	int i;

	memset(B, 0, m_b * m_b * sizeof(double));

	if( restart > 0 ) {
		F77_NAME(dlacpy)("N", &m, &restart, Vo, &n, V, &n);
		F77_NAME(dlacpy)("N", &m, &restart, Uo, &m, W, &m);
		for(i = 0; i < restart; ++i) {
			B[i * m_b + i] = so[i];
		}
	}

	if( init != NULL )
		F77_NAME(dcopy)(&n, init, &inc, V + restart * n, &inc);
	else
		randn(V + restart * n, n, data);

	/* normalize starting vector */
	S = F77_NAME(dnrm2)(&n, V + restart * n, &inc);
	if(S < eps())
		return -1;
	S = 1 / S;
	F77_NAME(dscal)(&n, &S, V + restart * n, &inc);

	info = irlba_(m, n, nu, m_b, maxit, restart, tol, so, Uo, Vo,
		V, W, F, T, B, BS, BU, BV, U1, V1, work, lwork, iwork,
		randn, matmul, data);

	free(work);
	free(iwork);
	free(V1);
	free(U1);
	free(BV);
	free(BU);
	free(BS);
	free(B);
	free(T);
	free(F);
	free(W);
	free(V);

	return info;
}

