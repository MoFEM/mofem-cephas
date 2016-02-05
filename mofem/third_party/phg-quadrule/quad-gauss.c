/* Parallel Hierarchical Grid -- an adaptive finite element library.
 *
 * Copyright (C) 2005-2010 State Key Laboratory of Scientific and
 * Engineering Computing, Chinese Academy of Sciences. */

/* This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA */

/* $Id: quad-gauss.c,v 1.20 2014/04/04 04:30:02 zlb Exp $
 *
 * Function for generating tensor product rules using 1D Gauss-Jacobi rules */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <quad.h>

#define DBL_EPSILON 2.2204460492503131E-16

// #include "sturm.h"

#if 0
static double
Gamma(double x)
/* Computes Gamma(x) for real x, x != 0, -1, -2, ...
 * converted from http://jin.ece.uiuc.edu/routines/mgamma.for */
{
    static double g[] = {1.0, 0.5772156649015329, -0.6558780715202538,
	-0.420026350340952e-1, 0.1665386113822915,-.421977345555443e-1,
	-.96219715278770e-2, .72189432466630e-2, -.11651675918591e-2,
	-.2152416741149e-3, .1280502823882e-3, -.201348547807e-4,
	-.12504934821e-5, .11330272320e-5, -.2056338417e-6, .61160950e-8,
	.50020075e-8, -.11812746e-8, .1043427e-9, .77823e-11,
	-.36968e-11, .51e-12, -.206e-13, -.54e-14, .14e-14, .1e-15};
    int k, m;
    double ga, gr, z, r;

    if (x == (double)(m = (int)x)) {
	if (x > 0) {
	    ga = 1.;
	    for (k = 2; k < m; k++)
		ga *= (double)k;
	}
	else {
	    ga = 1.e300;
	}
    }
    else {
	if ((z = fabs(x)) > 1.) {
	    m = (int)z;
	    r = 1.0;
	    for (k = 1; k <= m; k++)
		r *= z - k;
	    z -= m;
	}
	else {
	    z = x;
	}
	gr = g[25];
	for (k = 24; k >= 0; k--)
	    gr = gr * z + g[k];
	ga = 1. / (gr * z);
	if (fabs(x) > 1.) {
	    ga *= r;
	    if (x < 0.)
		ga = -M_PI / (x * ga * sin(M_PI * x));
	}
    }

    return ga;
}
#else
/* Note: the tgamma*() functions seem to be missing in math.h of gcc-4.1.x */
extern float tgammaf(float x);
extern double tgamma(double x);
extern long double tgammal(long double x);
#endif

static double
peval(int n, double *p, double x)
/* evaluates the value of p(x) */
{
    int i;
    double d;

    d = p[0];
    for (i = 1; i <= n; i++)
	d = d * x + p[i];
    return (d);
}

static int
comp_FLOAT(const void *p0, const void *p1)
{
    double d = *((const double *)p0) - *((const double *)p1);
    return  d > 0 ? 1 : (d < 0 ? -1 : 0);
}

static void
gauss_jacobi(double x[], double w[], int n, double alpha, double beta)
/* compute abscissas and weights of the n-point 1D Gauss-Jacobi rule */
{
    double p0[n + 1], p1[n + 1], p[n + 1];
    double c, d, e, f, temp0, temp;
    int i, j, iwk[n + 1];

    if (n == 0)
	return;

    /* evaluate the n-th order Gauss-Jacobi polynomial */
    p1[0] = 1.0;
    p[0] = (2.0 + alpha + beta) * 0.5;
    p[1] = (alpha - beta) * 0.5;
    if (n > 1) {
	temp0 = alpha * alpha - beta * beta;
	for (j = 1; j < n; j++) {
	    memcpy(p0, p1, j * sizeof(p[0]));
	    memcpy(p1, p, (j + 1) * sizeof(p[0]));
	    temp = 2. * j + alpha + beta + 1.;
	    c = 2. * (j + 1.) * (temp - j) * (temp - 1.);
	    d = temp * temp0 / c;
	    e = (temp * temp - 1.) * temp / c;
	    f = 2. * (j + alpha) * (j + beta) * (temp + 1.) / c;
	    /* P = (d + e x) * P1 - f * P0 */
	    p[0] = e * p1[0];
	    p[1] = e * p1[1] + d * p1[0];
	    p1[j + 1] = 0.0;
	    for (i = 2; i <= j + 1; i++)
		p[i] = e * p1[i] + d * p1[i - 1] - f * p0[i - 2];
	}
    }
    /* save p'(x) in p0 */
    for (i = 0; i < n; i++)
	p0[i] = (n - i) * p[i];
    /* compute the roots of p[] between (-1, 1) */
    if (alpha == beta) {
	/* compute last half of the roots, evaluate the others by symmetry */
	j = phgQuadSturm(p, n, -1., 0., DBL_EPSILON, NULL, iwk);
	assert(j == (n + 1) / 2);
	qsort(p, j, sizeof(p[0]), comp_FLOAT);
	for (i = 0; i < j; i++) {
	    p[n - 1 - i] = -p[i];
	    iwk[n - 1 - i] = iwk[i];
	}
    }
    else {
	j = phgQuadSturm(p, n, -1., 1., DBL_EPSILON, NULL, iwk);
	assert(j == n);
    }
    temp = Gamma(alpha + n) * Gamma(beta + n) /
	   (Gamma(n + 1.) * Gamma(n + alpha + beta + 1.)) *
	   (n + n + alpha + beta) * Pow((double)2., alpha + beta);
    for (i = 0; i < n; i++) {
	assert(iwk[i] == 1);
	x[i] = p[i];
	w[i] = temp / (peval(n - 1, p0, x[i]) * peval(n - 1, p1, x[i]));
    }
}

#if HAVE_LIBGMP	/*-------------------------------------------------------*/

#include <gmp.h>

static double
gmp_eval(int N, mpf_t X[], mpf_t A[], mpf_t B[], int alpha)
/* computes the residual (B) and the Jacobian (A), returns max(residual) */
{
    int n = N / 2;
    int i, k;
    mpf_t a, b, c, d, error;
    double err;

    mpf_init(a);
    mpf_init(b);
    mpf_init(c);
    mpf_init(d);
    mpf_init(error);

    /* TODO: avoid calling mpf_pow_ui by swapping the i and k loop */

    /* Loop over polynomials (1-x)^k, k = 0, ..., N - 1 */
    for (k = 0; k < N; k++) {
	/* Note: \int (1-x)^(alpha+k) = 2^(1+alpha+k) / (1+alpha+k) */
	/* compute d = 1 / \int (1-x)^(alpha+k) */
	mpf_set_ui(d, 2);
	mpf_pow_ui(d, d, alpha + k + 1);
	mpf_ui_div(d, alpha + k + 1, d);
	/* the residual */
	/* compute d * sum (1-x_i)^k * w_i, i = 0, ..., N/2-1 */
	mpf_set_ui(a, 0);
	for (i = 0; i < n; i++) {
	    mpf_ui_sub(c, 1.0, X[i]);	/* c = 1-x_i */
	    /* dF/dw_i = (1-x_i)^k * d */
	    mpf_pow_ui(b, c, k);	/* b = (1-x_i)^k */
	    mpf_mul(A[k * N + n + i], b, d);	/* J = d * (1-x_i)^k */
	    /* multiply by d */
	    mpf_mul(b, A[k * N + n + i], X[n + i]); /* b = d * (1-x_i)^k w_i */
	    mpf_add(a, a, b);			/* add to sum */
	    /* dF/dx_i = -d*k*(1-x_i)^(k-1) * w_i */
	    if (k == 0) {
		mpf_set_ui(A[k * N + i], 0);
	    }
	    else {
#if 0
		mpf_pow_ui(b, c, k - 1);
		mpf_set_si(c, -k);
		mpf_mul(b, b, c);
		mpf_mul(b, b, d);
		mpf_mul(A[k * N + i], b, X[n + i]);
#else
		mpf_div(b, b, c);	/* b = d * (1-x_i)^(k-1) * w_i */
		mpf_set_si(c, -k);
		mpf_mul(A[k * N + i], b, c);
#endif
	    }
	}
	mpf_sub_ui(B[k], a, 1);
	mpf_abs(b, B[k]);
	if (mpf_cmp(error, b) < 0)
	    mpf_set(error, b);
    }

    mpf_clear(a);
    mpf_clear(b);
    mpf_clear(c);
    mpf_clear(d);
    err = (double)mpf_get_d(error);
    mpf_clear(error);

    return err;
}

static BOOLEAN
gmp_solve(int N, mpf_t A[], mpf_t B[])
/* solves the linear system of equations AX = B by Gauss elimination */
{
    int i, j, k;
    mpf_t d, r;

    mpf_init(d);
    mpf_init(r);

    for (i = 0; i < N; i++) {
	mpf_abs(d, A[i * N + i]);
	k = i;
	for (j = i + 1; j < N; j++) {
	    mpf_abs(r, A[j * N + i]);
	    if (mpf_cmp(r, d) > 0) {
		mpf_set(d, r);
		k = j;
	    }
	}
	if (mpf_cmp_ui(d, 0) == 0) {
#if 0
	    printf("Dense solver failed: i = %d\n", i);
	    for (j = 0; j < N; j++) {
		for (k = 0; k < N; k++)
		    gmp_printf(" %9.2Fe", A[j * N + k]);
		gmp_printf(" = %9.2Fe\n", B[j]);
	    }
#endif
	    mpf_clear(d);
	    mpf_clear(r);
	    return FALSE;
	}
	if (k != i) {
	    /* exchange row i and row k */
	    for (j = i; j < N; j++)
		mpf_swap(A[i * N + j], A[k * N + j]);
	    mpf_swap(B[i], B[k]);
	}
	if (mpf_cmp_ui(A[i * N + i], 1) != 0) {
	    mpf_ui_div(d, 1, A[i * N + i]);
	    for (j = i + 1; j < N; j++)
		mpf_mul(A[i * N + j], A[i * N + j], d);
	    mpf_mul(B[i], B[i], d);
	}
	for (j = i + 1; j < N; j++) {
	    if (mpf_cmp_ui(A[j * N + i], 0) == 0)
		continue;
	    for (k = i + 1; k < N; k++) {
		mpf_mul(d, A[j * N + i], A[i * N + k]);
		mpf_sub(A[j * N + k], A[j * N + k], d);
	    }
	    mpf_mul(d, A[j * N + i], B[i]);
	    mpf_sub(B[j], B[j], d);
	}
    }

    for (i = N - 2; i >= 0; i--) {
	for (j = i + 1; j < N; j++) {
	    mpf_mul(d, A[i * N + j], B[j]);
	    mpf_sub(B[i], B[i], d);
	}
    }

    mpf_clear(d);
    mpf_clear(r);

    return TRUE;
}

static void
gmp_newton(int n, double x[], double w[], int alpha)
/* performs Newton iterations using high precision arithmetics to improve
 * precision of an n-point Jacobi-Gauss quadrature rule (beta = 0) */
{
    int N = n + n;
    mpf_t *A, *B, *X;
    int i, j, it;
    double d = 0.;

    A = phgAlloc(N * N * sizeof(*A));
    B = phgAlloc(N * sizeof(*B));
    X = phgAlloc(N * sizeof(*X));

    mpf_set_default_prec(sizeof(double) * 8 * 4);

    for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++)
	    mpf_init(A[i * N + j]);
	mpf_init(B[i]);
	mpf_init(X[i]);
	mpf_set_d(X[i], (double)(i < n ? x[i] : w[i - n]));
    }

    /* Newton iterations */
    if (phgRank == 0)
	phgInfo(2, "gmp_newton: alpha = %d\n", alpha);
    for (it = 0; it < 8; it++) {
	d = gmp_eval(N, X, A, B, alpha);
	if (phgRank == 0)
	    phgInfo(2, "gmp_newton: it = %d, error = %lg\n", it, (double)d);
	if (d < 1e-40 || d < (1e4 * FLOAT_MIN))
	    break;
	if (!gmp_solve(N, A, B))
	    phgError(1, "%s:%d: unexpected error!\n", __FILE__, __LINE__);
	/* xnew = xold - J^(-1)F(xold) (B = J^(-1)F(xold)) */
	for (i = 0; i < N; i++)
	    mpf_sub(X[i], X[i], B[i]);
    }

    if (d >= 1e-40 && (d >= 1e4 * FLOAT_MIN))
	phgError(1, "%s:%d: unexpected error!\n", __FILE__, __LINE__);

    for (i = 0; i < N; i++) {
	d = (double)mpf_get_d(X[i]);
	if (sizeof(double) > sizeof(double)) {
	    mpf_set_d(B[i], (double)d);
	    mpf_sub(X[i], X[i], B[i]);
	    d = ((double)((double)d)) + (double)mpf_get_d(X[i]);
	}
	i < n ? (x[i] = d) : (w[i - n] = d);
	for (j = 0; j < N; j++)
	    mpf_clear(A[i * N + j]);
	mpf_clear(B[i]);
	mpf_clear(X[i]);
    }

    phgFree(A);
    phgFree(B);
    phgFree(X);
}

#else	/*--------------------------- HAVE_LIBGMP ---------------------------*/

static void
gmp_newton(int n, double x[], double w[], int alpha)
{
    Unused(n);
    Unused(x);
    Unused(w);
    Unused(alpha);
}

#endif	/*--------------------------- HAVE_LIBGMP ---------------------------*/

QUAD *
phgQuadTensorProductRule(int dim, int order)
/* This code computes a 1D, 2D or 3D quadrature rule of order 'order'
 * (2D or 3D rules are constructed using tensor product of 1D
 * Gauss-Jacobi rules).
 *
 * Note: let y = (1-x) y', z = (1-x-y) z', then:
 *	\int_T f(x,y,z) dx =
 *		\int_0^1\int_0^1\int_0^1 (1-x)^2 (1-y') f(x,y',z') dx dy' dz'
 * thus the integral is approximated by:
 *		Gauss-Jacobi(2,0) + Gauss-Jacobi(1,0) + Gauss-Legendre */
{
    int i, j, k, n, npoints;
    double x, y, z;
    double *x0, *w0, *x1 = NULL, *w1 = NULL, *x2 = NULL, *w2 = NULL;
    double *lambda, *weight;
    QUAD *quad;
    char s[16];

    assert(order > 0 && dim > 0 && dim <= 3);

    n = order / 2 + 1;
    x0 = malloc(n * sizeof(*x0));
    w0 = malloc(n * sizeof(*w0));

    if (dim > 1) {
	x1 = malloc(n * sizeof(*x1));
	w1 = malloc(n * sizeof(*w1));
	if (dim > 2) {
	    x2 = malloc(n * sizeof(*x2));
	    w2 = malloc(n * sizeof(*w2));
	    npoints = n * n * n;
	}
	else {
	    npoints = n * n;
	}
    }
    else {
	npoints = n;
    }

    quad = phgAlloc(sizeof(*quad));
    sprintf(s, "%dD P%dt", dim, n + n - 1);
    quad->name = strdup(s);
    quad->dim = dim;
    quad->order = n + n - 1;
    quad->npoints = npoints;
    quad->points = phgAlloc(npoints * (dim + 1) * sizeof(*quad->points));
    quad->weights = phgAlloc(npoints * sizeof(*quad->weights));
    quad->id = -1;

    lambda = quad->points;
    weight = quad->weights;

    /* 1D Gauss-Legendre rule, which is the same as Gauss-Jacobi(0,0) */
    gauss_jacobi(x0, w0, n, 0.0, 0.0);
    /* perform a few Newton iterations using high precision arithmetics */
    gmp_newton(n, x0, w0, 0);
    /* scale to (0, 1) */
    for (i = 0; i < n; i++) {
	x0[i] = .5 * (x0[i] + 1.);
	w0[i] *= .5;
    }

    if (dim > 1) {
	/* 1D Gauss-Jacobi(1,0) rule
	 * Note: let x = 2y - 1, then
	 *		\int_{-1}^1 (1-x) f(x) dx = 4 \int (1-y) f(2y-1) dy
	 * Thus x[i] and w[i] should be scaled by:
	 *		x[i] := 0.5 * (x[i] + 1.0);
	 *		w[i] := 0.25 * w[i] */
	gauss_jacobi(x1, w1, n, 1.0, 0.0);
	gmp_newton(n, x1, w1, 1);
	for (i = 0; i < n; i++) {
	    x1[i] = .5 * (x1[i] + 1.);
	    w1[i] *= .25;
	}

	if (dim > 2) {
	    /* 1D Gauss-Jacobi(2,0) rule
	     * Note: let x = 2y - 1, then
	     *	\int_{-1}^1 (1-x)^2 f(x) dx = 8 \int (1-y)^2 f(2y-1) dy
	     * Thus x[i] and w[i] should be scaled by:
	     *		x[i] := 0.5 * (x[i] + 1.0);
	     *		w[i] := 0.125 * w[i] */
	    gauss_jacobi(x2, w2, n, 2.0, 0.0);
	    gmp_newton(n, x2, w2, 2);
	    for (i = 0; i < n; i++) {
		x2[i] = _F(.5) * (x2[i] + 1.);
		w2[i] *= _F(.125);
	    }
	}
    }

    if (dim == 3) {
	for (i = 0; i < n; i++) {
	    for (j = 0; j < n; j++) {
		for (k = 0; k < n; k++) {
		    x = x2[k];
		    y = x1[j] * (1. - x2[k]);
		    z = x0[i] * (1. - x2[k]) * (1. - x1[j]);
		    *(lambda++) = 1. - (x + y + z);
		    *(lambda++) = x;
		    *(lambda++) = y;
		    *(lambda++) = z;
		    *(weight++) = 6. * w0[i] * w1[j] * w2[k];
		}
	    }
	}
    }
    else if (dim == 2) {
	for (i = 0; i < n; i++) {
	    for (j = 0; j < n; j++) {
		x = x1[j];
		y = x0[i] * (1. - x1[j]);
		*(lambda++) = 1. - (x + y);
		*(lambda++) = x;
		*(lambda++) = y;
		*(weight++) = 2. * w0[i] * w1[j];
	    }
	}
    }
    else {
	for (i = 0; i < n; i++) {
		x = x0[i];
		*(lambda++) = 1. - x;
		*(lambda++) = x;
		*(weight++) = w0[i];
	}
    }

    if (dim > 1) {
	if (dim > 2) {
	    free(x2);
	    free(w2);
	}
	free(x1);
	free(w1);
    }
    free(x0);
    free(w0);

    return quad;
}

#ifdef TEST

int
main(int argc, char *argv[])
{
    int i, l, order, dim, npoints;
    QUAD *quad;
    double sum;

    if (argc != 3) {
usage:
	fprintf(stderr, "Usage: %s dim order\n", argv[0]);
	exit(1);
    }

    order = atoi(argv[2]);
    dim = atoi(argv[1]);

    if (order <= 0 || dim <= 0 || dim > 3) {
	fprintf(stderr, "Invalid arguments.\n");
	goto usage;
    }

    quad = phgQuadTensorProductRule(dim, order);
    npoints = quad->npoints;

    for (i = 0, sum = 0.; i < npoints; i++) {
	sum += quad->weights[i];
    }

    fprintf(stderr, "dimension: %d\n", dim);
    fprintf(stderr, "number of points: %d\n", quad->npoints);
    fprintf(stderr, "algebraic order: %d\n", order);
    fprintf(stderr, "sum of weights:  %0.16le\n", (double)sum);

#define LAMBDA(l,i)	(double)quad->points[l * (dim + 1) + i]
#define WEIGHTS(l)	(double)quad->weights[l]

    /* output the rule */
    printf("\nstatic double QUAD_%dD_P%d_wts[] = {\n", dim, order);

    if (dim == 1)
	npoints = (quad->npoints + 1) / 2;	/* symmetry */

    for (l = 0; l < npoints; l++) {
	if (dim > 1)
	    printf("    Dup0(%0.16le)%s\n",
			WEIGHTS(l), l + 1 == npoints ? "" : ",");
	else
	    printf(2 * l == quad->npoints - 1 ? "    Dup2(%0.16le)%s\n" :
				    "    Dup11(%0.16le)%s\n",
			WEIGHTS(l), l + 1 == npoints ? "" : ",");
    }
    printf("};\nstatic double QUAD_%dD_P%d_pts[Length(QUAD_%dD_P%d_wts) "
	   "* %d] = {\n", dim, order, dim, order, dim + 1);

    for (l = 0; l < npoints; l++) {
	switch (dim) {
	    case 3:
		printf("    Perm0(%0.16le,%0.16le,%0.16le)%s\n",
			LAMBDA(l,0), LAMBDA(l,1), LAMBDA(l,2),
			l + 1 == npoints ? "" : ",");
		break;
	    case 2:
		printf("    Perm0(%0.16le, %0.16le)%s\n",
			LAMBDA(l,0), LAMBDA(l,1),
			l + 1 == npoints ? "" : ",");
		break;
	    case 1:
		printf(2 * l == quad->npoints - 1 ? "    Perm2(%0.16le)%s\n" :
					"    Perm11(%0.16le)%s\n",
			LAMBDA(l,0), l + 1 == npoints ? "" : ",");
		break;
	}
    }
    printf("};\nQUAD QUAD_%dD_P%d_ = {\n", dim, order);
    printf("    \"%dD P%d\",\t\t\t/* name */\n", dim, order);
    printf("    %d,\t\t\t\t/* dim */\n", dim);
    printf("    %d,\t\t\t\t/* order */\n", order);
    printf("    Length(QUAD_1D_P%d_wts),\t/* npoints (%d) */\n",
		order, quad->npoints);
    printf("    QUAD_%dD_P%d_pts,\t\t/* points */\n", dim, order);
    printf("    QUAD_%dD_P%d_wts,\t\t/* weights */\n", dim, order);
    printf("    -1\t\t\t\t/* id */\n};\n\n");

    phgQuadFree(&quad);

    return 0;
}

#endif	/* defined(TEST) */
