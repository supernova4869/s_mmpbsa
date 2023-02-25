/**
 *  @file       vblas.h
 *  @brief      The primary header for the BLAS.
 *              (Basic Linear Algebra Subroutines.)
 *  @note       We provide this header whether or not we provide the BLAS
 *              library itself.  This gives us some compile-time type-checking
 *              even for architecture-dependent assembly-coded BLAS.
 *  @version    $Id: vblas.h,v 1.15 2010/08/12 05:52:32 fetk Exp $
 *  @author     Michael Holst
 *  
 *  @attention
 *  @verbatim
 *
 * PUNC = < Portable Understructure for Numerical Computing >
 * Copyright (C) 1994-- Michael Holst
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * 
 *  @endverbatim
 */


#ifndef _VBLAS_H_
#define _VBLAS_H_

#include <punc/punc_base.h>

#include <punc/vf2c.h>

/*
 * ***************************************************************************
 * Library VBLAS prototypes
 * ***************************************************************************
 */

/** @brief Library VBLAS prototypes 
    @author Michael Holst */
VEXTERNC logical lsame_(char *ca, char *cb);
/** @brief Library VBLAS prototypes 
    @author Michael Holst */
VEXTERNC int xerbla_(char *srname, integer *info);
/** @brief Library VBLAS prototypes 
    @author Michael Holst */
VEXTERNC int caxpy_(integer *n, realcomplex *ca, realcomplex *cx, integer *incx, realcomplex *cy, integer *incy);
/** @brief Library VBLAS prototypes 
    @author Michael Holst */
VEXTERNC int ccopy_(integer *n, realcomplex *cx, integer *incx, realcomplex *cy, integer *incy);
/** @brief Library VBLAS prototypes 
    @author Michael Holst */
VEXTERNC C_f cdotc_(realcomplex * ret_val, integer *n, realcomplex *cx, integer *incx, realcomplex *cy, integer *incy);
/** @brief Library VBLAS prototypes 
    @author Michael Holst */
VEXTERNC C_f cdotu_(realcomplex * ret_val, integer *n, realcomplex *cx, integer *incx, realcomplex *cy, integer *incy);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124. \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int cgbmv_(char *trans, integer *m, integer *n, integer *kl, integer *ku, realcomplex *alpha, realcomplex *a, integer *lda, realcomplex *x, integer *incx, realcomplex *beta, realcomplex *y, integer *incy, ftnlen trans_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int cgemm_(char *transa, char *transb, integer *m, integer *n, integer *k, realcomplex *alpha, realcomplex *a, integer *lda, realcomplex *b, integer *ldb, realcomplex *beta, realcomplex *c__, integer *ldc, ftnlen transa_len, ftnlen transb_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int cgemv_(char *trans, integer *m, integer *n, realcomplex *alpha, realcomplex *a, integer *lda, realcomplex *x, integer *incx, realcomplex *beta, realcomplex *y, integer *incy, ftnlen trans_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int cgerc_(integer *m, integer *n, realcomplex *alpha, realcomplex *x, integer *incx, realcomplex *y, integer *incy, realcomplex *a, integer *lda);
/** @brief Library VBLAS prototypes. \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int cgeru_(integer *m, integer *n, realcomplex *alpha, realcomplex *x, integer *incx, realcomplex *y, integer *incy, realcomplex *a, integer *lda);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int chbmv_(char *uplo, integer *n, integer *k, realcomplex *alpha, realcomplex *a, integer *lda, realcomplex *x, integer *incx, realcomplex *beta, realcomplex *y, integer *incy, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int chemm_(char *side, char *uplo, integer *m, integer *n, realcomplex *alpha, realcomplex *a, integer *lda, realcomplex *b, integer *ldb, realcomplex *beta, realcomplex *c__, integer *ldc, ftnlen side_len, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int chemv_(char *uplo, integer *n, realcomplex *alpha, realcomplex *a, integer *lda, realcomplex *x, integer *incx, realcomplex *beta, realcomplex *y, integer *incy, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int cher_(char *uplo, integer *n, real *alpha, realcomplex *x, integer *incx, realcomplex *a, integer *lda, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int cher2_(char *uplo, integer *n, realcomplex *alpha, realcomplex *x, integer *incx, realcomplex *y, integer *incy, realcomplex *a, integer *lda, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int cher2k_(char *uplo, char *trans, integer *n, integer *k, realcomplex *alpha, realcomplex *a, integer *lda, realcomplex *b, integer *ldb, real *beta, realcomplex *c__, integer *ldc, ftnlen uplo_len, ftnlen trans_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int cherk_(char *uplo, char *trans, integer *n, integer *k, real *alpha, realcomplex *a, integer *lda, real *beta, realcomplex *c__, integer *ldc, ftnlen uplo_len, ftnlen trans_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int chpmv_(char *uplo, integer *n, realcomplex *alpha, realcomplex *ap, realcomplex *x, integer *incx, realcomplex *beta, realcomplex *y, integer *incy, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int chpr_(char *uplo, integer *n, real *alpha, realcomplex *x, integer *incx, realcomplex *ap, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int chpr2_(char *uplo, integer *n, realcomplex *alpha, realcomplex *x, integer *incx, realcomplex *y, integer *incy, realcomplex *ap, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. 
    @author Michael Holst */
VEXTERNC int crotg_(realcomplex *ca, realcomplex *cb, real *c__, realcomplex *s);
/** @brief Library VBLAS prototypes. 
    @author Michael Holst */
VEXTERNC int cscal_(integer *n, realcomplex *ca, realcomplex *cx, integer *incx);
/** @brief Library VBLAS prototypes. 
    @author Michael Holst */
VEXTERNC int csscal_(integer *n, real *sa, realcomplex *cx, integer *incx);
/** @brief Library VBLAS prototypes. 
    @author Michael Holst */
VEXTERNC int cswap_(integer *n, realcomplex *cx, integer *incx, realcomplex *cy, integer *incy);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int csymm_(char *side, char *uplo, integer *m, integer *n, realcomplex *alpha, realcomplex *a, integer *lda, realcomplex *b, integer *ldb, realcomplex *beta, realcomplex *c__, integer *ldc, ftnlen side_len, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int csyr2k_(char *uplo, char *trans, integer *n, integer *k, realcomplex *alpha, realcomplex *a, integer *lda, realcomplex *b, integer *ldb, realcomplex *beta, realcomplex *c__, integer *ldc, ftnlen uplo_len, ftnlen trans_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int csyrk_(char *uplo, char *trans, integer *n, integer *k, realcomplex *alpha, realcomplex *a, integer *lda, realcomplex *beta, realcomplex *c__, integer *ldc, ftnlen uplo_len, ftnlen trans_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int ctbmv_(char *uplo, char *trans, char *diag, integer *n, integer *k, realcomplex *a, integer *lda, realcomplex *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int ctbsv_(char *uplo, char *trans, char *diag, integer *n, integer *k, realcomplex *a, integer *lda, realcomplex *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int ctpmv_(char *uplo, char *trans, char *diag, integer *n, realcomplex *ap, realcomplex *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int ctpsv_(char *uplo, char *trans, char *diag, integer *n, realcomplex *ap, realcomplex *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int ctrmm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, realcomplex *alpha, realcomplex *a, integer *lda, realcomplex *b, integer *ldb, ftnlen side_len, ftnlen uplo_len, ftnlen transa_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int ctrmv_(char *uplo, char *trans, char *diag, integer *n, realcomplex *a, integer *lda, realcomplex *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int ctrsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, realcomplex *alpha, realcomplex *a, integer *lda, realcomplex *b, integer *ldb, ftnlen side_len, ftnlen uplo_len, ftnlen transa_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int ctrsv_(char *uplo, char *trans, char *diag, integer *n, realcomplex *a, integer *lda, realcomplex *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. 
    @author Michael Holst */
VEXTERNC doublereal dasum_(integer *n, doublereal *dx, integer *incx);
/** @brief Library VBLAS prototypes. 
    @author Michael Holst */
VEXTERNC int daxpy_(integer *n, doublereal *da, doublereal *dx, integer *incx, doublereal *dy, integer *incy);
/** @brief Library VBLAS prototypes. 
    @author Michael Holst */
VEXTERNC doublereal dcabs1_(doublecomplex *z__);
/** @brief Library VBLAS prototypes. 
    @author Michael Holst */
VEXTERNC int dcopy_(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy);
/** @brief Library VBLAS prototypes. 
    @author Michael Holst */
VEXTERNC doublereal ddot_(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int dgbmv_(char *trans, integer *m, integer *n, integer *kl, integer *ku, doublereal *alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, doublereal *beta, doublereal *y, integer *incy, ftnlen trans_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int dgemm_(char *transa, char *transb, integer *m, integer *n, integer *k, doublereal *alpha, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *beta, doublereal *c__, integer *ldc, ftnlen transa_len, ftnlen transb_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int dgemv_(char *trans, integer *m, integer *n, doublereal *alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, doublereal *beta, doublereal *y, integer *incy, ftnlen trans_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int dger_(integer *m, integer *n, doublereal *alpha, doublereal *x, integer *incx, doublereal *y, integer *incy, doublereal *a, integer *lda);
/** @brief Library VBLAS prototypes. 
    @author Michael Holst */
VEXTERNC doublereal dnrm2_(integer *n, doublereal *x, integer *incx);
/** @brief Library VBLAS prototypes. 
    @author Michael Holst */
VEXTERNC int drot_(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy, doublereal *c__, doublereal *s);
/** @brief Library VBLAS prototypes. 
    @author Michael Holst */
VEXTERNC int drotg_(doublereal *da, doublereal *db, doublereal *c__, doublereal *s);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int dsbmv_(char *uplo, integer *n, integer *k, doublereal *alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, doublereal *beta, doublereal *y, integer *incy, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. 
    @author Michael Holst */
VEXTERNC int dscal_(integer *n, doublereal *da, doublereal *dx, integer *incx);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int dspmv_(char *uplo, integer *n, doublereal *alpha, doublereal *ap, doublereal *x, integer *incx, doublereal *beta, doublereal *y, integer *incy, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int dspr_(char *uplo, integer *n, doublereal *alpha, doublereal *x, integer *incx, doublereal *ap, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int dspr2_(char *uplo, integer *n, doublereal *alpha, doublereal *x, integer *incx, doublereal *y, integer *incy, doublereal *ap, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. 
    @author Michael Holst */
VEXTERNC int dswap_(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int dsymm_(char *side, char *uplo, integer *m, integer *n, doublereal *alpha, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *beta, doublereal *c__, integer *ldc, ftnlen side_len, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int dsymv_(char *uplo, integer *n, doublereal *alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, doublereal *beta, doublereal *y, integer *incy, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int dsyr_(char *uplo, integer *n, doublereal *alpha, doublereal *x, integer *incx, doublereal *a, integer *lda, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int dsyr2_(char *uplo, integer *n, doublereal *alpha, doublereal *x, integer *incx, doublereal *y, integer *incy, doublereal *a, integer *lda, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int dsyr2k_(char *uplo, char *trans, integer *n, integer *k, doublereal *alpha, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *beta, doublereal *c__, integer *ldc, ftnlen uplo_len, ftnlen trans_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int dsyrk_(char *uplo, char *trans, integer *n, integer *k, doublereal *alpha, doublereal *a, integer *lda, doublereal *beta, doublereal *c__, integer *ldc, ftnlen uplo_len, ftnlen trans_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int dtbmv_(char *uplo, char *trans, char *diag, integer *n, integer *k, doublereal *a, integer *lda, doublereal *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int dtbsv_(char *uplo, char *trans, char *diag, integer *n, integer *k, doublereal *a, integer *lda, doublereal *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int dtpmv_(char *uplo, char *trans, char *diag, integer *n, doublereal *ap, doublereal *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int dtpsv_(char *uplo, char *trans, char *diag, integer *n, doublereal *ap, doublereal *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int dtrmm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, doublereal *alpha, doublereal *a, integer *lda, doublereal *b, integer *ldb, ftnlen side_len, ftnlen uplo_len, ftnlen transa_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int dtrmv_(char *uplo, char *trans, char *diag, integer *n, doublereal *a, integer *lda, doublereal *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int dtrsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, doublereal *alpha, doublereal *a, integer *lda, doublereal *b, integer *ldb, ftnlen side_len, ftnlen uplo_len, ftnlen transa_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int dtrsv_(char *uplo, char *trans, char *diag, integer *n, doublereal *a, integer *lda, doublereal *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: dcabs1_ 7 1 9
    @author Michael Holst  */
VEXTERNC doublereal dzasum_(integer *n, doublecomplex *zx, integer *incx);
/** @brief Library VBLAS prototypes.
    @author Michael Holst  */
VEXTERNC doublereal dznrm2_(integer *n, doublecomplex *x, integer *incx);
/** @brief Library VBLAS prototypes.
    @author Michael Holst  */
VEXTERNC integer icamax_(integer *n, realcomplex *cx, integer *incx);
/** @brief Library VBLAS prototypes.
    @author Michael Holst  */
VEXTERNC integer idamax_(integer *n, doublereal *dx, integer *incx);
/** @brief Library VBLAS prototypes.
    @author Michael Holst  */
VEXTERNC integer isamax_(integer *n, real *sx, integer *incx);
/** @brief Library VBLAS prototypes. \n
 *:ref: dcabs1_ 7 1 9
    @author Michael Holst  */
VEXTERNC integer izamax_(integer *n, doublecomplex *zx, integer *incx);
/** @brief Library VBLAS prototypes.
    @author Michael Holst  */
VEXTERNC E_f sasum_(integer *n, real *sx, integer *incx);
/** @brief Library VBLAS prototypes.
    @author Michael Holst  */
VEXTERNC int saxpy_(integer *n, real *sa, real *sx, integer *incx, real *sy, integer *incy);
/** @brief Library VBLAS prototypes.
    @author Michael Holst  */
VEXTERNC E_f scasum_(integer *n, realcomplex *cx, integer *incx);
/** @brief Library VBLAS prototypes.
    @author Michael Holst  */
VEXTERNC E_f scnrm2_(integer *n, realcomplex *x, integer *incx);
/** @brief Library VBLAS prototypes.
    @author Michael Holst  */
VEXTERNC int scopy_(integer *n, real *sx, integer *incx, real *sy, integer *incy);
/** @brief Library VBLAS prototypes.
    @author Michael Holst  */
VEXTERNC E_f sdot_(integer *n, real *sx, integer *incx, real *sy, integer *incy);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124
    @author Michael Holst  */
VEXTERNC int sgbmv_(char *trans, integer *m, integer *n, integer *kl, integer *ku, real *alpha, real *a, integer *lda, real *x, integer *incx, real *beta, real *y, integer *incy, ftnlen trans_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124
    @author Michael Holst  */
VEXTERNC int sgemm_(char *transa, char *transb, integer *m, integer *n, integer *k, real *alpha, real *a, integer *lda, real *b, integer *ldb, real *beta, real *c__, integer *ldc, ftnlen transa_len, ftnlen transb_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124
    @author Michael Holst  */
VEXTERNC int sgemv_(char *trans, integer *m, integer *n, real *alpha, real *a, integer *lda, real *x, integer *incx, real *beta, real *y, integer *incy, ftnlen trans_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: xerbla_ 14 3 13 4 124
    @author Michael Holst  */
VEXTERNC int sger_(integer *m, integer *n, real *alpha, real *x, integer *incx, real *y, integer *incy, real *a, integer *lda);
/** @brief Library VBLAS prototypes.
    @author Michael Holst  */
VEXTERNC E_f snrm2_(integer *n, real *x, integer *incx);
/** @brief Library VBLAS prototypes.
    @author Michael Holst  */
VEXTERNC int srot_(integer *n, real *sx, integer *incx, real *sy, integer *incy, real *c__, real *s);
/** @brief Library VBLAS prototypes.
    @author Michael Holst  */
VEXTERNC int srotg_(real *sa, real *sb, real *c__, real *s);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124
    @author Michael Holst  */
VEXTERNC int ssbmv_(char *uplo, integer *n, integer *k, real *alpha, real *a, integer *lda, real *x, integer *incx, real *beta, real *y, integer *incy, ftnlen uplo_len);
/** @brief Library VBLAS prototypes.
    @author Michael Holst  */
VEXTERNC int sscal_(integer *n, real *sa, real *sx, integer *incx);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124
    @author Michael Holst  */
VEXTERNC int sspmv_(char *uplo, integer *n, real *alpha, real *ap, real *x, integer *incx, real *beta, real *y, integer *incy, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124
    @author Michael Holst  */
VEXTERNC int sspr_(char *uplo, integer *n, real *alpha, real *x, integer *incx, real *ap, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124
    @author Michael Holst  */
VEXTERNC int sspr2_(char *uplo, integer *n, real *alpha, real *x, integer *incx, real *y, integer *incy, real *ap, ftnlen uplo_len);
/** @brief Library VBLAS prototypes.
    @author Michael Holst  */
VEXTERNC int sswap_(integer *n, real *sx, integer *incx, real *sy, integer *incy);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int ssymm_(char *side, char *uplo, integer *m, integer *n, real *alpha, real *a, integer *lda, real *b, integer *ldb, real *beta, real *c__, integer *ldc, ftnlen side_len, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124     
    @author Michael Holst */
VEXTERNC int ssymv_(char *uplo, integer *n, real *alpha, real *a, integer *lda, real *x, integer *incx, real *beta, real *y, integer *incy, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124
    @author Michael Holst  */
VEXTERNC int ssyr_(char *uplo, integer *n, real *alpha, real *x, integer *incx, real *a, integer *lda, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int ssyr2_(char *uplo, integer *n, real *alpha, real *x, integer *incx, real *y, integer *incy, real *a, integer *lda, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int ssyr2k_(char *uplo, char *trans, integer *n, integer *k, real *alpha, real *a, integer *lda, real *b, integer *ldb, real *beta, real *c__, integer *ldc, ftnlen uplo_len, ftnlen trans_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int ssyrk_(char *uplo, char *trans, integer *n, integer *k, real *alpha, real *a, integer *lda, real *beta, real *c__, integer *ldc, ftnlen uplo_len, ftnlen trans_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int stbmv_(char *uplo, char *trans, char *diag, integer *n, integer *k, real *a, integer *lda, real *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int stbsv_(char *uplo, char *trans, char *diag, integer *n, integer *k, real *a, integer *lda, real *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int stpmv_(char *uplo, char *trans, char *diag, integer *n, real *ap, real *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int stpsv_(char *uplo, char *trans, char *diag, integer *n, real *ap, real *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int strmm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, real *alpha, real *a, integer *lda, real *b, integer *ldb, ftnlen side_len, ftnlen uplo_len, ftnlen transa_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int strmv_(char *uplo, char *trans, char *diag, integer *n, real *a, integer *lda, real *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int strsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, real *alpha, real *a, integer *lda, real *b, integer *ldb, ftnlen side_len, ftnlen uplo_len, ftnlen transa_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int strsv_(char *uplo, char *trans, char *diag, integer *n, real *a, integer *lda, real *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: dcabs1_ 7 1 9 
    @author Michael Holst */
VEXTERNC int zaxpy_(integer *n, doublecomplex *za, doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy);
/** @brief Library VBLAS prototypes. 
    @author Michael Holst */
VEXTERNC int zcopy_(integer *n, doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy);
/** @brief Library VBLAS prototypes. 
    @author Michael Holst */
VEXTERNC Z_f zdotc_(doublecomplex * ret_val, integer *n, doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy);
/** @brief Library VBLAS prototypes. 
    @author Michael Holst */
VEXTERNC Z_f zdotu_(doublecomplex * ret_val, integer *n, doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy);
/** @brief Library VBLAS prototypes. 
    @author Michael Holst */
VEXTERNC int zdscal_(integer *n, doublereal *da, doublecomplex *zx, integer *incx);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int zgbmv_(char *trans, integer *m, integer *n, integer *kl, integer *ku, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, doublecomplex *beta, doublecomplex *y, integer *incy, ftnlen trans_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int zgemm_(char *transa, char *transb, integer *m, integer *n, integer *k, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublecomplex *beta, doublecomplex *c__, integer *ldc, ftnlen transa_len, ftnlen transb_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int zgemv_(char *trans, integer *m, integer *n, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, doublecomplex *beta, doublecomplex *y, integer *incy, ftnlen trans_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int zgerc_(integer *m, integer *n, doublecomplex *alpha, doublecomplex *x, integer *incx, doublecomplex *y, integer *incy, doublecomplex *a, integer *lda);
/** @brief Library VBLAS prototypes. \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int zgeru_(integer *m, integer *n, doublecomplex *alpha, doublecomplex *x, integer *incx, doublecomplex *y, integer *incy, doublecomplex *a, integer *lda);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int zhbmv_(char *uplo, integer *n, integer *k, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, doublecomplex *beta, doublecomplex *y, integer *incy, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int zhemm_(char *side, char *uplo, integer *m, integer *n, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublecomplex *beta, doublecomplex *c__, integer *ldc, ftnlen side_len, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int zhemv_(char *uplo, integer *n, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, doublecomplex *beta, doublecomplex *y, integer *incy, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int zher_(char *uplo, integer *n, doublereal *alpha, doublecomplex *x, integer *incx, doublecomplex *a, integer *lda, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int zher2_(char *uplo, integer *n, doublecomplex *alpha, doublecomplex *x, integer *incx, doublecomplex *y, integer *incy, doublecomplex *a, integer *lda, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int zher2k_(char *uplo, char *trans, integer *n, integer *k, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublereal *beta, doublecomplex *c__, integer *ldc, ftnlen uplo_len, ftnlen trans_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int zherk_(char *uplo, char *trans, integer *n, integer *k, doublereal *alpha, doublecomplex *a, integer *lda, doublereal *beta, doublecomplex *c__, integer *ldc, ftnlen uplo_len, ftnlen trans_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int zhpmv_(char *uplo, integer *n, doublecomplex *alpha, doublecomplex *ap, doublecomplex *x, integer *incx, doublecomplex *beta, doublecomplex *y, integer *incy, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int zhpr_(char *uplo, integer *n, doublereal *alpha, doublecomplex *x, integer *incx, doublecomplex *ap, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int zhpr2_(char *uplo, integer *n, doublecomplex *alpha, doublecomplex *x, integer *incx, doublecomplex *y, integer *incy, doublecomplex *ap, ftnlen uplo_len);
/** @brief Library VBLAS prototypes.
    @author Michael Holst */
VEXTERNC int zrotg_(doublecomplex *ca, doublecomplex *cb, doublereal *c__, doublecomplex *s);
/** @brief Library VBLAS prototypes. 
    @author Michael Holst */
VEXTERNC int zscal_(integer *n, doublecomplex *za, doublecomplex *zx, integer *incx);
/** @brief Library VBLAS prototypes. 
    @author Michael Holst */
VEXTERNC int zswap_(integer *n, doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int zsymm_(char *side, char *uplo, integer *m, integer *n, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublecomplex *beta, doublecomplex *c__, integer *ldc, ftnlen side_len, ftnlen uplo_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int zsyr2k_(char *uplo, char *trans, integer *n, integer *k, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublecomplex *beta, doublecomplex *c__, integer *ldc, ftnlen uplo_len, ftnlen trans_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int zsyrk_(char *uplo, char *trans, integer *n, integer *k, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *beta, doublecomplex *c__, integer *ldc, ftnlen uplo_len, ftnlen trans_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int ztbmv_(char *uplo, char *trans, char *diag, integer *n, integer *k, doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int ztbsv_(char *uplo, char *trans, char *diag, integer *n, integer *k, doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int ztpmv_(char *uplo, char *trans, char *diag, integer *n, doublecomplex *ap, doublecomplex *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int ztpsv_(char *uplo, char *trans, char *diag, integer *n, doublecomplex *ap, doublecomplex *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int ztrmm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, ftnlen side_len, ftnlen uplo_len, ftnlen transa_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int ztrmv_(char *uplo, char *trans, char *diag, integer *n, doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int ztrsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, ftnlen side_len, ftnlen uplo_len, ftnlen transa_len, ftnlen diag_len);
/** @brief Library VBLAS prototypes. \n
 *:ref: lsame_ 12 4 13 13 124 124 \n
 *:ref: xerbla_ 14 3 13 4 124 
    @author Michael Holst */
VEXTERNC int ztrsv_(char *uplo, char *trans, char *diag, integer *n, doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);


#endif /* _VBLAS_H_ */

