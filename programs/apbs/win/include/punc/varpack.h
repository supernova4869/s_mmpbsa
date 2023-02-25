/**
 *  @file       varpack.h
 *  @brief      The primary header ARPACK.
 *              (Arnoldi Package for solving sparse eigenproblems.)
 *  @note       We provide this header whether or not we provide ARPACK
 *              library itself.  This gives us some compile-time type-checking
 *              even for architecture-dependent assembly-coded ARPACK.
 *  @version    $Id: varpack.h,v 1.11 2010/08/12 05:52:31 fetk Exp $
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

#ifndef _VARPACK_H_
#define _VARPACK_H_

#include <punc/punc_base.h>

#include <punc/vf2c.h>

/*
 * ***************************************************************************
 * Library VARPACK prototypes
 * ***************************************************************************
 */

/** @brief  Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note 
 *  comlen debug_ 96\n
 *  comlen timing_ 124\n
 *  :ref: second_ 14 1 6\n 
 *  :ref: clarnv_ 14 4 4 4 4 8\n
 *  :ref: ccopy_ 14 5 4 8 4 8 4\n
 *  :ref: cdotc_ 8 6 8 4 8 4 8 4\n
 *  :ref: slapy2_ 6 2 6 6\n
 *  :ref: scnrm2_ 6 3 4 8 4\n
 *  :ref: cgemv_ 14 12 13 4 4 8 8 4 8 4 8 8 4 124\n
 *  :ref: svout_ 14 6 4 4 6 4 13 124\n
 *  :ref: cvout_ 14 6 4 4 8 4 13 124 */
VEXTERNC int cgetv0_(integer *ido, char *bmat, integer *itry, logical *initv, integer *n, integer *j, realcomplex *v, integer *ldv, realcomplex *resid, real *rnorm, integer *ipntr, realcomplex *workd, integer *ierr, ftnlen bmat_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note  
 *  comlen debug_ 96 \n
 *  comlen timing_ 124 \n
 *  :ref: slamch_ 6 2 13 124 \n
 *  :ref: slabad_ 14 2 6 6 \n
 *  :ref: second_ 14 1 6 \n
 *  :ref: ivout_ 14 6 4 4 4 4 13 124 \n
 *  :ref: svout_ 14 6 4 4 6 4 13 124 \n
 *  :ref: cgetv0_ 14 14 4 13 4 12 4 4 8 4 8 6 4 8 4 124 \n
 *  :ref: ccopy_ 14 5 4 8 4 8 4 \n
 *  :ref: csscal_ 14 4 4 6 8 4 \n
 *  :ref: clascl_ 14 11 13 4 4 6 6 4 4 8 4 4 124 \n
 *  :ref: cdotc_ 8 6 8 4 8 4 8 4 \n
 *  :ref: slapy2_ 6 2 6 6 \n
 *  :ref: scnrm2_ 6 3 4 8 4 \n
 *  :ref: cgemv_ 14 12 13 4 4 8 8 4 8 4 8 8 4 124 \n
 *  :ref: cvout_ 14 6 4 4 8 4 13 124 \n
 *  :ref: caxpy_ 14 6 4 8 8 4 8 4 \n
 *  :ref: clanhs_ 6 6 13 4 8 4 8 124 \n
 *  :ref: cmout_ 14 8 4 4 4 8 4 4 13 124 */
VEXTERNC int cnaitr_(integer *ido, char *bmat, integer *n, integer *k, integer *np, integer *nb, realcomplex *resid, real *rnorm, realcomplex *v, integer *ldv, realcomplex *h__, integer *ldh, integer *ipntr, realcomplex *workd, integer *info, ftnlen bmat_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  comlen debug_ 96 \n
 *  comlen timing_ 124 \n
 *  :ref: slamch_ 6 2 13 124 \n
 *  :ref: slabad_ 14 2 6 6 \n
 *  :ref: second_ 14 1 6 \n
 *  :ref: claset_ 14 8 13 4 4 8 8 8 4 124 \n
 *  :ref: ivout_ 14 6 4 4 4 4 13 124 \n
 *  :ref: cvout_ 14 6 4 4 8 4 13 124 \n
 *  :ref: clanhs_ 6 6 13 4 8 4 8 124 \n
 *  :ref: clartg_ 14 5 8 8 6 8 8 \n
 *  :ref: slapy2_ 6 2 6 6 \n
 *  :ref: cscal_ 14 4 4 8 8 4 \n
 *  :ref: cgemv_ 14 12 13 4 4 8 8 4 8 4 8 8 4 124 \n
 *  :ref: ccopy_ 14 5 4 8 4 8 4 \n
 *  :ref: clacpy_ 14 8 13 4 4 8 4 8 4 124 \n
 *  :ref: caxpy_ 14 6 4 8 8 4 8 4 \n
 *  :ref: cmout_ 14 8 4 4 4 8 4 4 13 124 \n
 */
VEXTERNC int cnapps_(integer *n, integer *kev, integer *np, realcomplex *shift, realcomplex *v, integer *ldv, realcomplex *h__, integer *ldh, realcomplex *resid, realcomplex *q, integer *ldq, realcomplex *workl, realcomplex *workd);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  comlen debug_ 96 \n
 *  comlen timing_ 124 \n
 *  :ref: second_ 14 1 6 \n
 *  :ref: slamch_ 6 2 13 124 \n
 *  :ref: cgetv0_ 14 14 4 13 4 12 4 4 8 4 8 6 4 8 4 124 \n
 *  :ref: cnaitr_ 14 16 4 13 4 4 4 4 8 6 8 4 8 4 4 8 4 124 \n
 *  :ref: ivout_ 14 6 4 4 4 4 13 124 \n
 *  :ref: svout_ 14 6 4 4 6 4 13 124 \n
 *  :ref: cneigh_ 14 11 6 4 8 4 8 8 8 4 8 6 4 \n
 *  :ref: ccopy_ 14 5 4 8 4 8 4 \n
 *  :ref: cngets_ 14 7 4 13 4 4 8 8 124 \n
 *  :ref: slapy2_ 6 2 6 6 \n
 *  :ref: cvout_ 14 6 4 4 8 4 13 124 \n
 *  :ref: csortc_ 14 6 13 12 4 8 8 124 \n
 *  :ref: cnapps_ 14 13 4 4 4 8 8 4 8 4 8 8 4 8 8 \n
 *  :ref: cdotc_ 8 6 8 4 8 4 8 4 \n
 *  :ref: scnrm2_ 6 3 4 8 4 \n
 *  :ref: cmout_ 14 8 4 4 4 8 4 4 13 124 
 */
VEXTERNC int cnaup2_(integer *ido, char *bmat, integer *n, char *which, integer *nev, integer *np, real *tol, realcomplex *resid, integer *mode, integer *iupd, integer *ishift, integer *mxiter, realcomplex *v, integer *ldv, realcomplex *h__, integer *ldh, realcomplex *ritz, realcomplex *bounds, realcomplex *q, integer *ldq, realcomplex *workl, integer *ipntr, realcomplex *workd, real *rwork, integer *info, ftnlen bmat_len, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  comlen debug_ 96 \n
 *  comlen timing_ 124 \n
 *  :ref: cstatn_ 14 0 \n
 *  :ref: second_ 14 1 6 \n
 *  :ref: slamch_ 6 2 13 124 \n
 *  :ref: cnaup2_ 14 27 4 13 4 13 4 4 6 8 4 4 4 4 8 4 8 4 8 8 8 4 8 4 8 6 4 124 124 \n
 *  :ref: ivout_ 14 6 4 4 4 4 13 124 \n
 *  :ref: cvout_ 14 6 4 4 8 4 13 124 
 */
VEXTERNC int cnaupd_(integer *ido, char *bmat, integer *n, char *which, integer *nev, real *tol, realcomplex *resid, integer *ncv, realcomplex *v, integer *ldv, integer *iparam, integer *ipntr, realcomplex *workd, realcomplex *workl, integer *lworkl, real *rwork, integer *info, ftnlen bmat_len, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  comlen debug_ 96 \n
 *  comlen timing_ 124 \n
 *  :ref: second_ 14 1 6 \n
 *  :ref: cmout_ 14 8 4 4 4 8 4 4 13 124 \n
 *  :ref: clacpy_ 14 8 13 4 4 8 4 8 4 124 \n
 *  :ref: claset_ 14 8 13 4 4 8 8 8 4 124 \n
 *  :ref: clahqr_ 14 13 12 12 4 4 4 8 4 8 4 4 8 4 4 \n
 *  :ref: ccopy_ 14 5 4 8 4 8 4 \n
 *  :ref: cvout_ 14 6 4 4 8 4 13 124 \n
 *  :ref: ctrevc_ 14 17 13 13 12 4 8 4 8 4 8 4 4 4 8 6 4 124 124 \n
 *  :ref: scnrm2_ 6 3 4 8 4 \n
 *  :ref: csscal_ 14 4 4 6 8 4 \n
 */
VEXTERNC int cneigh_(real *rnorm, integer *n, realcomplex *h__, integer *ldh, realcomplex *ritz, realcomplex *bounds, realcomplex *q, integer *ldq, realcomplex *workl, real *rwork, integer *ierr);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  comlen debug_ 96 \n
 *  comlen timing_ 124 \n
 *  :ref: slamch_ 6 2 13 124 \n
 *  :ref: cvout_ 14 6 4 4 8 4 13 124 \n
 *  :ref: cngets_ 14 7 4 13 4 4 8 8 124 \n 
 *  :ref: slapy2_ 6 2 6 6 \n
 *  :ref: ivout_ 14 6 4 4 4 4 13 124 \n
 *  :ref: ccopy_ 14 5 4 8 4 8 4 \n
 *  :ref: claset_ 14 8 13 4 4 8 8 8 4 124 \n
 *  :ref: clahqr_ 14 13 12 12 4 4 4 8 4 8 4 4 8 4 4 \n
 *  :ref: cmout_ 14 8 4 4 4 8 4 4 13 124 \n
 *  :ref: ctrsen_ 14 17 13 13 12 4 8 4 8 4 8 4 6 6 8 4 4 124 124 \n
 *  :ref: cgeqr2_ 14 7 4 4 8 4 8 8 4 \n
 *  :ref: cunm2r_ 14 14 13 13 4 4 4 8 4 8 8 4 8 4 124 124 \n
 *  :ref: clacpy_ 14 8 13 4 4 8 4 8 4 124 \n
 *  :ref: cscal_ 14 4 4 8 8 4 \n
 *  :ref: ctrevc_ 14 17 13 13 12 4 8 4 8 4 8 4 4 4 8 6 4 124 124 \n
 *  :ref: scnrm2_ 6 3 4 8 4 \n
 *  :ref: csscal_ 14 4 4 6 8 4 \n
 *  :ref: cdotc_ 8 6 8 4 8 4 8 4 \n
 *  :ref: ctrmm_ 14 15 13 13 13 13 4 4 8 8 4 8 4 124 124 124 124 \n
 *  :ref: cgeru_ 14 9 4 4 8 8 4 8 4 8 4 
 */
VEXTERNC int cneupd_(logical *rvec, char *howmny, logical *select, realcomplex *d__, realcomplex *z__, integer *ldz, realcomplex *sigma, realcomplex *workev, char *bmat, integer *n, char *which, integer *nev, real *tol, realcomplex *resid, integer *ncv, realcomplex *v, integer *ldv, integer *iparam, integer *ipntr, realcomplex *workd, realcomplex *workl, integer *lworkl, real *rwork, integer *info, ftnlen howmny_len, ftnlen bmat_len, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  comlen debug_ 96 \n
 *  comlen timing_ 124 \n
 *  :ref: second_ 14 1 6 \n
 *  :ref: csortc_ 14 6 13 12 4 8 8 124 \n
 *  :ref: ivout_ 14 6 4 4 4 4 13 124 \n
 *  :ref: cvout_ 14 6 4 4 8 4 13 124 
 */
VEXTERNC int cngets_(integer *ishift, char *which, integer *kev, integer *np, realcomplex *ritz, realcomplex *bounds, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  :ref: slapy2_ 6 2 6 6 
 */
VEXTERNC int csortc_(char *which, logical *apply, integer *n, realcomplex *x, realcomplex *y, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  comlen timing_ 124 */
VEXTERNC int cstatn_(void);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  comlen debug_ 96 \n
 *  comlen timing_ 124 \n  
 *  :ref: second_ 14 1 6 \n  
 *  :ref: dlarnv_ 14 4 4 4 4 7 \n  
 *  :ref: dcopy_ 14 5 4 7 4 7 4 \n  
 *  :ref: ddot_ 7 5 4 7 4 7 4 \n  
 *  :ref: dnrm2_ 7 3 4 7 4 \n  
 *  :ref: dgemv_ 14 12 13 4 4 7 7 4 7 4 7 7 4 124 \n  
 *  :ref: dvout_ 14 6 4 4 7 4 13 124 
 */
VEXTERNC int dgetv0_(integer *ido, char *bmat, integer *itry, logical *initv, integer *n, integer *j, doublereal *v, integer *ldv, doublereal *resid, doublereal *rnorm, integer *ipntr, doublereal *workd, integer *ierr, ftnlen bmat_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  :ref: dlamch_ 7 2 13 124 \n
 *  :ref: dlabad_ 14 2 7 7 \n
 *  :ref: dlanhs_ 7 6 13 4 7 4 7 124 \n
 *  :ref: dcopy_ 14 5 4 7 4 7 4 \n
 *  :ref: dlarfg_ 14 5 4 7 7 4 7 \n
 *  :ref: dlanv2_ 14 10 7 7 7 7 7 7 7 7 7 7 \n
 *  :ref: drot_ 14 7 4 7 4 7 4 7 7 
 */
VEXTERNC int dlaqrb_(logical *wantt, integer *n, integer *ilo, integer *ihi, doublereal *h__, integer *ldh, doublereal *wr, doublereal *wi, doublereal *z__, integer *info);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  comlen debug_ 96 \n
 *  comlen timing_ 124 \n
 *  :ref: dlamch_ 7 2 13 124 \n
 *  :ref: dlabad_ 14 2 7 7 \n
 *  :ref: second_ 14 1 6 \n
 *  :ref: ivout_ 14 6 4 4 4 4 13 124 \n
 *  :ref: dvout_ 14 6 4 4 7 4 13 124 \n
 *  :ref: dgetv0_ 14 14 4 13 4 12 4 4 7 4 7 7 4 7 4 124 \n
 *  :ref: dcopy_ 14 5 4 7 4 7 4 \n
 *  :ref: dscal_ 14 4 4 7 7 4 \n
 *  :ref: dlascl_ 14 11 13 4 4 7 7 4 4 7 4 4 124 \n
 *  :ref: ddot_ 7 5 4 7 4 7 4 \n
 *  :ref: dnrm2_ 7 3 4 7 4 \n
 *  :ref: dgemv_ 14 12 13 4 4 7 7 4 7 4 7 7 4 124 \n
 *  :ref: daxpy_ 14 6 4 7 7 4 7 4 \n
 *  :ref: dlanhs_ 7 6 13 4 7 4 7 124 \n
 *  :ref: dmout_ 14 8 4 4 4 7 4 4 13 124 
 */
VEXTERNC int dnaitr_(integer *ido, char *bmat, integer *n, integer *k, integer *np, integer *nb, doublereal *resid, doublereal *rnorm, doublereal *v, integer *ldv, doublereal *h__, integer *ldh, integer *ipntr, doublereal *workd, integer *info, ftnlen bmat_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  comlen debug_ 96 \n
 *  comlen timing_ 124 \n 
 *  :ref: dlamch_ 7 2 13 124 \n 
 *  :ref: dlabad_ 14 2 7 7 \n 
 *  :ref: second_ 14 1 6 \n 
 *  :ref: dlaset_ 14 8 13 4 4 7 7 7 4 124 \n 
 *  :ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *  :ref: dvout_ 14 6 4 4 7 4 13 124 \n 
 *  :ref: dlanhs_ 7 6 13 4 7 4 7 124 \n 
 *  :ref: dlartg_ 14 5 7 7 7 7 7 \n 
 *  :ref: dlapy2_ 7 2 7 7 \n 
 *  :ref: dlarfg_ 14 5 4 7 7 4 7 \n 
 *  :ref: dlarf_ 14 10 13 4 4 7 4 7 7 4 7 124 \n 
 *  :ref: dscal_ 14 4 4 7 7 4 \n 
 *  :ref: dgemv_ 14 12 13 4 4 7 7 4 7 4 7 7 4 124 \n 
 *  :ref: dcopy_ 14 5 4 7 4 7 4 \n 
 *  :ref: dlacpy_ 14 8 13 4 4 7 4 7 4 124 \n 
 *  :ref: daxpy_ 14 6 4 7 7 4 7 4 \n 
 *  :ref: dmout_ 14 8 4 4 4 7 4 4 13 124 
 */
VEXTERNC int dnapps_(integer *n, integer *kev, integer *np, doublereal *shiftr, doublereal *shifti, doublereal *v, integer *ldv, doublereal *h__, integer *ldh, doublereal *resid, doublereal *q, integer *ldq, doublereal *workl, doublereal *workd);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  comlen debug_ 96 \n 
 *  comlen timing_ 124 \n 
 *  :ref: second_ 14 1 6 \n 
 *  :ref: dlamch_ 7 2 13 124 \n 
 *  :ref: dgetv0_ 14 14 4 13 4 12 4 4 7 4 7 7 4 7 4 124 \n 
 *  :ref: dnaitr_ 14 16 4 13 4 4 4 4 7 7 7 4 7 4 4 7 4 124 \n 
 *  :ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *  :ref: dvout_ 14 6 4 4 7 4 13 124 \n 
 *  :ref: dneigh_ 14 11 7 4 7 4 7 7 7 7 4 7 4 \n 
 *  :ref: dcopy_ 14 5 4 7 4 7 4 \n 
 *  :ref: dngets_ 14 10 4 13 4 4 7 7 7 7 7 124 \n 
 *  :ref: dnconv_ 14 6 4 7 7 7 7 4 \n 
 *  :ref: dsortc_ 14 7 13 12 4 7 7 7 124 \n 
 *  :ref: dlapy2_ 7 2 7 7 \n 
 *  :ref: dnapps_ 14 14 4 4 4 7 7 7 4 7 4 7 7 4 7 7 \n 
 *  :ref: ddot_ 7 5 4 7 4 7 4 \n 
 *  :ref: dnrm2_ 7 3 4 7 4 \n 
 *  :ref: dmout_ 14 8 4 4 4 7 4 4 13 124 \n 
 */
VEXTERNC int dnaup2_(integer *ido, char *bmat, integer *n, char *which, integer *nev, integer *np, doublereal *tol, doublereal *resid, integer *mode, integer *iupd, integer *ishift, integer *mxiter, doublereal *v, integer *ldv, doublereal *h__, integer *ldh, doublereal *ritzr, doublereal *ritzi, doublereal *bounds, doublereal *q, integer *ldq, doublereal *workl, integer *ipntr, doublereal *workd, integer *info, ftnlen bmat_len, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  comlen debug_ 96 \n 
 *  comlen timing_ 124 \n 
 *  :ref: dstatn_ 14 0 \n 
 *  :ref: second_ 14 1 6 \n 
 *  :ref: dlamch_ 7 2 13 124 \n 
 *  :ref: dnaup2_ 14 27 4 13 4 13 4 4 7 7 4 4 4 4 7 4 7 4 7 7 7 7 4 7 4 7 4 124 124 \n 
 *  :ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *  :ref: dvout_ 14 6 4 4 7 4 13 124
 */
VEXTERNC int dnaupd_(integer *ido, char *bmat, integer *n, char *which, integer *nev, doublereal *tol, doublereal *resid, integer *ncv, doublereal *v, integer *ldv, integer *iparam, integer *ipntr, doublereal *workd, doublereal *workl, integer *lworkl, integer *info, ftnlen bmat_len, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  comlen debug_ 96 \n 
 *  comlen timing_ 124 \n 
 *  :ref: second_ 14 1 6 \n 
 *  :ref: dlamch_ 7 2 13 124 \n 
 *  :ref: dlapy2_ 7 2 7 7 
 */
VEXTERNC int dnconv_(integer *n, doublereal *ritzr, doublereal *ritzi, doublereal *bounds, doublereal *tol, integer *nconv);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  comlen debug_ 96 \n 
 *  comlen timing_ 124 \n 
 *  :ref: second_ 14 1 6 \n 
 *  :ref: dmout_ 14 8 4 4 4 7 4 4 13 124 \n 
 *  :ref: dlacpy_ 14 8 13 4 4 7 4 7 4 124 \n 
 *  :ref: dlaqrb_ 14 10 12 4 4 4 7 4 7 7 7 4 \n 
 *  :ref: dvout_ 14 6 4 4 7 4 13 124 \n 
 *  :ref: dtrevc_ 14 16 13 13 12 4 7 4 7 4 7 4 4 4 7 4 124 124 \n 
 *  :ref: dnrm2_ 7 3 4 7 4 \n 
 *  :ref: dscal_ 14 4 4 7 7 4 \n 
 *  :ref: dlapy2_ 7 2 7 7 \n 
 *  :ref: dgemv_ 14 12 13 4 4 7 7 4 7 4 7 7 4 124 
 */
VEXTERNC int dneigh_(doublereal *rnorm, integer *n, doublereal *h__, integer *ldh, doublereal *ritzr, doublereal *ritzi, doublereal *bounds, doublereal *q, integer *ldq, doublereal *workl, integer *ierr);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  comlen debug_ 96 \n 
 *  comlen timing_ 124 \n 
 *  :ref: dlamch_ 7 2 13 124 \n 
 *  :ref: dvout_ 14 6 4 4 7 4 13 124 \n 
 *  :ref: dngets_ 14 10 4 13 4 4 7 7 7 7 7 124 \n
 *  :ref: dlapy2_ 7 2 7 7 \n 
 *  :ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *  :ref: dcopy_ 14 5 4 7 4 7 4 \n 
 *  :ref: dlaset_ 14 8 13 4 4 7 7 7 4 124 \n 
 *  :ref: dlahqr_ 14 14 12 12 4 4 4 7 4 7 7 4 4 7 4 4 \n 
 *  :ref: dmout_ 14 8 4 4 4 7 4 4 13 124 \n 
 *  :ref: dtrsen_ 14 20 13 13 12 4 7 4 7 4 7 7 4 7 7 7 4 4 4 4 124 124 \n 
 *  :ref: dgeqr2_ 14 7 4 4 7 4 7 7 4 \n 
 *  :ref: dorm2r_ 14 14 13 13 4 4 4 7 4 7 7 4 7 4 124 124 \n 
 *  :ref: dlacpy_ 14 8 13 4 4 7 4 7 4 124 \n 
 *  :ref: dscal_ 14 4 4 7 7 4 \n 
 *  :ref: dtrevc_ 14 16 13 13 12 4 7 4 7 4 7 4 4 4 7 4 124 124 \n 
 *  :ref: dnrm2_ 7 3 4 7 4 \n 
 *  :ref: dgemv_ 14 12 13 4 4 7 7 4 7 4 7 7 4 124 \n 
 *  :ref: dtrmm_ 14 15 13 13 13 13 4 4 7 7 4 7 4 124 124 124 124 \n 
 *  :ref: dger_ 14 9 4 4 7 7 4 7 4 7 4 
 */
VEXTERNC int dneupd_(logical *rvec, char *howmny, logical *select, doublereal *dr, doublereal *di, doublereal *z__, integer *ldz, doublereal *sigmar, doublereal *sigmai, doublereal *workev, char *bmat, integer *n, char *which, integer *nev, doublereal *tol, doublereal *resid, integer *ncv, doublereal *v, integer *ldv, integer *iparam, integer *ipntr, doublereal *workd, doublereal *workl, integer *lworkl, integer *info, ftnlen howmny_len, ftnlen bmat_len, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  comlen debug_ 96 \n 
 *  comlen timing_ 124 \n 
 *  :ref: second_ 14 1 6 \n 
 *  :ref: dsortc_ 14 7 13 12 4 7 7 7 124 \n 
 *  :ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *  :ref: dvout_ 14 6 4 4 7 4 13 124
 */
VEXTERNC int dngets_(integer *ishift, char *which, integer *kev, integer *np, doublereal *ritzr, doublereal *ritzi, doublereal *bounds, doublereal *shiftr, doublereal *shifti, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  comlen debug_ 96 \n 
 *  comlen timing_ 124 \n 
 *  :ref: dlamch_ 7 2 13 124 \n 
 *  :ref: second_ 14 1 6 \n 
 *  :ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *  :ref: dvout_ 14 6 4 4 7 4 13 124 \n 
 *  :ref: dgetv0_ 14 14 4 13 4 12 4 4 7 4 7 7 4 7 4 124 \n 
 *  :ref: dcopy_ 14 5 4 7 4 7 4 \n 
 *  :ref: dscal_ 14 4 4 7 7 4 \n 
 *  :ref: dlascl_ 14 11 13 4 4 7 7 4 4 7 4 4 124 \n 
 *  :ref: ddot_ 7 5 4 7 4 7 4 \n 
 *  :ref: dnrm2_ 7 3 4 7 4 \n 
 *  :ref: dgemv_ 14 12 13 4 4 7 7 4 7 4 7 7 4 124 
 */
VEXTERNC int dsaitr_(integer *ido, char *bmat, integer *n, integer *k, integer *np, integer *mode, doublereal *resid, doublereal *rnorm, doublereal *v, integer *ldv, doublereal *h__, integer *ldh, integer *ipntr, doublereal *workd, integer *info, ftnlen bmat_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  comlen debug_ 96 \n 
 *  comlen timing_ 124 \n 
 *  :ref: dlamch_ 7 2 13 124 \n
 *  :ref: second_ 14 1 6 \n 
 *  :ref: dlaset_ 14 8 13 4 4 7 7 7 4 124 \n 
 *  :ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *  :ref: dvout_ 14 6 4 4 7 4 13 124 \n 
 *  :ref: dlartg_ 14 5 7 7 7 7 7 \n 
 *  :ref: dscal_ 14 4 4 7 7 4 \n 
 *  :ref: dgemv_ 14 12 13 4 4 7 7 4 7 4 7 7 4 124 \n 
 *  :ref: dcopy_ 14 5 4 7 4 7 4 \n 
 *  :ref: dlacpy_ 14 8 13 4 4 7 4 7 4 124 \n 
 *  :ref: daxpy_ 14 6 4 7 7 4 7 4 
 */
VEXTERNC int dsapps_(integer *n, integer *kev, integer *np, doublereal *shift, doublereal *v, integer *ldv, doublereal *h__, integer *ldh, doublereal *resid, doublereal *q, integer *ldq, doublereal *workd);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  comlen debug_ 96 \n 
 *  comlen timing_ 124 \n 
 *  :ref: second_ 14 1 6 \n 
 *  :ref: dlamch_ 7 2 13 124 \n 
 *  :ref: dgetv0_ 14 14 4 13 4 12 4 4 7 4 7 7 4 7 4 124 \n 
 *  :ref: dsaitr_ 14 16 4 13 4 4 4 4 7 7 7 4 7 4 4 7 4 124 \n 
 *  :ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *  :ref: dvout_ 14 6 4 4 7 4 13 124 \n 
 *  :ref: dseigt_ 14 8 7 4 7 4 7 7 7 4 \n 
 *  :ref: dcopy_ 14 5 4 7 4 7 4 \n 
 *  :ref: dsgets_ 14 8 4 13 4 4 7 7 7 124 \n 
 *  :ref: dsconv_ 14 5 4 7 7 7 4 \n 
 *  :ref: dsortr_ 14 6 13 12 4 7 7 124 \n 
 *  :ref: dswap_ 14 5 4 7 4 7 4 \n 
 *  :ref: dsapps_ 14 12 4 4 4 7 7 4 7 4 7 7 4 7 \n 
 *  :ref: ddot_ 7 5 4 7 4 7 4 \n 
 *  :ref: dnrm2_ 7 3 4 7 4
 */
VEXTERNC int dsaup2_(integer *ido, char *bmat, integer *n, char *which, integer *nev, integer *np, doublereal *tol, doublereal *resid, integer *mode, integer *iupd, integer *ishift, integer *mxiter, doublereal *v, integer *ldv, doublereal *h__, integer *ldh, doublereal *ritz, doublereal *bounds, doublereal *q, integer *ldq, doublereal *workl, integer *ipntr, doublereal *workd, integer *info, ftnlen bmat_len, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  comlen debug_ 96 \n 
 *  comlen timing_ 124 \n 
 *  :ref: dstats_ 14 0 \n 
 *  :ref: second_ 14 1 6 \n 
 *  :ref: dlamch_ 7 2 13 124 \n 
 *  :ref: dsaup2_ 14 26 4 13 4 13 4 4 7 7 4 4 4 4 7 4 7 4 7 7 7 4 7 4 7 4 124 124 \n 
 *  :ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *  :ref: dvout_ 14 6 4 4 7 4 13 124 
 */
VEXTERNC int dsaupd_(integer *ido, char *bmat, integer *n, char *which, integer *nev, doublereal *tol, doublereal *resid, integer *ncv, doublereal *v, integer *ldv, integer *iparam, integer *ipntr, doublereal *workd, doublereal *workl, integer *lworkl, integer *info, ftnlen bmat_len, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  comlen debug_ 96 \n 
 *  comlen timing_ 124 \n 
 *  :ref: second_ 14 1 6 \n 
 *  :ref: dlamch_ 7 2 13 124 
 */
VEXTERNC int dsconv_(integer *n, doublereal *ritz, doublereal *bounds, doublereal *tol, integer *nconv);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  comlen debug_ 96 \n 
 *  comlen timing_ 124 \n 
 *  :ref: second_ 14 1 6 \n 
 *  :ref: dvout_ 14 6 4 4 7 4 13 124 \n 
 *  :ref: dcopy_ 14 5 4 7 4 7 4 \n 
 *  :ref: dstqrb_ 14 6 4 7 7 7 7 4 
 */
VEXTERNC int dseigt_(doublereal *rnorm, integer *n, doublereal *h__, integer *ldh, doublereal *eig, doublereal *bounds, doublereal *workl, integer *ierr);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  :ref: dswap_ 14 5 4 7 4 7 4 
 */
VEXTERNC int dsesrt_(char *which, logical *apply, integer *n, doublereal *x, integer *na, doublereal *a, integer *lda, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: dlamch_ 7 2 13 124 \n 
 *:ref: dnrm2_ 7 3 4 7 4 \n 
 *:ref: dvout_ 14 6 4 4 7 4 13 124 \n 
 *:ref: dsgets_ 14 8 4 13 4 4 7 7 7 124 \n 
 *:ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *:ref: dcopy_ 14 5 4 7 4 7 4 \n 
 *:ref: dsteqr_ 14 9 13 4 7 7 7 4 7 4 124 \n 
 *:ref: dsesrt_ 14 8 13 12 4 7 4 7 4 124 \n 
 *:ref: dsortr_ 14 6 13 12 4 7 7 124 \n 
 *:ref: dscal_ 14 4 4 7 7 4 \n 
 *:ref: dgeqr2_ 14 7 4 4 7 4 7 7 4 \n 
 *:ref: dorm2r_ 14 14 13 13 4 4 4 7 4 7 7 4 7 4 124 124 \n 
 *:ref: dlacpy_ 14 8 13 4 4 7 4 7 4 124 \n 
 *:ref: dger_ 14 9 4 4 7 7 4 7 4 7 4 
 */
VEXTERNC int dseupd_(logical *rvec, char *howmny, logical *select, doublereal *d__, doublereal *z__, integer *ldz, doublereal *sigma, char *bmat, integer *n, char *which, integer *nev, doublereal *tol, doublereal *resid, integer *ncv, doublereal *v, integer *ldv, integer *iparam, integer *ipntr, doublereal *workd, doublereal *workl, integer *lworkl, integer *info, ftnlen howmny_len, ftnlen bmat_len, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: second_ 14 1 6 \n 
 *:ref: dsortr_ 14 6 13 12 4 7 7 124 \n 
 *:ref: dswap_ 14 5 4 7 4 7 4 \n 
 *:ref: dcopy_ 14 5 4 7 4 7 4 \n 
 *:ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *:ref: dvout_ 14 6 4 4 7 4 13 124
 */
VEXTERNC int dsgets_(integer *ishift, char *which, integer *kev, integer *np, doublereal *ritz, doublereal *bounds, doublereal *shifts, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *:ref: dlapy2_ 7 2 7 7 
 */
VEXTERNC int dsortc_(char *which, logical *apply, integer *n, doublereal *xreal, doublereal *ximag, doublereal *y, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst 
 */
VEXTERNC int dsortr_(char *which, logical *apply, integer *n, doublereal *x1, doublereal *x2, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen timing_ 124 */
VEXTERNC int dstatn_(void);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen timing_ 124 */
VEXTERNC int dstats_(void);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *:ref: dlamch_ 7 2 13 124 \n 
 *:ref: dlanst_ 7 5 13 4 7 7 124 \n 
 *:ref: dlascl_ 14 11 13 4 4 7 7 4 4 7 4 4 124 \n 
 *:ref: dlaev2_ 14 7 7 7 7 7 7 7 7 \n 
 *:ref: dlae2_ 14 5 7 7 7 7 7 \n 
 *:ref: dlapy2_ 7 2 7 7 \n 
 *:ref: dlartg_ 14 5 7 7 7 7 7 \n 
 *:ref: dlasr_ 14 12 13 13 13 4 4 7 7 7 4 124 124 124 \n 
 *:ref: dlasrt_ 14 5 13 4 7 4 124 
 */
VEXTERNC int dstqrb_(integer *n, doublereal *d__, doublereal *e, doublereal *z__, doublereal *work, integer *info);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: second_ 14 1 6 \n 
 *:ref: slarnv_ 14 4 4 4 4 6 \n 
 *:ref: scopy_ 14 5 4 6 4 6 4 \n 
 *:ref: sdot_ 6 5 4 6 4 6 4 \n 
 *:ref: snrm2_ 6 3 4 6 4 \n 
 *:ref: sgemv_ 14 12 13 4 4 6 6 4 6 4 6 6 4 124 \n 
 *:ref: svout_ 14 6 4 4 6 4 13 124 
 */
VEXTERNC int sgetv0_(integer *ido, char *bmat, integer *itry, logical *initv, integer *n, integer *j, real *v, integer *ldv, real *resid, real *rnorm, integer *ipntr, real *workd, integer *ierr, ftnlen bmat_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *:ref: slamch_ 6 2 13 124 \n 
 *:ref: slabad_ 14 2 6 6 \n 
 *:ref: slanhs_ 6 6 13 4 6 4 6 124 \n 
 *:ref: scopy_ 14 5 4 6 4 6 4 \n 
 *:ref: slarfg_ 14 5 4 6 6 4 6 \n 
 *:ref: slanv2_ 14 10 6 6 6 6 6 6 6 6 6 6 \n 
 *:ref: srot_ 14 7 4 6 4 6 4 6 6 
 */
VEXTERNC int slaqrb_(logical *wantt, integer *n, integer *ilo, integer *ihi, real *h__, integer *ldh, real *wr, real *wi, real *z__, integer *info);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: slamch_ 6 2 13 124 \n 
 *:ref: slabad_ 14 2 6 6 \n 
 *:ref: second_ 14 1 6 \n 
 *:ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *:ref: svout_ 14 6 4 4 6 4 13 124 \n 
 *:ref: sgetv0_ 14 14 4 13 4 12 4 4 6 4 6 6 4 6 4 124 \n 
 *:ref: scopy_ 14 5 4 6 4 6 4 \n 
 *:ref: sscal_ 14 4 4 6 6 4 \n 
 *:ref: slascl_ 14 11 13 4 4 6 6 4 4 6 4 4 124 \n 
 *:ref: sdot_ 6 5 4 6 4 6 4 \n 
 *:ref: snrm2_ 6 3 4 6 4 \n 
 *:ref: sgemv_ 14 12 13 4 4 6 6 4 6 4 6 6 4 124 \n 
 *:ref: saxpy_ 14 6 4 6 6 4 6 4 \n 
 *:ref: slanhs_ 6 6 13 4 6 4 6 124 \n 
 *:ref: smout_ 14 8 4 4 4 6 4 4 13 124 
 */
VEXTERNC int snaitr_(integer *ido, char *bmat, integer *n, integer *k, integer *np, integer *nb, real *resid, real *rnorm, real *v, integer *ldv, real *h__, integer *ldh, integer *ipntr, real *workd, integer *info, ftnlen bmat_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: slamch_ 6 2 13 124 \n 
 *:ref: slabad_ 14 2 6 6 \n 
 *:ref: second_ 14 1 6 \n 
 *:ref: slaset_ 14 8 13 4 4 6 6 6 4 124 \n 
 *:ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *:ref: svout_ 14 6 4 4 6 4 13 124 \n 
 *:ref: slanhs_ 6 6 13 4 6 4 6 124 \n 
 *:ref: slartg_ 14 5 6 6 6 6 6 \n 
 *:ref: slapy2_ 6 2 6 6 \n 
 *:ref: slarfg_ 14 5 4 6 6 4 6 \n 
 *:ref: slarf_ 14 10 13 4 4 6 4 6 6 4 6 124 \n 
 *:ref: sscal_ 14 4 4 6 6 4 \n 
 *:ref: sgemv_ 14 12 13 4 4 6 6 4 6 4 6 6 4 124 \n 
 *:ref: scopy_ 14 5 4 6 4 6 4 \n 
 *:ref: slacpy_ 14 8 13 4 4 6 4 6 4 124 \n 
 *:ref: saxpy_ 14 6 4 6 6 4 6 4 \n 
 *:ref: smout_ 14 8 4 4 4 6 4 4 13 124 
 */
VEXTERNC int snapps_(integer *n, integer *kev, integer *np, real *shiftr, real *shifti, real *v, integer *ldv, real *h__, integer *ldh, real *resid, real *q, integer *ldq, real *workl, real *workd);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: second_ 14 1 6 \n 
 *:ref: slamch_ 6 2 13 124 \n 
 *:ref: sgetv0_ 14 14 4 13 4 12 4 4 6 4 6 6 4 6 4 124 \n 
 *:ref: snaitr_ 14 16 4 13 4 4 4 4 6 6 6 4 6 4 4 6 4 124 \n 
 *:ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *:ref: svout_ 14 6 4 4 6 4 13 124 \n 
 *:ref: sneigh_ 14 11 6 4 6 4 6 6 6 6 4 6 4 \n 
 *:ref: scopy_ 14 5 4 6 4 6 4 \n 
 *:ref: sngets_ 14 10 4 13 4 4 6 6 6 6 6 124 \n 
 *:ref: snconv_ 14 6 4 6 6 6 6 4 \n 
 *:ref: ssortc_ 14 7 13 12 4 6 6 6 124 \n 
 *:ref: slapy2_ 6 2 6 6 \n 
 *:ref: snapps_ 14 14 4 4 4 6 6 6 4 6 4 6 6 4 6 6 \n 
 *:ref: sdot_ 6 5 4 6 4 6 4 \n 
 *:ref: snrm2_ 6 3 4 6 4 \n 
 *:ref: smout_ 14 8 4 4 4 6 4 4 13 124 
 */
VEXTERNC int snaup2_(integer *ido, char *bmat, integer *n, char *which, integer *nev, integer *np, real *tol, real *resid, integer *mode, integer *iupd, integer *ishift, integer *mxiter, real *v, integer *ldv, real *h__, integer *ldh, real *ritzr, real *ritzi, real *bounds, real *q, integer *ldq, real *workl, integer *ipntr, real *workd, integer *info, ftnlen bmat_len, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: sstatn_ 14 0 \n 
 *:ref: second_ 14 1 6 \n 
 *:ref: slamch_ 6 2 13 124 \n 
 *:ref: snaup2_ 14 27 4 13 4 13 4 4 6 6 4 4 4 4 6 4 6 4 6 6 6 6 4 6 4 6 4 124 124 \n 
 *:ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *:ref: svout_ 14 6 4 4 6 4 13 124 
 */
VEXTERNC int snaupd_(integer *ido, char *bmat, integer *n, char *which, integer *nev, real *tol, real *resid, integer *ncv, real *v, integer *ldv, integer *iparam, integer *ipntr, real *workd, real *workl, integer *lworkl, integer *info, ftnlen bmat_len, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: second_ 14 1 6 \n 
 *:ref: slamch_ 6 2 13 124 \n 
 *:ref: slapy2_ 6 2 6 6 
 */
VEXTERNC int snconv_(integer *n, real *ritzr, real *ritzi, real *bounds, real *tol, integer *nconv);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: second_ 14 1 6 \n 
 *:ref: smout_ 14 8 4 4 4 6 4 4 13 124 \n 
 *:ref: slacpy_ 14 8 13 4 4 6 4 6 4 124 \n 
 *:ref: slaqrb_ 14 10 12 4 4 4 6 4 6 6 6 4 \n 
 *:ref: svout_ 14 6 4 4 6 4 13 124 \n 
 *:ref: strevc_ 14 16 13 13 12 4 6 4 6 4 6 4 4 4 6 4 124 124 \n 
 *:ref: snrm2_ 6 3 4 6 4 \n 
 *:ref: sscal_ 14 4 4 6 6 4 \n 
 *:ref: slapy2_ 6 2 6 6 \n 
 *:ref: sgemv_ 14 12 13 4 4 6 6 4 6 4 6 6 4 124 
 */
VEXTERNC int sneigh_(real *rnorm, integer *n, real *h__, integer *ldh, real *ritzr, real *ritzi, real *bounds, real *q, integer *ldq, real *workl, integer *ierr);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: slamch_ 6 2 13 124 \n 
 *:ref: svout_ 14 6 4 4 6 4 13 124 \n 
 *:ref: sngets_ 14 10 4 13 4 4 6 6 6 6 6 124 \n 
 *:ref: slapy2_ 6 2 6 6 \n 
 *:ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *:ref: scopy_ 14 5 4 6 4 6 4 \n 
 *:ref: slaset_ 14 8 13 4 4 6 6 6 4 124 \n 
 *:ref: slahqr_ 14 14 12 12 4 4 4 6 4 6 6 4 4 6 4 4 \n 
 *:ref: smout_ 14 8 4 4 4 6 4 4 13 124 \n 
 *:ref: strsen_ 14 20 13 13 12 4 6 4 6 4 6 6 4 6 6 6 4 4 4 4 124 124 \n 
 *:ref: sgeqr2_ 14 7 4 4 6 4 6 6 4 \n 
 *:ref: sorm2r_ 14 14 13 13 4 4 4 6 4 6 6 4 6 4 124 124 \n 
 *:ref: slacpy_ 14 8 13 4 4 6 4 6 4 124 \n 
 *:ref: sscal_ 14 4 4 6 6 4 \n 
 *:ref: strevc_ 14 16 13 13 12 4 6 4 6 4 6 4 4 4 6 4 124 124 \n 
 *:ref: snrm2_ 6 3 4 6 4 \n 
 *:ref: sgemv_ 14 12 13 4 4 6 6 4 6 4 6 6 4 124 \n 
 *:ref: strmm_ 14 15 13 13 13 13 4 4 6 6 4 6 4 124 124 124 124 \n 
 *:ref: sger_ 14 9 4 4 6 6 4 6 4 6 4 
 */
VEXTERNC int sneupd_(logical *rvec, char *howmny, logical *select, real *dr, real *di, real *z__, integer *ldz, real *sigmar, real *sigmai, real *workev, char *bmat, integer *n, char *which, integer *nev, real *tol, real *resid, integer *ncv, real *v, integer *ldv, integer *iparam, integer *ipntr, real *workd, real *workl, integer *lworkl, integer *info, ftnlen howmny_len, ftnlen bmat_len, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: second_ 14 1 6 \n 
 *:ref: ssortc_ 14 7 13 12 4 6 6 6 124 \n 
 *:ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *:ref: svout_ 14 6 4 4 6 4 13 124 
 */
VEXTERNC int sngets_(integer *ishift, char *which, integer *kev, integer *np, real *ritzr, real *ritzi, real *bounds, real *shiftr, real *shifti, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: slamch_ 6 2 13 124 \n 
 *:ref: second_ 14 1 6 \n 
 *:ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *:ref: svout_ 14 6 4 4 6 4 13 124 \n 
 *:ref: sgetv0_ 14 14 4 13 4 12 4 4 6 4 6 6 4 6 4 124 \n 
 *:ref: scopy_ 14 5 4 6 4 6 4 \n 
 *:ref: sscal_ 14 4 4 6 6 4 \n 
 *:ref: slascl_ 14 11 13 4 4 6 6 4 4 6 4 4 124 \n 
 *:ref: sdot_ 6 5 4 6 4 6 4 \n 
 *:ref: snrm2_ 6 3 4 6 4 \n 
 *:ref: sgemv_ 14 12 13 4 4 6 6 4 6 4 6 6 4 124 
 */
VEXTERNC int ssaitr_(integer *ido, char *bmat, integer *n, integer *k, integer *np, integer *mode, real *resid, real *rnorm, real *v, integer *ldv, real *h__, integer *ldh, integer *ipntr, real *workd, integer *info, ftnlen bmat_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: slamch_ 6 2 13 124 \n 
 *:ref: second_ 14 1 6 \n 
 *:ref: slaset_ 14 8 13 4 4 6 6 6 4 124 \n 
 *:ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *:ref: svout_ 14 6 4 4 6 4 13 124 \n 
 *:ref: slartg_ 14 5 6 6 6 6 6 \n 
 *:ref: sscal_ 14 4 4 6 6 4 \n 
 *:ref: sgemv_ 14 12 13 4 4 6 6 4 6 4 6 6 4 124 \n 
 *:ref: scopy_ 14 5 4 6 4 6 4 \n 
 *:ref: slacpy_ 14 8 13 4 4 6 4 6 4 124 \n 
 *:ref: saxpy_ 14 6 4 6 6 4 6 4 
 */
VEXTERNC int ssapps_(integer *n, integer *kev, integer *np, real *shift, real *v, integer *ldv, real *h__, integer *ldh, real *resid, real *q, integer *ldq, real *workd);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: second_ 14 1 6 \n 
 *:ref: slamch_ 6 2 13 124 \n 
 *:ref: sgetv0_ 14 14 4 13 4 12 4 4 6 4 6 6 4 6 4 124 \n 
 *:ref: ssaitr_ 14 16 4 13 4 4 4 4 6 6 6 4 6 4 4 6 4 124 \n 
 *:ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *:ref: svout_ 14 6 4 4 6 4 13 124 \n 
 *:ref: sseigt_ 14 8 6 4 6 4 6 6 6 4 \n 
 *:ref: scopy_ 14 5 4 6 4 6 4 \n 
 *:ref: ssgets_ 14 8 4 13 4 4 6 6 6 124 \n 
 *:ref: ssconv_ 14 5 4 6 6 6 4 \n 
 *:ref: ssortr_ 14 6 13 12 4 6 6 124 \n 
 *:ref: sswap_ 14 5 4 6 4 6 4 \n 
 *:ref: ssapps_ 14 12 4 4 4 6 6 4 6 4 6 6 4 6 \n 
 *:ref: sdot_ 6 5 4 6 4 6 4 \n 
 *:ref: snrm2_ 6 3 4 6 4 */
VEXTERNC int ssaup2_(integer *ido, char *bmat, integer *n, char *which, integer *nev, integer *np, real *tol, real *resid, integer *mode, integer *iupd, integer *ishift, integer *mxiter, real *v, integer *ldv, real *h__, integer *ldh, real *ritz, real *bounds, real *q, integer *ldq, real *workl, integer *ipntr, real *workd, integer *info, ftnlen bmat_len, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: sstats_ 14 0 \n 
 *:ref: second_ 14 1 6 \n 
 *:ref: slamch_ 6 2 13 124 \n 
 *:ref: ssaup2_ 14 26 4 13 4 13 4 4 6 6 4 4 4 4 6 4 6 4 6 6 6 4 6 4 6 4 124 124 \n 
 *:ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *:ref: svout_ 14 6 4 4 6 4 13 124 
 */
VEXTERNC int ssaupd_(integer *ido, char *bmat, integer *n, char *which, integer *nev, real *tol, real *resid, integer *ncv, real *v, integer *ldv, integer *iparam, integer *ipntr, real *workd, real *workl, integer *lworkl, integer *info, ftnlen bmat_len, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: second_ 14 1 6 \n 
 *:ref: slamch_ 6 2 13 124 
 */
VEXTERNC int ssconv_(integer *n, real *ritz, real *bounds, real *tol, integer *nconv);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: second_ 14 1 6 \n 
 *:ref: svout_ 14 6 4 4 6 4 13 124 \n 
 *:ref: scopy_ 14 5 4 6 4 6 4 \n 
 *:ref: sstqrb_ 14 6 4 6 6 6 6 4 
 */
VEXTERNC int sseigt_(real *rnorm, integer *n, real *h__, integer *ldh, real *eig, real *bounds, real *workl, integer *ierr);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *:ref: sswap_ 14 5 4 6 4 6 4 */
VEXTERNC int ssesrt_(char *which, logical *apply, integer *n, real *x, integer *na, real *a, integer *lda, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: slamch_ 6 2 13 124 \n 
 *:ref: snrm2_ 6 3 4 6 4 \n 
 *:ref: svout_ 14 6 4 4 6 4 13 124 \n 
 *:ref: ssgets_ 14 8 4 13 4 4 6 6 6 124 \n 
 *:ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *:ref: scopy_ 14 5 4 6 4 6 4 \n 
 *:ref: ssteqr_ 14 9 13 4 6 6 6 4 6 4 124 \n 
 *:ref: ssesrt_ 14 8 13 12 4 6 4 6 4 124 \n 
 *:ref: ssortr_ 14 6 13 12 4 6 6 124 \n 
 *:ref: sscal_ 14 4 4 6 6 4 \n 
 *:ref: sgeqr2_ 14 7 4 4 6 4 6 6 4 \n 
 *:ref: sorm2r_ 14 14 13 13 4 4 4 6 4 6 6 4 6 4 124 124 \n 
 *:ref: slacpy_ 14 8 13 4 4 6 4 6 4 124 \n 
 *:ref: sger_ 14 9 4 4 6 6 4 6 4 6 4 
 */
VEXTERNC int sseupd_(logical *rvec, char *howmny, logical *select, real *d__, real *z__, integer *ldz, real *sigma, char *bmat, integer *n, char *which, integer *nev, real *tol, real *resid, integer *ncv, real *v, integer *ldv, integer *iparam, integer *ipntr, real *workd, real *workl, integer *lworkl, integer *info, ftnlen howmny_len, ftnlen bmat_len, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: second_ 14 1 6 \n 
 *:ref: ssortr_ 14 6 13 12 4 6 6 124 \n 
 *:ref: sswap_ 14 5 4 6 4 6 4 \n 
 *:ref: scopy_ 14 5 4 6 4 6 4 \n 
 *:ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *:ref: svout_ 14 6 4 4 6 4 13 124
 */
VEXTERNC int ssgets_(integer *ishift, char *which, integer *kev, integer *np, real *ritz, real *bounds, real *shifts, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *:ref: slapy2_ 6 2 6 6 */
VEXTERNC int ssortc_(char *which, logical *apply, integer *n, real *xreal, real *ximag, real *y, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
*/
VEXTERNC int ssortr_(char *which, logical *apply, integer *n, real *x1, real *x2, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen timing_ 124 */
VEXTERNC int sstatn_(void);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen timing_ 124 */
VEXTERNC int sstats_(void);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *:ref: slamch_ 6 2 13 124 \n 
 *:ref: slanst_ 6 5 13 4 6 6 124 \n 
 *:ref: slascl_ 14 11 13 4 4 6 6 4 4 6 4 4 124 \n 
 *:ref: slaev2_ 14 7 6 6 6 6 6 6 6 \n 
 *:ref: slae2_ 14 5 6 6 6 6 6 \n 
 *:ref: slapy2_ 6 2 6 6 \n 
 *:ref: slartg_ 14 5 6 6 6 6 6 \n 
 *:ref: slasr_ 14 12 13 13 13 4 4 6 6 6 4 124 124 124 \n 
 *:ref: slasrt_ 14 5 13 4 6 4 124 
 */
VEXTERNC int sstqrb_(integer *n, real *d__, real *e, real *z__, real *work, integer *info);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: second_ 14 1 6 \n 
 *:ref: zlarnv_ 14 4 4 4 4 9 \n 
 *:ref: zcopy_ 14 5 4 9 4 9 4 \n 
 *:ref: zdotc_ 9 6 9 4 9 4 9 4 \n 
 *:ref: dlapy2_ 7 2 7 7 \n 
 *:ref: dznrm2_ 7 3 4 9 4 \n 
 *:ref: zgemv_ 14 12 13 4 4 9 9 4 9 4 9 9 4 124 \n 
 *:ref: dvout_ 14 6 4 4 7 4 13 124 \n 
 *:ref: zvout_ 14 6 4 4 9 4 13 124 
 */
VEXTERNC int zgetv0_(integer *ido, char *bmat, integer *itry, logical *initv, integer *n, integer *j, doublecomplex *v, integer *ldv, doublecomplex *resid, doublereal *rnorm, integer *ipntr, doublecomplex *workd, integer *ierr, ftnlen bmat_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: dlamch_ 7 2 13 124 \n 
 *:ref: dlabad_ 14 2 7 7 \n 
 *:ref: second_ 14 1 6 \n 
 *:ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *:ref: dvout_ 14 6 4 4 7 4 13 124 \n 
 *:ref: zgetv0_ 14 14 4 13 4 12 4 4 9 4 9 7 4 9 4 124 \n 
 *:ref: zcopy_ 14 5 4 9 4 9 4 \n 
 *:ref: zdscal_ 14 4 4 7 9 4 \n 
 *:ref: zlascl_ 14 11 13 4 4 7 7 4 4 9 4 4 124 \n 
 *:ref: zdotc_ 9 6 9 4 9 4 9 4 \n 
 *:ref: dlapy2_ 7 2 7 7 \n 
 *:ref: dznrm2_ 7 3 4 9 4 \n 
 *:ref: zgemv_ 14 12 13 4 4 9 9 4 9 4 9 9 4 124 \n 
 *:ref: zvout_ 14 6 4 4 9 4 13 124 \n 
 *:ref: zaxpy_ 14 6 4 9 9 4 9 4 \n 
 *:ref: zlanhs_ 7 6 13 4 9 4 9 124 \n 
 *:ref: zmout_ 14 8 4 4 4 9 4 4 13 124 
 */
VEXTERNC int znaitr_(integer *ido, char *bmat, integer *n, integer *k, integer *np, integer *nb, doublecomplex *resid, doublereal *rnorm, doublecomplex *v, integer *ldv, doublecomplex *h__, integer *ldh, integer *ipntr, doublecomplex *workd, integer *info, ftnlen bmat_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: dlamch_ 7 2 13 124 \n 
 *:ref: dlabad_ 14 2 7 7 \n 
 *:ref: second_ 14 1 6 \n 
 *:ref: zlaset_ 14 8 13 4 4 9 9 9 4 124 \n 
 *:ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *:ref: zvout_ 14 6 4 4 9 4 13 124 \n 
 *:ref: zlanhs_ 7 6 13 4 9 4 9 124 \n 
 *:ref: zlartg_ 14 5 9 9 7 9 9 \n 
 *:ref: dlapy2_ 7 2 7 7 \n 
 *:ref: zscal_ 14 4 4 9 9 4 \n 
 *:ref: zgemv_ 14 12 13 4 4 9 9 4 9 4 9 9 4 124 \n 
 *:ref: zcopy_ 14 5 4 9 4 9 4 \n 
 *:ref: zlacpy_ 14 8 13 4 4 9 4 9 4 124 \n 
 *:ref: zaxpy_ 14 6 4 9 9 4 9 4 \n 
 *:ref: zmout_ 14 8 4 4 4 9 4 4 13 124 
 */
VEXTERNC int znapps_(integer *n, integer *kev, integer *np, doublecomplex *shift, doublecomplex *v, integer *ldv, doublecomplex *h__, integer *ldh, doublecomplex *resid, doublecomplex *q, integer *ldq, doublecomplex *workl, doublecomplex *workd);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: second_ 14 1 6 \n 
 *:ref: dlamch_ 7 2 13 124 \n 
 *:ref: zgetv0_ 14 14 4 13 4 12 4 4 9 4 9 7 4 9 4 124 \n 
 *:ref: znaitr_ 14 16 4 13 4 4 4 4 9 7 9 4 9 4 4 9 4 124 \n 
 *:ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *:ref: dvout_ 14 6 4 4 7 4 13 124 \n 
 *:ref: zneigh_ 14 11 7 4 9 4 9 9 9 4 9 7 4 \n 
 *:ref: zcopy_ 14 5 4 9 4 9 4 \n 
 *:ref: zngets_ 14 7 4 13 4 4 9 9 124 \n 
 *:ref: dlapy2_ 7 2 7 7 \n 
 *:ref: zvout_ 14 6 4 4 9 4 13 124 \n 
 *:ref: zsortc_ 14 6 13 12 4 9 9 124 \n 
 *:ref: znapps_ 14 13 4 4 4 9 9 4 9 4 9 9 4 9 9 \n 
 *:ref: zdotc_ 9 6 9 4 9 4 9 4 \n 
 *:ref: dznrm2_ 7 3 4 9 4 \n 
 *:ref: zmout_ 14 8 4 4 4 9 4 4 13 124 
 */
VEXTERNC int znaup2_(integer *ido, char *bmat, integer *n, char *which, integer *nev, integer *np, doublereal *tol, doublecomplex *resid, integer *mode, integer *iupd, integer *ishift, integer *mxiter, doublecomplex *v, integer *ldv, doublecomplex *h__, integer *ldh, doublecomplex *ritz, doublecomplex *bounds, doublecomplex *q, integer *ldq, doublecomplex *workl, integer *ipntr, doublecomplex *workd, doublereal *rwork, integer *info, ftnlen bmat_len, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: zstatn_ 14 0 \n 
 *:ref: second_ 14 1 6 \n 
 *:ref: dlamch_ 7 2 13 124 \n 
 *:ref: znaup2_ 14 27 4 13 4 13 4 4 7 9 4 4 4 4 9 4 9 4 9 9 9 4 9 4 9 7 4 124 124 \n 
 *:ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *:ref: zvout_ 14 6 4 4 9 4 13 124 
 */
VEXTERNC int znaupd_(integer *ido, char *bmat, integer *n, char *which, integer *nev, doublereal *tol, doublecomplex *resid, integer *ncv, doublecomplex *v, integer *ldv, integer *iparam, integer *ipntr, doublecomplex *workd, doublecomplex *workl, integer *lworkl, doublereal *rwork, integer *info, ftnlen bmat_len, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: second_ 14 1 6 \n 
 *:ref: zmout_ 14 8 4 4 4 9 4 4 13 124 \n 
 *:ref: zlacpy_ 14 8 13 4 4 9 4 9 4 124 \n 
 *:ref: zlaset_ 14 8 13 4 4 9 9 9 4 124 \n 
 *:ref: zlahqr_ 14 13 12 12 4 4 4 9 4 9 4 4 9 4 4 \n 
 *:ref: zcopy_ 14 5 4 9 4 9 4 \n 
 *:ref: zvout_ 14 6 4 4 9 4 13 124 \n 
 *:ref: ztrevc_ 14 17 13 13 12 4 9 4 9 4 9 4 4 4 9 7 4 124 124 \n 
 *:ref: dznrm2_ 7 3 4 9 4 \n 
 *:ref: zdscal_ 14 4 4 7 9 4 
 */
VEXTERNC int zneigh_(doublereal *rnorm, integer *n, doublecomplex *h__, integer *ldh, doublecomplex *ritz, doublecomplex *bounds, doublecomplex *q, integer *ldq, doublecomplex *workl, doublereal *rwork, integer *ierr);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: dlamch_ 7 2 13 124 \n 
 *:ref: zvout_ 14 6 4 4 9 4 13 124 \n 
 *:ref: zngets_ 14 7 4 13 4 4 9 9 124 \n 
 *:ref: dlapy2_ 7 2 7 7 \n 
 *:ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *:ref: zcopy_ 14 5 4 9 4 9 4 \n 
 *:ref: zlaset_ 14 8 13 4 4 9 9 9 4 124 \n 
 *:ref: zlahqr_ 14 13 12 12 4 4 4 9 4 9 4 4 9 4 4 \n 
 *:ref: zmout_ 14 8 4 4 4 9 4 4 13 124 \n 
 *:ref: ztrsen_ 14 17 13 13 12 4 9 4 9 4 9 4 7 7 9 4 4 124 124 \n 
 *:ref: zgeqr2_ 14 7 4 4 9 4 9 9 4 \n 
 *:ref: zunm2r_ 14 14 13 13 4 4 4 9 4 9 9 4 9 4 124 124 \n 
 *:ref: zlacpy_ 14 8 13 4 4 9 4 9 4 124 \n 
 *:ref: zscal_ 14 4 4 9 9 4 \n 
 *:ref: ztrevc_ 14 17 13 13 12 4 9 4 9 4 9 4 4 4 9 7 4 124 124 \n 
 *:ref: dznrm2_ 7 3 4 9 4 \n 
 *:ref: zdscal_ 14 4 4 7 9 4 \n 
 *:ref: zdotc_ 9 6 9 4 9 4 9 4 \n 
 *:ref: ztrmm_ 14 15 13 13 13 13 4 4 9 9 4 9 4 124 124 124 124 \n 
 *:ref: zgeru_ 14 9 4 4 9 9 4 9 4 9 4 
 */
VEXTERNC int zneupd_(logical *rvec, char *howmny, logical *select, doublecomplex *d__, doublecomplex *z__, integer *ldz, doublecomplex *sigma, doublecomplex *workev, char *bmat, integer *n, char *which, integer *nev, doublereal *tol, doublecomplex *resid, integer *ncv, doublecomplex *v, integer *ldv, integer *iparam, integer *ipntr, doublecomplex *workd, doublecomplex *workl, integer *lworkl, doublereal *rwork, integer *info, ftnlen howmny_len, ftnlen bmat_len, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 * comlen debug_ 96 \n 
 * comlen timing_ 124 \n 
 *:ref: second_ 14 1 6 \n 
 *:ref: zsortc_ 14 6 13 12 4 9 9 124 \n 
 *:ref: ivout_ 14 6 4 4 4 4 13 124 \n 
 *:ref: zvout_ 14 6 4 4 9 4 13 124 
 */
VEXTERNC int zngets_(integer *ishift, char *which, integer *kev, integer *np, doublecomplex *ritz, doublecomplex *bounds, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *:ref: dlapy2_ 7 2 7 7 */
VEXTERNC int zsortc_(char *which, logical *apply, integer *n, doublecomplex *x, doublecomplex *y, ftnlen which_len);

/** @brief Library VARPACK prototypes 
 *  @author Michael Holst
 *  @note
 *  comlen timing_ 124 */
VEXTERNC int zstatn_(void);

#endif /* _VARPACK_H_ */

