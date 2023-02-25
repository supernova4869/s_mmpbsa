/*
 * ***************************************************************************
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
 * rcsid="$Id: vcgcode.h,v 1.7 2010/08/12 05:52:32 fetk Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     vcgcode.h
 *
 * Purpose:  The primary header for CgCode.
 *
 * Notes:    We provide this header whether or not we provide the BLAS
 *           library itself.  This gives us some compile-time type-checking
 *           even for architecture-dependent assembly-coded BLAS.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#ifndef _VCGCODE_H_
#define _VCGCODE_H_

#include <punc/punc_base.h>

#include <punc/vf2c.h>

/*
 * ***************************************************************************
 * Library VCGCODE prototypes
 * ***************************************************************************
 */

VEXTERNC doublereal d1mach_(integer *idum);
VEXTERNC int dcbfix_(S_fp matvec, S_fp pcondl, doublereal *a, integer *ia, doublereal *c__, integer *ic, integer *ipc, doublereal *aa, doublereal *bb, integer *nsteps, doublereal *b, doublereal *x, doublereal *r__, doublereal *dx, doublereal *work, integer *n);
/*:ref: dcopy_ 14 5 4 7 4 7 4 */
VEXTERNC int dcg_(S_fp matvec, doublereal *a, integer *ia, doublereal *x, doublereal *b, integer *n, integer *iparam, doublereal *rparam, integer *iwork, doublereal *r__, doublereal *ap, doublereal *d__, doublereal *e, doublereal *cndwk, integer *ierror);
/*:ref: dcgchk_ 14 3 4 7 4 */
/*:ref: d1mach_ 7 1 4 */
/*:ref: dnrm2_ 7 3 4 7 4 */
/*:ref: mdstop_ 4 16 4 4 4 7 7 4 7 7 7 4 7 7 7 7 7 4 */
/*:ref: dcopy_ 14 5 4 7 4 7 4 */
/*:ref: ddot_ 7 5 4 7 4 7 4 */
/*:ref: daxpy_ 14 6 4 7 7 4 7 4 */
/*:ref: donest_ 14 10 4 7 7 7 7 4 4 7 7 7 */
VEXTERNC int dcgchk_(integer *ipar, doublereal *rpar, integer *n);
/*:ref: d1mach_ 7 1 4 */
VEXTERNC int dcgdrv_(U_fp matvec, U_fp pcondl, U_fp pcondr, doublereal *a, integer *ia, doublereal *x, doublereal *b, integer *n, doublereal *q, integer *iq, doublereal *p, integer *ip, integer *iparam, doublereal *rparam, integer *iwork, doublereal *rwork, integer *ierror);
/*:ref: dcg_ 14 15 200 7 4 7 7 4 4 7 4 7 7 7 7 7 4 */
/*:ref: dcr_ 14 16 200 7 4 7 7 4 4 7 4 7 7 7 7 7 7 4 */
/*:ref: dcrind_ 14 18 200 7 4 7 7 4 4 7 4 7 7 7 7 7 7 7 7 4 */
/*:ref: dpcg_ 14 18 200 200 7 4 7 7 4 7 4 4 7 4 7 7 7 7 7 4 */
/*:ref: dcgnr_ 14 15 200 7 4 7 7 4 4 7 4 7 7 7 7 7 4 */
/*:ref: dcgne_ 14 15 200 7 4 7 7 4 4 7 4 7 7 7 7 7 4 */
/*:ref: dpcgnr_ 14 19 200 200 7 4 7 7 4 7 4 4 7 4 7 7 7 7 7 7 4 */
/*:ref: dpcgne_ 14 19 200 200 7 4 7 7 4 7 4 4 7 4 7 7 7 7 7 7 4 */
/*:ref: dppcg_ 14 21 200 200 7 4 7 7 4 7 4 4 7 4 7 7 7 7 7 7 7 7 4 */
/*:ref: dpcgca_ 14 21 200 200 7 4 7 7 4 7 4 4 7 4 7 7 7 7 7 7 7 7 4 */
VEXTERNC int dcgne_(S_fp matvec, doublereal *a, integer *ia, doublereal *x, doublereal *b, integer *n, integer *iparam, doublereal *rparam, integer *iwork, doublereal *r__, doublereal *h__, doublereal *d__, doublereal *e, doublereal *cndwk, integer *ierror);
/*:ref: dcgchk_ 14 3 4 7 4 */
/*:ref: d1mach_ 7 1 4 */
/*:ref: dnrm2_ 7 3 4 7 4 */
/*:ref: mdstop_ 4 16 4 4 4 7 7 4 7 7 7 4 7 7 7 7 7 4 */
/*:ref: daxpy_ 14 6 4 7 7 4 7 4 */
/*:ref: donest_ 14 10 4 7 7 7 7 4 4 7 7 7 */
VEXTERNC int dcgnr_(S_fp matvec, doublereal *a, integer *ia, doublereal *x, doublereal *b, integer *n, integer *iparam, doublereal *rparam, integer *iwork, doublereal *r__, doublereal *h__, doublereal *d__, doublereal *e, doublereal *cndwk, integer *ierror);
/*:ref: dcgchk_ 14 3 4 7 4 */
/*:ref: d1mach_ 7 1 4 */
/*:ref: dnrm2_ 7 3 4 7 4 */
/*:ref: mdstop_ 4 16 4 4 4 7 7 4 7 7 7 4 7 7 7 7 7 4 */
/*:ref: daxpy_ 14 6 4 7 7 4 7 4 */
/*:ref: donest_ 14 10 4 7 7 7 7 4 4 7 7 7 */
VEXTERNC int dckchb_(integer *ipar, doublereal *rpar, doublereal *pegmin, doublereal *pegmax, doublereal *condca);
/*:ref: d1mach_ 7 1 4 */
/*:ref: ddpchb_ 14 9 4 7 7 4 7 7 7 7 4 */
VEXTERNC int dcr_(S_fp matvec, doublereal *a, integer *ia, doublereal *x, doublereal *b, integer *n, integer *iparam, doublereal *rparam, integer *iwork, doublereal *r__, doublereal *ar, doublereal *ap, doublereal *d__, doublereal *e, doublereal *cndwk, integer *ierror);
/*:ref: dcgchk_ 14 3 4 7 4 */
/*:ref: d1mach_ 7 1 4 */
/*:ref: dnrm2_ 7 3 4 7 4 */
/*:ref: mdstop_ 4 16 4 4 4 7 7 4 7 7 7 4 7 7 7 7 7 4 */
/*:ref: ddot_ 7 5 4 7 4 7 4 */
/*:ref: dcopy_ 14 5 4 7 4 7 4 */
/*:ref: daxpy_ 14 6 4 7 7 4 7 4 */
/*:ref: donest_ 14 10 4 7 7 7 7 4 4 7 7 7 */
VEXTERNC int dcrind_(S_fp matvec, doublereal *a, integer *ia, doublereal *x, doublereal *b, integer *n, integer *iparam, doublereal *rparam, integer *iwork, doublereal *r__, doublereal *aap, doublereal *ap, doublereal *pold, doublereal *apold, doublereal *d__, doublereal *e, doublereal *cndwk, integer *ierror);
/*:ref: d1mach_ 7 1 4 */
/*:ref: dcgchk_ 14 3 4 7 4 */
/*:ref: dnrm2_ 7 3 4 7 4 */
/*:ref: mdstop_ 4 16 4 4 4 7 7 4 7 7 7 4 7 7 7 7 7 4 */
/*:ref: ddot_ 7 5 4 7 4 7 4 */
VEXTERNC int ddpchb_(integer *iounit, doublereal *aa, doublereal *bb, integer *n, doublereal *cl, doublereal *cu, doublereal *cond, doublereal *ereps, integer *iadapt);
/*:ref: d1mach_ 7 1 4 */
VEXTERNC doublereal depsln_(doublereal *x);
VEXTERNC int donest_(integer *iounit, doublereal *d__, doublereal *e, doublereal *w1, doublereal *w2, integer *ind, integer *nt, doublereal *eigmin, doublereal *eigmax, doublereal *cond);
/*:ref: sratqr_ 14 12 4 7 7 7 7 4 7 4 7 12 4 4 */
VEXTERNC int dpcg_(S_fp matvec, S_fp pcondl, doublereal *a, integer *ia, doublereal *x, doublereal *b, integer *n, doublereal *q, integer *iq, integer *iparam, doublereal *rparam, integer *iwork, doublereal *r__, doublereal *h__, doublereal *d__, doublereal *e, doublereal *cndwk, integer *ierror);
/*:ref: dcgchk_ 14 3 4 7 4 */
/*:ref: d1mach_ 7 1 4 */
/*:ref: dnrm2_ 7 3 4 7 4 */
/*:ref: mdstop_ 4 16 4 4 4 7 7 4 7 7 7 4 7 7 7 7 7 4 */
/*:ref: ddot_ 7 5 4 7 4 7 4 */
/*:ref: daxpy_ 14 6 4 7 7 4 7 4 */
/*:ref: donest_ 14 10 4 7 7 7 7 4 4 7 7 7 */
VEXTERNC int dpcgca_(S_fp matvec, S_fp pcondl, doublereal *a, integer *ia, doublereal *x, doublereal *b, integer *n, doublereal *q, integer *iq, integer *iparam, doublereal *rparam, integer *iwork, doublereal *p, doublereal *h__, doublereal *cap, doublereal *w1, doublereal *w2, doublereal *w3, doublereal *d__, doublereal *e, integer *ierror);
/*:ref: dcgchk_ 14 3 4 7 4 */
/*:ref: dckchb_ 14 5 4 7 7 7 7 */
/*:ref: dnrm2_ 7 3 4 7 4 */
/*:ref: dcbfix_ 14 16 214 214 7 4 7 4 4 7 7 4 7 7 7 7 7 4 */
/*:ref: ddot_ 7 5 4 7 4 7 4 */
/*:ref: dcopy_ 14 5 4 7 4 7 4 */
/*:ref: mdstop_ 4 16 4 4 4 7 7 4 7 7 7 4 7 7 7 7 7 4 */
/*:ref: donest_ 14 10 4 7 7 7 7 4 4 7 7 7 */
/*:ref: ddpchb_ 14 9 4 7 7 4 7 7 7 7 4 */
VEXTERNC int dpcgne_(S_fp matvec, S_fp pcondl, doublereal *a, integer *ia, doublereal *x, doublereal *b, integer *n, doublereal *q, integer *iq, integer *iparam, doublereal *rparam, integer *iwork, doublereal *r__, doublereal *h__, doublereal *ap, doublereal *d__, doublereal *e, doublereal *cndwk, integer *ierror);
/*:ref: dcgchk_ 14 3 4 7 4 */
/*:ref: d1mach_ 7 1 4 */
/*:ref: dnrm2_ 7 3 4 7 4 */
/*:ref: mdstop_ 4 16 4 4 4 7 7 4 7 7 7 4 7 7 7 7 7 4 */
/*:ref: daxpy_ 14 6 4 7 7 4 7 4 */
/*:ref: donest_ 14 10 4 7 7 7 7 4 4 7 7 7 */
VEXTERNC int dpcgnr_(S_fp matvec, S_fp pcondl, doublereal *a, integer *ia, doublereal *x, doublereal *b, integer *n, doublereal *q, integer *iq, integer *iparam, doublereal *rparam, integer *iwork, doublereal *r__, doublereal *h__, doublereal *ap, doublereal *d__, doublereal *e, doublereal *cndwk, integer *ierror);
/*:ref: dcgchk_ 14 3 4 7 4 */
/*:ref: d1mach_ 7 1 4 */
/*:ref: dnrm2_ 7 3 4 7 4 */
/*:ref: mdstop_ 4 16 4 4 4 7 7 4 7 7 7 4 7 7 7 7 7 4 */
/*:ref: daxpy_ 14 6 4 7 7 4 7 4 */
/*:ref: donest_ 14 10 4 7 7 7 7 4 4 7 7 7 */
VEXTERNC int dppcg_(S_fp matvec, S_fp pcondl, doublereal *a, integer *ia, doublereal *x, doublereal *b, integer *n, doublereal *q, integer *iq, integer *iparam, doublereal *rparam, integer *iwork, doublereal *r__, doublereal *h__, doublereal *w1, doublereal *w2, doublereal *w3, doublereal *w4, doublereal *d__, doublereal *e, integer *ierror);
/*:ref: dcgchk_ 14 3 4 7 4 */
/*:ref: dckchb_ 14 5 4 7 7 7 7 */
/*:ref: dnrm2_ 7 3 4 7 4 */
/*:ref: dcbfix_ 14 16 214 214 7 4 7 4 4 7 7 4 7 7 7 7 7 4 */
/*:ref: ddot_ 7 5 4 7 4 7 4 */
/*:ref: mdstop_ 4 16 4 4 4 7 7 4 7 7 7 4 7 7 7 7 7 4 */
/*:ref: daxpy_ 14 6 4 7 7 4 7 4 */
/*:ref: donest_ 14 10 4 7 7 7 7 4 4 7 7 7 */
/*:ref: ddpchb_ 14 9 4 7 7 4 7 7 7 7 4 */
VEXTERNC int dratqr_(integer *n, doublereal *eps1, doublereal *d__, doublereal *e, doublereal *e2, integer *m, doublereal *w, integer *ind, doublereal *bd, logical *type__, integer *idef, integer *ierr);
/*:ref: depsln_ 7 1 7 */
VEXTERNC integer mdstop_(integer *istop, integer *iters, integer *itmax, doublereal *errtol, doublereal *stptst, integer *ierror, doublereal *r__, doublereal *s, doublereal *z__, integer *n, doublereal *rnorm, doublereal *snorm, doublereal *znorm, doublereal *denom, doublereal *conda, integer *ido);
/*:ref: dnrm2_ 7 3 4 7 4 */
VEXTERNC integer msstop_(integer *istop, integer *iters, integer *itmax, real *errtol, real *stptst, integer *ierror, real *r__, real *s, real *z__, integer *n, real *rnorm, real *snorm, real *znorm, real *denom, real *conda, integer *ido);
/*:ref: snrm2_ 6 3 4 6 4 */
VEXTERNC E_f r1mach_(integer *idum);
VEXTERNC int scbfix_(S_fp matvec, S_fp pcondl, real *a, integer *ia, real *c__, integer *ic, integer *ipc, real *aa, real *bb, integer *nsteps, real *b, real *x, real *r__, real *dx, real *work, integer *n);
/*:ref: scopy_ 14 5 4 6 4 6 4 */
VEXTERNC int scg_(S_fp matvec, real *a, integer *ia, real *x, real *b, integer *n, integer *iparam, real *rparam, integer *iwork, real *r__, real *ap, real *d__, real *e, real *cndwk, integer *ierror);
/*:ref: scgchk_ 14 3 4 6 4 */
/*:ref: r1mach_ 6 1 4 */
/*:ref: snrm2_ 6 3 4 6 4 */
/*:ref: msstop_ 4 16 4 4 4 6 6 4 6 6 6 4 6 6 6 6 6 4 */
/*:ref: scopy_ 14 5 4 6 4 6 4 */
/*:ref: sdot_ 6 5 4 6 4 6 4 */
/*:ref: saxpy_ 14 6 4 6 6 4 6 4 */
/*:ref: sonest_ 14 10 4 6 6 6 6 4 4 6 6 6 */
VEXTERNC int scgchk_(integer *ipar, real *rpar, integer *n);
/*:ref: r1mach_ 6 1 4 */
VEXTERNC int scgdrv_(U_fp matvec, U_fp pcondl, U_fp pcondr, real *a, integer *ia, real *x, real *b, integer *n, real *q, integer *iq, real *p, integer *ip, integer *iparam, real *rparam, integer *iwork, real *rwork, integer *ierror);
/*:ref: scg_ 14 15 200 6 4 6 6 4 4 6 4 6 6 6 6 6 4 */
/*:ref: scr_ 14 16 200 6 4 6 6 4 4 6 4 6 6 6 6 6 6 4 */
/*:ref: scrind_ 14 18 200 6 4 6 6 4 4 6 4 6 6 6 6 6 6 6 6 4 */
/*:ref: spcg_ 14 18 200 200 6 4 6 6 4 6 4 4 6 4 6 6 6 6 6 4 */
/*:ref: scgnr_ 14 15 200 6 4 6 6 4 4 6 4 6 6 6 6 6 4 */
/*:ref: scgne_ 14 15 200 6 4 6 6 4 4 6 4 6 6 6 6 6 4 */
/*:ref: spcgnr_ 14 19 200 200 6 4 6 6 4 6 4 4 6 4 6 6 6 6 6 6 4 */
/*:ref: spcgne_ 14 19 200 200 6 4 6 6 4 6 4 4 6 4 6 6 6 6 6 6 4 */
/*:ref: sppcg_ 14 21 200 200 6 4 6 6 4 6 4 4 6 4 6 6 6 6 6 6 6 6 4 */
/*:ref: spcgca_ 14 21 200 200 6 4 6 6 4 6 4 4 6 4 6 6 6 6 6 6 6 6 4 */
VEXTERNC int scgne_(S_fp matvec, real *a, integer *ia, real *x, real *b, integer *n, integer *iparam, real *rparam, integer *iwork, real *r__, real *h__, real *d__, real *e, real *cndwk, integer *ierror);
/*:ref: scgchk_ 14 3 4 6 4 */
/*:ref: r1mach_ 6 1 4 */
/*:ref: snrm2_ 6 3 4 6 4 */
/*:ref: msstop_ 4 16 4 4 4 6 6 4 6 6 6 4 6 6 6 6 6 4 */
/*:ref: saxpy_ 14 6 4 6 6 4 6 4 */
/*:ref: sonest_ 14 10 4 6 6 6 6 4 4 6 6 6 */
VEXTERNC int scgnr_(S_fp matvec, real *a, integer *ia, real *x, real *b, integer *n, integer *iparam, real *rparam, integer *iwork, real *r__, real *h__, real *d__, real *e, real *cndwk, integer *ierror);
/*:ref: scgchk_ 14 3 4 6 4 */
/*:ref: r1mach_ 6 1 4 */
/*:ref: snrm2_ 6 3 4 6 4 */
/*:ref: msstop_ 4 16 4 4 4 6 6 4 6 6 6 4 6 6 6 6 6 4 */
/*:ref: saxpy_ 14 6 4 6 6 4 6 4 */
/*:ref: sonest_ 14 10 4 6 6 6 6 4 4 6 6 6 */
VEXTERNC int sckchb_(integer *ipar, real *rpar, real *pegmin, real *pegmax, real *condca);
/*:ref: r1mach_ 6 1 4 */
/*:ref: sdpchb_ 14 9 4 6 6 4 6 6 6 6 4 */
VEXTERNC int scr_(S_fp matvec, real *a, integer *ia, real *x, real *b, integer *n, integer *iparam, real *rparam, integer *iwork, real *r__, real *ar, real *ap, real *d__, real *e, real *cndwk, integer *ierror);
/*:ref: scgchk_ 14 3 4 6 4 */
/*:ref: r1mach_ 6 1 4 */
/*:ref: snrm2_ 6 3 4 6 4 */
/*:ref: msstop_ 4 16 4 4 4 6 6 4 6 6 6 4 6 6 6 6 6 4 */
/*:ref: sdot_ 6 5 4 6 4 6 4 */
/*:ref: scopy_ 14 5 4 6 4 6 4 */
/*:ref: saxpy_ 14 6 4 6 6 4 6 4 */
/*:ref: sonest_ 14 10 4 6 6 6 6 4 4 6 6 6 */
VEXTERNC int scrind_(S_fp matvec, real *a, integer *ia, real *x, real *b, integer *n, integer *iparam, real *rparam, integer *iwork, real *r__, real *aap, real *ap, real *pold, real *apold, real *d__, real *e, real *cndwk, integer *ierror);
/*:ref: r1mach_ 6 1 4 */
/*:ref: scgchk_ 14 3 4 6 4 */
/*:ref: snrm2_ 6 3 4 6 4 */
/*:ref: msstop_ 4 16 4 4 4 6 6 4 6 6 6 4 6 6 6 6 6 4 */
/*:ref: sdot_ 6 5 4 6 4 6 4 */
VEXTERNC int sdpchb_(integer *iounit, real *aa, real *bb, integer *n, real *cl, real *cu, real *cond, real *ereps, integer *iadapt);
/*:ref: r1mach_ 6 1 4 */
VEXTERNC int sonest_(integer *iounit, real *d__, real *e, real *w1, real *w2, integer *ind, integer *nt, real *eigmin, real *eigmax, real *cond);
/*:ref: sratqr_ 14 12 4 6 6 6 6 4 6 4 6 12 4 4 */
VEXTERNC int spcg_(S_fp matvec, S_fp pcondl, real *a, integer *ia, real *x, real *b, integer *n, real *q, integer *iq, integer *iparam, real *rparam, integer *iwork, real *r__, real *h__, real *d__, real *e, real *cndwk, integer *ierror);
/*:ref: scgchk_ 14 3 4 6 4 */
/*:ref: r1mach_ 6 1 4 */
/*:ref: snrm2_ 6 3 4 6 4 */
/*:ref: msstop_ 4 16 4 4 4 6 6 4 6 6 6 4 6 6 6 6 6 4 */
/*:ref: sdot_ 6 5 4 6 4 6 4 */
/*:ref: saxpy_ 14 6 4 6 6 4 6 4 */
/*:ref: sonest_ 14 10 4 6 6 6 6 4 4 6 6 6 */
VEXTERNC int spcgca_(S_fp matvec, S_fp pcondl, real *a, integer *ia, real *x, real *b, integer *n, real *q, integer *iq, integer *iparam, real *rparam, integer *iwork, real *p, real *h__, real *cap, real *w1, real *w2, real *w3, real *d__, real *e, integer *ierror);
/*:ref: scgchk_ 14 3 4 6 4 */
/*:ref: sckchb_ 14 5 4 6 6 6 6 */
/*:ref: snrm2_ 6 3 4 6 4 */
/*:ref: scbfix_ 14 16 214 214 6 4 6 4 4 6 6 4 6 6 6 6 6 4 */
/*:ref: sdot_ 6 5 4 6 4 6 4 */
/*:ref: scopy_ 14 5 4 6 4 6 4 */
/*:ref: msstop_ 4 16 4 4 4 6 6 4 6 6 6 4 6 6 6 6 6 4 */
/*:ref: sonest_ 14 10 4 6 6 6 6 4 4 6 6 6 */
/*:ref: sdpchb_ 14 9 4 6 6 4 6 6 6 6 4 */
VEXTERNC int spcgne_(S_fp matvec, S_fp pcondl, real *a, integer *ia, real *x, real *b, integer *n, real *q, integer *iq, integer *iparam, real *rparam, integer *iwork, real *r__, real *h__, real *ap, real *d__, real *e, real *cndwk, integer *ierror);
/*:ref: scgchk_ 14 3 4 6 4 */
/*:ref: r1mach_ 6 1 4 */
/*:ref: snrm2_ 6 3 4 6 4 */
/*:ref: msstop_ 4 16 4 4 4 6 6 4 6 6 6 4 6 6 6 6 6 4 */
/*:ref: saxpy_ 14 6 4 6 6 4 6 4 */
/*:ref: sonest_ 14 10 4 6 6 6 6 4 4 6 6 6 */
VEXTERNC int spcgnr_(S_fp matvec, S_fp pcondl, real *a, integer *ia, real *x, real *b, integer *n, real *q, integer *iq, integer *iparam, real *rparam, integer *iwork, real *r__, real *h__, real *ap, real *d__, real *e, real *cndwk, integer *ierror);
/*:ref: scgchk_ 14 3 4 6 4 */
/*:ref: r1mach_ 6 1 4 */
/*:ref: snrm2_ 6 3 4 6 4 */
/*:ref: msstop_ 4 16 4 4 4 6 6 4 6 6 6 4 6 6 6 6 6 4 */
/*:ref: saxpy_ 14 6 4 6 6 4 6 4 */
/*:ref: sonest_ 14 10 4 6 6 6 6 4 4 6 6 6 */
VEXTERNC int sppcg_(S_fp matvec, S_fp pcondl, real *a, integer *ia, real *x, real *b, integer *n, real *q, integer *iq, integer *iparam, real *rparam, integer *iwork, real *r__, real *h__, real *w1, real *w2, real *w3, real *w4, real *d__, real *e, integer *ierror);
/*:ref: scgchk_ 14 3 4 6 4 */
/*:ref: sckchb_ 14 5 4 6 6 6 6 */
/*:ref: snrm2_ 6 3 4 6 4 */
/*:ref: scbfix_ 14 16 214 214 6 4 6 4 4 6 6 4 6 6 6 6 6 4 */
/*:ref: sdot_ 6 5 4 6 4 6 4 */
/*:ref: msstop_ 4 16 4 4 4 6 6 4 6 6 6 4 6 6 6 6 6 4 */
/*:ref: saxpy_ 14 6 4 6 6 4 6 4 */
/*:ref: sonest_ 14 10 4 6 6 6 6 4 4 6 6 6 */
/*:ref: sdpchb_ 14 9 4 6 6 4 6 6 6 6 4 */
VEXTERNC int sratqr_(integer *n, real *eps1, real *d__, real *e, real *e2, integer *m, real *w, integer *ind, real *bd, logical *type__, integer *idef, integer *ierr);
/*:ref: r1mach_ 6 1 4 */

#endif /* _VCGCODE_H_ */

