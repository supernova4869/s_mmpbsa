/**
 *  @file       vpmg.h
 *  @brief      The primary header for PMG.
 *  @note       We provide this header whether or not we provide the BLAS
 *              library itself.  This gives us some compile-time type-checking
 *              even for architecture-dependent assembly-coded BLAS.
 *  @version    $Id: vpmg.h,v 1.7 2010/08/12 05:52:34 fetk Exp $
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

#ifndef _VPMG_H_
#define _VPMG_H_

#include <punc/punc_base.h>

#include <punc/vf2c.h>

/*
 * ***************************************************************************
 * Library VPMG prototypes
 * ***************************************************************************
 */

/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int builda_(integer *nx, integer *ny, integer *nz, integer *ipkey, integer *mgdisc, integer *numdia, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *xf, doublereal *yf, doublereal *zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal *fcf);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int builda_fv__(integer *nx, integer *ny, integer *nz, integer *ipkey, integer *numdia, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *fc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *xf, doublereal *yf, doublereal *zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal *fcf);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int builda_fe__(integer *nx, integer *ny, integer *nz, integer *ipkey, integer *numdia, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *fc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *one, doublereal *onw, doublereal *ue, doublereal *uw, doublereal *un, doublereal *us, doublereal *une, doublereal *unw, doublereal *use, doublereal *usw, doublereal *xf, doublereal *yf, doublereal *zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal *fcf);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int buildband_(integer *key, integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, integer *ipcb, doublereal *rpcb, doublereal *acb);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int buildband1_7__(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *oe, doublereal *on, doublereal *uc, integer *ipcb, doublereal *rpcb, doublereal *acb, integer *n, integer *m, integer *lda);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int buildband1_27__(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *one, doublereal *onw, doublereal *ue, doublereal *uw, doublereal *un, doublereal *us, doublereal *une, doublereal *unw, doublereal *use, doublereal *usw, integer *ipcb, doublereal *rpcb, doublereal *acb, integer *n, integer *m, integer *lda);
/** @brief  Library VPMG prototypes *
 *  @note:ref: dpbfa_ 14 5 7 4 4 4 4
    @author Michael Holst  */
VEXTERNC int buildg_(integer *nxf, integer *nyf, integer *nzf, integer *nxc, integer *nyc, integer *nzc, integer *numdia, doublereal *pcff, doublereal *acff, doublereal *ac);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int buildg_1__(integer *nxf, integer *nyf, integer *nzf, integer *nx, integer *ny, integer *nz, doublereal *opc, doublereal *opn, doublereal *ops, doublereal *ope, doublereal *opw, doublereal *opne, doublereal *opnw, doublereal *opse, doublereal *opsw, doublereal *upc, doublereal *upn, doublereal *ups, doublereal *upe, doublereal *upw, doublereal *upne, doublereal *upnw, doublereal *upse, doublereal *upsw, doublereal *dpc, doublereal *dpn, doublereal *dps, doublereal *dpe, doublereal *dpw, doublereal *dpne, doublereal *dpnw, doublereal *dpse, doublereal *dpsw, doublereal *oc, doublereal *xoc, doublereal *xoe, doublereal *xon, doublereal *xuc, doublereal *xone, doublereal *xonw, doublereal *xue, doublereal *xuw, doublereal *xun, doublereal *xus, doublereal *xune, doublereal *xunw, doublereal *xuse, doublereal *xusw);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int buildg_7__(integer *nxf, integer *nyf, integer *nzf, integer *nx, integer *ny, integer *nz, doublereal *opc, doublereal *opn, doublereal *ops, doublereal *ope, doublereal *opw, doublereal *opne, doublereal *opnw, doublereal *opse, doublereal *opsw, doublereal *upc, doublereal *upn, doublereal *ups, doublereal *upe, doublereal *upw, doublereal *upne, doublereal *upnw, doublereal *upse, doublereal *upsw, doublereal *dpc, doublereal *dpn, doublereal *dps, doublereal *dpe, doublereal *dpw, doublereal *dpne, doublereal *dpnw, doublereal *dpse, doublereal *dpsw, doublereal *oc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *xoc, doublereal *xoe, doublereal *xon, doublereal *xuc, doublereal *xone, doublereal *xonw, doublereal *xue, doublereal *xuw, doublereal *xun, doublereal *xus, doublereal *xune, doublereal *xunw, doublereal *xuse, doublereal *xusw);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int buildg_27__(integer *nxf, integer *nyf, integer *nzf, integer *nx, integer *ny, integer *nz, doublereal *opc, doublereal *opn, doublereal *ops, doublereal *ope, doublereal *opw, doublereal *opne, doublereal *opnw, doublereal *opse, doublereal *opsw, doublereal *upc, doublereal *upn, doublereal *ups, doublereal *upe, doublereal *upw, doublereal *upne, doublereal *upnw, doublereal *upse, doublereal *upsw, doublereal *dpc, doublereal *dpn, doublereal *dps, doublereal *dpe, doublereal *dpw, doublereal *dpne, doublereal *dpnw, doublereal *dpse, doublereal *dpsw, doublereal *oc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *one, doublereal *onw, doublereal *ue, doublereal *uw, doublereal *un, doublereal *us, doublereal *une, doublereal *unw, doublereal *use, doublereal *usw, doublereal *xoc, doublereal *xoe, doublereal *xon, doublereal *xuc, doublereal *xone, doublereal *xonw, doublereal *xue, doublereal *xuw, doublereal *xun, doublereal *xus, doublereal *xune, doublereal *xunw, doublereal *xuse, doublereal *xusw);
/** @brief  Library VPMG prototypes
    @author Michael Holst */
VEXTERNC int buildp_(integer *nxf, integer *nyf, integer *nzf, integer *nxc, integer *nyc, integer *nzc, integer *mgprol, integer *ipc, doublereal *rpc, doublereal *pc, doublereal *ac, doublereal *xf, doublereal *yf, doublereal *zf);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int buildp_trilin__(integer *nxf, integer *nyf, integer *nzf, integer *nxc, integer *nyc, integer *nzc, doublereal *pc, doublereal *xf, doublereal *yf, doublereal *zf);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int buildpb_trilin__(integer *nxf, integer *nyf, integer *nzf, integer *nxc, integer *nyc, integer *nzc, doublereal *opc, doublereal *opn, doublereal *ops, doublereal *ope, doublereal *opw, doublereal *opne, doublereal *opnw, doublereal *opse, doublereal *opsw, doublereal *upc, doublereal *upn, doublereal *ups, doublereal *upe, doublereal *upw, doublereal *upne, doublereal *upnw, doublereal *upse, doublereal *upsw, doublereal *dpc, doublereal *dpn, doublereal *dps, doublereal *dpe, doublereal *dpw, doublereal *dpne, doublereal *dpnw, doublereal *dpse, doublereal *dpsw, doublereal *xf, doublereal *yf, doublereal *zf);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int buildp_op7__(integer *nxf, integer *nyf, integer *nzf, integer *nxc, integer *nyc, integer *nzc, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *pc);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int buildpb_op7__(integer *nxf, integer *nyf, integer *nzf, integer *nxc, integer *nyc, integer *nzc, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *opc, doublereal *opn, doublereal *ops, doublereal *ope, doublereal *opw, doublereal *opne, doublereal *opnw, doublereal *opse, doublereal *opsw, doublereal *upc, doublereal *upn, doublereal *ups, doublereal *upe, doublereal *upw, doublereal *upne, doublereal *upnw, doublereal *upse, doublereal *upsw, doublereal *dpc, doublereal *dpn, doublereal *dps, doublereal *dpe, doublereal *dpw, doublereal *dpne, doublereal *dpnw, doublereal *dpse, doublereal *dpsw);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int buildp_op27__(integer *nxf, integer *nyf, integer *nzf, integer *nxc, integer *nyc, integer *nzc, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *pc);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int buildpb_op27__(integer *nxf, integer *nyf, integer *nzf, integer *nxc, integer *nyc, integer *nzc, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *one, doublereal *onw, doublereal *ue, doublereal *uw, doublereal *un, doublereal *us, doublereal *une, doublereal *unw, doublereal *use, doublereal *usw, doublereal *opc, doublereal *opn, doublereal *ops, doublereal *ope, doublereal *opw, doublereal *opne, doublereal *opnw, doublereal *opse, doublereal *opsw, doublereal *upc, doublereal *upn, doublereal *ups, doublereal *upe, doublereal *upw, doublereal *upne, doublereal *upnw, doublereal *upse, doublereal *upsw, doublereal *dpc, doublereal *dpn, doublereal *dps, doublereal *dpe, doublereal *dpw, doublereal *dpne, doublereal *dpnw, doublereal *dpse, doublereal *dpsw);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int buildp_modop7__(integer *nxf, integer *nyf, integer *nzf, integer *nxc, integer *nyc, integer *nzc, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *pc);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int buildpb_modop7__(integer *nxf, integer *nyf, integer *nzf, integer *nxc, integer *nyc, integer *nzc, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *opc, doublereal *opn, doublereal *ops, doublereal *ope, doublereal *opw, doublereal *opne, doublereal *opnw, doublereal *opse, doublereal *opsw, doublereal *upc, doublereal *upn, doublereal *ups, doublereal *upe, doublereal *upw, doublereal *upne, doublereal *upnw, doublereal *upse, doublereal *upsw, doublereal *dpc, doublereal *dpn, doublereal *dps, doublereal *dpe, doublereal *dpw, doublereal *dpne, doublereal *dpnw, doublereal *dpse, doublereal *dpsw);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int buildp_modop27__(integer *nxf, integer *nyf, integer *nzf, integer *nxc, integer *nyc, integer *nzc, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *pc);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int buildpb_modop27__(integer *nxf, integer *nyf, integer *nzf, integer *nxc, integer *nyc, integer *nzc, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *one, doublereal *onw, doublereal *ue, doublereal *uw, doublereal *un, doublereal *us, doublereal *une, doublereal *unw, doublereal *use, doublereal *usw, doublereal *opc, doublereal *opn, doublereal *ops, doublereal *ope, doublereal *opw, doublereal *opne, doublereal *opnw, doublereal *opse, doublereal *opsw, doublereal *upc, doublereal *upn, doublereal *ups, doublereal *upe, doublereal *upw, doublereal *upne, doublereal *upnw, doublereal *upse, doublereal *upsw, doublereal *dpc, doublereal *dpn, doublereal *dps, doublereal *dpe, doublereal *dpw, doublereal *dpne, doublereal *dpnw, doublereal *dpse, doublereal *dpsw);
/** @brief  Library VPMG prototypes \n
 *:ref: mresid_ 14 10 4 4 4 4 7 7 7 7 7 7 \n
 *:ref: xnrm2_ 7 4 4 4 4 7 \n
 *:ref: xdot_ 7 5 4 4 4 7 7 \n
 *:ref: xcopy_ 14 5 4 4 4 7 7 \n
 *:ref: xaxpy_ 14 6 4 4 4 7 7 7 \n
 *:ref: xscal_ 14 5 4 4 4 7 7 \n
 *:ref: matvec_ 14 9 4 4 4 4 7 7 7 7 7 
    @author Michael Holst */
VEXTERNC int cghs_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *x, doublereal *p, doublereal *ap, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int fcgmg_(integer *nx, integer *ny, integer *nz, doublereal *x, integer *iz, doublereal *w0, doublereal *w1, doublereal *w2, doublereal *w3, integer *istop, integer *itmax, integer *iters, integer *ierror, integer *nlev, integer *ilev, integer *nlev_real__, integer *mgsolv, integer *iok, integer *iinfo, doublereal *epsiln, doublereal *errtol, doublereal *omega, integer *nu1, integer *nu2, integer *mgsmoo, doublereal *w4, doublereal *w5, doublereal *w6, integer *ipc, doublereal *rpc, doublereal *pc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *tru);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int cgmg_(integer *nx, integer *ny, integer *nz, doublereal *x, integer *iz, doublereal *w0, doublereal *w1, doublereal *w2, doublereal *w3, integer *istop, integer *itmax, integer *iters, integer *ierror, integer *nlev, integer *ilev, integer *nlev_real__, integer *mgsolv, integer *iok, integer *iinfo, doublereal *epsiln, doublereal *errtol, doublereal *omega, integer *nu1, integer *nu2, integer *mgsmoo, doublereal *rr, doublereal *zz, doublereal *pp, integer *ipc, doublereal *rpc, doublereal *pc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *tru);
/** @brief  Library VPMG prototypes 
 *:ref: mkcors_ 14 7 4 4 4 4 4 4 4 \n
 *:ref: mkfine_ 14 7 4 4 4 4 4 4 4 \n
 *:ref: interp_ 14 9 4 4 4 4 4 4 7 7 7 \n
 *:ref: prtini_ 14 1 4 \n
 *:ref: prtstp_ 14 5 4 4 7 7 7 \n
 *:ref: xnrm1_ 7 4 4 4 4 7 \n
 *:ref: xnrm2_ 7 4 4 4 4 7 \n
 *:ref: matvec_ 14 9 4 4 4 4 7 7 7 7 7 \n
 *:ref: xdot_ 7 5 4 4 4 7 7 \n
 *:ref: mresid_ 14 10 4 4 4 4 7 7 7 7 7 7 \n
 *:ref: azeros_ 14 4 4 4 4 7 *\n
 *:ref: mvcs_ 14 32 4 4 4 7 4 7 7 7 7 4 4 4 4 4 4 4 4 4 4 7 7 7 4 4 4 4 7 7 7 7 7 7 \n
 *:ref: xcopy_ 14 5 4 4 4 7 7 \n
 *:ref: xaxpy_ 14 6 4 4 4 7 7 7 \n
 *:ref: xscal_ 14 5 4 4 4 7 7 \n
 *:ref: restrc_ 14 9 4 4 4 4 4 4 7 7 7 
    @author Michael Holst */
VEXTERNC int getpre_(integer *nx, integer *ny, integer *nz, integer *iz, integer *lev, integer *nlev_real__, doublereal *r__, doublereal *pc);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int cgmgdriv_(integer *iparm, doublereal *rparm, integer *iwork, doublereal *rwork, doublereal *u, doublereal *xf, doublereal *yf, doublereal *zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal *fcf, doublereal *tcf);
/** @brief  Library VPMG prototypes \n
 *:ref: maxlev_ 4 3 4 4 4 \n
 *:ref: mgsz_ 14 19 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 \n
 *:ref: prtstp_ 14 5 4 4 7 7 7 \n
 *:ref: buildstr_ 14 5 4 4 4 4 4 \n
 *:ref: tstart_ 14 2 7 7 \n
 *:ref: buildops_ 14 30 4 4 4 4 4 4 4 4 4 4 4 4 4 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 \n
 *:ref: tstop_ 14 3 7 7 7 \n
 *:ref: buildalg_ 14 14 4 4 4 4 4 4 4 7 7 7 7 7 7 7 \n
 *:ref: epsmac_ 7 1 4 \n
 *:ref: fbound00_ 14 4 4 4 4 7 \n
 *:ref: cgmg_ 14 35 4 4 4 7 4 7 7 7 7 4 4 4 4 4 4 4 4 4 4 7 7 7 4 4 4 7 7 7 4 7 7 7 7 7 7 \n
 *:ref: fcgmg_ 14 35 4 4 4 7 4 7 7 7 7 4 4 4 4 4 4 4 4 4 4 7 7 7 4 4 4 7 7 7 4 7 7 7 7 7 7 \n
 *:ref: fbound_ 14 8 4 4 4 4 7 7 7 7 
    @author Michael Holst */
VEXTERNC int cgmgdriv2_(integer *iparm, doublereal *rparm, integer *nx, integer *ny, integer *nz, doublereal *u, integer *iz, doublereal *w1, doublereal *w2, integer *ipc, doublereal *rpc, doublereal *pc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *xf, doublereal *yf, doublereal *zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal *fcf, doublereal *tcf);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int delget_(integer *nx, integer *ny, integer *nz, doublereal *xmin, doublereal *xmax, doublereal *ymin, doublereal *ymax, doublereal *zmin, doublereal *zmax, doublereal *epsin, doublereal *epsout, doublereal *rionst, doublereal *temper, integer *ncrgpt, integer *iepsmap, integer *idebmap, integer *icrgpos, doublereal *crg, doublereal *phi);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int gsrb_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int gsrb7_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *fc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int gsrb27_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *fc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *one, doublereal *onw, doublereal *ue, doublereal *uw, doublereal *un, doublereal *us, doublereal *une, doublereal *unw, doublereal *use, doublereal *usw, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes \n
 *:ref: mresid7_1s__ 14 13 4 4 4 4 7 7 7 7 7 7 7 7 7 \n
 *:ref: mresid27_1s__ 14 23 4 4 4 4 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 
    @author Michael Holst */
VEXTERNC int gsrb7x_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *fc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int readit_(integer *iparm, doublereal *rparm, integer *nx, integer *ny, integer *nz, integer *nlev, integer *nrwk, integer *niwk, integer *key, integer *meth);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int packmg_(integer *iparm, doublereal *rparm, integer *nrwk, integer *niwk, integer *nx, integer *ny, integer *nz, integer *nlev, integer *nu1, integer *nu2, integer *mgkey, integer *itmax, integer *istop, integer *ipcon, integer *nonlin, integer *mgsmoo, integer *mgprol, integer *mgcoar, integer *mgsolv, integer *mgdisc, integer *iinfo, doublereal *errtol, integer *ipkey, doublereal *omegal, doublereal *omegan, integer *irite, integer *iperf);
/** @brief  Library VPMG prototypes \n
 *:ref: fillco_ 14 17 4 7 4 4 4 7 7 7 7 7 7 7 7 7 7 7 7 
    @author Michael Holst */
VEXTERNC int writit_(integer *iparm, doublereal *rparm, integer *nx, integer *ny, integer *nz, doublereal *u, doublereal *xf, doublereal *yf, doublereal *zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal *fcf, doublereal *tcf, integer *key);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int matvec_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *x, doublereal *y);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int matvec7_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *x, doublereal *y);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int matvec7_1s__(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *x, doublereal *y);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int matvec27_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *x, doublereal *y);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int matvec27_1s__(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *one, doublereal *onw, doublereal *ue, doublereal *uw, doublereal *un, doublereal *us, doublereal *une, doublereal *unw, doublereal *use, doublereal *usw, doublereal *x, doublereal *y);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int mresid_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *x, doublereal *r__);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int mresid7_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *x, doublereal *r__);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int mresid7_1s__(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *fc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *x, doublereal *r__);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int mresid27_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *x, doublereal *r__);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int mresid27_1s__(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *fc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *one, doublereal *onw, doublereal *ue, doublereal *uw, doublereal *un, doublereal *us, doublereal *une, doublereal *unw, doublereal *use, doublereal *usw, doublereal *x, doublereal *r__);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nmatvec_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *x, doublereal *y, doublereal *w1);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nmatvec7_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *x, doublereal *y, doublereal *w1);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nmatvecd7_1s__(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *x, doublereal *y, doublereal *w1);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nmatvec27_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *x, doublereal *y, doublereal *w1);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nmatvecd27_1s__(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *one, doublereal *onw, doublereal *ue, doublereal *uw, doublereal *un, doublereal *us, doublereal *une, doublereal *unw, doublereal *use, doublereal *usw, doublereal *x, doublereal *y, doublereal *w1);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nmresid_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *x, doublereal *r__, doublereal *w1);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nmresid7_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *x, doublereal *r__, doublereal *w1);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nmresid7_1s__(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *fc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *x, doublereal *r__, doublereal *w1);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nmresid27_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *x, doublereal *r__, doublereal *w1);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nmresid27_1s__(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *fc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *one, doublereal *onw, doublereal *ue, doublereal *uw, doublereal *un, doublereal *us, doublereal *une, doublereal *unw, doublereal *use, doublereal *usw, doublereal *x, doublereal *r__, doublereal *w1);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int interp_(integer *nxc, integer *nyc, integer *nzc, integer *nxf, integer *nyf, integer *nzf, doublereal *xin, doublereal *xout, doublereal *pc);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int interp2_(integer *nxc, integer *nyc, integer *nzc, integer *nxf, integer *nyf, integer *nzf, doublereal *xin, doublereal *xout, doublereal *opc, doublereal *opn, doublereal *ops, doublereal *ope, doublereal *opw, doublereal *opne, doublereal *opnw, doublereal *opse, doublereal *opsw, doublereal *upc, doublereal *upn, doublereal *ups, doublereal *upe, doublereal *upw, doublereal *upne, doublereal *upnw, doublereal *upse, doublereal *upsw, doublereal *dpc, doublereal *dpn, doublereal *dps, doublereal *dpe, doublereal *dpw, doublereal *dpne, doublereal *dpnw, doublereal *dpse, doublereal *dpsw);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int restrc_(integer *nxf, integer *nyf, integer *nzf, integer *nxc, integer *nyc, integer *nzc, doublereal *xin, doublereal *xout, doublereal *pc);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int restrc2_(integer *nxf, integer *nyf, integer *nzf, integer *nxc, integer *nyc, integer *nzc, doublereal *xin, doublereal *xout, doublereal *opc, doublereal *opn, doublereal *ops, doublereal *ope, doublereal *opw, doublereal *opne, doublereal *opnw, doublereal *opse, doublereal *opsw, doublereal *upc, doublereal *upn, doublereal *ups, doublereal *upe, doublereal *upw, doublereal *upne, doublereal *upnw, doublereal *upse, doublereal *upsw, doublereal *dpc, doublereal *dpn, doublereal *dps, doublereal *dpe, doublereal *dpw, doublereal *dpne, doublereal *dpnw, doublereal *dpse, doublereal *dpsw);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int extrac_(integer *nxf, integer *nyf, integer *nzf, integer *nxc, integer *nyc, integer *nzc, doublereal *xin, doublereal *xout);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int interpz_(integer *nxc, integer *nyc, integer *nzc, integer *nxf, integer *nyf, integer *nzf, doublereal *xin, doublereal *xout);
/** @brief  Library VPMG prototypes \n
 *:ref: c_vec__ 14 7 7 7 7 4 4 4 4 \n
 *:ref: fbound00_ 14 4 4 4 4 7 
    @author Michael Holst */
VEXTERNC int restrcz_(integer *nxf, integer *nyf, integer *nzf, integer *nxc, integer *nyc, integer *nzc, doublereal *xin, doublereal *xout);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int daxpy_(integer *n, doublereal *alpha, doublereal *x, integer *istep, doublereal *y, integer *jstep);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int dcopy_(integer *n, doublereal *x, integer *istep, doublereal *y, integer *jstep);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC doublereal dasum_(integer *n, doublereal *x, integer *istep);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC doublereal dnrm1_(integer *n, doublereal *x, integer *istep);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC doublereal dnrm2_(integer *n, doublereal *x, integer *istep);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC doublereal dnrm8_(integer *n, doublereal *x, integer *istep);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int dscal_(integer *n, doublereal *fac, doublereal *x, integer *istep);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC doublereal ddot_(integer *n, doublereal *x, integer *istep, doublereal *y, integer *jstep);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC integer idamax_(integer *n, doublereal *sx, integer *incx);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int fmvcs_(integer *nx, integer *ny, integer *nz, doublereal *x, integer *iz, doublereal *w0, doublereal *w1, doublereal *w2, doublereal *w3, integer *istop, integer *itmax, integer *iters, integer *ierror, integer *nlev, integer *ilev, integer *nlev_real__, integer *mgsolv, integer *iok, integer *iinfo, doublereal *epsiln, doublereal *errtol, doublereal *omega, integer *nu1, integer *nu2, integer *mgsmoo, integer *ipc, doublereal *rpc, doublereal *pc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *tru);
/** @brief  Library VPMG prototypes \n
 *:ref: mkcors_ 14 7 4 4 4 4 4 4 4 \n
 *:ref: mkfine_ 14 7 4 4 4 4 4 4 4 \n
 *:ref: interp_ 14 9 4 4 4 4 4 4 7 7 7 \n
 *:ref: prtini_ 14 1 4 \n
 *:ref: prtstp_ 14 5 4 4 7 7 7 \n
 *:ref: xnrm1_ 7 4 4 4 4 7 \n
 *:ref: xnrm2_ 7 4 4 4 4 7 \n
 *:ref: matvec_ 14 9 4 4 4 4 7 7 7 7 7 \n
 *:ref: xdot_ 7 5 4 4 4 7 7 \n
 *:ref: azeros_ 14 4 4 4 4 7 \n
 *:ref: smooth_ 14 19 4 4 4 4 7 7 7 7 7 7 7 7 4 4 7 7 4 4 4 \n
 *:ref: xcopy_small__ 14 5 4 4 4 7 7 \n
 *:ref: dpbsl_ 14 5 7 4 4 4 7 \n
 *:ref: xcopy_large__ 14 5 4 4 4 7 7 \n
 *:ref: fbound00_ 14 4 4 4 4 7 \n
 *:ref: mresid_ 14 10 4 4 4 4 7 7 7 7 7 7 \n
 *:ref: xcopy_ 14 5 4 4 4 7 7 \n
 *:ref: xaxpy_ 14 6 4 4 4 7 7 7 \n
 *:ref: ivariv_ 4 2 4 4 \n
 *:ref: restrc_ 14 9 4 4 4 4 4 4 7 7 7 
    @author Michael Holst */
VEXTERNC int mvcs_(integer *nx, integer *ny, integer *nz, doublereal *x, integer *iz, doublereal *w0, doublereal *w1, doublereal *w2, doublereal *w3, integer *istop, integer *itmax, integer *iters, integer *ierror, integer *nlev, integer *ilev, integer *nlev_real__, integer *mgsolv, integer *iok, integer *iinfo, doublereal *epsiln, doublereal *errtol, doublereal *omega, integer *nu1, integer *nu2, integer *mgsmoo, integer *ipc, doublereal *rpc, doublereal *pc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *tru);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int mgdriv_(integer *iparm, doublereal *rparm, integer *iwork, doublereal *rwork, doublereal *u, doublereal *xf, doublereal *yf, doublereal *zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal *fcf, doublereal *tcf);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int mgdriv2_(integer *iparm, doublereal *rparm, integer *nx, integer *ny, integer *nz, doublereal *u, integer *iz, integer *ipc, doublereal *rpc, doublereal *pc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *xf, doublereal *yf, doublereal *zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal *fcf, doublereal *tcf);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int mgsz_(integer *mgcoar, integer *mgdisc, integer *mgsolv, integer *nx, integer *ny, integer *nz, integer *nlev, integer *nxc, integer *nyc, integer *nzc, integer *nf, integer *nc, integer *narr, integer *narrc, integer *n_rpc__, integer *n_iz__, integer *n_ipc__, integer *iretot, integer *iintot);
/** @brief  Library VPMG prototypes \n
 *:ref: maxlev_ 4 3 4 4 4 \n
 *:ref: prtstp_ 14 5 4 4 7 7 7 \n
 *:ref: buildstr_ 14 5 4 4 4 4 4 \n
 *:ref: tstart_ 14 2 7 7 \n
 *:ref: buildops_ 14 30 4 4 4 4 4 4 4 4 4 4 4 4 4 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 \n
 *:ref: tstop_ 14 3 7 7 7 \n
 *:ref: epsmac_ 7 1 4 \n
 *:ref: mkcors_ 14 7 4 4 4 4 4 4 4 \n
 *:ref: power_ 14 19 4 4 4 4 4 4 7 7 7 7 7 7 7 7 7 7 4 4 4 \n
 *:ref: azeros_ 14 4 4 4 4 7 \n
 *:ref: ipower_ 14 33 4 4 4 7 4 7 7 7 7 7 7 7 7 4 4 4 4 4 4 4 4 7 7 7 4 4 4 4 7 7 7 7 7 \n
 *:ref: mpower_ 14 33 4 4 4 7 4 7 7 7 7 7 7 7 4 4 4 4 4 4 4 4 7 7 7 4 4 4 4 7 7 7 7 7 7 \n
 *:ref: buildalg_ 14 14 4 4 4 4 4 4 4 7 7 7 7 7 7 7 \n
 *:ref: fbound00_ 14 4 4 4 4 7 \n
 *:ref: mvcs_ 14 32 4 4 4 7 4 7 7 7 7 4 4 4 4 4 4 4 4 4 4 7 7 7 4 4 4 4 7 7 7 7 7 7 \n
 *:ref: fmvcs_ 14 32 4 4 4 7 4 7 7 7 7 4 4 4 4 4 4 4 4 4 4 7 7 7 4 4 4 4 7 7 7 7 7 7 \n
 *:ref: mvfas_ 14 33 4 4 4 7 4 7 7 7 7 7 4 4 4 4 4 4 4 4 4 4 7 7 7 4 4 4 4 7 7 7 7 7 7 \n
 *:ref: fmvfas_ 14 33 4 4 4 7 4 7 7 7 7 7 4 4 4 4 4 4 4 4 4 4 7 7 7 4 4 4 4 7 7 7 7 7 7 \n
 *:ref: fbound_ 14 8 4 4 4 4 7 7 7 7 
    @author Michael Holst */
VEXTERNC int mgsize_(integer *mgcoar, integer *mgdisc, integer *mgsolv, integer *nx, integer *ny, integer *nz, integer *nlev);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int fmvfas_(integer *nx, integer *ny, integer *nz, doublereal *x, integer *iz, doublereal *w0, doublereal *w1, doublereal *w2, doublereal *w3, doublereal *w4, integer *istop, integer *itmax, integer *iters, integer *ierror, integer *nlev, integer *ilev, integer *nlev_real__, integer *mgsolv, integer *iok, integer *iinfo, doublereal *epsiln, doublereal *errtol, doublereal *omega, integer *nu1, integer *nu2, integer *mgsmoo, integer *ipc, doublereal *rpc, doublereal *pc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *tru);
/** @brief  Library VPMG prototypes \n
 *:ref: mkcors_ 14 7 4 4 4 4 4 4 4 \n
 *:ref: mkfine_ 14 7 4 4 4 4 4 4 4 \n
 *:ref: interp_ 14 9 4 4 4 4 4 4 7 7 7 \n
 *:ref: prtini_ 14 1 4 \n
 *:ref: prtstp_ 14 5 4 4 7 7 7 \n
 *:ref: azeros_ 14 4 4 4 4 7 \n
 *:ref: nmresid_ 14 11 4 4 4 4 7 7 7 7 7 7 7 \n
 *:ref: xnrm1_ 7 4 4 4 4 7 \n
 *:ref: xnrm2_ 7 4 4 4 4 7 \n
 *:ref: nmatvec_ 14 10 4 4 4 4 7 7 7 7 7 7 \n
 *:ref: xdot_ 7 5 4 4 4 7 7 \n
 *:ref: nsmooth_ 14 19 4 4 4 4 7 7 7 7 7 7 7 7 4 4 7 7 4 4 4 \n
 *:ref: xcopy_ 14 5 4 4 4 7 7 \n
 *:ref: xaxpy_ 14 6 4 4 4 7 7 7 \n
 *:ref: ivariv_ 4 2 4 4 \n
 *:ref: restrc_ 14 9 4 4 4 4 4 4 7 7 7 \n
 *:ref: extrac_ 14 8 4 4 4 4 4 4 7 7 \n
 *:ref: linesearch_ 14 15 4 4 4 7 4 7 7 7 7 7 7 7 7 7 7 
    @author Michael Holst */
VEXTERNC int mvfas_(integer *nx, integer *ny, integer *nz, doublereal *x, integer *iz, doublereal *w0, doublereal *w1, doublereal *w2, doublereal *w3, doublereal *w4, integer *istop, integer *itmax, integer *iters, integer *ierror, integer *nlev, integer *ilev, integer *nlev_real__, integer *mgsolv, integer *iok, integer *iinfo, doublereal *epsiln, doublereal *errtol, doublereal *omega, integer *nu1, integer *nu2, integer *mgsmoo, integer *ipc, doublereal *rpc, doublereal *pc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *tru);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC integer maxlev_(integer *n1, integer *n2, integer *n3);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int mkcors_(integer *numlev, integer *nxold, integer *nyold, integer *nzold, integer *nxnew, integer *nynew, integer *nznew);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int corsr_(integer *nold, integer *nnew);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int mkfine_(integer *numlev, integer *nxold, integer *nyold, integer *nzold, integer *nxnew, integer *nynew, integer *nznew);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int finer_(integer *nold, integer *nnew);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC integer ivariv_(integer *nu, integer *level);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int prtini_(integer *istop);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int prtstp_(integer *iok, integer *iters, doublereal *rsnrm, doublereal *rsden, doublereal *orsnrm);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int buildstr_(integer *nx, integer *ny, integer *nz, integer *nlev, integer *iz);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int buildops_(integer *nx, integer *ny, integer *nz, integer *nlev, integer *ipkey, integer *iinfo, integer *ido, integer *iz, integer *mgprol, integer *mgcoar, integer *mgsolv, integer *mgdisc, integer *ipc, doublereal *rpc, doublereal *pc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *xf, doublereal *yf, doublereal *zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal *fcf, doublereal *tcf);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int buildcopy0_(integer *nx, integer *ny, integer *nz, integer *nxf, integer *nyf, integer *nzf, doublereal *xc, doublereal *yc, doublereal *zc, doublereal *gxc, doublereal *gyc, doublereal *gzc, doublereal *a1c, doublereal *a2c, doublereal *a3c, doublereal *cc, doublereal *fc, doublereal *tc, doublereal *xf, doublereal *yf, doublereal *zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal *fcf, doublereal *tcf);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int buildharm0_(integer *nx, integer *ny, integer *nz, integer *nxf, integer *nyf, integer *nzf, doublereal *xc, doublereal *yc, doublereal *zc, doublereal *gxc, doublereal *gyc, doublereal *gzc, doublereal *a1c, doublereal *a2c, doublereal *a3c, doublereal *cc, doublereal *fc, doublereal *tc, doublereal *xf, doublereal *yf, doublereal *zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal *fcf, doublereal *tcf);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int buildgaler0_(integer *nxf, integer *nyf, integer *nzf, integer *nxc, integer *nyc, integer *nzc, integer *ipkey, integer *numdia, doublereal *pcff, integer *ipcff, doublereal *rpcff, doublereal *acff, doublereal *ccff, doublereal *fcff, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int buildgaler1_(integer *nxf, integer *nyf, integer *nzf, integer *nxc, integer *nyc, integer *nzc, integer *numdia, doublereal *pcff, integer *ipcff, doublereal *rpcff, doublereal *ccff, integer *ipc, doublereal *rpc, doublereal *cc);
/** @brief  Library VPMG prototypes \n
 *:ref: tstart_ 14 2 7 7 \n
 *:ref: tstop_ 14 3 7 7 7 \n
 *:ref: builda_ 14 22 4 4 4 4 4 4 4 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 \n
 *:ref: prtmatd_ 14 6 4 4 4 4 7 7 \n
 *:ref: buildp_ 14 14 4 4 4 4 4 4 4 4 7 7 7 7 7 7 \n
 *:ref: extrac_ 14 8 4 4 4 4 4 4 7 7 \n
 *:ref: buildband_ 14 10 4 4 4 4 4 7 7 4 7 7 \n
 *:ref: buildg_ 14 10 4 4 4 4 4 4 4 7 7 7 \n
 *:ref: restrc_ 14 9 4 4 4 4 4 4 7 7 7 \n
 *:ref: nmatvec_ 14 11 4 4 4 4 7 7 7 7 7 7 7 \n
 *:ref: matvec_ 14 9 4 4 4 4 7 7 7 7 7 
    @author Michael Holst */
VEXTERNC int buildalg_(integer *nx, integer *ny, integer *nz, integer *mode, integer *nlev, integer *iz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *x, doublereal *y, doublereal *tmp);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC doublereal epsmac_(integer *idum);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int tstart_(doublereal *before, doublereal *overhd);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int tstop_(doublereal *before, doublereal *overhd, doublereal *cputme);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int xaxpy_(integer *nx, integer *ny, integer *nz, doublereal *alpha, doublereal *x, doublereal *y);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int xcopy_(integer *nx, integer *ny, integer *nz, doublereal *x, doublereal *y);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC doublereal xnrm1_(integer *nx, integer *ny, integer *nz, doublereal *x);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC doublereal xnrm2_(integer *nx, integer *ny, integer *nz, doublereal *x);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC doublereal xnrm8_(integer *nx, integer *ny, integer *nz, doublereal *x);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int xscal_(integer *nx, integer *ny, integer *nz, doublereal *fac, doublereal *x);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC doublereal xdot_(integer *nx, integer *ny, integer *nz, doublereal *x, doublereal *y);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC doublereal xdot3_(integer *nx, integer *ny, integer *nz, doublereal *x, doublereal *y);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int zeros_(integer *nx, integer *ny, integer *nz, doublereal *x);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int xrand_(integer *nx, integer *ny, integer *nz, doublereal *x);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int cinit_(integer *nx, integer *ny, integer *nz, doublereal *x, doublereal *value);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int azeros_(integer *nx, integer *ny, integer *nz, doublereal *x);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int axrand_(integer *nx, integer *ny, integer *nz, doublereal *x);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int xcopy_small__(integer *nx, integer *ny, integer *nz, doublereal *x, doublereal *y);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int xcopy_large__(integer *nx, integer *ny, integer *nz, doublereal *x, doublereal *y);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int fbound_(integer *ibound, integer *nx, integer *ny, integer *nz, doublereal *x, doublereal *gxc, doublereal *gyc, doublereal *gzc);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int fbound00_(integer *nx, integer *ny, integer *nz, doublereal *x);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int xprint_(integer *nx, integer *ny, integer *nz, doublereal *x);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int prtmatd_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int prtmatd7_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *oe, doublereal *on, doublereal *uc);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int prtmatd27_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *one, doublereal *onw, doublereal *ue, doublereal *uw, doublereal *un, doublereal *us, doublereal *une, doublereal *unw, doublereal *use, doublereal *usw);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int prtmatb_(doublereal *a, integer *n, integer *m, integer *lda);
/** @brief  Library VPMG prototypes \n
 *:ref: tsecnd_ 7 0 \n
 *:ref: rand_ 6 1 4 \n
 *:ref: matvec_ 14 9 4 4 4 4 7 7 7 7 7 \n
 *:ref: c_vec__ 14 7 7 7 7 4 4 4 4 \n
 *:ref: dc_vec__ 14 7 7 7 7 4 4 4 4 
    @author Michael Holst */
VEXTERNC int linesearch_(integer *nx, integer *ny, integer *nz, doublereal *alpha, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *p, doublereal *x, doublereal *r__, doublereal *ap, doublereal *zk, doublereal *zkp1);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int dpbco_(doublereal *abd, integer *lda, integer *n, integer *m, doublereal *rcond, doublereal *z__, integer *info);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int dpbfa_(doublereal *abd, integer *lda, integer *n, integer *m, integer *info);
/** @brief  Library VPMG prototypes \n
 *:ref: dasum_ 7 3 4 7 4 \n
 *:ref: dscal_ 14 4 4 7 7 4 \n
 *:ref: daxpy_ 14 6 4 7 7 4 7 4 \n
 *:ref: ddot_ 7 5 4 7 4 7 4 
    @author Michael Holst */
VEXTERNC int dpbsl_(doublereal *abd, integer *lda, integer *n, integer *m, doublereal *b);
/** @brief  Library VPMG prototypes \n
 *:ref: nmresid_ 14 11 4 4 4 4 7 7 7 7 7 7 7 \n
 *:ref: xdot_ 7 5 4 4 4 7 7 \n
 *:ref: xcopy_ 14 5 4 4 4 7 7 \n
 *:ref: xaxpy_ 14 6 4 4 4 7 7 7 \n
 *:ref: xscal_ 14 5 4 4 4 7 7 \n
 *:ref: linesearch_ 14 15 4 4 4 7 4 7 7 7 7 7 7 7 7 7 7 
    @author Michael Holst */
VEXTERNC int ncghs_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *x, doublereal *p, doublereal *ap, doublereal *r__, doublereal *zk, doublereal *zkp1, doublereal *tmp, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int ncghsdriv_(integer *iparm, doublereal *rparm, integer *iwork, doublereal *rwork, doublereal *u, doublereal *xf, doublereal *yf, doublereal *zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal *fcf, doublereal *tcf);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */VEXTERNC int ncghsdriv2_(integer *iparm, doublereal *rparm, integer *nx, integer *ny, integer *nz, doublereal *u, integer *iz, doublereal *w0, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *xf, doublereal *yf, doublereal *zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal *fcf, doublereal *tcf);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int ncghsgo_(integer *nx, integer *ny, integer *nz, doublereal *x, doublereal *r__, doublereal *p, doublereal *ap, doublereal *zk, doublereal *zkp1, doublereal *tmp, integer *istop, integer *itmax, integer *iters, integer *ierror, integer *iok, integer *iinfo, doublereal *epsiln, doublereal *errtol, doublereal *omega, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *tru);
/** @brief  Library VPMG prototypes \n
 *:ref: prtstp_ 14 5 4 4 7 7 7 \n
 *:ref: buildstr_ 14 5 4 4 4 4 4 \n
 *:ref: tstart_ 14 2 7 7 \n
 *:ref: buildops_ 14 30 4 4 4 4 4 4 4 4 4 4 4 4 4 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 \n
 *:ref: tstop_ 14 3 7 7 7 \n
 *:ref: nmatvec_ 14 10 4 4 4 4 7 7 7 7 7 7 \n
 *:ref: matvec_ 14 9 4 4 4 4 7 7 7 7 7 \n
 *:ref: epsmac_ 7 1 4 \n
 *:ref: fbound00_ 14 4 4 4 4 7 \n
 *:ref: fbound_ 14 8 4 4 4 4 7 7 7 7 \n
 *:ref: prtini_ 14 1 4 \n
 *:ref: azeros_ 14 4 4 4 4 7 \n
 *:ref: nmresid_ 14 11 4 4 4 4 7 7 7 7 7 7 7 \n
 *:ref: xnrm1_ 7 4 4 4 4 7 \n
 *:ref: xnrm2_ 7 4 4 4 4 7 \n
 *:ref: xdot_ 7 5 4 4 4 7 7 \n
 *:ref: xcopy_ 14 5 4 4 4 7 7 \n
 *:ref: xaxpy_ 14 6 4 4 4 7 7 7 \n
 *:ref: xscal_ 14 5 4 4 4 7 7 \n
 *:ref: linesearch_ 14 15 4 4 4 7 4 7 7 7 7 7 7 7 7 7 7 \n
 *:ref: mresid_ 14 10 4 4 4 4 7 7 7 7 7 7 
    @author Michael Holst */
VEXTERNC int cghsgo_(integer *nx, integer *ny, integer *nz, doublereal *x, doublereal *r__, doublereal *p, doublereal *ap, doublereal *zk, doublereal *zkp1, doublereal *tmp, integer *istop, integer *itmax, integer *iters, integer *ierror, integer *iok, integer *iinfo, doublereal *epsiln, doublereal *errtol, doublereal *omega, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *tru);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int newdriv_(integer *iparm, doublereal *rparm, integer *iwork, doublereal *rwork, doublereal *u, doublereal *xf, doublereal *yf, doublereal *zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal *fcf, doublereal *tcf);
/** @brief  Library VPMG prototypes \n
 *:ref: maxlev_ 4 3 4 4 4 \n
 *:ref: mgsz_ 14 19 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 \n
 *:ref: prtstp_ 14 5 4 4 7 7 7 \n
 *:ref: buildstr_ 14 5 4 4 4 4 4 \n
 *:ref: tstart_ 14 2 7 7 \n
 *:ref: buildops_ 14 30 4 4 4 4 4 4 4 4 4 4 4 4 4 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 \n
 *:ref: tstop_ 14 3 7 7 7 \n
 *:ref: buildalg_ 14 14 4 4 4 4 4 4 4 7 7 7 7 7 7 7 \n
 *:ref: epsmac_ 7 1 4 \n
 *:ref: fbound00_ 14 4 4 4 4 7 \n
 *:ref: newton_ 14 35 4 4 4 7 4 7 7 7 7 4 4 4 4 4 4 4 4 4 4 7 7 7 4 4 4 7 7 7 4 7 7 7 7 7 7 \n
 *:ref: fnewton_ 14 35 4 4 4 7 4 7 7 7 7 4 4 4 4 4 4 4 4 4 4 7 7 7 4 4 4 7 7 7 4 7 7 7 7 7 7 \n
 *:ref: fbound_ 14 8 4 4 4 4 7 7 7 7 
    @author Michael Holst */
VEXTERNC int newdriv2_(integer *iparm, doublereal *rparm, integer *nx, integer *ny, integer *nz, doublereal *u, integer *iz, doublereal *w1, doublereal *w2, integer *ipc, doublereal *rpc, doublereal *pc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *xf, doublereal *yf, doublereal *zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal *fcf, doublereal *tcf);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int fnewton_(integer *nx, integer *ny, integer *nz, doublereal *x, integer *iz, doublereal *w0, doublereal *w1, doublereal *w2, doublereal *w3, integer *istop, integer *itmax, integer *iters, integer *ierror, integer *nlev, integer *ilev, integer *nlev_real__, integer *mgsolv, integer *iok, integer *iinfo, doublereal *epsiln, doublereal *errtol, doublereal *omega, integer *nu1, integer *nu2, integer *mgsmoo, doublereal *cprime, doublereal *rhs, doublereal *xtmp, integer *ipc, doublereal *rpc, doublereal *pc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *tru);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int newton_(integer *nx, integer *ny, integer *nz, doublereal *x, integer *iz, doublereal *w0, doublereal *w1, doublereal *w2, doublereal *w3, integer *istop, integer *itmax, integer *iters, integer *ierror, integer *nlev, integer *ilev, integer *nlev_real__, integer *mgsolv, integer *iok, integer *iinfo, doublereal *epsiln, doublereal *errtol, doublereal *omega, integer *nu1, integer *nu2, integer *mgsmoo, doublereal *cprime, doublereal *rhs, doublereal *xtmp, integer *ipc, doublereal *rpc, doublereal *pc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *tru);
/** @brief  Library VPMG prototypes \n
 *:ref: mkcors_ 14 7 4 4 4 4 4 4 4 \n
 *:ref: mkfine_ 14 7 4 4 4 4 4 4 4 \n
 *:ref: interp_ 14 9 4 4 4 4 4 4 7 7 7 \n
 *:ref: prtini_ 14 1 4 \n
 *:ref: prtstp_ 14 5 4 4 7 7 7 \n
 *:ref: azeros_ 14 4 4 4 4 7 \n
 *:ref: nmresid_ 14 11 4 4 4 4 7 7 7 7 7 7 7 \n
 *:ref: xnrm1_ 7 4 4 4 4 7 \n
 *:ref: xnrm2_ 7 4 4 4 4 7 \n
 *:ref: nmatvec_ 14 10 4 4 4 4 7 7 7 7 7 7 \n
 *:ref: xdot_ 7 5 4 4 4 7 7 \n
 *:ref: xcopy_ 14 5 4 4 4 7 7 \n
 *:ref: mvcs_ 14 32 4 4 4 7 4 7 7 7 7 4 4 4 4 4 4 4 4 4 4 7 7 7 4 4 4 4 7 7 7 7 7 7 \n
 *:ref: xaxpy_ 14 6 4 4 4 7 7 7 \n
 *:ref: power_ 14 19 4 4 4 4 4 4 7 7 7 7 7 7 7 7 7 7 4 4 4 \n
 *:ref: ipower_ 14 33 4 4 4 7 4 7 7 7 7 7 7 7 7 4 4 4 4 4 4 4 4 7 7 7 4 4 4 4 7 7 7 7 7 \n
 *:ref: dc_vec__ 14 7 7 7 7 4 4 4 4 \n
 *:ref: restrc_ 14 9 4 4 4 4 4 4 7 7 7 
    @author Michael Holst */
VEXTERNC int getjac_(integer *nx, integer *ny, integer *nz, integer *nlev_real__, integer *iz, integer *lev, integer *ipkey, doublereal *x, doublereal *r__, doublereal *cprime, doublereal *rhs, doublereal *cc, doublereal *pc);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int ngsrb_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int ngsrb7_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *fc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int ngsrb27_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *fc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *one, doublereal *onw, doublereal *ue, doublereal *uw, doublereal *un, doublereal *us, doublereal *une, doublereal *unw, doublereal *use, doublereal *usw, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes \n
 *:ref: c_scal__ 7 3 7 7 4 \n
 *:ref: dc_scal__ 7 3 7 7 4 \n
 *:ref: nmresid7_1s__ 14 14 4 4 4 4 7 7 7 7 7 7 7 7 7 7 \n
 *:ref: nmresid27_1s__ 14 24 4 4 4 4 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 
    @author Michael Holst */
VEXTERNC int ngsrb7x_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *fc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int ngsrbdriv_(integer *iparm, doublereal *rparm, integer *iwork, doublereal *rwork, doublereal *u, doublereal *xf, doublereal *yf, doublereal *zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal *fcf, doublereal *tcf);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int ngsrbdriv2_(integer *iparm, doublereal *rparm, integer *nx, integer *ny, integer *nz, doublereal *u, integer *iz, doublereal *w0, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *xf, doublereal *yf, doublereal *zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal *fcf, doublereal *tcf);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int ngsrbgo_(integer *nx, integer *ny, integer *nz, doublereal *x, doublereal *r__, doublereal *w1, doublereal *w2, integer *istop, integer *itmax, integer *iters, integer *ierror, integer *iok, integer *iinfo, doublereal *epsiln, doublereal *errtol, doublereal *omega, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *tru);
/** @brief  Library VPMG prototypes \n
 *:ref: prtstp_ 14 5 4 4 7 7 7 \n
 *:ref: buildstr_ 14 5 4 4 4 4 4 \n
 *:ref: tstart_ 14 2 7 7 \n
 *:ref: buildops_ 14 30 4 4 4 4 4 4 4 4 4 4 4 4 4 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 \n
 *:ref: tstop_ 14 3 7 7 7 \n
 *:ref: nmatvec_ 14 10 4 4 4 4 7 7 7 7 7 7 \n
 *:ref: matvec_ 14 9 4 4 4 4 7 7 7 7 7 \n
 *:ref: epsmac_ 7 1 4 \n
 *:ref: fbound00_ 14 4 4 4 4 7 \n
 *:ref: fbound_ 14 8 4 4 4 4 7 7 7 7 \n
 *:ref: prtini_ 14 1 4 \n
 *:ref: azeros_ 14 4 4 4 4 7 \n
 *:ref: nmresid_ 14 11 4 4 4 4 7 7 7 7 7 7 7 \n
 *:ref: xnrm1_ 7 4 4 4 4 7 \n
 *:ref: xnrm2_ 7 4 4 4 4 7 \n
 *:ref: xdot_ 7 5 4 4 4 7 7 \n
 *:ref: xcopy_ 14 5 4 4 4 7 7 \n
 *:ref: ngsrb_ 14 18 4 4 4 4 7 7 7 7 7 7 7 7 4 4 7 7 4 4 \n
 *:ref: xaxpy_ 14 6 4 4 4 7 7 7 \n
 *:ref: gsrb_ 14 18 4 4 4 4 7 7 7 7 7 7 7 7 4 4 7 7 4 4 
    @author Michael Holst */
VEXTERNC int gsrbgo_(integer *nx, integer *ny, integer *nz, doublereal *x, doublereal *r__, doublereal *w1, doublereal *w2, integer *istop, integer *itmax, integer *iters, integer *ierror, integer *iok, integer *iinfo, doublereal *epsiln, doublereal *errtol, doublereal *omega, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *tru);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int ninterp_(integer *nxc, integer *nyc, integer *nzc, integer *nxf, integer *nyf, integer *nzf, doublereal *xin, doublereal *xout, doublereal *pc, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int ninterp7_(integer *nxc, integer *nyc, integer *nzc, integer *nxf, integer *nyf, integer *nzf, doublereal *xin, doublereal *xout, doublereal *opc, doublereal *opn, doublereal *ops, doublereal *ope, doublereal *opw, doublereal *opne, doublereal *opnw, doublereal *opse, doublereal *opsw, doublereal *upc, doublereal *upn, doublereal *ups, doublereal *upe, doublereal *upw, doublereal *upne, doublereal *upnw, doublereal *upse, doublereal *upsw, doublereal *dpc, doublereal *dpn, doublereal *dps, doublereal *dpe, doublereal *dpw, doublereal *dpne, doublereal *dpnw, doublereal *dpse, doublereal *dpsw, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *cc, doublereal *fc);
/** @brief  Library VPMG prototypes \n
 *:ref: fbound00_ 14 4 4 4 4 7 \n
 *:ref: c_scal__ 7 3 7 7 4 \n
 *:ref: dc_scal__ 7 3 7 7 4 
    @author Michael Holst */
VEXTERNC int ninterp27_(integer *nxc, integer *nyc, integer *nzc, integer *nxf, integer *nyf, integer *nzf, doublereal *xin, doublereal *xout, doublereal *opc, doublereal *opn, doublereal *ops, doublereal *ope, doublereal *opw, doublereal *opne, doublereal *opnw, doublereal *opse, doublereal *opsw, doublereal *upc, doublereal *upn, doublereal *ups, doublereal *upe, doublereal *upw, doublereal *upne, doublereal *upnw, doublereal *upse, doublereal *upsw, doublereal *dpc, doublereal *dpn, doublereal *dps, doublereal *dpe, doublereal *dpw, doublereal *dpne, doublereal *dpnw, doublereal *dpse, doublereal *dpsw, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *one, doublereal *onw, doublereal *ue, doublereal *uw, doublereal *un, doublereal *us, doublereal *une, doublereal *unw, doublereal *use, doublereal *usw, doublereal *cc, doublereal *fc);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nrich_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nrich7_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *fc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes \n
 *:ref: nmresid7_1s__ 14 14 4 4 4 4 7 7 7 7 7 7 7 7 7 7 \n
 *:ref: mresid27_1s__ 14 23 4 4 4 4 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 \n
 *:ref: nmresid27_1s__ 14 24 4 4 4 4 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 
    @author Michael Holst */
VEXTERNC int nrich27_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *fc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *one, doublereal *onw, doublereal *ue, doublereal *uw, doublereal *un, doublereal *us, doublereal *une, doublereal *unw, doublereal *use, doublereal *usw, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nrichdriv_(integer *iparm, doublereal *rparm, integer *iwork, doublereal *rwork, doublereal *u, doublereal *xf, doublereal *yf, doublereal *zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal *fcf, doublereal *tcf);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nrichdriv2_(integer *iparm, doublereal *rparm, integer *nx, integer *ny, integer *nz, doublereal *u, integer *iz, doublereal *w0, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *xf, doublereal *yf, doublereal *zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal *fcf, doublereal *tcf);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nrichgo_(integer *nx, integer *ny, integer *nz, doublereal *x, doublereal *r__, doublereal *w1, doublereal *w2, integer *istop, integer *itmax, integer *iters, integer *ierror, integer *iok, integer *iinfo, doublereal *epsiln, doublereal *errtol, doublereal *omega, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *tru);
/** @brief  Library VPMG prototypes \n
 *:ref: prtstp_ 14 5 4 4 7 7 7 \n
 *:ref: buildstr_ 14 5 4 4 4 4 4 \n
 *:ref: tstart_ 14 2 7 7 \n
 *:ref: buildops_ 14 30 4 4 4 4 4 4 4 4 4 4 4 4 4 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 \n
 *:ref: tstop_ 14 3 7 7 7 \n
 *:ref: nmatvec_ 14 10 4 4 4 4 7 7 7 7 7 7 \n
 *:ref: matvec_ 14 9 4 4 4 4 7 7 7 7 7 \n
 *:ref: epsmac_ 7 1 4 \n
 *:ref: fbound00_ 14 4 4 4 4 7 \n
 *:ref: fbound_ 14 8 4 4 4 4 7 7 7 7 \n
 *:ref: prtini_ 14 1 4 \n
 *:ref: azeros_ 14 4 4 4 4 7 \n
 *:ref: nmresid_ 14 11 4 4 4 4 7 7 7 7 7 7 7 \n
 *:ref: xnrm1_ 7 4 4 4 4 7 \n
 *:ref: xnrm2_ 7 4 4 4 4 7 \n
 *:ref: xdot_ 7 5 4 4 4 7 7 \n
 *:ref: xcopy_ 14 5 4 4 4 7 7 \n
 *:ref: nrich_ 14 18 4 4 4 4 7 7 7 7 7 7 7 7 4 4 7 7 4 4 \n
 *:ref: xaxpy_ 14 6 4 4 4 7 7 7 \n
 *:ref: rich_ 14 18 4 4 4 4 7 7 7 7 7 7 7 7 4 4 7 7 4 4 
    @author Michael Holst */
VEXTERNC int richgo_(integer *nx, integer *ny, integer *nz, doublereal *x, doublereal *r__, doublereal *w1, doublereal *w2, integer *istop, integer *itmax, integer *iters, integer *ierror, integer *iok, integer *iinfo, doublereal *epsiln, doublereal *errtol, doublereal *omega, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *tru);
/** @brief  Library VPMG prototypes \n
 *:ref: nwjac_ 14 18 4 4 4 4 7 7 7 7 7 7 7 7 4 4 7 7 4 4 \n
 *:ref: ngsrb_ 14 18 4 4 4 4 7 7 7 7 7 7 7 7 4 4 7 7 4 4 \n
 *:ref: nsor_ 14 18 4 4 4 4 7 7 7 7 7 7 7 7 4 4 7 7 4 4 \n
 *:ref: nrich_ 14 18 4 4 4 4 7 7 7 7 7 7 7 7 4 4 7 7 4 4 
    @author Michael Holst */
VEXTERNC int nsmooth_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint, integer *meth);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nsor_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nsor7_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *fc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes \n
/*:ref: c_scal__ 7 3 7 7 4 \n
/*:ref: dc_scal__ 7 3 7 7 4 \n
/*:ref: nmresid7_1s__ 14 14 4 4 4 4 7 7 7 7 7 7 7 7 7 7 \n
/*:ref: nmresid27_1s__ 14 24 4 4 4 4 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 
    @author Michael Holst */
VEXTERNC int nsor27_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *fc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *one, doublereal *onw, doublereal *ue, doublereal *uw, doublereal *un, doublereal *us, doublereal *une, doublereal *unw, doublereal *use, doublereal *usw, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nsordriv_(integer *iparm, doublereal *rparm, integer *iwork, doublereal *rwork, doublereal *u, doublereal *xf, doublereal *yf, doublereal *zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal *fcf, doublereal *tcf);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nsordriv2_(integer *iparm, doublereal *rparm, integer *nx, integer *ny, integer *nz, doublereal *u, integer *iz, doublereal *w0, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *xf, doublereal *yf, doublereal *zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal *fcf, doublereal *tcf);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nsorgo_(integer *nx, integer *ny, integer *nz, doublereal *x, doublereal *r__, doublereal *w1, doublereal *w2, integer *istop, integer *itmax, integer *iters, integer *ierror, integer *iok, integer *iinfo, doublereal *epsiln, doublereal *errtol, doublereal *omega, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *tru);
/** @brief  Library VPMG prototypes \n
 *:ref: prtstp_ 14 5 4 4 7 7 7 \n
 *:ref: buildstr_ 14 5 4 4 4 4 4 \n
 *:ref: tstart_ 14 2 7 7 \n
 *:ref: buildops_ 14 30 4 4 4 4 4 4 4 4 4 4 4 4 4 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 \n
 *:ref: tstop_ 14 3 7 7 7 \n
 *:ref: nmatvec_ 14 10 4 4 4 4 7 7 7 7 7 7 \n
 *:ref: matvec_ 14 9 4 4 4 4 7 7 7 7 7 \n
 *:ref: epsmac_ 7 1 4 \n
 *:ref: fbound00_ 14 4 4 4 4 7 \n
 *:ref: fbound_ 14 8 4 4 4 4 7 7 7 7 \n
 *:ref: prtini_ 14 1 4 \n
 *:ref: azeros_ 14 4 4 4 4 7 \n
 *:ref: nmresid_ 14 11 4 4 4 4 7 7 7 7 7 7 7 \n
 *:ref: xnrm1_ 7 4 4 4 4 7 \n
 *:ref: xnrm2_ 7 4 4 4 4 7 \n
 *:ref: xdot_ 7 5 4 4 4 7 7 \n
 *:ref: xcopy_ 14 5 4 4 4 7 7 \n
 *:ref: nsor_ 14 18 4 4 4 4 7 7 7 7 7 7 7 7 4 4 7 7 4 4 \n
 *:ref: xaxpy_ 14 6 4 4 4 7 7 7 \n
 *:ref: sor_ 14 18 4 4 4 4 7 7 7 7 7 7 7 7 4 4 7 7 4 4 
    @author Michael Holst */
VEXTERNC int sorgo_(integer *nx, integer *ny, integer *nz, doublereal *x, doublereal *r__, doublereal *w1, doublereal *w2, integer *istop, integer *itmax, integer *iters, integer *ierror, integer *iok, integer *iinfo, doublereal *epsiln, doublereal *errtol, doublereal *omega, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *tru);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nwjac_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nwjac7_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *fc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes \n
 *:ref: c_scal__ 7 3 7 7 4 \n
 *:ref: dc_scal__ 7 3 7 7 4 \n
 *:ref: nmresid7_1s__ 14 14 4 4 4 4 7 7 7 7 7 7 7 7 7 7 \n
 *:ref: nmresid27_1s__ 14 24 4 4 4 4 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 
    @author Michael Holst */
VEXTERNC int nwjac27_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *fc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *one, doublereal *onw, doublereal *ue, doublereal *uw, doublereal *un, doublereal *us, doublereal *une, doublereal *unw, doublereal *use, doublereal *usw, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nwjacdriv_(integer *iparm, doublereal *rparm, integer *iwork, doublereal *rwork, doublereal *u, doublereal *xf, doublereal *yf, doublereal *zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal *fcf, doublereal *tcf);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nwjacdriv2_(integer *iparm, doublereal *rparm, integer *nx, integer *ny, integer *nz, doublereal *u, integer *iz, doublereal *w0, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *xf, doublereal *yf, doublereal *zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal *fcf, doublereal *tcf);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int nwjacgo_(integer *nx, integer *ny, integer *nz, doublereal *x, doublereal *r__, doublereal *w1, doublereal *w2, integer *istop, integer *itmax, integer *iters, integer *ierror, integer *iok, integer *iinfo, doublereal *epsiln, doublereal *errtol, doublereal *omega, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *tru);
/** @brief  Library VPMG prototypes \n
 *:ref: prtstp_ 14 5 4 4 7 7 7 \n
 *:ref: buildstr_ 14 5 4 4 4 4 4 \n
 *:ref: tstart_ 14 2 7 7 \n
 *:ref: buildops_ 14 30 4 4 4 4 4 4 4 4 4 4 4 4 4 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 \n
 *:ref: tstop_ 14 3 7 7 7 \n
 *:ref: nmatvec_ 14 10 4 4 4 4 7 7 7 7 7 7 \n
 *:ref: matvec_ 14 9 4 4 4 4 7 7 7 7 7 \n
 *:ref: epsmac_ 7 1 4 \n
 *:ref: fbound00_ 14 4 4 4 4 7 \n
 *:ref: fbound_ 14 8 4 4 4 4 7 7 7 7 \n
 *:ref: prtini_ 14 1 4 \n
 *:ref: azeros_ 14 4 4 4 4 7 \n
 *:ref: nmresid_ 14 11 4 4 4 4 7 7 7 7 7 7 7 \n
 *:ref: xnrm1_ 7 4 4 4 4 7 \n
 *:ref: xnrm2_ 7 4 4 4 4 7 \n
 *:ref: xdot_ 7 5 4 4 4 7 7 \n
 *:ref: xcopy_ 14 5 4 4 4 7 7 \n
 *:ref: nwjac_ 14 18 4 4 4 4 7 7 7 7 7 7 7 7 4 4 7 7 4 4 \n
 *:ref: xaxpy_ 14 6 4 4 4 7 7 7 \n
 *:ref: wjac_ 14 18 4 4 4 4 7 7 7 7 7 7 7 7 4 4 7 7 4 4 
    @author Michael Holst */
VEXTERNC int wjacgo_(integer *nx, integer *ny, integer *nz, doublereal *x, doublereal *r__, doublereal *w1, doublereal *w2, integer *istop, integer *itmax, integer *iters, integer *ierror, integer *iok, integer *iinfo, doublereal *epsiln, doublereal *errtol, doublereal *omega, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *tru);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC doublereal c_scal__(doublereal *coef, doublereal *u, integer *ipkey);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC doublereal dc_scal__(doublereal *coef, doublereal *u, integer *ipkey);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int c_vec__(doublereal *coef, doublereal *uin, doublereal *uout, integer *nx, integer *ny, integer *nz, integer *ipkey);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int dc_vec__(doublereal *coef, doublereal *uin, doublereal *uout, integer *nx, integer *ny, integer *nz, integer *ipkey);
/** @brief  Library VPMG prototypes
    @author Michael Holst */
VEXTERNC int fillco_(integer *iparm, doublereal *rparm, integer *nx, integer *ny, integer *nz, doublereal *xf, doublereal *yf, doublereal *zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal *fcf, doublereal *tcf);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int power_(integer *nx, integer *ny, integer *nz, integer *iz, integer *ilev, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *w1, doublereal *w2, doublereal *w3, doublereal *w4, doublereal *eigmax, doublereal *eigmax_model__, doublereal *tol, integer *itmax, integer *iters, integer *iinfo);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int ipower_(integer *nx, integer *ny, integer *nz, doublereal *u, integer *iz, doublereal *w0, doublereal *w1, doublereal *w2, doublereal *w3, doublereal *w4, doublereal *eigmin, doublereal *eigmin_model__, doublereal *tol, integer *itmax, integer *iters, integer *nlev, integer *ilev, integer *nlev_real__, integer *mgsolv, integer *iok, integer *iinfo, doublereal *epsiln, doublereal *errtol, doublereal *omega, integer *nu1, integer *nu2, integer *mgsmoo, integer *ipc, doublereal *rpc, doublereal *pc, doublereal *ac, doublereal *cc, doublereal *tru);
/** @brief  Library VPMG prototypes \n
 *:ref: axrand_ 14 4 4 4 4 7 \n
 *:ref: azeros_ 14 4 4 4 4 7 \n
 *:ref: xnrm2_ 7 4 4 4 4 7 \n
 *:ref: xscal_ 14 5 4 4 4 7 7 \n
 *:ref: matvec_ 14 9 4 4 4 4 7 7 7 7 7 \n
 *:ref: xdot_ 7 5 4 4 4 7 7 \n
 *:ref: xcopy_ 14 5 4 4 4 7 7 \n
 *:ref: xaxpy_ 14 6 4 4 4 7 7 7 \n
 *:ref: mvcs_ 14 32 4 4 4 7 4 7 7 7 7 4 4 4 4 4 4 4 4 4 4 7 7 7 4 4 4 4 7 7 7 7 7 7 
    @author Michael Holst */
VEXTERNC int mpower_(integer *nx, integer *ny, integer *nz, doublereal *u, integer *iz, doublereal *w0, doublereal *w1, doublereal *w2, doublereal *w3, doublereal *w4, doublereal *eigmax, doublereal *tol, integer *itmax, integer *iters, integer *nlev, integer *ilev, integer *nlev_real__, integer *mgsolv, integer *iok, integer *iinfo, doublereal *epsiln, doublereal *errtol, doublereal *omega, integer *nu1, integer *nu2, integer *mgsmoo, integer *ipc, doublereal *rpc, doublereal *pc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *tru);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int rich_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int rich7_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *fc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes \n
 *:ref: mresid7_1s__ 14 13 4 4 4 4 7 7 7 7 7 7 7 7 7 \n
 *:ref: mresid27_1s__ 14 23 4 4 4 4 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 
    @author Michael Holst */
VEXTERNC int rich27_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *fc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *one, doublereal *onw, doublereal *ue, doublereal *uw, doublereal *un, doublereal *us, doublereal *une, doublereal *unw, doublereal *use, doublereal *usw, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes \n
 *:ref: wjac_ 14 18 4 4 4 4 7 7 7 7 7 7 7 7 4 4 7 7 4 4 \n
 *:ref: gsrb_ 14 18 4 4 4 4 7 7 7 7 7 7 7 7 4 4 7 7 4 4 \n
 *:ref: sor_ 14 18 4 4 4 4 7 7 7 7 7 7 7 7 4 4 7 7 4 4 \n
 *:ref: rich_ 14 18 4 4 4 4 7 7 7 7 7 7 7 7 4 4 7 7 4 4 \n
 *:ref: cghs_ 14 18 4 4 4 4 7 7 7 7 7 7 7 7 4 4 7 7 4 4 
    @author Michael Holst */
VEXTERNC int smooth_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint, integer *meth);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int sor_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int sor7_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *fc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes \n
 *:ref: mresid7_1s__ 14 13 4 4 4 4 7 7 7 7 7 7 7 7 7 \n
 *:ref: mresid27_1s__ 14 23 4 4 4 4 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 
    @author Michael Holst */
VEXTERNC int sor27_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *fc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *one, doublereal *onw, doublereal *ue, doublereal *uw, doublereal *un, doublereal *us, doublereal *une, doublereal *unw, doublereal *use, doublereal *usw, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int wjac_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *ac, doublereal *cc, doublereal *fc, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes 
    @author Michael Holst */
VEXTERNC int wjac7_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *fc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);
/** @brief  Library VPMG prototypes \n
 *:ref: mresid7_1s__ 14 13 4 4 4 4 7 7 7 7 7 7 7 7 7 \n
 *:ref: mresid27_1s__ 14 23 4 4 4 4 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 
    @author Michael Holst */
VEXTERNC int wjac27_(integer *nx, integer *ny, integer *nz, integer *ipc, doublereal *rpc, doublereal *oc, doublereal *cc, doublereal *fc, doublereal *oe, doublereal *on, doublereal *uc, doublereal *one, doublereal *onw, doublereal *ue, doublereal *uw, doublereal *un, doublereal *us, doublereal *une, doublereal *unw, doublereal *use, doublereal *usw, doublereal *x, doublereal *w1, doublereal *w2, doublereal *r__, integer *itmax, integer *iters, doublereal *errtol, doublereal *omega, integer *iresid, integer *iadjoint);


#endif /* _VPMG_H_ */

