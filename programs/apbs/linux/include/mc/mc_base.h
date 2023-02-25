/**
 * @defgroup global_mc global_mc class
 * @brief    Global group for MC
 */

/**
 *  @file       mc_base.h
 *  @ingroup    global_mc
 *  @brief      The base (or foundation) header for MC.
 *  @author     Michael Holst
 *  @note       None 
 *  @version    $Id: mc_base.h,v 1.30 2010/08/12 05:18:40 fetk Exp $ 
 *
 *  @attention
 *  @verbatim
 *
 * MC = < Manifold Code >
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

#ifndef _MC_BASE_H_
#define _MC_BASE_H_

#include <maloc/maloc.h>
#include <punc/punc.h>

/*
 * ***************************************************************************
 * Some parameters 
 * ***************************************************************************
 */

/** @brief Some parameters */
#define SPARSE_CUTOFF 1000

/* define SPARSE_CUTOFF 60000 */
/* define SPARSE_CUTOFF 10000 */
/* define SPARSE_CUTOFF 1000 */
/* define SPARSE_CUTOFF 40 */
/* define SPARSE_CUTOFF 5 */

/** @brief max num of matrix blocks allowed    */
#define MAXV                      4 
/** @brief max num of boundary types allowed   */
#define VMAX_BDTYPE             512  
/** @brief max num of links per row in lnkmats */
#define LN_MAX_ENTRIES_PER_ROW 1000  

/** @brief  3=max quadrature order (3=exact for cubics)        */
#define VMAXO   10  
/** @brief 20=max number of polys/coefs (20=up to 3D cubics)   */ 
#define VMAXP   20 
/** @brief 20=max number of quadrature points allowed          */ 
#define VMAXQ   20  
/** @brief 10=max degrees of freedom per v/e/f/s allowed       */
#define VMAXDF  10 

/** @brief permutations of (0,1,2) to label rows in Re::spmt   */
#define VPMT_012 0 
/** @brief and Re::spmthi arrays, which describe automorphisms */
#define VPMT_021 1 
/** @brief of the set of surface quadrature points under       */
#define VPMT_102 2  
/** @brief permutations of vertices of a triangular face       */
#define VPMT_120 3  
/** @brief permutations of surface quadrature points */
#define VPMT_201 4
/** @brief permutations of volumn quadrature points */
#define VPMT_210 5

/** @brief struct quadInfo */
typedef struct quadInfo {
  /** @brief master element integration weight           */
    double w;    
  /** @brief master element integration point            */       
    double x[3]; 
  /** @brief values of local basis funcs at integ pt     */
    double phi[VMAXP];
  /** @brief x-derivs of local basis funcs at integ pt   */      
    double phix[VMAXP];
  /** @brief y-derivs of local basis funcs at integ pt   */
    double phiy[VMAXP]; 
  /** @brief z-derivs of local basis funcs at integ pt   */
    double phiz[VMAXP];  
} quadInfo;

/** @brief struct simHelper */
typedef struct simHelper {
  /** @brief current color                                     */
    int color;
  /** @brief row sum for the matrix diagonal entry             */
    int diag; 
  /** @brief global face numbers for global faces we keep      */
    int faceId[4];
  /** @brief mass for the element                              */ 
    double mass;  
  /** @brief error estimate for the element                    */
    double error;
  /** @brief baricenter                                        */ 
    double bc[3];
} simHelper;

/** @brief struct Emat */
typedef struct Emat {
    /** @brief  global index for the node */
    int NI[MAXV][VMAXP]; 
    /** @brief boundary type for the node */
    int TP[MAXV][VMAXP];
    /** @brief the element of the force block vector*/
    double F[MAXV][VMAXP];
    /** @brief the element of the stiffness matrix */
    double A[MAXV][MAXV][VMAXP][VMAXP];
    /** @brief the element of the mass matrix */
    double M[MAXV][MAXV][VMAXP][VMAXP];
    /** @brief the energy */
    double J;
} Emat;

/*
 * ***************************************************************************
 * Inlining via macros for speed
 * ***************************************************************************
 */

#if 1
/** @brief Inlining via macros for speed */
#   define VINLINE_APRX
/** @brief Inlining via macros for speed */
#   define VINLINE_BAM
/** @brief Inlining via macros for speed */
#   define VINLINE_GEM
/** @brief Inlining via macros for speed */
#   define VINLINE_MCSH
/** @brief Inlining via macros for speed */
#   define VINLINE_NAM
/** @brief Inlining via macros for speed */
#   define VINLINE_PDE
/** @brief Inlining via macros for speed */
#   define VINLINE_ZBLAS
/** @brief Inlining via macros for speed */
#   define VINLINE_ZSLU
#endif

/*
 * ***************************************************************************
 * General element support
 * ***************************************************************************
 */

#if 1
/** @brief General element support */
#   define VG_ELEMENT
#endif

/** @brief Some macros  */
#define VBOUNDARY(x)  (x)
/** @brief Some macros  */
#define VINTERIOR(x)  (!(x))
/** @brief Some macros  */
#define VDIRICHLET(x) VODD(x)
/** @brief Some macros  */
#define VNEUMANN(x)   VEVENP(x)

#endif /* _MC_BASE_H_ */

