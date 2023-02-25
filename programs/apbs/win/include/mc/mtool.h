/**
 *  @file       mtool.h
 *  @ingroup    global_mc
 *  @brief      Some sparse matrix tools.
 *  @author     Michael Holst
 *  @note       None
 *  @version    $Id: mtool.h,v 1.30 2010/08/12 05:18:35 fetk Exp $ 
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


#ifndef _MTOOL_H_
#define _MTOOL_H_

#include <mc/mc_base.h>

/*
 * ***************************************************************************
 * Class MTOOL: Parameters and datatypes                  
 * ***************************************************************************
 */

/**
 * @ingroup global_mc
 * @brief   A sparse matrix tool 
 * @author  Michael Holst
 * @note    Class MTOOL: Parameters and datatypes 
 */

typedef struct Link {
    /** @brief the index of array */
    int         idx;
    /** @brief Point to the next element in the list */
    struct Link *next;
} Link;


/**
 * @ingroup global_mc
 * @brief   A sparse matrix tool
 * @author  Michael Holst
 * @note    Class MTOOL: Parameters and datatypes 
 */

typedef struct LinkA {
    /** @brief the index of array */
    int         idx;
    /** @brief value of the sparse matrix component */
    double      val;
    /** @brief Point to the next element in the list */
    struct LinkA *next;
} LinkA;

/**
 * @ingroup global_mc
 * @brief   A sparse matrix tool 
 * @author  Michael Holst
 * @note    Class MTOOL: Parameters and datatypes 
 */

typedef struct LinkRC {
    /** @brief the index of array */
    int         idx;
    /** @brief Point to the next element in the list */
    struct LinkRC *next;

    /** @brief the index of trans array */
    int         idxT;
    /** @brief Point to the next element in the trans list */
    struct LinkRC *nextT;

} LinkRC;

/**
 * @ingroup global_mc
 * @brief   A sparse matrix tool 
 * @author  Michael Holst
 * @note    Class MTOOL: Parameters and datatypes 
 */

typedef struct LinkRCS {
    /** @brief the index of array */
    int         idx;
    /** @brief Point to the next element in the list */
    struct LinkRC *next;

    /** @brief the index of trans array */
    int         idxT;
    /** @brief Point to the next element in the trans list */
    struct LinkRC *nxtT;

    /** @brief value of the sparse matrix component */
    double      val;

} LinkRCS;

/**
 * @ingroup global_mc
 * @brief   A sparse matrix tool
 * @author  Michael Holst
 * @note    Class MTOOL: Parameters and datatypes 
 */

typedef struct LinkRCN {
    /** @brief the index of array */
    int         idx;
    /** @brief Point to the next element in the list */
    struct LinkRC *next;

    /** @brief the index of trans array */
    int         idxT;
    /** @brief Point to the next element in the trans list */
    struct LinkRC *nxtT;

    /** @brief value of the sparse matrix component */
    double      val;
    /** @brief value of the trans sparse matrix component */
    double      valT;

} LinkRCN;

/**
 * @ingroup global_mc
 * @brief   A sparse matrix tool
 * @author  Michael Holst
 * @note    Class MTOOL: Parameters and datatypes 
 */

typedef enum MATformat {
    ZERO_FORMAT,
    DRC_FORMAT,
    ROW_FORMAT,
    COL_FORMAT,
    SLU_FORMAT,
    RLN_FORMAT,
    CLN_FORMAT,
    XLN_FORMAT,
    RFL_FORMAT,
    CFL_FORMAT
} MATformat;

/**
 * @ingroup global_mc
 * @brief   matrix symmantec options
 * @author  Michael Holst
 * @note    Class MTOOL: Parameters and datatypes 
 */

typedef enum MATsym {
    ISNOT_SYM,
    IS_SYM,
    STRUC_SYM
} MATsym;

/**
 * @ingroup global_mc
 * @brief   matrix mirrow options 
 * @author  Michael Holst
 * @note    Class MTOOL: Parameters and datatypes 
 */

typedef enum MATmirror {
    ISNOT_MIRROR,
    IS_MIRROR
} MATmirror;

/**
 * @ingroup global_mc
 * @brief   the sparse matrix impl
 * @author  Michael Holst
 * @note    Class MTOOL: Parameters and datatypes 
 */

typedef enum MATimpl {
    ISNOT_IMPL,
    IS_IMPL
} MATimpl;

/**
 * @ingroup global_mc
 * @brief   different states of the sparse matrix 
 * @author  Michael Holst
 * @note    Class MTOOL: Parameters and datatypes 
 */

typedef enum MATstate {
    NULL_STATE,
    ZERO_STATE,
    ASSEMBLED_STATE,
    FACTORED_STATE
} MATstate;

/*
 * ***************************************************************************
 * Class MTOOL: Inlineable methods (mtool.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_BAM)
#else /* if defined(VINLINE_BAM) */
#endif /* if !defined(VINLINE_BAM) */


/**
 * @ingroup global_mc
 * @brief   Set or add a value to a matrix row/column. 
 * @author  Michael Holst
 * @note    Class MTOOL: Non-inlineable methods (mtool.c)
 * @return  None
 * @param   IA   pos in JA/offU/offL for row/col start
 * @param   JA   row/col indices for nonzeros in col/row
 * @param   A    packed nozeros:
 *               DRC: [ diag ; offU ; offL ]
 *               ROW: [ offU ]
 *               COL: [ offL ] 
 *               RFL: [ everything; stored row-wise ]
 *               CFL: [ everything; stored col-wise ]       
 * @param   key  0 ==> Set the value, 1 ==> Add the value
 * @param   i    an index
 * @param   j    an index 
 * @param   val  the value of the matrix component
 */
VEXTERNC void mPlaceit(int *IA, int *JA, double *A,
    int key, int i, int j, double val);

/**
 * @ingroup global_mc
 * @brief   Set or add a value to a linked matrix entry list.
 * @author  Michael Holst
 * @note    Class MTOOL: Non-inlineable methods (mtool.c)
 * @return  None
 * @param   mtpool  Pointer to a Class Vset
 * @param   key     0 ==> Set the value, 1 ==> Add the value
 * @param   count   Pointer to the inserted position
 * @param   i    an index
 * @param   j    an index 
 * @param   val  the value of the matrix component
 */
VEXTERNC void mContrib(Vset *mtpool,
    int key, int *count, int i, int j, double val);

/**
 * @ingroup global_mc
 * @brief   Add a link to a linked graph entry list.     
 * @authors Michael Holst and Stephen Bond
 * @note    Class MTOOL: Non-inlineable methods (mtool.c)
 * @return  None
 * @param   mtpool  Pointer to a Class Vset
 * @param   count   Pointer to the inserted position
 * @param   i       an index
 * @param   j       an index
 */
VEXTERNC void lContrib(Vset *mtpool, int *count, int i, int j);

/**
 * @ingroup global_mc
 * @brief   Builds an index (and value) array in transposed format given   
 *          the index and value array of a ROW or COL matrix. 
 * @author  Stephen Bond
 * @note    Class MTOOL: Non-inlineable methods (mtool.c)\n
 *          (flag == 0)  Build index array only, don't use matrix values.\n
 *          (flag == 1)  Build index array only, guard against zero values.\n
 *          (flag == 2)  Build index and val array, don't guard for zeros.\n
 *          (flag == 3)  Build index and val arrays, guard for zero values.\n\n
 *          If the original A is ROW (COL) it must be of size nxm (mxn).\n
 *          Hence, IJA, A, and work are n+1+numO, numO, and m respectively.    
 * @return  None
 * @param   vmem    Memory management object
 * @param   IJAT    Pointer to the index array
 * @param   AT      Pointer to the trasposed A array
 * @param   ATnumO  number of nonzeros in each row (col)
 * @param   IJA     integer structure [ IA ; JA ]
 * @param   A    packed nozeros:
 *               DRC: [ diag ; offU ; offL ]
 *               ROW: [ offU ]
 *               COL: [ offL ] 
 *               RFL: [ everything; stored row-wise ]
 *               CFL: [ everything; stored col-wise ]       
 * @param   n       an index
 * @param   m       an index
 * @param   flag    0 ==> Build index array only, don't use matrix values.\n
 *                  1 ==> Build index array only, guard against zero values.\n
 *                  2 ==> Build index and val array, don't guard for zeros.\n
 *                  3 ==>  Build index and val arrays, guard for zero values.
 * @param   work    Pointer to the work array with zeros
 */
VEXTERNC void mBuildGraphT(Vmem *vmem, int **IJAT, double **AT, int *ATnumO,
    int *IJA, double *A, int n, int m, int flag, int *work);

/**
 * @ingroup global_mc
 * @brief   Produce the sparse triple matrix product:   B = R*A*P,
 *          where A is a sparse NxN square matrix, R is a sparse 
 *          MxN retangular matrix, and P is a sparse NxM rectangular 
 *          matrix, with M < N (possibly M << N).\n
 *          The result is a smaller (but possibly much more dense) 
 *          square MxM matrix B.  
 * @author  Michael Holst
 * @note    Class MTOOL: Non-inlineable methods (mtool.c)\n
 * @verbatim
 *           The input matrices R,A,P are assumed to have the following 
 *           storage formats:  R is stored column-wise (ROW-format), P is
 *           stored row-wise (ROW-format), and A is stored in one of three
 *           forms, namely, either row-wise (ROW), col-wise (COL), or by
 *           diagonal followed by upper-triangle row-wise and then by lower
 *           triangle columne-wise (DRC).  The resulting B is stored in the
 *           same format as the input matrix A.  I.e., the matrices have the
 *           following possible storage format combinations:
 *
 *           The matrices A and B in DRC format:
 *
 *                B  =       R       *    A    *  P
 *
 *               \--   | | | | | | |   \------   ---
 *               |\- = | | | | | | | * |\----- * ---
 *               ||\   | | | | | | |   ||\----   ---
 *                                     |||\---   ---
 *                                     ||||\--   ---
 *                                     |||||\-   ---
 *                                     ||||||\   ---
 *
 *           Or, the matrices A and B in ROW format:
 *
 *                B  =       R       *    A    *  P
 *
 *               ---   | | | | | | |   -------   ---
 *               --- = | | | | | | | * ------- * ---
 *               ---   | | | | | | |   -------   ---
 *                                     -------   ---
 *                                     -------   ---
 *                                     -------   ---
 *                                     -------   ---
 *
 *           Or, finally, the matrices A and B in COL format:
 *
 *                B  =       R       *    A    *  P
 *
 *               |||   | | | | | | |   |||||||   ---
 *               ||| = | | | | | | | * ||||||| * ---
 *               |||   | | | | | | |   |||||||   ---
 *                                     |||||||   ---
 *                                     |||||||   ---
 *                                     |||||||   ---
 *                                     |||||||   ---
 *
 *           We compute the product B = R*A*P by splitting up A into its
 *           three structures A=D+L+U, and then by splitting up the sum:
 *
 *               B = R*D*P + R*L*P + R*U*P
 *
 *           In the case of A in ROW form, only the last term is
 *           present.  In the case of COL form, only the middle
 *           term is present.  In the DRC case, all terms are present.
 *           Note that in the ROW and COL forms, U and L are not
 *           actually upper and lower triangular, but general rectangular
 *           matrices; this is all that we need here.  In the DRC form,
 *           they are strict triangles, but we never use this fact directly.
 *
 *           The first and last terms are easy to compute using the 
 *           datastructures used for R,D,U,P.  However, the middle term
 *           R*L*P is not straightforward to compute directly, since the
 *           lower-triangle L of A is stored column-wise rather than
 *           row-wise.  However, this product can be easily formed by 
 *           computing its transpose instead:
 *
 *               B = R*D*P + (P^T*L^T*R^T)^T + R*U*P
 *
 *           The product P^T*L^T*R^T is easy to compute using the column-wise
 *           structure of L, hence row-wise structure of L^T; the result is 
 *           then transposed and added into B.
 *
 *           The product is formed component-wise; no temporary matrices
 *           are stored, except for a linked-list structure to build the
 *           non-zero structure of B, which is not possible to predict
 *           without actually doing the product (at least symbolically).
 * @endverbatim
 * @return  None
 * @param   vmem    Memory management object
 * @param   frmt    possible format types of this matrix
 * @param   sym     symmetry keys for the matrix
 * @param   m       an index
 * @param   n       an index
 * @param   numO    num of nonzeros we are actually storing in the strict upper-triangle of
 *                  matrix. (DRC only)
 * @param   numA    num of nonzeros we are actually storing, counting the diagonal, the strict
 *                  upper-triangle, and also the strict lower-triangle if we are actually
 *                  storing the lower-triangle (sym=0).
 * @param   ijb     Pointer to memory location of the arry ijb
 * @param   b       Pointer to memory location of the arry b
 * @param   ib      Pointer to memory location of the arry ib
 * @param   jb      Pointer to memory location of the arry jb
 * @param   diagb   Pointer to memory location of the arry diagb
 * @param   offUb   Pointer to memory location of the arry offUb
 * @param   offLb   Pointer to memory location of the arry offLb
 * @param   ir      Pointer to the arry ir
 * @param   jr      Pointer to the arry jr
 * @param   r       Pointer to the arry r
 * @param   ia      Pointer to the arry ia
 * @param   ja      Pointer to the arry ja
 * @param   a       Pointer to the arry a
 * @param   diag    Pointer to the arry diag
 * @param   offU    Pointer to the arry offU
 * @param   offL    Pointer to the arry offL
 * @param   ip      Pointer to the arry ip
 * @param   jp      Pointer to the arry jp
 * @param   p       Pointer to the arry p
 */
VEXTERNC void buildG(Vmem *vmem, MATformat frmt, MATsym sym,
    int m, int n, int *numO, int *numA,
    int **ijb, double **b,
    int **ib, int **jb, double **diagb, double **offUb, double **offLb,
    int *ir, int *jr, double *r,
    int *ia, int *ja, double *a,
    double *diag, double *offU, double *offL,
    int *ip, int *jp, double *p);

#endif /* _MTOOL_H_ */

