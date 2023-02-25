/**
 * @defgroup Slu Slu class
 * @brief    Wrapper class for a generic sparse direct solver.
 */

/**
 *  @file       slu.h
 *  @ingroup    Slu
 *  @brief      Class Slu: Wrapper class for a generic sparse direct solver.
 *  @author     Michael Holst and Stephen Bond
 *  @note       None
 *  @version    $Id: slu.h,v 1.9 2010/08/12 05:18:36 fetk Exp $
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

#ifndef _SLU_H_
#define _SLU_H_

#include <mc/mc_base.h>

/**
 * @ingroup Slu
 * @brief   Contains public data members for Slu class
 * @author  Michael Holst and Stephen Bond
 */

struct sSlu {

  /** @brief the memory manager */
    Vmem *vmem;
  /** @brief did i make vmem or was it inherited */
    int  iMadeVmem;

  /** @brief status of the factors */
    int statLU;

  /** @brief storage format (0=row-wise YSMP, 1=col-wise YSMP)    */
    int skey;
  /** @brief number of rows                                       */
    int m;
  /** @brief number of cols                                       */
    int n;
  /** @brief number of nonzeros                                   */
    int nnz;

  /** @brief row-start (skey=0) or col-start (skey=1) in < j a >& < a > */
    int    *ia;
  /** @brief col indices (skey=0) or row-indices (skey=1) for < a > */
    int    *ja;
  /** @brief the nonzeros in the matrix                           */
    double *a;

  /** @brief solver-dependent structure for factor work area      */
    void *work;
};

/**
 * @ingroup Slu
 * @brief   Declaration of the Slu class as the Slu structure
 * @author  Michael Holst and Stephen Bond
 * @return  None
 */
typedef struct sSlu Slu;

/*
 * ***************************************************************************
 * Class Slu: Inlineable methods (slu.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_SLU)
#else /* if defined(VINLINE_SLU) */
#endif /* if !defined(VINLINE_SLU) */

/**
 * @ingroup Slu
 * @brief   The Slu constructor
 * @author  Michael Holst and Stephen Bond
 * @note    Class Slu: Non-inlineable methods (slu.c)
 * @return  Pointer to a newly allocated (empty) Slu class
 * @param   vmem  Memory management object
 * @param   skey  index for storage format (0=row-wise YSMP, 1=col-wise YSMP)
 * @param   m     number of rows
 * @param   n     number of cols
 * @param   nnz   number of nonzeros
 * @param   ia    row-start (skey=0) or col-start (skey=1) in < j a >& < a >
 * @param   ja    col indices (skey=0) or row-indices (skey=1) for < a >
 * @param   a     the nonzeros in the matrix
 */
VEXTERNC Slu* Slu_ctor(Vmem *vmem, int skey, int m, int n, int nnz,
    int *ia, int *ja, double *a);

/**
 * @ingroup Slu
 * @brief   The Slu destructor
 * @author  Michael Holst and Stephen Bond
 * @note    Class Slu: Non-inlineable methods (slu.c)
 * @return  None
 * @param   thee  Pointer to the Slu class
 */
VEXTERNC void Slu_dtor(Slu **thee);

/**
 * @ingroup Slu
 * @brief   Sparse LU factor the system.
 * @author  Michael Holst and Stephen Bond
 * @note    Class Slu: Non-inlineable methods (slu.c)
 * @return  Sparse LU factor the system
 * @param   thee  Pointer to the Slu class
 */
VEXTERNC int Slu_factor(Slu *thee);

/**
 * @ingroup Slu
 * @brief   Use sparse LU factors to back/forward solve a linear system.
 * @author  Michael Holst and Stephen Bond
 * @note    Class Slu: Non-inlineable methods (slu.c)
 * @return  Success enumeration
 * @param   thee  Pointer to the Slu class
 * @param   key   index for different solution methods
 * @param   b     Pointer to the array b
 * @param   x     Pointer to the array x
 */
VEXTERNC int Slu_solve(Slu *thee, int key, double *b, double *x);

/**
 * @ingroup Slu
 * @brief   Calculate the log of the determinant of a factored matrix.
 * @author  Stephen Bond
 * @note    Class Slu: Non-inlineable methods (slu.c)
 * @return  the log of the determinant of a factored matrix
 * @param   thee  Pointer to the Slu class
 */
VEXTERNC double Slu_lnDet(Slu *thee);

#endif /* _SLU_H_ */

