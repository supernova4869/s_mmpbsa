/**
 * @defgroup Zslu Zslu class
 * @brief    Wrapper class for a generic sparse direct solver.
 */

/**
 *  @file       zslu.h
 *  @ingroup    Zslu
 *  @brief      Class Zslu: Wrapper class for a generic sparse direct solver.
 *  @author     Michael Holst
 *  @note       None
 *  @version    $Id: zslu.h,v 1.3 2010/08/12 05:18:37 fetk Exp $
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

#ifndef _ZSLU_H_
#define _ZSLU_H_

#include <mc/mc_base.h>

/**
 * @ingroup Zslu
 * @brief   Contains public data members for Zslu class
 * @author  Michael Holst
 */

struct sZslu {

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
  /** @brief number of rhs vectors                                */
    int nrhs;
  /** @brief work area size                                       */
    int lwork;
  /** @brief equilibration/scaling flag                           */
    char equed[1];

  /** @brief row-start (skey=0) or col-start (skey=1) in < j a >& < a > */
    int    *ia;
  /** @brief col indices (skey=0) or row-indices (skey=1) for < a > */
    int    *ja;
  /** @brief the nonzeros in the matrix                           */
    double *a;

  /** @brief solver-dependent structure for the matrix            */
    void *A;
  /** @brief solver-dependent structure for factored lower-triang */
    void *L;
  /** @brief solver-dependent structure for factored upper-triang */
    void *U;
  /** @brief solver-dependent structure for rhs(s)                */
    void *B;
  /** @brief solver-dependent structure for soln(s)               */
    void *X;
  /** @brief solver-dependent structure for status of affairs     */
    void *stat;
  /** @brief solver-dependent structure for options               */
    void *opts;
  /** @brief solver-dependent structure for factor work area      */
    void *work;

  /** @brief solver-dependent stuff                               */
    double *R;
  /** @brief solver-dependent stuff                               */
    double *C;
  /** @brief solver-dependent stuff                               */
    double *rhsb;
  /** @brief solver-dependent stuff                               */
    double *rhsx;
  /** @brief solver-dependent stuff                               */
    double *ferr;
  /** @brief solver-dependent stuff                               */
    double *berr;

  /** @brief row permutation vector                               */
    int *pr;
  /** @brief col permutation vector                               */
    int *pc;
  /** @brief factorization work area                              */
    int *et;

};

/**
 * @ingroup Zslu
 * @brief   Declaration of the Zslu class as the Zslu structure
 * @author  Michael Holst
 * @return  None
 */
typedef struct sZslu Zslu;

/*
 * ***************************************************************************
 * Class Zslu: Inlineable methods (zslu.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_ZSLU)
#else /* if defined(VINLINE_ZSLU) */
#endif /* if !defined(VINLINE_ZSLU) */

/**
 * @ingroup Zslu
 * @brief   The Zslu constructor
 * @author  Michael Holst
 * @note    Class Zslu: Non-inlineable methods (zslu.c)
 * @return  Pointer to a newly allocated (empty) Zslu class
 * @param   vmem  Memory management object
 * @param   skey  index for storage format (0=row-wise YSMP, 1=col-wise YSMP)
 * @param   m     number of rows
 * @param   n     number of cols
 * @param   nnz   number of nonzeros
 * @param   ia    row-start (skey=0) or col-start (skey=1) in < j a >& < a >
 * @param   ja    col indices (skey=0) or row-indices (skey=1) for < a >
 * @param   a     the nonzeros in the matrix
 */
VEXTERNC Zslu* Zslu_ctor(Vmem *vmem, int skey, int m, int n, int nnz,
    int *ia, int *ja, double *a);

/**
 * @ingroup Zslu
 * @brief   The Zslu destructor
 * @author  Michael Holst
 * @note    Class Zslu: Non-inlineable methods (zslu.c)
 * @return  None
 * @param   thee  Pointer to Zslu class
 */
VEXTERNC void Zslu_dtor(Zslu **thee);

/**
 * @ingroup Zslu
 * @brief   Sparse LU factor the system.
 * @author  Michael Holst
 * @note    Class Zslu: Non-inlineable methods (zslu.c)
 * @return  Success enumeration
 * @param   thee  Pointer to Zslu class
 */
VEXTERNC int Zslu_factor(Zslu *thee);

/**
 * @ingroup Zslu
 * @brief   Use sparse LU factors to back/forward solve a linear system.
 * @author  Michael Holst
 * @note    Class Zslu: Non-inlineable methods (zslu.c)
 * @return  Success enumeration
 * @param   thee  Pointer to Zslu class
 * @param   key   index for different solution methods
 * @param   b     Pointer to the array b
 * @param   x     Pointer to the array x
 */
VEXTERNC int Zslu_solve(Zslu *thee, int key, double *b, double *x);

/**
 * @ingroup Zslu
 * @brief   Print the exact current malloc usage.
 * @author  Michael Holst
 * @note    Class Zslu: Non-inlineable methods (zslu.c)
 * @return  None
 * @param   thee  Pointer to Zslu class
 */
VEXTERNC void Zslu_memChk(Zslu *thee);

/**
 * @ingroup Zslu
 * @brief   Calculate the log of the determinant of a SuperLU matrix.
 * @author  Nathan Baker (some restructuring by Michael Holst)
 * @note    Class Zslu: Non-inlineable methods (zslu.c) \n
 *          Uses SLU factorization and will be very slow for large matrices.\n
 * According to Sherry Li, author of SuperLU:\n
 *   The diagonal blocks of both L and U are stored in the L matrix,
 *   which is returned from dgstrf().  The L matrix is a supernodal matrix,
 *   its structure is called SCformat in supermatrix.h.  This is also
 *   illustrated by a small 5x5 example in Section 2.3 of the Users' Guide,
 *   see Figures 2.1 and 2.3.   This example is in the code
 *   EXAMPLE/superlu.c.  Since L is unit-diagonal, so the ones are not
 *   stored. Instead, the diagonal stored in L is really the diagonal for U.
 *   Therefore, you only need to extract those diagonal elements.  One
 *   routine that you can hack to get the diagonal is
 *   dPrint_SuperNode_Matrix() in dutil.c.  Another tricky part is the sign
 *   of the determinant. Since I am doing the following factorization Pr*A*Pc
 *   = LU, i.e., both row and column permutations may be applied, they are
 *   called perm_r and perm_c in the code. Their determinants will be 1 or
 *   -1, but you need to find out the sign by going through these
 *   permutations.
 * @return  the log of the determinant of a SuperLU matrix
 * @param   thee  Pointer to Zslu class
 */
VEXTERNC double Zslu_lnDet(Zslu *thee);

#endif /* _ZSLU_H_ */

