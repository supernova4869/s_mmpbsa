/**
 * @defgroup Vec Vec class
 * @brief    an Nx1 vector class.
 */

/**
 *  @file       vec.h
 *  @ingroup    Vec
 *  @brief      Class Vec: an Nx1 vector class.
 *  @version    $Id: vec.h,v 1.22 2010/08/12 05:18:36 fetk Exp $ 
 *  @author     Michael Holst
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


#ifndef _VEC_H_
#define _VEC_H_

#include <mc/mc_base.h>

#include <mc/mat.h>

/**
 * @ingroup Vec
 * @author  Michael Holst
 * @brief   Contains public data memebers for Vec class
 */
struct sVec {

    Vmem   *vmem;          /**< @brief the memory manager                   */
    int    iMadeVmem;      /**< @brief did i make vmem or was it inherited  */

    char   name[10];       /**< @brief character string name for this vector*/
    int    n;              /**< @brief number of components in the vector   */
    double *u;             /**< @brief the vector components                */
    int    iMallocU;       /**< @brief did i malloc u or was it passed in   */

};

/**
 * @brief   Declaration of the Vec class as the Vec structure
 * @ingroup Vec
 * @author  Michael Holst
 * @return  None
 */
typedef struct sVec Vec;

/*
 * ***************************************************************************
 * Class Vec: Inlineable methods (vec.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_BAM)
#else /* if defined(VINLINE_BAM) */
#endif /* if !defined(VINLINE_BAM) */

/**
 * @ingroup Vec 
 * @brief   The vector constructor (data array malloc'd by us).
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  Pointer to a newly allocated (empty) vector
 * @param   vmem    Memory management object
 * @param   name    character string name for this vector
 * @param   length  the length of the vector
 */
VEXTERNC Vec* Vec_ctor(Vmem *vmem, const char *name, int length);

/**
 * @ingroup Vec 
 * @brief   The vector constructor (data array passed in). 
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  Pointer to a newly allocated (empty) vector
 * @param   vmem    Memory management object
 * @param   name    character string name for this vector
 * @param   length  the length of the vector
 * @param   data    pointer to the initial value of the vector component
 */
VEXTERNC Vec* Vec_ctor2(Vmem *vmem, const char *name, int length, double *data);

/**
 * @ingroup Vec 
 * @brief   The vector destructor.
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  None
 * @param   thee  Pointer to the vector
 */
VEXTERNC void Vec_dtor(Vec **thee);

/**
 * @ingroup Vec 
 * @brief   Length of a vector.    
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  Length of a vector
 * @param   thee  Pointer to the vector
 */
VEXTERNC int Vec_len(Vec *thee);

/**
 * @ingroup Vec 
 * @brief   Address of a vector data.  
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  Address of a vector data
 * @param   thee  Pointer to the vector
 */
VEXTERNC double *Vec_addr(Vec *thee);

/**
 * @ingroup Vec 
 * @brief   Return value of component of vector.   
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  value of component of vector.   
 * @param   thee  Pointer to the vector
 * @param   i     index of the component of vector
 */
VEXTERNC double Vec_val(Vec *thee, int i);

/**
 * @ingroup Vec 
 * @brief   Set value of component of vector.   
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  None
 * @param   thee  Pointer to the vector
 * @param   i     index of the component of vector
 * @param   val   value of the component of vector
 */
VEXTERNC void Vec_set(Vec *thee, int i, double val);

/**
 * @ingroup Vec 
 * @brief   1-norm of a vector.  
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  1-norm of a vector
 * @param   thee  Pointer to the vector
 */
VEXTERNC double Vec_nrm1(Vec *thee);

/**
 * @ingroup Vec 
 * @brief   2-norm of a vector.  
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  2-norm of a vector
 * @param   thee  Pointer to the vector
 */
VEXTERNC double Vec_nrm2(Vec *thee);

/**
 * @ingroup Vec 
 * @brief   oo-norm of a vector.  
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  oo-norm of a vector.
 * @param   thee  Pointer to the vector
 */
VEXTERNC double Vec_nrm8(Vec *thee);

/**
 * @ingroup Vec 
 * @brief   1-norm of a vector.  
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  1-norm of the difference of two vectors
 * @param   thee  Pointer to the vector
 * @param   s     the source vector
 */
VEXTERNC double Vec_dif1(Vec *thee, Vec *s);

/**
 * @ingroup Vec 
 * @brief   2-norm of a vector.  
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  2-norm of the difference of two vectors
 * @param   thee  Pointer to the vector
 * @param   s     the source vector
 */
VEXTERNC double Vec_dif2(Vec *thee, Vec *s);

/**
 * @ingroup Vec 
 * @brief   oo-norm of a vector.  
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  oo-norm of the difference of two vectors
 * @param   thee  Pointer to the vector
 * @param   s     the source vector
 */
VEXTERNC double Vec_dif8(Vec *thee, Vec *s);

/**
 * @ingroup Vec 
 * @brief   dot product of two vectors.
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  dot product of two vectors
 * @param   thee  Pointer to the vector
 * @param   s     the source vector
 */
VEXTERNC double Vec_dot(Vec *thee, Vec *s);

/**
 * @ingroup Vec 
 * @brief   initialize a vector to be a constant.   
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  None
 * @param   thee  Pointer to the vector
 * @param   val   value to be initialized for the vector
 */
VEXTERNC void Vec_init(Vec *thee, double val);

/**
 * @ingroup Vec 
 * @brief   vector scale.
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  None
 * @param   thee  Pointer to the vector
 * @param   val   value to be scaled to the vector
 */
VEXTERNC void Vec_scal(Vec *thee, double val);

/**
 * @ingroup Vec 
 * @brief   vector copy.
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  None
 * @param   thee  Pointer to the vector
 * @param   s     the source vector
 */
VEXTERNC void Vec_copy(Vec *thee, Vec *s);

/**
 * @ingroup Vec 
 * @brief   scalar times vector plus vector.  
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  None
 * @param   thee  Pointer to the vector
 * @param   s     the source vector
 * @param   val   coeficient for scaling the source vector
 */
VEXTERNC void Vec_axpy(Vec *thee, Vec *s, double val);

/**
 * @ingroup Vec 
 * @brief   print a vector.
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  None
 * @param   thee  Pointer to the vector
 */
VEXTERNC void Vec_print(Vec *thee);

/**
 * @ingroup Vec 
 * @brief   print a vector.
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  None
 * @param   thee   Pointer to the vector
 * @param   fname  the output file name
 * @param   pflag  0 ==> write, 1 ==> append
 */
VEXTERNC void Vec_printSp(Vec *thee, char *fname, int pflag);

/**
 * @ingroup Vec 
 * @brief   Apply inverse of diagonal to a vector.           
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  None
 * @param   thee  Pointer to the vector
 * @param   amat  system matrix
 * @param   f     source vector slot
 */
VEXTERNC void Vec_diagScale(Vec *thee, Mat *amat, Vec *f);

/**
 * @ingroup Vec 
 * @brief   Matrix times vector.
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  None
 * @param   thee  Pointer to the vector
 * @param   amat  system matrix
 * @param   v     the source vector
 * @param   key   which of A or A' to use for the matvec \n
 *                (0=A, 1=A', 2=A(accumulate), 3=A'(accumulate))\n
 *                (NOTE: 0=prolongate, 1=restrict, for PRO)
 * @param   part  which part of matrix to use for the matvec \n
 *                (0=A=D+L+U, 1=D, 2=D+L, 3=D+U, 4=L+U) \n
 *                (NOTE: part!=0 only well-defined for DRC) 
 */
VEXTERNC void Vec_matvec(Vec *thee, Mat *amat, Vec *v, int key, int part);

/**
 * @ingroup Vec 
 * @brief   Generic smoothing operator.      
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  None
 * @param   thee    Pointer to the vector
 * @param   amat    system matrix
 * @param   f       source vector slot
 * @param   w       the work block vector
 * @param   key     smooth with A or A' (0=A, 1=A') 
 * @param   ioflag  debug output level (0=normal, 1=none, 2=lots, ... )
 * @param   meth    smoother choice (0=jac, 1=gs, ... ) 
 * @param   adj     adjoint choice (0=normal, 1=adjoint)
 * @param   itmax   number of iterations to do (the maximum allowed) 
 * @param   etol    error tolerance (currently ignored) 
 * @param   omega   SOR parameter (currently ignored)
 */
VEXTERNC void Vec_smooth(Vec *thee, Mat *amat, Vec *f, Vec *w,
    int key, int ioflag, int meth, int adj, int itmax, double etol,
    double omega);

/**
 * @ingroup Vec 
 * @brief   Jacobi iteration. 
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  None
 * @param   thee    Pointer to the vector
 * @param   amat    system matrix
 * @param   f       source vector slot
 * @param   w       the work block vector
 * @param   key     smooth with A or A' (0=A, 1=A') 
 * @param   ioflag  debug output level (0=normal, 1=none, 2=lots, ... )
 * @param   itmax   number of iterations to do (the maximum allowed) 
 * @param   etol    error tolerance (currently ignored) 
 * @param   omega   iteration parameter
 */
VEXTERNC void Vec_jac(Vec *thee, Mat *amat, Vec *f, Vec *w,
    int key, int ioflag, int itmax, double etol, double omega);

/**
 * @ingroup Vec 
 * @brief   Gauss-Seidel iteration.   
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  None
 * @param   thee    Pointer to the vector
 * @param   amat    system matrix
 * @param   f       source vector slot
 * @param   w       the work block vector
 * @param   key     smooth with A or A' (0=A, 1=A') 
 * @param   ioflag  debug output level (0=normal, 1=none, 2=lots, ... )
 * @param   adj     adjoint choice (0=normal, 1=adjoint) 
 * @param   itmax   number of iterations to do (the maximum allowed) 
 * @param   etol    error tolerance (currently ignored) 
 */
VEXTERNC void Vec_gs(Vec *thee, Mat *amat, Vec *f, Vec *w,
    int key, int ioflag, int adj, int itmax, double etol);

/**
 * @ingroup Vec 
 * @brief   Apply boundary condition.
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  None
 * @param   thee    Pointer to the vector
 * @param   mat     system matrix
 * @param   key     index for setting zeroes in the row/column components
 */
VEXTERNC void Vec_bnd(Vec *thee, Mat *mat, int key);

/**
 * @ingroup Vec 
 * @brief   Log of abs value of entries of block vector, shifted by < val >
 * @author  Michael Holst
 * @note    Class Vec: Non-Inlineable methods (vec.c) 
 * @return  None
 * @param   thee    Pointer to the vector
 * @param   val     value to be shifted to abs value of entries of block vector
 */
VEXTERNC void Vec_absLog(Vec *thee, double val);

#endif /* _VEC_H_ */

