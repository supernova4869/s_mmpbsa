/** 
 * @defgroup Bvec Bvec class
 * @brief    A block vector object.
 */

/**
 *  @file       bvec.h
 *  @ingroup    Bvec
 *  @brief      Class Bvec: a block vector object.
 *  @author     Michael Holst
 *  @note       None
 *  @version    $Id: bvec.h,v 1.40 2010/08/12 05:18:34 fetk Exp $ 
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


#ifndef _BVEC_H_
#define _BVEC_H_

#include <mc/mc_base.h>

#include <mc/vec.h>
#include <mc/bmat.h>

/**
 * @ingroup Bvec
 * @brief   Contains public data members for Bvec class
 * @author  Michael Holst
 */

struct sBvec {

    Vmem   *vmem;          /**< @brief the memory manager                   */
    int    iMadeVmem;      /**< @brief did i make vmem or was it inherited  */

    char   name[10];       /**< @brief character string name for this vector*/
    int    n;              /**< @brief number of components in the vector   */
    double *u;             /**< @brief the vector components                */
    int    iMallocU;       /**< @brief did i malloc u or was it passed in   */

    struct sBvec *coarse;   /**< @brief next coarser object in the hierarchy */
    struct sBvec *fine;     /**< @brief next finer object in the hierarchy   */

    int    numB;           /**< @brief num vector blocks                    */
    int    numR[MAXV];     /**< @brief num of rows in each block vector     */
    Vec    *VC[MAXV];      /**< @brief vector block substructure            */

};

/**
 * @brief   Declaration of the Bvec class as the Bvec structure
 * @ingroup Bvec
 * @author  Michael Holst
 * @return  None
 */
typedef struct sBvec Bvec;

/*
 * ***************************************************************************
 * Class Bvec: Inlineable methods (bvec.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_BAM)
    /**
     * @ingroup Bvec
     * @brief   Return the number of blocks.      
     * @author  Michael Holst
     * @note    Class Bvec: Non-inlineable methods (bvec.c) 
     * @return  the number of blocks.
     * @param   thee  pointer to the block vector
     */    
    VEXTERNC int Bvec_numB(Bvec *thee);
#else /* if defined(VINLINE_BAM) */

    /**
     * @ingroup Bvec
     * @brief   Return the number of blocks.      
     * @author  Michael Holst
     * @note    Class Bvec: Non-inlineable methods (bvec.c) 
     * @return  the number of blocks.
     * @param   thee  pointer to the block vector
     */    
#   define Bvec_numB(thee)              ((thee)->numB)
#endif /* if !defined(VINLINE_BAM) */

/**
 * @ingroup Bvec
 * @brief   The block vector constructor (data array malloc'd by us).
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  Pointer to a newly allocated (empty) block vector
 * @param   vmem  Memory management object
 * @param   name  character string name for this vector
 * @param   pnumB num vector blocks
 * @param   pnumR num of rows in each block vector
 */
VEXTERNC Bvec* Bvec_ctor(Vmem *vmem,
    const char *name, int pnumB, int pnumR[MAXV]);

/**
 * @ingroup Bvec
 * @brief   The block vector constructor (data array passed in). 
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  Pointer to a newly allocated (empty) block vector
 * @param   vmem  Memory management object
 * @param   name  character string name for this vector
 * @param   pnumB num vector blocks
 * @param   pnumR num of rows in each block vector
 * @param   data  the vector components
 */
VEXTERNC Bvec* Bvec_ctor2(Vmem *vmem,
    const char *name, int pnumB, int pnumR[MAXV], double *data);

/**
 * @ingroup Bvec
 * @brief   The vector constructor (data array malloc'd by us). 
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  Pointer to a newly allocated (empty) block vector
 * @param   vmem    Memory management object
 * @param   name    character string name for this vector
 * @param   length  number of components in the vector
 */
VEXTERNC Bvec* Bvec_ctor3(Vmem *vmem,
    const char *name, int length);

/**
 * @ingroup Bvec
 * @brief   The vector constructor (data array malloc'd by us). 
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  Pointer to a newly allocated (empty) block vector
 * @param   vmem    Memory management object
 * @param   name    character string name for this vector
 * @param   length  number of components in the vector
 * @param   data    the vector components
 */
VEXTERNC Bvec* Bvec_ctor4(Vmem *vmem,
    const char *name, int length, double *data);

/**
 * @ingroup Bvec
 * @brief   The block vector destructor.   
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  None
 * @param   thee  pointer to the block vector
 */
VEXTERNC void Bvec_dtor(Bvec **thee);

/**
 * @ingroup Bvec
 * @brief   Create a set of vectors with duplicate attributes.  
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) \n
 *          We actually do this using contiguous memory, i.e., each
 *          vector follows contiguously from the preceeding one.  
 * @return  None
 * @param   thee  pointer to the block vector
 * @param   vecs  the created block vector
 * @param   num   the number of created vectors
 */
VEXTERNC void Bvec_createVectors(Bvec *thee, Bvec *vecs[], int num);

/**
 * @ingroup Bvec
 * @brief   Destroy set of vectors previously created by Bvec_createVectors.  
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  None
 * @param   thee  pointer to the block vector
 * @param   vecs  the created block vector
 * @param   num   the number of created vectors
 */
VEXTERNC void Bvec_destroyVectors(Bvec *thee, Bvec *vecs[], int num);

/**
 * @ingroup Bvec
 * @brief   Create a matrix from columns of vectors.  
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  None
 * @param   thee  pointer to the block vector
 * @param   vecs  the created block vector
 * @param   num   the number of created vectors
 * @param   mat   Pointer to the sparse matrix
 */
VEXTERNC void Bvec_createVecMat(Bvec *thee, Bvec *vecs[], int num, Mat **mat);

/**
 * @ingroup Bvec
 * @brief   Destroy matrix previously created by Bvec_createVecMat.       
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  None
 * @param   thee  pointer to the block vector
 * @param   vecs  the created block vector
 * @param   num   the number of created vectors
 * @param   mat   Pointer to the sparse matrix
 */
VEXTERNC void Bvec_destroyVecMat(Bvec *thee, Bvec *vecs[], int num, Mat **mat);

/**
 * @ingroup Bvec
 * @brief   Return the total number of rows.  
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  the total number of rows.
 * @param   thee  pointer to the block vector
 */
VEXTERNC int Bvec_len(Bvec *thee);

/**
 * @ingroup Bvec
 * @brief   Return a pointer to the total data array.       
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  a pointer to the total data arary
 * @param   thee  pointer to the block vector
 */
VEXTERNC double* Bvec_addr(Bvec *thee);

/**
 * @ingroup Bvec
 * @brief   Return value of component of vector.
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  value of component of vector
 * @param   thee  pointer to the block vector
 * @param   i     index of the block vector
 */
VEXTERNC double Bvec_val(Bvec *thee, int i);

/**
 * @ingroup Bvec
 * @brief   Set value of component of vector.   
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  None
 * @param   thee  pointer to the block vector
 * @param   i     index of the block vector
 * @param   val   value of the i-th component of the block vector
 */
VEXTERNC void Bvec_set(Bvec *thee, int i, double val);

/**
 * @ingroup Bvec
 * @brief   1-norm of a block vector.  
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  1-norm of a block vector.
 * @param   thee  pointer to the block vector
 */
VEXTERNC double Bvec_nrm1(Bvec *thee);

/**
 * @ingroup Bvec
 * @brief   2-norm of a block vector. 
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  2-norm of a block vector
 * @param   thee  pointer to the block vector
 */
VEXTERNC double Bvec_nrm2(Bvec *thee);

/**
 * @ingroup Bvec
 * @brief   oo-norm of a block vector.   
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  oo-norm of a block vector.
 * @param   thee  pointer to the block vector
 */
VEXTERNC double Bvec_nrm8(Bvec *thee);

/**
 * @ingroup Bvec
 * @brief   1-norm of the difference of two block vectors.              
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  1-norm of the difference of two block vectors.
 * @param   thee  pointer to the block vector
 * @param   v     the source block vector
 */
VEXTERNC double Bvec_dif1(Bvec *thee, Bvec *v);

/**
 * @ingroup Bvec
 * @brief   2-norm of the difference of two block vectors.    
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  2-norm of the difference of two block vectors.
 * @param   thee  pointer to the block vector
 * @param   v     the source block vector
 */
VEXTERNC double Bvec_dif2(Bvec *thee, Bvec *v);

/**
 * @ingroup Bvec
 * @brief   oo-norm of the difference of two block vectors.      
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  oo-norm of the difference of two block vectors.
 * @param   thee  pointer to the block vector
 * @param   v     the source block vector
 */
VEXTERNC double Bvec_dif8(Bvec *thee, Bvec *v);

/**
 * @ingroup Bvec
 * @brief   Dot product of two block vectors.      
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  Dot product of two block vectors.
 * @param   thee  pointer to the block vector
 * @param   v     the source block vector
 */
VEXTERNC double Bvec_dot(Bvec *thee, Bvec *v);

/**
 * @ingroup Bvec
 * @brief   Initialize a block vector.     
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  None
 * @param   thee  pointer to the block vector
 * @param   val   value to be initialized for the block vector
 */
VEXTERNC void Bvec_init(Bvec *thee, double val);

/**
 * @ingroup Bvec
 * @brief   Scale a block vector.  
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  None
 * @param   thee  pointer to the block vector
 * @param   val   value of the i-th component of the block vector
 */
VEXTERNC void Bvec_scal(Bvec *thee, double val);

/**
 * @ingroup Bvec
 * @brief   Copy a block vector.    
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  None
 * @param   thee  pointer to the block vector
 * @param   v     the source block vector
 */
VEXTERNC void Bvec_copy(Bvec *thee, Bvec *v);

/**
 * @ingroup Bvec
 * @brief   Saxpy for a block vector.   
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  None
 * @param   thee  pointer to the block vector
 * @param   v     the source block vector
 * @param   val   coeficient for scaling the source block vector
 */
VEXTERNC void Bvec_axpy(Bvec *thee, Bvec *v, double val);

/**
 * @ingroup Bvec
 * @brief   Print the block vector
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  None
 * @param   thee  pointer to the block vector
 */
VEXTERNC void Bvec_print(Bvec *thee);

/**
 * @ingroup Bvec
 * @brief   Print the block vector
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  None
 * @param   thee  pointer to the block vector
 * @param   fname the file name of the block vector
 */
VEXTERNC void Bvec_printSp(Bvec *thee, char *fname);

/**
 * @ingroup Bvec
 * @brief   Apply inverse of diagonal to a vector. 
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  None
 * @param   thee   pointer to the block vector
 * @param   A      system matrix
 * @param   f      source vector slot
 */
VEXTERNC void Bvec_diagScale(Bvec *thee, Bmat *A, Bvec *f);

/**
 * @ingroup Bvec
 * @brief   Block matrix times vector.
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  None
 * @param   thee  pointer to the block vector
 * @param   A     system matrix
 * @param   v     the source block vector
 * @param   key   index for different accumulation options
 * @param   part  index for future use
 */
VEXTERNC void Bvec_matvec(Bvec *thee, Bmat *A, Bvec *v, int key, int part);

/**
 * @ingroup Bvec
 * @brief   Generic smoothing operator (for a general YSMP-based matrix).
 * @authors Michael Holst and Stephen Bond
 * @note    Class Bvec: Non-inlineable methods (bvec.c) \n
 *          meth==1, adj==0  : Forward Gauss-Seidel\n
 *          u = (D + L)^{-1}*(f - U*u)\n
 *          meth==1, adj==1  : Adjoint Gauss-Seidel\n
 *          u = (D + U)^{-1}*(f - L*u)\n
 *          meth==2 : Self-adjoint or Symmetric Gauss-Seidel\n
 *          Apply (meth==1, adj==0) then (meth==1, adj==1)          
 * @return  None
 * @param   thee    pointer to the block vector
 * @param   A       system matrix
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
VEXTERNC void Bvec_smooth(Bvec *thee, Bmat *A, Bvec *f, Bvec *w,
    int key, int ioflag, int meth, int adj, int itmax, double etol,
    double omega);

/**
 * @ingroup Bvec
 * @brief   Apply boundary condition.
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  None
 * @param   thee  pointer to the block vector    
 * @param   bmat  the block sparse matrix
 * @param   key   index for setting zeroes in the row/column components.    
 */
VEXTERNC void Bvec_bnd(Bvec *thee, Bmat *bmat, int key);

/**
 * @ingroup Bvec
 * @brief   Print the exact current malloc usage.
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bvec.c) 
 * @return  None
 * @param   thee  pointer to the block vector
 */
VEXTERNC void Bvec_memChk(Bvec *thee);

/** 
 * @ingroup Bvec
 * @brief   Return the total number of rows.
 * @author  Michael Holst
 * @note    BLOCK-ONLY STUFF 
 * @return  the total number of rows
 * @param   thee  pointer to the block vector
 */
VEXTERNC int Bvec_numRT(Bvec *thee);

/** 
 * @ingroup Bvec
 * @brief   Return the size of a block.
 * @author  Michael Holst
 * @note    BLOCK-ONLY STUFF 
 * @return  the size of a block
 * @param   thee  pointer to the block vector
 * @param   i     index of the block vector
 */
VEXTERNC int Bvec_numRB(Bvec *thee, int i);

/** 
 * @ingroup Bvec
 * @brief   Return a pointer to the data in a block.   
 * @author  Michael Holst
 * @note    BLOCK-ONLY STUFF 
 * @return  a pointer to the data in a block.   
 * @param   thee  pointer to the block vector
 * @param   i     index of the block vector
 */
VEXTERNC double* Bvec_addrB(Bvec *thee, int i);

/** 
 * @ingroup Bvec
 * @brief   Return value of component of a vector.
 * @author  Michael Holst
 * @note    BLOCK-ONLY STUFF 
 * @return  value of component of a vector.
 * @param   thee  pointer to the block vector
 * @param   i     index of the block vector
 * @param   which index of the component in the vector
 */
VEXTERNC double Bvec_valB(Bvec *thee, int i, int which);

/** 
 * @ingroup Bvec
 * @brief   Set value of component of a vector.
 * @author  Michael Holst
 * @note    BLOCK-ONLY STUFF 
 * @return  None
 * @param   thee  pointer to the block vector
 * @param   i     index of the block vector
 * @param   which index of the component in the vector
 * @param   val   value of the i-th component of the block vector
 */
VEXTERNC void Bvec_setB(Bvec *thee, int i, int which, double val);

/** 
 * @ingroup Bvec
 * @brief   Add value to component of a vector. 
 * @author  Michael Holst
 * @note    BLOCK-ONLY STUFF 
 * @return  None
 * @param   thee  pointer to the block vector
 * @param   i     index of the block vector
 * @param   which index of the component in the vector
 * @param   val   value of the i-th component of the block vector
 */
VEXTERNC void Bvec_addToB(Bvec *thee, int i, int which, double val);

/** 
 * @ingroup Bvec
 * @brief   Initialize a particular block of a block vector. 
 * @author  Michael Holst
 * @note    BLOCK-ONLY STUFF 
 * @return  None
 * @param   thee  pointer to the block vector
 * @param   i     index of the block vector
 * @param   val   value of the i-th component of the block vector
 */
VEXTERNC void Bvec_initB(Bvec *thee, int i, double val);

/** 
 * @ingroup Bvec
 * @brief   Log of abs value of entries of block vector, shifted by < val >  
 * @author  Michael Holst
 * @note    BLOCK-ONLY STUFF 
 * @return  None
 * @param   thee  pointer to the block vector
 * @param   val   value of the i-th component of the block vector
 */
VEXTERNC void Bvec_absLog(Bvec *thee, double val);

/**
 * @ingroup Bvec
 * @brief   Executes one of a number of linear solvers. 
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bsolv.c) 
 * @return  None
 * @param   thee   pointer to the block vector
 * @param   A      system matrix
 * @param   f      source vector slot
 * @param   r      residual slot
 * @param   ut     block NODAL analytical solution
 * @param   key    index for different solution methods: 0=(Au=f), 1=(A'u=f)
 * @param   flag   index for i/o control
 * @param   itmax  number of iterations to do (the maximum allowed)
 * @param   etol   error tolerance
 * @param   prec   index for different preconditioners
 * @param   cycle  index for symmetric/nonsymmetric multigrids
 * @param   P      prolongation matrix maintained by aprx
 * @param   meth   method choice (0=slu,1=mg,2=cg,3=bcg,4=pcg,5=pbcg)
 */
VEXTERNC void Bvec_lmethod(Bvec *thee, Bmat *A, Bvec *f, Bvec *r, Bvec *ut,
    int key, int flag, int itmax, double etol, int prec, int cycle, Bmat *P,
    int meth);

/**
 * @ingroup Bvec
 * @brief   Sparse LU direct solver (requires only a nonsingular A).
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bsolv.c) 
 * @return  None
 * @param   thee   pointer to the block vector
 * @param   A      system matrix
 * @param   f      source vector slot
 * @param   r      residual slot
 * @param   ut     block NODAL analytical solution
 * @param   key    key ==0   --> Solve: A  u = f\n
 *                 key ==1   --> Solve: A' u = f  
 * @param   flag   Determines which "mode" we run in:\n
 *                 flag==0,2 --> Normal: normal i/o\n
 *                 flag==1,3 --> Silent: no i/o
 */
VEXTERNC void Bvec_slu(Bvec *thee, Bmat *A, Bvec *f, Bvec *r, Bvec *ut,
    int key, int flag);

/**
 * @ingroup Bvec
 * @brief   MG algorithm front end (NOT recursively callable).   
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bsolv.c)\n
 *          We use P to build an algebraic multilevel hierarchy in
 *          A,u(thee),f,r,ut from the fine-level-only representations,
 *          and then make the first call to the recursively-callable  
 *          core MG algorithm Bvec_mgCore(). 
 * @return  None
 * @param   thee   solution slot
 * @param   A      system matrix
 * @param   f      source vector slot
 * @param   r      residual slot
 * @param   ut     true solution slot
 * @param   key    just passed along to Bvec_mgCore()    
 * @param   flag   just passed along to Bvec_mgCore()    
 * @param   itmax  just passed along to Bvec_mgCore()    
 * @param   etol   just passed along to Bvec_mgCore()    
 * @param   prec   just passed along to Bvec_mgCore()    
 * @param   cycle  just passed along to Bvec_mgCore()    
 * @param   P      prolongation operators (linked list) 
 */
VEXTERNC void Bvec_mg(Bvec *thee, Bmat *A, Bvec *f, Bvec *r, Bvec *ut,
    int key, int flag, int itmax, double etol, int prec, int cycle, Bmat *P);

/**
 * @ingroup Bvec
 * @brief   Initialize the algebraic multilevel hierarchy for Bvec_mg
 * @authors Michael Holst and Stephan Bond
 * @note    Class Bvec: Non-inlineable methods (bsolv.c)\n
 *          We use P to build an algebraic multilevel hierarchy in
 *          A,u(thee),f,r,ut from the fine-level-only representations. 
 *          The vectors are created after the matrix to reduce the total 
 *          memory usage since Bmat_galerkin uses temporary arrays.  
 * @return  None
 * @param   thee   solution slot
 * @param   A      system matrix
 * @param   f      source vector slot
 * @param   r      residual slot
 * @param   P      prolongation operators (linked list)
 */
VEXTERNC void Bvec_mgInit(Bvec *thee, Bmat *A, Bvec *f, Bvec *r, Bmat *P);

/**
 * @ingroup Bvec
 * @brief   Destroy algebraic multilevel hierarchy created by Bvec_mgInit 
 * @authors Michael Holst and Stephan Bond
 * @note    Class Bvec: Non-inlineable methods (bsolv.c)\n
 *          Only the coarser level matrices and vectors are destructed. 
 * @return  None
 * @param   thee   solution slot
 * @param   A      system matrix
 * @param   f      source vector slot
 * @param   r      residual slot
 */
VEXTERNC void Bvec_mgDestroy(Bvec *thee, Bmat *A, Bvec *f, Bvec *r);

/**
 * @ingroup Bvec
 * @brief   CG method (ODIR algorithm, requires SPD operators A, B, and M).
 * @authors Michael Holst and Stephan Bond
 * @note    Class Bvec: Non-inlineable methods (bsolv.c)\n
 * @verbatim
 *           ODIR(M,B,A) (for Au=f, preconditioner B, inner-product M)
 *           ---------------------------------------------------------
 *
 *           Given u^0
 *           r^0 = f - A u^0
 *           p^0 = B r^0 / ||B r^0||
 *           for i=0,1,2,....
 *               alpha_i = <M e^i, p^i> / <M p^i, p^i>
 *               u^{i+1} = u^i + alpha_i p^i
 *               r^{i+1} = r^i - alpha_i A p^i
 *               gamma_i = <M B A p^i, p^i> / <M p^i, p^i>
 *               sigma_i = <M B A p^i, p^{i-1}> / <M p^{i-1}, p^{i-1}>
 *               v       = B A p^i - gamma_i p^i - sigma_i p^{i-1}
 *               p^{i+1} = v / ||v||
 *           end for
 *
 *           CG=ODIR(M=A,B=B,A=A)
 *           --------------------
 *
 *           Given u^0
 *           r^0 = f - A u^0
 *           p^0 = B r^0 / ||B r^0||
 *           for i=0,1,2,....
 *               alpha_i = <r^i, p^i> / <A p^i, p^i>
 *               u^{i+1} = u^i + alpha_i p^i
 *               r^{i+1} = r^i - alpha_i A p^i
 *               gamma_i = <B A p^i, A p^i> / <A p^i, p^i>
 *               sigma_i = <B A p^i, A p^{i-1}> / <A p^{i-1}, p^{i-1}>
 *               v       = B A p^i - gamma_i p^i - sigma_i p^{i-1}
 *               p^{i+1} = v / ||v||
 *           end for
 * @endverbatim
 * @return  None
 * @param   thee   pointer to the block vector
 * @param   A      system matrix
 * @param   f      source vector slot
 * @param   r      residual slot
 * @param   ut     block NODAL analytical solution
 * @param   key    index for different solution methods: 0=(Au=f), 1=(A'u=f)
 * @param   flag   0 --> Normal: check itmax, error tolerance; normal i/o\n
 *                 1 --> Silent: check only itmax; no i/o\n
 *                 2 --> Subcycle: check itmax, error tolerance; subcycle i/o\n
 *                 3 --> Subcycle: check itmax, error tolerance; no i/o
 * @param   itmax  number of iterations to do (the maximum allowed)
 * @param   etol   error tolerance
 * @param   prec   index for different preconditioners
 * @param   cycle  index for symmetric/nonsymmetric multigrids
 * @param   P      prolongation matrix maintained by aprx
 * @param   p      the block vector p
 * @param   ap     another block vector ap
 * @param   bap    third block vector bap
 * @param   po     the block vector po
 * @param   apo    the block vector apo
 * @param   tp     temporary block vector tp
 */
VEXTERNC void Bvec_cg(Bvec *thee, Bmat *A, Bvec *f, Bvec *r, Bvec *ut,
    int key, int flag, int itmax, double etol, int prec, int cycle, Bmat *P,
    Bvec *p, Bvec *ap, Bvec *bap, Bvec *po, Bvec *apo, Bvec *tp);

/**
 * @ingroup Bvec
 * @brief   BCG method (pure bi-orthogonal extension of CG). 
 * @author  Michael Holst
 * @note    Class Bvec: Non-inlineable methods (bsolv.c)\n
 *          This is a pure Petrov-Galerkin formulation of the 
 *          ODIR implementation of CG, and reduces to (exactly) 
 *          the ODIR implementation of CG in the case of a    
 *          symmetric (possibly indefinite) system matrix A and a  
 *          symmetric (possibly indefinite) preconditioner B.   
 *          (However, the operation count of BCG is rougly double
 *          that of CG.)    
 * @verbatim
 *           PG-ODIR(M,B,A) (for Au=f, preconditioner B, inner-product M)
 *           ------------------------------------------------------------
 *
 *           Given u^0
 *           r^0 = f - A  u^0
 *           s^0 = f - A' u^0
 *           p^0 = B r^0 / ||B r^0||
 *           q^0 = B's^0 / ||B's^0||
 *           for i=0,1,2,....
 *               alpha_i = <M e^i, q^i> / <M p^i, q^i>
 *               u^{i+1} = u^i + alpha_i p^i
 *               r^{i+1} = r^i - alpha_i A p^i
 *               gamma_i = <M B A p^i, q^i> / <M p^i, q^i>
 *               sigma_i = <M B A p^i, q^{i-1}> / <M p^{i-1}, q^{i-1}>
 *               delta_i = <M B A p^{i-1}, q^i> / <M p^{i-1}, q^{i-1}>
 *               v       = B A  p^i - gamma_i p^i - sigma_i p^{i-1}
 *               w       = B'A' q^i - gamma_i q^i - delta_i q^{i-1}
 *               p^{i+1} = v / ||v||
 *               q^{i+1} = w / ||w||
 *           end for
 *
 *           BCG=PG-ODIR(M=A,B=B,A=A)
 *           ------------------------
 *           Given u^0
 *           r^0 = f - A  u^0
 *           s^0 = f - A' u^0
 *           p^0 = B r^0 / ||B r^0||
 *           q^0 = B's^0 / ||B's^0||
 *           for i=0,1,2,....
 *               alpha_i = <r^i, q^i> / <A p^i, q^i>
 *               u^{i+1} = u^i + alpha_i p^i
 *               r^{i+1} = r^i - alpha_i A p^i
 *               gamma_i = <B A p^i, A' q^i> / <A p^i, q^i>
 *               sigma_i = <B A p^i, A' q^{i-1}> / <A p^{i-1}, q^{i-1}>
 *               delta_i = <A p^{i-1}, B'A' q^i> / <A p^{i-1}, q^{i-1}>
 *               v       = B A  p^i - gamma_i p^i - sigma_i p^{i-1}
 *               w       = B'A' q^i - gamma_i q^i - delta_i q^{i-1}
 *               p^{i+1} = v / ||v||
 *               q^{i+1} = w / ||w||
 *           end for
 * @endverbatim
 * @return  None
 * @param   thee   pointer to the block vector
 * @param   A      system matrix
 * @param   f      source vector slot
 * @param   r      residual slot
 * @param   ut     block NODAL analytical solution
 * @param   key    smooth with A or A' (0=A, 1=A')
 * @param   flag   determines which "mode" we run in:\n
 *                 0 --> Normal: check itmax, error tolerance; normal i/o\n
 *                 1 --> Silent: check only itmax; no i/o\n
 *                 2 --> Subcycle: check itmax, error tolerance; subcycle i/o\n
 *                 3 --> Subcycle: check itmax, error tolerance; no i/o
 * @param   itmax  number of iterations to do (the maximum allowed)
 * @param   etol   error tolerance
 * @param   prec   index for different preconditioners
 * @param   cycle  index for symmetric/nonsymmetric multigrids
 * @param   P      prolongation matrix maintained by aprx
 * @param   s      residual block vector s: f - A'u
 * @param   p      the block vector p
 * @param   ap     another block vector ap
 * @param   bap    third block vector bap
 * @param   po     the block vector po
 * @param   apo    the block vector apo
 * @param   q      the block vector q
 * @param   atq    the block vector atq
 * @param   btatq  the block vector btatq
 * @param   qo     the block vector qo
 * @param   atqo   the block vector atqo
 */
VEXTERNC void Bvec_bcg(Bvec *thee, Bmat *A, Bvec *f, Bvec *r, Bvec *ut,
    int key, int flag, int itmax, double etol, int prec, int cycle, Bmat *P,
    Bvec *s,
    Bvec *p, Bvec *ap,  Bvec *bap,  Bvec *po,  Bvec *apo,
    Bvec *q, Bvec *atq, Bvec *btatq, Bvec *qo, Bvec *atqo);

/**
 * @ingroup Bvec
 * @brief   Calculate the eigenvector corresponding to the second smallest 
 *          eigenvalue of the given (symmetric adjacency) matrix.  
 * @author  Michael Holst 
 * @note    Class Bvec: Non-inlineable methods (bsolv.c)\n\n
 *          We make the following assumptions about the matrix A, which
 *          are valid when A is the adjaceny matrix of a simplex mesh
 *          (e.g., the laplace operator for the dual mesh):  \n\n
 *          (1) A is symmetric positive semi-definite\n
 *          (2) A is singular with null space of dimension 1\n
 *          (3) A has a unique smallest eigenvalue "0" and a   
 *              corresponding unique (up to constant multiplier)
 *              eigenvector (or null vector) representing the set of           
 *              constant vectors.\n\n
 *          Note that since A is symmetric, it has a complete set of   
 *          orthonormal eigenvectors. \n\n
 *          The algorithm we employ here is a simple inverse iteration  
 *          (or Rayleigh quotient iteration) which pulls out the  
 *          second smallest eigenvalue (the smallest positive eigenvalue)    
 *          and the corresponding eigenvector.  The algorithm simply  
 *          does a sequence of inexact inverse iterations, while    
 *          repeatedly orthogonalizing the computed eigenvector against
 *          the unique null vector of constants.   The inverse step  
 *          damps out all eigendirections other than the smallest, and
 *          the orthgonalizing step removes the smallest eigen direction.\n\n
 *          On convergence (change in eigenvector meets a tolerance) then 
 *          we generate the corresponding eigenvalue with a final   
 *          Rayleigh quotient.  
 * @return  None
 * @param   thee   pointer to the block vector
 * @param   A      system matrix
 * @param   litmax number of iterations to do (the maximum allowed)
 * @param   letol  error tolerance
 * @param   lambda Point to the eigenvalue
 * @param   key    smooth with A or A' (0=A, 1=A')
 * @param   flag   determines which "mode" we run in:\n
 *                 0 --> Normal: check itmax, error tolerance; normal i/o\n
 *                 1 --> Silent: check only itmax; no i/o\n
 *                 2 --> Subcycle: check itmax, error tolerance; subcycle i/o\n
 *                 3 --> Subcycle: check itmax, error tolerance; no i/o
 * @param   itmax  number of iterations to do (the maximum allowed)
 * @param   etol   error tolerance
 */
VEXTERNC void Bvec_eig(Bvec *thee, Bmat *A,
    int litmax, double letol, double *lambda,
    int key, int flag, int itmax, double etol);

#endif /* _BVEC_H_ */

