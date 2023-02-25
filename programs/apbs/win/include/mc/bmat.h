/**
 * @defgroup Bmat Bmat class
 * @brief    A block sparse matrix object.
 */

/**
 *  @file       bmat.h
 *  @ingroup    Bmat
 *  @brief      Class Bmat: a block sparse matrix object.
 *  @author     Michael Holst
 *  @note       None
 *  @verbatim
 * This class extends the Mat class to block matrices built as
 * retangular arrays of Mat objects.
 *
 * The Bmat class is very efficient in memory and operation
 * complexity for large sparse block matrices having a few
 * (e.g. 1-20) blocks, each of which are themselves large
 * (e.g. > 1000 rows) sparse matrices.  This class IS NOT very
 * efficient for large sparse matrices consisting of many
 * (e.g. > 20) small blocks (e.g. < 1000 rows).  This is because
 * each block is a Mat object which contains some overhead.
 * Note that If some of the matrix blocks are symmetric images
 * of each other, we represent the block only once.
 *
 * Blocks:   The global matrix from a petrov-galerkin FEM discretization of 
 *           a PDE will be stored in blocks; for example, an 18x18 matrix
 *           consisting of 3 blocks of varying sizes (as in the case of a
 *           3-component system, with different number of degrees of freedom
 *           on different components) will be stored as:
 *
 *           A = \----- --- ---------    The block sizes in this example are:
 *               |\---- --- ---------
 *               ||\--- --- ---------        block[0][0] = 6x6
 *               |||\-- --- ---------        block[1][1] = 3x3
 *               ||||\- --- ---------        block[2][2] = 9x9
 *               |||||\ --- ---------        block[0][1] = 6x3
 *                                           block[0][2] = 6x9
 *               |||||| \-- ---------        block[1][2] = 3x9
 *               |||||| |\- ---------        block[1][0] = 3x6
 *               |||||| ||\ ---------        block[2][0] = 9x6
 *                                           block[2][1] = 9x3
 *               |||||| ||| \--------
 *               |||||| ||| |\-------
 *               |||||| ||| ||\------
 *               |||||| ||| |||\-----
 *               |||||| ||| ||||\----
 *               |||||| ||| |||||\---
 *               |||||| ||| ||||||\--
 *               |||||| ||| |||||||\-
 *               |||||| ||| ||||||||\
 *
 *           A prolongation matrix will be stored row-wise in blocks.
 *           For example, a prolongation matrix for a 2-component system,
 *           mapping 3-vectors in the first component to 9-vectors, and
 *           mapping 6-vectors in the second component to 10-vectors, will
 *           be stored in blocks as:
 *
 *           P = --- oooooo    The blocks have the shapes:
 *               --- oooooo
 *               --- oooooo        block[0][0] = 9x3
 *               --- oooooo        block[0][1] = 9x6
 *               --- oooooo        block[1][0] = 10x3
 *               --- oooooo        block[1][1] = 10x6
 *               --- oooooo
 *               --- oooooo
 *               --- oooooo
 *               ooo ------
 *               ooo ------
 *               ooo ------
 *               ooo ------
 *               ooo ------
 *               ooo ------
 *               ooo ------
 *               ooo ------
 *               ooo ------
 *               ooo ------
 *
 *           The adjoint would likely be used as the restriction matrix,
 *           which would then be stored columnwise by simply adjusting
 *           some pointers:
 *
 *           R = ||||||||| oooooooooo    The blocks have the shapes:
 *               ||||||||| oooooooooo
 *               ||||||||| oooooooooo        block[0][0] = 3x9
 *               ooooooooo ||||||||||        block[0][1] = 3x10
 *               ooooooooo ||||||||||        block[1][0] = 6x9
 *               ooooooooo ||||||||||        block[1][1] = 6x10
 *               ooooooooo ||||||||||
 *               ooooooooo ||||||||||
 *               ooooooooo ||||||||||
 *  @endverbatim
 *  @version    $Id: bmat.h,v 1.35 2010/08/12 05:18:34 fetk Exp $ 
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

#ifndef _BMAT_H_
#define _BMAT_H_

#include <mc/mc_base.h>

#include <mc/mat.h>

/**
 * @ingroup Bmat
 * @brief   Contains public data memebers for Bmat class
 * @author  Michael Holst
 */
struct sBmat {

    /** @brief the memory manager */
    Vmem *vmem;    
    /** @brief did i make vmem or was it inherited  */
    int  iMadeVmem;            

    /** @brief Character string name for this matrix */
    char name[10]; 
    /** @brief num of blocks (both row and col) */
    int  numB;

    /** @brief block mirror keys for the matrix:\n
     *  0 => ISNOT (block is stored as a MAT)\n
     *  1 => IS    (block is mirror of a MAT) */
    MATmirror mirror[MAXV][MAXV];  

    /** @brief blocks of this block matrix:\n      
     *  mirror[p][q]=0 ==> AD[p][q] exists\n
     *  mirror[p][q]=1 ==> AD[p][q]->AD[q][p] */
    Mat  *AD[MAXV][MAXV];   

    /** @brief global sparse matrix */
    Mat *AG;                      

    /** @brief next coarser object in the hierarchy */
    struct sBmat *coarse;          
    /** @brief next finer object in the hierarchy */
    struct sBmat *fine; 

};

/**
 * @ingroup Bmat
 * @brief   Declaration of the Bmat class as the Bmat structure
 * @author  Michael Holst
 * @return  None
 */
typedef struct sBmat Bmat;
   
/*
 * ***************************************************************************
 * Class Bmat: Inlineable methods (bmat.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_BAM)
#else /* if defined(VINLINE_BAM) */
#endif /* if !defined(VINLINE_BAM) */

/** 
 * @ingroup Bmat
 * @brief   The block sparse matrix constructor.      
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)\n
 *          This constructor only does minimal initialization of a Bmat 
 *          object after creating the object itself; it doesn't create   
 *          any storage for either the integer structure arrays or the 
 *          nonzero storage arrays.\n
 *          This constructor only fixes the number of blocks and the 
 *          numbers of rows and columns in each block; the nonzero  
 *          structure is not set yet.   
 * @return  Pointer to a newly allocated (empty) block sparse matrix
 * @param   vmem     Memory management object 
 * @param   name     character string name for this matrix
 * @param   pnumB    number of blocks
 * @param   pnumR    num of rows in the matrix
 * @param   pnumC    num of cols (DRC REQUIRES numC=numR)
 * @param   pmirror  the mirror-ness of the block
 */
VEXTERNC Bmat* Bmat_ctor(Vmem *vmem, const char *name,
    int pnumB, int pnumR[MAXV], int pnumC[MAXV],
    MATmirror pmirror[MAXV][MAXV]);

/** 
 * @ingroup Bmat
 * @brief   The block sparse matrix destructor.
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)\n
 *          This destructor does the reverse of Bmat_ctor, and if 
 *          necessary first reverses Bmat_initStructure  
 *          (or Bmat_copyStructure).  I.e., if necessary, 
 *          it first frees the large integer and real arrays created by 
 *          Bmat_initStructure or Bmat_copyStructure, and then frees the 
 *          Bmat object itself at the last moment.  
 * @return  None
 * @param   thee  Pointer to the block sparse matrix 
 */
VEXTERNC void Bmat_dtor(Bmat **thee);

/** 
 * @ingroup Bmat
 * @brief   Initialize the nonzero structure given structure information.
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)\n
 *          This routine actually does the storage creation for both the
 *          integer structure information arrays and the nonzero value  
 *          arrays.
 * @return  None
 * @param   thee     Pointer to the block sparse matrix 
 * @param   pfrmt    format types of the block sparse matrix 
 * @param   psym     symmetric types of the block sparse matrix
 * @param   pnumO    num of nonzeros we are actually storing in the
 *                   strict upper-triangle of matrix. (DRC only)
 * @param   IJA      integer structure [ IA ; JA ]  
 */
VEXTERNC void Bmat_initStructure(Bmat *thee,
    MATformat pfrmt[MAXV][MAXV], MATsym psym[MAXV][MAXV],
    int pnumO[MAXV][MAXV], int *IJA[MAXV][MAXV]);

/** 
 * @ingroup Bmat
 * @brief   Copy the nonzero structure from an input model.
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  None
 * @param   thee     Pointer to the block sparse matrix
 * @param   model    an input model

 */
VEXTERNC void Bmat_copyStructure(Bmat *thee, Bmat *model);

/** 
 * @ingroup Bmat
 * @brief   Kill the nonzero structure.    
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  None
 * @param   thee     Pointer to the block sparse matrix
 */
VEXTERNC void Bmat_killStructure(Bmat *thee);

/** 
 * @ingroup Bmat
 * @brief   Return the number of blocks.     
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the number of blocks
 * @param   thee     Pointer to the block sparse matrix
 */
VEXTERNC int Bmat_numB(Bmat *thee);

/** 
 * @ingroup Bmat
 * @brief   Return the number of rows.     
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the number of rows
 * @param   thee  Pointer to the block sparse matrix
 * @param   p     the index of the blocks
 * @param   q     the index of the blocks 
 */
VEXTERNC int Bmat_numR(Bmat *thee, int p, int q);

/** 
 * @ingroup Bmat
 * @brief   Return the number of columns.               
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the number of columns.
 * @param   thee  Pointer to the block sparse matrix
 * @param   p     the index of the blocks
 * @param   q     the index of the blocks 
 */
VEXTERNC int Bmat_numC(Bmat *thee, int p, int q);

/** 
 * @ingroup Bmat
 * @brief   Return the number of nonzeros we actually store in a block.   
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the number of nonzeros we actually store in a block.   
 * @param   thee  Pointer to the block sparse matrix
 * @param   p     the index of the blocks
 * @param   q     the index of the blocks 
 */
VEXTERNC int Bmat_numA(Bmat *thee, int p, int q);

/** 
 * @ingroup Bmat
 * @brief   Return the number of nonzeros we actually store in a block 
 *          which are actually in the strict upper-triangle of the 
 *          block (DRC only).
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the number of nonzeros we actually store in a block 
 *          which are actually in the strict upper-triangle of the 
 *          block (DRC only).
 * @param   thee  Pointer to the block sparse matrix
 * @param   p     the index of the blocks
 * @param   q     the index of the blocks 
 */
VEXTERNC int Bmat_numO(Bmat *thee, int p, int q);

/** 
 * @ingroup Bmat
 * @brief   Return the number of nonzeros we WOULD actually store in a block
 *          if we ignored symmetry and stored everything.  
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the number of nonzeros we WOULD actually store in a block
 *          if we ignored symmetry and stored everything.  
 * @param   thee  Pointer to the block sparse matrix
 * @param   p     the index of the blocks
 * @param   q     the index of the blocks 
 */
VEXTERNC int Bmat_numZ(Bmat *thee, int p, int q);

/** 
 * @ingroup Bmat
 * @brief   Return the total number of rows.    
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the total number of rows
 * @param   thee  Pointer to the block sparse matrix
 */
VEXTERNC int Bmat_numRT(Bmat *thee);

/** 
 * @ingroup Bmat
 * @brief   Return the total number of columns.                 
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the total number of columns
 * @param   thee  Pointer to the block sparse matrix
 */
VEXTERNC int Bmat_numCT(Bmat *thee);

/** 
 * @ingroup Bmat
 * @brief   Return the total number of nonzeros we are actually storing. 
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the total number of nonzeros we are actually storing. 
 * @param   thee  Pointer to the block sparse matrix
 */
VEXTERNC int Bmat_numAT(Bmat *thee);

/** 
 * @ingroup Bmat
 * @brief   Return the total number of nonzeros we are actually storing
 *          which are located in the strict upper-triangle. 
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the total number of nonzeros we are actually storing
 *          which are located in the strict upper-triangle. 
 * @param   thee  Pointer to the block sparse matrix
 */
VEXTERNC int Bmat_numOT(Bmat *thee);

/** 
 * @ingroup Bmat
 * @brief   Return the total number of nonzeros we WOULD be storing if we 
 *          ignored all symmetry in all blocks and stored everything.  
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the total number of nonzeros we WOULD be storing if we 
 *          ignored all symmetry in all blocks and stored everything.  
 * @param   thee  Pointer to the block sparse matrix
 */
VEXTERNC int Bmat_numZT(Bmat *thee);

/** 
 * @ingroup Bmat
 * @brief   Return the format of the block.       
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the format of the block.
 * @param   thee  Pointer to the block sparse matrix
 * @param   p     the index of the blocks
 * @param   q     the index of the blocks 
 */
VEXTERNC MATformat Bmat_format(Bmat *thee, int p, int q);

/** 
 * @ingroup Bmat
 * @brief   Return the symmetry of the block.     
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the symmetry of the block.     
 * @param   thee  Pointer to the block sparse matrix
 * @param   p     the index of the blocks
 * @param   q     the index of the blocks 
 */
VEXTERNC MATsym Bmat_sym(Bmat *thee, int p, int q);

/** 
 * @ingroup Bmat
 * @brief   Return the state of the block.    
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the state of the block.    
 * @param   thee  Pointer to the block sparse matrix
 * @param   p     the index of the blocks
 * @param   q     the index of the blocks 
 */
VEXTERNC MATstate Bmat_state(Bmat *thee, int p, int q);

/** 
 * @ingroup Bmat
 * @brief   Return the implicitness of the block.   
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the implicitness of the block.   
 * @param   thee  Pointer to the block sparse matrix
 * @param   p     the index of the blocks
 * @param   q     the index of the blocks 
 */
VEXTERNC MATimpl Bmat_impl(Bmat *thee, int p, int q);

/** 
 * @ingroup Bmat
 * @brief   Return the mirror-ness of the block.          
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the mirror-ness of the block.          
 * @param   thee  Pointer to the block sparse matrix
 * @param   p     the index of the blocks
 * @param   q     the index of the blocks 
 */
VEXTERNC MATmirror Bmat_mirror(Bmat *thee, int p, int q);

/** 
 * @ingroup Bmat
 * @brief   Return the number of nonzeros in all blocks.       
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the number of nonzeros in all blocks.
 * @param   thee  Pointer to the block sparse matrix
 */
VEXTERNC int Bmat_sizeA(Bmat *thee);

/** 
 * @ingroup Bmat
 * @brief   Return the numer of integer storage locations used.   
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the numer of integer storage locations used.   
 * @param   thee  Pointer to the block sparse matrix
 */
VEXTERNC int Bmat_sizeIJA(Bmat *thee);

/** 
 * @ingroup Bmat
 * @brief   Return the integer structure IJA.  
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the integer structure IJA.  
 * @param   thee  Pointer to the block sparse matrix
 * @param   p     the index of the blocks
 * @param   q     the index of the blocks 
 */
VEXTERNC int *Bmat_IJA(Bmat *thee, int p, int q);

/** 
 * @ingroup Bmat
 * @brief   Return the integer structure IA.        
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the integer structure IA.  
 * @param   thee  Pointer to the block sparse matrix
 * @param   p     the index of the blocks
 * @param   q     the index of the blocks 
 */
VEXTERNC int *Bmat_IA(Bmat *thee, int p, int q);

/** 
 * @ingroup Bmat
 * @brief   Return the integer structure JA.        
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the integer structure JA.  
 * @param   thee  Pointer to the block sparse matrix
 * @param   p     the index of the blocks
 * @param   q     the index of the blocks 
 */
VEXTERNC int *Bmat_JA(Bmat *thee, int p, int q);

/** 
 * @ingroup Bmat
 * @brief   Return the real structure A.  
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the real structure 
 * @param   thee  Pointer to the block sparse matrix
 * @param   p     the index of the blocks
 * @param   q     the index of the blocks 
 */
VEXTERNC double *Bmat_A(Bmat *thee, int p, int q);

/** 
 * @ingroup Bmat
 * @brief   Return the diagonal of A (DRC only).  
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the diagonal of A (DRC only).  
 * @param   thee  Pointer to the block sparse matrix
 * @param   p     the index of the blocks
 * @param   q     the index of the blocks 
 */
VEXTERNC double *Bmat_diag(Bmat *thee, int p, int q);

/** 
 * @ingroup Bmat
 * @brief   Return the strict upper triangle of A (DRC only).  
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the strict upper triangle of A (DRC only).  
 * @param   thee  Pointer to the block sparse matrix
 * @param   p     the index of the blocks
 * @param   q     the index of the blocks 
 */
VEXTERNC double *Bmat_offU(Bmat *thee, int p, int q);

/** 
 * @ingroup Bmat
 * @brief   Return the strict lower triangle of A (DRC only).  
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  the strict lower triangle of A (DRC only).  
 * @param   thee  Pointer to the block sparse matrix
 * @param   p     the index of the blocks
 * @param   q     the index of the blocks
 */
VEXTERNC double *Bmat_offL(Bmat *thee, int p, int q);

/** 
 * @ingroup Bmat
 * @brief   Print the prolongation matrix blocks.   
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  None
 * @param   thee  Pointer to the block sparse matrix
 */
VEXTERNC void Bmat_print(Bmat *thee);

/** 
 * @ingroup Bmat
 * @brief   Print the prolongation matrix in MATLAB sparse form.   
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  None
 * @param   thee    Pointer to the block sparse matrix
 * @param   fname   output file/buff/unix/inet name
 * @param   pflag   flag==0 ==> write, flag==1 ==> append
 */
VEXTERNC void Bmat_printSp(Bmat *thee, char *fname, int pflag);

/** 
 * @ingroup Bmat
 * @brief   Print the matrix as a DENSE matrix in MATLAB format, 
 *          but first zero out any rows/cols corresponding to  
 *          Dirichlet boundary points.     
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)\n
 *          This routine is useful for e.g. checking that Galerkin
 *          conditions hold for stiffness matrices.  Removing the
 *          dirichlet equations is crucial; otherwise the Galerkin 
 *          condition cannot hold.  Note that the matrix (and the  
 *          Galerkin coarse matrix) are then of course singular.  
 * @return  None
 * @param   thee    Pointer to the block sparse matrix
 */
VEXTERNC void Bmat_printNoD(Bmat *thee);

/** 
 * @ingroup Bmat
 * @brief   Print the matrix as a DENSE matrix in MATLAB format, 
 *          but first zero out any rows/cols corresponding to  
 *          Dirichlet boundary points.     
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)\n
 *          This routine is useful for e.g. checking that Galerkin   
 *          conditions hold for stiffness matrices.  Removing the
 *          dirichlet equations is crucial; otherwise the Galerkin
 *          condition cannot hold.  Note that the matrix (and the
 *          Galerkin coarse matrix) are then of course singular. 
 * @return  None
 * @param   thee    Pointer to the block sparse matrix
 * @param   fname   output file/buff/unix/inet name
 * @param   pflag   index for write/append
 */
VEXTERNC void Bmat_printSpNoD(Bmat *thee, char *fname, int pflag);

/** 
 * @ingroup Bmat
 * @brief   Clear the floating point storage for the sparse matrix.
 *          Also clear any sparse factorization storage.  
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)\n
 *          This is basically done in preparation for an accumulation as  
 *          part of a matrix assembly, and before a new sparse factorization.
 * @return  None
 * @param   thee    Pointer to the block sparse matrix
 */
VEXTERNC void Bmat_zero(Bmat *thee);

/** 
 * @ingroup Bmat
 * @brief   Setup the dirichlet equations. 
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  None
 * @param   thee    Pointer to the block sparse matrix
 */
VEXTERNC void Bmat_diri(Bmat *thee);

/** 
 * @ingroup Bmat
 * @brief   Set the (i,j)-th entry of the (p,q)-th block to < val >      
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  None
 * @param   thee  Pointer to the block sparse matrix
 * @param   p     the index of the blocks
 * @param   q     the index of the blocks 
 * @param   i     the index of the matrix
 * @param   j     the index of the matrix
 * @param   val   the value of the matrix element
 */
VEXTERNC void Bmat_set(Bmat *thee, int p, int q, int i, int j, double val);

/** 
 * @ingroup Bmat
 * @brief   Contribute < val > to the (i,j)-th entry of the (p,q)-th block.
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  None
 * @param   thee  Pointer to the block sparse matrix
 * @param   p     the index of the blocks
 * @param   q     the index of the blocks 
 * @param   i     the index of the matrix
 * @param   j     the index of the matrix
 * @param   val   the value of the matrix element
 */
VEXTERNC void Bmat_addTo(Bmat *thee, int p, int q, int i, int j, double val);

/** 
 * @ingroup Bmat
 * @brief   Enforce the Galerkin conditions algebraically.       
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  None
 * @param   thee    Pointer to the block sparse matrix
 * @param   rmat    R matrix which is stored column-wise (ROW-format)
 * @param   amat    A matrix which is stored in one of three
 *           forms, namely, either row-wise (ROW), col-wise (COL), or by
 *           diagonal followed by upper-triangle row-wise and then by lower
 *           triangle columne-wise (DRC).
 * @param   pmat    P matrix which is stored row-wise (ROW-format)
 */
VEXTERNC void Bmat_galerkin(Bmat *thee, Bmat *rmat, Bmat *amat, Bmat *pmat);

/** 
 * @ingroup Bmat
 * @brief   Make a decision about whether or not a sparse direct solver  
 *          should be used in place of an iterative solver, based on the 
 *          size of the system.
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)\n
 *          This is obviously heuristic in nature; in general the cutoff 
 *          size where iterative methods start to win is smaller in 3D.
 * @return  the decision about whether or not a sparse direct solver
 *           should be used in place of an iterative solver, based on the
 *           size of the system.
 * @param   thee    Pointer to the block sparse matrix
 */
VEXTERNC int Bmat_sluDirect(Bmat *thee);

/** 
 * @ingroup Bmat
 * @brief   Create the global matrix from the blocks in the modified      
 *          ROW or COL storage format.  This is useful for preparing a       
 *          single global matrix for input to e.g. a sparse direct solver. 
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  None
 * @param   thee    Pointer to the block sparse matrix
 */
VEXTERNC void Bmat_sluCreate(Bmat *thee);

/** 
 * @ingroup Bmat
 * @brief   Create the sparse LU factors for global matrix. 
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  Sparse LU factor
 * @param   thee    Pointer to the block sparse matrix
 */
VEXTERNC int Bmat_sluFactor(Bmat *thee);

/** 
 * @ingroup Bmat
 * @brief   Forward/backward solve using sparse LU factors of global matrix.  
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)\n
 *          This requires that Bmat_sluFactor has been previously called. 
 * @return  Success enumeration
 * @param   thee  Pointer to the block sparse matrix
 * @param   key   index for choosing NOTRANS or TRANS
 * @param   f     the number of right-hand sides
 * @param   u     solution 
 */
VEXTERNC int Bmat_sluSolve(Bmat *thee, int key, double *f, double *u);

/** 
 * @ingroup Bmat
 * @brief   Destroy the sparse LU factors for the system matrix.      
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)\n
 *          This frees our < ia,ja,a > storage, and also the internal
 *          storage that was malloc'd by the sparse direct library.     
 * @return  None
 * @param   thee  Pointer to the block sparse matrix
 */
VEXTERNC void Bmat_sluDestroy(Bmat *thee);

/** 
 * @ingroup Bmat
 * @brief   Print the exact current malloc usage.         
 * @author  Michael Holst
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  None
 * @param   thee  Pointer to the block sparse matrix
 */
VEXTERNC void Bmat_memChk(Bmat *thee);

/** 
 * @ingroup Bmat
 * @brief   Construct a clone of a Bmat with the same structure.
 * @author  Stephen Bond
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  None
 * @param   vmem  Memory management object
 * @param   name  character string name for this matrix
 * @param   X     The matrix which is cloned
 */
VEXTERNC Bmat *Bmat_clone(Vmem *vmem, char *name, Bmat *X);

/** 
 * @ingroup Bmat
 * @brief   Copy a block matrix.
 * @author  Stephen Bond
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  None
 * @param   Y     The new matrix
 * @param   X     The old matrix
 */
VEXTERNC void Bmat_copy(Bmat *Y, Bmat *X);

/** 
 * @ingroup Bmat
 * @brief   Remove the boundary rows or columns from a block matrix.
 * @author  Stephen Bond
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  None
 * @param   thee  Pointer to the block sparse matrix
 * @param   key   index for removing the boundary rows or columns from a matrix
 */
VEXTERNC void Bmat_squeezeBRC(Bmat *thee, int key);

/** 
 * @ingroup Bmat
 * @brief   Copy a block matrix.
 * @author  Stephen Bond
 * @note    Class Bmat: Non-inlineable methods (bmat.c)
 * @return  None
 * @param   Y     The new matrix
 * @param   X     The old matrix
 */
VEXTERNC void Bmat_copy2(Bmat *Y, Bmat *X);

/** 
 * @ingroup Bmat
 * @brief   Scalar times a Bmat plus a Bmat:  Y += val*X. 
 * @author  Stephen Bond
 * @note    Class Bmat: Non-inlineable methods (bmat.c)\n
 *          The function of this routine can be controlled using "key"   
 * @return  None
 * @param   Y     The new matrix
 * @param   X     The old matrix
 * @param   val   the coeficient for timing X matrix
 * @param   key   index for different X and Y matrices.
 */
VEXTERNC void Bmat_axpy(Bmat *Y, Bmat *X, double val, int key);

#endif /* _BMAT_H_ */

