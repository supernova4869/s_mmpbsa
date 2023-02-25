/** 
 * @defgroup Mat Mat class
 * @brief    A Spare Matrix object.
 */

/**
 *  @file       mat.h
 *  @ingroup    Mat
 *  @brief      Class Mat: a sparse matrix object.
 *  @author     Michael Holst
 *  @note       None
 *  @verbatim
 * This class support several datastructures for sparse matrices.
 *           The following formats are supported (see below for descriptions):
 *
 *             ZERO (the zero matrix; no storage at all)
 *             DRC  (diag-row-col-YSMP variant)
 *             ROW  (row-YSMP)
 *             COL  (col-YSMP)
 *             SLU  (sparse-LU)
 *             RLN  (row-linked-list)
 *             CLN  (col-linked-list)
 *             XLN  (row-linked-list AND col-linked-list)
 *             RFL  (row-wise dense matrix; i.e., C-style 2D array)
 *             CFL  (col-wise dense matrix; i.e., FORTRAN-style 2D array)
 *
 *           NOTE: This class is very efficient in both memory and operation
 *           complexity for LARGE sparse matrices.  It also supports dense
 *           matrices (RFL and CFL formats).
 *
 *           ANOTHER NOTE:  Only ONE format is supported at a time; i.e.,
 *           you cannot simultaneously have a matrix in multiple formats
 *           maintained WITHIN that Mat datastructure.  You must create a
 *           new one of the desired type, and copy the old one into it, if
 *           you want to have an existing matrix represented in a different
 *           format.  (This is a change from the previous version of this
 *           library.)
 *
 * Formats:  Here is a brief description of the suppored matrix formats:
 *
 *           ZERO:       .                  This is a zero matrix (no storage).
 *
 *           DRC:        \-----             This is a symmetric storage
 *                       |\----             format; only the upper triangle
 *                       ||\---             need be stored in the case of
 *                       |||\--             symmetry.  However, we must assume
 *                       ||||\-             that the upper and lower triangles
 *                       |||||\             have identical nonzero structures;
 *                                          the matrix MUST BE SQUARE.
 *                                          The diagonal entries are stored
 *                                          separately from the triangles.
 *
 *                                          NOTE: In the case of symmetry,
 *                                          we simply point the lower triangle
 *                                          nonzeros A to the upper, as well as
 *                                          the IA and JA pointers.
 *
 *           ROW:        ---------          This is a completely nonsymmetric
 *                       ---------          storage format; no provision is
 *                       ---------          made to handle storage savings
 *                       ---------          in the case of symmetry.
 *                                          No assumptions are made about the
 *                                          nonzero-structure of the matrix;
 *                                          the matrix can be non-square.
 *                                          Row-start pointers are kept in IA,
 *                                          and column indices are kept in JA.
 *                                          The diagonal entriecs are treated
 *                                          like any other row entry.
 *
 *           COL:        ||||               This is a column-wise variant
 *                       ||||               of the ROW format.  Col-start
 *                       ||||               pointers are kept in IA, and
 *                       ||||               row indices are kept in JA.
 *                       ||||
 *                       ||||
 *                       ||||
 *                       ||||
 *                       ||||
 *
 *           SLU:        [LU]               This format is determined by
 *                                          the particular sparse direct
 *                                          solver which we use to factor
 *                                          the matrix, and we do not use
 *                                          any information about the
 *                                          particular storage format.
 *
 *           RLN:        [linked-list]      This is a row-wise linked list
 *                                          representation of the matrix.
 *                                          It is usually used to accumulate
 *                                          a matrix product for which there
 *                                          is no a priori knowledge about
 *                                          the resulting nonzero structure.
 *                                          It is usually converted into one
 *                                          of the other matrix formats before
 *                                          it is used for anything else, since
 *                                          linked-list implementations of
 *                                          operations such as matrix-vector
 *                                          products tend to be inefficient.
 *
 *           CLN:        [linked-list]      This is a column-wise variant
 *                                          of the RLN format.
 *
 *           XLN:        [linked-list]      Simultaneous RLN and CLN.
 *
 *           RFL:        ---------          This is a dense row-wise storage
 *                       ---------          format.  All nonzeros are stored;
 *                       ---------          no integer storage is used.
 *
 *           CFL:        ||||               This is a dense col-wise storage
 *                       ||||               format.  All nonzeros are stored;
 *                       ||||               no integer storage is used.
 *
 * Storage:  Details of the DRC/ROW/COL format storage are as follows.
 *           The integer part IJA of the structure has the following layout:
 *
 *               IJA = [ IA ; JA ]
 *
 *               length(IA)    = N+1      row(col) start pointers into JA/A
 *               length(JA)    = NZ       col(row) indices for each row(col)
 *               ------------------------
 *               length(IJA)   = N+1+NZ
 *
 *           which is a DRC, ROW, or COL pointer structure for a matrix.
 *           The IA part of the array points into the JA portion of the array
 *           for the beginning of each row (or column) left-to-right
 *           (or top-to-bottom).
 *
 *           The JA portion then contains the column (or row) indices for the
 *           corresponding entries in that row (or column), ordered
 *           left-to-right (or top to bottom).
 *
 *           The corresponding array of the actual nonzeros in the case of
 *           the DRC style format has the form:
 *
 *               A = [ diag ; offU ; offL ]
 *
 *               length(diag) = N         (diag)
 *               length(offU) = NZ        (upper-triang)
 *               length(offL) = NZ        (lower-triang, or null if symmetric)
 *               -----------------------
 *               length(A)    = N+2*NZ    (or N+NZ if symmetric)
 *
 *           In the case of ROW or COL, the matrix A is simply a row-wise or
 *           col-wise ordering of all of the nozeros, with no special role
 *           played by the diagonal entries.
 *  @endverbatim
 *  @version    $Id: mat.h,v 1.47 2010/08/12 05:18:35 fetk Exp $ 
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


#ifndef _MAT_H_
#define _MAT_H_

#include <mc/mc_base.h>

#include <mc/mtool.h>
#include <mc/slu.h>


/**
 * @ingroup Mat
 * @brief   Contains public data memebers for Mat class
 * @author  Michael Holst
 */
struct sMat {

    /** @brief SETUP (name, memory management, etc) 
     *         character string name for this matrix */
    char   name[10];  
    /** @brief the memory manager */ 
    Vmem   *vmem;     
    /** @brief did i make vmem or was it inherited */
    int    iMadeVmem;  

    /**
     *  @brief PARAMETERS (storage format, symmetry, etc) 
     *  @note  possible format types of this matrix:             \n
     *         0 => ZERO (zero matrix; no structure at all)      \n
     *         1 => DRC  (diag-row-col-YSMP variant)             \n
     *         2 => ROW  (row-YSMP)                              \n
     *         3 => COL  (col-YSMP)                              \n
     *         4 => SLU  (sparse-LU)                             \n
     *         5 => RLN  (row-linked-list)                       \n
     *         6 => CLN  (col-linked-list)                       \n
     *         7 => XLN  (row-linked-list AND col-linked-list)   \n
     *         9 => RFL  (row-wise dense matrix)                 \n
     *         10 => CFL  (col-wise dense matrix) 
     */
    MATformat format;
 
    /** 
     * @brief possible states of this matrix format     
     * @note  0 => NULL (exists; shape fixed; no storage)       \n
     *        1 => ZERO (matrix storage available and zeroed)   \n
     *        2 => ASSEMBLED (matrix ready to use)              \n
     *        3 => FACTORED (matrix factored) (SLU only) 
     */
    MATstate state;    

    /**
     *  @brief symmetry keys for the matrix
     *  @note  0 => ISNOT (store all)                           \n
     *         1 => IS    (store upper tri) (DRC only)          \n
     *         2 => STRUC (store all, but can reuse int structs)
     */
    MATsym sym;        

    /**
     *  @brief implicit diagonal block not stored.
     *  @note  0 => ISNOT (nothing implicit; everything stored) \n
     *         1 => IS    (implicit diagonal block not stored)  \n
     *                                                          \n\n
     *         If impl==1, then an implicit diagonal block is   
     *         present but not stored; this is to avoid the     
     *         extra log N-like storage requirement for storing  
     *         the identity block in all prolongation matrices,  
     *         but yet account for their impact when multiplying 
     *         by prolongation and restriction matrices.         
     *                                                          \n\n
     *         The setting (impl==1) is valid ONLY for ROW and   
     *         COL formats.  Moreover, there are restrictions in
     *         the dimensions of these two formats to support   
     *         the implicit identity; these are:                 
     *                                                          \n
     *             ROW:   numR >= numC  (i.e., tall and skinny)  \n
     *             COL:   numR <= numC  (i.e., short and fat)   
     */
    MATimpl impl;   

    /** @brief DIMENSIONS (row and col dimensions, nonzeros, etc)  
     *         num of rows in the matrix */
    int    numR;   
    /** @brief DIMENSIONS (row and col dimensions, nonzeros, etc)  
     *         num of cols (DRC REQUIRES numC=numR) */
    int    numC;  
    /** @brief DIMENSIONS (row and col dimensions, nonzeros, etc)  
     *         num of nonzeros we are actually storing, counting 
     *         the diagonal, the strict upper-triangle, and  
     *         also the strict lower-triangle if we are       
     *         actually storing the lower-triangle (sym=0).   */
    int    numA;       
    /** @brief DIMENSIONS (row and col dimensions, nonzeros, etc)  
     *         num of nonzeros we are actually storing in the   
     *         strict upper-triangle of matrix. (DRC only)    */
    int    numO;      
    /**
     *  @brief DIMENSIONS (row and col dimensions, nonzeros, etc)  
     *         num of nonzeros we WOULD be storing if we ignored 
     *         symmetry. (DRC only).                           
     *  @note  The relationships between numZ/numA and numO are: \n
     *         non-DRC:           numO = numZ = numA             \n
     *  DRC-symmetric:     numA = numR + numO, numZ = numA + numO\n
     *  DRC-non-symmetric: numA = numR + 2*numO,  numZ = numA    \n
     */
    int    numZ; 
    /** @brief DIMENSIONS (row and col dimensions, nonzeros, etc)  
     *         num of boundary rows */
    int    numBR;    
    /** @brief DIMENSIONS (row and col dimensions, nonzeros, etc)  
     *          num of boundary cols */
    int    numBC;      

     /** @brief MALLOC AREAS (high-order storage). Did I malloc IJA?  */
    int    iMallocIJA;
    /** @brief MALLOC AREAS (high-order storage).  Did I malloc A? */
    int    iMallocA;   
    /** @brief MALLOC AREAS (high-order storage)
     *         integer structure [ IA ; JA ] */
    int    *IJA;       
    /** @brief MALLOC AREAS (high-order storage)    \n
     *         packed nozeros:                      \n
     *         DRC: [ diag ; offU ; offL ]          \n        
     *         ROW: [ offU ]                        \n         
     *         COL: [ offL ]                        \n         
     *         RFL: [ everything; stored row-wise ] \n         
     *         CFL: [ everything; stored col-wise ] */
    double *A;    
    /** @brief MALLOC AREAS (high-order storage) 
     *         boundary rows (optionally used) */
    int    *BR;   
    /** @brief MALLOC AREAS (high-order storage) 
     *         boundary cols (optionally used) */
    int    *BC; 

    /** @brief ALIASES (pointers into the above malloc areas; low-order storage)
     *         pos in JA/offU/offL for row/col start */
    int    *IA;  
    /** @brief ALIASES (pointers into the above malloc areas; low-order storage)
    *         row/col indices for nonzeros in col/row */
    int    *JA;   
    /** @brief ALIASES (pointers into the above malloc areas; low-order storage)
     *         diagonal of the matrix */
    double *diag;  
    /** @brief ALIASES (pointers into the above malloc areas; low-order storage)
     *         upper-triangle off-diag row-wise nonzeros  */
    double *offU; 
    /** @brief ALIASES (pointers into the above malloc areas; low-order storage)
     *         lower-triangle off-diag col-wise nonzeros */
    double *offL; 

    /** @brief EXTERNAL SUPPORT (handles to other complex objects)  
     *         Sparse LU factorization container object */
    Slu    *slu;
    /** @brief EXTERNAL SUPPORT (handles to other complex objects)
     *         Vset object for linked list utilities and RLN */
    Vset   *lnkL;   
    /** @brief EXTERNAL SUPPORT (handles to other complex objects)
     *         Vset object for linked list utilities and CLN  */
    Vset   *lnkU;
    /** @brief EXTERNAL SUPPORT (handles to other complex objects)       
     *         Support for XLN */
    void   *xln;   
    /** @brief EXTERNAL SUPPORT (handles to other complex objects)  
     *         Support for XLN */
    void   *xlnt;   

};

/**
 * @ingroup Mat
 * @brief   Declaration of the Mat class as the Mat structure
 * @author  Michael Holst
 * @return  None
 */
typedef struct sMat Mat;

/*
 * ***************************************************************************
 * Class Mat: Inlineable methods (mat.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_BAM)
#else /* if defined(VINLINE_BAM) */
#endif /* if !defined(VINLINE_BAM) */

/** 
 * @ingroup Mat
 * @brief   The sparse matrix constructor.   
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c) \n
 *          This constructor only fixes the number of rows and columns 
 *          in the matrix; the nonzero structure is not set.  
 * @return  Pointer to a newly allocated (empty) sparse matrix
 * @param   vmem   Memory management object
 * @param   name   character string name for this matrix
 * @param   pnumR  num of rows in the matrix
 * @param   pnumC  num of cols (DRC REQUIRES numC=numR)
 */
VEXTERNC Mat* Mat_ctor(Vmem *vmem, const char *name, int pnumR, int pnumC);

/** 
 * @ingroup Mat
 * @brief   The sparse matrix destructor. 
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c) \n
 *          This destructor does the reverse of Mat_ctor, and if 
 *          necessary first reverses Mat_initStructure 
 *          (or Mat_copyStructure).  I.e., if necessary,
 *          it first frees the large integer and real arrays created   
 *          by Mat_initStructure or Mat_copyStructure, and then frees 
 *          the Mat object itself at the last moment.     
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 */
VEXTERNC void Mat_dtor(Mat **thee);

/** 
 * @ingroup Mat
 * @brief   Initialize the nonzero structure given structure information. 
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c) \n
 *          This routine actually does the storage creation for both the 
 *          integer structure information arrays and the nonzero value  
 *          arrays.
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   frmt  format types of the sparse matrix
 * @param   sym   symmetric types of the sparse matrix
 * @param   numO  num of nonzeros we are actually storing in the
 *                   strict upper-triangle of matrix. (DRC only)
 * @param   IJA   integer structure [ IA ; JA ]
 * @param   A     packed nozeros
 */
VEXTERNC void Mat_initStructure(Mat *thee,
    MATformat frmt, MATsym sym, int numO, int *IJA, double *A);

/** 
 * @ingroup Mat
 * @brief   Initialize the nonzero structure given structure information.  
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c) \n
 *          This routine actually does the storage creation for both the
 *          integer structure information arrays and the nonzero value    
 *          arrays.
 * @return  None
 * @param   thee   Pointer to the sparse matrix
 * @param   model  an input matrix
 */
VEXTERNC void Mat_copyStructure(Mat *thee, Mat *model);

/** 
 * @ingroup Mat
 * @brief   Kill the nonzero structure and structure information.           
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c) \n
 *          This routine does the reverse of Mat_initStructure 
 *          (or Mat_copyStructure).  It leaves only the information 
 *          about the number of blocks, number of rows, and number of
 *          columns per block.  I.e., what is left is only what was 
 *          present after the initial call to Mat_ctor.                
 * @return  None
 * @param   thee   Pointer to the sparse matrix
 */
VEXTERNC void Mat_killStructure(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Return number of rows in the matrix
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  number of rows in the matrix
 * @param   thee   Pointer to the sparse matrix
 */
VEXTERNC int Mat_numR(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Return number of columns in the matrix
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  number of columns in the matrix
 * @param   thee   Pointer to the sparse matrix
 */
VEXTERNC int Mat_numC(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Return number of nonzeros we are actually storing.           
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  number of nonzeros we are actually storing.           
 * @param   thee   Pointer to the sparse matrix
 */
VEXTERNC int Mat_numA(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Return number of nonzeros we are actually storing 
 *          which are located in upper (or lower) triangle. 
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  number of nonzeros we are actually storing 
 *          which are located in upper (or lower) triangle. 
 * @param   thee   Pointer to the sparse matrix
 */
VEXTERNC int Mat_numO(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Return number of nonzeros we WOULD be storing if we were 
 *          ignoring symmetry and storing all nonzeros.      
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  number of nonzeros we WOULD be storing if we were 
 *          ignoring symmetry and storing all nonzeros.      
 * @param   thee   Pointer to the sparse matrix
 */
VEXTERNC int Mat_numZ(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Return the format
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  the format
 * @param   thee   Pointer to the sparse matrix
 */
VEXTERNC MATformat Mat_format(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Return the symmetry.
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  the symmetry.
 * @param   thee   Pointer to the sparse matrix
 */
VEXTERNC MATsym Mat_sym(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Return the state
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  the state
 * @param   thee   Pointer to the sparse matrix
 */
VEXTERNC MATstate Mat_state(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Return the impl.
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  the impl
 * @param   thee   Pointer to the sparse matrix
 */
VEXTERNC MATimpl Mat_impl(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Set the format
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  None
 * @param   thee    Pointer to the sparse matrix
 * @param   format  the sparse matrix format
 */
VEXTERNC void Mat_setFormat(Mat *thee, MATformat format);

/** 
 * @ingroup Mat
 * @brief   Set the symmetry
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   sym   symmetric types of the sparse matrix
 */
VEXTERNC void Mat_setSym(Mat *thee, MATsym sym);

/** 
 * @ingroup Mat
 * @brief   Set the state
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  None
 * @param   thee   Pointer to the sparse matrix
 * @param   state  possible states of the sparse matrix format
 */
VEXTERNC void Mat_setState(Mat *thee, MATstate state);

/** 
 * @ingroup Mat
 * @brief   Set the impl
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   impl  implicit diagonal block not stored
 */
VEXTERNC void Mat_setImpl(Mat *thee, MATimpl impl);

/** 
 * @ingroup Mat
 * @brief   Return total number of INTEGER STORAGE LOCATIONS in the matrix.
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  total number of INTEGER STORAGE LOCATIONS in the matrix.
 * @param   thee  Pointer to the sparse matrix
 */
VEXTERNC int Mat_sizeIJA(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Return total number of REAL STORAGE LOCATIONS in the matrix.  
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  total number of REAL STORAGE LOCATIONS in the matrix.  
 * @param   thee  Pointer to the sparse matrix
 */
VEXTERNC int Mat_sizeA(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Return the integer structure IJA.              
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  the integer structure IJA.              
 * @param   thee  Pointer to the sparse matrix
 */
VEXTERNC int *Mat_IJA(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Return the integer structure IA.    
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  the integer structure IA.    
 * @param   thee  Pointer to the sparse matrix
 */
VEXTERNC int *Mat_IA(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Return the integer structure JA.    
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  the integer structure JA.    
 * @param   thee  Pointer to the sparse matrix
 */
VEXTERNC int *Mat_JA(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Return the real structure A.    
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  the real structure A
 * @param   thee  Pointer to the sparse matrix
 */
VEXTERNC double *Mat_A(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Return the diagonal structure A. 
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  the diagonal structure A. 
 * @param   thee  Pointer to the sparse matrix
 */
VEXTERNC double *Mat_diag(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Return the upper-triangle structure A.              
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  the upper-triangle structure A.              
 * @param   thee  Pointer to the sparse matrix
 */
VEXTERNC double *Mat_offU(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Return the lower-triangle structure A.              
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  the lower-triangle structure A.              
 * @param   thee  Pointer to the sparse matrix
 */
VEXTERNC double *Mat_offL(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Print the matrix as a DENSE matrix in MATLAB format.  
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 */
VEXTERNC void Mat_print(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Print the matrix as a SPARSE matrix in MATLAB format.
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  None
 * @param   thee   Pointer to the sparse matrix
 * @param   fname  the output sparse matrix name
 * @param   pflag  0 ==> write, 1 ==> append
 */
VEXTERNC void Mat_printSp(Mat *thee, char *fname, int pflag);

/** 
 * @ingroup Mat
 * @brief   Print the matrix as a DENSE matrix in MATLAB format, 
 *          but first zero out any rows/cols corresponding to 
 *          Dirichlet boundary points. 
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)\n
 *          This routine is useful for e.g. checking that Galerkin
 *          conditions hold for stiffness matrices.  Removing the 
 *          dirichlet equations is crucial; otherwise the Galerkin
 *          condition cannot hold.  Note that the matrix (and the
 *          Galerkin coarse matrix) are then of course singular.
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 */
VEXTERNC void Mat_printNoD(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Print the matrix as a DENSE matrix in MATLAB format, 
 *          but first zero out any rows/cols corresponding to 
 *          Dirichlet boundary points. 
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)\n
 *          This routine is useful for e.g. checking that Galerkin
 *          conditions hold for stiffness matrices.  Removing the 
 *          dirichlet equations is crucial; otherwise the Galerkin
 *          condition cannot hold.  Note that the matrix (and the
 *          Galerkin coarse matrix) are then of course singular.
 * @return  None
 * @param   thee   Pointer to the sparse matrix
 * @param   fname  the output sparse matrix name
 * @param   pflag  index for write/append
 */
VEXTERNC void Mat_printSpNoD(Mat *thee, char *fname, int pflag);

/** 
 * @ingroup Mat
 * @brief   Clear the floating point storage for the sparse matrix. 
 *          Also clear any sparse factorization storage.  
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)\n
 *          This is basically done in preparation for an accumulation as
 *          part of a matrix assembly, and before a new sparse factorization. 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 */
VEXTERNC void Mat_zero(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Set a value in a matrix.
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   i     the index for row
 * @param   j     the index for column
 * @param   val   the value of the sparse matrix element
 */
VEXTERNC void Mat_set(Mat *thee, int i, int j, double val);

/** 
 * @ingroup Mat
 * @brief   Add to a value in a matrix.
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   i     the index for row
 * @param   j     the index for column
 * @param   val   the value to be added on the sparse matrix element
 */
VEXTERNC void Mat_addTo(Mat *thee, int i, int j, double val);

/** 
 * @ingroup Mat
 * @brief   Set the boundary row and column information.    
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)\n
 *          key=0 ==> set pointers, key=1 ==> do a copy
 * @return  None
 * @param   thee   Pointer to the sparse matrix
 * @param   numBR  num of boundary rows
 * @param   numBC  num of boundary columns
 * @param   BR     boundary rows (optionally used)
 * @param   BC     boundary columns
 */
VEXTERNC void Mat_buildBRC(Mat *thee, int numBR, int numBC, int *BR, int *BC);

/** 
 * @ingroup Mat
 * @brief   Apply the boundary row and column information.   
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 */
VEXTERNC void Mat_zeroBRC(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Place identity entries on diagonal of boundary row/col.
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 */
VEXTERNC void Mat_diagBRC(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Enforce the Galerkin conditions algebraically.
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  None
 * @param   thee    Pointer to the sparse matrix
 * @param   rmat    R matrix which is stored column-wise (ROW-format)
 * @param   amat    A matrix which is stored in one of three
 *           forms, namely, either row-wise (ROW), col-wise (COL), or by
 *           diagonal followed by upper-triangle row-wise and then by
 *           lower
 *           triangle columne-wise (DRC).
 * @param   pmat    P matrix which is stored row-wise (ROW-format)
 */
VEXTERNC void Mat_galerkin(Mat *thee, Mat *rmat, Mat *amat, Mat *pmat);

/** 
 * @ingroup Mat
 * @brief   Make a decision about whether or not a sparse direct solver   
 *          should be used in place of an iterative solver, based on the 
 *          size of the system.
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)\n
 *          This is obviously heuristic in nature; in general the cutoff  
 *          size where iterative methods start to win is smaller in 3D.    
 * @return  the decision about whether or not a sparse direct solver
 *           should be used in place of an iterative solver, based on
 *           the size of the system.
 * @param   thee  Pointer to the sparse matrix
 */
VEXTERNC int Mat_sluDirect(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Setup for a sparse LU factorization of matrix. 
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)\n
 *          Creates internal <ia,ja,a> storage which is later freed   
 *          by Mat_sluDestroy.  Also initializes the sparse direct
 *          library.
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 */
VEXTERNC void Mat_sluCreate(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Create the sparse LU factors for the system matrix.     
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 */
VEXTERNC int Mat_sluFactor(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Performs a forward/backward solve using the sparse LU factors.   
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)\n
 *          This requires that Mat_sluFactor has been previously called.   
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   key   index for choosing NOTRANS or TRANS
 * @param   f     the number of right-hand sides
 * @param   u     solution pointer
 */
VEXTERNC int Mat_sluSolve(Mat *thee, int key, double *f, double *u);

/** 
 * @ingroup Mat
 * @brief   Destroy the sparse LU factors for the system matrix.  
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)\n
 *          This frees our <ia,ja,a> storage, and also the internal  
 *          storage that was malloc'd by the sparse direct library. 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 */
VEXTERNC void Mat_sluDestroy(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Print the exact current malloc usage      
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (mat.c)
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 */
VEXTERNC void Mat_memChk(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Copy a matrix.
 * @author  Stephen Bond
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   smat  The source matrix
 */
VEXTERNC void Mat_copy(Mat *thee, Mat *smat);

/** 
 * @ingroup Mat
 * @brief   Copy a matrix from ROW to COL.
 * @author  Stephen Bond
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   smat  the source matrix
 */
VEXTERNC void Mat_copyROW2COL(Mat *thee, Mat *smat);

/** 
 * @ingroup Mat
 * @brief   Copy a matrix from COL to ROW.
 * @author  Stephen Bond
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   smat  the source matrix
 */
VEXTERNC void Mat_copyCOL2ROW(Mat *thee, Mat *smat);

/** 
 * @ingroup Mat
 * @brief   Copy a matrix from DRC to RLN.       
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   smat  the source matrix
 */
VEXTERNC void Mat_copyDRC2RLN(Mat *thee, Mat *smat);

/** 
 * @ingroup Mat
 * @brief   Copy a matrix from ROW to RLN.       
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   smat  the source matrix
 */
VEXTERNC void Mat_copyROW2RLN(Mat *thee, Mat *smat);

/** 
 * @ingroup Mat
 * @brief   Copy a matrix from COL to RLN.     
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   smat  the source matrix
 */
VEXTERNC void Mat_copyCOL2RLN(Mat *thee, Mat *smat);

/** 
 * @ingroup Mat
 * @brief   Copy a matrix from RLN to ROW     
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   smat  the source matrix
 */
VEXTERNC void Mat_copyRLN2ROW(Mat *thee, Mat *smat);

/** 
 * @ingroup Mat
 * @brief   Copy a matrix from DRC to CLN
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   smat  the source matrix
 */
VEXTERNC void Mat_copyDRC2CLN(Mat *thee, Mat *smat);

/** 
 * @ingroup Mat
 * @brief   Copy a matrix from ROW to CLN
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   smat  the source matrix
 */
VEXTERNC void Mat_copyROW2CLN(Mat *thee, Mat *smat);

/** 
 * @ingroup Mat
 * @brief   Copy a matrix from COL to CLN
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   smat  the source matrix
 */
VEXTERNC void Mat_copyCOL2CLN(Mat *thee, Mat *smat);

/** 
 * @ingroup Mat
 * @brief   Copy a matrix from CLN to COL
 * @author  Michael Holst
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   smat  the source matrix
 */
VEXTERNC void Mat_copyCLN2COL(Mat *thee, Mat *smat);

/** 
 * @ingroup Mat
 * @brief   Copy a matrix from CLN to RLN
 * @authors Michael Holst and Stephen Bond
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   smat  the source matrix
 */
VEXTERNC void Mat_copyCLN2RLN(Mat *thee, Mat *smat);

/** 
 * @ingroup Mat
 * @brief   Copy a matrix from RLN to CLN
 * @authors Michael Holst and Stephen Bond
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   smat  the source matrix
 */
VEXTERNC void Mat_copyRLN2CLN(Mat *thee, Mat *smat);

/** 
 * @ingroup Mat
 * @brief   Copy a matrix from DRC to XLN
 * @author  Stephen Bond
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   smat  the source matrix
 */
VEXTERNC void Mat_copyDRC2XLN(Mat *thee, Mat *smat);

/** 
 * @ingroup Mat
 * @brief   Copy a matrix from ROW to XLN
 * @author  Stephen Bond
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   smat  the source matrix
 */
VEXTERNC void Mat_copyROW2XLN(Mat *thee, Mat *smat);

/** 
 * @ingroup Mat
 * @brief   Copy a matrix from COL to XLN
 * @author  Stephen Bond
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   smat  the source matrix
 */
VEXTERNC void Mat_copyCOL2XLN(Mat *thee, Mat *smat);

/** 
 * @ingroup Mat
 * @brief   Copy a matrix from RLN to XLN
 * @author  Stephen Bond
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   smat  the source matrix
 */
VEXTERNC void Mat_copyRLN2XLN(Mat *thee, Mat *smat);

/** 
 * @ingroup Mat
 * @brief   Copy a matrix from CLN to XLN
 * @author  Stephen Bond
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   smat  the source matrix
 */
VEXTERNC void Mat_copyCLN2XLN(Mat *thee, Mat *smat);

/** 
 * @ingroup Mat
 * @brief   Copy a matrix from XLN to DRC
 * @author  Stephen Bond
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   smat  the source matrix
 */
VEXTERNC void Mat_copyXLN2DRC(Mat *thee, Mat *smat);

/** 
 * @ingroup Mat
 * @brief   Copy a matrix from XLN to ROW
 * @author  Stephen Bond
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   smat  the source matrix
 */
VEXTERNC void Mat_copyXLN2ROW(Mat *thee, Mat *smat);

/** 
 * @ingroup Mat
 * @brief   Copy a matrix from XLN to COL
 * @author  Stephen Bond
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   smat  the source matrix
 */
VEXTERNC void Mat_copyXLN2COL(Mat *thee, Mat *smat);

/** 
 * @ingroup Mat
 * @brief   Copy a matrix from XLN to RLN
 * @author  Stephen Bond
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   smat  the source matrix
 */
VEXTERNC void Mat_copyXLN2RLN(Mat *thee, Mat *smat);

/** 
 * @ingroup Mat
 * @brief   Copy a matrix from XLN to CLN
 * @author  Stephen Bond
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   smat  the source matrix
 */
VEXTERNC void Mat_copyXLN2CLN(Mat *thee, Mat *smat);

/** 
 * @ingroup Mat
 * @brief   Remove the boundary rows or columns from a matrix.
 * @author  Stephen Bond
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   key   index for removing the boundary rows or columns from a matrix
 */
VEXTERNC void Mat_squeezeBRC(Mat *thee, int key);

/** 
 * @ingroup Mat
 * @brief   Raw copy of the nonzeros of X into Y.
 * @author  Stephen Bond
 * @note    Class Mat: Non-inlineable methods (matcopy.c) 
 * @return  None
 * @param   Y  the object matrix
 * @param   X  the source matrix
 */
VEXTERNC void Mat_copy2(Mat *Y, Mat *X);

/** 
 * @ingroup Mat
 * @brief   scalar times a Mat plus a Mat:  Y += val*X.   
 * @author  Stephen Bond
 * @note    Class Mat: Non-inlineable methods (matcopy.c) \n
 *          The function of this routine can be controlled using "key"
 * @return  None
 * @param   Y    the output matrix
 * @param   X    the source matrix
 * @param   val  the coeficient for scaling matrix X
 * @param   key  0 ==> X and Y have the EXACT SAME nonzero structure (nZ).
 *                     Very fast with no checking or temporary matrices\n
 *               1 ==> The nZ of Y is a SUBSET of the nZ of X.    
 *                     Still fast but requires a few extra checks.\n
 *               2 ==> X and Y have arbitrary nZ structure. 
 *                      Slowest requiring creation of a temporary link list.
 */
VEXTERNC void Mat_axpy(Mat *Y, Mat *X, double val, int key);


/** 
 * @ingroup Mat
 * @brief   Initialize the nonzero structure given structure information. 
 * @authors Stephen Bond and Michael Holst
 * @note    Class Mat: Non-inlineable methods (matln.c) \n
 *          This routine actually does the storage creation for all
 *          internal Vset, link arrays, and link pointer arrays. 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   frmt  format types of the sparse matrix
 * @param   sym   symmetric types of the sparse matrix
 */
VEXTERNC void Mat_initStructureLN(Mat *thee, MATformat frmt, MATsym sym);

/** 
 * @ingroup Mat
 * @brief   Kill the nonzero structure and structure information.       
 * @authors Stephen Bond and Michael Holst
 * @note    Class Mat: Non-inlineable methods (matln.c) \n
 *          This routine does the reverse of Mat_initStructureLN 
 *          (or Mat_copyStructureLN).  It leaves only the information 
 *          about the number of blocks, number of rows, and number of  
 *          columns per block.  I.e., what is left is only what was
 *          present after the initial call to Mat_ctor.
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 */
VEXTERNC void Mat_killStructureLN(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Access the first element in the ROW or COL of an XLN.
 * @author  Stephen Bond 
 * @note    Class Mat: Non-inlineable methods (matln.c) \n
 *          In the symmetric cases, we are just returning a pointer 
 *          to the diagonal.  In the nonsymmetric case, we are returning 
 *          a pointer to the first element in the row or column.  In  
 *          the symmetric cases the columns are linked in reverse index
 *          ordering to save the storage of an additional pointer array.  
 * @return  pointer to the diagonal element
 * @param   thee  Pointer to the sparse matrix
 * @param   idx   the index of array
 * @param   key   0 ==> return a ROW pointer\n
 *                1 ==> return a COL pointer
 */
VEXTERNC LinkRC* Mat_accessXLN(Mat *thee, int idx, int key);

/** 
 * @ingroup Mat
 * @brief   Set or add a value to a doubly linked matrix entry array.    
 * @author  Stephen Bond 
 * @note    Class Mat: Non-inlineable methods (matln.c) \n
 *          This is a doubly linked variant of mContrib.
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   key   0 ==> Set the value\n
 *                1 ==> Add the value
 * @param   i     the location when traversing rowwise
 * @param   j     the new inserted location when traversing rowwise
 * @param   val   the contribution of the new inserted position
 */
VEXTERNC void Mat_contribXLN(Mat *thee, int key, int i, int j, double val);

/** 
 * @ingroup Mat
 * @brief   Set or add a value to a NOTSYM XLN matrix.         
 * @author  Stephen Bond 
 * @note    Class Mat: Non-inlineable methods (matln.c) \n
 *          This is a doubly linked variant of mContrib.
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   key   0 ==> Set the value\n
 *                1 ==> Add the value
 * @param   i     the location when traversing rowwise
 * @param   j     the new inserted location when traversing rowwise
 * @param   val   the contribution of the new inserted position
 */
VEXTERNC void Mat_contribNSYMXLN(Mat *thee, int key, int i, int j, double val);

/** 
 * @ingroup Mat
 * @brief   Set or add a value to a STRUC_SYM XLN matrix. 
 * @author  Stephen Bond 
 * @note    Class Mat: Non-inlineable methods (matln.c) \n
 *          This is a doubly linked variant of mContrib.
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   key   0 ==> Set the value\n
 *                1 ==> Add the value
 * @param   i     the location when traversing rowwise
 * @param   j     the new inserted location when traversing rowwise
 * @param   val   the contribution of the new inserted position
 */
VEXTERNC void Mat_contribSSYMXLN(Mat *thee, int key, int i, int j, double val);

/** 
 * @ingroup Mat
 * @brief   Set or add a value to a SYM XLN Matrix.         
 * @author  Stephen Bond 
 * @note    Class Mat: Non-inlineable methods (matln.c) \n
 *          This is a doubly linked variant of mContrib.
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 * @param   key   0 ==> Set the value\n
 *                1 ==> Add the value
 * @param   i     the location when traversing rowwise
 * @param   j     the new inserted location when traversing rowwise
 * @param   val   the contribution of the new inserted position
 */
VEXTERNC void Mat_contribSYMXLN(Mat *thee, int key, int i, int j, double val);

/** 
 * @ingroup Mat
 * @brief   Print an LN format matrix as a DENSE matrix in MATLAB format.   
 * @author  Stephen Bond 
 * @note    Class Mat: Non-inlineable methods (matln.c) 
 * @return  None
 * @param   thee  Pointer to the sparse matrix
 */
VEXTERNC void Mat_printLN(Mat *thee);

/** 
 * @ingroup Mat
 * @brief   Print an LN format matrix as a SPARSE matrix in MATLAB format.     
 * @authors Michael Holst and Stephen Bond 
 * @note    Class Mat: Non-inlineable methods (matln.c) 
 * @return  None
 * @param   thee   Pointer to the sparse matrix
 * @param   fname  the output matrix name
 * @param   pflag  0 ==> write, 1 ==> append
 */
VEXTERNC void Mat_printLNSp(Mat *thee, char *fname, int pflag);

#endif /* _MAT_H_ */

