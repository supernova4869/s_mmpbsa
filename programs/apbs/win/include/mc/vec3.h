/**
 * @defgroup Vec3 Vec3 class 
 * @brief    3x1 vector object
 * @defgroup Mat3 Mat3 class 
 * @brief    3x3 dense matrix object
 */

/**
 *  @file       vec3.h
 *  @ingroup    Vec3 Mat3
 *  @brief      Classes Vec3,Mat3: 3x1 vector and 3x3 dense matrix objects.
 *  @author     Michael Holst
 *  @note       None
 *  @version    $Id: vec3.h,v 1.5 2010/08/12 05:18:36 fetk Exp $ 
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

#ifndef _VEC3_H_
#define _VEC3_H_

#include <mc/mc_base.h>

/*
 * ***************************************************************************
 * Class Vec3,Mat3: Parameters and datatypes
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * Class Vec3,Mat3: Definition
 * ***************************************************************************
 */

/**
 * @ingroup Vec3
 * @brief   Vec3 definition. 3x1 vector object
 * @author  Michael Holst
 */
typedef double Vec3[3];

/**
 * @ingroup Mat3
 * @brief   Mat3 definition. 3x3 dense matrix object
 * @author  Michael Holst
 */
typedef double Mat3[3][3];

/*
 * ***************************************************************************
 * Class Vec3,Mat3: Inlineable methods (vec3.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_BAM)
#else /* if defined(VINLINE_BAM) */
#endif /* if !defined(VINLINE_BAM) */

/**
 * @ingroup Vec3
 * @brief   1-norm of a 3-vector. 
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  1-norm of a 3-vector
 * @param   u  3-vector
 */
VEXTERNC double Vec3_nrm1(Vec3 u);

/**
 * @ingroup Vec3
 * @brief   2-norm of a 3-vector. 
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  2-norm of a 3-vector
 * @param   u  3-vector
 */
VEXTERNC double Vec3_nrm2(Vec3 u);

/**
 * @ingroup Vec3
 * @brief   oo-norm of a 3-vector. 
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  oo-norm of a 3-vector
 * @param   u  3-vector
 */
VEXTERNC double Vec3_nrm8(Vec3 u);

/**
 * @ingroup Vec3
 * @brief   1-norm of a 3-vector. 
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  1-norm of a 3-vector
 * @param   u  3-vector u
 * @param   v  3-vector v
 */
VEXTERNC double Vec3_dif1(Vec3 u, Vec3 v);

/**
 * @ingroup Vec3
 * @brief   2-norm of a 3-vector. 
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  2-norm of a 3-vector
 * @param   u  3-vector u
 * @param   v  3-vector v
 */
VEXTERNC double Vec3_dif2(Vec3 u, Vec3 v);

/**
 * @ingroup Vec3
 * @brief   oo-norm of a 3-vector. 
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  oo-norm of a 3-vector
 * @param   u  3-vector u
 * @param   v  3-vector v
 */
VEXTERNC double Vec3_dif8(Vec3 u, Vec3 v);

/**
 * @ingroup Vec3
 * @brief   dot product of two 3-vectors. 
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  dot product of two 3-vectors
 * @param   u  3-vector u
 * @param   v  3-vector v
 */
VEXTERNC double Vec3_dot(Vec3 u, Vec3 v);

/**
 * @ingroup Vec3
 * @brief   Initialize a 3-vector to be a constant. 
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  None
 * @param   u    3-vector u
 * @param   val  the constant to be initialized for a 3-vector
 */
VEXTERNC void Vec3_init(Vec3 u, double val);

/**
 * @ingroup Vec3
 * @brief   3-vector scale.
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  None
 * @param   u    3-vector u
 * @param   val  the constant to be initialized for a 3-vector
 */
VEXTERNC void Vec3_scal(Vec3 u, double val);

/**
 * @ingroup Vec3
 * @brief   3-vector copy.
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  None
 * @param   u  the destiny 3-vector
 * @param   v  the source 3-vector
 */
VEXTERNC void Vec3_copy(Vec3 u, Vec3 v);

/**
 * @ingroup Vec3
 * @brief   scalar times 3-vector plus 3-vector.   
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  None
 * @param   u    the output 3-vector
 * @param   v    the source 3-vector
 * @param   val  the coeficient for scaling 3-vector v
 */
VEXTERNC void Vec3_axpy(Vec3 u, Vec3 v, double val);

/**
 * @ingroup Vec3
 * @brief   cross-product of two 3-vectors. 
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  None
 * @param   u    the output 3-vector
 * @param   v    the source 3-vector
 * @param   w    another source 3-vector
 */
VEXTERNC void Vec3_xcry(Vec3 u, Vec3 v, double *w);

/**
 * @ingroup Vec3
 * @brief   normalize a 3-vector.  
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  None
 * @param   u      the output 3-vector
 * @param   scale  scale value of the 3-vector u
 */
VEXTERNC void Vec3_nrmlize(Vec3 u, double scale);

/**
 * @ingroup Vec3
 * @brief   normalize a 3-vector (no error check).
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  None
 * @param   u      the output 3-vector
 * @param   scale  scale value of the 3-vector u
 */
VEXTERNC void Vec3_nrmlizeNE(Vec3 u, double scale);

/**
 * @ingroup Vec3
 * @brief   print a 3-vector.
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  None
 * @param   u     the output 3-vector
 * @param   name  the output file name
 */
VEXTERNC void Vec3_print(Vec3 u, const char *name);

/**
 * @ingroup Vec3
 * @brief   multiply a 3-matrix and a 3-vector.
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  None
 * @param   u  the output 3-vector
 * @param   A  the 3-matrix
 * @param   v  the source 3-vector
 */
VEXTERNC void Vec3_mult(Vec3 u, Mat3 A, Vec3 v);

/**
 * @ingroup Vec3
 * @brief   get 3-vector column of a 3-matrix.   
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  None
 * @param   u    the output 3-vector
 * @param   A    the 3-matrix
 * @param   col  the index of column
 */
VEXTERNC void Vec3_getCol(Vec3 u, Mat3 A, int col);

/**
 * @ingroup Vec3
 * @brief   get 3-vector row of a 3-matrix.   
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  None
 * @param   u    the output 3-vector
 * @param   A    the 3-matrix
 * @param   row  the index of row
 */
VEXTERNC void Vec3_getRow(Vec3 u, Mat3 A, int row);

/**
 * @ingroup Mat3
 * @brief   1-norm of a 3-matrix.
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  1-norm of a 3-matrix
 * @param   A  the 3-matrix
 */
VEXTERNC double Mat3_nrm1(Mat3 A);

/**
 * @ingroup Mat3
 * @brief   2-norm of a 3-matrix.
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  2-norm of a 3-matrix
 * @param   A  the 3-matrix
 */
VEXTERNC double Mat3_nrm2(Mat3 A);

/**
 * @ingroup Mat3
 * @brief   oo-norm of a 3-matrix.
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  oo-norm of a 3-matrix
 * @param   A  the 3-matrix
 */
VEXTERNC double Mat3_nrm8(Mat3 A);

/**
 * @ingroup Mat3
 * @brief   1-norm of difference of two 3-matrices.     
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  1-norm of difference of two 3-matrices
 * @param   A  the 3-matrix A
 * @param   B  the 3-matrix B
 */
VEXTERNC double Mat3_dif1(Mat3 A, Mat3 B);

/**
 * @ingroup Mat3
 * @brief   2-norm of difference of two 3-matrices.     
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  2-norm of difference of two 3-matrices
 * @param   A  the 3-matrix A
 * @param   B  the 3-matrix B
 */
VEXTERNC double Mat3_dif2(Mat3 A, Mat3 B);

/**
 * @ingroup Mat3
 * @brief   oo-norm of difference of two 3-matrices.     
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  oo-norm of difference of two 3-matrices
 * @param   A  the 3-matrix A
 * @param   B  the 3-matrix B
 */
VEXTERNC double Mat3_dif8(Mat3 A, Mat3 B);

/**
 * @ingroup Mat3
 * @brief   oo-norm of lower-triangle of a 3-matrix. 
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  oo-norm of lower-triangle of a 3-matrix
 * @param   A  the 3-matrix
 */
VEXTERNC double Mat3_nrm8Low(Mat3 A);

/**
 * @ingroup Mat3
 * @brief   identity 3-matrix.    
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  None
 * @param   A  the 3-matrix
 */
VEXTERNC void Mat3_eye(Mat3 A);

/**
 * @ingroup Mat3
 * @brief   initialize a 3-matrix.     
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  None
 * @param   A    the 3-matrix
 * @param   val  initialized value of the 3-matrix component
 */
VEXTERNC void Mat3_init(Mat3 A, double val);

/**
 * @ingroup Mat3
 * @brief   normalize a 3-matrix.      
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  None
 * @param   A    the 3-matrix
 * @param   val  the coeficient of the 3-matrix components
 */
VEXTERNC void Mat3_scal(Mat3 A, double val);

/**
 * @ingroup Mat3
 * @brief   copy a 3-matrix. 
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  None
 * @param   A  the output 3-matrix
 * @param   B  the source 3-matrix
 */
VEXTERNC void Mat3_copy(Mat3 A, Mat3 B);

/**
 * @ingroup Mat3
 * @brief   Saxpy for 3-matrices.   
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  None
 * @param   A    the output 3-matrix
 * @param   B    the source 3-matrix
 * @param   val  the coeficient of the 3-matrix B components
 */
VEXTERNC void Mat3_axpy(Mat3 A, Mat3 B, double val);

/**
 * @ingroup Mat3
 * @brief   multiply 3-matrices.
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  None
 * @param   C  the output 3-matrix
 * @param   A  the source 3-matrix A
 * @param   B  the source 3-matrix B
 */
VEXTERNC void Mat3_mult(Mat3 C, Mat3 A, Mat3 B);

/**
 * @ingroup Mat3
 * @brief   put a 3-vector column in a 3-matrix. 
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  None
 * @param   A    the 3-matrix
 * @param   u    the 3-vector
 * @param   col  the index of the 3-vector column in a 3-matrix
 */
VEXTERNC void Mat3_putCol(Mat3 A, Vec3 u, int col);

/**
 * @ingroup Mat3
 * @brief   put a 3-vector row in a 3-matrix. 
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  None
 * @param   A    the 3-matrix
 * @param   u    the 3-vector
 * @param   row  the index of the 3-vector row in a 3-matrix
 */
VEXTERNC void Mat3_putRow(Mat3 A, Vec3 u, int row);

/**
 * @ingroup Mat3
 * @brief   print a 3-matrix.
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  None
 * @param   A     the 3-matrix
 * @param   name  the output name of the 3-matrix
 */
VEXTERNC void Mat3_print(Mat3 A, const char *name);

/**
 * @ingroup Mat3
 * @brief   QR iteration for 3-matrices.
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  QR iteration for 3-matrices
 * @param   V  the output 3-matrix
 * @param   D  the diagonal 3-matrix
 * @param   A  the source 3-matrix
 */
VEXTERNC double Mat3_qri(Mat3 V, Mat3 D, Mat3 A);

/**
 * @ingroup Mat3
 * @brief   QR factorization of 3-matrix (via modified Graham-Schmidt). 
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  None
 * @param   Q  the Q 3-matrix
 * @param   R  the R 3-matrix
 * @param   A  the source 3-matrix
 */
VEXTERNC void Mat3_gramSch(Mat3 Q, Mat3 R, Mat3 A);

/**
 * @ingroup Mat3
 * @brief   a single QR iteration for a 3-matrix.  
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  None
 * @param   Q  the Q 3-matrix
 * @param   R  the R 3-matrix
 * @param   A  the source 3-matrix
 */
VEXTERNC void Mat3_qr(Mat3 Q, Mat3 R, Mat3 A);

/**
 * @ingroup Mat3
 * @brief   determinant of a 3-matrix.    
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  determinant of a 3-matrix
 * @param   A  the 3-matrix
 */
VEXTERNC double Mat3_det(Mat3 A);

/**
 * @ingroup Mat3
 * @brief   inverse of a 3-matrix.
 * @author  Michael Holst
 * @note    Class Vec3,Mat3: Non-inlineable methods (vec3.c) 
 * @return  None
 * @param   A     the 3-matrix
 * @param   Ainv  inverse of a 3-matrix A
 */
VEXTERNC void Mat3_inverse(Mat3 A, Mat3 Ainv);

#endif /* _VEC3_H_ */

