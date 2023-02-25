/**
 *  @file       vpred.h
 *  @brief      Header file for the Geometric Predicates.
 *  @version    $Id: vpred.h,v 1.4 2010/08/12 05:40:37 fetk Exp $ 
 *  @author     Michael Holst
 *
 *  @attention
 *  @verbatim
 *
 * MALOC = < Minimal Abstraction Layer for Object-oriented C >
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

#ifndef _VPRED_H_
#define _VPRED_H_

#include <maloc/maloc_base.h>

/* random() prototype seems to be missing in <stdlib.h> */
/*
 * if !defined(VOSF1)
 *     extern long int random(void);
 * endif
 */

/* On some machines, the exact arithmetic routines might be defeated by the  */
/*   use of internal extended precision floating-point registers.  Sometimes */
/*   this problem can be fixed by defining certain values to be volatile,    */
/*   thus forcing them to be stored to memory and rounded off.  This isn't   */
/*   a great solution, though, as it slows the arithmetic down.              */
/*                                                                           */
/* To try this out, write "#define INEXACT volatile" below.  Normally,       */
/*   however, INEXACT should be defined to be nothing.  ("#define INEXACT".) */

/** @brief  Parameters and constants "INEXACT" */
#define INEXACT /* Nothing */
/* #define INEXACT volatile */

/** @brief  float or double */
#define REAL double

/** @brief  Print the bit representation of a double */
#define REALPRINT doubleprint

/** @brief  Generate a double with random 53-bit significand and a
 random exponent in [0, 511]. */
#define REALRAND doublerand

/** @brief Generate a double with random 53-bit significand 
 and a random exponent in [0, 7]. */
#define NARROWRAND narrowdoublerand

/** @brief Generate a double with random 53-bit significand. */
#define UNIFORMRAND uniformdoublerand

/** 
 * @brief Initialize the variables used for exact arithmetic. 
 * @note  `epsilon' is the largest power of two such that 1.0 + epsilon = 1.0
 * in floating-point arithmetic.  `epsilon' bounds the relative roundoff
 * error.  It is used for floating-point error analysis.
 * `splitter' is used to split floating-point numbers into two half-
 * length significands for exact multiplication.
 * I imagine that a highly optimizing compiler might be too smart for
 * its own good, and somehow cause this routine to fail, if it pretends
 * that floating-point arithmetic is too much like real arithmetic.
 * Don't change this routine unless you fully understand it.
*/
void Vpred_exactinit(void);

/** 
 * @brief Adaptive exact 2D orientation test.  Robust.
 * @return a positive value if the points pa, pb, and pc occur in
 *         counterclockwise order; a negative value if they occur
 *         in clockwise order; and zero if they are collinear. The
 *         result is also a rough approximation of twice the signed
 *         area of the triangle defined by the three points.
 * @param pa  Pointer to a real parameter
 * @param pb  Pointer to a real parameter
 * @param pc  Pointer to a real parameter
 */      
REAL Vpred_orient2d(REAL *pa, REAL *pb, REAL *pc);

/** 
 * @brief Approximate 2D orientation test.  Nonrobust.
 * @return a positive value if the points pa, pb, and pc occur in
 *         counterclockwise order; a negative value if they occur
 *         in clockwise order; and zero if they are collinear. The
 *         result is also a rough approximation of twice the signed
 *         area of the triangle defined by the three points.
 * @param pa  Pointer to a real parameter
 * @param pb  Pointer to a real parameter
 * @param pc  Pointer to a real parameter
 */    
REAL Vpred_orient2dfast(REAL *pa, REAL *pb, REAL *pc);

/** 
 * @brief Exact 2D orientation test.  Robust.
 * @return a positive value if the points pa, pb, and pc occur in
 *         counterclockwise order; a negative value if they occur
 *         in clockwise order; and zero if they are collinear. The
 *         result is also a rough approximation of twice the signed
 *         area of the triangle defined by the three points.
 * @param pa  Pointer to a real parameter
 * @param pb  Pointer to a real parameter
 * @param pc  Pointer to a real parameter
 */
REAL Vpred_orient2dexact(REAL *pa, REAL *pb, REAL *pc);

/** 
 * @brief Adaptive exact 3D orientation test.  Robust.
 * @return a positive value if the point pd lies below the plane passing
 *         through pa, pb, and pc; "below" is defined so that pa, pb, and pc
 *         appear in counterclockwise order when viewed from above the plane.
 *         Returns a negative value if pd lies above the plane.  Returns zero if
 *         the points are coplanar.  The result is also a rough approximation of
 *         six times the signed volume of the tetrahedron defined by the four 
 *         points.
 * @param pa  Pointer to a real parameter
 * @param pb  Pointer to a real parameter
 * @param pc  Pointer to a real parameter
 * @param pd  Pointer to a real parameter
 */
REAL Vpred_orient3d(REAL *pa, REAL *pb, REAL *pc, REAL *pd);

/** 
 * @brief  Approximate 3D orientation test.  Nonrobust.
 * @return a positive value if the point pd lies below the plane passing
 *         through pa, pb, and pc; "below" is defined so that pa, pb, and pc
 *         appear in counterclockwise order when viewed from above the plane.
 *         Returns a negative value if pd lies above the plane.  Returns zero if
 *         the points are coplanar.  The result is also a rough approximation of
 *         six times the signed volume of the tetrahedron defined by the four 
 *         points.
 * @param pa  Pointer to a real parameter
 * @param pb  Pointer to a real parameter
 * @param pc  Pointer to a real parameter
 * @param pd  Pointer to a real parameter
 */
REAL Vpred_orient3dfast(REAL *pa, REAL *pb, REAL *pc, REAL *pd);

/** 
 * @brief  Exact 3D orientation test.  Robust. 
 * @return a positive value if the point pd lies below the plane passing
 *         through pa, pb, and pc; "below" is defined so that pa, pb, and pc
 *         appear in counterclockwise order when viewed from above the plane.
 *         Returns a negative value if pd lies above the plane.  Returns zero if
 *         the points are coplanar.  The result is also a rough approximation of
 *         six times the signed volume of the tetrahedron defined by the four 
 *         points.
 * @param pa  Pointer to a real parameter
 * @param pb  Pointer to a real parameter
 * @param pc  Pointer to a real parameter
 * @param pd  Pointer to a real parameter
 */
REAL Vpred_orient3dexact(REAL *pa, REAL *pb, REAL *pc, REAL *pd);

/**
 * @brief  Adaptive exact 2D incircle test.  Robust.
 * @return a positive value if the point pd lies inside the
 *         circle passing through pa, pb, and pc; a negative value if
 *         it lies outside; and zero if the four points are cocircular.
 *         The points pa, pb, and pc must be in counterclockwise
 *         order, or the sign of the result will be reversed.
 * @param pa  Pointer to a real parameter
 * @param pb  Pointer to a real parameter
 * @param pc  Pointer to a real parameter
 * @param pd  Pointer to a real parameter 
 */
REAL Vpred_incircle(REAL *pa, REAL *pb, REAL *pc, REAL *pd);

/**
 * @brief  Approximate 2D incircle test.  Nonrobust. 
 * @return a positive value if the point pd lies inside the
 *         circle passing through pa, pb, and pc; a negative value if
 *         it lies outside; and zero if the four points are cocircular.
 *         The points pa, pb, and pc must be in counterclockwise
 *         order, or the sign of the result will be reversed.
 * @param pa  Pointer to a real parameter
 * @param pb  Pointer to a real parameter
 * @param pc  Pointer to a real parameter
 * @param pd  Pointer to a real parameter 
 */
REAL Vpred_incirclefast(REAL *pa, REAL *pb, REAL *pc, REAL *pd);

/**
 * @brief  Exact 2D incircle test.  Robust.
 * @return a positive value if the point pd lies inside the
 *         circle passing through pa, pb, and pc; a negative value if
 *         it lies outside; and zero if the four points are cocircular.
 *         The points pa, pb, and pc must be in counterclockwise
 *         order, or the sign of the result will be reversed.
 * @param pa  Pointer to a real parameter
 * @param pb  Pointer to a real parameter
 * @param pc  Pointer to a real parameter
 * @param pd  Pointer to a real parameter 
 */
REAL Vpred_incircleexact(REAL *pa, REAL *pb, REAL *pc, REAL *pd);

/**
 * @brief  Adaptive exact 3D insphere test.  Robust.
 * @return a positive value if the point pe lies inside the sphere passing through
 *         pa, pb, pc, and pd; a negative value if it lies outside; and zero if the
 *         five points are cospherical.  The points pa, pb, pc, and pd must be
 *         ordered so that they have a positive orientation (as defined by
 *         orient3d()), or the sign of the result will be reversed.
 * @param pa  Pointer to a real parameter
 * @param pb  Pointer to a real parameter
 * @param pc  Pointer to a real parameter
 * @param pd  Pointer to a real parameter
 * @param pe  Pointer to a real parameter 
 */
REAL Vpred_insphere(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe);

/**
 * @brief  Approximate 3D insphere test.  Nonrobust.
 * @return a positive value if the point pe lies inside the sphere passing through
 *         pa, pb, pc, and pd; a negative value if it lies outside; and zero if the
 *         five points are cospherical.  The points pa, pb, pc, and pd must be
 *         ordered so that they have a positive orientation (as defined by
 *         orient3d()), or the sign of the result will be reversed.
 * @param pa  Pointer to a real parameter
 * @param pb  Pointer to a real parameter
 * @param pc  Pointer to a real parameter
 * @param pd  Pointer to a real parameter
 * @param pe  Pointer to a real parameter 
 */
REAL Vpred_inspherefast(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe);

/**
 * @brief  Exact 3D insphere test.  Robust.
 * @return a positive value if the point pe lies inside the sphere passing through
 *         pa, pb, pc, and pd; a negative value if it lies outside; and zero if the
 *         five points are cospherical.  The points pa, pb, pc, and pd must be
 *         ordered so that they have a positive orientation (as defined by
 *         orient3d()), or the sign of the result will be reversed.
 * @param pa  Pointer to a real parameter
 * @param pb  Pointer to a real parameter
 * @param pc  Pointer to a real parameter
 * @param pd  Pointer to a real parameter
 * @param pe  Pointer to a real parameter 
 */
REAL Vpred_insphereexact(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe);

#endif /* _VPRED_H_ */

