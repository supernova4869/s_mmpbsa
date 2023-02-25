/**
 *  @defgroup Re Re class
 *  @brief    A reference simplex element object.
 */

/**
 *  @file       re.h
 *  @ingroup    Re
 *  @brief      Class Re: a reference simplex element object.
 *  @author     Michael Holst
 *  @note       None
 *  @version    $Id: re.h,v 1.21 2010/08/12 05:18:23 fetk Exp $ 
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

#ifndef _RE_H_
#define _RE_H_

#include <mc/mc_base.h>

#include <mc/bam.h>

/**
 * @ingroup Re
 * @brief   Class Re: Parameters and datatypes. Class Re: Definition 
 * @author  Michael Holst
 */
struct sRe {

  /** @brief spatial dimension of the simplex */
    int dim;                  

  /** @brief number of i-simplices in element (v/e/f/s) */
    int dimIS[4];             
  /** @brief  degrees of freedom per i-simplex (v/e/f/s) */
    int numDF[4];             

  /** @brief order = requested quadrature order <= VMAXO */
    int qorder;              
  /** @brief numP  =  number of polys in volume element 
       =  numV * numVDF +  numE * numEDF +  numF * numFDF +  numS * numSDF
       <= VMAXP
   */
    int numP;  
  /** @brief numPS =  number of polys in surface elem
      =  numV * numVDF +  numE * numEDF +  numF * numFDF <= VMXP 
  */
    int numPS;   
  /** @brief i-simplex offsets in global numbering */
    int offIS[4];     

  /** @brief coefficients of the polynomials */
    double c[VMAXP][VMAXP];  
  /** @brief coefficients of the x-partial of polys */
    double cx[VMAXP][VMAXP]; 
  /** @brief coefficients of the y-partial of polys */
    double cy[VMAXP][VMAXP];  
  /** @brief coefficients of the y-partial of polys */
    double cz[VMAXP][VMAXP];  

  /** @brief number of volume quad pts */
    int q;  
  /** @brief number of surface quad pts */
    int qs; 

  /** @brief number of volume quad pts -- high accuracy */
    int qhi;  
  /** @brief number of surface quad pts -- high accuracy */     
    int qshi; 

  /** @brief volume quadrature point data */
    quadInfo v[VMAXQ];    
  /** @brief surface quadrature point data */    
    quadInfo s[VMAXQ][4]; 
  /** @brief permutations of surface quadrature points  
      describes how s. quad. pts transform with 
      permutations of face vertices (0,1,2)      
  */
    int spmt[6][VMAXQ]; 

  /** @brief volume quad pt data -- high accuracy */
    quadInfo vhi[VMAXQ];   
  /** @brief surface quad pt data -- high accuracy */
    quadInfo shi[VMAXQ][4];
  /** @brief permutations of hi-acc surface quad.points */ 
    int spmthi[6][VMAXQ];  

  /** @brief This is a function-pointer for a user-defined function */
    int (*simplexBasisInit)(int key, int dim, int comp,
        int *ndof, int dof[]);

  /** @brief This is a function-pointer for a user-defined function */
    void (*simplexBasisForm)(int key, int dim, int comp,
        int pdkey, double xq[], double basis[]);

};

/**
 * @ingroup Re
 * @brief   Declaration of the Re class as the Re structure 
 * @author  Michael Holst
 * @return None
 */
typedef struct sRe Re;

/*
 * ***************************************************************************
 * Class Re: Inlineable methods (re.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_APRX)
#else /* if defined(VINLINE_APRX) */
#endif /* if !defined(VINLINE_APRX) */

/**
 * @brief   Class Re: Non-inlineable methods (Re.c) 
 * @ingroup Re
 * @author  Michael Holst
 */
VEXTERNC Re* Re_ctor(int key, int dim,
    int (*simplexBasisInit)(int key, int dim, int comp,
        int *ndof, int dof[]),
    void (*simplexBasisForm)(int key, int dim, int comp,
        int pdkey, double xq[], double basis[]),int qorder);

/**
 * @ingroup Re
 * @brief   Class Re: Non-inlineable methods (Re.c)
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  None
 * @param   thee Pointer to a Re allocated memory location
 */
VEXTERNC void Re_dtor(Re **thee);

/**
 * @ingroup Re
 * @brief   the dimension of Re class
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  the dimension of Re class
 * @param   thee Pointer to a Re allocated memory location
 */
VEXTERNC int Re_dim(Re *thee);

/**
 * @ingroup Re
 * @brief   Vertex number of i-simplices in element
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  Vertex number of i-simplices in element
 * @param   thee Pointer to a Re allocated memory location
 */
VEXTERNC int Re_dimV(Re *thee);

/**
 * @ingroup Re
 * @brief   Edge number of i-simplices in element
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  Edge number of i-simplices in element
 * @param   thee Pointer to a Re allocated memory location
 */
VEXTERNC int Re_dimE(Re *thee);

/**
 * @ingroup Re
 * @brief   Face number of i-simplices in element
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  Face number of i-simplices in element
 * @param   thee Pointer to a Re allocated memory location
 */
VEXTERNC int Re_dimF(Re *thee);

/**
 * @ingroup Re
 * @brief   Simplex number of i-simplices in element
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  Simplex number of i-simplices in element
 * @param   thee Pointer to a Re allocated memory location
 */
VEXTERNC int Re_dimS(Re *thee);

/**
 * @ingroup Re
 * @brief   Vertex degrees of freedom per i-simplex
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  Vertex degrees of freedom per i-simplex
 * @param   thee Pointer to a Re allocated memory location
 */
VEXTERNC int Re_numVDF(Re *thee);

/**
 * @ingroup Re
 * @brief   Edge degrees of freedom per i-simplex
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  Edge degrees of freedom per i-simplex
 * @param   thee Pointer to a Re allocated memory location
 */
VEXTERNC int Re_numEDF(Re *thee);

/**
 * @ingroup Re
 * @brief   Face degrees of freedom per i-simplex 
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  Face degrees of freedom per i-simplex
 * @param   thee Pointer to a Re allocated memory location
 */
VEXTERNC int Re_numFDF(Re *thee);

/**
 * @ingroup Re
 * @brief   Simplex degrees of freedom per i-simplex
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  Simplex degrees of freedom per i-simplex
 * @param   thee Pointer to a Re allocated memory location
 */
VEXTERNC int Re_numSDF(Re *thee);

/**
 * @ingroup Re
 * @brief   Requested quadrature order
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  Requested quadrature order
 * @param   thee Pointer to a Re allocated memory location
 */
VEXTERNC int Re_qorder(Re *thee);

/**
 * @ingroup Re
 * @brief   number of polys in volume element
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  number of polys in volume element
 *          numP  =  numV * numVDF
 *                +  numE * numEDF                 
 *                +  numF * numFDF                    
 *                +  numS * numSDF                     
 *                <= VMAXP                            
 * @param   thee Pointer to a Re allocated memory location
 * @param   f    index for the face types
 */
VEXTERNC int Re_numP(Re *thee, int f);

/**
 * @ingroup Re
 * @brief   number of volume/surface quad pts
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  number of volume/surface quad pts
 * @param   thee Pointer to a Re allocated memory location
 * @param   f    index for the face types
 */
VEXTERNC int Re_numQ(Re *thee, int f);

/**
 * @ingroup Re
 * @brief   master element integration weight of volume/surface quadrature point data
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  master element integration weight of volume/surface quadrature point data
 * @param   thee Pointer to a Re allocated memory location
 * @param   m    index of the quadrature points
 * @param   f    index for the face types
 */
VEXTERNC double Re_w(Re *thee, int m, int f);

/**
 * @ingroup Re
 * @brief   master element integration point coordinates of volume/surface quadrature
 *          point data
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  master element integration point coordinates of volume/surface quadrature
 *          point data
 * @param   thee Pointer to a Re allocated memory location
 * @param   m    index of the quadrature points
 * @param   i    index of polys/coefs
 * @param   f    index for the face types
 */
VEXTERNC double Re_x(Re *thee, int m, int i, int f);

/**
 * @ingroup Re
 * @brief   values of local basis funcs at integ pt 
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  values of local basis funcs at integ pt
 * @param   thee Pointer to a Re allocated memory location
 * @param   m    index of the quadrature points
 * @param   i    index of polys/coefs
 * @param   f    index for the face types
 */
VEXTERNC double Re_phi(Re *thee, int m, int i, int f);

/**
 * @ingroup Re
 * @brief   x-derivs of local basis funcs at integ pt
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  x-derivs of local basis funcs at integ pt
 * @param   thee Pointer to a Re allocated memory location
 * @param   m    index of the quadrature points
 * @param   i    index of polys/coefs
 * @param   f    index for the face types
 */
VEXTERNC double Re_phix(Re *thee, int m, int i, int f);

/**
 * @ingroup Re
 * @brief   y-derivs of local basis funcs at integ pt
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  y-derivs of local basis funcs at integ pt
 * @param   thee Pointer to a Re allocated memory location
 * @param   m    index of the quadrature points
 * @param   i    index of polys/coefs
 * @param   f    index for the face types
 */
VEXTERNC double Re_phiy(Re *thee, int m, int i, int f);

/**
 * @ingroup Re
 * @brief   z-derivs of local basis funcs at integ pt
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  z-derivs of local basis funcs at integ pt
 * @param   thee Pointer to a Re allocated memory location
 * @param   m    index of the quadrature points
 * @param   i    index of the xyz coordinates
 * @param   f    index for the face types
 */
VEXTERNC double Re_phiz(Re *thee, int m, int i, int f);

/**
 * @ingroup Re
 * @brief   Derivs of local basis funcs at integ pt
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  Derivs of local basis funcs at integ pt
 * @param   thee Pointer to a Re allocated memory location
 * @param   m    index of the quadrature points
 * @param   i    index of polys/coefs
 * @param   j    index of the xyz coordinates
 * @param   f    index for the face types
 */
VEXTERNC double Re_phix2(Re *thee, int m, int i, int j, int f);

/**
 * @ingroup Re
 * @brief   hi-acc number of volume/surface quad pts
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  hi-acc number of volume/surface quad pts
 * @param   thee Pointer to a Re allocated memory location
 * @param   f    index for the face types
 */
VEXTERNC int Re_numQ_hi(Re *thee, int f);

/**
 * @ingroup Re
 * @brief   high defination master element integration weight of volume/surface 
 *          quadrature point data
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  high defination master element integration weight of volume/surface 
 *          quadrature point data
 * @param   thee Pointer to a Re allocated memory location
 * @param   m    index of the quadrature points
 * @param   f    index for the face types
 */
VEXTERNC double Re_w_hi(Re *thee, int m, int f);

/**
 * @ingroup Re
 * @brief   hi-acc master element integration point coordinates of volume/surface quadrature
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  hi-acc master element integration point coordinates of volume/surface quadrature
 * @param   thee Pointer to a Re allocated memory location
 * @param   m    index of the quadrature points
 * @param   i    index of polys/coefs
 * @param   f    index for the face types
 */
VEXTERNC double Re_x_hi(Re *thee, int m, int i, int f);

/**
 * @ingroup Re
 * @brief   
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  hi-acc values of local basis funcs at integ pt
 * @param   thee Pointer to a Re allocated memory location
 * @param   m    index of the quadrature points
 * @param   i    index of polys/coefs
 * @param   f    index for the face types
 */
VEXTERNC double Re_phi_hi(Re *thee, int m, int i, int f);

/**
 * @ingroup Re
 * @brief   hi-acc x-derivs of local basis funcs at integ pt
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  hi-acc x-derivs of local basis funcs at integ pt
 * @param   thee Pointer to a Re allocated memory location
 * @param   m    index of the quadrature points
 * @param   i    index of polys/coefs
 * @param   f    index for the face types
 */
VEXTERNC double Re_phix_hi(Re *thee, int m, int i, int f);

/**
 * @ingroup Re
 * @brief   hi-acc y-derivs of local basis funcs at integ pt
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  hi-acc y-derivs of local basis funcs at integ pt
 * @param   thee Pointer to a Re allocated memory location
 * @param   m    index of the quadrature points
 * @param   i    index of polys/coefs
 * @param   f    index for the face types
 */
VEXTERNC double Re_phiy_hi(Re *thee, int m, int i, int f);

/**
 * @ingroup Re
 * @brief   hi-acc z-derivs of local basis funcs at integ pt
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  hi-acc z-derivs of local basis funcs at integ pt
 * @param   thee Pointer to a Re allocated memory location
 * @param   m    index of the quadrature points
 * @param   i    index of polys/coefs
 * @param   f    index for the face types
 */
VEXTERNC double Re_phiz_hi(Re *thee, int m, int i, int f);

/**
 * @ingroup Re
 * @brief   hi-acc derivs of local basis funcs at integ pt
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  hi-acc derivs of local basis funcs at integ pt
 * @param   thee Pointer to a Re allocated memory location
 * @param   m    index of the quadrature points
 * @param   i    index of polys/coefs
 * @param   j    index of the xyz coordinates
 * @param   f    index for the face types
 */
VEXTERNC double Re_phix2_hi(Re *thee, int m, int i, int j, int f);

/**
 * @ingroup Re
 * @brief   permutations of surface quadrature points
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  permutations of surface quadrature points
 * @param   thee     Pointer to a Re allocated memory location
 * @param   vxn_pmt  permutations of (0,1,2) to label rows in Re:
 * @param   m        index of the quadrature points
 */
VEXTERNC int Re_sqPmt(Re *thee, int vxn_pmt, int m);

/**
 * @ingroup Re
 * @brief   permutations of hi-acc surface quad.points
 * @author  Michael Holst
 * @note    Class Re: Non-inlineable methods (Re.c) 
 * @return  permutations of hi-acc surface quad.points
 * @param   thee     Pointer to a Re allocated memory location
 * @param   vxn_pmt  permutations of (0,1,2) to label rows in Re:
 * @param   m        index of the quadrature points
 */
VEXTERNC int Re_sqPmt_hi(Re *thee, int vxn_pmt, int m);



/**
 * RYAN -- NEED TO FILL THIS IN
 * A class for capturing a (vector) finite element space.
 */

struct sFES {

  /** @brief function domain dimension */
    int dim;

  /** @brief function range dimension */
    int dimV;

  /** @brief total DOF in an i-simplex */
    int numDF;

  /** @brief the reference elements */
    Re *re[MAXV];

};

/**
 * RYAN -- NEED TO FILL THIS IN
 */
typedef struct sFES FES;


#endif /* _RE_H_ */

