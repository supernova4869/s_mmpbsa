/**
 * @defgroup PDE PDE class
 * @brief    the differential equation object.
 */

/**
 *  @file       pde.h
 *  @ingroup    PDE
 *  @brief      Class PDE: the differential equation object.
 *  @note
 *  @verbatim
 *           The differential operator information is provided by
 *           specifying the following:
 *
 *       vec               Product dim (unknowns per spatial point)
 *       sym[MAXV][MAXV];  Symmetry/nonsymmetry of bilinear form
 *
 *       initAssemble      Once-per-assembly initialization
 *       initElement       Once-per-element initialization
 *       initFace          Once-per-face initialization
 *       initPoint         Once-per-point (e.g. quadrature pt) initialization
 *
 *       Fu                The nonlinear strong form defining the PDE
 *       Ju                The nonlinear energy (PDE is the Euler condition)
 *       Fu_v              The nonlinear weak form defining the PDE
 *       DFu_wv            The bilinear linearization weak form
 *       delta             The delta function source term (if present)
 *       u_D               The dirichlet boundary condition function
 *       u_T               Optional analytical solution function (for testing)
 *
 *       bisectEdge        The rule used to spit edges of simplices
 *       mapBoundary       The rule used to recover boundary surfaces exactly
 *       markSimplex       The rule used to mark simplices for refinement
 *       oneChart          Coordinate transformations
 *
 *       simplexBasisInit  Trial and test space bases initialization
 *       simplexBasisForm  Trial and test space bases formation
 *
 *          The following two additional parameters are determined by the
 *          input manifold, and are then set appropriately to be available
 *          to the user's PDE definition routines above:
 *
 *       dim               Manifold intrinsic dim (determined by mesh input)
 *       dimII             Manifold imbedding dim (determined by mesh input)
 *
 * Notes on the forms: Fu_v/DFu_wv
 *
 *           At the single given point x, we must evaluate the functionals
 *
 *               Fu_v   = F(u)(v)    = F (x, u(x),grad u(x), v(x),grad v(x))
 *
 *               DFu_wv = DF(u)(w,v) = DF(x, u(x),grad u(x), v(x),grad v(x),
 *                                           w(x),grad w(x))
 *
 *           The functional F(u)(v):HxH-->R is assumed to be linear in 
 *           the argument v, but possibly nonlinear in the argument u.
 *
 *           The functional DF(u)(w,v):HxHxH-->R is assumed to be linear in
 *           the arguments v and w, but possibly nonlinear in the argument u.
 *
 *           Here, H is some Banach space of functions, such as the
 *           Sobolev space W^{1,p}(\Omega), 1<=p<\infty, or some product
 *           space based on these spaces.
 *
 * The forms Fu_v/DFu_wv and their use in Newton iterations:
 *
 *           The functional F(u)(v) is a weak formulation of some nonlinear 
 *           partial differential equation.  This functional represents a 
 *           nonlinear operator in an equation of the form:
 *
 *               (P1) Find u in H such that F(u)(v)=0 for all v in H.
 *
 *           As a function of its first argument only, the functional F(u)(v) 
 *           can be viewed as a mapping from the space H into the space of 
 *           linear functionals on H, or F(u)(v):H->L(H,R).  Similarly, the 
 *           Frechet derivative of F(u)(v) with respect to u represents a 
 *           mapping DF(u)(w,v):H->L(HxH,R), which is a mapping from H into
 *           the space of bilinear forms on H.  At a point z in H, such a
 *           bilinear form will be denoted:
 *          
 *               DF(z)(u,v).
 *
 *           A Newton iteration for solving (P1) takes the form:
 *
 *               (0) Given some initial guess of u in H.
 *                   (For linear problems, we will always start with u=0.)
 *               (1) Find w in H such that DF(u)(w,v)=-F(u)(v), for all v in H.
 *               (2) Update u:  u = u + w
 *               (3) If (F(u)(v)>TOL) for some v in H, goto (1)
 *               (4) Quit
 *
 *           If this is a linear problem, in that F(u)(v) is only affine in the
 *           argument u rather than truly nonlinear, then the operator has 
 *           the general form:
 *
 *               F(u)(v) = A(u,v) - (f,v),
 *
 *           where A(u,v) is bilinear.  Therefore, problem (P1) has the form:
 *
 *               (P2) Find u in H such that A(u,v)=(f,v) for all v in H.
 *
 *           The Frechet derivative of F(u)(v) with respect to u at a point z
 *           is simply the bilinear form (independent of z)
 *
 *               DF(z)(u,v) = A(u,v).
 *
 *           Now, since we will require u=0 to start a Newton iteration for a 
 *           linear problem, this will produce the linear source term as 
 *           required to solve the linear problem:
 *
 *               (Aw,v) = DF(0)(w,v) = -F(0)(v) = (f,v)
 *
 *           So, Step (1) in our Newton iteration will be exactly the
 *           solution of the linear problem (P2), which is the proper
 *           linearization of a linear (affine) problem.  Step (2) of the
 *           Newton iteration them forms u as u=0+w, and Step (3) then sees
 *           F(u)(v)=A(w,v)-(f,v)=0, which leads to completion in Step (4).
 *
 * Sub-form implementation of Fu_v/DFu_wv and computational details:
 *
 *           For the sake of efficiency, these forms may be implemented in
 *           terms of sub-forms on orthogonal spaces.  For systems of partial
 *           differential equations, the functions u,v,w above have the form:
 *
 *               u=(u_1, u_2, ... , u_m),
 *
 *           so that u then lies in the product space H=[H_1xH_2x...xH_m],
 *           with each component function u_i in the space H_i.
 *
 *           While the argument appearing nonlinearly in each form (namely u)
 *           may have all components u_i nonzero, it is assumed that the other
 *           two arguments which appear linearly always have only ONE component
 *           (v_i != 0) for a given i; thus,
 *
 *               v=(0, 0, ... , 0, v_i, 0, ... ,0), for some i,
 *               w=(0, 0, ... , 0, w_j, 0, ... ,0), for some j.
 *
 *           There may be (quite) significant reduction in computational costs
 *           if this is known apriori.  The arguments v and w will ALWAYS have
 *           this form.  Due to their linear appearance in the forms and the
 *           orthogonality of spaces H_i and H_j for (i!=j) in the setting of
 *           a product of spaces, the above forms can be viewed as:
 *
 *               F(u)(v)    = F_1(u)(v_1) + ... + F_m(u)(v_m).
 *
 *               DF(u)(w,v) = DF_11(u)(w_1,v_1) + ... + DF_1m(u)(w_1,v_m)
 *                          + DF_21(u)(w_2,v_1) + ... + DF_2m(u)(w_2,v_m)
 *                                                   .
 *                                                   .
 *                                                   .
 *                          + DF_m1(u)(w_m,v_1) + ... + DF_mm(u)(w_m,v_m)
 *
 *           Therefore, it may be more efficient to represent the forms
 *           F(u)(v) and DF(u)(w,v) as additive collections of the sub-forms:
 *
 *               F_i(u)(v_i),        i=1,...,m
 *
 *               DF_ij(u)(w_j,v_i),  i=1,...,m,  j=1,...,m
 *
 *           In fact, in a finite element discretization, somewhere we will
 *           have to split the forms into the above sub-forms.  We can attempt
 *           to exploit this known structure and save on some computation by
 *           designing the F(u)(v) and DF(u)(w,v) routines to return ARRAYS
 *           of values, GAMMA[m][1] and THETA[m][m], containing the values of
 *           the above sub-forms at the given input point.  If it is difficult
 *           to split the forms by hand, one can of course simply use the full
 *           forms F(u)(v) and DF(u)(w,v) to fill the arrays GAMMA and THETA;
 *           you just need to plug in the basis functions having the form:
 *
 *               v=(0, 0, ... , 0, v_i, 0, ... ,0),
 *               w=(0, 0, ... , 0, w_j, 0, ... ,0),
 *
 *           taking i=1,...,m, and also j=1,...,m.
 *
 *  @endverbatim
 *  @version    $Id: pde.h,v 1.23 2010/08/12 05:19:18 fetk Exp $ 
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

#ifndef _PDE_H_
#define _PDE_H_

#include <mc/mc_base.h>

/*
 * ***************************************************************************
 * Class PDE: Parameters and datatypes
 * ***************************************************************************
 */

/** 
 * @ingroup PDE
 * @brief   Class PDE: Definition 
 * @author  Michael Holst
 */

struct sPDE {

   /*
    * ************************************************************************
    * User-supplied problem information -- ALL must be supplied by the user.
    * ************************************************************************
    */

  /** @brief User-supplied - Manifold intrinsic dim (determined by mesh)   */
    int dim;         
  /** @brief User-supplied - Manifold imbedding dim (determined by mesh)   */
    int dimII;            

  /** @brief User-supplied - Product dim (unknowns per spatial point)      */
    int vec;              
  /** 
   * @brief User-supplied symmetry/nonsymmetry of bilinear form:\n
   *        sym[i][j]=0 ==> nonsym bilinear form\n
   *        sym[i][j]=1 ==> sym bilinear form\n      
   *        sym[i][j]=2 ==> sym image of block in upper triangle
   *                        (2 requires i>j)             
   */
    int sym[MAXV][MAXV];  

  /** @brief User-supplied - error estimator weights                */
    double est[MAXV];             
  /** @brief User-supplied - boundary type (re)mapping              */
    int bmap[MAXV][VMAX_BDTYPE];  

  /** @brief Once-per-assembly/element/face/point initialization */
    void (*initAssemble)(struct sPDE *thee, int ip[], double rp[]);
  /** @brief Once-per-assembly/element/face/point initialization */
    void (*initElement)(struct sPDE *thee, int elementType,
        int chart, double tvx[][3], void *data);
  /** @brief Once-per-assembly/element/face/point initialization */
    void (*initFace)(struct sPDE *thee, int faceType,
        int chart, double tnvec[]);
  /** @brief Once-per-assembly/element/face/point initialization */
    void (*initPoint)(struct sPDE *thee, int pointType,
        int chart, double txq[],
        double U[], double dU[][3]);

  /** @brief Nonlinear/bilinear strong/weak forms defining PDE and 
   *linearization */
    void (*Fu)(struct sPDE *thee, int key, double F[]);
  /** @brief Nonlinear/bilinear strong/weak forms defining PDE and 
   *linearization */
    double (*Ju)(struct sPDE *thee, int key);
  /** @brief Nonlinear/bilinear strong/weak forms defining PDE and 
   *linearization */
    double (*Fu_v)(struct sPDE *thee, int key,
        double V[], double dV[][3]);
  /** @brief Nonlinear/bilinear strong/weak forms defining PDE and 
   *linearization */
    double (*DFu_wv)(struct sPDE *thee, int key,
        double W[], double dW[][3],
        double V[], double dV[][3]);
  /** @brief Nonlinear/bilinear strong/weak forms defining PDE and 
   *linearization */
    double (*p_wv)(struct sPDE *thee,
        int key, double W[], double V[]);

  /** @brief Delta functions, Dirichlet boundary condition, analytical 
   * solution */
    void (*delta)(struct sPDE *thee, int type,
        int chart, double txq[], void *data, double F[]);
  /** @brief Delta functions, Dirichlet boundary condition, analytical 
   * solution */
    void (*u_D)(struct sPDE *thee, int type,
        int chart, double txq[], double F[]);
  /** @brief Delta functions, Dirichlet boundary condition, analytical 
   * solution */
    void (*u_T)(struct sPDE *thee, int type,
		int chart, double txq[], double F[], double dF[][3]);

  /** @brief Bisect edge rule (for refinement) and the chart unification 
   * routine */
    void (*bisectEdge)(int dim, int dimII,
        int edgeType, int chart[], double vx[][3]);
  /** @brief Bisect edge rule (for refinement) and the chart unification 
   * routine */
    void (*mapBoundary)(int dim, int dimII, int vertexType,
        int chart, double vx[3]);
  /** @brief Bisect edge rule (for refinement) and the chart unification 
   * routine */
    int (*markSimplex)(int dim, int dimII,
        int simplexType, int faceType[4], int vertexType[4],
        int chart[], double vx[][3], void *data);
  /** @brief Bisect edge rule (for refinement) and the chart unification 
   * routine */
    void (*oneChart)(int dim, int dimII,
        int objType, int chart[], double vx[][3], int dimV);

  /** @brief Trial and test space bases evaluation routine */
    int (*simplexBasisInit)(int key, int dim, int comp,
        int *ndof, int dof[]);
  /** @brief Trial and test space bases evaluation routine */
    void (*simplexBasisForm)(int key, int dim, int comp,
        int pdkey, double xq[], double basis[]);

  /** @brief To allow user to hang problem-specific variables off PDE
   * structure */ 
    void *user;

};

/**
 * @ingroup PDE
 * @brief   Declaration of the PDE class as the PDE structure
 * @author  Michael Holst
 */
typedef struct sPDE PDE;

/*
 * ***************************************************************************
 * Class PDE: Inlineable methods (pde.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_PDE)
    /**   
     * @ingroup PDE 
     * @brief   Set the extrinsic spatial dimension. 
     * @author  Michael Holst
     * @note    Class PDE: Inlineable methods (pde.c) 
     * @return  None
     * @param   thee Pointer to the PDE object
     * @param   d    the extrinsic spatial dimension
     */
    VEXTERNC void PDE_setDim(PDE *thee, int d);

    /**   
     * @ingroup PDE 
     * @brief   Set the extrinsic spatial dimension. 
     * @author  Michael Holst
     * @note    Class PDE: Inlineable methods (pde.c) 
     * @return  None
     * @param   thee Pointer to the PDE object
     * @param   d    the intrinsic spatial dimension
     */
    VEXTERNC void PDE_setDimII(PDE *thee, int d);

    /**   
     * @ingroup PDE 
     * @brief   Return the extrinsic spatial dimension. 
     * @author  Michael Holst
     * @note    Class PDE: Inlineable methods (pde.c) 
     * @return  the extrinsic spatial dimension. 
     * @param   thee Pointer to the PDE object
     */
    VEXTERNC int PDE_dim(PDE *thee);

    /**   
     * @ingroup PDE 
     * @brief   Return the extrinsic spatial dimension. 
     * @author  Michael Holst
     * @note    Class PDE: Inlineable methods (pde.c) 
     * @return  the extrinsic spatial dimension. 
     * @param   thee Pointer to the PDE object
     */
    VEXTERNC int PDE_dimII(PDE *thee);

    /**   
     * @ingroup PDE 
     * @brief   Return the PDE product space dimension.   
     * @author  Michael Holst
     * @note    Class PDE: Inlineable methods (pde.c) 
     * @return  the PDE product space dimension.   
     * @param   thee Pointer to the PDE object
     */
    VEXTERNC int PDE_vec(PDE *thee);
#else /* if defined(VINLINE_PDE) */
#   define PDE_setDim(thee,d)         ( (thee)->dim = (d) )
/** @brief Class PDE: Inlineable methods (pde.c) if defined(VINLINE_PDE) */
#   define PDE_setDimII(thee,d)       ( (thee)->dimII = (d) )
/** @brief Class PDE: Inlineable methods (pde.c) if defined(VINLINE_PDE) */
#   define PDE_dim(thee)              ( (thee)->dim )
/** @brief Class PDE: Inlineable methods (pde.c) if defined(VINLINE_PDE) */
#   define PDE_dimII(thee)            ( (thee)->dimII )
/** @brief Class PDE: Inlineable methods (pde.c) if defined(VINLINE_PDE) */
#   define PDE_vec(thee)              ( (thee)->vec )
#endif /* if !defined(VINLINE_PDE) */

/*
 * ***************************************************************************
 * Class PDE: Non-inlineable methods (pde.c)
 * ***************************************************************************
 */

/**   
 * @ingroup PDE 
 * @brief   Construct a fake differential equation object in case there  
 *          is not one provided.
 * @author  Michael Holst
 * @note    Class PDE: Non-Inlineable methods (pde.c) 
 * @return  Pointer to a new allocated PDE object
 */
VEXTERNC PDE* PDE_ctor_default(void);

/**   
 * @ingroup PDE 
 * @brief   Destroy a fake differential equation object. 
 * @author  Michael Holst
 * @note    Class PDE: Non-Inlineable methods (pde.c) 
 * @return  None
 * @param   thee Pointer to the PDE object
 */
VEXTERNC void PDE_dtor_default(PDE **thee);

/**   
 * @ingroup PDE 
 * @brief   Return the bilinear form minor symmetry (i,j).         
 * @author  Michael Holst
 * @note    Class PDE: Non-Inlineable methods (pde.c) 
 * @return  the bilinear form minor symmetry (i,j).         
 * @param   thee Pointer to the PDE object
 * @param   i    index for the symmetry test
 * @param   j    index for the symmetry test
 */
VEXTERNC int PDE_sym(PDE *thee, int i, int j);

#endif /* _PDE_H_ */

