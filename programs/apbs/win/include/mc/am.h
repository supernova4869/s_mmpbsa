/**
 * @defgroup AM AM class
 * @brief    A multilevel finite element algebra object.
 */

/**
 *  @file       am.h
 *  @ingroup    AM
 *  @brief      Class AM: a multilevel finite element algebra object. 
 *  @author     Michael Holst
 *  @note       None
 *  @version    $Id: am.h,v 1.57 2010/08/12 05:19:14 fetk Exp $ 
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


#ifndef _AM_H_
#define _AM_H_

#include <mc/mc_base.h>

#include <mc/aprx.h>

/*
 * ***************************************************************************
 * Class AM: Parameters and datatypes
 * ***************************************************************************
 */

/** 
 * @ingroup AM
 * @brief   Class AM: Definition 
 * @author  Michael Holst
 */
struct sAM {

    /* ------- Objects out of our control (always exist) ------------------ */

  /** @brief Objects out of our control (always exist) the memory manager   */
    Vmem       *vmem;        
  /** @brief Objects out of our control (always exist) 
   * did i make vmem or was it inherited */
    int        iMadeVmem;   

  /** @brief Objects out of our control (always exist) 
   * ptr to the geometry manager object          */
    Aprx       *aprx;       
  /** @brief Objects out of our control (always exist) 
   * prolongation matrix maintained by aprx      */
    Bmat       *P;       

    /* ------- Objects under our control ---------------------------------- */

  /** @brief Objects under our control.
   * have the objects below been constructed     */
    int        algExist;    

  /** @brief Objects under our control.
   * block NODAL tangent matrix                  */
    Bmat       *A;   
  /** @brief Objects under our control.
   * block NODAL mass matrix                     */
    Bmat       *M;   

  /** @brief Objects under our control.
   * block NODAL load vector                     */
    Bvec       *f;  
  /** @brief Objects under our control.
   * block NODAL solution vector                 */
    Bvec       *u;     
  /** @brief Objects under our control.
   * block NODAL dirichlet vector                */
    Bvec       *ud;  
  /** @brief Objects under our control.
   * block NODAL interior dirichlet vector       */
    Bvec       *ui;   
  /** @brief Objects under our control.
   * block NODAL analytical solution             */
    Bvec       *ut;   
  /** @brief Objects under our control.
   * block NODAL residual vector                 */
    Bvec       *r; 
  /** @brief Objects under our control.
   * block NODAL work vector                     */
    Bvec       *w0;     

};

/**
 * @ingroup AM
 * @brief   Declaraction of the AM class as the AM structure 
 * @author  Michael Holst
 */
typedef struct sAM AM;

/*
 * ***************************************************************************
 * Class AM: Inlineable methods (am.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_NAM)
#else /* if defined(VINLINE_NAM) */
#endif /* if !defined(VINLINE_NAM) */

/*
 * ***************************************************************************
 * Class AM: Non-inlineable methods (am.c)
 * ***************************************************************************
 */

/**
 * @ingroup AM
 * @brief   The AM constructor.
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c) 
 * @return  Pointer to a newly allocated (empty) AM class
 * @param   vmem  Memory management object
 * @param   taprx Pointer to linear Approximation object
 */
VEXTERNC AM* AM_ctor(Vmem *vmem, Aprx *taprx);

/**
 * @ingroup AM
 * @brief   The AM destructor.
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c) 
 * @return  None
 * @param   thee Pointer to Class AM
 */
VEXTERNC void AM_dtor(AM **thee);

/**
 * @ingroup AM
 * @brief   Create the following internal Alg structures:\n
 *          A     ==> The tangent matrix (linearization operator)\n
 *          M     ==> The mass matrix\n
 *          W[]   ==> The node vectors
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c) 
 * @return  None
 * @param   thee Pointer to Class AM
 */
VEXTERNC void AM_create(AM *thee);

/**
 * @ingroup AM
 * @brief   Destroy the following internal Alg structures:\n
 *          A     ==> The tangent matrix (linearization operator)\n
 *          M     ==> The mass matrix\n
 *          W[]   ==> The node vectors
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c) 
 * @return  None
 * @param   thee Pointer to Class AM
 */
VEXTERNC void AM_destroy(AM *thee);

/**
 * @ingroup AM
 * @brief   Mark a given mesh for refinement.
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c) 
 * @return  number of marked simplices
 * @param   thee   Pointer to class AM
 * @param   key     index of different types of a posteriori error estimators
 * @param   color   chart type of marked simplices
 * @param   bkey    Bisection type
 * @param   elevel  level to determine the fraction of simplices
 */
VEXTERNC int AM_markRefine(AM *thee, int key, int color,
    int bkey, double elevel);

/**
 * @ingroup AM
 * @brief   Refine the mesh.
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c) 
 * @return  Success enumeration
 * @param   thee   Pointer to class AM
 * @param   rkey  Boolean sets the rkey type for refining the mesh.
 * Input:    If (rkey==0) Perform recursive simplex bisection until conformity
 *
 *           If (rkey==1) Perform first quadra-[octa-]-section, followed by
 *                        recursive simplex bisection until conformity
 *
 *                        IMPORTANT NOTE: In 2D, (rkey==1) WILL generate
 *                        a conforming mesh.  However, in 3D, this procedure
 *                        will in general produce nonconforming simplices.
 *                        To produce a conforming mesh in 3D would require an
 *                        implementation covering all possible face refinement
 *                        combinations (something like 169 cases).  This has
 *                        been done e.g. by Jurgen Bey in AGM, but we are
 *                        not that patient; use (rkey==0) above if you want
 *                        a conforming mesh...
 *
 *           If (rkey==2) As a test of the conformity procedure, perform
 *                        quadra-[octa-]-section until conformity, which
 *                        should produce a uniformly regularly refined mesh.
 *                        (In 2D, each triangle should be divided into four
 *                        children, and in 3D each tetrahedron should be
 *                        divided into eight children.)
 * @param   bkey  Boolean sets the bkey type for bisecting the mesh
 *           If (bkey==0) Bisection type: Longest Edge
 *           If (bkey==1) Bisection type: Newest Vertex
 *           If (bkey==2) Bisection type: Newest Pair
 * @param   pkey  Boolean sets the pkey type to prolongate a vector 
 */
VEXTERNC int AM_refine(AM *thee, int rkey, int bkey, int pkey);

/**
 * @ingroup AM
 * @brief   Un-refine the mesh.
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c) 
 * @return  Success enumeration
 * @param   thee   Pointer to class AM
 * @param   rkey  Boolean sets the rkey type for un-refining the mesh.
 * @param   pkey  Boolean sets the pkey type to prolongate a vector
 */
VEXTERNC int AM_unRefine(AM *thee, int rkey, int pkey);

/**
 * @ingroup AM
 * @brief   Deform the mesh.
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c) 
 * @return  Success enumeration for deforming the mesh.
 * @param   thee   Pointer to class AM
 */
VEXTERNC int AM_deform(AM *thee);

/**
 * @ingroup AM
 * @brief   Read in the user-specified initial mesh given in the
 *          "MCSF" or "MCEF" format, and transform into our internal
 *          datastructures.\n
 *          Do a little more than a "Aprx_read", in that we also              
 *          initialize the extrinsic and intrinsic spatial dimensions
 *          corresponding to the input mesh, and we also then build the
 *          reference elements. 
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c) \n
 *          See the documentation to Aprx_read for a description of the
 *          mesh input data file format.
 * @return  Success enumeration
 * @param   thee   Pointer to class AM
 * @param   key   input format type
 *                    key=0 ==> simplex format
 *                    key=1 ==> edge format
 *                    key=2 ==> simplex-nabor format
 * @param   sock socket for reading the external mesh data (NULL otherwise) 
 */
VEXTERNC int AM_read(AM *thee, int key, Vio *sock);

/**
 * @ingroup AM
 * @brief   Assemble the linearized problem at a given level.  
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c) 
 * @return  the assembled energy
 * @param   thee       Pointer to class AM
 * @param   evalKey    == 0 ==> Primal problem: evaluation at z=[0+ud].\n
 *                     == 1 ==> Primal problem: evaluation at z=[u+ud].\n
 *                     == 2 ==> Dual problem: evaluation at z=[0+ud].\n
 *                     == 3 ==> Dual problem: evaluation at z=[u+ud].\n
 *                     == ? ==> Same as 0.   
 * @param   energyKey  == 0 ==> DON'T assemble an energy.\n
 *                     == 1 ==> Assemble the energy J(z).\n
 *                     == 2 ==> Assemble the energy 0.\n
 *                     == ? ==> Same as 0.
 * @param   residKey   == 0 ==> DON'T assemble a residual.\n
 *                     == 1 ==> Assemble the residual F(z)(v).\n
 *                     == 2 ==> Assemble the residual 0.\n
 *                     == ? ==> Same as 0. 
 * @param   tangKey    == 0 ==> DON'T assemble a tangent matrix.\n
 *                     == 1 ==> Assemble the tangent matrix DF(z)(w,v).\n
 *                     == 2 ==> Assemble the tangent matrix 0.\n
 *                     == ? ==> Same as 0.  
 * @param   massKey    == 0 ==> DON'T assemble a mass matrix.\n
 *                     == 1 ==> Assemble the mass matrix p(w,v).\n
 *                     == 2 ==> Assemble the mass matrix 0.\n
 *                     == ? ==> Same as 0.
 * @param   bumpKey    == 0 ==> Do not assemble with bumps.\n
 *                     == 1 ==> Assemble bilinear form with bumps.\n
 *                     == 2 ==> Assemble residual form with bumps.\n
 *                     == ? ==> Same as 1.\n
 * @param   u          block NODAL solution vector
 * @param   ud         block NODAL dirichlet vector
 * @param   f          block NODAL residual vector
 * @param   ip         index for assembled energies
 * @param   rp         parameter for initially assembled PDE
 */
VEXTERNC double AM_assem(AM *thee,
    int evalKey, int energyKey, int residKey, int tangKey, int massKey,
    int bumpKey,
    Bvec *u, Bvec *ud, Bvec *f, int ip[], double rp[]);

/**
 * @ingroup AM
 * @brief   Assemble the energy functional at the current solution. 
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c) 
 * @return  the assembled energy
 * @param   thee Pointer to Class AM
 */
VEXTERNC double AM_evalJ(AM *thee);

/**
 * @ingroup AM
 * @brief   Evaluate a finite element function.  
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c) 
 * @return  None
 * @param   thee   Pointer to class AM
 * @param   number the value of variable in the environment as an integer
 * @param   block   index for the block
 * @param   numPts  number of all interpolation pts for one simplex
 * @param   pts     coordinates of the pt
 * @param   vals    solution for the interpolated pts
 * @param   marks   index for marked or not
 */
VEXTERNC void AM_evalFunc(AM *thee,
    int number, int block, int numPts, double *pts,
    double *vals, int *marks);

/**
 * @ingroup AM
 * @brief   Perform a boundary integral.      
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c) 
 * @return  None
 * @param   thee Pointer to Class AM
 */
VEXTERNC void AM_bndIntegral(AM *thee);

/**
 * @ingroup AM
 * @brief   Evaluate error in the current solution.   
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c) 
 * @return  error in solution in one of several norms.
 * @param   thee   Pointer to class AM
 * @param   pcolor  simplex chart type
 * @param   key     If (key == 0) ==> L^2 norm of the error.\n
 *                  If (key == 1) ==> L^{\\infty} norm of the error.\n
 *                  If (key == 2) ==> H^1 norm of the error. 
 */
VEXTERNC double AM_evalError(AM *thee, int pcolor, int key);

/**
 * @ingroup AM
 * @brief   Apply zero dirichlet condition at a given level.  
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c) 
 * @return  None
 * @param   thee  Pointer to class AM
 * @param   which the block vector
 */
VEXTERNC void AM_applyDiriZero(AM *thee, Bvec *which);

/**
 * @ingroup AM
 * @brief   Setup an initial guess at a given level.     
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c) 
 * @return  None
 * @param   thee  Pointer to class AM
 * @param   which the block vector
 */
VEXTERNC void AM_iniGuess(AM *thee, Bvec *which);

/**
 * @ingroup AM
 * @brief   Partition the mesh using the matching Alg level.      
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c)
 * @return  Success enumeration
 * @param   thee Pointer to Class AM
 * @param   pkey  index for different partition option
 * @param   pwht  index for weighted partitioning
 * @param   ppow  partitioning steps
 */
VEXTERNC int AM_part(AM *thee, int pkey, int pwht, int ppow);

/**
 * @ingroup AM
 * @brief   Set the partition color.
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c)
 * @return  Success enumeration
 * @param   thee   Pointer to class AM
 * @param   pcolor simplex chart type
 */
VEXTERNC int AM_partSet(AM *thee, int pcolor);

/**
 * @ingroup AM
 * @brief   Do a partition smoothing.
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c)
 * @return  counted number of partitions
 * @param   thee Pointer to Class AM
 */
VEXTERNC int AM_partSmooth(AM *thee);

/**
 * @ingroup AM
 * @brief   Print the energy.
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c)
 * @return  None
 * @param   thee Pointer to Class AM
 */
VEXTERNC void AM_printJ(AM *thee);

/**
 * @ingroup AM
 * @brief   Print the system matrix.    
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c)
 * @return  None
 * @param   thee Pointer to Class AM
 */
VEXTERNC void AM_printA(AM *thee);

/**
 * @ingroup AM
 * @brief   Print the system matrix with Dirichlet rows/cols zeroed.  
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c)
 * @return  None
 * @param   thee Pointer to Class AM
 */
VEXTERNC void AM_printAnoD(AM *thee);

/**
 * @ingroup AM
 * @brief   Print the system matrix in MATLAB sparse format. 
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c)
 * @return  None
 * @param   thee  Pointer to class AM
 * @param   fname the system matrix file name
 */
VEXTERNC void AM_printAsp(AM *thee, char *fname);

/**
 * @ingroup AM
 * @brief   Print the system matrix in MATLAB sparse format with   
 *          Dirichlet rows/cols zeroed.   
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c)
 * @return  None
 * @param   thee  Pointer to class AM
 * @param   fname the system matrix file name
 */
VEXTERNC void AM_printAspNoD(AM *thee, char *fname);

/**
 * @ingroup AM
 * @brief   Print the prolongation matrix. 
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c)
 * @return  None
 * @param   thee  Pointer to class AM
 */
VEXTERNC void AM_printP(AM *thee);

/**
 * @ingroup AM
 * @brief   Print the prolongation matrix in MATLAB sparse format. 
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c)
 * @return  None
 * @param   thee  Pointer to class AM
 * @param   fname the system matrix file name
 */
VEXTERNC void AM_printPsp(AM *thee, char *fname);

/**
 * @ingroup AM
 * @brief   Print a specified vector. 
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c)
 * @return  None
 * @param   thee Pointer to Class AM
 * @param   num  the number of block work vector
 */
VEXTERNC void AM_printV(AM *thee, int num);

/**
 * @ingroup AM
 * @brief   Print a vector in MATLAB sparse format. 
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c)
 * @return  None
 * @param   thee  Pointer to class AM
 * @param   num   the number of block work vector
 * @param   fname the output file name
 */
VEXTERNC void AM_printVsp(AM *thee, int num, char *fname);

/**
 * @ingroup AM
 * @brief   Write out a mesh in some format.  
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c)
 * @return  None
 * @param   thee      Pointer to class AM
 * @param   sock      socket for reading the external mesh data (NULL otherwise)
 * @param   defKey    defKey == 0  ==> draw mesh as it is
 *                    defKey == 1  ==> use "def??" as new vertex coords (deformation)
 *                    defKey == 2  ==> add "def??" to old vertex coords (displacement)
 *
 * @param   colKey    colKey == 0  ==> color simplices all same default color
 *                    colKey == 1  ==> color simplices based on their chart
 *                    colKey == 2  ==> color boundary simplices based on type
 *
 * @param   chartKey  chartKey <  0  ==> draw all simplices
 *                    chartKey >= 0  ==> draw only simplices with chart chartKey
 *
 * @param   gluVal    gluVal == 1. ==> draw all simplices glued together
 *               0. < gluVal <  1. ==> draw simplices with some separation
 *
 * @param   fkey      fkey   == 0  ==> draw simplices
 *                    fkey   == 1  ==> draw only simplex boundary faces
 *                    fkey   == 2  ==> draw only simplices with a boundary face
 * @param   number    the number of block work vector
 * @param   format    GV/MATH format
 */
VEXTERNC void AM_writeGEOM(AM *thee, Vio *sock,
    int defKey, int colKey, int chartKey, double gluVal, int fkey,
    int number, char *format);

/**
 * @ingroup AM
 * @brief   Write out a solution in some format.  
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c)
 * @return  None
 * @param   thee     Pointer to class AM
 * @param   sock     socket for reading the external mesh data (NULL otherwise)
 * @param   number   the number of block work vector
 * @param   format   Pointer to GV/MATH format
 */
VEXTERNC void AM_writeSOL(AM *thee, Vio *sock, int number, char *format);

/**
 * @ingroup AM
 * @brief   Print the exact current malloc usage.
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (am.c)
 * @return  None
 * @param   thee Pointer to Class AM
 */
VEXTERNC void AM_memChk(AM *thee);

/*
 * ***************************************************************************
 * Class AM: Non-inlineable methods (lsolv.c)
 * ***************************************************************************
 */

/**
 * @ingroup AM
 * @brief   Linear solver.
 * @author  Michael Holst
 * @note    Class AM: Non-inlineable methods (lsolv.c)
 * @return  None
 * @param   thee  Pointer to class AM
 * @param   prob  index for primal or dual problem
 * @param   meth  method choice (0=slu,1=mg,2=cg,3=bcg,4=pcg,5=pbcg)
 * @param   itmax number of iterations to do (the maximum allowed)
 * @param   etol  error tolerance (currently ignored)
 * @param   prec  index for different preconditioners
 * @param   gues  index for initial guess
 * @param   pjac  index for printing the system matrix
 */
VEXTERNC void AM_lSolve(AM *thee, int prob,
    int meth, int itmax, double etol, int prec, int gues, int pjac);

/**
 * @ingroup AM
 * @brief   Hierarchical linear solver.  
 * @authors Burak Aksoylu, Stephen Bond, and Michael Holst 
 * @note    Class AM: Non-inlineable methods (lsolv.c)
 * @return  None
 * @param   thee  Pointer to class AM
 * @param   prob  index for primal or dual problem
 * @param   meth  method choice (0=slu,1=mg,2=cg,3=bcg,4=pcg,5=pbcg)
 * @param   itmax number of iterations to do (the maximum allowed)
 * @param   etol  error tolerance (currently ignored)
 * @param   prec  index for different preconditioners
 * @param   gues  index for initial guess
 * @param   pjac  index for printing the system matrix
 */
VEXTERNC void AM_hlSolve(AM *thee, int prob,
    int meth, int itmax, double etol, int prec, int gues, int pjac);
#if 0

/**
 * @ingroup AM
 * @brief   PBCG method (pure bi-orthogonal extension of PCG). 
 * @authors Stephen Bond and Michael Holst
 * @note    Class AM: Non-inlineable methods (lsolv.c)\n
 * @verbatim
 * Notes:    This is a pure Petrov-Galerkin formulation of the
 *           ODIR implementation of CG, and reduces to (exactly)
 *           the ODIR implementation of CG in the case of a
 *           symmetric (possibly indefinite) system matrix A and a
 *           symmetric (possibly indefinite) preconditioner B.
 *           (However, the operation count of BCG is rougly double
 *           that of CG.)
 *
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
 *               delta_i = <B A p^{i-1}, A' q^i> / <A p^{i-1}, q^{i-1}>
 *               v       = B A  p^i - gamma_i p^i - sigma_i p^{i-1}
 *               w       = B'A' q^i - gamma_i q^i - delta_i q^{i-1}
 *               p^{i+1} = v / ||v||
 *               q^{i+1} = w / ||w||
 *           end for
 * @endverbatim
 * @return  None
 * @param   thee  Pointer to class AM
 * @param   lev   parameter to be determined
 * @param   key   coeficient for scaling the block vector
 * @param   flag  The "flag" variable determines which "mode" we run in:\n
 *                flag==0 --> Normal: check itmax, error tolerance; 
 *                            normal i/o\n
 *                flag==1 --> Silent: check only itmax; no i/o\n
 *                flag==2 --> Subcycle: check itmax, error tolerance; 
 *                            subcycle i/o\n
 *                flag==3 --> Subcycle: check itmax, error tolerance; no i/o
 * @param   itmax number of iterations to do (the maximum allowed)
 * @param   etol  error tolerance (currently ignored)
 * @param   meth  method choice (0=slu,1=mg,2=cg,3=bcg,4=pcg,5=pbcg)
 * @param   uu    index of the work block vector
 * @param   ff    index of the work block vector
 * @param   rr    index of the work block vector
 * @param   dd    index of the work block vector
 * @param   ut    index to be determined
 */
VEXTERNC void AM_hPbcg(AM *thee,
    int lev, int key, int flag, int itmax, double etol, int meth,
    int uu, int ff, int rr, int dd, int ut);
#endif

/*
 * ***************************************************************************
 * Class AM: Non-inlineable methods (nsolv.c)
 * ***************************************************************************
 */

/**
 * @ingroup AM
 * @brief   
 * @authors 
 * @note    Class AM: Non-inlineable methods (nsolv.c)
 * @return  None
 * @param   thee   Pointer to class AM
 * @param   meth   method choice (0=slu,1=mg,2=cg,3=bcg,4=pcg,5=pbcg)
 * @param   itmax  number of iterations to do (the maximum allowed)
 * @param   etol   error tolerance (currently ignored)
 * @param   lmeth  index for specified linear solver
 * @param   litmax iteration max for linear solver
 * @param   letol  error tolerance for linear solver
 * @param   lprec  preconditioner for linear solver
 * @param   gues   index for initial guess
 * @param   pjac   index for printing the system matrix
 */
VEXTERNC void AM_nSolve(AM *thee,
    int meth, int itmax, double etol,
    int lmeth, int litmax, double letol, int lprec, int gues, int pjac);

/**
 * @ingroup AM
 * @brief   Damped-Inexact Newton iteration.   
 * @authors Michael Holst
 * @note    Class AM: Non-inlineable methods (nsolv.c)
 * @return  None
 * @param   thee     Pointer to class AM
 * @param   itmax    number of iterations to do (the maximum allowed)
 * @param   etol     error tolerance (currently ignored)
 * @param   lmeth    index for specified linear solver
 * @param   litmax   iteration max for linear solver
 * @param   letol    error tolerance for linear solver
 * @param   lprec    preconditioner for linear solver
 * @param   pjac     index for printing the system matrix
 * @param   loadParm the incoming load parameter
 */
VEXTERNC void AM_newton(AM *thee,
    int itmax, double etol,
    int lmeth, int litmax, double letol, int lprec, int pjac, double loadParm);

/**
 * @ingroup AM
 * @brief   Homotopy and the related "incremental loading" iterations.
 * @authors Michael Holst
 * @note    Class AM: Non-inlineable methods (nsolv.c)
 * @return  None
 * @param   thee     Pointer to class AM
 * @param   itmax    number of iterations to do (the maximum allowed)
 * @param   etol     error tolerance (currently ignored)
 * @param   lmeth    index for specified linear solver
 * @param   litmax   iteration max for linear solver
 * @param   letol    error tolerance for linear solver
 * @param   lprec    preconditioner for linear solver
 * @param   pjac     index for printing the system matrix
 */
VPUBLIC void AM_homotopy(AM *thee,
    int itmax, double etol,
    int lmeth, int litmax, double letol, int lprec, int pjac);

#endif /* _AM_H_ */

