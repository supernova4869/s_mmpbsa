/**
 * @defgroup Aprx Aprx class
 * @brief    The APpRoXimation library object. 
 */

/**
 *  @file       aprx.h
 *  @ingroup    Aprx
 *  @brief      Class Aprx: the APpRoXimation library object.
 *  @author     Michael Holst
 *  @note       None
 *  @version    $Id: aprx.h,v 1.37 2010/08/12 05:18:21 fetk Exp $ 
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

#ifndef _APRX_H_
#define _APRX_H_

#include <mc/mc_base.h>

#include <mc/gem.h>
#include <mc/bam.h>
#include <mc/whb.h>

#include <mc/re.h>
#include <mc/node.h>
#include <mc/bnode.h>

/**
 * @ingroup Aprx
 * @brief Class Aprx: Parameters and datatypes Class Aprx: Definition. 
 * @author  Michael Holst 
 */

struct sAprx {

    /** @brief the memory manager                          */
    Vmem     *vmem;         
    /** @brief did i make vmem or was it inherited         */
    int      iMadeVmem;     

    /** @brief geometry machine                            */
    Gem      *gm;       
    /** @brief container for user-provided functions       */
    PDE      *pde;      

    /** @brief ptr to the Vset object for interactions     */
    Vset     *lnkG;      

    /** @brief ptr to the master elements                  */
    Re       *re[MAXV];   
    /** @brief ptr to the master bump elements             */  
    Re       *reB[MAXV];  

    /** @brief Global Number of 0-simplices (vertices)     */
    int      numV;
    /** @brief Global Number of 1-simplices (edges)        */
    int      numE;
    /** @brief Global Number of 2-simplices (faces)        */
    int      numF;
    /** @brief Global Number of 3-simplices (simplices)    */
    int      numS;

    /** @brief Global Vertex-to-Edge map                   */
    int      *I_V2E;     
    /** @brief Global Vertex-to-Face map                   */
    int      *I_V2F;     
    /** @brief Global Vertex-to-Simplex map                */
    int      *I_V2S;     
    /** @brief Global Edge-to-Face map                     */
    int      *I_E2F;     
    /** @brief Global Edge-to-Simplex map                  */
    int      *I_E2S;     
    /** @brief Global Face-to-Simplex map                  */
    int      *I_F2S;     

    /** @brief complete node information                   */
    Bnode    *B;         

    /** @brief number of dirichlet boundary nodes          */
    int      numBV[MAXV]; 
    /** @brief dirichlet boundary node numbers             */
    int      *ibv[MAXV];  
    /** @brief dirichlet boundary node values              */
    double   *bv[MAXV];    

    /** @brief ELEMENT work vectors                        */
    Bvec     *wev;        
    /** @brief current a posteriori global error over mesh */
    double   gerror;        

    /** @brief number of blocks                            */
    int      numB;       
    /** @brief block prolongation operator TO this mesh    */
    Bmat     *P;           

    /** @brief order of approximation */
    int order;

};

/**
 * @ingroup Aprx
 * @brief   Declaration of the Aprx class as the Aprx structure
 * @author  Michael Holst
 * @return  None
 */

typedef struct sAprx Aprx;

/*
 * ***************************************************************************
 * Class Aprx: Inlineable methods (aprx.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_APRX)
#else /* if defined(VINLINE_APRX) */
#endif /* if !defined(VINLINE_APRX) */

/**
 * @ingroup Aprx
 * @brief   Construct a linear Approximation object.
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (aprx.c)
 * @return  Pointer to a newly allocated (empty) linear Approximation object
 * @param   vmem  Memory management object
 * @param   tgm   Geometry manager source
 * @param   tpde  PDE object
 */
VEXTERNC Aprx* Aprx_ctor(Vmem *vmem, Gem *tgm, PDE *tpde);

/**
 * @ingroup Aprx
 * @brief   Approximation object constructor.
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (aprx.c)
 * @return  thee  Pointer to an Aprx allocated memory location
 * @param   vmem  Memory management object
 * @param   tgm   Geometry manager source
 * @param   tpde  PDE object
 * @param   order Order of approximation to be constructed
 */
VEXTERNC Aprx* Aprx_ctor2(Vmem *vmem, Gem *tgm, PDE *tpde,int order);

/**
 * @ingroup Aprx
 * @brief   Approximation object destructor. 
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (aprx.c)
 * @param   thee Pointer to object to be destroyed
 */
VEXTERNC void Aprx_dtor(Aprx **thee);

/**
 * @ingroup Aprx
 * @brief   Return the number of blocks.
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (aprx.c)
 * @return  the number of blocks
 * @param   thee Pointer to an Aprx allocated memory location
 */
VEXTERNC int Aprx_numB(Aprx *thee);

/**
 * @ingroup Aprx
 * @brief   Return a pointer to the current prolongation matrix.  
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (aprx.c)
 * @return  a pointer to the current prolongation matrix.  
 * @param   thee Pointer to an Aprx allocated memory location
 */
VEXTERNC Bmat *Aprx_P(Aprx *thee);

/**
 * @ingroup Aprx
 * @brief   Read in the user-specified initial mesh given in the 
 *          "MCSF" or "MCEF" format, and transform into our internal 
 *          datastructures.\n
 *          Do a little more than a "Gem_read", in that we also 
 *          initialize the extrinsic and intrinsic spatial dimensions 
 *          corresponding to the input mesh, and we also then build the 
 *          reference elements.  
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (aprx.c)\n
 *          See the documentation to Gem_read for a description of the
 *          mesh input data file format.  
 * @return  Success enumeration
 * @param   thee Pointer to an Aprx allocated memory location
 * @param   key   input format type
 *                    key=0 ==> simplex format
 *                    key=1 ==> edge format
 *                    key=2 ==> simplex-nabor format
 * @param   sock socket for reading the external mesh data (NULL otherwise) 
 */
VEXTERNC int Aprx_read(Aprx *thee, int key, Vio *sock);

/**
 * @ingroup Aprx
 * @brief   Reset all of the datastructures. 
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (aprx.c)
 * @return  None
 * @param   thee Pointer to an Aprx allocated memory location
 */
VEXTERNC void Aprx_reset(Aprx *thee);

/**
 * @ingroup Aprx
 * @brief   Make node information.
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (aprx.c)
 * @return  None
 * @param   thee Pointer to an Aprx allocated memory location
 */
VEXTERNC void Aprx_initNodes(Aprx *thee);

/**
 * @ingroup Aprx
 * @brief   Determine the node types via a single fast traversal of
 *          the elements by looking at the element faces.     
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (aprx.c)
 * @return  None
 * @param   thee Pointer to an Aprx allocated memory location
 */
VEXTERNC void Aprx_buildNodes(Aprx *thee);

/**
 * @ingroup Aprx
 * @brief   Evaluate the trace information.   
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (aprx.c)
 * @return  None
 * @param   thee  Pointer to an Aprx allocated memory location
 * @param   Wud   block NODAL dirichlet work vector
 * @param   Wui   block NODAL interior dirichlet work vector
 * @param   Wut   block NODAL analytical work vector
 */
VEXTERNC void Aprx_evalTrace(Aprx *thee, Bvec *Wud, Bvec *Wui, Bvec *Wut);

/**
 * @ingroup Aprx
 * @brief   Refine the mesh.
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (aprx.c)
 * @return  Success enumeration
 * @param   thee Pointer to an Aprx allocated memory location
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
VEXTERNC int Aprx_refine(Aprx *thee, int rkey, int bkey, int pkey);

/**
 * @ingroup Aprx
 * @brief   Unrefine the mesh.
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (aprx.c)
 * @return  Success enumeration
 * @param   thee  Pointer to an Aprx allocated memory location
 * @param   rkey  Boolean sets the rkey type for un-refining the mesh.
 * @param   pkey  Boolean sets the pkey type to prolongate a vector
 */
VEXTERNC int Aprx_unRefine(Aprx *thee, int rkey, int pkey);

/**
 * @ingroup Aprx
 * @brief   Deform the mesh
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (aprx.c)
 * @return  Success enumeration for deforming the mesh.
 * @param   thee Pointer to an Aprx allocated memory location
 * @param   def  The corresponding deforming block vector
 */
VEXTERNC int Aprx_deform(Aprx *thee, Bvec *def);

/**
 * @ingroup Aprx
 * @brief   Create storage for and build prolongation operator; the operator 
 *          prolongates a vector FROM a coarser level TO this level.
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (aprx.c)
 * @return  None
 * @param   thee   Pointer to an Aprx allocated memory location
 * @param   numRc  number of vertices before refinement
 * @param   numR   number of vertices after refinement
 * @param   pkey   Boolean sets the pkey type to prolongate a vector
 */
VEXTERNC void Aprx_buildProlong(Aprx *thee, int numRc, int numR, int pkey);

/**
 * @ingroup Aprx
 * @brief   Create block version of prolongation operator; also 
 *          updates the Bnode object.    
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (aprx.c)
 * @return  None
 * @param   thee Pointer to an Aprx allocated memory location
 * @param   p    Pointer to Bnode matrix 
 */
VEXTERNC void Aprx_buildProlongBmat(Aprx *thee, Mat *p);

/**
 * @ingroup Aprx
 * @brief   Count up the total number of basis functions (nodes). 
 *          This is basically an a priori count of the number of rows in
 *          the stiffness matrix.  
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (aprx.c) \n
 *          This routine can handle COMPLETELY GENERAL ELEMENTS; see the
 *          comments in Aprx_interact().  
 * @return  the total number of basis functions (nodes).
 * @param   thee  Pointer to an Aprx allocated memory location
 * @param   re    a reference simplex element object.
 */
VEXTERNC int Aprx_nodeCount(Aprx *thee, Re *re);

/**
 * @ingroup Aprx
 * @brief   Count the number of basis (linear/quadratic/etc) functions which 
 *          interact with a given basis function.  This is basically an 
 *          a priori count of the number of nonzeros per row in the 
 *          stiffness matrix, and the number of nozeros per row above the
 *          diagonal (and thus the number of nonzeros per column below the
 *          diagonal if we have symmetric nonzero structure).  By doing this  
 *          count, we can one-time malloc storage for matrices before   
 *          assembly, and avoid extra copies as well as avoid the use  
 *          of linked list structures for representing matrices with 
 *          a priori unknown nonzero structure. 
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (aprx.c)\n
 *          This routine can handle COMPLETELY GENERAL ELEMENTS; all it needs
 *          to do is simply compute all interactions of all nodes in the 
 *          element.  The nodes may be any combination of vertex-based,
 *          edge-based, face-based (3D only), or simplex-based (interior)
 *          degrees of freedom.\n
 *          This calculation is made possible because all of the node numbers 
 *          for all of the vertex/edge/face/simplex-nodes in each simplex are 
 *          available in O(1) time for each element, using the geometry
 *          manager routine Gem_simplexInfo().
 * @return  None
 * @param   thee    Pointer to an Aprx allocated memory location
 * @param   re      a reference simplex element object.
 * @param   reT     a reference simplex element object.
 * @param   numR    num of rows in the matrix
 * @param   numO    num of nonzeros we are actually storing in the
 *                   strict upper-triangle of matrix. (DRC only)
 * @param   numOYR  YSMP-row(col)
 * @param   IJA     integer structure [ IA ; JA ] 
 * @param   IJAYR   lower-triangle block from upper triangle block
 */
VEXTERNC void Aprx_interact(Aprx *thee, Re *re, Re *reT,
    int *numR, int *numO, int *numOYR, int **IJA, int **IJAYR);

/**
 * @ingroup Aprx
 * @brief   Determine all basis function interations.
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (aprx.c)
 * @return  None
 * @param   thee    Pointer to an Aprx allocated memory location
 * @param   re      array of reference simplex element objects
 * @param   reT     array of reference simplex element objects
 * @param   sym     symmetry keys for the matrix
 * @param   mirror  block mirror keys for the matrix
 * @param   frmt    Matrix format
 * @param   numR    num of rows in the matrix
 * @param   numO    num of nonzeros we are actually storing in the
 *                   strict upper-triangle of matrix. (DRC only)
 * @param   IJA     integer structure [ IA ; JA ] 
 */
VEXTERNC void Aprx_interactBlock(Aprx *thee,  Re *re[MAXV], Re *reT[MAXV],
    MATsym sym[MAXV][MAXV], MATmirror mirror[MAXV][MAXV], MATformat frmt[MAXV][MAXV],
    int numR[MAXV], int numO[MAXV][MAXV], int *IJA[MAXV][MAXV]);

/**
 * @ingroup Aprx
 * @brief   Compute global index of node of dimension dim in the matrix   
 * @author  Michael Holst, Oleg Korobkin
 * @note    Class Aprx: Non-inlineable methods (aprx.c)
 * @return  global index of node of dimension dim in the matrix
 * @param   thee  Pointer to an Aprx allocated memory location
 * @param   re    Pointer to the reference elemented to be used
 * @param   i     i-simplex global number
 * @param   dim   interactions index
 * @param   q     number of volume quadrature points
 */
VEXTERNC int Aprx_nodeIndex (Aprx *thee, Re *re, int i, int dim, int q);

/**
 * @ingroup Aprx
 * @brief   Given an index of node DF within the simplex and  
 *          its TT structure, compute global index for that node 
 * @authors Michael Holst and Oleg Korobkin  
 * @note    Class Aprx: Non-inlineable methods (aprx.c)
 * @return  global index for the node by given an index of node DF within
 *          the simplex
 * @param   thee  Pointer to an Aprx allocated memory location
 * @param   re    Pointer to the reference elemented to be used
 * @param   t     Pointer to class T
 * @param   idf   index of node DF within the simplex
 */
VEXTERNC int Aprx_nodeIndexTT (Aprx *thee, Re *re, TT *t, int idf);

/**
 * @ingroup Aprx
 * @brief   Build boundary Bnode information into two Bmats.
 * @authors Michael Holst
 * @note    Class Aprx: Non-inlineable methods (aprx.c)\n
 *          Either/both Bmat objects can be VNULL without error;
 *          the result for that VNULL Bmat is a No-Op.
 * @return  None
 * @param   thee Pointer to an Aprx allocated memory location
 * @param   A    Pointer to stiffness matrix
 * @param   M    Pointer to Mass matrix
 */
VEXTERNC void Aprx_buildBRC(Aprx *thee, Bmat *A, Bmat *M);

/**
 * @ingroup Aprx
 * @brief   Initialize vector test function "i" and its gradient, 
 *          generated by scalar test function "r" for this element.
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (aprx.c)
 * @return  None
 * @param   thee  Pointer to an Aprx allocated memory location
 * @param   i     vector test function index 
 * @param   r     scalar test function index
 * @param   phi   values of local basis funcs at integ pt
 * @param   phix  first derivs of local basis funcs at integ pt
 * @param   V     vector test function
 * @param   dV    first derivs of vector test function
 */
VEXTERNC void Aprx_initSpace(Aprx *thee, int i, int r,
    double phi[], double phix[][3], double V[], double dV[][3]);


#if 0
/**
 * @ingroup Aprx
 * @brief   Build the scalar basis functions and their gradients   
 *          on this element.  Also build an interpolant of the
 *          solution and its gradient.\n
 *          evalKey==0 ==> use U=[0+ud] as evaluation point \n
 *          evalKey>=1 ==> use U=[u+ud] as evaluation point 
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (assem.c) 
 * @return  None
 * @param   thee     Pointer to an Aprx allocated memory location
 * @param   re       a reference simplex element object.
 * @param   evalKey  evalKey==0 ==> use U=[0+ud] as evaluation point \n
 *                   evalKey>=1 ==> use U=[u+ud] as evaluation point 
 * @param   qp       index of quad pts
 * @param   face     index of face types
 * @param   t        pointer to class T
 * @param   xq       coordinates of quad pts
 * @param   phi      values of local basis funcs at integ pt
 * @param   phix     first derivs of local basis funcs at integ pt
 * @param   U        block NODAL solution vector 
 * @param   dU       first derivs of block NODAL solution vector
 * @param   Wu       block NODAL solution work vector
 * @param   Wud      block NODAL dirichlet work vector
 */
VEXTERNC void Aprx_buildBasis(Aprx *thee, Re *re,
    int evalKey, int qp, int face, TT *t,
    double xq[], double phi[], double phix[][3],
    double U[], double dU[][3],
    Bvec *Wu, Bvec *Wud);
#endif

/**
 * @ingroup Aprx
 * @brief   Initialize an element matrix.  
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (assem.c) 
 * @param   thee  Pointer to an Aprx allocated memory location
 * @param   bumpKey    == 0 ==> Do not assemble with bumps.\n
 *                     == 1 ==> Assemble bilinear form with bumps.\n
 *                     == 2 ==> Assemble residual form with bumps.\n
 *                     == ? ==> Same as 1.\n
 * @param   t     Pointer to class T
 * @param   em    Pointer to element matrix
 */
VEXTERNC void Aprx_initEmat(Aprx *thee, int bumpKey, TT *t, Emat *em);

/**
 * @ingroup Aprx
 * @brief   Build an element matrix and an element load vector for a 
 *          single given element, using quadrature. 
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (assem.c) \n
 *          We handle both the volume and surface cases here, using
 *          "face" and "t" as follows: \n
 *          face <  0 ==> do volume quadrature; we take the node numbers 
 *          as just the identity mapping.\n
 *          face >= 0 ==> do face quadrature for face "face", where the node
 *          numbers for face are taken from t->loc[face][].
 * @return  None
 * @param   thee       Pointer to an Aprx allocated memory location
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
 * @param   face       face <  0 ==> do volume quadrature; we take the node 
 *                     numbers  as just the identity mapping.\n
 *                     face >= 0 ==> do face quadrature for face "face", where
 *                     the node numbers for face are taken from t->loc[face][].
 * @param   t          Pointer to class T
 * @param   em         Pointer to element matrix
 * @param   Wu         block NODAL solution work vector
 * @param   Wud        block NODAL dirichlet work vector
 */
VEXTERNC void Aprx_quadEmat(Aprx *thee,
    int evalKey, int energyKey, int residKey, int tangKey, int massKey,
    int bumpKey,
    int face, TT *t, Emat *em,
    Bvec *Wu, Bvec *Wud);

/**
 * @ingroup Aprx
 * @brief   Fan out an element matrix to the global matrix,
 *          and fan out an element vector to the global vector. 
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (assem.c) 
 * @return  None
 * @param   thee       Pointer to an Aprx allocated memory location
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
 * @param   em         Pointer to element matrix
 * @param   A          Pointer to stiffness matrix
 * @param   M          Pointer to Mass matrix
 * @param   F          Pointer to residual vector
 */
VEXTERNC void Aprx_fanEmat(Aprx *thee,
    int evalKey, int energyKey, int residKey, int tangKey, int massKey,
    int bumpKey,
    Emat *em,
    Bmat *A, Bmat *M, Bvec *F);

/**
 * @ingroup Aprx
 * @brief   Assemble one element matrix.  
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (assem.c) 
 * @return  None
 * @param   thee       Pointer to an Aprx allocated memory location
 * @param   sm         Pointer to a given simplex
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
 * @param   em         Pointer to element matrix
 * @param   Wu         block NODAL solution work vector
 * @param   Wud        block NODAL dirichlet work vector
 */
VEXTERNC void Aprx_assemEmat(Aprx *thee, SS *sm,
    int evalKey, int energyKey, int residKey, int tangKey, int massKey,
    int bumpKey,
    Emat *em,
    Bvec *Wu, Bvec *Wud);

/**
 * @ingroup Aprx
 * @brief   Assemble the tangent matrix and nonlinear residual    
 *          vector (which reduces to the normal stiffness matrix and load   
 *          vector in the case of a linear problem) for a general nonlinear 
 *          system of m components in spatial dimension d.
 * @note    Class Aprx: Non-inlineable methods (assem.c) \n
 *          The nonlinear residual should always be ZERO at dirichlet nodes!
 * @author  Michael Holst
 * @return  the assembled energy
 * @param   thee       Pointer to an Aprx allocated memory location
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
 * @param   ip         index for assembled energies 
 * @param   rp         parameter for initially assembled PDE
 * @param   A          Pointer to stiffness matrix
 * @param   M          Pointer to Mass matrix
 * @param   F          block NODAL residual vector
 * @param   Wu         block NODAL solution work vector
 * @param   Wud        block NODAL dirichlet work vector
 */
VEXTERNC double Aprx_assem(Aprx *thee,
    int evalKey, int energyKey, int residKey, int tangKey, int massKey,
    int bumpKey,
    int ip[], double rp[],
    Bmat *A, Bmat *M, Bvec *F, Bvec *Wu, Bvec *Wud);

/**
 * @ingroup Aprx
 * @brief   Mark simplices to be refined
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (estim.c)\n
 *          If (  (key == -1) || (key ==  0) || (key ==  1)) ==>\n
 *          Let geometry manager deal with it.\n
 *          If (key == 2) ==> Nonlinear residual indicator.\n
 *          If (key == 3) ==> Dual-weighted residual indicator.\n
 *          If (key == 4) ==> Local problem indicator.
 * @return  number of marked simplices
 * @param   thee    Pointer to an Aprx allocated memory location
 * @param   key     If (  (key == -1) || (key ==  0) || (key ==  1)) ==>\n
 *                  Let geometry manager deal with it.\n
 *                  If (key == 2) ==> Nonlinear residual indicator.\n
 *                  If (key == 3) ==> Dual-weighted residual indicator.\n
 *                  If (key == 4) ==> Local problem indicator.
 * @param   color   chart type of marked simplices
 * @param   bkey    Bisection type
 * @param   elevel  level to determine the fraction of simplices 
 * @param   u       block NODAL solution vector
 * @param   ud      block NODAL dirichlet vector
 * @param   f       block NODAL residual vector
 * @param   r       block NODAL temporary vector
 */
VEXTERNC int Aprx_markRefine(Aprx *thee, int key, int color,
    int bkey, double elevel, Bvec *u, Bvec *ud, Bvec *f, Bvec *r);

/**
 * @ingroup Aprx
 * @brief   A posteriori error estimation via several different approaches.
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (estim.c)
 * @verbatim
   Apply one of several a posteriori error estimators:   

        If (key == 2) ==> Nonlinear residual indicator.            
        If (key == 3) ==> Dual-weighted residual indicator.            
        If (key == 4) ==> Local problem indicator.                     
        If (key == 5) ==> Recovered gradient indicator.                

   The error tolerance can be used in several ways:                   

        If (bkey == 0) ==> mark if:  [error/S] > TOL                   
        If (bkey == 1) ==> mark if:  [error/S] > (TOL^2/numS)^{1/2}    

        If (bkey == 2) ==> mark fixed number/fraction (dep. on ETOL): 
            if elevel < 1, mark fraction of simplices                  
            if elevel >=1, mark that many simplices (rounded to        
            the smaller integer)                                           
   @endverbatim
 * @return  number of marked simplices
 * @param   thee    Pointer to an Aprx allocated memory location
 * @param   key     index of different types of a posteriori error estimators
 * @param   color   chart type of marked simplices
 * @param   bkey    Bisection type
 * @param   elevel  level to determine the fraction of simplices
 * @param   u       block NODAL solution vector
 * @param   ud      block NODAL dirichlet vector
 * @param   f       block NODAL residual vector
 * @param   r       block NODAL temporary vector
 */
VEXTERNC int Aprx_estRefine(Aprx *thee, int key, int color,
    int bkey, double elevel, Bvec *u, Bvec *ud, Bvec *f, Bvec *r);

/**
 * @ingroup Aprx
 * @brief   Selects a fixed number/fraction of simplicies to be refined.
 * @author  Olen Korobkin, Michael Holst
 * @note    Class Aprx: Non-inlineable methods (estim.c)
 * @verbatim
   This function is called by Aprx_markRefine if fixed number / 
   fraction of simplices <ilevel> has to be refined. Returns 
   number of marked elements.
   @endverbatim
 * @return  Number of marked elements
 * @param   thee    Pointer to an Aprx allocated memory location
 * @param   num2ref index for the refined simplex
 * @param   color   chart type of marked simplices
 */
VEXTERNC int Aprx_markRefineFixed (Aprx *thee, int num2ref, int color); 


/**
 * @ingroup Aprx
 * @brief   Selects a fixed number/fraction of simplicies to be refined.
 * @author  Olen Korobkin, Michael Holst
 * @note    Class Aprx: Non-inlineable methods (estim.c)
 * @verbatim
   This function is called by Aprx_markRefine. Returns 
   number of marked elements.
   @endverbatim
 * @return  Number of marked elements
 * @param   thee    Pointer to an Aprx allocated memory location
 * @param   percentToRefine Percentage to refine
 * @param   color   chart type of marked simplices
 */
VEXTERNC int Aprx_markRefineDorfler (Aprx *thee, double percentToRefine, int color); 


/**
 * @ingroup Aprx
 * @brief   Evaluate a finite element solution at a set of arbitrary points.
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (eval.c)
 * @return  None
 * @param   thee    Pointer to an Aprx allocated memory location
 * @param   u       solution vector
 * @param   block   index for the block
 * @param   numPts  number of all interpolation pts for one simplex
 * @param   pts     coordinates of the pt
 * @param   vals    solution for the interpolated pts
 * @param   marks   index for marked or not
 */
VEXTERNC void Aprx_evalFunc(Aprx *thee,
    Bvec *u, int block, int numPts, double *pts,
    double *vals, int *marks);

/**
 * @ingroup Aprx
 * @brief   Interpolates data onto finite difference grid 
 * @author  Michael Holst, Doug Arnold and Karen Camarda
 * @note    Class Aprx: Non-inlineable methods (eval.c)\n
 *          Much is lifted from Aprx_evalFunc,
 *          with some ideas from Doug Arnold.
 *          Contributed by Karen Camarda.
 * @verbatim
 * Example:  Below is a matlab example courtesy of Doug Arnold.
 *
 *    %-- BEGIN MATLAB CODE --------------------------------------------------
 *    
 *    function colormesh
 *    
 *    % This file demonstrates how to go from an unstructured grid to
 *    % a structured one.  Each point on the structured grid is colored
 *    % according to which triangle in the unstructured grid it belongs to.
 *    
 *    % Douglas N. Arnold, 3/4/98
 *    
 *    % set up a mesh
 *    npts=10;
 *    x=[ 1.13 1.66 0.50 0.46 0.26 0.20 0.70 1.93 1.01 0.65 ];
 *    y=[ 1.81 1.23 0.99 1.30 1.29 1.50 1.10 0.68 0.60 0.07 ];
 *    tri = delaunay(x,y);
 *    ntri=size(tri,1);
 *    
 *    % grid spacing
 *    h=.03;
 *    k=.04;
 *    
 *    clf
 *    [z,i]=sort(rand(64,1));
 *    colormap hsv;
 *    cmap = colormap;
 *    cmap=cmap(i,:,:);
 *    
 *    % loop over triangles
 *    for m = 1:ntri
 *      % get vertices of triangle i
 *      v1=[x(tri(m,1)),y(tri(m,1))];
 *      v2=[x(tri(m,2)),y(tri(m,2))];
 *      v3=[x(tri(m,3)),y(tri(m,3))];
 *      
 *      % plot the triangle
 *      ttt=[v1;v2;v3;v1];
 *      line(ttt(:,1),ttt(:,2),'color','green') 
 *    
 *      % find bounding box of triangle
 *      xmin = min([v1(1),v2(1),v3(1)]);
 *      xmax = max([v1(1),v2(1),v3(1)]);
 *      ymin = min([v1(2),v2(2),v3(2)]);
 *      ymax = max([v1(2),v2(2),v3(2)]);
 *    
 *      % loop through grid points in bounding box, plotting those in triangle
 *      for i = ceil(xmin/h):floor(xmax/h)
 *        for j = ceil(ymin/k):floor(ymax/k)
 *          if intriangle(v1,v2,v3,[i*h,j*k])
 *            % set color according to triangle number
 *        color = cmap(m,:);
 *            line(i*h,j*k,'marker','.','markersize',10,'color',color)
 *          end
 *        end
 *      end
 *    
 *    end
 *    
 *    function i = intriangle(v1,v2,v3,x)
 *    % INTRIANGLE -- check if a point belongs to a triangle
 *    %                
 *    % f = intriangle(v1,v2,v3,x)
 *    %
 *    % input:
 *    %   v1,v2,v3   2-vectors giving the vertices of the triangle
 *    %   x          2-vector giving the point
 *    %
 *    % output:
 *    %   i          1 if x belongs to the close triangle, 0 else
 *    
 *    % map point under affine mapping from given triangle to ref triangle
 *    xhat = (x-v1)/[v2-v1;v3-v1];
 *    % check if mapped point belongs to reference triangle
 *    i = xhat(1) >= 0 & xhat(2) >= 0 & 1-xhat(1)-xhat(2) >= 0;
 *    %-- END   MATLAB CODE --------------------------------------------------
 *
 * ***************************************************************************
   @endverbatim
 * @return  None   
 * @param   thee        Pointer to an Aprx allocated memory location
 * @param   u           solution vector
 * @param   block       index of block
 * @param   nx          number of grids in X coordinates
 * @param   ny          number of grids in y coordinates
 * @param   nz          number of grids in z coordinates
 * @param   x0          x coordinate of lower grid corner
 * @param   y0          y coordinate of lower grid corner
 * @param   z0          z coordinate of lower grid corner
 * @param   dx          Grid spacing in x direction
 * @param   dy          Grid spacing in y direction
 * @param   dz          Grid spacing in z direction
 * @param   derivs      index for different nfunc types
 * @param   outTypeKey  ASCII output or Binary output
 */
VEXTERNC void Aprx_fe2fd(Aprx *thee, Bvec *u, int block,
    int nx, int ny, int nz,
    double x0, double y0, double z0, double dx, double dy, double dz,
    int derivs, int outTypeKey);

/**
 * @ingroup Aprx
 * @brief   Evaluate a given boundary integral. 
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (eval.c)
 * @return  None
 * @param   thee  Pointer to an Aprx allocated memory location
 * @param   u     block NODAL solution vector
 * @param   ud    block NODAL dirichlet vector
 * @param   ut    block NODAL analytical vector
 */
VEXTERNC void Aprx_bndIntegral(Aprx *thee, Bvec *u, Bvec *ud, Bvec *ut);

/**
 * @ingroup Aprx
 * @brief   Evaluate ADM mass from conformal factor.        
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (eval.c)
 * @return  None
 * @param   thee   Pointer to an Aprx allocated memory location
 * @param   u      block NODAL solution vector
 * @param   ud     block NODAL dirichlet vector
 * @param   ut     block NODAL analytical vector
 * @param   block  index of a block
 */
VEXTERNC void Aprx_admMass(Aprx *thee, Bvec *u, Bvec *ud, Bvec *ut, int block);

/**
 * @ingroup Aprx
 * @brief   Evaluate error in solution in one of several norms. 
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (eval.c) \n
 *          We evaluate the error in the solution (in the case that 
 *          an exact analytical solution is available) in one of several 
 *          norms:\n
 *          If (key == 0) ==> L^2 norm of the error.\n
 *          If (key == 1) ==> L^{\\infty} norm of the error.\n
 *          If (key == 2) ==> H^1 norm of the error. 
 * @return  error in solution in one of several norms.
 * @param   thee    Pointer to an Aprx allocated memory location
 * @param   pcolor  simplex chart type
 * @param   key     If (key == 0) ==> L^2 norm of the error.\n
 *                  If (key == 1) ==> L^{\\infty} norm of the error.\n
 *                  If (key == 2) ==> H^1 norm of the error. 
 * @param   u       block NODAL solution vector
 * @param   ud      block NODAL dirichlet vector
 * @param   ut      block NODAL analytical vector
 */
VEXTERNC double Aprx_evalError(Aprx *thee, int pcolor, int key,
    Bvec *u, Bvec *ud, Bvec *ut);

/**
 * @ingroup Aprx
 * @brief   Evaluate difference of two functions in H1 norm in one element.
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (eval.c) \n
 *          We evaluate difference of two functions in H1 norm in one element.
 * @return  Difference of two functions in H1 norm in one element.
 * @param   thee    Pointer to an Aprx allocated memory location
 * @param   sm      simplex
 * @param   u       block NODAL solution vector
 * @param   ud      block NODAL dirichlet vector
 * @param   ut      block NODAL analytical vector
 */
VEXTERNC double Aprx_evalH1(Aprx *thee, SS *sm, Bvec *u, Bvec *ud, Bvec *ut);

/**
 * @ingroup Aprx
 * @brief   Write a finite element mesh or mesh function in some format. 
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (io.c)
 * @return  None
 * @param   thee      Pointer to an Aprx allocated memory location
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
 * @param   u         block NODAL solution vector
 * @param   format    GV/MATH format
 */
VEXTERNC void Aprx_writeGEOM(Aprx *thee, Vio *sock,
    int defKey, int colKey, int chartKey, double gluVal, int fkey,
    Bvec *u, char *format);

/**
 * @ingroup Aprx
 * @brief   Write a finite element mesh or mesh function in some format. 
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (io.c)
 * @return  None
 * @param   thee      Pointer to an Aprx allocated memory location
 * @param   sock      socket for reading the external mesh data (NULL otherwise)
 * @param   u         block NODAL solution vector
 * @param   format    GV/MATH format
 */
VEXTERNC void Aprx_writeSOL(Aprx *thee, Vio *sock, Bvec *u, char *format);

/**
 * @ingroup Aprx
 * @brief   Clear the partitioning (set all colors to pcolor).
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (parti.c)
 * @return  Success enumeration
 * @param   thee     Pointer to an Aprx allocated memory location
 * @param   pcolor   simplex chart type
 */
VEXTERNC int Aprx_partSet(Aprx *thee, int pcolor);

/**
 * @ingroup Aprx
 * @brief   Smooth the partitioning. 
 * @author  Michael Holst 
 * @note    Class Aprx: Non-inlineable methods (parti.c)
 * @return  counted number of partitions
 * @param   thee Pointer to an Aprx allocated memory location
 */
VEXTERNC int Aprx_partSmooth(Aprx *thee);

/**
 * @ingroup Aprx
 * @brief   Partition the domain using inertial bisection.
 *          Partition sets of points in R^d (d=2 or d=3) by viewing them 
 *          as point masses of a rigid body, and by then employing the
 *          classical mechanics ideas of inertia and Euler axes. 
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (parti.c)\n
 *          We first locate the center of mass, then change the coordinate
 *          system so that the center of mass is located at the origin.
 *          We then form the (symmetric) dxd inertia tensor, and then find 
 *          the set of (real) eigenvalues and (orthogonal) eigenvectors.
 *          The eigenvectors represent the principle inertial rotation axes,
 *          and the eigenvalues represent the inertial strength in those  
 *          principle directions.  The smallest inerial component along an
 *          axis represents a direction along which the rigid body is most 
 *          "line-like" (assuming all the points have the same mass). \n
 *          For our purposes, it makes sense to using the axis (eigenvector)  
 *          corresponding to the smallest inertia (eigenvalue) as the line to 
 *          bisect with a line (d=2) or a plane (d=3).  We know the center of 
 *          mass, and once we also have this particular eigenvector, we can 
 *          effectively bisect the point set into the two regions separated 
 *          by the line/plane simply by taking an inner-product of the    
 *          eigenvector with each point (or rather the 2- or 3-vector  
 *          representing the point).  A positive inner-product represents one 
 *          side of the cutting line/plane, and a negative inner-product 
 *          represents the other side (a zero inner-product is right on the 
 *          cutting line/plane, so we arbitrarily assign it to one region or
 *          the other). 
 * @return  Success enumeration
 * @param   thee    Pointer to an Aprx allocated memory location
 * @param   pcolor  simplex chart type
 * @param   numC    number of columns in the matrix
 * @param   evec    eigenvector 
 * @param   simH    simHelper class
 */
VEXTERNC int Aprx_partInert(Aprx *thee, int pcolor,
    int numC, double *evec, simHelper *simH);

/**
 * @ingroup Aprx
 * @brief   Partition the domain using spectral bisection.  
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (parti.c)\n
 * @verbatim
 *          We solve the following generalized eigenvalue problem:
 *
 *                Ax = lambda Bx
 *
 *           for second smallest eigenpair.  We then return the eigenvector
 *           from the pair.  We make this happen by turning it into a
 *           regular eigenvalue problem:
 *
 *               B^{-1/2} A B^{-1/2} ( B^{1/2} x ) = lambda ( B^{1/2} x )
 *
 *           or rather
 *
 *               C y = lambda y,   where C=B^{-1/2}AB^{-1/2}, y=B^{1/2}x.
 *
 *           The matrix "B" is simply a diagonal matrix with a (positive)
 *           error estimate for the element on the diagonal.  Therefore, the
 *           matrix B^{-1/2} is a well-defined positive diagonal matrix.
 *
 *           We explicitly form the matrix C and then send it to the inverse
 *           Rayleigh-quotient iteration to recover the second smallest
 *           eigenpair.  On return, we scale the eigenvector y by B^{-1/2} to
 *           recover the actual eigenvector x = B^{-1/2} y.
 *
 *           To handle the possibility that an element has zero error, in
 *           which case the B matrix had a zero on the diagonal, we set the
 *           corresponding entry in B^{-1/2} to be 1.
 * @endverbatim
 * @return  Success enumeration
 * @param   thee     Pointer to an Aprx allocated memory location
 * @param   pcolor   simplex chart type
 * @param   numC     number of columns in the matrix
 * @param   evec     eigenvector 
 * @param   simH     simHelper class
 * @param   ford     temporary storage ford array
 * @param   rord     temporary storage rord array
 * @param   general  index for creating different scaling matrix 
 */
VEXTERNC int Aprx_partSpect(Aprx *thee, int pcolor,
    int numC, double *evec, simHelper *simH, int *ford, int *rord,
    int general);

/**
 * @ingroup Aprx
 * @brief   Do recursive bisection. 
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (parti.c)
 * @return  Success enumeration
 * @param   thee  Pointer to an Aprx allocated memory location
 * @param   pkey  index for different partition option
 * @param   pwht  index for weighted partitioning
 * @param   ppow  partitioning steps
 */
VEXTERNC int Aprx_part(Aprx *thee, int pkey, int pwht, int ppow);

/**
 * @ingroup Aprx
 * @brief   Do a single bisection step.
 * @author  Michael Holst
 * @note    Class Aprx: Non-inlineable methods (parti.c)
 * @param   thee    Pointer to an Aprx allocated memory location
 * @param   pkey    partition method 
 * @param   pwht    partitioning weighting
 * @param   pcolor  simplex chart type
 * @param   poff    an index
 */
VEXTERNC int Aprx_partOne(Aprx *thee, int pkey, int pwht, int pcolor, int poff);

#endif /* _APRX_H_ */

