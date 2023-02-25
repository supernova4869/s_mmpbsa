/** 
 * @defgroup HBvec HBvec class
 * @brief    HBvec class
 */

/** 
 * @defgroup HBmat HBmat class
 * @brief    HBmat class
 */

/** 
 * @defgroup Bchar Bchar class
 * @brief    Bchar class
 */

/**
 *  @file       whb.h
 *  @ingroup    HBvec HBmat global_mc
 *  @brief      Class Whb: stabilized hierarchical basis library
 *  @author     Burak Aksoylu, Stephen Bond, and Michael Holst
 *  @note       None
 *  @version    $Id: whb.h,v 1.21 2010/08/12 05:19:26 fetk Exp $ 
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

#ifndef _WHB_H
#define _WHB_H

#include <mc/mc_base.h>

#include <mc/bam.h>

/**
 * @ingroup global_mc
 * @brief   Class Whb: Parameters and datatypes 
 */

typedef enum HBMATtype {
    ZERO_TYPE,
    AMIN_TYPE,
    ANOR_TYPE,
    GHB_TYPE,
    GWM_TYPE
} HBMATtype;

/* Class Whb: Definition */

/**
 * @ingroup HBmat
 * @brief   struct HBmat Definition
 * @author  Michael Holst
 */

struct sHBmat {

  /** @brief the memory manager                          */
    Vmem       *vmem;        
  /** @brief did i make vmem or was it inherited         */
    int        iMadeVmem;   

  /** @brief indicates which data structure is employed \n
   * ZERO: zero matrix, all blocks are null.\n
   * AMIN: coarsest level A, with A22 ptr = A11.\n
   * ANOR: multilevel A block, with A12/A21/A22.\n
   * GHB:  G matrix for HB, with only G21 block.\n
   * GWM:  G matrix for WMHB, with G12/G21/G22   */
    HBMATtype  type;

  /** @brief number of blocks in block system            */
    int        numB;   

  /** @brief (1,2)-block of a block matrix               */
    Bmat       *A12; 
  /** @brief (2,1)-block of a block matrix               */
    Bmat       *A21;  
  /** @brief (2,2)-block of a block matrix               */  
    Bmat       *A22;  

  /** @brief our coarser parent HBmat                    */
    struct sHBmat *next;    

};

/**
 * @ingroup HBmat
 * @brief   Declaration of the HBmat class as the HBmat structure
 * @author  Michael Holst
 */
typedef struct sHBmat HBmat;

/**
 * @ingroup HBvec
 * @brief   struct HBvec Definition
 * @author  Michael Holst
 */

struct sHBvec {

  /** @brief the memory manager                          */
    Vmem       *vmem;     
  /** @brief did i make vmem or was it inherited         */
    int        iMadeVmem; 
  /** @brief matrix state                                */
    int        state;

  /** @brief number of blocks in block system            */
    int        numB;  

  /** @brief pointer to the full block vector            */
    Bvec       *bv;   
  /** @brief pointer to the tail of a block vector       */
    Bvec       *bv2;  

  /** @brief our coarser parent HBvec                    */
    struct sHBvec *next;  

};

/**
 * @ingroup HBvec
 * @brief   Declaration of the HBvec class as the HBvec structure
 * @author  Michael Holst
 */
typedef struct sHBvec HBvec;

/**
 * @ingroup HBvec
 * @brief   struct Bchar Definition
 * @author  Michael Holst
 */

struct sBchar {
  /** @brief the memory manager                          */
    Vmem       *vmem;     
  /** @brief did i make vmem or was it inherited         */
    int        iMadeVmem;  

  /** @brief character string name for this array        */
    char       name[10];    
  /** @brief total number of rows                        */
    int        n;     
  /** @brief the character array components              */
    char       *uflat;  

  /** @brief next coarser object in the hierarchy        */
    struct sBchar *coarse;  

  /** @brief num character array blocks                  */
    int        numB;       
  /** @brief num of rows in each array block             */
    int        numR[MAXV];
  /** @brief array of pointers to blocks of characters   */ 
    char       *u[MAXV];  

};

/**
 * @ingroup Bchar
 * @brief   Declaration of the Bchar class as the Bchar structure
 * @author  Michael Holst
 */
typedef struct sBchar Bchar;

/*
 * ***************************************************************************
 * Class Whb: Inlineable methods (whb.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_WHB)
#else /* if defined(VINLINE_WHB) */
#endif /* if !defined(VINLINE_WHB) */

/*
 * ***************************************************************************
 * Class Whb: Non-Inlineable methods (whb.c)
 * ***************************************************************************
 */

/**
 * @ingroup HBvec
 * @brief   Hierarchical Basis MG Vcycle (Bank, Dupont, and Yserentant). 
 * @authors Burak Aksoylu and Stephen Bond 2002/01/09  
 * @note    Class Whb: Non-Inlineable methods (whb.c)\n
 * @verbatim
 * Warning:  This version is uses RESIDUAL/ERROR form!  The initial guess is
 *           always zero, and the righthand side is always the residual from
 *           the previous iteration.
 *
 *   Notes:  By using a residual/error form for the Vcycle, the code inside
 *           the Vcycle is much cleaner.  In addition, it's much simpler to
 *           decide when the residual tolerance has been met.  If a full
 *           approximation scheme is used, or if the residual is calculated
 *           "on the fly", the magnitude of the residual is not known until
 *           the very bottom of the Vcycle.  This makes Vcycle code a bit
 *           more complicated.
 *
 * HB crash course:
 *
 * A typical multigrid algorithm can be implemented in two different ways:
 * either as a full approximation scheme (FAS) or as a coarse grid correction
 * scheme (CGCS). CGCS seems to be more popular than FAS. Here, Alg_hbVcycle
 * is implemented as CGCS.  As far as we can tell, they both require the same
 * number of operations.
 *
 * Here is the MATLAB equivalent of HBvec_hbVcyc.
 *
 * function [u]=vhbmg(rhs,level) ;
 *
 * %%% global variables: prolongation and stiffness matrices
 * global P_12 P_23 P_34 P_45 P_56 P_67 P_78 P_89;
 * global A_1 A_2 A_3 A_4 A_5 A_6 A_7 A_8 A_9;
 *
 *   %%% get the stiffness matrix on this level
 *   A = eval(['A_' num2str(level) ]);
 *
 *  if (level == 1)
 *   u = A \ rhs;
 *  else
 *    %%%%%%% recover the dimensions
 *    P = eval(['P_' num2str(level-1) num2str(level)]);
 *    [N,M] = size(P);
 *    n = N-M;
 *    P2 = [zeros(M,n);eye(n,n)];
 *    myzeros = zeros(n,1);
 *
 *    %%% get the level-2 block of the stiffness matrix
 *    A_22 = A(N-n+1:N, N-n+1:N);
 *    %%% grab the fine part of rhs
 *    rhs_2 = P2'*rhs;
 *    %%% pre-smoothing by block symmetric Gauss-Seidel
 *    u_2 = smooth_block(A_22,myzeros,rhs_2);
 *    %%% transform basis
 *    w = P2*u_2;
 *    %%% restriction of the residual for the coarse grid
 *    res = P'*(rhs - A*w);
 *
 *        u = vhbmg(res,level-1);
 *
 *    %%% Prolongate the error
 *    %%% Be careful u is the error not the solution because of CGCR
 *    u = P*u;
 *    %%% update fine-grid pseudo residual. Pseudo residual because of
 *    %%% CGCS style. In fact, pseudoRes_2 = residual_2 + A_22*u_2.
 *    pseudoRes_2 = P2'*(rhs - A*u);
 *    %%% post-smoothing by block symmetric Gauss-Seidel
 *    %%% Rather than u_2 = smooth_block(A_22,myzeros,residual_2), we use
 *    u_2 = smooth_block(A_22,u_2,pseudoRes_2);
 *
 *    %%% embed the updated error to the solution by inclusion
 *    error = P2*u_2;
 *    u = u + error;
 *  end;
 * @endverbatim
 * @return  None
 * @param   dd    Pointer to the HBvec object
 * @param   Ahb   nodal stiffness matrix in doubly linked format
 * @param   rr    Pointer to the HBvec object
 * @param   Ghb   change of basis matrices for each level
 * @param   ww    Pointer to the HBvec object
 * @param   key   index for different solution methods: 0=(Au=f), 1=(A'u=f)
 * @param   csolv index for different HB options
 *          csolv==0 --> nonsymmetric HB\n
 *          csolv==1 --> adjoint of nonsymmetric HB\n
 *          csolv==2 --> symmetric HB\n
 */
VEXTERNC void HBvec_hbVcyc(HBvec *dd, HBmat *Ahb, HBvec *rr,
    HBmat *Ghb, HBvec *ww, int key, int csolv);

/**
 * @ingroup HBvec
 * @brief   Creates a chain of HBvec pointers to a Bvec.    
 * @authors Burak Aksoylu and Stephen Bond
 * @note    Class Whb: Non-Inlineable methods (whb.c)\n
 * @verbatim
 *   Notes:  This routine takes a null pointer to an HBvec array, and
 *           returns a linked (hierarchical) chain of fully assembled
 *           HBvec's.  Each of the links contains two Bvec pointers,
 *           which correspond to logical subsections of a global Bvec.
 *           The fully assemble HBvec array consist of only pointers,
 *           sharing the data space used by the specified Bvec.
 *
 *           Basic Inputs:
 *
 *           Bmat *P0, *P1, ...  the prologation matrices for each level
 *           Bvec *u             a block vector on the finest level
 *                               (passed using the HBvec bv pointer)
 *
 *           Basic Outputs (inside the HBvec chain):
 *
 *           Bvec *bv, *bv2     pointers to logical subsections of u.
 * @endverbatim
 * @return  None
 * @param   thee Pointer to the HBvec object
 * @param   Ahb  nodal stiffness matrix in doubly linked format
 */
VEXTERNC void HBvec_initStructure(HBvec *thee, HBmat *Ahb);

/**
 * @ingroup HBvec
 * @brief   Kills a chain of HBvec pointers.
 * @authors Burak Aksoylu and Stephen Bond
 * @note    Class Whb: Non-Inlineable methods (whb.c)
 * @return  None
 * @param   thee Pointer to the HBvec object
 */
VEXTERNC void HBvec_killStructure(HBvec *thee);

/**
 * @ingroup HBmat
 * @brief   Assembles the change of basis matrices for the HB Vcycle.
 * @author  Burak Aksoylu and Stephen Bond
 * @note    Class Whb: Non-Inlineable methods (whb.c)\n
 * @verbatim
 *   Notes:  This routine takes a pointer to an HBmat array, and
 *           returns a linked (hierarchical) chain of fully assembled
 *           HBmat's.  Each of the links contains several block matrices,
 *           which correspond to the logical subsections of the change
 *           of basis matrix:
 *
 *           [ S11 S12 ] = [  I  0  ] + [ G11 G12 ]
 *           [ S21 S22 ]   [  0  I  ]   [ G21 G22 ]
 *
 *           To keep optimal storage complexity, only G is stored, and
 *           any blocks which are entirely zero have VNULL pointers.
 *
 *           Basic Inputs:
 *
 *           BXLN *M   nodal mass matrix in doubly linked format
 *           Bmat *R     tails of the prolongation matrices (inside Algs)
 *           int minLev  coarsest level in the hierarchy
 *
 *           Basic Outputs (inside the HBmat chain):
 *
 *           Bmat *G21, *G22, *G12   logical blocks of G = S - I.
 *
 *           The mass matrix can be safely passed as VNULL in the HB
 *           case, since it is only touched in the WMHB case.
 * @endverbatim
 * @return  None
 * @param   thee  Pointer to the HBmat object
 * @param   Ppro  tails of the prolongation matrices (inside Algs)
 *                3 --> Standard HB\n
 *                5 --> Wavelet Modified HB
 * @param   Mlink nodal mass matrix in doubly linked format
 * @param   meth  method choice (0=Standard HB,1=Wavelet Modified HB)
 */
VEXTERNC void HBmat_initG(HBmat *thee, Bmat *Ppro, Bmat *Mlink, int meth);

/**
 * @ingroup HBmat
 * @brief   Assembles the hierarchical matrix for the HB Vcycle.
 * @author  Burak Aksoylu and Stephen Bond
 * @note    Class Whb: Non-Inlineable methods (whb.c)\n
 * @verbatim
 *   Notes:  This routine takes a pointer to an HBmat array, and
 *           returns a linked (hierarchical) chain of fully assembled
 *           HBmat's.  Each of the links contains several block matrices,
 *           which correspond to the stabilized subsections of the global
 *           block (stiffness) matrix.  All together, these stabilized
 *           blocks form the entire stiffness matrix, with respect to the
 *           recursively defined "hierarchical basis".
 *
 *           Basic Inputs:
 *
 *           HBmat *G    change of basis matrices for each level
 *           BXLN *A   nodal stiffness matrix in doubly linked format
 *           int minLev  coarsest level in the hierarchy
 *
 *           Basic Outputs (inside the HBmat chain):
 *
 *           Bmat *A21, *A22, *A12   stabilized Mat blocks on each level
 * @endverbatim
 * @return  None
 * @param   Ahb   nodal stiffness matrix in doubly linked format
 * @param   Ghb   change of basis matrices for each level
 * @param   Alink coarsest level in the hierarchy
 */
VEXTERNC void HBmat_initA(HBmat *Ahb, HBmat *Ghb, Bmat *Alink);

/**
 * @ingroup HBmat
 * @brief   Kills a chain of HBmat pointers.
 * @author  Burak Aksoylu and Stephen Bond
 * @note    Class Whb: Non-Inlineable methods (whb.c)
 * @return  None
 * @param   thee  Pointer to the HBmat object
 */
VEXTERNC void HBmat_killStructure(HBmat *thee);
/*
 * ***************************************************************************
 * Class Whb: Non-Inlineable methods (bvechb.c)
 * ***************************************************************************
 */

/**
 * @ingroup Bvec
 * @brief   Executes one of a number of hierarchical linear solvers. 
 * @author  Stephen Bond 2007/09/24   
 * @note    Class Whb: Non-Inlineable methods (bvechb.c)
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
 * @param   M      the mass matrix
 * @param   meth   method choice (0=Standard HB,1=Wavelet Modified HB)
 */
VEXTERNC void Bvec_hlmethod(Bvec *thee, Bmat *A, Bvec *f, Bvec *r, Bvec *ut,
    int key, int flag, int itmax, double etol, int prec, int cycle, Bmat *P,
    Bmat *M, int meth);

/**
 * @ingroup Bvec
 * @brief   HBMG algorithm due to Bank, Dupont, and Yserentant     
 * @author  Burak Aksoylu and Stephen Bond 
 * @note    Class Whb: Non-Inlineable methods (bvechb.c)\n
 *          csolv==0 --> nonsymmetric HB\n
 *          csolv==1 --> adjoint of nonsymmetric HB\n
 *          csolv==2 --> symmetric HB\n
 *          prec==3  --> Standard HB\n
 *          prec==5  --> Wavelet Modified HB
 * @return  None
 * @param   thee  pointer to the block vector
 * @param   A     system matrix
 * @param   f     source vector slot
 * @param   r     residual slot
 * @param   ut    block NODAL analytical solution
 * @param   key   0 --> solve Au=f\n
 *                1 --> solve A'u=f
 * @param   flag  0  --> Normal: check itmax, error tolerance; normal i/o\n
 *                1  --> Normal: check itmax; no error tolerance; no i/o
 * @param   itmax number of iterations to do (the maximum allowed)
 * @param   etol  error tolerance
 * @param   meth  method choice (0=Standard HB,1=Wavelet Modified HB)
 * @param   Ppro  tails of the prolongation matrices (inside Algs)
 *                3 --> Standard HB\n
 *                5 --> Wavelet Modified HB
 * @param   M     the mass matrix
 */
VEXTERNC void Bvec_hb(Bvec *thee, Bmat *A, Bvec *f, Bvec *r, Bvec *ut,
    int key, int flag, int itmax, double etol, int meth, Bmat *Ppro, Bmat *M);

/**
 * @ingroup Bvec
 * @brief   CG method with Bvec_hb as an optional preconditioner.
 * @author  Michael Holst and Stephen Bond 
 * @note    Class Whb: Non-Inlineable methods (bvechb.c)\n
 *          This is an exact copy of Bvec_cg with Bvec_mg -> Bvec_hb
 * @return  None
 * @param   thee  pointer to the block vector
 * @param   A     system matrix
 * @param   f     source vector slot
 * @param   r     residual slot
 * @param   ut    block NODAL analytical solution
 * @param   key   index for different solution methods: 0=(Au=f), 1=(A'u=f)
 * @param   flag  Determines which "mode" we run in:\n
 *                0 --> Normal: check itmax, error tolerance; normal i/o\n
 *                1 --> Silent: check only itmax; no i/o\n
 *                2 --> Subcycle: check itmax, error tolerance; subcycle i/o\n
 *                3 --> Subcycle: check itmax, error tolerance; no i/o
 * @param   itmax number of iterations to do (the maximum allowed)
 * @param   etol  error tolerance
 * @param   prec  index for different preconditioners
 * @param   P     prolongation matrix maintained by aprx
 * @param   M     the mass matrix
 * @param   p     the block vector p
 * @param   ap    another block vector ap
 * @param   bap   third block vector bap
 * @param   po    the block vector po
 * @param   apo   the block vector apo
 * @param   tp    temporary block vector tp
 */
VEXTERNC void Bvec_hcg(Bvec *thee, Bmat *A, Bvec *f, Bvec *r, Bvec *ut,
    int key, int flag, int itmax, double etol, int prec, Bmat *P, Bmat *M,
    Bvec *p, Bvec *ap, Bvec *bap, Bvec *po, Bvec *apo, Bvec *tp);

/**
 * @ingroup Bvec
 * @brief   Bvec_hbcg with Bvec_hb as an optional preconditioner.
 * @author  Michael Holst and Stephen Bond 
 * @note    Class Whb: Non-Inlineable methods (bvechb.c)\n
 *          This is an exact copy of Bvec_bcg with Bvec_mg -> Bvec_hb
 * @return  None
 * @param   thee  pointer to the block vector
 * @param   A     system matrix
 * @param   f     source vector slot
 * @param   r     residual slot
 * @param   ut    block NODAL analytical solution
 * @param   key   smooth with A or A' (0=A, 1=A')
 * @param   flag  Determines which "mode" we run in:\n
 *                0 --> Normal: check itmax, error tolerance; normal i/o\n
 *                1 --> Silent: check only itmax; no i/o\n
 *                2 --> Subcycle: check itmax, error tolerance; subcycle i/o\n
 *                3 --> Subcycle: check itmax, error tolerance; no i/o
 * @param   itmax number of iterations to do (the maximum allowed)
 * @param   etol  error tolerance
 * @param   prec  index for different preconditioners
 * @param   P     prolongation matrix maintained by aprx
 * @param   M     the mass matrix
 * @param   s     residual block vector s: f - A'u
 * @param   p     the block vector p
 * @param   ap    another block vector ap
 * @param   bap   third block vector bap
 * @param   po    the block vector po
 * @param   apo   the block vector apo
 * @param   q     the block vector q
 * @param   atq   the block vector atq
 * @param   btatq the block vector btatq
 * @param   qo    the block vector qo
 * @param   atqo  the block vector atqo
 */
VEXTERNC void Bvec_hbcg(Bvec *thee, Bmat *A, Bvec *f, Bvec *r, Bvec *ut,
    int key, int flag, int itmax, double etol, int prec, Bmat *P, Bmat *M,
    Bvec *s,
    Bvec *p, Bvec *ap, Bvec *bap, Bvec *po, Bvec *apo,
    Bvec *q, Bvec *atq, Bvec *btatq, Bvec *qo, Bvec *atqo);

/**
 * @ingroup HBmat
 * @brief   Assembles the hierarchical matrix structures for the HB Vcycle. 
 * @authors Burak Aksoylu and Stephen Bond
 * @note    Class Whb: Non-Inlineable methods (bvechb.c)\n
 * @verbatim
 *   Notes:  This routine takes a null pointer to two HBmat pointers and
 *           returns a linked (hierarchical) chain of fully assembled
 *           HBmat's.
 *
 *           Basic Inputs:
 *
 *           Bmat *A             the finest level stiffness matrix
 *           Bmat *P0, *P1, ...  the prologation matrices for each level
 *
 *           Basic Outputs (inside the HBmat chain):
 *
 *           Bmat *A21, *A22, *A12   stabilized Mat blocks on each level
 *           Bmat *A11               stabilized A11 block on coarse level
 *           Bvec *d1, *d2           defect section pointers
 *           Bvec *r1, *r2           residual section pointers
 *
 * Example:  Suppose we have a two-level system, resulting from a
 *           mesh refinement.  Let P be the prologation (restriction)
 *           matrix which prolongates (restricts) a trial solution
 *           from level 0 to 1 (1 to 0).  The stiffness matrix, A,
 *           can be naturally decomposed into four blocks, if the nodes
 *           added during refinement (fine nodes) are logically
 *           "appended" to the end of the solution vector:
 *
 *           A = [ A11 A12 ] ; P = [ I ] ; d = [ d1 ] ; r = [ r1 ]
 *               [ A21 A22 ]       [ R ]       [ d2 ]       [ r2 ]
 *
 *           Hierarchical Basis Multigrid (HB), doesn't use the standard
 *           nodal basis, using instead the "hierarchical basis".  The
 *           HB-basis stiffness matrix is related to the nodal basis
 *           stiffness matrix through a "change-of-basis" transform:
 *
 *           A_hb = S'*A*S ; d = S*d_hb ; r_hb = S'*r ; S = [ I 0 ]
 *                                                          [ R I ]
 *           Written out in blocks, A_hb has the form:
 *
 *           A_hb = [ A11 + A12*R + R'*A21 + R'*A22*R ; A12 + R'*A22 ]
 *                  [ A21 + A22*R                     ; A22          ]
 *
 *           The resulting HBmat and HBvec chains will have two levels,
 *           0 and 1, with a pointer to the chains returned as output.
 *           Level 0 contains the only non-null pointer to A11.
 * @endverbatim
 * @return  None
 * @param   Ahb  nodal stiffness matrix in doubly linked format
 * @param   Ghb  change of basis matrices for each level
 * @param   A    system matrix
 * @param   M    mass matrix
 * @param   Ppro tails of the prolongation matrices (inside Algs)
 *                3 --> Standard HB\n
 *                5 --> Wavelet Modified HB
 * @param   meth method choice (0=Standard HB,1=Wavelet Modified HB)
 */
VEXTERNC void HBmat_initMulti( HBmat **Ahb, HBmat **Ghb,
    Bmat *A, Bmat *M, Bmat *Ppro, int meth );

/**
 * @ingroup HBmat
 * @brief   Destruct the hierarchical structures
 * @authors Burak Aksoylu and Stephen Bond
 * @note    Class Whb: Non-Inlineable methods (bvechb.c)
 * @return  None
 * @param   Ahb   nodal stiffness matrix in doubly linked format
 * @param   Ghb   change of basis matrices for each level
 */
VEXTERNC void HBmat_killMulti( HBmat **Ahb, HBmat **Ghb );

/*
 * ***************************************************************************
 * Class Whb: Non-Inlineable methods (hbtools.c)
 * ***************************************************************************
 */

/**
 * @ingroup HBmat
 * @brief   The HBmat constructor. 
 * @author  Stephen Bond 2002/01/07   
 * @note    Class Whb: Non-Inlineable methods (hbtools.c)
 * @return  Pointer to the HBmat object
 * @param   vmem  Memory management object
 * @param   bname character string name for the matrix
 * @param   pnumB num vector blocks
 */
VEXTERNC HBmat* HBmat_ctor(Vmem *vmem, const char *bname, int pnumB);

/**
 * @ingroup HBmat
 * @brief   The HBmat destructor. 
 * @author  Stephen Bond 2002/01/07   
 * @note    Class Whb: Non-Inlineable methods (hbtools.c)
 * @return  None
 * @param   thee Pointer to the HBmat object
 */
VEXTERNC void HBmat_dtor(HBmat **thee);

/**
 * @ingroup HBvec
 * @brief   The HBvec constructor.
 * @author  Stephen Bond 2002/01/07   
 * @note    Class Whb: Non-Inlineable methods (hbtools.c)
 * @return  Pointer to the HBvec object
 * @param   vmem  Memory management object
 * @param   bname character string name for the block vector
 * @param   pnumB num vector blocks
 */
VEXTERNC HBvec* HBvec_ctor(Vmem *vmem, const char *bname, int pnumB);

/**
 * @ingroup HBvec
 * @brief   The HBvec destructor.
 * @author  Stephen Bond 2002/01/07   
 * @note    Class Whb: Non-Inlineable methods (hbtools.c)
 * @return  None
 * @param   thee Pointer to the HBvec object
 */
VEXTERNC void HBvec_dtor(HBvec **thee);

/**
 * @ingroup Bvec
 * @brief   Create a Bvec pointer to the tail of a Bvec.
 * @author  Stephen Bond 
 * @note    Class Whb: Non-Inlineable methods (hbtools.c)
 * @return  Pointer to the Bvec object
 * @param   v     Pointer to the Bvec object
 * @param   name  character string name for the block vector
 * @param   ibase pointer to size array
 * @param   numR  num of rows in each block vector
 */
VEXTERNC Bvec * Bvec_ctorPoint(Bvec *v, const char *name, int *ibase, int *numR);

/**
 * @ingroup Bvec
 * @brief   The block vec pointer destructor.
 * @author  Stephen Bond 
 * @note    Class Whb: Non-Inlineable methods (hbtools.c)
 * @return  None
 * @param   thee Pointer to Pointer of the Bvec object
 */
VEXTERNC void Bvec_dtorPoint(Bvec **thee);

/**
 * @ingroup Vec
 * @brief   Create a Vec pointer to the tail of a Vec.
 * @author  Stephen Bond 
 * @note    Class Whb: Non-Inlineable methods (hbtools.c)
 * @return  None
 * @param   v    Pointer to the vec object
 * @param   name character string name for the vector
 * @param   sep  pointer to size array
 * @param   numR num of rows in each vector
 */
VEXTERNC Vec * Vec_ctorPoint(Vec *v, const char *name, int sep, int numR);

/**
 * @ingroup Vec
 * @brief   The vec pointer destructor.
 * @author  Stephen Bond 
 * @note    Class Whb: Non-Inlineable methods (hbtools.c)
 * @return  None
 * @param   thee Pointer to Pointer of the Bvec object
 */
VEXTERNC void Vec_dtorPoint(Vec **thee);

/**
 * @ingroup Bmat
 * @brief   Create a Bmat pointer to the tail of a Bmat.  
 * @author  Stephen Bond 
 * @note    Class Whb: Non-Inlineable methods (hbtools.c)
 * @return  Bmat  pointer to the tail of a Bmat.
 * @param   P     prolongation matrix maintained by aprx
 * @param   name  character string name for the block matrix
 * @param   ibase pointer to size array
 * @param   numR  num of rows in each array block
 */
VEXTERNC Bmat *Bmat_ctorPoint(Bmat *P,
    const char *name, int *ibase, int *numR);

/**
 * @ingroup Bmat
 * @brief   The Bmat pointer destructor.
 * @author  Stephen Bond 
 * @note    Class Whb: Non-Inlineable methods (hbtools.c)
 * @return  None
 * @param   thee Pointer to Pointer of the Bmat object
 */
VEXTERNC void Bmat_dtorPoint(Bmat **thee);

/**
 * @ingroup Mat
 * @brief   Create a Mat pointer to the tail of a Pro.  
 * @author  Stephen Bond 
 * @note    Class Whb: Non-Inlineable methods (hbtools.c)
 * @return  Pointer to the Mat object
 * @param   P     prolongation matrix maintained by aprx
 * @param   name  character string name for the matrix
 * @param   ibase pointer to size array
 * @param   numR  num of rows in each array
 */
VEXTERNC Mat *Mat_ctorPoint(Mat *P, const char *name, int ibase, int numR);

/**
 * @ingroup Mat
 * @brief   The Mat pointer destructor.
 * @author  Stephen Bond 
 * @note    Class Whb: Non-Inlineable methods (hbtools.c)
 * @return  None
 * @param   thee Pointer to Pointer of the Mat object
 */
VEXTERNC void Mat_dtorPoint(Mat **thee);

/**
 * @ingroup HBmat
 * @brief   Print an HBmat in MATLAB sparse form.
 * @author  Stephen Bond 2002/08/11
 * @note    Class Whb: Non-Inlineable methods (hbtools.c)
 * @return  None
 * @param   thee  Pointer to the HBmat object
 * @param   fname character string name for the HBmatrix
 * @param   pflag index for write/append
 */
VEXTERNC void HBmat_printSp(HBmat *thee, char *fname, int pflag);

/**
 * @ingroup Bmat
 * @brief   Print the Dirichlet info of a Bmat in MATLAB sparse form. 
 * @author  Stephen Bond 2002/08/24
 * @note    Class Whb: Non-Inlineable methods (hbtools.c)
 * @return  None
 * @param   thee  Pointer to the Bmat object
 * @param   fname character string name for the block matrix
 * @param   pflag index for write/append
 */
VEXTERNC void Bmat_printDiriSp(Bmat *thee, char *fname, int pflag);

/**
 * @ingroup Mat
 * @brief   Print the Dirichlet info of a Bmat in MATLAB sparse form. 
 * @author  Stephen Bond 2002/08/24
 * @note    Class Whb: Non-Inlineable methods (hbtools.c)
 * @return  None
 * @param   thee  Pointer to the Mat object
 * @param   fname character string name for the matrix
 * @param   pflag index for write/append
 */
VEXTERNC void Mat_printDiriSp(Mat *thee, char *fname, int pflag);

/**
 * @ingroup Bmat
 * @brief   Return the Mat structure AD.  
 * @author  Michael Holst
 * @note    Class Whb: Non-Inlineable methods (hbtools.c)
 * @return  the Mat structure AD
 * @param   thee Pointer to the Bmat object
 * @param   p    index for the row
 * @param   q    index for the column
 */
VEXTERNC Mat *Bmat_AD(Bmat *thee, int p, int q);

/**
 * @ingroup HBvec
 * @brief   Special matrix-vector multiplication for HB.
 * @author  Michael Holst
 * @note    Class Whb: Non-Inlineable methods (hbtools.c)\n
 *          Blocks which are entirely zero should be specified using
 *          NULL pointers.  The work vector should be of the same length
 *          as the number of rows in G21. This vector is only used if 
 *          A12 != 0 or A22 != 0, otherwise it can be safely set to VNULL.
 * @return  None
 * @param   thee Pointer to the HBvec object
 * @param   Gmat Pointer to the source HBmat
 * @param   key  0   standard product, u += G *u\n
 *               1  transpose product, u += G'*u
 * @param   work Pointer to the block vector
 */
VEXTERNC void HBvec_matvec(HBvec *thee, HBmat *Gmat, int key, Bvec *work);

/**
 * @ingroup Bmat
 * @brief   Calculate the log of the determinant of a SPD Bmat matrix. 
 * @author  Stephen Bond
 * @note    Class Whb: Non-Inlineable methods (cholesky.c)\n
 *          Uses Cholesky factorization (slow for large matrices).
 * @return  None
 * @param   thee Pointer to the Bmat object
 */
VEXTERNC double Bmat_lnDet(Bmat *thee);

/*
 * ***************************************************************************
 * Class Whb: Non-Inlineable methods (brcmat.c)
 * ***************************************************************************
 */

/**
 * @ingroup BXLN
 * @brief   Create an BXLN copy of a Bmat.  
 * @author  Stephen Bond 2002/01/28
 * @note    Class Whb: Non-Inlineable methods (brcmat.c)\n
 *          Only the upper triangular part of any symmetric part is copied. 
 * @return  None
 * @param   thee Pointer to the Bmat object
 * @param   Amat Pointer to the Bmat object
 */
VEXTERNC void BXLN_copyBmat(Bmat *thee, Bmat *Amat);

/**
 * @ingroup Bmat
 * @brief   Produces a Bmat copy of a BXLN.
 * @author  Stephen Bond 2002/01/26
 * @note    Class Whb: Non-Inlineable methods (brcmat.c)
 * @return  None
 * @param   thee  Pointer to the Bmat object
 * @param   Alink Pointer to the Bmat object
 */
VEXTERNC void Bmat_copyBXLN(Bmat *thee, Bmat *Alink);

/**
 * @ingroup BXLN
 * @brief   Produce the HB sparse triple matrix product:
 * @verbatim
 *           A(p,q)_hb = (I + G(p)') A(p,q) (I + G(q))
 *
 *                  = [  I     Gp21'    ] [ Apq11 Apq12 ] [  I      Gq12   ]
 *                    [ Gp12' (I+Gp22') ] [ Apq21 Apq22 ] [ Gq21  (I+Gq22) ]
 * @endverbatim
 *          for every (p,q) block of A_hb.
 * @authors Stephen Bond and Burak Aksoylu 2002/01/26
 * @note    Class Whb: Non-Inlineable methods (brcmat.c)\n
 * @verbatim
 *   Notes:  The G21, G22, and G12 blocks must be ROW, DRC, and COL
 *           formats respectively.  Any NULL blocks are treated as
 *           zero matrices.  If the subblock of the BXLN is of symmetric
 *           type, only the upper triangular portion is used/modified.
 *
 *           The matrix A is overwritten with the A11 result.
 *
 *  Method:  CASE 1: Heirarchical Basis (HB)
 *
 *           [A11 A12] += [ A12 Gq21 + Gp21'(A21 + A22 Gq21); Gp21'A22 ]
 *           [A21 A22]    [                        A22 Gq21 ;  0       ]
 *
 *           CASE 2: Wavelet Modified Heirarchical Basis (WMHB)
 *
 *   [A11 A12] += [ Gp21'(A21+A22*Gq21) + A12*Gq21 ; 0 ]
 *   [A21 A22]    [ Gp22'(A21+A22*Gq21) + A22*Gq21 ; 0 ]
 *
 *              + [ 0                   ; Gp21'(A22+A22*Gq22+A21*Gq12) ]
 *                [ Gp21'(A11+A12*Gq21) ; Gp22'(A22+A22*Gq22+A21*Gq12) ]
 *
 *              + [ 0 ; A11*Gq12 + A12*Gq22                                ]
 *                [ 0 ; A21*Gq12 + A22*Gq22 + Gp12'(A12+A12*Gq22+A11*Gq12) ]
 * @endverbatim
 * @return  None
 * @param   thee Pointer to the Bmat object
 * @param   G    Pointer to the HBmat object
 */
VEXTERNC void BXLN_hbTriple(Bmat *thee, HBmat *G);

/**
 * @ingroup BLXN
 * @brief   Copies the A12, A21, and A22 blocks of a XLN, storing them
 *          in Mats of type COL, ROW, and DRC respectively. 
 * @authors Burak Aksoylu and Stephen Bond 2002/01/26
 * @note    Class Whb: Non-Inlineable methods (brcmat.c)\n
 * @verbatim
 *   Notes:  The logical position of eacg block is determined from the
 *           sizes of the other blocks.
 *
 *                    [   A11    |   A12   ]
 *                    [ sepxsep  | sepxn   ]
 *           A_hb  =  [----------|---------]  where numR = sep + n
 *                    [   A21    |   A22   ]
 *                    [   nxsep  |   nxn   ]
 *
 *           A12, A21, A22 must be passed w/ NULL_STATE blocks.
 * @endverbatim
 * @return  None
 * @param   A   Pointer to the Bmat object
 * @param   A12 Pointer to the Bmat object
 * @param   A21 Pointer to the Bmat object
 * @param   A22 Pointer to the Bmat object
 */
VEXTERNC void BXLN_copyBlocks(Bmat *A, Bmat *A12, Bmat *A21, Bmat *A22);

/**
 * @ingroup BLXN
 * @brief   Sets the sizes of the logical blocks inside a BXLN using
 *          the number of rows and cols in two different Bmats.
 * @authors Burak Aksoylu and Stephen Bond 2002/01/26
 * @note    Class Whb: Non-Inlineable methods (brcmat.c)\n
 * @verbatim
 *   Notes:  This method is used exclusively by HB methods to "replace"
 *           the matrix A with its A11 subblock.  The size of the A11
 *           sublock can be uniquely determined using the sizes of the A12
 *           and A21 sublocks:
 *
 *                 [   A11    |   A12   ]
 *                 [ sepxsep  | sepxn   ]
 *           A  =  [----------|---------]  where numR = sep + n
 *                 [   A21    |         ]
 *                 [   nxsep  |         ]
 *
 *           In the future, this routine may actually shrink the matrix,
 *           but for now it only sets the internal numR and numC arrays
 *           to reflect the size of the logical A11 block.
 * @endverbatim
 * @return  None
 * @param   A   Pointer to the Bmat object
 * @param   A12 Pointer to the Bmat object
 * @param   A21 Pointer to the Bmat object
 */
VEXTERNC void BXLN_shrinkLogical(Bmat *A, Bmat *A12, Bmat *A21);

/**
 * @ingroup HBmat
 * @brief   Constructs G12/G22 for WMHB using Mhb and G21.
 * @author  Stephen Bond 2002/02/03
 * @note    Class Whb: Non-Inlineable methods (brcmat.c)\n
 * @verbatim
 *   Notes:  This routine is the second step in the construction of the
 *           change of basis matrix for WMHB.  It constructs the G12 and
 *           G22 blocks using the mass matrix in the HB basis:
 *
 *           [ G11 ] = [ 0 ]  ;   [ G12 ] =  [ I ] inv(-Mhb_11) Mhb_12
 *           [ G21 ]   [ R ]      [ G22 ]    [ R ]
 *
 *           where Mhb is the previously formed mass matrix in the HB basis:
 *
 *           Mhb = [ Mhb_11 Mhb_12 ] = [ I R'] [ M_11 M_12 ] [ I 0 ]
 *                 [ Mhb_21 Mhb_22 ]   [ 0 I ] [ M_21 M_22 ] [ R I ]
 *
 *           Here inv( . ) represents any approximation to the inverse.
 *           To keep optimal complexity this should be chosen very carefully.
 *           Some possible choices include special matrix polynomials,
 *           mass lumping, or simply the inverse of the diagonal.
 *
 *           The G12 and G22 blocks should be passed in with NULL_STATE
 *           subblocks, with only their shape determined (i.e. just after
 *           ctor).  The G21 block must already exist, and contain the tail
 *           of the prolongation matrix, i.e. P = [I ; R].
 * @endverbatim
 * @return  None
 * @param   thee Pointer to the HBmat object
 * @param   M    Pointer to the Bmat object
 */
VEXTERNC void HBmat_GHB2WMHB(HBmat *thee, Bmat *M);

/*
 * ***************************************************************************
 * Class Whb: Non-Inlineable methods (rcmat.c)
 * ***************************************************************************
 */

/**
 * @ingroup XLN
 * @brief   Produce the HB sparse triple matrix product:
 * @verbatim
 *           A_hb = (I + GL') A (I + GR)
 *
 *                  = [  I     GL21'    ] [ A11 A12 ] [  I      GR12   ]
 *                    [ GL12' (I+GL22') ] [ A21 A22 ] [ GR21  (I+GR22) ]
 * @endverbatim
 *          where GL/GR are the left/right change of basis matrices.
 * @authors Stephen Bond and Burak Aksoylu 2002/01/26
 * @note    Class Whb: Non-Inlineable methods (rcmat.c)
 * @verbatim
 *   Notes:  The G21, G22, and G12 blocks must be ROW, DRC, and COL
 *           formats respectively.  Any NULL blocks are treated as
 *           zero matrices.  If the XLN is of symmetric type, only the
 *           upper triangular portion is used/modified.
 *
 *           The matrix A is overwritten with the A11 result.
 *
 *  Method:  CASE 1: Heirarchical Basis (HB)
 *
 *           [A11 A12] += [ A12 GL21 + GL21'(A21 + A22 GR21); GL21'A22 ]
 *           [A21 A22]    [                        A22 GR21 ;  0       ]
 *
 *           CASE 2: Wavelet Modified Heirarchical Basis (WMHB)
 *
 *   [A11 A12] += [ GL21'(A21+A22*GR21) + A12*GR21 ; 0 ]
 *   [A21 A22]    [ GL22'(A21+A22*GR21) + A22*GR21 ; 0 ]
 *
 *              + [ 0                   ; GL21'(A22+A22*GR22+A21*GR12) ]
 *                [ GL21'(A11+A12*GR21) ; GL22'(A22+A22*GR22+A21*GR12) ]
 *
 *              + [ 0 ; A11*GR12 + A12*GR22                                ]
 *                [ 0 ; A21*GR12 + A22*GR22 + GL12'(A12+A12*GR22+A11*GR12) ]
 * @endverbatim
 * @return  None
 * @param   thee Pointer to the Mat object
 * @param   GL21 Pointer to the Mat object
 * @param   GL12 Pointer to the Mat object
 * @param   GL22 Pointer to the Mat object
 * @param   GR21 Pointer to the Mat object
 * @param   GR12 Pointer to the Mat object
 * @param   GR22 Pointer to the Mat object
 */
VEXTERNC void XLN_hbTriple(Mat *thee,
    Mat *GL21, Mat *GL12, Mat *GL22, Mat *GR21, Mat *GR12, Mat *GR22 );

/**
 * @ingroup XLN
 * @brief   Multiplies two Mat's contributing the result to an XLN. 
 * @author  Stephen Bond
 * @note    Class Whb: Non-Inlineable methods (rcmat.c)\n
 *          The offset of the contrib is controlled by ibase, jbase.
 *          If flag1 (flag2) is nonzero, the first (second) matrix 
 *          is applied using its transpose.\n
 *          This routine is used by the WMHB method while forming
 *          the triple matrix product.
 * @return  None
 * @param   thee    Pointer to the Mat object
 * @param   Ablock1 Pointer to the Mat object
 * @param   ibase   pointer to size array
 * @param   flag1   index for transpose options
 * @param   Ablock2 Pointer to the Mat object
 * @param   jbase   pointer to size array
 * @param   flag2   index for transpose options
 */
VEXTERNC void XLN_matmatContrib(Mat *thee,
    Mat *Ablock1, int ibase, int flag1, Mat *Ablock2, int jbase, int flag2);

/**
 * @ingroup XLN
 * @brief   Copies the A12, A21, or A22 subblock of a XLN, storing it
 *          in a Mat of type COL, ROW, or DRC respectively.
 * @authors Burak Aksoylu and Stephen Bond 2002/01/26 
 * @note    Class Whb: Non-Inlineable methods (rcmat.c)
 * @verbatim
 *   Notes:  The logical position of each block is determined from the
 *           sizes of the other blocks.
 *
 *                    [   A11    |   A12   ]
 *                    [ sepxsep  | sepxn   ]
 *           A_hb  =  [----------|---------]  where numR = sep + n
 *                    [   A21    |   A22   ]
 *                    [   nxsep  |   nxn   ]
 *
 *           The subblock must be passed w/ NULL_STATE (i.e. after ctor).
 * @endverbatim
 * @return  None
 * @param   thee Pointer to the Mat object
 * @param   Amat Pointer to the Mat object
 * @param   flag index for different copy options
 */
VEXTERNC void XLN_copySubblock(Mat *thee, Mat *Amat, int flag);

/**
 * @ingroup XLN
 * @brief   Copies the A12, A21, and A22 blocks of a XLN, storing them
 *          in Mats of type COL, ROW, and DRC respectively.
 * @authors Burak Aksoylu and Stephen Bond 2002/01/26 
 * @note    Class Whb: Non-Inlineable methods (rcmat.c)
 * @verbatim
 *   Notes:  The logical position of each block is determined from the
 *           sizes of the other blocks.
 *
 *                    [   A11    |   A12   ]
 *                    [ sepxsep  | sepxn   ]
 *           A_hb  =  [----------|---------]  where numR = sep + n
 *                    [   A21    |   A22   ]
 *                    [   nxsep  |   nxn   ]
 *
 *           A12, A21, A22 must be passed w/ NULL_STATE (i.e. after ctor).
 * @endverbatim
 * @return  None
 * @param   Ablock Pointer to the Mat object
 * @param   A12    Pointer to the Mat object
 * @param   A21    Pointer to the Mat object
 * @param   A22    Pointer to the Mat object
 */
VEXTERNC void XLN_copyBlocks(Mat *Ablock, Mat *A12, Mat *A21, Mat *A22);

/**
 * @ingroup XLN
 * @brief   Sets the sizes of the logical blocks inside a XLN using
 *          the number of rows and cols in two different Mats.
 * @authors Burak Aksoylu and Stephen Bond 2002/01/26 
 * @note    Class Whb: Non-Inlineable methods (rcmat.c)
 * @verbatim
 *   Notes:  This method is used exclusively by HB methods to "replace"
 *           the matrix A with its A11 subblock.  The size of the A11
 *           sublock can be uniquely determined using the sizes of the A12
 *           and A21 sublocks:
 *
 *                 [   A11    |   A12   ]
 *                 [ sepxsep  | sepxn   ]
 *           A  =  [----------|---------]  where numR = sep + n
 *                 [   A21    |         ]
 *                 [   nxsep  |         ]
 *
 *           In the future, this routine may actually shrink the matrix,
 *           but for now it only sets the internal numR and numC values
 *           to reflect the size of the logical A11 block.
 * @endverbatim
 * @return  None
 * @param   A      Pointer to the Mat object
 * @param   A12    Pointer to the Mat object
 * @param   A21    Pointer to the Mat object
 */
VEXTERNC void XLN_shrinkLogical(Mat *A, Mat *A12, Mat *A21);

/**
 * @ingroup Mat
 * @brief   Constructs G12/G22 for WMHB using Mhb and G21. 
 * @author  Stephen Bond 2002/02/03
 * @note    Class Whb: Non-Inlineable methods (rcmat.c)
 * @verbatim
 *   Notes:  This routine is the second step in the construction of the
 *           change of basis matrix for WMHB.  It constructs the G12 and
 *           G22 blocks using the mass matrix in the HB basis:
 *
 *           [ G11 ] = [ 0 ]  ;   [ G12 ] =  [ I ] inv(-Mhb_11) Mhb_12
 *           [ G21 ]   [ R ]      [ G22 ]    [ R ]
 *
 *           where Mhb is the previously formed mass matrix in the HB basis:
 *
 *           Mhb = [ Mhb_11 Mhb_12 ] = [ I R'] [ M_11 M_12 ] [ I 0 ]
 *                 [ Mhb_21 Mhb_22 ]   [ 0 I ] [ M_21 M_22 ] [ R I ]
 *
 *           Here inv( . ) represents any approximation to the inverse.
 *           To keep optimal complexity this should be chosen very carefully.
 *           Some possible choices include special matrix polynomials,
 *           mass lumping, or simply the inverse of the diagonal.
 *
 *           The G12 and G22 blocks should be passed in with zero state,
 *           with only their shape determined (i.e. just after ctor).  The
 *           G21 block must already exist, and contain the tail of the
 *           prolongation matrix, R.
 * @endverbatim
 * @return  None
 * @param   G12    Pointer to the Mat object
 * @param   G22    Pointer to the Mat object
 * @param   G21    Pointer to the Mat object
 * @param   Mblock Pointer to the Mat object
 */
VEXTERNC void Mat_initGWMHB(Mat *G12, Mat *G22, Mat *G21, Mat *Mblock);

/*
 * ***************************************************************************
 * Class Whb: Non-Inlineable methods (subopt.c)
 * ***************************************************************************
 */

/**
 * @ingroup Bvec
 * @brief   Executes one of a number of linear solvers.
 * @authors Stephen Bond and Michael Holst
 * @note    Class Whb: Non-Inlineable methods (subopt.c)
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
 * @param   meth   method choice (0=Standard HB,1=Wavelet Modified HB)
 */
VEXTERNC void Bvec_submethod(Bvec *thee, Bmat *A, Bvec *f, Bvec *r, Bvec *ut,
    int key, int flag, int itmax, double etol, int prec, int cycle, Bmat *P,
    int meth);

/**
 * @ingroup Bvec
 * @brief   F-point smoothing operator
 * @authors Stephen Bond and Michael Holst
 * @note    Class Whb: Non-Inlineable methods (subopt.c)
 * @verbatim
 * Notes:    Consider the following partitioned linear system
 *
 *             [ Acc Acf ] [ xc ] = [ bc ]
 *             [ Afc Aff ] [ xf ]   [ bf ]
 *
 *           where xc is a vector of "coarse" or c-points and
 *           xf is a vector of "fine" or f-points.  This partitioning
 *           induces an analogous partitioning of the matrix A and the
 *           vector b as illustrated above. Given such a partitioning,
 *           an f-point smoother applies a Gauss-Seidel (or Jacobi)
 *           smoother to the f-points, leaving the c-points unchanged.
 *
 *           The forward Gauss-Seidel f-point smoother can be written as
 *               xc = xc
 *               xf = ( Dff + Lff )^(-1) ( bf - Afc xc - Uff xf )
 *
 *           The adjoint Gauss-Seidel f-point smoother can be written as
 *               xc = xc
 *               xf = ( Dff + Uff )^(-1) ( bf - Afc xc - Lff xf )
 *
 *           where, Aff = Uff + Lff + Dff is the standard gs splitting.
 *
 *           This method is invoked in the same manner as Bvec_smooth
 *           with one additional input.  This additional input is a char
 *           array of point types, with coarse and fine points indicated
 *           by 'c' and 'f' respectively.
 *
 *           For example, suppose we have 3 blocks with fpoints (1,5)
 *           in the first block, (0,3) in the second block, and (2,3,4)
 *           in the third block.  In this case fc should be
 *
 *           fc = {{'c', 'f', 'c', 'c', 'c', 'f'},
 *                    {'f', 'c', 'c', 'f', 'c', 'c'},
 *                    {'c', 'c', 'f', 'f', 'f', 'c'}}
 *
 *           This algorithm is not fast and its only practical purpose is
 *           for understanding and debugging HB and BPX style preconditioners.
 *           The HBMG algorithm is output-equivalent to standard
 *           multigrid with f-point smoothing where the f-points are the
 *           degrees of freedom introduced on the current level.  For
 *           multiplicative BPX, the list of f-points is expanded to include
 *           the one-ring of the f-points used by HB.
 *
 *           The optimal versions of multiplicative HB and BPX require
 *           the use of a change of basis so the f-point smoother does
 *           not need to touch the c-points explicitly.  However, the
 *           optimal version is mathematically "output-equivalent" to
 *           the suboptimal version described above.
 * @endverbatim
 * @return  None
 * @param   thee   pointer to the block vector                                               
 * @param   amat   system matrix
 * @param   f      source vector slot
 * @param   w      the work block vector
 * @param   key    smooth with A or A' (0=A, 1=A')
 * @param   ioflag debug output level (0=normal, 1=none, 2=lots, ... )  
 * @param   meth   smoother choice (0=jac, 1=gs, ... )   
 * @param   adj    adjoint choice (0=normal, 1=adjoint) 
 * @param   itmax  number of iterations to do (the maximum allowed)
 * @param   fc     Pointer to the block character                
 */
VEXTERNC void Bvec_fsmooth(Bvec *thee, Bmat *amat, Bvec *f, Bvec *w,
    int key, int ioflag, int meth, int adj, int itmax, Bchar *fc);

/**
 * @ingroup Bchar
 * @brief   The block character constructor
 * @authors Stephen Bond and Michael Holst
 * @note    Class Whb: Non-Inlineable methods (subopt.c)
 * @return  Pointer to a new allocated block character
 * @param   vmem  Memory management object
 * @param   name  character string name for the matrix
 * @param   pnumB num vector blocks
 * @param   pnumR num of rows in each block vector
 */
VEXTERNC Bchar* Bchar_ctor(Vmem *vmem,
    const char *name, int pnumB, int pnumR[MAXV]);

/**
 * @ingroup Bchar
 * @brief   The block character destructor
 * @authors Stephen Bond and Michael Holst
 * @note    Class Whb: Non-Inlineable methods (subopt.c)
 * @return  None
 * @param   thee Pointer to pointer of the block character
 */
VEXTERNC void Bchar_dtor(Bchar **thee);

/**
 * @ingroup Bchar
 * @brief   Build a block list of F and C points
 * @author  Stephen Bond 
 * @note    Class Whb: Non-Inlineable methods (subopt.c)
 * @return  None
 * @param   thee Pointer to the block character
 * @param   key  0 :  MG    -- All points are F-points\n
 *               1 :  HBMG  -- All new nodes are F-points\n
 *               2 :  BPXMG -- All new nodes plus 1-ring are F-points
 * @param   Ppro tails of the prolongation matrices (inside Algs)
 *                3 --> Standard HB\n
 *                5 --> Wavelet Modified HB
 */
VEXTERNC void Bchar_assem(Bchar *thee, int key, Bmat *Ppro);


/**
 * @ingroup Bchar
 * @brief   Build a block list of F and C points
 * @author  Stephen Bond 
 * @note    Class Whb: Non-Inlineable methods (subopt.c)
 * @return  None
 * @param   thee Pointer to the block character
 * @param   key  0 :  MG    -- All points are F-points\n
 *               1 :  HBMG  -- All new nodes are F-points\n
 *               2 :  BPXMG -- All new nodes plus 1-ring are F-points
 * @param   amat system matrix 
 * @param   Ppro tails of the prolongation matrices (inside Algs)
 *                3 --> Standard HB\n
 *                5 --> Wavelet Modified HB
 */
VEXTERNC void Bchar_assem2(Bchar *thee, int key, Bmat *amat, Bmat *Ppro);

#endif /* _WHB_H_ */
