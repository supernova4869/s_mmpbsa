/**
 * @defgroup VV VV class
 * @brief    VV class
 */

/**
 * @defgroup EE EE class
 * @brief    EE class
 */

/**
 * @defgroup FF FF class
 * @brief    FF class
 */

/**
 * @defgroup SS SS class
 * @brief    SS class
 */

/**
 *  @file       ves.h
 *  @ingroup    VV EE FF SS global_mc
 *  @brief      Class V,E,S: the fundamental simplex geometry objects.
 *  @author     Michael Holst
 *  @note       This is an implementation of the RIVER datastructure.
 *              RIVER = RInged-VERtex data Structure
 *  @version    $Id: ves.h,v 1.24 2010/08/12 05:19:07 fetk Exp $ 
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

#ifndef _VES_H_
#define _VES_H_

#include <mc/mc_base.h>

#include <mc/vel.h>

/*
 * ***************************************************************************
 * Class VV,VE,VS: Parameters and datatypes
 * ***************************************************************************
 */

/**
 * @ingroup VV
 * @brief   Class VV: Definition (8+32=40 bytes) 
 * @author  Michael Holst
 */
struct sVV {
  /** @brief --------> Geometric information package           */
    Gip      g;
  /** @brief --------> Vertex information package              */  
    Vip      d;
};

/**
 * @ingroup VV
 * @brief   Declaration of the VV class as the VV structure
 * @author  Michael Holst
 */
typedef struct sVV VV;

/**
 * @ingroup EE
 * @brief   Class EE: Definition (8+20=28 bytes) 
 * @author  Michael Holst
 */
struct sEE {
  /** @brief --------> Geometric information package           */
    Gip      g;
  /** @brief --------> Edge information package                */   
    Eip      d;
};

/**
 * @ingroup EE
 * @brief   Declaration of the EE class as the EE structure
 * @author  Michael Holst
 */
typedef struct sEE EE;

/**
 * @ingroup FF
 * @brief   Class FF: Definition (8+20=28 bytes) 
 * @author  Michael Holst
 */
struct sFF {
  /** @brief --------> Geometric information package           */
    Gip      g;
  /** @brief --------> Face information package                */  
    Fip      d;
};

/**
 * @ingroup FF
 * @brief   Declaration of the FF class as the FF structure
 * @author  Michael Holst
 */
typedef struct sFF FF; 

/**
 * @ingroup SS
 * @brief   Class SS: Definition (8+40=48 bytes OR 8+80=88 bytes) 
 * @author  Michael Holst
 */
struct sSS {
  /** @brief --------> Geometric information package           */
    Gip      g;   
  /** @brief --------> Simplex information package             */  
    Sip      d;   
};

/**
 * @ingroup SS
 * @brief   Declaration of the SS class as the SS structure
 * @author  Michael Holst
 */
typedef struct sSS SS; 


/**
 * @ingroup global_mc
 * @brief   Class T: Definition 
 * @author  Michael Holst
 */
typedef struct TT {

  /** @brief common chart for vertex coordinates        */
    int gchart;  
  /** @brief charts for vertex coordinates              */
    int chart[4]; 

  /** @brief jacobian determinant of the transformation */
    double D;  
  /** @brief jacobian determinant of cooked trans       */
    double Dcook;
  /** @brief face jacobian determinants                 */
    double faceD[4]; 

  /** @brief affine trans from master to arbitrary el   */
    double ff[3][3];
  /** @brief affine trans from master to arbitrary el   */
    double bb[3];
  /** @brief affine trans from arbitrary el to master   */ 
    double gg[3][3];
  /** @brief affine trans from arbitrary el to master   */ 
    double  cc[3];

  /** @brief local ordering of vertices for each face   */
    int loc[4][3];
  /** @brief vertex coordinate labels                   */
    double vx[4][3]; 
  /** @brief normal vectors to the faces                */
    double nvec[4][3];
  /** @brief edge vectors                               */ 
    double evec[6][3];
  /** @brief edge vector lengths                        */
    double elen[6];
  /** @brief baricenter                                 */   
    double bc[4];  

  /** @brief number of vertices in the d-simplex        */
    int dimV;
  /** @brief number of edges in the d-simplex           */
    int dimE;   
  /** @brief number of faces in the d-simplex           */ 
    int dimF; 
  /** @brief number of simplices in the d-simplex (=1)  */  
    int dimS; 

  /** @brief global simplex ID                          */
    int sid;
  /** @brief global vertex IDs                          */
    int vid[4];
  /** @brief global face IDs                            */
    int fid[4]; 
  /** @brief global edge IDs                            */ 
    int eid[6]; 

  /** @brief global simplex type                        */
    int stype;  
  /** @brief global vertex types                        */
    int vtype[4];  
  /** @brief global face types                          */
    int ftype[4];  
  /** @brief LOCAL edge types                           */ 
    int etype[6];  

  /** @brief pointer to the simplex                     */
    SS *s;  
  /** @brief pointers to vertices of the simplex        */
    VV *v[4];

  /** @brief local basis function of the current simplex */
    double phi[4], phix[4][3];

  /** @brief face and/or basis functions that are current focus */
    int focusFace, focusIV, focusIW;
} TT;

/*
 * ***************************************************************************
 * Class V: Inlineable methods (ves.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_GEM)
    /**
     * @ingroup VV
     * @brief   Intialize a vertex.
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a vertex
     * @param   dim  the dimension of the mesh
     * @param   id   the ID bits
     */
    VEXTERNC void VV_init(VV *thee, int dim, int id);

    /**
     * @ingroup VV
     * @brief   Re-intialize a vertex.
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a vertex
     */
    VEXTERNC void VV_reinit(VV *thee);

    /**
     * @ingroup VV
     * @brief   Set the reality
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a vertex
     * @param   type the reality bits
     */
    VEXTERNC void VV_setReality(VV *thee, int type);

    /**
     * @ingroup VV
     * @brief   Set the dimension
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a vertex
     * @param   dim  the dimension bits
     */
    VEXTERNC void VV_setDim(VV *thee, int dim);

    /**
     * @ingroup VV
     * @brief   Set the class.
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a vertex
     * @param   type the class bits
     */
    VEXTERNC void VV_setClass(VV *thee, int type);

    /**
     * @ingroup VV
     * @brief   Set the type
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a vertex
     * @param   type the type bits
     */
    VEXTERNC void VV_setType(VV *thee, int type);

    /**
     * @ingroup VV
     * @brief   Set the chart.
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee  Pointer to a vertex
     * @param   chart the chart bits
     */
    VEXTERNC void VV_setChart(VV *thee, int chart);

    /**
     * @ingroup VV
     * @brief   Set the ID
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a vertex
     * @param   id   the ID bits
     */
    VEXTERNC void VV_setId(VV *thee, int id);

    /**
     * @ingroup VV
     * @brief   Return the reality
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  the reality
     * @param   thee Pointer to a vertex
     */
    VEXTERNC unsigned int VV_reality(VV *thee);

    /**
     * @ingroup VV
     * @brief   Return the dimension
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  the dimension
     * @param   thee Pointer to a vertex
     */
    VEXTERNC unsigned int VV_dim(VV *thee);

    /**
     * @ingroup VV
     * @brief   Return the number of vertices in a simplex. 
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  the number of vertices in a simplex. 
     * @param   thee Pointer to a vertex
     */
    VEXTERNC unsigned int VV_dimVV(VV *thee);

    /**
     * @ingroup VV
     * @brief   Return the number of edges in a simplex.
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  the number of edges in a simplex.
     * @param   thee Pointer to a vertex
     */
    VEXTERNC unsigned int VV_dimEE(VV *thee);

    /**
     * @ingroup VV
     * @brief   Return the number of faces in a simplex.          
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  the number of faces in a simplex.          
     * @param   thee Pointer to a vertex
     */
    VEXTERNC unsigned int VV_dimFF(VV *thee);

    /**
     * @ingroup VV
     * @brief   Return the class
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  the class
     * @param   thee Pointer to a vertex
     */
    VEXTERNC unsigned int VV_class(VV *thee);

    /**
     * @ingroup VV
     * @brief   Return the type
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  the type
     * @param   thee Pointer to a vertex
     */
    VEXTERNC unsigned int VV_type(VV *thee);

    /**
     * @ingroup VV
     * @brief   Return the chart
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  the chart
     * @param   thee Pointer to a vertex
     */
    VEXTERNC unsigned int VV_chart(VV *thee);

    /**
     * @ingroup VV
     * @brief   Return the ID
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  the ID
     * @param   thee Pointer to a vertex
     */
    VEXTERNC unsigned int VV_id(VV *thee);

    /**
     * @ingroup VV
     * @brief   Set the coordinate
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a vertex
     * @param   i    index of the coordinate 
     * @param   val  the value of the coordinate
     */
    VEXTERNC void VV_setCoord(VV *thee, int i, double val);

    /**
     * @ingroup VV
     * @brief   Set the first edge in an edge ring.          
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a vertex
     * @param   eg   Pointer to an edge
     */
    VEXTERNC void VV_setFirstEE(VV *thee, EE *eg);

    /**
     * @ingroup VV
     * @brief   Set the first simplex in a simplex ring.     
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a vertex
     * @param   sm   Pointer to a simplex
     */
    VEXTERNC void VV_setFirstSS(VV *thee, SS *sm);

    /**
     * @ingroup VV
     * @brief   Set the parent edge pointer. 
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a vertex
     * @param   eg   Pointer to an edge
     */
    VEXTERNC void VV_setParentEE(VV *thee, EE *eg);

    /**
     * @ingroup VV
     * @brief   Return the coordinate
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  the coordinate
     * @param   thee Pointer to a vertex
     * @param   i    Index for a vertex in one simplex
     */
    VEXTERNC double VV_coord(VV *thee, int i);

    /**
     * @ingroup VV
     * @brief   Set the first edge in an edge ring. 
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  Pointer to the first edge in an edge ring.
     * @param   thee Pointer to a vertex
     */
    VEXTERNC EE* VV_firstEE(VV *thee);

    /**
     * @ingroup VV
     * @brief   Set the first simplex in a simplex ring.
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  Pointer to the first simplex in a simplex ring.
     * @param   thee Pointer to a vertex
     */
    VEXTERNC SS* VV_firstSS(VV *thee);

    /**
     * @ingroup VV
     * @brief   the parent edge pointer.
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  Pointer to the parent edge
     * @param   thee Pointer to a vertex
     */
    VEXTERNC EE* VV_parentEE(VV *thee);

    /**
     * @ingroup VV
     * @brief   Add an edge to an edge ring. 
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a vertex
     * @param   eg   Pointer to an edge
     */
    VEXTERNC void VV_addEdgeToRing(VV *thee, EE *eg);

    /**
     * @ingroup VV
     * @brief   Add a simplex to a simplex ring.   
     * @author  Michael Holst
     * @note    Class V: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a vertex
     * @param   sm   Pointer to a simplex
     */
    VEXTERNC void VV_addSimplexToRing(VV *thee, SS *sm);
#else /* if defined(VINLINE_GEM) */
/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_init(thee,dim,id) ( \
        (thee)->d.ePtr   = VNULL, \
        (thee)->d.sPtr   = VNULL, \
        (thee)->d.x[0]   = 0.0, \
        (thee)->d.x[1]   = 0.0, \
        (thee)->d.x[2]   = 0.0, \
        Vel_init((Vel*)(thee), dim, id) \
    )
/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_reinit(thee) ( \
        (thee)->d.ePtr   = VNULL, \
        (thee)->d.sPtr   = VNULL, \
        (thee)->d.x[0]   = 0.0, \
        (thee)->d.x[1]   = 0.0, \
        (thee)->d.x[2]   = 0.0, \
        Vel_reinit((Vel*)(thee)) \
    )

/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_setReality(thee,reel) Vel_setReality( (Vel*)(thee), (reel)  )
/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_setDim(thee,dim)      Vel_setDim(     (Vel*)(thee), (dim)   )
/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_setClass(thee,clas)   Vel_setClass(   (Vel*)(thee), (clas)  )
/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_setType(thee,type)    Vel_setType(    (Vel*)(thee), (type)  )
/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_setChart(thee,chart)  Vel_setChart(   (Vel*)(thee), (chart) )
/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_setId(thee,id)        Vel_setId(      (Vel*)(thee), (id)    )

/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_reality(thee)         Vel_reality(    (Vel*)(thee)          )
/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_dim(thee)             Vel_dim(        (Vel*)(thee)          )
/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_dimVV(thee)           Vel_dimVV(      (Vel*)(thee)          )
/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_dimEE(thee)           Vel_dimEE(      (Vel*)(thee)          )
/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_dimFF(thee)           Vel_dimFF(      (Vel*)(thee)          )
/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_class(thee)           Vel_class(      (Vel*)(thee)          )
/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_type(thee)            Vel_type(       (Vel*)(thee)          )
/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_chart(thee)           Vel_chart(      (Vel*)(thee)          )
/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_id(thee)              Vel_id(         (Vel*)(thee)          )

/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_setCoord(thee,i,val)  ((thee)->d.x[(i)] = (val))
/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_setFirstEE(thee,eg)   ((thee)->d.ePtr = (eg))
/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_setFirstSS(thee,eg)   ((thee)->d.sPtr = (eg))
/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_setParentEE(thee,eg)  ((thee)->d.eParent = (eg))

/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_coord(thee,i)  ((thee)->d.x[(i)])
/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_firstEE(thee)  (EE*)((thee)->d.ePtr)
/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_firstSS(thee)  (SS*)((thee)->d.sPtr)
/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_parentEE(thee) (EE*)((thee)->d.eParent)

/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_addEdgeToRing(thee,eg) ( \
        EE_setLink( (eg), (thee), VV_firstEE(thee) ), \
        VV_setFirstEE( (thee), (eg) ) \
    )
/** @brief Class V: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define VV_addSimplexToRing(thee,sm) ( \
        SS_setLink( (sm), (thee), VV_firstSS(thee) ), \
        VV_setFirstSS( (thee), (sm) ) \
    )
#endif /* if !defined(VINLINE_GEM) */

/*
 * ***************************************************************************
 * Class V: Non-inlineable methods (ves.c)
 * ***************************************************************************
 */

/**
 * @ingroup VV
 * @brief   The vertex constructor.
 * @author  Michael Holst
 * @note    Class V: Inlineable methods (ves.c) 
 * @return  Pointer to a new allocated vertex
 * @param   dim  the dimension bits
 * @param   myid the ID bits
 */
VEXTERNC VV* VV_ctor(int dim, int myid);

/**
 * @ingroup VV
 * @brief   The vertex destructor.
 * @author  Michael Holst
 * @note    Class V: Inlineable methods (ves.c) 
 * @return  None
 * @param   thee Pointer to a vertex
 */
VEXTERNC void VV_dtor(VV **thee);

/**
 * @ingroup VV
 * @brief   Remove an edge from an edge ring.   
 * @author  Michael Holst
 * @note    Class V: Inlineable methods (ves.c) 
 * @return  None
 * @param   thee Pointer to a vertex
 * @param   eg   Pointer to an edge
 */
VEXTERNC void VV_removeEdgeFromRing(VV *thee, EE *eg);

/**
 * @ingroup VV
 * @brief   Remove a simplex from a simplex ring.  
 * @author  Michael Holst
 * @note    Class V: Inlineable methods (ves.c) 
 * @return  None
 * @param   thee Pointer to a vertex
 * @param   sm   Pointer to a simplex
 */
VEXTERNC void VV_removeSimplexFromRing(VV *thee, SS *sm);

/**
 * @ingroup VV
 * @brief   Is an edge in the edge ring.
 * @author  Michael Holst
 * @note    Class V: Inlineable methods (ves.c) 
 * @return  Success enumeration
 * @param   thee Pointer to a vertex
 * @param   eg   Pointer to an edge
 */
VEXTERNC int VV_edgeInRing(VV *thee, EE *eg);

/**
 * @ingroup VV
 * @brief   Is a simplex in a simplex ring.    
 * @author  Michael Holst
 * @note    Class V: Inlineable methods (ves.c) 
 * @return  Success enumeration
 * @param   thee Pointer to a vertex
 * @param   sm   Pointer to a simplex
 */
VEXTERNC int VV_simplexInRing(VV *thee, SS *sm);

/**
 * @ingroup VV
 * @brief   Return the common edge between two vertices.
 * @author  Michael Holst
 * @note    Class V: Inlineable methods (ves.c) 
 * @return  the common edge between two vertices.
 * @param   thee Pointer to a vertex
 * @param   v0   Pointer to a vertex
 */
VEXTERNC EE* VV_commonEdge(VV *thee, VV *v0);

/**
 * @ingroup VV
 * @brief   Return parent edge having one of two vertices as midpoint. 
 * @author  Michael Holst
 * @note    Class V: Inlineable methods (ves.c) 
 * @return  parent edge having one of two vertices as midpoint. 
 * @param   thee Pointer to a vertex
 * @param   v0   Pointer to a vertex
 */
VEXTERNC EE* VV_parentEdge(VV *thee, VV *v0);

/**
 * @ingroup VV
 * @brief   Return a 2-simplex that is shared by 2 vertices, 
 *          and which is different from the input simplex sm.
 * @author  Michael Holst
 * @note    Class V: Inlineable methods (ves.c) 
 * @return  2-simplex that is shared by 2 vertices, 
 *          and which is different from the input simplex sm.
 * @param   thee Pointer to a vertex
 * @param   v0   Pointer to a vertex
 * @param   sm   Pointer to a simplex
 */
VEXTERNC SS* VV_commonSimplex2(VV *thee, VV *v0, SS *sm);

/**
 * @ingroup VV
 * @brief   Return a 3-simplex that is shared by 3 vertices, 
 *          and which is different from the input simplex sm.
 * @author  Michael Holst
 * @note    Class V: Inlineable methods (ves.c) 
 * @return  3-simplex that is shared by 3 vertices, 
 *          and which is different from the input simplex sm.
 * @param   thee Pointer to a vertex
 * @param   v0   Pointer to a vertex
 * @param   v1   Pointer to a vertex
 * @param   sm   Pointer to a smplex
 */
VEXTERNC SS* VV_commonSimplex3(VV *thee, VV *v0, VV *v1, SS *sm);

/**
 * @ingroup VV
 * @brief   Return a 3-simplex that is shared by 4 vertices, 
 *          and which is different from the input simplex sm.
 * @author  Michael Holst
 * @note    Class V: Inlineable methods (ves.c) 
 * @return  3-simplex that is shared by 4 vertices, 
 *          and which is different from the input simplex sm.
 * @param   thee Pointer to a vertex
 * @param   v0   Pointer to a vertex
 * @param   v1   Pointer to a vertex
 * @param   v2   Pointer to a vertex
 * @param   sm   Pointer to a simplex
 */
VEXTERNC SS* VV_commonSimplex4(VV *thee, VV *v0, VV *v1, VV *v2, SS *sm);

/**
 * @ingroup VV
 * @brief   the common vertex pointer
 * @author  Michael Holst
 * @note    Class V: Inlineable methods (ves.c) 
 * @return  Pointer to the common vertex
 * @param   thee Pointer to a vertex
 * @param   v0   Pointer to a vertex
 * @param   v1   Pointer to a vertex
 */
VEXTERNC VV* VV_commonVertex3(VV *thee, VV *v0, VV *v1);

/**
 * @ingroup VV
 * @brief   the common vertex pointer
 * @author  Michael Holst
 * @note    Class V: Inlineable methods (ves.c) 
 * @return  Pointer to the common vertex
 * @param   thee Pointer to a vertex
 * @param   v0   Pointer to a vertex
 * @param   v1   Pointer to a vertex
 * @param   v2   Pointer to a vertex
 */
VEXTERNC VV* VV_commonVertex4(VV *thee, VV *v0, VV *v1, VV *v2);


#if !defined(VINLINE_GEM)
    /**
     * @ingroup EE
     * @brief   Initialize an edge. 
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to an edge
     * @param   dim  the dimension bits
     * @param   id   the ID bits 
     */
    VEXTERNC void EE_init(EE *thee, int dim, int id);

    /**
     * @ingroup EE
     * @brief   Re-initialize an edge. 
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to an edge
     */
    VEXTERNC void EE_reinit(EE *thee);

    /**
     * @ingroup EE
     * @brief   Set the reality.
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to an edge
     * @param   type the reality bits
     */
    VEXTERNC void EE_setReality(EE *thee, int type);

    /**
     * @ingroup EE
     * @brief   Set the dimension. 
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to an edge
     * @param   dim  the dimension bits
     */
    VEXTERNC void EE_setDim(EE *thee, int dim);

    /**
     * @ingroup EE
     * @brief   Set the class.
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to an edge
     * @param   type the class bits
     */
    VEXTERNC void EE_setClass(EE *thee, int type);

    /**
     * @ingroup EE
     * @brief   Set the type.
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to an edge
     * @param   type the type bits
     */
    VEXTERNC void EE_setType(EE *thee, int type);

    /**
     * @ingroup EE
     * @brief   Set the chart.
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee  Pointer to a vertex
     * @param   chart the chart bits
     */
    VEXTERNC void EE_setChart(EE *thee, int chart);

    /**
     * @ingroup EE
     * @brief   Set the ID.
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to an edge
     * @param   id   the ID bits
     */
    VEXTERNC void EE_setId(EE *thee, int id);

    /**
     * @ingroup EE
     * @brief   Return the dimension.
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  the dimension.
     * @param   thee Pointer to an edge
     */
    VEXTERNC unsigned int EE_dim(EE *thee);

    /**
     * @ingroup EE
     * @brief   Return the number of vertices in a simplex. 
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  the number of vertices in a simplex. 
     * @param   thee Pointer to an edge
     */
    VEXTERNC unsigned int EE_dimVV(EE *thee);

    /**
     * @ingroup EE
     * @brief   Return the number of edges in a simplex.     
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  the number of edges in a simplex.     
     * @param   thee Pointer to an edge
     */
    VEXTERNC unsigned int EE_dimEE(EE *thee);

    /**
     * @ingroup EE
     * @brief   Return the number of faces in a simplex. 
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  the number of faces in a simplex.     
     * @param   thee Pointer to an edge
     */
    VEXTERNC unsigned int EE_dimFF(EE *thee);

    /**
     * @ingroup EE
     * @brief   Return the reality
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  the reality
     * @param   thee Pointer to an edge
     */
    VEXTERNC unsigned int EE_reality(EE *thee);

    /**
     * @ingroup EE
     * @brief   Return the class
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  the class
     * @param   thee Pointer to an edge
     */
    VEXTERNC unsigned int EE_class(EE *thee);

    /**
     * @ingroup EE
     * @brief   Return the type
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  the type
     * @param   thee Pointer to an edge
     */
    VEXTERNC unsigned int EE_type(EE *thee);

    /**
     * @ingroup EE
     * @brief   Return the chart
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  the chart
     * @param   thee Pointer to an edge
     */
    VEXTERNC unsigned int EE_chart(EE *thee);

    /**
     * @ingroup EE
     * @brief   Return the ID
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  the ID
     * @param   thee Pointer to an edge
     */
    VEXTERNC unsigned int EE_id(EE *thee);

    /**
     * @ingroup EE
     * @brief   Set one of the two vertices.   
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to an edge
     * @param   i    index for a vertex
     * @param   vx   Pointer to a vertex
     */
    VEXTERNC void EE_setVertex(EE *thee, int i, VV *vx);

    /**
     * @ingroup EE
     * @brief   Set the midpoint of the edge. 
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to an edge
     * @param   vx   Pointer to a vertex
     */
    VEXTERNC void EE_setMidPoint(EE *thee, VV *vx);

    /**
     * @ingroup EE
     * @brief   Set the parent edge pointer.
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to an edge
     * @param   eg   Pointer to an edge
     */
    VEXTERNC void EE_setParent(EE *thee, EE *eg);

    /**
     * @ingroup EE
     * @brief   Set the vertex order.  
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee  Pointer to an edge
     * @param   vxTmp Pointer to a vertex
     */
    VEXTERNC void EE_setVertexOrder(EE *thee, VV* vxTmp);

    /**
     * @ingroup EE
     * @brief   Set the next edge pointer associated with one of the vertices. 
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to an edge
     * @param   vx   Pointer to a vertex
     * @param   eg   Pointer to an edge
     */
    VEXTERNC void EE_setLink(EE *thee, VV *vx, EE *eg);

    /**
     * @ingroup EE
     * @brief   Return one of the vertices.
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  one of the vertcies.
     * @param   thee Pointer to an edge
     * @param   i    index for a vertex
     */
    VEXTERNC VV* EE_vertex(EE *thee, int i);

    /**
     * @ingroup EE
     * @brief   Return the midpoint of the edge.    
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  the midpoint of the edge.
     * @param   thee Pointer to an edge
     */
    VEXTERNC VV* EE_midPoint(EE *thee);

    /**
     * @ingroup EE
     * @brief   Return the parent edge pointer.  
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  the parent edge pointer.
     * @param   thee Pointer to an edge
     */
    VEXTERNC EE* EE_parent(EE *thee);

    /**
     * @ingroup EE
     * @brief   Return the vertex not matching the input vertex. 
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  the vertex not matching the input vertex. 
     * @param   thee Pointer to an edge
     * @param   vx   Pointer to a vertex
     */
    VEXTERNC VV* EE_otherVertex(EE *thee, VV *vx);

    /**
     * @ingroup EE
     * @brief   Return the next edge in the ring around a given vertex.    
     * @author  Michael Holst
     * @note    Class E: Inlineable methods (ves.c) 
     * @return  the next edge in the ring around a given vertex.    
     * @param   thee Pointer to an edge
     * @param   vx   Pointer to a vertex
     */
    VEXTERNC EE* EE_link(EE *thee, VV *vx);
#else /* if defined(VINLINE_GEM) */
/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_init(thee,dim,id) ( \
        (thee)->d.vPtr[0] = VNULL, \
        (thee)->d.vPtr[1] = VNULL, \
        (thee)->d.ePtr[0] = VNULL, \
        (thee)->d.ePtr[1] = VNULL, \
        (thee)->d.midPtr  = VNULL, \
        Vel_init((Vel*)(thee), dim, id) \
    )
/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_reinit(thee) ( \
        (thee)->d.vPtr[0] = VNULL, \
        (thee)->d.vPtr[1] = VNULL, \
        (thee)->d.ePtr[0] = VNULL, \
        (thee)->d.ePtr[1] = VNULL, \
        (thee)->d.midPtr  = VNULL, \
        Vel_reinit((Vel*)(thee)) \
    )

/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_setReality(thee,reel) Vel_setReality( (Vel*)(thee), (reel)  )
/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_setDim(thee,dim)      Vel_setDim(     (Vel*)(thee), (dim)   )
/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_setClass(thee,clas)   Vel_setClass(   (Vel*)(thee), (clas)  )
/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_setType(thee,type)    Vel_setType(    (Vel*)(thee), (type)  )
/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_setChart(thee,chart)  Vel_setChart(   (Vel*)(thee), (chart) )
/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_setId(thee,id)        Vel_setId(      (Vel*)(thee), (id)    )

/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_reality(thee)         Vel_reality(    (Vel*)(thee)          )
/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_dim(thee)             Vel_dim(        (Vel*)(thee)          )
/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_dimVV(thee)           Vel_dimVV(      (Vel*)(thee)          )
/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_dimEE(thee)           Vel_dimEE(      (Vel*)(thee)          )
/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_dimFF(thee)           Vel_dimFF(      (Vel*)(thee)          )
/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_class(thee)           Vel_class(      (Vel*)(thee)          )
/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_type(thee)            Vel_type(       (Vel*)(thee)          )
/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_chart(thee)           Vel_chart(      (Vel*)(thee)          )
/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_id(thee)              Vel_id(         (Vel*)(thee)          )

/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_setVertex(thee,i,vx) ((thee)->d.vPtr[(i)]=(vx))
/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_setMidPoint(thee,vx) ((thee)->d.midPtr=(vx))
/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_setParent(thee,eg)  ((thee)->d.eParent = (eg))
/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_setVertexOrder(thee,vxTmp) ( \
        (VV_id((thee)->d.vPtr[0]) > VV_id((thee)->d.vPtr[1])) \
            ? ( (vxTmp) = (thee)->d.vPtr[1], \
                (thee)->d.vPtr[1] = (thee)->d.vPtr[0], \
                (thee)->d.vPtr[0] = (vxTmp) ) \
            : VNULL \
    )
/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_setLink(thee,vx,eg) ( \
           ((vx) == (VV*)(thee)->d.vPtr[0]) ? (thee)->d.ePtr[0] = eg \
        : (((vx) == (VV*)(thee)->d.vPtr[1]) ? (thee)->d.ePtr[1] = eg : VNULL) \
    )
/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_vertex(thee,i) (VV*)((thee)->d.vPtr[(i)])
/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_midPoint(thee) (VV*)((thee)->d.midPtr)
/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_parent(thee)   (EE*)((thee)->d.eParent)
/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_otherVertex(thee,vx) (VV*)( \
           ((vx) == (VV*)(thee)->d.vPtr[0]) ? (thee)->d.vPtr[1] \
        : (((vx) == (VV*)(thee)->d.vPtr[1]) ? (thee)->d.vPtr[0] : VNULL) \
    )
/** @brief Class E: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#   define EE_link(thee,vx) (EE*)( \
           ((vx) == (VV*)(thee)->d.vPtr[0]) ? (thee)->d.ePtr[0] \
        : (((vx) == (VV*)(thee)->d.vPtr[1]) ? (thee)->d.ePtr[1] : VNULL) \
    )
#endif /* if !defined(VINLINE_GEM) */

/*
 * ***************************************************************************
 * Class E: Non-inlineable methods (ves.c)
 * ***************************************************************************
 */

/**
 * @ingroup EE
 * @brief   The edge constructor
 * @author  Michael Holst
 * @note    Class E: Inlineable methods (ves.c) 
 * @return  Pointer to a new allocated edge
 * @param   dim  the dimension bits
 * @param   myid the ID bits
 */
VEXTERNC EE* EE_ctor(int dim, int myid);

/**
 * @ingroup EE
 * @brief   The edge destructor.     
 * @author  Michael Holst
 * @note    Class E: Inlineable methods (ves.c) 
 * @return  None
 * @param   thee Pointer to an edge
 */
VEXTERNC void EE_dtor(EE **thee);

/**
 * @ingroup EE
 * @brief   Initialize the edge rings.
 * @author  Michael Holst
 * @note    Class E: Inlineable methods (ves.c) 
 * @return  None
 * @param   thee Pointer to an edge
 */
VEXTERNC void EE_initRing(EE *thee);

/**
 * @ingroup EE
 * @brief   Destroy the edge rings.   
 * @author  Michael Holst
 * @note    Class E: Inlineable methods (ves.c) 
 * @return  None
 * @param   thee Pointer to an edge
 */
VEXTERNC void EE_meltRing(EE *thee);

/**
 * @ingroup EE
 * @brief   Build the edge rings.      
 * @author  Michael Holst
 * @note    Class E: Inlineable methods (ves.c) 
 * @return  None
 * @param   thee Pointer to an edge
 */
VEXTERNC void EE_buildRing(EE *thee);


#if !defined(VINLINE_GEM)
    /**
     * @ingroup SS
     * @brief   Initialize a simplex.   
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a simplex
     * @param   dim  the dimension bits
     * @param   id   the ID bits
     */
    VEXTERNC void SS_init(SS *thee, int dim, int id);

    /**
     * @ingroup SS
     * @brief   Re-Initialize a simplex.   
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a simplex
     */
    VEXTERNC void SS_reinit(SS *thee);

    /**
     * @ingroup SS
     * @brief   Set the reality
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a simplex
     * @param   type the reality ID
     */
    VEXTERNC void SS_setReality(SS *thee, int type);

    /**
     * @ingroup SS
     * @brief   Set the dimension
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a simplex
     * @param   dim  the dimension bits
     */
    VEXTERNC void SS_setDim(SS *thee, int dim);

    /**
     * @ingroup SS
     * @brief   Set the class
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a simplex
     * @param   type the class bits
     */
    VEXTERNC void SS_setClass(SS *thee, int type);

    /**
     * @ingroup SS
     * @brief   Set the type
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a simplex
     * @param   type the type value of a simplex
     */
    VEXTERNC void SS_setType(SS *thee, int type);

    /**
     * @ingroup SS
     * @brief   Set the chart
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee  Pointer to a vertex
     * @param   chart index for a chart in a simplex
     */
    VEXTERNC void SS_setChart(SS *thee, int chart);

    /**
     * @ingroup SS
     * @brief   Set the ID
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a simplex
     * @param   id   index for the simplex ID
     */
    VEXTERNC void SS_setId(SS *thee, int id);

    /**
     * @ingroup SS
     * @brief   Return the reality
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  the reality
     * @param   thee Pointer to a simplex
     */
    VEXTERNC unsigned int SS_reality(SS *thee);

    /**
     * @ingroup SS
     * @brief   Return the dimension
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  the dimension
     * @param   thee Pointer to a simplex
     */
    VEXTERNC unsigned int SS_dim(SS *thee);

    /**
     * @ingroup SS
     * @brief   Return the number of vertices in a simplex.  
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  the number of vertices in a simplex.  
     * @param   thee Pointer to a simplex
     */
    VEXTERNC unsigned int SS_dimVV(SS *thee);

    /**
     * @ingroup SS
     * @brief   Return the number of edges in a simplex.  
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  the number of edges in a simplex.  
     * @param   thee Pointer to a simplex
     */
    VEXTERNC unsigned int SS_dimEE(SS *thee);

    /**
     * @ingroup SS
     * @brief   Return the number of faces in a simplex.  
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  the number of faces in a simplex.  
     * @param   thee Pointer to a simplex
     */
    VEXTERNC unsigned int SS_dimFF(SS *thee);

    /**
     * @ingroup SS
     * @brief   Return the class
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  the class
     * @param   thee Pointer to a simplex
     */
    VEXTERNC unsigned int SS_class(SS *thee);

    /**
     * @ingroup SS
     * @brief   Return the type
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  the type
     * @param   thee Pointer to a simplex
     */
    VEXTERNC unsigned int SS_type(SS *thee);

    /**
     * @ingroup SS
     * @brief   Return the chart
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  the chart
     * @param   thee Pointer to a simplex
     */
    VEXTERNC unsigned int SS_chart(SS *thee);

    /**
     * @ingroup SS
     * @brief   Return the ID
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  the ID
     * @param   thee Pointer to a simplex
     */
    VEXTERNC unsigned int SS_id(SS *thee);

    /**
     * @ingroup SS
     * @brief   Set the face type
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a simplex
     * @param   f    index for a face
     * @param   type number of the face type
     */
    VEXTERNC void SS_setFaceType(SS *thee, int f, int type);

    /**
     * @ingroup SS
     * @brief   Set the refinement edge.
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a simplex
     * @param   i    the refinement edge
     */
    VEXTERNC void SS_setRefinementEdge(SS *thee, int i);

    /**
     * @ingroup SS
     * @brief   Set the first marked edge.
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a simplex
     * @param   i    the first marked edge
     */
    VEXTERNC void SS_setMarkedEdge1(SS *thee, int i);

    /**
     * @ingroup SS
     * @brief   Set the second marked edge.
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a simplex
     * @param   i    the second marked edge
     */
    VEXTERNC void SS_setMarkedEdge2(SS *thee, int i);

    /**
     * @ingroup SS
     * @brief   Set the third marked edge.
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a simplex
     * @param   i    the third marked edge.
     */
    VEXTERNC void SS_setMarkedEdge3(SS *thee, int i);

    /**
     * @ingroup SS
     * @brief   Set degenerate marking flag.
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a simplex
     * @param   i    degenerate marking flag
     */
    VEXTERNC void SS_setDegen(SS *thee, int i);

    /**
     * @ingroup SS
     * @brief   Set refinement count.
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a simplex
     * @param   i    index for refinement count
     */
    VEXTERNC void SS_setRefinementCount(SS *thee, int i);

    /**
     * @ingroup SS
     * @brief   Set refinement key.
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee  Pointer to a simplex
     * @param   which index for the refinement
     * @param   key   the value for the refinement
     */
    VEXTERNC void SS_setRefineKey(SS *thee, int which, int key);
#   if defined(VG_ELEMENT)
        /**
	 * @ingroup SS
	 * @brief   the ID numer of a face
	 * @author  Michael Holst
	 * @note    Class S: Inlineable methods (ves.c) 
	 * @return  None
	 * @param   thee  Pointer to a simplex
	 * @param   i     index of a face
	 * @param   fn    the ID numer of a face
	 */
        VEXTERNC void SS_setFaceNumber(SS *thee, int i, int fn);

        /**
	 * @ingroup SS
	 * @brief   the ID numer of an edge
	 * @author  Michael Holst
	 * @note    Class S: Inlineable methods (ves.c) 
	 * @return  None
	 * @param   thee  Pointer to a simplex
	 * @param   i     index for an edge
	 * @param   fn    the ID numer of an edge
	 */
        VEXTERNC void SS_setEdgeNumber(SS *thee, int i, int en);
#   endif
    /**
     * @ingroup SS
     * @brief   Set a given vertex.
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a simplex
     * @param   i    index of a vertex
     * @param   vx   Pointer to a vertex
     */
    VEXTERNC void SS_setVertex(SS *thee, int i, VV *vx);

    /**
     * @ingroup SS
     * @brief   Set the next simplex pointer for a given simplex ring.
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a simplex
     * @param   vx   Pointer to a vertex
     * @param   sm   Pointer to a simplex
     */
    VEXTERNC void SS_setLink(SS *thee, VV *vx, SS *sm);

    /**
     * @ingroup SS
     * @brief   Set the type of a given face.
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  None
     * @param   thee Pointer to a simplex
     * @param   f    index of a face
     */
    VEXTERNC unsigned int SS_faceType(SS *thee, int f);

    /**
     * @ingroup SS
     * @brief   Return the refinement edge. 
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  the refinement edge.
     * @param   thee Pointer to a simplex
     */
    VEXTERNC unsigned int SS_refinementEdge(SS *thee);

    /**
     * @ingroup SS
     * @brief   Return the first marked edge.    
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  the first marked edge.
     * @param   thee Pointer to a simplex
     */
    VEXTERNC unsigned int SS_markedEdge1(SS *thee);

    /**
     * @ingroup SS
     * @brief   Return the second marked edge.   
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  the second marked edge.
     * @param   thee Pointer to a simplex
     */
    VEXTERNC unsigned int SS_markedEdge2(SS *thee);

    /**
     * @ingroup SS
     * @brief   Return the third marked edge.   
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  the third marked edge.
     * @param   thee Pointer to a simplex
     */
    VEXTERNC unsigned int SS_markedEdge3(SS *thee);

    /**
     * @ingroup SS
     * @brief   Return the degenerate edge marker.
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  the degenerate edge marker.
     * @param   thee Pointer to a simplex
     */
    VEXTERNC unsigned int SS_degen(SS *thee);

    /**
     * @ingroup SS
     * @brief   Return the refinement count.
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  the refinement count.
     * @param   thee Pointer to a simplex
     */
    VEXTERNC unsigned int SS_refinementCount(SS *thee);

    /**
     * @ingroup SS
     * @brief   Return the refinement key.
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  the refinement key.
     * @param   thee  Pointer to a simplex
     * @param   which index for refinement nor not
     */
    VEXTERNC unsigned int SS_refineKey(SS *thee, int which);

#   if defined(VG_ELEMENT)
       /**
        * @ingroup SS
        * @brief   the ID number of a face
        * @author  Michael Holst
        * @note    Class S: Inlineable methods (ves.c) 
        * @return  the ID number of a face
        * @param   thee  Pointer to a simplex
        * @param   i     index for a face
        */
        VEXTERNC int SS_faceNumber(SS *thee, int i);

       /**
        * @ingroup SS
        * @brief   the ID number of an edge
        * @author  Michael Holst
        * @note    Class S: Inlineable methods (ves.c) 
        * @return  the ID number of an edge
        * @param   thee  Pointer to a simplex
        * @param   i     index for an edge
        */
        VEXTERNC int SS_edgeNumber(SS *thee, int i);
#   endif

    /**
     * @ingroup SS
     * @brief   Return a given vertex.
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  a given vertex.
     * @param   thee Pointer to a simplex
     * @param   i    index for a vertex in a simplex
     */
    VEXTERNC VV* SS_vertex(SS *thee, int i);

    /**
     * @ingroup SS
     * @brief   Return the next simplex in a given simplex ring. 
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  the next simplex in a given simplex ring. 
     * @param   thee Pointer to a simplex
     * @param   vx   Pointer to a vertex
     */
    VEXTERNC SS* SS_link(SS *thee, VV *vx);

    /**
     * @ingroup SS
     * @brief   Return local vertex number i for face f.  
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  local vertex number i for face f.  
     * @param   thee Pointer to a simplex
     * @param   f    index of the face
     * @param   i    local vertex number
     */
    VEXTERNC int SS_faceVertexNumber(SS *thee, int f, int i);

    /**
     * @ingroup SS
     * @brief   Return the local vertex number associated with a given vertex.
     * @author  Michael Holst
     * @note    Class S: Inlineable methods (ves.c) 
     * @return  the local vertex number associated with a given vertex.
     * @param   thee Pointer to a simplex
     * @param   vx   Pointer to a vertex
     */
    VEXTERNC int SS_vptr2localVnum(SS *thee, VV *vx);
#else /* if defined(VINLINE_GEM) */
#   if defined(VG_ELEMENT)
/** @brief Class S: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#       define SS_init(thee,dim,id) ( \
            (thee)->d.vPtr[0] = VNULL, \
            (thee)->d.vPtr[1] = VNULL, \
            (thee)->d.vPtr[2] = VNULL, \
            (thee)->d.vPtr[3] = VNULL, \
            (thee)->d.sPtr[0] = VNULL, \
            (thee)->d.sPtr[1] = VNULL, \
            (thee)->d.sPtr[2] = VNULL, \
            (thee)->d.sPtr[3] = VNULL, \
            (thee)->d.fNum[0] = -1, \
            (thee)->d.fNum[1] = -1, \
            (thee)->d.fNum[2] = -1, \
            (thee)->d.fNum[3] = -1, \
            (thee)->d.eNum[0] = -1, \
            (thee)->d.eNum[1] = -1, \
            (thee)->d.eNum[2] = -1, \
            (thee)->d.eNum[3] = -1, \
            (thee)->d.eNum[4] = -1, \
            (thee)->d.eNum[5] = -1, \
            (thee)->d.tags    = MASK_00000000000000000000000000000000, \
            Vel_init((Vel*)(thee), dim, id) \
        )
/** @brief Class S: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#       define SS_reinit(thee) ( \
            (thee)->d.vPtr[0] = VNULL, \
            (thee)->d.vPtr[1] = VNULL, \
            (thee)->d.vPtr[2] = VNULL, \
            (thee)->d.vPtr[3] = VNULL, \
            (thee)->d.sPtr[0] = VNULL, \
            (thee)->d.sPtr[1] = VNULL, \
            (thee)->d.sPtr[2] = VNULL, \
            (thee)->d.sPtr[3] = VNULL, \
            (thee)->d.fNum[0] = -1, \
            (thee)->d.fNum[1] = -1, \
            (thee)->d.fNum[2] = -1, \
            (thee)->d.fNum[3] = -1, \
            (thee)->d.eNum[0] = -1, \
            (thee)->d.eNum[1] = -1, \
            (thee)->d.eNum[2] = -1, \
            (thee)->d.eNum[3] = -1, \
            (thee)->d.eNum[4] = -1, \
            (thee)->d.eNum[5] = -1, \
            (thee)->d.tags    = MASK_00000000000000000000000000000000, \
            Vel_reinit((Vel*)(thee)) \
        )
#   else
/** @brief Class S: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#       define SS_init(thee,dim,id) ( \
            (thee)->d.vPtr[0] = VNULL, \
            (thee)->d.vPtr[1] = VNULL, \
            (thee)->d.vPtr[2] = VNULL, \
            (thee)->d.vPtr[3] = VNULL, \
            (thee)->d.sPtr[0] = VNULL, \
            (thee)->d.sPtr[1] = VNULL, \
            (thee)->d.sPtr[2] = VNULL, \
            (thee)->d.sPtr[3] = VNULL, \
            (thee)->d.tags    = MASK_00000000000000000000000000000000, \
            Vel_init((Vel*)(thee), dim, id) \
        )
/** @brief Class S: Inlineable methods (ves.c) if defined(VINLINE_GEM) */
#       define SS_reinit(thee) ( \
            (thee)->d.vPtr[0] = VNULL, \
            (thee)->d.vPtr[1] = VNULL, \
            (thee)->d.vPtr[2] = VNULL, \
            (thee)->d.vPtr[3] = VNULL, \
            (thee)->d.sPtr[0] = VNULL, \
            (thee)->d.sPtr[1] = VNULL, \
            (thee)->d.sPtr[2] = VNULL, \
            (thee)->d.sPtr[3] = VNULL, \
            (thee)->d.tags    = MASK_00000000000000000000000000000000, \
            Vel_reinit((Vel*)(thee)) \
        )
#   endif

/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_setReality(thee,reel) Vel_setReality( (Vel*)(thee), (reel)  )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_setDim(thee,dim)      Vel_setDim(     (Vel*)(thee), (dim)   )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_setClass(thee,clas)   Vel_setClass(   (Vel*)(thee), (clas)  )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_setType(thee,type)    Vel_setType(    (Vel*)(thee), (type)  )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_setChart(thee,chart)  Vel_setChart(   (Vel*)(thee), (chart) )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_setId(thee,id)        Vel_setId(      (Vel*)(thee), (id)    )

/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_reality(thee)         Vel_reality(    (Vel*)(thee)          )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_dim(thee)             Vel_dim(        (Vel*)(thee)          )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_dimVV(thee)           Vel_dimVV(      (Vel*)(thee)          )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_dimEE(thee)           Vel_dimEE(      (Vel*)(thee)          )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_dimFF(thee)           Vel_dimFF(      (Vel*)(thee)          )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_class(thee)           Vel_class(      (Vel*)(thee)          )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_type(thee)            Vel_type(       (Vel*)(thee)          )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_chart(thee)           Vel_chart(      (Vel*)(thee)          )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_id(thee)              Vel_id(         (Vel*)(thee)          )

/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_setFaceType(thee,f,type) (void)( \
           ((f) == 0) ? \
              ( (thee)->d.faces &= MASK_11111111111111111111111100000000, \
                (thee)->d.faces |= (type) ) \
        : (((f) == 1) ? \
              ( (thee)->d.faces &= MASK_11111111111111110000000011111111, \
                (thee)->d.faces |= ((type) << 8) ) \
        : (((f) == 2) ? \
              ( (thee)->d.faces &= MASK_11111111000000001111111111111111, \
                (thee)->d.faces |= ((type) << 16) ) \
        : (((f) == 3) ? \
              ( (thee)->d.faces &= MASK_00000000111111111111111111111111, \
                (thee)->d.faces |= ((type) << 24) ) \
        : 0))) \
    )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_setRefinementEdge(thee,i) ( \
        (thee)->d.tags &= MASK_11111111111111111111111111111000, \
        (thee)->d.tags |=  (i) \
    )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_setMarkedEdge1(thee,i) ( \
        (thee)->d.tags &= MASK_11111111111111111111111111000111, \
        (thee)->d.tags |= ((i) << 3) \
    )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_setMarkedEdge2(thee,i) ( \
        (thee)->d.tags &= MASK_11111111111111111111111000111111, \
        (thee)->d.tags |= ((i) << 6) \
    )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_setMarkedEdge3(thee,i) ( \
        (thee)->d.tags &= MASK_11111111111111111111000111111111, \
        (thee)->d.tags |= ((i) << 9) \
    )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_setDegen(thee,i) ( \
        (thee)->d.tags &= MASK_11111111111111111000111111111111, \
        (thee)->d.tags |= ((i) << 12) \
    )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_setRefinementCount(thee,i) ( \
        (thee)->d.tags &= MASK_11111111111111000111111111111111, \
        (thee)->d.tags |= ((i) << 15) \
    )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_setRefineKey(thee,which,key) (void)( \
        ((which) == 0) \
            ? ( (thee)->d.tags &= MASK_11111111111110111111111111111111, \
                (thee)->d.tags |= ((key) << 18) ) \
            : ( (thee)->d.tags &= MASK_11111111111101111111111111111111, \
                (thee)->d.tags |= ((key) << 19) ) \
    )
#   if defined(VG_ELEMENT)
/** @brief Class S: Inlineable methods (ves.c) */
#       define SS_setFaceNumber(thee,i,fn) ((thee)->d.fNum[(i)] = (fn))
/** @brief Class S: Inlineable methods (ves.c) */
#       define SS_setEdgeNumber(thee,i,en) ((thee)->d.eNum[(i)] = (en))
#   endif
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_setVertex(thee,i,vx) ((thee)->d.vPtr[(i)] = (vx))
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_setLink(thee,vx,sm) ( \
           ((vx) == (VV*)(thee)->d.vPtr[0]) ? (thee)->d.sPtr[0] = sm \
        : (((vx) == (VV*)(thee)->d.vPtr[1]) ? (thee)->d.sPtr[1] = sm \
        : (((vx) == (VV*)(thee)->d.vPtr[2]) ? (thee)->d.sPtr[2] = sm \
        : (((vx) == (VV*)(thee)->d.vPtr[3]) ? (thee)->d.sPtr[3] = sm \
        : VNULL))) \
    )

/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_faceType(thee,f) ( \
           ((f) == 0) ? \
               ((thee)->d.faces & MASK_00000000000000000000000011111111) \
        : (((f) == 1) ? \
               ((thee)->d.faces & MASK_00000000000000001111111100000000) >> 8 \
        : (((f) == 2) ? \
               ((thee)->d.faces & MASK_00000000111111110000000000000000) >> 16 \
        : (((f) == 3) ? \
               ((thee)->d.faces & MASK_11111111000000000000000000000000) >> 24 \
        : 0))) \
    )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_refinementEdge(thee) ( \
         ((thee)->d.tags & MASK_00000000000000000000000000000111) \
    )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_markedEdge1(thee) ( \
        (((thee)->d.tags & MASK_00000000000000000000000000111000) >> 3) \
    )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_markedEdge2(thee) ( \
        (((thee)->d.tags & MASK_00000000000000000000000111000000) >> 6) \
    )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_markedEdge3(thee) ( \
        (((thee)->d.tags & MASK_00000000000000000000111000000000) >> 9) \
    )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_degen(thee) ( \
        (((thee)->d.tags & MASK_00000000000000000111000000000000) >> 12) \
    )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_refinementCount(thee) ( \
        (((thee)->d.tags & MASK_00000000000000111000000000000000) >> 15) \
    )
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_refineKey(thee,which) ( \
        (which == 0) \
          ? (((thee)->d.tags & MASK_00000000000001000000000000000000) >> 18) \
          : (((thee)->d.tags & MASK_00000000000010000000000000000000) >> 19) \
    )
#   if defined(VG_ELEMENT)
/** @brief Class S: Inlineable methods (ves.c) */
#       define SS_faceNumber(thee,i) ((thee)->d.fNum[(i)])
/** @brief Class S: Inlineable methods (ves.c) */
#       define SS_edgeNumber(thee,i) ((thee)->d.eNum[(i)])
#   endif
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_vertex(thee,i) (VV*)((thee)->d.vPtr[(i)])
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_link(thee,vx) (SS*)( \
           ((vx) == (VV*)((thee)->d.vPtr[0])) \
                  ?      ((thee)->d.sPtr[0]) \
        : (((vx) == (VV*)((thee)->d.vPtr[1])) \
                  ?      ((thee)->d.sPtr[1]) \
        : (((vx) == (VV*)((thee)->d.vPtr[2])) \
                  ?      ((thee)->d.sPtr[2]) \
        : (((vx) == (VV*)((thee)->d.vPtr[3])) \
                  ?      ((thee)->d.sPtr[3]) : VNULL))) \
    )

/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_faceVertexNumber(thee,f,i) (vmapF[(f)][(i)])
/** @brief Class S: Inlineable methods (ves.c) */
#   define SS_vptr2localVnum(thee,vx) ( \
           ((vx) == (VV*)(thee)->d.vPtr[0]) ? 0 \
        : (((vx) == (VV*)(thee)->d.vPtr[1]) ? 1 \
        : (((vx) == (VV*)(thee)->d.vPtr[2]) ? 2 \
        : (((vx) == (VV*)(thee)->d.vPtr[3]) ? 3 \
        : 0))) \
    )
#endif /* if !defined(VINLINE_GEM) */

/*
 * ***************************************************************************
 * Class S: Non-inlineable methods (ves.c)
 * ***************************************************************************
 */

/**
 * @ingroup SS
 * @brief   The simplex constructor.
 * @author  Michael Holst
 * @note    Class S: Inlineable methods (ves.c) 
 * @return  Pointer to a new allocated simplex object
 * @param   dim  the dimension bits
 * @param   myid the ID bits
 */
VEXTERNC SS* SS_ctor(int dim, int myid);

/**
 * @ingroup SS
 * @brief   The simplex destructor.
 * @author  Michael Holst
 * @note    Class S: Inlineable methods (ves.c) 
 * @return  None
 * @param   thee Pointer to a simplex
 */
VEXTERNC void SS_dtor(SS **thee);

/**
 * @ingroup SS
 * @brief   Initialize the simplex rings.
 * @author  Michael Holst
 * @note    Class S: Inlineable methods (ves.c)
 * @return  None
 * @param   thee Pointer to a simplex
 */
VEXTERNC void SS_initRing(SS *thee);

/**
 * @ingroup SS
 * @brief   Destroy the simplex rings.  
 * @author  Michael Holst
 * @note    Class S: Inlineable methods (ves.c) 
 * @return  None
 * @param   thee Pointer to a simplex
 */
VEXTERNC void SS_meltRing(SS *thee);

/**
 * @ingroup SS
 * @brief   Build the simplex rings.     
 * @author  Michael Holst
 * @note    Class S: Inlineable methods (ves.c) 
 * @return  None
 * @param   thee Pointer to a simplex
 */
VEXTERNC void SS_buildRing(SS *thee);

/**
 * @ingroup SS
 * @brief   Return the simplex sharing face i (opposite vertex i).   
 * @author  Michael Holst
 * @note    Class S: Inlineable methods (ves.c) 
 * @return  the simplex sharing face i (opposite vertex i).   
 * @param   thee Pointer to a simplex
 * @param   i    index for the face
 */
VEXTERNC SS* SS_nabor(SS *thee, int i);

/**
 * @ingroup SS
 * @brief   Return the edge connecting vertices i and j.      
 * @author  Michael Holst
 * @note    Class S: Inlineable methods (ves.c) 
 * @return  the edge connecting vertices i and j.      
 * @param   thee Pointer to a simplex
 * @param   i    index for a vertex
 * @param   j    index for a vertex
 */
VEXTERNC EE* SS_edge(SS *thee, int i, int j);

/**
 * @ingroup SS
 * @brief   Return local face number of face shared between two simplices. 
 * @author  Michael Holst
 * @note    Class S: Inlineable methods (ves.c) 
 * @return  local face number of face shared between two simplices. 
 * @param   thee Pointer to a simplex
 * @param   sm   Pointer to a simplex
 */
VEXTERNC int SS_sharedFaceLocNum(SS *thee, SS *sm);

/**
 * @ingroup SS
 * @brief   Return the proper edge type (interior/boundary) of an edge
 *          connecting the two given vertices. 
 * @author  Michael Holst
 * @note    Class S: Inlineable methods (ves.c)\n
 *          ASSUME that the edge DOES NOT exist or does not have a correct
 *          and consistent type, and calculate its type from the face types.
 *          \n
 *          The type we calculate is based only on local info; in 3D it may
 *          NOT be the correct edge type for the edge as far as the global
 *          mesh is concerned.  However, the local information about all
 *          simplices using the edge can be combined to determine the 
 *          correct global type.
 *          \n
 *          We must not be tempted to look for an existing edge and then 
 *          grab its type; the edge may have just been created (with 
 *          interior type by default), and in fact WE are being called as 
 *          part of an attempt to calculate its global type...          
 * @return  the proper edge type (interior/boundary) of an edge connecting the two
 *          given vertices
 * @param   thee Pointer to a simplex
 * @param   v0   Pointer to a vertex
 * @param   v1   Pointer to a vertex
 */
VEXTERNC int SS_localEdgeType(SS *thee, VV *v0, VV *v1);

/**
 * @ingroup SS
 * @brief   Reverse the orientation of a simplex by swapping vertices.    
 * @author  Michael Holst
 * @note    Class S: Inlineable methods (ves.c)
 * @return  None
 * @param   thee Pointer to a vertex
 */
VEXTERNC void SS_reverse(SS *thee);

/**
 * @ingroup FF
 * @brief   The FF constructor.
 * @author  Michael Holst
 * @note    Class S: Inlineable methods (ves.c)
 * @return  Pointer to a new allocated FF class
 * @param   grandParentFace ID number for gip sub-structure
 * @param   sPtr            Pointer to a simplex 
 * @param   color           index for a face
 */
VPUBLIC FF* FF_ctor(int grandParentFace, SS *sPtr, int color);

#endif /* _VES_H_ */

