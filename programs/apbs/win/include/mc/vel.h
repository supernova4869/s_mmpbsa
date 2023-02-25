/**
 *  @file       vel.h
 *  @ingroup    global_mc
 *  @brief      Class Vel: a canonical set element object.
 *  @author     Michael Holst
 *  @note       None 
 *  @version    $Id: vel.h,v 1.21 2010/08/12 05:19:07 fetk Exp $ 
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

#ifndef _VEL_H_
#define _VEL_H_

#include <mc/mc_base.h>

/*
 * ***************************************************************************
 * Class Vel: Parameters and datatypes
 * ***************************************************************************
 */

/** @brief Gip.bits masks (implicit "0/1"s in front that pad out the word) 
 *         none */
#define MASK_00000000000000000000000000000000  000000000000 
/** @brief Gip.bits masks (implicit "0/1"s in front that pad out the word) 
 *         reality */
#define MASK_00000000000000000000000000000011  000000000003 
/** @brief Gip.bits masks (implicit "0/1"s in front that pad out the word) 
 *         dim */
#define MASK_00000000000000000000000000001100  000000000014
/** @brief Gip.bits masks (implicit "0/1"s in front that pad out the word) 
 *         class */
#define MASK_00000000000000000000000000110000  000000000060 
/** @brief Gip.bits masks (implicit "0/1"s in front that pad out the word) 
 *         type */
#define MASK_00000000000000000011111111000000  000000037700
/** @brief Gip.bits masks (implicit "0/1"s in front that pad out the word) 
 *         chart */
#define MASK_11111111111111111100000000000000  007777740000 
/* ------------------------------------------------------------------------ */
/** @brief Gip.bits masks (implicit "0/1"s in front that pad out the word) 
 *         none */
#define MASK_11111111111111111111111111111111 ~000000000000 
/** @brief Gip.bits masks (implicit "0/1"s in front that pad out the word) 
 *         reality */
#define MASK_11111111111111111111111111111100 ~000000000003 
/** @brief Gip.bits masks (implicit "0/1"s in front that pad out the word) 
 *         dim */
#define MASK_11111111111111111111111111110011 ~000000000014 
/** @brief Gip.bits masks (implicit "0/1"s in front that pad out the word) 
 *         class */
#define MASK_11111111111111111111111111001111 ~000000000060 
/** @brief Gip.bits masks (implicit "0/1"s in front that pad out the word) 
 *         type */
#define MASK_11111111111111111100000000111111 ~000000037700 
/** @brief Gip.bits masks (implicit "0/1"s in front that pad out the word) 
 *         chart */
#define MASK_00000000000000000011111111111111 ~007777740000 
/* ------------------------------------------------------------------------ */

/** @brief Sip.tags masks (implicit "0/1"s in front that pad out the word)
 *         refine edge0 */
#define MASK_00000000000000000000000000000111  000000000007 
/** @brief Sip.tags masks (implicit "0/1"s in front that pad out the word)
 *         marked edge1 */
#define MASK_00000000000000000000000000111000  000000000070 
/** @brief Sip.tags masks (implicit "0/1"s in front that pad out the word)
 *         marked edge2 */
#define MASK_00000000000000000000000111000000  000000000700 
/** @brief Sip.tags masks (implicit "0/1"s in front that pad out the word)
 *         marked edge3 */
#define MASK_00000000000000000000111000000000  000000007000 
/** @brief Sip.tags masks (implicit "0/1"s in front that pad out the word)
 *         degen marks */
#define MASK_00000000000000000111000000000000  000000070000
/** @brief Sip.tags masks (implicit "0/1"s in front that pad out the word)
 *         refine count */
#define MASK_00000000000000111000000000000000  000000700000 
/** @brief Sip.tags masks (implicit "0/1"s in front that pad out the word)
 *         refine Q0 */
#define MASK_00000000000001000000000000000000  000001000000
/** @brief Sip.tags masks (implicit "0/1"s in front that pad out the word)
 *         refine Q1 */
#define MASK_00000000000010000000000000000000  000002000000
/* ------------------------------------------------------------------------ */
/** @brief Sip.tags masks (implicit "0/1"s in front that pad out the word)
 *         refine edge0 */
#define MASK_11111111111111111111111111111000 ~000000000007
/** @brief Sip.tags masks (implicit "0/1"s in front that pad out the word)
 *         refine edge1 */
#define MASK_11111111111111111111111111000111 ~000000000070
/** @brief Sip.tags masks (implicit "0/1"s in front that pad out the word)
 *         refine edge2 */
#define MASK_11111111111111111111111000111111 ~000000000700 
/** @brief Sip.tags masks (implicit "0/1"s in front that pad out the word)
 *         refine edge3 */
#define MASK_11111111111111111111000111111111 ~000000007000 
/** @brief Sip.tags masks (implicit "0/1"s in front that pad out the word)
 *         degen marks */
#define MASK_11111111111111111000111111111111 ~000000070000 
/** @brief Sip.tags masks (implicit "0/1"s in front that pad out the word)
 *         refine count */
#define MASK_11111111111111000111111111111111 ~000000700000 
/** @brief Sip.tags masks (implicit "0/1"s in front that pad out the word)
 *         refine Q0 */
#define MASK_11111111111110111111111111111111 ~000001000000 
/** @brief Sip.tags masks (implicit "0/1"s in front that pad out the word)
 *         refine Q1 */
#define MASK_11111111111101111111111111111111 ~000002000000 
/* ------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------ */
/** @brief Sip.faces masks (implicit "0/1"s in front that pad out the word) 
 *         face 0 type */
#define MASK_00000000000000000000000011111111  000000000377 
/** @brief Sip.faces masks (implicit "0/1"s in front that pad out the word) 
 *         face 1 type */
#define MASK_00000000000000001111111100000000  000000177400 
/** @brief Sip.faces masks (implicit "0/1"s in front that pad out the word) 
 *         face 2 type */
#define MASK_00000000111111110000000000000000  000077600000 
/** @brief Sip.faces masks (implicit "0/1"s in front that pad out the word) 
 *         face 3 type */
#define MASK_11111111000000000000000000000000  037700000000 
/* ------------------------------------------------------------------------ */
/** @brief Sip.faces masks (implicit "0/1"s in front that pad out the word) 
 *         face 0 type */
#define MASK_11111111111111111111111100000000 ~000000000377
/** @brief Sip.faces masks (implicit "0/1"s in front that pad out the word) 
 *         face 1 type */
#define MASK_11111111111111110000000011111111 ~000000177400
/** @brief Sip.faces masks (implicit "0/1"s in front that pad out the word) 
 *         face 2 type */
#define MASK_11111111000000001111111111111111 ~000077600000 
/** @brief Sip.faces masks (implicit "0/1"s in front that pad out the word) 
 *         face 3 type */
#define MASK_00000000111111111111111111111111 ~037700000000 
/* ------------------------------------------------------------------------ */

/**
 * @ingroup global_mc
 * @brief   struct Gip (2*4=8 bytes) 
 * @author  Michael Holst
 */
typedef struct Gip {
  /**
   * @brief See note
   * @note
   * @verbatim 
   * (0-1)---> Reality (0--3=[2^2-1])           
   * -------->     0 = real                          
   * -------->     1 = ghost                          
   * -------->     2 = deleted                         
   * -------->     3 = <unused>                        
   * (2-3)---> Intrinsic Dimension (0--3=[2^2-1])     
   * -------->     0 = <unused>                       
   * -------->     1 = 1D                              
   * -------->     2 = 2D                              
   * -------->     3 = 3D                              
   * (4-5)---> Simplex Class (0--3=[2^2-1])            
   * -------->     0 = 0-simplex (vertex)             
   * -------->     1 = 1-simplex (edge)               
   * -------->     2 = 2-simplex (triangle)           
   * -------->     3 = 3-simplex (tetrahedron)        
   * (6-13)--> Material Type (0--256=[2^8-1])         
   * -------->     Simplex interpretation:            
   * -------->       0,1,2,3,...       = Material tag 
   * -------->     Vertex and Edge interpretation:    
   * -------->       0                 = Interior     
   * -------->       ODD (1,3,5,...)   = Dirichlet    
   * -------->       EVENP (2,4,6,...) = Neumann      
   * (14-31)-> Chart (in the atlas) (0--262K=[2^18-1])
   *  @endverbatim */
    unsigned int bits; 
  /** @brief --------> Unit id (0--4G=[2^32-1]) */
    unsigned int uid;
} Gip;

/**
 * @ingroup global_mc
 * @brief   struct Vip (3*8+3*4=36 bytes) 
 * @author  Michael Holst
 */
typedef struct Vip {
  /** @brief --------> My coordinates (wrt my chart)           */
    double x[3];  
  /** @brief --------> Ptr to first of my edges                */  
    void   *ePtr; 
  /** @brief --------> Ptr to first of my simplices            */
    void   *sPtr; 
  /** @brief --------> Ptr to parent edge of which I am midpt  */
    void   *eParent;  
} Vip;

/** 
 * @ingroup global_mc
 * @brief   struct Eip (6*4=24 bytes) 
 * @author  Michael Holst
 */
typedef struct Eip { 
  /** @brief --------> Ptr to vertices forming my edge         */
    void   *vPtr[2]; 
  /** @brief --------> Ptr to next edge sharing my verts       */
    void   *ePtr[2]; 
  /** @brief --------> Ptr to my refined midpoint vertex       */
    void   *midPtr;  
  /** @brief --------> Ptr to parent edge (1 of 2 siblings)    */
    void   *eParent;   
} Eip;

/** 
 * @ingroup global_mc
 * @brief   struct Fip (3*4=12 bytes) 
 * @author  Michael Holst
 */
typedef struct Fip { 
  /** @brief --------> Group ID for simplex with mirror face   */
    int    gID;    
  /** @brief --------> Ptr to my (single) owning simplex       */  
    void   *sPtr;  
  /** @brief -------> child ID containing refinement history   */
    unsigned int rhist;
  /** @brief -------> refinement level for this face           */
    unsigned int rlevel;
} Fip;

/**
 * @ingroup global_mc
 * @brief   (14*4=56 bytes OR 24*4=96 bytes) 
 * @author  Michael Holst
 */
typedef struct Sip {
  /** @brief --------> Ptr to my dim+1 verts                   */
    void   *vPtr[4]; 
  /** @brief --------> dim+1 simplex ring ptrs                 */
    void   *sPtr[4];
  /** @brief --------> dim+1 potential boundary faces          */ 
    void   *fPtr[4];
#if defined(VG_ELEMENT)
  /** @brief --------> Global face numbers                     */
    int    fNum[4];  
  /** @brief --------> Global edge numbers                     */
    int    eNum[6];   
#endif
  /**
   * @brief See note
   * @note
   * @verbatim
   * (0-2)---> Marked edge 0 (ref edge) (0--5<[2^3]-1)
   * (3-5)---> Marked edge 1            (0--5<[2^3]-1)
   * (6-8)---> Marked edge 2            (0--5<[2^3]-1)
   * (9-11)--> Marked edge 3            (0--5<[2^3]-1) 
   * (12-14)-> Degenerate markings      (0--5<[2^3]-1) 
   * (15-17)-> Refinement counter       (0--7=[2^3]-1) 
   * (18-18)-> I'm on refinement Q0     (0=no, 1=yes)  
   * (19-19)-> I'm on refinement Q1     (0=no, 1=yes)  
   * (20-31)-> <unused>                                
   * @endverbatim
   */ 
    unsigned int tags; 
  /** @verbatim
   * (0-7)---> Face 0 (0=Interior,ODD=Diri,EVENP=Neum) 
   * (8-15)--> Face 1 (0=Interior,ODD=Diri,EVENP=Neum)
   * (16-23)-> Face 2 (0=Interior,ODD=Diri,EVENP=Neum)
   * (24-31)-> Face 3 (0=Interior,ODD=Diri,EVENP=Neum)
   * @endverbatim
   */
    unsigned int faces;
} Sip;

/**
 * @ingroup global_mc
 * @brief   typedef union VES 
 *          (40 bytes=max[36,24,40] OR 80 bytes=max[36,24,80])
 * @author  Michael Holst
 */
typedef union VES {
  /** @brief Vertex information package */
    Vip v;  
  /** @brief Edge information package */
    Eip e;
  /** @brief Simplex information package */
    Sip s;
} VES;

/*
 * ***************************************************************************
 * Class Vel: Global objects
 * ***************************************************************************
 */

/** @brief Useful simplex geometry mappings */
VEXTERNC const int vmapE[4][4];
/** @brief Useful simplex geometry mappings */
VEXTERNC const int vmapEI[6][2];
/** @brief Useful simplex geometry mappings */
VEXTERNC const int vmapV[6][4];
/** @brief Useful simplex geometry mappings */
VEXTERNC const int vmapF[4][3];
/** @brief Useful simplex geometry mappings */
VEXTERNC const int vmapFE[4][3];
/** @brief Useful simplex geometry mappings */
VEXTERNC const int vmapOV1[4][4][4];
/** @brief Useful simplex geometry mappings */
VEXTERNC const int vmapOV2[2][4][4];
/** @brief Useful simplex geometry mappings */
VEXTERNC const int vmapOV3[4][3];
/** @brief Useful simplex geometry mappings */
VEXTERNC const int vmapOE[4][3];

/**
 * @ingroup global_mc
 * @brief   Class Vel: Definition (8+80 bytes) 
 * @author  Michael Holst
 */
typedef struct Vel {
    Gip      g;        /**< @brief geometric information package            */
    VES      d;        /**< @brief vertex,edge,simplex data                 */
} Vel;

/*
 * ***************************************************************************
 * Class Vel: Inlineable methods (vel.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_GEM)
    /** 
     *  @ingroup Vel
     *  @brief   Initialize the element. 
     *  @author  Michael Holst
     *  @note    Class Vel: Inlineable methods (vel.c) 
     *  @return  None
     *  @param   thee  Pointer to class Vel
     *  @param   tdim  the dimension bits
     *  @param   tid   the ID bits
     */
    VEXTERNC void Vel_init(Vel *thee, int tdim, int tid);

    /** 
     *  @ingroup Vel
     *  @brief   Re-Initialize the element. 
     *  @author  Michael Holst
     *  @note    Class Vel: Inlineable methods (vel.c) 
     *  @return  None
     *  @param   thee  Pointer to class Vel
     */
    VEXTERNC void Vel_reinit(Vel *thee);

    /** 
     *  @ingroup Vel
     *  @brief   Set the reality bits.
     *  @author  Michael Holst
     *  @note    Class Vel: Inlineable methods (vel.c) 
     *  @return  None
     *  @param   thee  Pointer to class Vel
     *  @param   reel  the reality bits
     */
    VEXTERNC void Vel_setReality(Vel *thee, int reel);

    /** 
     *  @ingroup Vel
     *  @brief   Set the dimension bits.
     *  @author  Michael Holst
     *  @note    Class Vel: Inlineable methods (vel.c) 
     *  @return  None
     *  @param   thee  Pointer to class Vel
     *  @param   dim   the dimension bits
     */
    VEXTERNC void Vel_setDim(Vel *thee, int dim);

    /** 
     *  @ingroup Vel
     *  @brief   Set the class bits.
     *  @author  Michael Holst
     *  @note    Class Vel: Inlineable methods (vel.c) 
     *  @return  None
     *  @param   thee  Pointer to class Vel
     *  @param   clas  the class bits
     */
    VEXTERNC void Vel_setClass(Vel *thee, int clas);

    /** 
     *  @ingroup Vel
     *  @brief   Set the type bits.
     *  @author  Michael Holst
     *  @note    Class Vel: Inlineable methods (vel.c) 
     *  @return  None
     *  @param   thee  Pointer to class Vel
     *  @param   type  the type bits
     */
    VEXTERNC void Vel_setType(Vel *thee, int type);

    /** 
     *  @ingroup Vel
     *  @brief   Set the chart bits.
     *  @author  Michael Holst
     *  @note    Class Vel: Inlineable methods (vel.c) 
     *  @return  None
     *  @param   thee   Pointer to class Vel
     *  @param   chart  the chart bits
     */
    VEXTERNC void Vel_setChart(Vel *thee, int chart);

    /** 
     *  @ingroup Vel
     *  @brief   Set the ID bits.
     *  @author  Michael Holst
     *  @note    Class Vel: Inlineable methods (vel.c) 
     *  @return  None
     *  @param   thee  Pointer to class Vel
     *  @param   id    the ID bits
     */
    VEXTERNC void Vel_setId(Vel *thee, int id);

    /** 
     *  @ingroup Vel
     *  @brief   Return the reality.
     *  @author  Michael Holst
     *  @note    Class Vel: Inlineable methods (vel.c) 
     *  @return  the reality
     *  @param   thee  Pointer to class Vel
     */
    VEXTERNC unsigned int Vel_reality(Vel *thee);

    /** 
     *  @ingroup Vel
     *  @brief   Return the dimension.
     *  @author  Michael Holst
     *  @note    Class Vel: Inlineable methods (vel.c) 
     *  @return  the dimension
     *  @param   thee  Pointer to class Vel
     */
    VEXTERNC unsigned int Vel_dim(Vel *thee);

    /** 
     *  @ingroup Vel
     *  @brief   Return the number of vertices in a simplex.  
     *  @author  Michael Holst
     *  @note    Class Vel: Inlineable methods (vel.c) 
     *  @return  the number of vertices in a simplex.  
     *  @param   thee  Pointer to class Vel
     */
    VEXTERNC unsigned int Vel_dimVV(Vel *thee);

    /** 
     *  @ingroup Vel
     *  @brief   Return the number of edges in a simplex.     
     *  @author  Michael Holst
     *  @note    Class Vel: Inlineable methods (vel.c) 
     *  @return  the number of edges in a simplex.     
     *  @param   thee  Pointer to class Vel
     */
    VEXTERNC unsigned int Vel_dimEE(Vel *thee);

    /** 
     *  @ingroup Vel
     *  @brief   Return the number of faces in a simplex.   
     *  @author  Michael Holst
     *  @note    Class Vel: Inlineable methods (vel.c) 
     *  @return  the number of faces in a simplex.     
     *  @param   thee  Pointer to class Vel
     */
    VEXTERNC unsigned int Vel_dimFF(Vel *thee);

    /** 
     *  @ingroup Vel
     *  @brief   Return the class
     *  @author  Michael Holst
     *  @note    Class Vel: Inlineable methods (vel.c) 
     *  @return  the class
     *  @param   thee  Pointer to class Vel
     */
    VEXTERNC unsigned int Vel_class(Vel *thee);

    /** 
     *  @ingroup Vel
     *  @brief   Return the type
     *  @author  Michael Holst
     *  @note    Class Vel: Inlineable methods (vel.c) 
     *  @return  the type
     *  @param   thee  Pointer to class Vel
     */
    VEXTERNC unsigned int Vel_type(Vel *thee);

    /** 
     *  @ingroup Vel
     *  @brief   Return the chart
     *  @author  Michael Holst
     *  @note    Class Vel: Inlineable methods (vel.c) 
     *  @return  the chart
     *  @param   thee  Pointer to class Vel
     */
    VEXTERNC unsigned int Vel_chart(Vel *thee);

    /** 
     *  @ingroup Vel
     *  @brief   Return the ID
     *  @author  Michael Holst
     *  @note    Class Vel: Inlineable methods (vel.c) 
     *  @return  the ID
     *  @param   thee  Pointer to class Vel
     */
    VEXTERNC unsigned int Vel_id(Vel *thee);
#else /* if defined(VINLINE_GEM) */
/** @brief   Initialize the element
 *  @ingroup Vel
 *  @brief   Initialize the element. 
 *  @author  Michael Holst
 *  @note    Class Vel: Inlineable methods (vel.c) if defined(VINLINE_GEM) 
 *  @return  None
 *  @param   thee  Pointer to class Vel
 *  @param   tdim  the dimension bits
 *  @param   tid   the ID bits
 */
#   define Vel_init(thee,tdim,tid) ( \
        (thee)->g.bits = MASK_00000000000000000000000000000000, \
        (thee)->g.uid  = MASK_00000000000000000000000000000000, \
        Vel_setDim((thee),(tdim)), \
        Vel_setId((thee),(tid)) \
    )

/** 
 *  @ingroup Vel
 *  @brief   Re-Initialize the element. 
 *  @author  Michael Holst
 *  @note    Class Vel: Inlineable methods (vel.c) if defined(VINLINE_GEM) 
 *  @return  None
 *  @param   thee  Pointer to class Vel
 */
#   define Vel_reinit(thee) ( \
        (thee)->g.bits &= MASK_11111111111111111100000000111111, \
        (thee)->g.uid  &= MASK_11111111111111111111111111111111 \
    )

/** 
 *  @ingroup Vel
 *  @brief   Set the reality bits.
 *  @author  Michael Holst
 *  @note    Class Vel: Inlineable methods (vel.c) if defined(VINLINE_GEM) 
 *  @return  None
 *  @param   thee  Pointer to class Vel
 *  @param   reel  the reality bits
 */
#   define Vel_setReality(thee,reel) ( \
        (thee)->g.bits &= MASK_11111111111111111111111111111100, \
        (thee)->g.bits |= (reel) \
    )

/** 
 *  @ingroup Vel
 *  @brief   Set the dimension bits.
 *  @author  Michael Holst
 *  @note    Class Vel: Inlineable methods (vel.c) if defined(VINLINE_GEM) 
 *  @return  None
 *  @param   thee  Pointer to class Vel
 *  @param   dim   the dimension bits
 */
#   define Vel_setDim(thee,dim) ( \
        (thee)->g.bits &= MASK_11111111111111111111111111110011, \
        (thee)->g.bits |= ((dim) << 2) \
    )

/** 
 *  @ingroup Vel
 *  @brief   Set the class bits.
 *  @author  Michael Holst
 *  @note    Class Vel: Inlineable methods (vel.c) if defined(VINLINE_GEM) 
 *  @return  None
 *  @param   thee  Pointer to class Vel
 *  @param   clas  the class bits
 */
#   define Vel_setClass(thee,clas) ( \
        (thee)->g.bits &= MASK_11111111111111111111111111001111, \
        (thee)->g.bits |= ((clas) << 4) \
    )

/** 
 *  @ingroup Vel
 *  @brief   Set the type bits.
 *  @author  Michael Holst
 *  @note    Class Vel: Inlineable methods (vel.c) if defined(VINLINE_GEM) 
 *  @return  None
 *  @param   thee  Pointer to class Vel
 *  @param   type  the type bits
 */
#   define Vel_setType(thee,type) ( \
        (thee)->g.bits &= MASK_11111111111111111100000000111111, \
        (thee)->g.bits |= ((type) << 6) \
    )

/** 
 *  @ingroup Vel
 *  @brief   Set the chart bits.
 *  @author  Michael Holst
 *  @note    Class Vel: Inlineable methods (vel.c) if defined(VINLINE_GEM) 
 *  @return  None
 *  @param   thee   Pointer to class Vel
 *  @param   chart  the chart bits
 */
#   define Vel_setChart(thee,chart) ( \
        (thee)->g.bits &= MASK_00000000000000000011111111111111, \
        (thee)->g.bits |= ((chart) << 14) \
    )

/** 
 *  @ingroup Vel
 *  @brief   Set the ID bits.
 *  @author  Michael Holst
 *  @note    Class Vel: Inlineable methods (vel.c) if defined(VINLINE_GEM) 
 *  @return  None
 *  @param   thee  Pointer to class Vel
 *  @param   id    the ID bits
 */
#   define Vel_setId(thee,id) ((thee)->g.uid = (id))

/** 
 *  @ingroup Vel
 *  @brief   Return the reality.
 *  @author  Michael Holst
 *  @note    Class Vel: Inlineable methods (vel.c) if defined(VINLINE_GEM) 
 *  @return  the reality
 *  @param   thee  Pointer to class Vel
 */
#   define Vel_reality(thee) ( \
        ((thee)->g.bits & MASK_00000000000000000000000000000011) \
    )

/** 
 *  @ingroup Vel
 *  @brief   Return the dimension.
 *  @author  Michael Holst
 *  @note    Class Vel: Inlineable methods (vel.c) if defined(VINLINE_GEM) 
 *  @return  the dimension
 *  @param   thee  Pointer to class Vel
 */
#   define Vel_dim(thee) ( \
        ((thee)->g.bits & MASK_00000000000000000000000000001100) >> 2 \
    )

/** 
 *  @ingroup Vel
 *  @brief   Return the number of vertices in a simplex.  
 *  @author  Michael Holst
 *  @note    Class Vel: Inlineable methods (vel.c) if defined(VINLINE_GEM) 
 *  @return  the number of vertices in a simplex.  
 *  @param   thee  Pointer to class Vel
 */
#   define Vel_dimVV(thee) ( \
        Vel_dim((thee))+1 \
    )

/** 
 *  @ingroup Vel
 *  @brief   Return the number of edges in a simplex.     
 *  @author  Michael Holst
 *  @note    Class Vel: Inlineable methods (vel.c) if defined(VINLINE_GEM) 
 *  @return  the number of edges in a simplex.     
 *  @param   thee  Pointer to class Vel
 */
#   define Vel_dimEE(thee) ( \
        3*(Vel_dim((thee))-1) \
    )

/** 
 *  @ingroup Vel
 *  @brief   Return the number of faces in a simplex.   
 *  @author  Michael Holst
 *  @note    Class Vel: Inlineable methods (vel.c) if defined(VINLINE_GEM) 
 *  @return  the number of faces in a simplex.     
 *  @param   thee  Pointer to class Vel
 */
#   define Vel_dimFF(thee) ( \
        Vel_dim((thee))+1 \
    )

/** 
 *  @ingroup Vel
 *  @brief   Return the class
 *  @author  Michael Holst
 *  @note    Class Vel: Inlineable methods (vel.c) if defined(VINLINE_GEM) 
 *  @return  the class
 *  @param   thee  Pointer to class Vel
 */
#   define Vel_class(thee) ( \
        ((thee)->g.bits & MASK_00000000000000000000000000110000) >> 4 \
    )

/** 
 *  @ingroup Vel
 *  @brief   Return the type
 *  @author  Michael Holst
 *  @note    Class Vel: Inlineable methods (vel.c) if defined(VINLINE_GEM) 
 *  @return  the type
 *  @param   thee  Pointer to class Vel
 */
#   define Vel_type(thee) ( \
        ((thee)->g.bits & MASK_00000000000000000011111111000000) >> 6 \
    )

/** 
 *  @ingroup Vel
 *  @brief   Return the chart
 *  @author  Michael Holst
 *  @note    Class Vel: Inlineable methods (vel.c) if defined(VINLINE_GEM) 
 *  @return  the chart
 *  @param   thee  Pointer to class Vel
 */
#   define Vel_chart(thee) ( \
        ((thee)->g.bits & MASK_11111111111111111100000000000000) >> 14 \
    )

/** 
 *  @ingroup Vel
 *  @brief   Return the ID
 *  @author  Michael Holst
 *  @note    Class Vel: Inlineable methods (vel.c) if defined(VINLINE_GEM) 
 *  @return  the ID
 *  @param   thee  Pointer to class Vel
 */
#   define Vel_id(thee) ((thee)->g.uid)
#endif /* if !defined(VINLINE_GEM) */


/**
 * @ingroup Vel
 * @brief   The Vel constructor.
 * @author  Michael Holst 
 * @note    Class Vel: Non-Inlineable methods (vel.c) 
 * @return  Pointer to a new allocated Vel class
 * @param   dim  the dimension bits
 * @param   id   the ID bits
 */
VEXTERNC Vel* Vel_ctor(int dim, int id);

/**
 * @ingroup Vel
 * @brief   The Vel destructor.
 * @author  Michael Holst 
 * @note    Class Vel: Non-Inlineable methods (vel.c) 
 * @return  None
 * @param   thee  Pointer to class Vel
 */
VEXTERNC void Vel_dtor(Vel **thee);

#endif /* _VEL_H_ */
