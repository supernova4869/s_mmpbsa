/**
 * @defgroup Gem Gem class
 * @brief    The geometry manager object.
 */

/**
 *  @file       gem.h
 *  @ingroup    Gem
 *  @brief      Class Gem: the geometry manager object.
 *  @author     Michael Holst
 *  @note       None
 *  @version    $Id: gem.h,v 1.43 2010/08/12 05:19:05 fetk Exp $ 
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


#ifndef _GEM_H_
#define _GEM_H_

#include <mc/mc_base.h>

#include <mc/pde.h>

#include <mc/vel.h>
#include <mc/ves.h>

#include <mc/bam.h>

/** @brief Class Gem: Parameters and datatypes */

#define VMAXSQ 2


/**
 * @ingroup Gem
 * @brief   Contains public data memebers for Gem class
 * @author  Michael Holst
 */

struct sGem {

  /** @brief Intrinsic spatial dim (2 for sphere)         */
    int     dim;     
  /** @brief Imbedded spatial dim (3 for sphere)          */       
    int     dimII;  
  /** @brief Number of vertices in a dim-simplex          */ 
    int     dimVV; 
  /** @brief Number of edges in a dim-simplex             */ 
    int     dimEE; 

  /** @brief Initial count of the number of vertices      */
    int     numVV0; 
  /** @brief Last count of the number of vertices         */
    int     numVV;  
  /** @brief Last count of the number of edges            */ 
    int     numEE; 
  /** @brief Last count of the number of faces            */ 
    int     numFF; 
  /** @brief Last count of the number of simplices        */
    int     numSS; 
  /** @brief Last count of boundary vertices              */
    int     numBV;  
  /** @brief Last count of boundary faces                 */
    int     numBF;  

  /** @brief the memory manager                           */
    Vmem     *vmem;   
  /** @brief did i make vmem or was it inherited          */
    int      iMadeVmem;    

  /** @brief the set of vertices                          */
    Vset     *vertices;
  /** @brief the set of edges                             */ 
    Vset     *edges;   
  /** @brief the set of simplices                         */
    Vset     *simplices; 

  /** @brief refinement/conformity/flipping simplex Qs    */
    Vset     *sQueM[VMAXSQ];

  /** @brief did I have to make a fake PDE object?        */
    int      iMadePDE; 
  /** @brief container for various user-provided functions*/
    PDE      *pde;   

  /** @brief Hook for external structure updating */
    int xUpFlag;
  /** @brief Hook for external structure updating */
    void (*xUp)(SS **sms, int numS);

};

/**
 * @brief   Declaration of the Gem class as the Gem structure
 * @ingroup Gem
 * @author  Michael Holst
 * @return  None
 */
typedef struct sGem Gem;

/*
 * ***************************************************************************
 * Class Gem: Inlineable methods (gem.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_GEM)
    /**
     * @ingroup Gem
     * @brief   Return the extrinsic spatial dimension.   
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the extrinsic spatial dimension.   
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC int Gem_dim(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return the extrinsic spatial dimension.   
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the extrinsic spatial dimension.   
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC int Gem_dimII(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return the number of vertices in a simplex.
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the number of vertices in a simplex
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC int Gem_dimVV(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return the number of edges in a simplex.
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the number of edges in a simplex
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC int Gem_dimEE(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return the logical number of vertices in the mesh.     
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the logical number of vertices in the mesh.     
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC int Gem_numVirtVV(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return the logical number of edges in the mesh. 
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the logical number of edges in the mesh. 
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC int Gem_numVirtEE(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return the logical number of faces in the mesh. 
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the logical number of faces in the mesh. 
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC int Gem_numVirtFF(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return the logical number of simplices in the mesh. 
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the logical number of simplices in the mesh. 
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC int Gem_numVirtSS(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Set the logical number of vertices in the mesh.     
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  None
     * @param   thee  Pointer to class Gem
     * @param   i     index for the logical number of vertices
     */
    VEXTERNC void Gem_setNumVirtVV(Gem *thee, int i);

    /**
     * @ingroup Gem
     * @brief   Set the logical number of edges in the mesh.     
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  None
     * @param   thee  Pointer to class Gem
     * @param   i     index for the logical number of edges
     */
    VEXTERNC void Gem_setNumVirtEE(Gem *thee, int i);

    /**
     * @ingroup Gem
     * @brief   Geometry manager constructor. 
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  None
     * @param   thee  Pointer to class Gem
     * @param   i     index for the logical number of faces
     */
    VEXTERNC void Gem_setNumVirtFF(Gem *thee, int i);

    /**
     * @ingroup Gem
     * @brief   Set the logical number of simplices in the mesh. 
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  None
     * @param   thee  Pointer to class Gem
     * @param   i     index for the logical number of simplices
     */
    VEXTERNC void Gem_setNumVirtSS(Gem *thee, int i);

    /**
     * @ingroup Gem
     * @brief   Return the number of vertices in the mesh.        
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the number of vertices in the mesh.        
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC int Gem_numVV(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return a given vertex.  
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  a given vertex.
     * @param   thee  Pointer to class Gem
     * @param   i     index for a given vertex
     */
    VEXTERNC VV* Gem_VV(Gem *thee, int i);

    /**
     * @ingroup Gem
     * @brief   Create a new vertex (becoming the new last vertex).   
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  Pointer to the new created vertex
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC VV* Gem_createVV(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return the first vertex.   
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the first vertex
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC VV* Gem_firstVV(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return the last vertex.   
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the last vertex
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC VV* Gem_lastVV(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return the next vertex.   
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the next vertex
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC VV* Gem_nextVV(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return the previous vertex.   
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the previous vertex
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC VV* Gem_prevVV(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Peek at the first vertex.   
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  Pointer to the first vertex in the list
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC VV* Gem_peekFirstVV(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Peek at the last vertex.   
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  Pointer to the last vertex in the list
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC VV* Gem_peekLastVV(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Destroy the last vertex.          
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  None
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC void Gem_destroyVV(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Destroy all of the vertices.     
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  None
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC void Gem_resetVV(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return the number of edge.
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the number of edges.
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC int Gem_numEE(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return a given edge.
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  a given edge
     * @param   thee  Pointer to class Gem
     * @param   i     index for a given edge
     */
    VEXTERNC EE* Gem_EE(Gem *thee, int i);

    /**
     * @ingroup Gem
     * @brief   Create a new edge (becoming the new last edge).
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  Pointer to the new created edge
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC EE* Gem_createEE(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return the first edge.   
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the first edge
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC EE* Gem_firstEE(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return the last edge.   
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the last edge
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC EE* Gem_lastEE(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return the next edge.   
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the next edge
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC EE* Gem_nextEE(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return the previous edge.   
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the previous edge
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC EE* Gem_prevEE(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Peek at the first edge. 
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  Pointer to the first edge in the list
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC EE* Gem_peekFirstEE(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Peek at the last edge. 
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  Pointer to the last edge in the list
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC EE* Gem_peekLastEE(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Destroy the last edge.                                
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  None
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC void Gem_destroyEE(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Destroy all of the edges.
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  None
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC void Gem_resetEE(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return the number of simplices.    
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the number of simplicies
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC int Gem_numSS(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return a given simplex.
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  a given simplex
     * @param   thee  Pointer to class Gem
     * @param   i     Pointer to a given simplex
     */
    VEXTERNC SS* Gem_SS(Gem *thee, int i);

    /**
     * @ingroup Gem
     * @brief   Create a new simplex (becoming the new last simplex). 
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  a new simplex (becoming the new last simplex). 
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC SS* Gem_createSS(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return the first simplex.
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the first simplex
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC SS* Gem_firstSS(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return the last simplex.
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the last simplex
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC SS* Gem_lastSS(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return the next simplex.
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the next simplex
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC SS* Gem_nextSS(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return the previous simplex.
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the previous simplex
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC SS* Gem_prevSS(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Peek at the first simplex.   
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  Pointer to the first simplex in the list
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC SS* Gem_peekFirstSS(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Peek at the last simplex.   
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  Pointer to the last simplex in the list
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC SS* Gem_peekLastSS(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Destroy the last simplex.  
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  None
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC void Gem_destroySS(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Destroy all of the simplices.      
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  None
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC void Gem_resetSS(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return the number of simplices in a given queue.  
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the number of simplices in a given queue.  
     * @param   thee      Pointer to class Gem
     * @param   currentQ  index of a given queue
     */
    VEXTERNC int Gem_numSQ(Gem *thee, int currentQ);

    /**
     * @ingroup Gem
     * @brief   Release all of the simplices in a given queue.   
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  None
     * @param   thee      Pointer to class Gem
     * @param   currentQ  index of a given queue
     */
    VEXTERNC void Gem_resetSQ(Gem *thee, int currentQ);

    /**
     * @ingroup Gem
     * @brief   Return the number of boundary faces.        
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the number of boundary faces.        
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC int Gem_numBF(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Return the number of boundary vertices.    
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  the number of boundary vertices
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC int Gem_numBV(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Set the number of boundary faces.        
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  None
     * @param   thee  Pointer to class Gem
     * @param   val   value for boundary faces
     */
    VEXTERNC void Gem_setNumBF(Gem *thee, int val);

    /**
     * @ingroup Gem
     * @brief   Set the number of boundary vertices.
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  None
     * @param   thee  Pointer to class Gem
     * @param   val   value for boundary vertices
     */
    VEXTERNC void Gem_setNumBV(Gem *thee, int val);

    /**
     * @ingroup Gem
     * @brief   Increment the number of boundary faces by a given integer.
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  None
     * @param   thee  Pointer to class Gem
     * @param   val   Value for incrementing the number of boundary faces
     */
    VEXTERNC void Gem_addToNumBF(Gem *thee, int val);

    /**
     * @ingroup Gem
     * @brief   Increment the number of boundary vertices by a given integer. 
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  None
     * @param   thee  Pointer to class Gem
     * @param   val   value for incrementing the number of boundary vertices
     */
    VEXTERNC void Gem_addToNumBV(Gem *thee, int val);

    /**
     * @ingroup Gem
     * @brief   Increment the number of boundary faces by 1.     
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  None
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC void Gem_numBFpp(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Increment the number of boundary vertices by 1.
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  None
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC void Gem_numBVpp(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Decrement the number of boundary faces by 1.    
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  None
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC void Gem_numBFmm(Gem *thee);

    /**
     * @ingroup Gem
     * @brief   Decrement the number of boundary vertices by 1. 
     * @author  Michael Holst
     * @note    Class Gem: Inlineable methods (gem.c) if !defined(VINLINE_GEM) 
     * @return  None
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC void Gem_numBVmm(Gem *thee);
#else /* if defined(VINLINE_GEM) */

/**
 * @ingroup Gem
 * @brief   Return the extrinsic spatial dimension.   
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the extrinsic spatial dimension.   
 * @param   thee  Pointer to class Gem
 */
#   define Gem_dim(thee) ((thee)->dim)

/**
 * @ingroup Gem
 * @brief   Return the extrinsic spatial dimension.   
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the extrinsic spatial dimension.   
 * @param   thee  Pointer to class Gem
 */
#   define Gem_dimII(thee) ((thee)->dimII)

/**
 * @ingroup Gem
 * @brief   Return the number of vertices in a simplex.
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the number of vertices in a simplex
 * @param   thee  Pointer to class Gem
 */
#   define Gem_dimVV(thee) ((thee)->dimVV)

/**
 * @ingroup Gem
 * @brief   Return the number of edges in a simplex.
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the number of edges in a simplex
 * @param   thee  Pointer to class Gem
 */
#   define Gem_dimEE(thee) ((thee)->dimEE)

/**
 * @ingroup Gem
 * @brief   Return the logical number of vertices in the mesh.     
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the logical number of vertices in the mesh.     
 * @param   thee  Pointer to class Gem
 */
#   define Gem_numVirtVV(thee) ((thee)->numVV)

/**
 * @ingroup Gem
 * @brief   Return the logical number of edges in the mesh. 
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the logical number of edges in the mesh. 
 * @param   thee  Pointer to class Gem
 */
#   define Gem_numVirtEE(thee) ((thee)->numEE)

/**
 * @ingroup Gem
 * @brief   Return the logical number of faces in the mesh. 
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the logical number of faces in the mesh. 
 * @param   thee  Pointer to class Gem
 */
#   define Gem_numVirtFF(thee) ((thee)->numFF)

/**
 * @ingroup Gem
 * @brief   Return the logical number of simplices in the mesh. 
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the logical number of simplices in the mesh. 
 * @param   thee  Pointer to class Gem
 */
#   define Gem_numVirtSS(thee) ((thee)->numSS)

/**
 * @ingroup Gem
 * @brief   Set the logical number of vertices in the mesh.     
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  None
 * @param   thee  Pointer to class Gem
 * @param   i     index for the logical number of vertices
 */
#   define Gem_setNumVirtVV(thee,i) ((thee)->numVV = (i))

/**
 * @ingroup Gem
 * @brief   Set the logical number of edges in the mesh.     
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  None
 * @param   thee  Pointer to class Gem
 * @param   i     index for the logical number of edges
 */
#   define Gem_setNumVirtEE(thee,i) ((thee)->numEE = (i))

/**
 * @ingroup Gem
 * @brief   Geometry manager constructor. 
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  None
 * @param   thee  Pointer to class Gem
 * @param   i     index for the logical number of faces
 */
#   define Gem_setNumVirtFF(thee,i) ((thee)->numFF = (i))

/**
 * @ingroup Gem
 * @brief   Set the logical number of simplices in the mesh. 
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  None
 * @param   thee  Pointer to class Gem
 * @param   i     index for the logical number of simplices
 */
#   define Gem_setNumVirtSS(thee,i) ((thee)->numSS = (i))

/**
 * @ingroup Gem
 * @brief   Return the number of vertices in the mesh.        
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the number of vertices in the mesh.        
 * @param   thee  Pointer to class Gem
 */
#   define Gem_numVV(thee) (Vset_num((thee)->vertices))

/**
 * @ingroup Gem
 * @brief   Return a given vertex.  
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  a given vertex.
 * @param   thee  Pointer to class Gem
 * @param   i     index for a given vertex
 */
#   define Gem_VV(thee,i) ((VV*)Vset_access((thee)->vertices,(i)))

/**
 * @ingroup Gem
 * @brief   Create a new vertex (becoming the new last vertex).   
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  Pointer to the new created vertex
 * @param   thee  Pointer to class Gem
 */
#   define Gem_createVV(thee) ((VV*)Vset_create((thee)->vertices))

/**
 * @ingroup Gem
 * @brief   Return the first vertex.   
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the first vertex
 * @param   thee  Pointer to class Gem
 */
#   define Gem_firstVV(thee) ((VV*)Vset_first((thee)->vertices))

/**
 * @ingroup Gem
 * @brief   Return the last vertex.   
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the last vertex
 * @param   thee  Pointer to class Gem
 */
#   define Gem_lastVV(thee) ((VV*)Vset_last((thee)->vertices))

/**
 * @ingroup Gem
 * @brief   Return the next vertex.   
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the next vertex
 * @param   thee  Pointer to class Gem
 */
#   define Gem_nextVV(thee) ((VV*)Vset_next((thee)->vertices))

/**
 * @ingroup Gem
 * @brief   Return the previous vertex.   
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the previous vertex
 * @param   thee  Pointer to class Gem
 */
#   define Gem_prevVV(thee) ((VV*)Vset_prev((thee)->vertices))

/**
 * @ingroup Gem
 * @brief   Peek at the first vertex.   
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  Pointer to the first vertex in the list
 * @param   thee  Pointer to class Gem
 */
#   define Gem_peekFirstVV(thee) ((VV*)Vset_peekFirst((thee)->vertices))

/**
 * @ingroup Gem
 * @brief   Peek at the last vertex.   
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  Pointer to the last vertex in the list
 * @param   thee  Pointer to class Gem
 */
#   define Gem_peekLastVV(thee) ((VV*)Vset_peekLast((thee)->vertices))

/**
 * @ingroup Gem
 * @brief   Destroy the last vertex.          
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  None
 * @param   thee  Pointer to class Gem
 */
#   define Gem_destroyVV(thee) (Vset_destroy((thee)->vertices))

/**
 * @ingroup Gem
 * @brief   Destroy all of the vertices.     
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  None
 * @param   thee  Pointer to class Gem
 */
#   define Gem_resetVV(thee) (Vset_reset((thee)->vertices))

/**
 * @ingroup Gem
 * @brief   Return the number of edge.
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the number of edges.
 * @param   thee  Pointer to class Gem
 */
#   define Gem_numEE(thee) (Vset_num((thee)->edges))

/**
 * @ingroup Gem
 * @brief   Return a given edge.
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  a given edge
 * @param   thee  Pointer to class Gem
 * @param   i     index for a given edge
 */
#   define Gem_EE(thee,i) ((EE*)Vset_access((thee)->edges,(i)))

/**
 * @ingroup Gem
 * @brief   Create a new edge (becoming the new last edge).
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  Pointer to the new created edge
 * @param   thee  Pointer to class Gem
 */
#   define Gem_createEE(thee) ((EE*)Vset_create((thee)->edges))

/**
 * @ingroup Gem
 * @brief   Return the first edge.   
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the first edge
 * @param   thee  Pointer to class Gem
 */
#   define Gem_firstEE(thee) ((EE*)Vset_first((thee)->edges))

/**
 * @ingroup Gem
 * @brief   Return the last edge.   
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the last edge
 * @param   thee  Pointer to class Gem
 */
#   define Gem_lastEE(thee) ((EE*)Vset_last((thee)->edges))

/**
 * @ingroup Gem
 * @brief   Return the next edge.   
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the next edge
 * @param   thee  Pointer to class Gem
 */
#   define Gem_nextEE(thee) ((EE*)Vset_next((thee)->edges))

/**
 * @ingroup Gem
 * @brief   Return the previous edge.   
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the previous edge
 * @param   thee  Pointer to class Gem
 */
#   define Gem_prevEE(thee) ((EE*)Vset_prev((thee)->edges))

/**
 * @ingroup Gem
 * @brief   Peek at the first edge. 
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  Pointer to the first edge in the list
 * @param   thee  Pointer to class Gem
 */
#   define Gem_peekFirstEE(thee) ((EE*)Vset_peekFirst((thee)->edges))

/**
 * @ingroup Gem
 * @brief   Peek at the last edge. 
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  Pointer to the last edge in the list
 * @param   thee  Pointer to class Gem
 */
#   define Gem_peekLastEE(thee) ((EE*)Vset_peekLast((thee)->edges))

/**
 * @ingroup Gem
 * @brief   Destroy the last edge.                                
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  None
 * @param   thee  Pointer to class Gem
 */
#   define Gem_destroyEE(thee) (Vset_destroy((thee)->edges))

/**
 * @ingroup Gem
 * @brief   Destroy all of the edges.
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  None
 * @param   thee  Pointer to class Gem
 */
#   define Gem_resetEE(thee) (Vset_reset((thee)->edges))

/**
 * @ingroup Gem
 * @brief   Return the number of simplices.    
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the number of simplicies
 * @param   thee  Pointer to class Gem
 */
#   define Gem_numSS(thee) (Vset_num((thee)->simplices))

/**
 * @ingroup Gem
 * @brief   Return a given simplex.
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  a given simplex
 * @param   thee  Pointer to class Gem
 * @param   i     Pointer to a given simplex
 */
#   define Gem_SS(thee,i) ((SS*)Vset_access((thee)->simplices,(i)))

/**
 * @ingroup Gem
 * @brief   Create a new simplex (becoming the new last simplex). 
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  a new simplex (becoming the new last simplex). 
 * @param   thee  Pointer to class Gem
 */
#   define Gem_createSS(thee) ((SS*)Vset_create((thee)->simplices))

/**
 * @ingroup Gem
 * @brief   Return the first simplex.
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the first simplex
 * @param   thee  Pointer to class Gem
 */
#   define Gem_firstSS(thee) ((SS*)Vset_first((thee)->simplices))

/**
 * @ingroup Gem
 * @brief   Return the last simplex.
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the last simplex
 * @param   thee  Pointer to class Gem
 */
#   define Gem_lastSS(thee) ((SS*)Vset_last((thee)->simplices))

/**
 * @ingroup Gem
 * @brief   Return the next simplex.
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the next simplex
 * @param   thee  Pointer to class Gem
 */
#   define Gem_nextSS(thee) ((SS*)Vset_next((thee)->simplices))

/**
 * @ingroup Gem
 * @brief   Return the previous simplex.
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the previous simplex
 * @param   thee  Pointer to class Gem
 */
#   define Gem_prevSS(thee) ((SS*)Vset_prev((thee)->simplices))

/**
 * @ingroup Gem
 * @brief   Peek at the first simplex.   
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  Pointer to the first simplex in the list
 * @param   thee  Pointer to class Gem
 */
#   define Gem_peekFirstSS(thee) ((SS*)Vset_peekFirst((thee)->simplices))

/**
 * @ingroup Gem
 * @brief   Peek at the last simplex.   
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  Pointer to the last simplex in the list
 * @param   thee  Pointer to class Gem
 */
#   define Gem_peekLastSS(thee) ((SS*)Vset_peekLast((thee)->simplices))

/**
 * @ingroup Gem
 * @brief   Destroy the last simplex.  
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  None
 * @param   thee  Pointer to class Gem
 */
#   define Gem_destroySS(thee) (Vset_destroy((thee)->simplices))

/**
 * @ingroup Gem
 * @brief   Destroy all of the simplices.      
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  None
 * @param   thee  Pointer to class Gem
 */
#   define Gem_resetSS(thee) (Vset_reset((thee)->simplices))

/**
 * @ingroup Gem
 * @brief   Return the number of simplices in a given queue.  
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the number of simplices in a given queue.  
 * @param   thee      Pointer to class Gem
 * @param   currentQ  index of a given queue
 */
#   define Gem_numSQ(thee,currentQ) (Vset_num((thee)->sQueM[(currentQ)]))

/**
 * @ingroup Gem
 * @brief   Release all of the simplices in a given queue.   
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  None
 * @param   thee      Pointer to class Gem
 * @param   currentQ  index of a given queue
 */
#   define Gem_resetSQ(thee,currentQ) (Vset_reset((thee)->sQueM[(currentQ)]))

/**
 * @ingroup Gem
 * @brief   Return the number of boundary faces.        
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the number of boundary faces.        
 * @param   thee  Pointer to class Gem
 */
#   define Gem_numBF(thee)  ((thee)->numBF)

/**
 * @ingroup Gem
 * @brief   Return the number of boundary vertices.    
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  the number of boundary vertices
 * @param   thee  Pointer to class Gem
 */
#   define Gem_numBV(thee)  ((thee)->numBV)

/**
 * @ingroup Gem
 * @brief   Set the number of boundary faces.        
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  None
 * @param   thee  Pointer to class Gem
 * @param   val   value for boundary faces
 */
#   define Gem_setNumBF(thee,val) ((thee)->numBF = (val))

/**
 * @ingroup Gem
 * @brief   Set the number of boundary vertices.
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  None
 * @param   thee  Pointer to class Gem
 * @param   val   value for boundary vertices
 */
#   define Gem_setNumBV(thee,val) ((thee)->numBV = (val))

/**
 * @ingroup Gem
 * @brief   Increment the number of boundary faces by a given integer.
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  None
 * @param   thee  Pointer to class Gem
 * @param   val   Value for incrementing the number of boundary faces
 */
#   define Gem_addToNumBF(thee,val) ((thee)->numBF += (val))

/**
 * @ingroup Gem
 * @brief   Increment the number of boundary vertices by a given integer. 
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  None
 * @param   thee  Pointer to class Gem
 * @param   val   value for incrementing the number of boundary vertices
 */
#   define Gem_addToNumBV(thee,val) ((thee)->numBV += (val))

/**
 * @ingroup Gem
 * @brief   Increment the number of boundary faces by 1.     
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  None
 * @param   thee  Pointer to class Gem
 */
#   define Gem_numBFpp(thee) ((thee)->numBF++)

/**
 * @ingroup Gem
 * @brief   Increment the number of boundary vertices by 1.
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  None
 * @param   thee  Pointer to class Gem
 */
#   define Gem_numBVpp(thee) ((thee)->numBV++)

/**
 * @ingroup Gem
 * @brief   Decrement the number of boundary faces by 1.    
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  None
 * @param   thee  Pointer to class Gem
 */
#   define Gem_numBFmm(thee) ((thee)->numBF--)

/**
 * @ingroup Gem
 * @brief   Decrement the number of boundary vertices by 1. 
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c) if defined(VINLINE_GEM) 
 * @return  None
 * @param   thee  Pointer to class Gem
 */
#   define Gem_numBVmm(thee) ((thee)->numBV--)
#endif /* if !defined(VINLINE_GEM) */


/**
 * @ingroup Gem
 * @brief   Geometry manager constructor.
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c)
 * @return  Pointer to a newly allocated (empty) Gem class
 * @param   vmem  Memory management object
 * @param   tpde  Pointer to the PDE object
 */
VEXTERNC Gem* Gem_ctor(Vmem *vmem, PDE *tpde);

/**
 * @ingroup Gem
 * @brief   Geometry manager destructor.
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c)
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_dtor(Gem **thee);

/**
 * @ingroup Gem
 * @brief   Reset all of the geometry datastructures.    
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c)
 * @return  None
 * @param   thee   Pointer to class Gem
 * @param   dim    the extrinsic spatial dimension
 * @param   dimII  the intrinsic spatial dimension
 */
VEXTERNC void Gem_reset(Gem *thee, int dim, int dimII);

/**
 * @ingroup Gem
 * @brief   Create and initialize a new vertex; return a point to it.   
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c)
 * @return  a point to the newly created and initialized vertex
 * @param   thee   Pointer to class Gem
 */
VEXTERNC VV* Gem_createAndInitVV(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Create and initialize a new edge; return a point to it. 
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c)
 * @return  a point to the newly created and initialized edge
 * @param   thee   Pointer to class Gem
 */
VEXTERNC EE* Gem_createAndInitEE(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Create and initialize a new simplex; return a point to it. 
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c)
 * @return  a point to the newly created and initialized simplex
 * @param   thee   Pointer to class Gem
 */
VEXTERNC SS* Gem_createAndInitSS(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Return the simplex at a particular location in the simplex Q.   
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c)
 * @return  the simplex at a particular location in the simplex Q.   
 * @param   thee      Pointer to class Gem
 * @param   currentQ  index of a given queue
 * @param   i         index for the number of simplices in a given queue
 */
VEXTERNC SS* Gem_SQ(Gem *thee, int currentQ, int i);

/**
 * @ingroup Gem
 * @brief   Append a simplex to the end of a simplex Q.
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c)
 * @return  None
 * @param   thee      Pointer to class Gem
 * @param   currentQ  index of a given queue
 * @param   qsm       Pointer to the simplex
 */
VEXTERNC void Gem_appendSQ(Gem *thee, int currentQ, SS *qsm);

/**
 * @ingroup Gem
 * @brief   Create all of the simplex rings.
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c)
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_createSimplexRings(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Destroy all of the simplex rings.    
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c)
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_destroySimplexRings(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Look for a common edges between two vertices.  
 *          If it doesn't yet exist, we create it, and then note that
 *          we did so.
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c)
 * @return  Pointer to the common edges between two vertices.
 * @param   thee  Pointer to class Gem
 * @param   v0    Pointer to the first vertex
 * @param   v1    Pointer to the second vertex
 * @param   iDid  index for creating an edge
 */
VEXTERNC EE* Gem_findOrCreateEdge(Gem *thee, VV *v0, VV *v1, int *iDid);

/**
 * @ingroup Gem
 * @brief   Create all of the edges. 
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c)\n
 *          Based on a simplex traversal. 
 *          We also set the edge numbers in the simplices while we we 
 *          are doing the edge creation.
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_createEdges(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Destroy all of the edges. 
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c)
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_destroyEdges(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Count all vertices, edges, faces, simplices, and do it   
 *          is cheaply as possible.  Also do a form check. 
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c)
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_countChk(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Count all of the faces without actually creating them. 
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c)\n
 *          Keep track of the global face numbers in each element.   
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_countFaces(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Clear all of the edge numbers in each simplex.    
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c)
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_clearEdges(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Count up all of the edges without actually creating them, and  
 *          keep track of the global edge numbers in each element. 
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c)\n
 *          Based on a simplex traversal.    
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_countEdges(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Make some specified hacked fix to a given mesh. 
 * @author  Michael Holst
 * @note    Class Gem: Inlineable methods (gem.c)
 * @return  None
 * @param   thee  Pointer to class Gem
 * @param   key   0 --> ?
 */
VEXTERNC void Gem_formFix(Gem *thee, int key);

#if 0
    /** 
     * @ingroup Gem
     * @brief   Hook for external structure updating 
     * @author  Michael Holst
     * @note    We intentionally do not define these three prototypes.    
     *          These three routines may be present in the library depending 
     *          on how it was compiled.  Users of these three functions must
     *          provide their own prototypes if they have built the library
     *          to enable them.
     * @return  index for hooking for external structure updating
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC int Gem_externalUpdateFlag(Gem *thee);
    /** 
     * @ingroup Gem
     * @brief   Hook for external structure updating 
     * @author  Michael Holst
     * @note    We intentionally do not define these three prototypes.   
     *          These three routines may be present in the library depending 
     *          on how it was compiled.  Users of these three functions must 
     *          provide their own prototypes if they have built the library 
     *          to enable them.      
     * @return  None
     * @param   thee  Pointer to class Gem
     * @param   fl    index for hooking for external structure updating
     */
    VEXTERNC void Gem_setExternalUpdateFlag(Gem *thee, int fl);
    /** 
     * @ingroup Gem
     * @brief   Hook for external structure updating 
     * @author  Michael Holst
     * @note    We intentionally do not define these three prototypes.   
     *          These three routines may be present in the library depending 
     *          on how it was compiled.  Users of these three functions must 
     *          provide their own prototypes if they have built the library 
     *          to enable them.      
     * @return  None
     * @param   thee  Pointer to class Gem
     */
    VEXTERNC void Gem_setExternalUpdateFunction(Gem *thee,
        void (*xUp)(SS **sms, int numS));
#endif


/**
 * @ingroup Gem
 * @brief   Build the basic master-to-element transformation information.
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemg.c) \n\n
 *          gchart             ==> unified common chart for vertex coordinates\n
 *          chart[4]           ==> individual charts for vertex coordinates\n
 *          \n
 *          D                  ==> jacobian determinant of the transformation\n
 *          Dcook              ==> jacobian determinant of cooked trans\n
 *          faceD[4]           ==> face jacobian determinants\n
 *          \n
 *          ff[3][3], bb[3]    ==> affine trans from master to arbitrary el\n
 *          gg[3][3], cc[3]    ==> affine trans from arbitrary el to master\n
 *          \n
 *          loc[4][3]          ==> local ordering of vertices for each face\n
 *          vx[4][3]           ==> vertex coordinate labels\n
 *          nvec[4][3]         ==> normal vectors to the faces\n
 *          evec[6][3]         ==> edge vectors\n
 *          elen[6]            ==> edge vector lengths\n
 *          \n
 *          dimV               ==> number of vertices in the d-simplex\n
 *          dimE               ==> number of edges in the d-simplex\n
 *          dimF               ==> number of faces in the d-simplex\n
 *          dimS               ==> number of simplices in the d-simplex (=1)\n
 *          \n
 *          sid                ==> global simplex ID\n
 *          vid[4]             ==> global vertex IDs\n
 *          fid[4]             ==> global face IDs\n
 *          eid[6]             ==> global edge IDs\n
 *          \n
 *          stype              ==> global simplex type\n
 *          vtype[4]           ==> global vertex types\n
 *          ftype[4]           ==> global face types\n
 *          etype[6]           ==> LOCAL edge types\n
 *          \n
 *          *s                 ==> pointer to the simplex\n
 *          *v[4]              ==> pointers to vertices of the simplex
 * @return  None
 * @param   thee  Pointer to class Gem
 * @param   sm    Pointer to a simplex
 * @param   t     Pointer to class TT
 */
VEXTERNC void Gem_simplexInfo(Gem *thee, SS *sm, TT *t);

/**
 * @ingroup Gem
 * @brief   Build the complete master-to-element transformatio.
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemg.c) \n
 *          We just call Gem_simplexInfo to build the basic information,
 *          and then we compute the transformation and some additional
 *          things like the inverse transformations and various determinants.
 * @return  None
 * @param   thee  Pointer to class Gem
 * @param   sm    Pointer to a simplex
 * @param   t     Pointer to class TT
 */
VEXTERNC void Gem_buildVolumeTrans(Gem *thee, SS *sm, TT *t);

/**
 * @ingroup Gem
 * @brief   Build the complete masterFace-to-elementFace transformation.
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemg.c) \n
 *          We ASSUME that Gem_simplexInfo has already been called to build
 *          the basic information in the TT structure we are given.  
 *          We then compute some additional things here for a SINGLE face 
 *          specified by "iface", such as the inverse transformations and     
 *          various determinants.                                              
 *          \n
 *          We do not actually have to assume that Gem_buildVolumeTrans
 *          has been previously called, only that Gem_simplexInfo has 
 *          been called.  However, it seems that all situations that 
 *          occur result in Gem_buildVolumeTrans being called for an
 *          element before Gem_buildSurfTrans is called on any face.
 * @return  None
 * @param   thee  Pointer to class Gem
 * @param   iface index for the faces in a simplex
 * @param   t     index for class TT
 */
VEXTERNC void Gem_buildSurfaceTrans(Gem *thee, int iface, TT *t);

/**
 * @ingroup Gem
 * @brief   Calculate the edge lengths of a simplex.  
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemg.c) 
 * @return  the edge lengths of a simplex
 * @param   thee  Pointer to class Gem
 * @param   v0    Pointer to the first vertex
 * @param   v1    Pointer to the second vertex
 */
VEXTERNC double Gem_edgeLength(Gem *thee, VV *v0, VV *v1);

/**
 * @ingroup Gem
 * @brief   Determine the edge of a simplex opposite the newest vertex. 
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemg.c) 
 * @return  The permutation map of the edge of a simplex opposite the newest 
 *          vertex
 * @param   thee  Pointer to class Gem
 * @param   sm    Pointer to the simplex
 * @param   face  index for the face
 * @param   len   Pointer to the current longest edge length
 */
VEXTERNC int Gem_newestVertex(Gem *thee, SS *sm, int face, double *len);

/**
 * @ingroup Gem
 * @brief   Determine the longest edge of a simplex or a simplex face. 
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemg.c) \n
 *          It is critical to have a consistent tie-breaking rule in order 
 *          to guarantee that recursive refinement procedures: \n
 *          (1) produce conforming meshes (two simplices will refine the  
 *              same edge of a shared face) \n
 *          (2) terminate in finite steps (due to (1)).           
 * @return  The permutation map of the longest edge of a simplex or a simplex face.
 * @param   thee  Pointer to class Gem
 * @param   sm    Pointer to the simplex
 * @param   face  index for the face
 * @param   len   Pointer to the current longest edge length
 */
VEXTERNC int Gem_longestEdge(Gem *thee, SS *sm, int face, double *len);

/**
 * @ingroup Gem
 * @brief   Determine the shortest edge of a simplex or a simplex face. 
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemg.c) \n
 *          It is critical to have a consistent tie-breaking rule in order 
 *          to guarantee that recursive refinement procedures: \n
 *          (1) produce conforming meshes (two simplices will refine the  
 *              same edge of a shared face) \n
 *          (2) terminate in finite steps (due to (1)).           
 * @return  The permutation map of the shortest edge of a simplex or a simplex face.
 * @param   thee  Pointer to class Gem
 * @param   sm    Pointer to the simplex
 * @param   face  index for the face
 * @param   len   Pointer to the current longest edge length
 */
VEXTERNC int Gem_shortestEdge(Gem *thee, SS *sm, int face, double *len);

/**
 * @ingroup Gem
 * @brief   Calculate the shape quality measure for this simplex.  
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemg.c) \n
 * @verbatim
 *           Let |s| denote the volume (which may be negative) of a given
 *           d-simplex "s", let v_i (i=0,...,d) denote the vertices of s,
 *           and let e_{ij} denote the d-vectors representing the 3 or 6 edges
 *           of s that connect v_i to v_j.  We compute the following shape
 *           quality measure for the simplex s:
 *
 *                           f(s,d)   2^{2(1-1/d)} * 3^{(d-1)/2} * |s|^{2/d}
 *               meas(s,d) = ------ = --------------------------------------
 *                           g(s,d)       \sum_{0<=i<j<=d} |e_{ij}|^2
 *
 * 2D Notes: The shape function meas(s,2) is (nearly) the same one used by
 *           Randy Bank and Kent Smith in their joint paper on mesh smoothing:
 *
 *                           f(s,2)         2 * 3^{1/2} * |s|
 *               meas(s,2) = ------ = ---------------------------
 *                           g(s,2)   \sum_{0<=i<j<=2} |e_{ij}|^2
 *
 *           It has the property that its maximal value of 1 is obtained
 *           for an equilateral triangle, and it is scaling invariant, i.e.,
 *           we don't have to worry about the size of the triangle.
 *           Their original function was normalized (their numerator was
 *           2*f(s,2) above) to yield a value of 1 for a equalateral triangle
 *           (with volume 1); this is modified to yield a maximal value
 *           of 1 for the unit triangle (with volume 1/2) to work better with
 *           the unit triangle code.  To effect this slightly different
 *           normalization, the numerator of the quality function was changed
 *           to 2(3)^{1/2}|s|.
 *
 * 3D Notes: The shape function meas(s,3) is (nearly) the same one used by
 *           Joe and Liu in their paper on quality measures for tetrahedra:
 *
 *                           f(s,3)     2^{4/3} * 3 * |s|^{2/3}
 *               meas(s,3) = ------ = ---------------------------
 *                           g(s,3)   \sum_{0<=i<j<=3} |e_{ij}|^2
 *
 *           It is also scaling invariant, so we don't have to worry about
 *           the size of the tetrahedron.  Their original function was 
 *           normalized (their numerator was 12(3|s|)^{2/3}) to yield a
 *           maximal value of 1 for a regular tetrahedron (with volume 1);
 *           this was modified to yield a maximal value of 1 for the unit
 *           tetrahedron (with volume 1/6) to work better with the unit
 *           tetrahedron code.  To effect this slightly different
 *           normalization, the numerator of the quality function was changed
 *           to f(s,3) = 12(3|s|/6)^{2/3} = 2^{4/3}*3*|s|^{2/3}, as above.
 * @endverbatim
 * @return  the shape quality measure for this simplex
 * @param   thee  Pointer to class Gem
 * @param   sm    Pointer to the simplex
 * @param   f     Pointer to volume scaling
 * @param   g     Pointer to the sum of d-vectors representing the 3 or 6 edges
 */
VEXTERNC double Gem_shapeMeasure(Gem *thee, SS *sm, double *f, double *g);

/**
 * @ingroup Gem
 * @brief   Calculate gradient of the shape quality measure for this simplex,  
 *          where the last vertex (vertex d) is treated as the set of
 *          independent variables.   
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemg.c)
 * @return  None
 * @param   thee  Pointer to class Gem
 * @param   sm    Pointer to the simplex
 * @param   vmap  the array of the map
 * @param   dm    gradient of the shape quality measure for this simplex
 */
VEXTERNC void Gem_shapeMeasureGrad(Gem *thee, SS *sm, int vmap[], double dm[]);

/**
 * @ingroup Gem
 * @brief   Calculate ratio of longest-to-shortest edge of a simplex.       
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemg.c) 
 * @return  ratio of longest-to-shortest edge of a simplex.
 * @param   thee  Pointer to class Gem
 * @param   sm    Pointer to the simplex
 */
VEXTERNC double Gem_edgeRatio(Gem *thee, SS *sm);

/**
 * @ingroup Gem
 * @brief   Calculate the determinant of the transformation from the   
 *          master element to this element.   
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemg.c) \n
 *          We just call Gem_simplexInfo to build the basic information,  
 *          and then we build the transformation and compute the determinant 
 *          here. \n
 *          For a manifold, we need to know the orientation of the manifold 
 *          in order to decide what is counter-clock-wise, and what is 
 *          clockwise, in terms of vertex orderings.  One will lead to a   
 *          positive volume, and the other to a negative volume, when we  
 *          compute volume using the determinant of the jacobian of the 
 *          affine transformation to the master element.  We assume here
 *          that the uniform chart computed by Gem_simplexInfo is such that
 *          the orientation in the chart reflects the correct orientation
 *          of the manifold (locally).
 * @return  the determinant of the transformation from the master element to this
 *          element.
 * @param   thee  Pointer to class Gem
 * @param   sm    Pointer to the simplex
 */
VEXTERNC double Gem_simplexVolume(Gem *thee, SS *sm);

/**
 * @ingroup Gem
 * @brief   Traverse the simplices and check their shapes.   
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemg.c)
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_shapeChk(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Produce the initial edge markings in all of the simplices.
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemg.c)
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_markEdges(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Go through simplices and enforce a vertex ordering that
 *          will produce a positive determinant.
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemg.c)\n
 *          This relies on an embedding of R2 into R3 (or R3 into R4)   
 *          and breaks e.g. if this is a non-orientable 2-manifold, etc.
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_reorderSV(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Smooth the mesh using simple Laplace smoothing.
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemg.c)\n
 *          We don't need the simplex rings here, but we need the edge rings.
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_smoothMeshLaplace(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Smooth the mesh using simple Laplace smoothing.  
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemg.c)\n
 *          We don't need the simplex rings here, but we need the edge rings.
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_smoothMesh(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Smooth the boundary mesh. 
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemg.c)
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_smoothMeshBnd(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Smooth the mesh using a volume optimization approach.   
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemg.c)
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_smoothMeshOpt(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Unify the charts of vertices.
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemg.c)
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_buildCharts(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Unify the charts of vertices.  
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemg.c)
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_clearCharts(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Print the exact current malloc usage.  
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemchk.c)
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_memChk(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Print the exact current malloc usage: vertices.
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemchk.c)
 * @return  None
 * @param   thee     Pointer to class Gem
 * @param   num      the global "T" counter -- how many "T"s in list
 * @param   size     size of the object in bytes
 * @param   vecUse   total object size in the the global "T" counter
 * @param   vecMal   total size of allocated blocks
 * @param   vecOhd   max size of blocks
 */
VEXTERNC void Gem_memChkVV(Gem *thee, int *num,
    int *size, int *vecUse, int *vecMal, int *vecOhd);

/**
 * @ingroup Gem
 * @brief   Print the exact current malloc usage: edges.   
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemchk.c)
 * @return  None
 * @param   thee     Pointer to class Gem
 * @param   tnum     the global "T" counter -- how many "T"s in list
 * @param   size     size of the object in bytes
 * @param   vecUse   total object size in the the global "T" counter
 * @param   vecMal   total size of allocated blocks
 * @param   vecOhd   max size of blocks
 */
VEXTERNC void Gem_memChkEE(Gem *thee, int *tnum,
    int *size, int *vecUse, int *vecMal, int *vecOhd);

/**
 * @ingroup Gem
 * @brief   Print the exact current malloc usage: simplices. 
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemchk.c)
 * @return  None
 * @param   thee     Pointer to class Gem
 * @param   tnum     the global "T" counter -- how many "T"s in list
 * @param   size     size of the object in bytes
 * @param   vecUse   total object size in the the global "T" counter
 * @param   vecMal   total size of allocated blocks
 * @param   vecOhd   max size of blocks
 */
VEXTERNC void Gem_memChkSS(Gem *thee, int *tnum,
    int *size, int *vecUse, int *vecMal, int *vecOhd);

/**
 * @ingroup Gem
 * @brief   Print the exact current malloc usage: simplex queues. 
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemchk.c)
 * @return  None
 * @param   thee     Pointer to class Gem
 * @param   currentQ index of a given queue
 * @param   tnum     the global "T" counter -- how many "T"s in list
 * @param   tsize    size of the object in bytes
 * @param   tVecUse  total object size in the the global "T" counter
 * @param   tVecMal  total size of allocated blocks
 * @param   tVecOhd  max size of blocks
 */
VEXTERNC void Gem_memChkSQ(Gem *thee, int currentQ,
    int *tnum, int *tsize, int *tVecUse, int *tVecMal, int *tVecOhd);

/**
 * @ingroup Gem
 * @brief   Estimate the current RAM usage.
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemchk.c)
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_memChkMore(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Calculate the cost to traverse the various structures.     
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemchk.c)
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_speedChk(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Check the self-consistency of the geometry datastructures.      
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemchk.c)
 * @return  None
 * @param   thee  Pointer to class Gem
 * @param   key   0 --> check: min (just vertices and simplices) \n
 *                1 --> check: min + simplex ring \n
 *                2 --> check: min + simplex ring + edge ring \n
 *                3 --> check: min + simplex ring + edge ring + conform
 */
VEXTERNC void Gem_formChk(Gem *thee, int key);

/**
 * @ingroup Gem
 * @brief   Print out contents of all geometry structures.  
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemchk.c)
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_contentChk(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Check some structures.
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemchk.c)
 * @return  None
 * @param   thee  Pointer to class Gem
 * @param   key   0 --> check: min (just vertices and simplices) \n
 *                1 --> check: min + simplex ring \n
 *                2 --> check: min + simplex ring + edge ring \n
 */
VEXTERNC void Gem_ramClear(Gem *thee, int key);

/**
 * @ingroup Gem
 * @brief   Force naborless faces to become boundary faces.    
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemchk.c)
 * @return  None
 * @param   thee  Pointer to class Gem
 * @param   btype 0 --> create interior boundary faces \n
 *                1 --> create boundary faces of type "1"\n
 *                2 --> create boundary faces of type "2"
 */
VEXTERNC void Gem_makeBnd(Gem *thee, int btype);

/**
 * @ingroup Gem
 * @brief   Mark selected boundary faces in a special way.   
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemchk.c)
 * @return  None
 * @param   thee  Pointer to class Gem
 * @param   key   index for selected boundary faces
 *          key==0 --> check: min (just vertices and simplices)
 *          key==1 --> check: min + simplex ring
 *          key==2 --> check: min + simplex ring + edge ring
 *          key==3 --> check: min + simplex ring + edge ring + conform
 */
VEXTERNC void Gem_makeBndExt(Gem *thee, int key);

/**
 * @ingroup Gem
 * @brief   Incremental flip Delaunay generator. 
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemgen.c)\n
 * @verbatim
 *           We use an incremental flip Delaunay mesh generator.
 *           In 2D, this is based on the standard edge-flipping.
 *           In 3D, this is based on 2-to-3 flips (shared face-to-edge flip),
 *           and 3-to-2 flips (a restricted edge-to-face flip).
 *
 *           The 2D version of the algorithm here is similar to several
 *           edge-flip-based Delaunay algorithms in the literature.
 *
 *           The 3D version of the algorithm here is similar to Barry Joe's
 *           and Ernst Mucke's.  However, the ringed-vertex datastructure used
 *           here leads to a more concise implementation, as compared to
 *           implementations using e.g. Mucke's edge-facet datastructure.
 *          
 *           Both this algorithm and similar incremental flip algorithms in
 *           the 3D case are based on several recent theoretical results
 *           due to Barry Joe, which guarantee the following:
 *
 *               (1) The flipping algorithm to re-establish Delaunay-ness
 *                   always terminates in a finite number of steps.
 *
 *               (2) The algorithm works regardless of the flipping order.
 *
 *               (3) Only the exterior faces of the star region of the new
 *                   vertex (what Mucke calls "link facets", because these
 *                   "facets" link two triangles or tets together), and in the
 *                   3D case the three edges which make up those faces, need
 *                   be tested and then possibly flipped.
 *
 *           The algorithm is as follows:
 *
 *               (1) Given N inputs points, a single enclosing simplex is
 *                   formed by adding d+1 additional points, and then forming
 *                   that single simplex.
 *
 *               (2) The N points are then added to the mesh one at a time
 *                   by locating the simplex containing each point, adding
 *                   the point, and splitting the containing simplex into
 *                   d+1 children (note that unless the points are in
 *                   "general" position, this may give rise to degenerate
 *                   simplices).
 *
 *               (3) The link-facets of the newly added vertex are checked
 *                   for Delaunayness.  The link-facets are the faces 
 *                   opposite the new vertex in each simplex that uses
 *                   the new vertex.  (This is the boundary of the support
 *                   region for a finite element basis function, for example.)
 *                   If a non-Delaunay face is located, we attempt to flip
 *                   one of the three edges of the face using a 3-to-2
 *                   (edge-to-face) flip.  We only flip such as edge if the
 *                   simplex ring about the edge has length three.  If no
 *                   such edge flips are possible, we do a 2-to-3 flip
 *                   (face-to-edge) if this is possible (it is only possible
 *                   if the two tets sharing the face form a convex region).
 *                   If no flipping is possible, we temporarily ignore this
 *                   non-Delaunay face and move on to the next face.
 *
 *               (4) Step (3) above terminates in a finite number of steps
 *                   thanks to the theoretical results of Barry Joe.
 * @endverbatim
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_delaunay(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Edge or face flip for the incremental flip algorithm.
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemgen.c)
 * @return  None
 * @param   thee  Pointer to class Gem
 * @param   vx    Pointer to the vertex
 */
VEXTERNC void Gem_flip(Gem *thee, VV *vx);

/**
 * @ingroup Gem
 * @brief   Find a simplex containing a given vertex.  
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemgen.c)
 * @return  Pointer to the simplex
 * @param   thee  Pointer to class Gem
 * @param   vx    Pointer to the vertex
 */
VEXTERNC SS* Gem_findSimplex(Gem *thee, VV *vx);

/**
 * @ingroup Gem
 * @brief   Initialize the geometric predicates.   
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemgen.c)
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_predinit(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Determine the orientation of the vertices in a simplex. 
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemgen.c)
 * @return  Success enumeration
 * @param   thee  Pointer to class Gem
 * @param   sm    Pointer to the simplex
 */
VEXTERNC int Gem_orient(Gem *thee, SS *sm);

/**
 * @ingroup Gem
 * @brief   Determine if a vertex lies in a sphere of other vertices. 
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemgen.c)
 * @return  Success enumeration
 * @param   thee      Pointer to class Gem
 * @param   sm        pointer to the simplex
 * @param   sm_facet  index for the face
 * @param   vx        pointer to the vertex
 * @param   vxnb      pointer to the given vertex 
 */
VEXTERNC int Gem_inSphere(Gem *thee, SS *sm, int sm_facet, VV *vx, VV *vxnb);

/**
 * @ingroup Gem
 * @brief   Determine whether or a not a point is in a simplex.
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemgen.c)
 * @return  Success enumeration
 * @param   thee  Pointer to class Gem
 * @param   sm    Pointer to the simplex
 * @param   x     arrary of the point position
 */
VEXTERNC int Gem_pointInSimplex(Gem *thee, SS *sm, double x[]);

/**
 * @ingroup Gem
 * @brief   Evaluate basis functions at a point in a simplex.      
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemgen.c)
 * @return  Success enumeration
 * @param   thee  Pointer to class Gem
 * @param   sm    Pointer to the simplex
 * @param   x     arrary of the point position
 * @param   phi   basis functions at a point in a simplex
 * @param   phix  derivs of basis functions at a point in a simplex
 */
VEXTERNC int Gem_pointInSimplexVal(Gem *thee, SS *sm, double x[],
    double phi[], double phix[][3]);

/**
 * @ingroup Gem
 * @brief   Edge or face flip for the incremental flip algorithm.        
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemgen.c)
 * @return  Success enumeration
 * @param   thee  Pointer to class Gem
 * @param   dimX  the intrinsic spatial dimension
 * @param   defX  Pointer to vertex deformation or displacement values
 */
VEXTERNC int Gem_deform(Gem *thee, int dimX, double *defX[MAXV]);

/**
 * @ingroup Gem
 * @brief   Mark simplices to be refined.    
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemref.c)
 * @return  number of marked simplices
 * @param   thee   Pointer to class Gem
 * @param   key    If (key == -1)  Clear all simplex refinement flags.\n
 *                 If (key == 0)   Mark all simplices for refinement.\n
 *                 If (key == 1)   Mark special simplices for a testcase 
 *                                 refinement.
 * @param   color  the chart of the simplex
 */
VEXTERNC int Gem_markRefine(Gem *thee, int key, int color);

/**
 * @ingroup Gem
 * @brief   Refine the manifold and also build a prolongation operator that
 *          can interpolate functions from the original manifold to the      
 *          new manifold.
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemref.c)\n
 *           Longest Edge, Newest Vertex, or Newest Pair, is used to bisect
 *           a single simplex in an asymptotically non-degenerate way in the
 *           "bisect[LE,NV,NP]" routines which are called from this routine.
 *           Marked simplices are subdivided into 2/4/8 child simplices.
 *           A closure algorithm is performed which continues subdivision
 *           until a conforming mesh is produced.  Boundary
 *           nodes/edges[faces] are correctly refined.
 *           \n\n
 *           We purposely do the following trick in order to facilitate the
 *           construction of a prolongation operator after refinement,
 *           if it is so desired.  We begin the refinement with no edges;
 *           only the conforming mesh as described by the list of vertices
 *           and the list of simplices using the vertices.  When a simplex
 *           is to be subdivided, the refinement edge (or edges) is then
 *           identified, and then created on the fly.  The new vertex which
 *           is created by the refinement, namely the midpoint of the edge,
 *           is then stored with the newly created edge.  The edge is added
 *           to the ring of edges around each of its two vertices.  This
 *           allows our refinement algorithm to easily detect whether or
 *           not an edge has already been refined by a naboring simplex by
 *           simply traversing the edge lists of the two vertices on the
 *           refinement edge.  If the edge exists, then it must have already
 *           been refined, since it is only created in order to refine it.
 *           Moreover, the new vertex at the midpoint of the edge is then
 *           also directly available from the edge structure for use in
 *           building the children simplices, without having to search for it.
 *           (The edge datastructure can be viewed as simply a holding cell
 *           for the newly created vertices so that they can be found without
 *           any searching.) 
 *           \n\n
 *           How does this help us to later build an appropriate prolongation
 *           operator between the original mesh and the final refined mesh?
 *           While we begin the refinement with no edges, we end with a
 *           list of edges that is precisely the set of edges that were
 *           refined.  Let us order the function values of a mesh function
 *           on the fine mesh (with function values at vertices) in the
 *           following order: vertices common to the coarse mesh in the same
 *           order as the coarse mesh, followed by vertices at the midpoints
 *           of the refined mesh, in the order of the edges in the edge list.
 *           The linear prolongation operator (for example) which would
 *           linearly interpolate a function from the original coarse mesh
 *           to the refined mesh is then a block 2x1 matrix.  The upper block
 *           is a square identity matrix with number of rows/columns equal
 *           to the number of vertices in the original mesh; it is completely
 *           clear how to build this upper block.  The lower block is a
 *           (generally) rectangular matrix, with number of rows equal to the
 *           number of edges that were refined.  Since we finish refinement
 *           with precisely the refined edges in the edge list, he can
 *           simply traverse the edge list to build the lower block of the
 *           prolongation matrix.  In particular, in the linear interpolation
 *           case, each row of the lower block will be zero, except for two
 *           columns, corresponding to the vertex numbers of the two coarse
 *           mesh vertices which lie on the ends of the edge that was refined.
 *           A value of 0.5 is then placed in those two columns.
 *           \n\n
 *           In the case of linear prolongation, the lower block of the
 *           prolongation matrix has exactly one row for each edge that was
 *           refined, with zeros as every entry except for the two columns
 *           corresponding to the vertex numbers that were on each end of the
 *           edge that was bisected.
 *           \n\n
 *           There is an opportunity for a problem with this approach; if an
 *           edge is multiply refined, then we must keep track of all of the
 *           resulting edges and their parent-child relationships, in order
 *           to build the correct interpolation.  The lower block of the
 *           prolongation matrix will now be slightly more complicated than
 *           described above.
 * @return  number of refined simplices
 * @param   thee  Pointer to class Gem
 * @param   rkey  If (rkey==0) Perform recursive simplex bisection until 
 *                conformity\n
 *                If (rkey==1) Perform first quadra-[octa-]-section, followed 
 *                by recursive simplex bisection until conformity.\n
 *                IMPORTANT NOTE: In 2D, (rkey==1) WILL generate        
 *                a conforming mesh.  However, in 3D, this procedure    
 *                will in general produce nonconforming simplices.      
 *                To produce a conforming mesh in 3D would require an   
 *                implementation covering all possible face refinement  
 *                combinations (something like 169 cases).  This has    
 *                been done e.g. by Jurgen Bey in AGM, but we are       
 *                not that patient; use (rkey==0) above if you want     
 *                a conforming mesh...\n
 *                If (rkey==2) As a test of the conformity procedure, 
 *                perform quadra-[octa-]-section until conformity, which        
 *                should produce a uniformly regularly refined mesh.    
 *                (In 2D, each triangle should be divided into four     
 *                children, and in 3D each tetrahedron should be        
 *                divided into eight children.)                   
 * @param   bkey  Boolean sets the bkey type for bisecting the mesh
 *           If (bkey==0) Bisection type: Longest Edge
 *           If (bkey==1) Bisection type: Newest Vertex
 *           If (bkey==2) Bisection type: Newest Pair
 * @param   pkey  Boolean sets the pkey type to prolongate a vector 
 */
VEXTERNC int Gem_refine(Gem *thee, int rkey, int bkey, int pkey);

/**
 * @ingroup Gem
 * @brief   We do three things in this routine:\n
 *          (1) Find the midpoint of existing edge, or create it\n
 *          (2) Tell simplices using edge that they are now nonconforming\n
 *          (3) Determine the "type" of the new point by using the type  
 *              of the edge.  The edge type must itself be calculated   
 *              on the fly, because we allow the use of lazy edge creation. 
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemref.c)\n\n
 *           Note that we MUST determine the type of the new point
 *           (interior or boundary) based effectively on EDGE types
 *           rather than vertex or face types.  It is easy to construct
 *           examples where typing based on vertex type can mark a new 
 *           interior point falsely as a boundary point.  Note that typing by
 *           the faces of a single simplex which uses the bisection edge
 *           can also be fooled in 3D (it is foolproof in 2D).
 *           For example, it may be the case in 3D that an edge of a tet 
 *           touches a boundary, but none of its faces are boundary faces.
 *           In this case, the new point would be marked (incorrectly)
 *           as an interior point.
 *           \n\n
 *           The solution to this problem is to determine the type of the
 *           new point by using the type of the edge.  The only problem we
 *           then face is how to do this without actually having edges around.
 *           In other words, we must determine the correct edge type on the 
 *           fly.  This can be handled by looking at all faces of all simplices 
 *           which use the edge, and applying the following rules:
 *           (1) If all faces are interior, the edge is interior
 *           (3) If at least one face is boundary, the edge is boundary
 *           \n\n
 *           Note that since we must look at all simplices on the ring
 *           around the edge anyway to handle the conformity situation,
 *           we don't have to do any additional work to determine the
 *           correct edge type.  Therefore, we will take this approach
 *           at determining the correct edge type on the fly, EVEN IF
 *           the edges are always around and their correct types are 
 *           recorded correctly once and for all when a mesh is built.
 *           \n\n
 *           This way we can also do lazy edge creation; i.e., create an 
 *           edge only when it needs to be refined.  The lazy edge is then
 *           in principle simply a holder for the new point, allowing O(1)
 *           access to the new point by other simplices, through the edge
 *           rings around their vertices.
 *           \n\n
 *           Note that lazy edge creation has a serious performance benefit
 *           to this routine in particular: if all edges are around, then to
 *           find a particular edge, we then always have to search both edge
 *           rings associated with each vertex for a common edge.  This means
 *           we look for the intersection of two sets of five elements on
 *           average in 2D, and two sets of fifteen elements on average in 3D.  
 *           With lazy edge creation, we search only through lists of edges
 *           that were created for refinement; these lists are usually a
 *           much smaller.
 *           \n\n
 *           A final advantage of lazy edge creation is that having a list
 *           of only the refined edges allows us to efficiently build a
 *           prolongation operator between the original mesh and the refined
 *           mesh.
 * @return  None
 * @param   thee      Pointer to class Gem
 * @param   currentQ  index of a given queue
 * @param   sm        pointer to the simplex
 * @param   v         pointer to the vertex
 * @param   vAB       pointer to the midpoint of the edge
 * @param   A         index for a vertex in a simplex
 * @param   B         index for a vertex in a simplex
 */
VEXTERNC void Gem_refineEdge(Gem *thee, int currentQ,
    SS *sm, VV *v[4], VV **vAB, int A, int B);

/**
 * @ingroup Gem
 * @brief   Uniform regular (quadrasection) refinement of a single simplex.  
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemref.c)\n
 *          Boundary nodes/edges[faces] are correctly refined.  
 * @return  None
 * @param   thee      Pointer to class Gem
 * @param   sm        pointer to the simplex
 * @param   currentQ  index of a given queue
 */
VEXTERNC void Gem_octsect(Gem *thee, SS *sm, int currentQ);

/**
 * @ingroup Gem
 * @brief   Bisection refinement of a single simplex by longest edge. 
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemref.c)\n
 *          Boundary nodes/edges[faces] are correctly refined.\n
 *          "Face Type rule": Boundary faces rule over Interior faces.\n
 *          In other words, if an edge is shared between a boundary face and
 *          an interior face, and the edge gets refined, the new point will 
 *          be a boundary point. 
 * @return  None
 * @param   thee      Pointer to class Gem
 * @param   sm        pointer to the simplex
 * @param   currentQ  index of a given queue
 */
VEXTERNC void Gem_bisectLE(Gem *thee, SS *sm, int currentQ);

/**
 * @ingroup Gem
 * @brief   Bisection refinement of a single simplex by newest vertex.
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemref.c)\n\n
 * @verbatim
 *           Boundary nodes/edges[faces] are correctly refined.
 *
 *           "Face Type rule": Boundary faces rule over Interior faces.
 *
 *           In other words, if an edge is shared between a boundary face and
 *           an interior face, and the edge gets refined, the new point will
 *           be a boundary point.
 *
 *           We use "newest vertex" approach to selecting the refinement edge
 *           choice generated by the Arnold-Mukherjee marking procedure.
 *           Below is a mathematica program that demonstrates the marking
 *           algorithm, courtesy of Doug Arnold and Arup Mukerjee.
 *
 *    (* -- BEGIN MATHEMATICA CODE ----------------------------------------- *)
 *    (*
 *     * NOTE: The following mathematica program, courtesy of doug arnold
 *     * and arup mukherjee, illustrates their particular edge choice,
 *     * which is provably non-degenerate.  The small comment at the
 *     * beginning, just following my comment here, is from Arup.
 *     * We implement their marking procedure in MC with our MC geometry
 *     * datastructures (as an alternative to longest edge choice) by 
 *     * allocating a few bits to indicate both the marked edges on each face
 *     * and then the refinement edge choice.  -mike
 *     *)
 *    (*
 *     * Arup's Comment:  I am also including a mathematica script for the 
 *     * bisection.  Given a tet and its "type" it does the bisection upto a 
 *     * required bisection level and lists the total number of similarity 
 *     * classes produced in the process (it uses a Liu-Joe indicator given on
 *     * one of their papers ... the example at the end of Chap 3 in the 
 *     * thesis has the specific reference on this). The tetrahedron specified
 *     * in the example attains the upper bound of 36 similarity classes.
 *     *)
 *    
 *    (* SET the value of nlevels (number of bisection levels)
 *    and the TYPE of t[0] p-uf, p-f, np-aa, etc ...
 *     before input
 *    to mathematica *)
 *    
 *    nlevels=8
 *    
 *    SetAttributes[bisect,Listable]
 *    (* these are the "bisection" rules 
 *      p-uf -- planar with flag off 
 *        v0-v1 is refinement edge and v0-v2 and v1-v2 are marked edges
 *          the plane containing all marked edges is v0 v1 v2
 *      p-f -- planar with flag on 
 *        v0-v1 is refinement edge and v0-v2 and v1-v2 are marked edges
 *          the plane containing all marked edges is v0 v1 v2
 *      np-aa -- non planar case with adjacent markings
 *        v0-v1 is refinement edge and v0-v2 and v1-v3 are marked edges 
 *      np-oo -- opp-opp  
 *        v0-v1 is refinement edge and v2-v3 is the other marked edge
 *      np-ao -- adj-opp
 *        v0-v1 is the refinement edge and v0-v2 and v2-v3 are marked
 *    *)
 *    
 *    bisect[Tetra[v0_,v1_,v2_,v3_,p-uf]]:= {
 *        Tetra[v0,v2,v3,(v0+v1)/2,p-f], Tetra[v1,v2,v3,(v0+v1)/2,p-f]  }
 *    
 *    bisect[Tetra[v0_,v1_,v2_,v3_,p-f]]:= {
 *        Tetra[v0,v2,v3,(v0+v1)/2,np-aa], Tetra[v1,v2,v3,(v0+v1)/2,np-aa] }
 *    
 *    bisect[Tetra[v0_,v1_,v2_,v3_,np-aa]]:={
 *        Tetra[v0,v2,v3,(v0+v1)/2,p-uf], Tetra[v1,v3,v2,(v0+v1)/2,p-uf] }
 *    
 *    bisect[Tetra[v0_,v1_,v2_,v3_,np-oo]]:={
 *        Tetra[v2,v3,v0,(v0+v1)/2,p-uf], Tetra[v2,v3,v1,(v0+v1)/2,p-uf] }
 *    
 *    bisect[Tetra[v0_,v1_,v2_,v3_,np-ao]]:={
 *        Tetra[v0,v2,v3,(v0+v1)/2,p-uf], Tetra[v2,v3,v1,(v0+v1)/2,p-uf] }
 *    
 *    (* a "generic" tetrahedron *)
 *    v0={0,0,0};
 *    v1={23,0,0};
 *    v2={7,0,11};
 *    v3={17,5,13};
 *    
 *    Clear[t]
 *    (* set t[0] to be pf, pt, or np-** as the case may be *)
 *    t[0]={Tetra[v0,v1,v2,v3,p-uf]};
 *    t[n_]:=t[n]=bisect[t[n-1]]
 *    
 *    (* The Liu-Joe quality indicator *)
 *    dist2[v0_,v1_] := (v0-v1).(v0-v1)
 *    dist[v0_,v1_] := Sqrt[dist[v0,v1]]
 *    SetAttributes[qual,Listable]
 *    qual[Tetra[v0_,v1_,v2_,v3_,any_]] :=
 *    12 Abs[Det[{v1-v0,v2-v0,v3-v0}]/2]^(2/3)/
 *    (dist2[v0,v1]+dist2[v0,v2]+dist2[v0,v3]
 *    +dist2[v1,v2]+dist2[v1,v3]+dist2[v2,v3])
 *    
 *    (* discretized Liu-Joe quality indicator *)
 *    (* (scaled to [0,100000] and rounded to an integer *)
 *    SetAttributes[dqual,Listable]
 *    dqual[t_] := Round[100000 N[qual[t],10]]
 *    
 *    (* q[i] --- list of qualities for level i
 *       qq[i] -- list of all qualities upto level i
 *       newout[i] -- the "number" of different similarity classes at 
 *                    level i
 *       totout[i] -- the "number" of different similarity classes at
 *                    or below level i    *)
 *    Do[q[i]=dqual[t[i]],{i,0,nlevels}]
 *    Do[qq[i]= Union@@Table[q[j],{j,0,i}],{i,0,nlevels}]
 *    Do[newout[i]=Dimensions[Union[Flatten[q[i]]]],{i,0,nlevels}]
 *    Do[totout[i]=Dimensions[Union[Flatten[qq[i]]]],{i,0,nlevels}]
 *    Table[{i,newout[i],totout[i]},{i,0,nlevels}]//TableForm
 *    (* -- END   MATHEMATICA CODE ---------------------------------------- *)
 * @endverbatim
 * @return  None
 * @param   thee      Pointer to class Gem
 * @param   sm        pointer to the simplex
 * @param   currentQ  index of a given queue
 */
VEXTERNC void Gem_bisectNV(Gem *thee, SS *sm, int currentQ);

/**
 * @ingroup Gem
 * @brief   Bisection refinement by pairs.    
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemref.c)\n
 *          Boundary nodes/edges[faces] are correctly refined.\n
 *          "Face Type rule": Boundary faces rule over Interior faces.\n
 *          In other words, if an edge is shared between a boundary face and
 *          an interior face, and the edge gets refined, the new point will 
 *          be a boundary point. 
 * @return  None
 * @param   thee      Pointer to class Gem
 * @param   sm        pointer to the simplex
 * @param   currentQ  index of a given queue
 */
VEXTERNC void Gem_bisectNP(Gem *thee, SS *sm, int currentQ);

/**
 * @ingroup Gem
 * @brief   Un-refine the mesh. 
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemunref.c)\n
 *          If (key==0) Simply toss out the marked simplices; this leaves
 *          "holes", and we mark all neighbor faces that were
 *          revealed as boundary faces. 
 * @return  number of unrefined mesh elements
 * @param   thee  Pointer to class Gem
 * @param   rkey  If (rkey==0) Perform recursive simplex bisection until 
 *                conformity\n
 *                If (rkey==1) Perform first quadra-[octa-]-section, followed 
 *                by recursive simplex bisection until conformity.\n
 *                IMPORTANT NOTE: In 2D, (rkey==1) WILL generate        
 *                a conforming mesh.  However, in 3D, this procedure    
 *                will in general produce nonconforming simplices.      
 *                To produce a conforming mesh in 3D would require an   
 *                implementation covering all possible face refinement  
 *                combinations (something like 169 cases).  This has    
 *                been done e.g. by Jurgen Bey in AGM, but we are       
 *                not that patient; use (rkey==0) above if you want     
 *                a conforming mesh...\n
 *                If (rkey==2) As a test of the conformity procedure, 
 *                perform quadra-[octa-]-section until conformity, which        
 *                should produce a uniformly regularly refined mesh.    
 *                (In 2D, each triangle should be divided into four     
 *                children, and in 3D each tetrahedron should be        
 *                divided into eight children.)                   
 * @param   pkey  Boolean sets the pkey type to prolongate a vector 
 */
VEXTERNC int Gem_unRefine(Gem *thee, int rkey, int pkey);

/**
 * @ingroup Gem
 * @brief   Delete a simplex cleanly, maintaining a consecutive list 
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemunref.c)
 * @return  None
 * @param   thee      Pointer to class Gem
 * @param   sm        pointer to the simplex
 * @param   currentQ  index of a given queue
 */
VEXTERNC void Gem_delSimplex(Gem *thee, SS *sm, int currentQ);

/**
 * @ingroup Gem
 * @brief   Toss out any hanging vertices. 
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemunref.c)
 * @return  None
 * @param   thee   Pointer to class Gem
 */
VEXTERNC void Gem_unHangVertices(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Read in the user-specified initial vertex-simplex mesh. 
 *          provided to us in MC-Simplex-Format (MCSF), and transform   
 *          into our internal datastructures. 
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemio.c)
 * @verbatim
 *        The user provides the following information about a domain and an
 *        initial simplex-triangulation:
 *
 *        D = DIMENSION  = spatial dimension of problem (1, 2, or 3)
 *        N = NVERTICES  = total number of vertices in mesh
 *        L = NSIMPLICES = total number of simplices in the mesh
 *
 *        double vertex* = List of all vertex coordinates in the form:
 *
 *        vertex[ N * (3+D) ] = {
 *            { id1, chart1, v1_x1, ..., v1_xD },
 *            { id2, chart2, v2_x1, ..., v2_xD },
 *                          ...
 *            { idN, chartN, vN_x1, ..., vN_xD }
 *        }
 *
 *        Here, vj_xk is the kth component (float or double) of the vertex 
 *        with global vertex number j.  The "idj" flag is a 32-bit integer 
 *        "name" for vertex j, "chartj" is a 32-bit number representing the
 *        "chart" with which to interpret the coordinates vj_xk, in the sense
 *        of the charts of an atlas of manifold domain.
 *
 *        int simplex* = List of all simplices by vertex number in the form:
 *
 *        simplex[ L * (2 + 2*(D+1)) ] = {
 *            { id1, g1, m1, f1_1, ..., f1_{D+1}, s1_v1, ..., s1_v{D+1} },
 *            { id2, g2, m2, f2_1, ..., f2_{D+1}, s2_v1, ..., s2_v{D+1} },
 *                                  ...
 *            { idL, gL, mL, fL_1, ..., f2_{D+1}, sL_v1, ..., sL_v{D+1} }
 *        }
 *
 *        Here, sj_vk is a 32-bit integer giving the global vertex number
 *        making up the kth vertex of the simplex with global number j.  The
 *        "idj" flag is a 32-bit integer "name" for simplex j, "gj" is a
 *        32-bit group number associated with simplex j (for grouping subsets
 *        of simplices together for various reasons), "mj" is a 32-bit integer
 *        containing the information about the simplex j such as its material
 *        type, and "fj_k" is a 32-bit integer containing the information 
 *        about each face k of simplex j such as their boundary types (each 
 *        face opposite vertex k in simplex j).
 *
 *        Thus, in 2D, a simplex (triangle) is specified by 3 consecutive 
 *        vertex numbers, and in 3D a simplex (tetrahedra) is specified by 4 
 *        consecutive vertex numbers.  The physical coordinates of any vertex
 *        k in the simplex array are given in the vertex array in the
 *        appropriate row of the array.
 *
 *        NOTE: The ordering of the vertices in a simplex is *extremely* 
 *        important here; see the note below.
 *
 * Ordering of the vertices in a simplex:
 *
 *    1D: Well, this is pretty straight-forward; we will order the vertices
 *        from left-to-right in each simplex; this will produce the correct
 *        sign in integration by parts.
 *
 *    2D: All closed triangles must be specified by three consecutive vertices
 *        in simplex and must be counter-clockwise-ordered by their vertices,
 *        as seen from the "up" side of the 2D body/shell/surface.  This 
 *        produces the correct surface-normals from the right-hand-rule for 
 *        surface (line in 2D) integrals.
 *
 *    3D: All closed tetrahedra must be specified by four consecutive vertices
 *        in "simplex" in the following way:  The first three vertices must 
 *        represent any one of the four faces as a counter-clockwise-ordered 
 *        triangle, as seen from INSIDE the tetrahedra.  I.e., you can think
 *        about this first triangle as lying in the plane, and you are 
 *        standing on the plane looking down at it. The fourth vertex 
 *        specified following the first three in the simplex is then some 
 *        height above the plane containing the first counter-clockwise 
 *        ordered triangle.  This specification allows the remaining three
 *        (of the four) triangles making up any tetrahedra to be correctly 
 *        specified in a counter-clockwise, inward-facing manner, so that the
 *        correct surface-normals can be calculated consistently.
 * @endverbatim
 * @return  Success enumeration
 * @param   thee  Pointer to class Gem
 * @param   key   input format type\n
 *                0 ==> simplex format\n
 *                1 ==> edge format\n
 *                2 ==> simplex-nabor format
 * @param   sock  socket for reading/writing a finite element mesh or mesh function
 */
VEXTERNC int Gem_read(Gem *thee, int key, Vio *sock);

/**
 * @ingroup Gem
 * @brief   Toss out any hanging vertices. 
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemio.c)
 * @return  None
 * @param   thee   Pointer to class Gem
 * @param   key    output format type\n
 *                 0 ==> simplex format\n
 *                 1 ==> edge format\n
 *                 2 ==> simplex-nabor format
 * @param   sock   socket for reading/writing a finite element mesh or mesh function
 * @param   fkey   simplex write option\n
 *                 0 ==> write simplices\n
 *                 1 ==> write only sipmlex boundary faces\n
 *                 2 ==> write only simplices that have
 *                       at least one boundary face 
 *                       (NOT IMPLEMENTED HERE)
 */
VEXTERNC void Gem_write(Gem *thee, int key, Vio *sock, int fkey);

/**
 * @ingroup Gem
 * @brief   Write out the faces of a 3-simplex mesh as a complete and legal
 *          2-simplex mesh in "MCSF" format (described above for Gem_read).
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemio.c)
 * @return  None
 * @param   thee   Pointer to class Gem
 * @param   sock   socket for reading/writing a finite element mesh or mesh function
 */
VEXTERNC void Gem_writeFace3d(Gem *thee, Vio *sock);

/**
 * @ingroup Gem
 * @brief   Read in the user-specified initial vertex-edge mesh
 *          provided to us in MC-Edge-Format (MCEF), and transform
 *          into our internal datastructures.
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemio.c)\n
 * @verbatim
 * Notes: The user provides the following information about a domain and an
 *        initial edge-based triangulation:
 *
 *        D = DIMENSION  = spatial dimension of problem (1, 2, or 3)
 *        N = NVERTICES  = total number of vertices in mesh
 *        L = NEDGES     = total number of edges in the mesh
 *
 *        double vertex* = List of all vertex coordinates in the form:
 *
 *        vertex[ N * (3+D) ] = {
 *            { id1, proc1, info1, v1_x1, ..., v1_xD },
 *            { id2, proc2, info2, v2_x1, ..., v2_xD },
 *                            ...
 *            { idN, procN, infoN, vN_x1, ..., vN_xD }
 *        }
 *
 *        Here, vj_xk is the kth component (float or double) of the vertex 
 *        with global vertex number j.  The "idj" flag is a 32-bit integer 
 *        "name" for vertex j, "procj" is a 32-bit color or processor number 
 *        associated with vertex j, and "infoj" is a 32-bit integer containing
 *        information about vertex j.
 *
 *        int edge* = List of all edges by vertex number in the form:
 *
 *        edge[ L * (2 + 2*(D+1)) ] = {
 *            { id1, proc1, info1, e1_v1, e1_v2 },
 *            { id2, proc2, info2, e2_v1, e2_v2 },
 *                            ...
 *            { idL, procL, infoL, eL_v1, eL_v2 }
 *        }
 *
 *        Here, ej_vk is a 32-bit integer giving the global vertex number
 *        making up the kth vertex of the edge with global number j.  The
 *        "idj" flag is a 32-bit integer "name" for edge j, "procj" is a
 *        32-bit color or processor number associated with edge j,
 *        "infoj" is a 32-bit integer containing the information about the
 *        edge j.
 *
 * Notes: To recover a simplex (triangle) mesh from the vertices and edges,
 *        we must traverse them and rebuild the simplex structure, and do it
 *        all in linear time.
 *
 *        The three-step simplex recovery algorithm is as follows:
 *
 *        (1) Definite vertices and edges from the input data.
 *
 *        (2) Create all possible simplices as follows:
 *
 *            For (vx=firstVV; vx!=lastVV; vx=nextVV) {
 *            | Build all simplices which use vx as a vertex as follows:
 *            |   2D CASE: For each distinct pair of vertices connected
 *            |            by an edge with vx, if this pair are also
 *            |            connected by an edge, then the three make a
 *            |            triangle.  If the triangle has not already
 *            |            been created, then do so and add to the
 *            |            simplex rings for each of the three vertices.
 *            |   3D CASE: For each distinct trio of vertices connected
 *            |            by an edge with vx, if this trio forms a triangle,
 *            |            then the foursome make a tetrahedron.  If the
 *            |            tetrahedron has not already been created, then do
 *            |            so and add to the simplex rings for each of the
 *            |            four vertices.
 *            EndFor
 *
 *        (3) Remove a few "bad" simplices which were created incorrectly
 *            in Step (2) above as follows:
 *
 *            For (sm=firstSS; sm!=lastSS; sm=nextSS) {
 *            | Remove all "bad" simplices, defined to be those satisfying:
 *            |   2D CASE: All interior edges have multiple nabors, and
 *            |            all boundary edges have at least one nabor.
 *            |   3D CASE: All interior faces have multiple nabors, and
 *            |            all boundary faces have at least one nabor.
 *            EndFor
 *
 *       NOTE: The crucial Step (3) of the above algorithm only works
 *       correctly if we correctly identify the boundary faces of the
 *       simplices, which is only possible if the input edge-based mesh
 *       was given with boundary edge information.  (We can construct
 *       correct face types from the types of the one or three edges
 *       forming the face in 2D or 3D respectively.)
 *
 *       We attempt to recover simplex face types from the given
 *       edge types.  However, (at least) one degenerate cases is not 
 *       recoverable: an isolated Neumann face surounded by Dirichlet faces
 *       will appear simply as three Dirichlet edges in the edge-based mesh,
 *       and we turn this into a Dirichlet face.  Note that if a Neumann
 *       face consists of more than one simplex face, then it will be
 *       recovered correctly from the edge types.  
 *
 *       Note also that if a non-simplex mesh is provided as input in
 *       as an edge-based mesh, we will build the edges as specified,
 *       but our attempt to build simplices will only recover those
 *       simplices which actually exist in the edge mesh.  Any other
 *       non-simplex polyedra will appear as holes in the final
 *       simplex mesh.
 * @endverbatim
 * @return  Success enumeration
 * @param   thee   Pointer to class Gem
 * @param   sock   socket for reading/writing a finite element mesh or mesh function
 */
VEXTERNC int Gem_readEdge(Gem *thee, Vio *sock);

/**
 * @ingroup Gem
 * @brief   Write out the edges of a 2- or 3-simplex mesh 
 *          in "MCEF" (MC-Edge-Format).
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemio.c)
 * @return  None
 * @param   thee   Pointer to class Gem
 * @param   sock   socket for reading/writing a finite element mesh or mesh function
 */
VEXTERNC void Gem_writeEdge(Gem *thee, Vio *sock);

/**
 * @ingroup Gem
 * @brief   Build a vertex-simplex mesh from a vertex-edge mesh.
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemio.c)\n
 * @verbatim
 * Notes:    Below is an email I sent to R. Bank outlining the idea of the
 *           algorithm.  I wanted to do this completely topologically, using
 *           no geometry information, so that it would work for abstract
 *           simplex manifolds.
 *
 * ---------------------------------------------------------------------------
 * Randy,
 * 
 * I realized on the drive that your 2-manifold example is not actually
 * a problem after all.
 * 
 * Consider first the planar situation we were worrying about, e.g., take a tet
 * and push the fourth vertex down into the plane of the other three vertices,
 * and you have the three (good) triangles inside one (bad) triangle.
 * As we agreed, any situation like this in the plane is managable (e.g. we can
 * detect the "bad" triangle) as long as we known what the boundary edges are
 * for the entire mesh.
 * 
 * I conjecture that the "bad" simplices in this case and other 2D cases, as
 * well as the 3D case, are precisely those whose faces (edges in 2D) satisfy
 * both of the following conditions.  (If this is a theorem, then everything
 * in this note is rigorous.)
 *    (a) All interior faces are shared with >1 other simplex
 *    (b) All boundary faces are shared with >0 other simplex
 * 
 * Consider now your example, e.g. the four surface triangles of a tetrahedron
 * as a 2-manifold without boundary.  We would like to be able to recover all
 * four surface triangles from the six edges, and we don't want to toss out one
 * of the triangles as we did above.  The key difference is that above, the
 * problem triangle either has one or more boundary edges, OR it is imbedded in
 * the interior of a mesh, so its interior edges have neighoring triangles.
 * 
 * That doesn't happen for the tet surface example.  We would first build all
 * possible triangles from triples of edges.  According to the above
 * "bad simplex" rule, all four of the triangles are "good", since they all
 * have exactly one naboring triangle sharing each edge.  So we get the correct
 * triangulation of the surface of the tet.
 * 
 * I conjecture that the following algorithm will rebuild d-simplices
 * (d=2 or d=3) from edges, using only topological information, with no
 * geometry information (and no floating point arithmetic at all, for that
 * matter).  The edge-based input mesh has to satisfy three properties:
 * 
 *    1. The dimension "d" of the mesh is specified (either d=2 or d=3)
 *    2. The edge-based mesh was built from an underlying d-simplex mesh
 *    3. The boundary edges are marked as such
 * 
 * The two-step simplex reconstruction algorithm is then as follows:
 * 
 *    1. For each vertex "v" do:
 *          For all vertices connected to "v" by a single edge do:
 *             Build every possible d-simplex
 *          EndFor
 *       EndFor
 *    2. For each simplex "s" that was built in step 1 do:
 *          3.  If d=3, calculate all "face" types of "s" based on edge types
 *          4.  Remove simplex "s" from list of simplices if (a) AND (b) hold:
 *              (a) ALL inter faces (edges if d=2) shared by >1 other simplices
 *              (b) ALL bndry faces (edges if d=2) shared by >0 other simplices
 *       EndFor
 * 
 * If you can come up with a triangulation of any 2-manifold, with or without
 * boundary, with or without holes, etc, which can't be turned into an edge
 * mesh and then recovered correctly using the above algorithm, then I'll
 * give you a dollar...  -mike
 * @endverbatim
 * @return  None
 * @param   thee   Pointer to class Gem
 */
VEXTERNC void Gem_buildSfromE(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Write out the boundary edges or boundary faces of a 
 *          2-simplex or 3-simplex mesh in a "BREP" format for input
 *          into Barry Joe's 2D/3D Geompak.
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemio.c)
 * @return  None
 * @param   thee   Pointer to class Gem
 * @param   sock   socket for reading/writing a finite element mesh or mesh function
 */
VEXTERNC void Gem_writeBrep(Gem *thee, Vio *sock);

/**
 * @ingroup Gem
 * @brief   Write out boundary edges of a 2-simplex mesh 
 *          in a "BREP" format for input into Barry Joe's 2D Geompak.
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemio.c)
 * @return  None
 * @param   thee   Pointer to class Gem
 * @param   sock   socket for reading/writing a finite element mesh or mesh function
 */
VEXTERNC void Gem_writeBrep2(Gem *thee, Vio *sock);

/**
 * @ingroup Gem
 * @brief   Write out boundary faces of a 3-simplex mesh
 *          in a "BREP" format for input into Barry Joe's 3D Geompak.
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemio.c)
 * @return  None
 * @param   thee   Pointer to class Gem
 * @param   sock   socket for reading/writing a finite element mesh or mesh function
 */
VEXTERNC void Gem_writeBrep3(Gem *thee, Vio *sock);

/**
 * @ingroup Gem
 * @brief   Write out triangles of a 2-manifold as a 3-simplex boundary mesh
 *          in a "BREP" format for input into Barry Joe's 3D Geompak.
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemio.c)
 * @return  None
 * @param   thee   Pointer to class Gem
 * @param   sock   socket for reading/writing a finite element mesh or mesh function
 */
VEXTERNC void Gem_writeBrep2to3(Gem *thee, Vio *sock);

/**
 * @ingroup Gem
 * @brief   Write a finite element mesh or mesh function in some format. 
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemdisp.c)
 * @return  None
 * @param   thee      Pointer to class Gem
 * @param   sock      socket for reading/writing a finite element mesh or mesh function
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
 * @param   defX      Pointer to vertex deformation or displacement values
 * @param   format    Pointer to GV/MATH format
 */
VEXTERNC void Gem_writeGEOM(Gem *thee, Vio *sock,
    int defKey, int colKey, int chartKey, double gluVal, int fkey,
    double *defX[MAXV], char *format);

/**
 * @ingroup Gem
 * @brief   Write a finite element mesh or mesh function in some format.  
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemdisp.c)
 * @return  None
 * @param   thee      Pointer to class Gem
 * @param   sock      socket for reading/writing a finite element mesh or mesh function
 * @param   fldKey    index for drawing the mesh
 *                    fldKey == 0  ==> draw mesh as it is
 *                    fldKey == 1  ==> write field[] as a scalar field over the mesh
 *                    fldKey == 2  ==> write field[] as a 2-vector field over the mesh
 *                    fldKey == 3  ==> write field[] as a 3-vector field over the mesh
 *                    fldKey == 4  ==> etc
 * @param   defX      Pointer to vertex deformation or displacement values
 * @param   format    Pointer to GV/MATH format
 */
VEXTERNC void Gem_writeSOL(Gem *thee, Vio *sock,
    int fldKey,
    double *defX[MAXV], char *format);

/**
 * @ingroup Gem
 * @brief   Produce an OFF file header for a volume mesh.
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemdisp.c)
 * @return  None
 * @param   thee      Pointer to class Gem
 * @param   sock      socket for reading/writing a finite element mesh or mesh function
 */
VEXTERNC void Gem_writeVolHeaderOFF(Gem *thee, Vio *sock);

/**
 * @ingroup Gem
 * @brief   Produce an OFF file header for a boundary mesh.  
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemdisp.c)
 * @return  None
 * @param   thee      Pointer to class Gem
 * @param   sock      socket for reading/writing a finite element mesh or mesh function
 */
VEXTERNC void Gem_writeBndHeaderOFF(Gem *thee, Vio *sock);

/**
 * @ingroup Gem
 * @brief   Produce an OFF file trailer. 
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemdisp.c)
 * @return  None
 * @param   thee      Pointer to class Gem
 * @param   sock      socket for reading/writing a finite element mesh or mesh function
 */
VEXTERNC void Gem_writeTrailerOFF(Gem *thee, Vio *sock);

/**
 * @ingroup Gem
 * @brief   Write out a mesh in "Geomview OFF" format to a file or socket.  
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemdisp.c)
 * @return  None
 * @param   thee      Pointer to class Gem
 * @param   sock      socket for reading/writing a finite element mesh or mesh function
 * @param   defKey    0  ==> draw mesh as it is\n
 *                    1  ==> use "def??" as new vertex coords (deformation)\n
 *                    2  ==> add "def??" to old vertex coords (displacement)
 * @param   colKey    0  ==> color simplices all same default color\n
 *                    1  ==> color simplices based on their chart\n 
 *                    2  ==> color boundary simplices based on type
 * @param   chartKey  <  0  ==> draw all simplices\n
 *                    >= 0  ==> draw only simplices with chart chartKey  
 * @param   gluVal    gluVal == 1. ==> draw all simplices glued together\n
 *                    0. < gluVal <  1. ==> draw simplices with some separation
 * @param   fkey      0  ==> draw simplices\n
 *                    1  ==> draw only simplex boundary faces\n
 *                    2  ==> draw only simplices with a boundary face
 * @param   defX      defX[][MAXV] ==> vertex deformation or displacement values
 */
VEXTERNC void Gem_writeGV(Gem *thee, Vio *sock,
    int defKey, int colKey, int chartKey, double gluVal, int fkey,
    double *defX[MAXV]);

/**
 * @ingroup Gem
 * @brief   Write out the faces of a 2-simplex mesh as a complete and legal
 *          1-simplex mesh in "Geomview SKEL" format. 
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemdisp.c)
 * @return  None
 * @param   thee      Pointer to class Gem
 * @param   sock      socket for reading/writing a finite element mesh or mesh function
 * @param   defKey    0  ==> draw mesh as it is\n
 *                    1  ==> use "def??" as new vertex coords (deformation)\n
 *                    2  ==> add "def??" to old vertex coords (displacement)
 * @param   colKey    0  ==> color simplices all same default color\n
 *                    1  ==> color simplices based on their chart\n 
 *                    2  ==> color boundary simplices based on type
 * @param   chartKey  <  0  ==> draw all simplices\n
 *                    >= 0  ==> draw only simplices with chart chartKey  
 * @param   gluVal    gluVal == 1. ==> draw all simplices glued together\n
 *                    0. < gluVal <  1. ==> draw simplices with some separation
 * @param   defX      defX[][MAXV] ==> vertex deformation or displacement values
 */
VEXTERNC void Gem_writeFace2dGV(Gem *thee, Vio *sock,
    int defKey, int colKey, int chartKey, double gluVal,
    double *defX[MAXV]);

/**
 * @ingroup Gem
 * @brief   Write out the faces of a 3-simplex mesh as a complete and legal
 *          2-simplex mesh in "Geomview OFF" format.   
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemdisp.c)
 * @return  None
 * @param   thee      Pointer to class Gem
 * @param   sock      socket for reading/writing a finite element mesh or mesh function
 * @param   defKey    0  ==> draw mesh as it is\n
 *                    1  ==> use "def??" as new vertex coords (deformation)\n
 *                    2  ==> add "def??" to old vertex coords (displacement)
 * @param   colKey    0  ==> color simplices all same default color\n
 *                    1  ==> color simplices based on their chart\n 
 *                    2  ==> color boundary simplices based on type
 * @param   chartKey  <  0  ==> draw all simplices\n
 *                    >= 0  ==> draw only simplices with chart chartKey  
 * @param   gluVal    gluVal == 1. ==> draw all simplices glued together\n
 *                    0. < gluVal <  1. ==> draw simplices with some separation
 * @param   defX      defX[][MAXV] ==> vertex deformation or displacement values
 */
VEXTERNC void Gem_writeFace3dGV(Gem *thee, Vio *sock,
    int defKey, int colKey, int chartKey, double gluVal,
    double *defX[MAXV]);

/**
 * @ingroup Gem
 * @brief   Produce an MATH file header for a volume mesh.     
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemdisp.c)
 * @return  None
 * @param   thee      Pointer to class Gem
 * @param   sock      socket for reading/writing a finite element mesh or mesh function
 */
VEXTERNC void Gem_writeHeaderMATH(Gem *thee, Vio *sock);

/**
 * @ingroup Gem
 * @brief   Produce an MATH file trailer.  
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemdisp.c)
 * @return  None
 * @param   thee      Pointer to class Gem
 * @param   sock      socket for reading/writing a finite element mesh or mesh function
 */
VEXTERNC void Gem_writeTrailerMATH(Gem *thee, Vio *sock);

/**
 * @ingroup Gem
 * @brief   Write out a mesh in "Mathematica" format to a file or socket.
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemdisp.c)
 * @return  None
 * @param   thee      Pointer to class Gem
 * @param   sock      socket for reading/writing a finite element mesh or mesh function
 * @param   defKey    0  ==> draw mesh as it is\n
 *                    1  ==> use "def??" as new vertex coords (deformation)\n
 *                    2  ==> add "def??" to old vertex coords (displacement)
 * @param   colKey    0  ==> color simplices all same default color\n
 *                    1  ==> color simplices based on their chart\n 
 *                    2  ==> color boundary simplices based on type
 * @param   chartKey  <  0  ==> draw all simplices\n
 *                    >= 0  ==> draw only simplices with chart chartKey  
 * @param   gluVal    gluVal == 1. ==> draw all simplices glued together\n
 *                    0. < gluVal <  1. ==> draw simplices with some separation
 * @param   fkey      0  ==> draw simplices\n
 *                    1  ==> draw only simplex boundary faces\n
 *                    2  ==> draw only simplices with a boundary face 
 * @param   defX      defX[][MAXV] ==> vertex deformation or displacement values
 */
VEXTERNC void Gem_writeMATH(Gem *thee, Vio *sock,
    int defKey, int colKey, int chartKey, double gluVal, int fkey,
    double *defX[MAXV]);

/**
 * @ingroup Gem
 * @brief   Write out a scalar or vector function over a simplex mesh
 *          in "GMV" (General Mesh Viewer) format.  
 * @author  Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemdisp.c)
 * @return  None
 * @param   thee      Pointer to class Gem
 * @param   sock      socket for reading/writing a finite element mesh or mesh function
 * @param   fldKey    0  ==> draw mesh as it is \n
 *                    1  ==> write field[] as a scalar field over the mesh \n
 *                    2  ==> write field[] as a 2-vector field over the mesh\n
 *                    3  ==> write field[] as a 3-vector field over the mesh\n
 *                    4  ==> etc
 * @param   defX      defX[][MAXV] ==> possible scalar or vector field
 */
VEXTERNC void Gem_writeGMV(Gem *thee, Vio *sock,
    int fldKey, double *defX[MAXV]);

/**
 * @ingroup Gem
 * @brief   Write out a scalar or vector function over a simplex mesh  
 *          in "UCD" (Unstructured Cell Data) format for AVS 5.0
 * @authors Amit Majumdar and Stephen Bond 
 *          (created using Mike Holst's Gem_writeGMV as template)
 * @note    Class Gem: Non-inlineable methods (gemdisp.c)\n
 *          Format:   Our particular use of the UCD format is as follows\n
 *          \n
 *          PART 1.1: [# NODES] [# CELLS] [DIM NODEDAT] [DIM CELLDAT] [0]\n
 *          PART 1.2: [NODE ID] [X COORD] [Y COORD] [Z COORD]\n
 *                    (REPEATED FOR EACH NODE)\n
 *          PART 1.3: [CELL ID] [MATERIAL] [CELL TYPE] [LIST OF NODE IDS]\n
 *                    (REPEATED FOR EACH CELL)\n     
 *          PART 2.1: [NUM NODE FIELD COMPONENTS] [LIST OF COMPONENT SIZES]\n
 *          PART 2.2: [COMPONENT LABEL], [COMPONENT UNITS]\n
 *                    (REPEATED FOR EACH NODE DATA COMPONENT)\n
 *          PART 2.3: [NODE ID] [DATA VALUES]\n
 *                    (REPEATED FOR EACH NODE)\n
 *          \n
 *          PART 3.1: [NUM CELL FIELD COMPONENTS] [LIST OF COMPONENT SIZES]\n
 *          PART 3.2: [COMPONENT LABEL], [COMPONENT UNITS]\n
 *                    (REPEATED FOR EACH CELL DATA COMPONENT)\n
 *          PART 3.3: [CELL ID] [DATA VALUES]\n
 *                    (REPEATED FOR EACH CELL)
 * @return  None
 * @param   thee      Pointer to class Gem
 * @param   sock      socket for reading/writing a finite element mesh or mesh function
 * @param   fldKey    0  ==> draw mesh as it is \n
 *                    1  ==> write field[] as a scalar field over the mesh \n
 *                    2  ==> write field[] as a 2-vector field over the mesh\n
 *                    3  ==> write field[] as a 3-vector field over the mesh\n
 *                    4  ==> etc
 * @param   defX      defX[][MAXV] ==> possible scalar or vector field
 */
VEXTERNC void Gem_writeUCD(Gem *thee, Vio *sock,
    int fldKey, double *defX[MAXV]);


/**
 * @ingroup Gem
 * @brief   Write out a scalar or vector function over a simplex mesh 
 *          in "DX" (www.opendx.org) format.  
 * @authors Nathan Baker, Stephen Bond, and Michael Holst  
 *          (created using Mike Holst's Gem_writeGMV as template)
 * @note    Class Gem: Non-inlineable methods (gemdisp.c)\n
 * @return  None
 * @param   thee      Pointer to class Gem
 * @param   sock      socket for reading/writing a finite element mesh or mesh function
 * @param   fldKey    0  ==> draw mesh as it is \n
 *                    1  ==> write field[] as a scalar field over the mesh \n
 *                    2  ==> write field[] as a 2-vector field over the mesh\n
 *                    3  ==> write field[] as a 3-vector field over the mesh\n
 *                    4  ==> etc
 * @param   defX      defX[][MAXV] ==> possible scalar or vector field
 */
VEXTERNC void Gem_writeDX(Gem *thee, Vio *sock,
    int fldKey, double *defX[MAXV]);

/**
 * @ingroup Gem
 * @brief   Write out a scalar or vector function over a simplex mesh 
 *          in "TEC" (www.opendx.org) format. 
 * @authors Nathan Baker, Jason Suen, and Michael Holst  
 *          (created using Mike Holst's Gem_writeGMV as template)
 * @note    Class Gem: Non-inlineable methods (gemdisp.c)\n
 * @return  None
 * @param   thee      Pointer to class Gem
 * @param   sock      socket for reading/writing a finite element mesh or mesh function
 * @param   fldKey    0  ==> draw mesh as it is \n
 *                    1  ==> write field[] as a scalar field over the mesh \n
 *                    2  ==> write field[] as a 2-vector field over the mesh\n
 *                    3  ==> write field[] as a 3-vector field over the mesh\n
 *                    4  ==> etc
 * @param   defX      defX[][MAXV] ==> possible scalar or vector field
 */
VEXTERNC void Gem_writeTEC(Gem *thee, Vio *sock,
    int fldKey, double *defX[MAXV]);

/*
 * ***************************************************************************
 * Class Gem: Non-inlineable methods (gemext.c)
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * Class Gem: Non-inlineable methods (gemcube.c)
 * ***************************************************************************
 */

/**
 * @ingroup Gem
 * @brief   Generate a unit cube domain. 
 * @authors Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemcube.c)\n
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_makeCube(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Generate a unit octahedron domain.  
 * @authors Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemcube.c)\n
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_makeOctahedron(Gem *thee);

/**
 * @ingroup Gem
 * @brief   Generate a unit icosahedron domain.   
 * @authors Michael Holst
 * @note    Class Gem: Non-inlineable methods (gemcube.c)\n
 * @return  None
 * @param   thee  Pointer to class Gem
 */
VEXTERNC void Gem_makeIcosahedron(Gem *thee);

#endif /* _GEM_H_ */

