/**
 * @defgroup Bnode Bnode class
 * @brief    A block node object.
 */

/**
 *  @file       bnode.h
 *  @ingroup    Bnode
 *  @brief      Class Bnode: a block node object.
 *  @version    $Id: bnode.h,v 1.21 2010/08/12 05:18:21 fetk Exp $ 
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

#ifndef _BNODE_H_
#define _BNODE_H_

#include <mc/mc_base.h>

#include <mc/bam.h>

#include <mc/node.h>

/** 
 * @ingroup Bnode
 * @brief Class Bnode: Parameters and datatypes. Class Bnode: Definition.
 * @author  Michael Holst
 */
struct sBnode {

    /** @brief the memory manager                           */
    Vmem   *vmem;        
    /** @brief did i make vmem or was it inherited          */ 
    int    iMadeVmem;    
  
    /** @brief Character string name for this node vector   */
    char   name[10];     
    /** @brief num vector blocks                            */
    int    numB;        
    /** @brief the nodes in each block                      */ 
    Node   *ND[MAXV];   

};

/**
 * @ingroup Bnode
 * @brief   Declaration of the Bnode class as the Vio structure
 * @author  Michael Holst
 * @return  None
 */
typedef struct sBnode Bnode;


#if !defined(VINLINE_APRX)
   /**
   * @ingroup Bnode
   * @brief   Return the number of blocks.  
   * @author  Michael Holst
   * @note    Class Bnode: Non-inlineable methods (bnode.c)\n
   *          if VINLINE_APRX is not defined
   * @return  the number of blocks
   * @param   thee the Bnode class
   */
    VEXTERNC int Bnode_numB(Bnode *thee);
#else 
   /**
   * @ingroup Bnode
   * @brief   Return the number of blocks.
   * @author  Michael Holst
   * @note    Class Bnode: Non-inlineable methods (bnode.c)\n
   *          if VINLINE_APRX is defined
   * @return  the number of blocks
   * @param   thee Pointer to the Bnode class 
   */
#   define Bnode_numB(thee)              ((thee)->numB)
#endif 


/**
 * @ingroup Bnode
 * @brief   The Block node destructor.
 * @author  Michael Holst
 * @note    Class Bnode: Non-inlineable methods (bnode.c)
 * @return  Pointer to a bnode allocated memory location 
 * @param   vmem   Memory management object
 * @param   pnumB  number of blocks
 * @param   pnumR  number of rows
 */
VEXTERNC Bnode* Bnode_ctor(Vmem *vmem, int pnumB, int pnumR[MAXV]);

/**
 * @ingroup Bnode
 * @brief   The Block node destructor.
 * @author  Michael Holst
 * @note    Class Bnode: Non-inlineable methods (bnode.c)
 * @return  None
 * @param   thee  Pointer to a bnode allocated memory location
 */
VEXTERNC void Bnode_dtor(Bnode **thee);

/**
 * @ingroup Bnode
 * @brief   Return the total number of nodes.         
 * @author  Michael Holst
 * @note    Class Bnode: Non-inlineable methods (bnode.c)
 * @return  the total number of nodes
 * @param   thee  Pointer to a bnode allocated memory location
 */
VEXTERNC int Bnode_numRT(Bnode *thee);

/**
 * @ingroup Bnode
 * @brief   Return the number of nodes in a block.
 * @author  Michael Holst
 * @note    Class Bnode: Non-inlineable methods (bnode.c)
 * @return  the number of nodes in a block.
 * @param   thee  Pointer to a bnode allocated memory location
 * @param   i     index of blocks
 */
VEXTERNC int Bnode_numR(Bnode *thee, int i);

/**
 * @ingroup Bnode
 * @brief   Return the nodes in one block.            
 * @author  Michael Holst
 * @note    Class Bnode: Non-inlineable methods (bnode.c)
 * @return  the nodes in one block
 * @param   thee  Pointer to a bnode allocated memory location
 * @param   i     index of blocks
 */
VEXTERNC Node *Bnode_nodes(Bnode *thee, int i);

/**
 * @ingroup Bnode
 * @brief   Return the actual nodes in one block.  
 * @author  Michael Holst
 * @note    Class Bnode: Non-inlineable methods (bnode.c)
 * @return  the actual nodes in one block
 * @param   thee  Pointer to a bnode allocated memory location
 * @param   i     index of blocks
 */
VEXTERNC aNode *Bnode_data(Bnode *thee, int i);

/**
 * @ingroup Bnode
 * @brief   Print the exact current malloc usage. 
 * @author  Michael Holst
 * @note    Class Bnode: Non-inlineable methods (bnode.c)
 * @return  None
 * @param   thee  Pointer to a bnode allocated memory location
 */
VEXTERNC void Bnode_memChk(Bnode *thee);

#endif /* _BNODE_H_ */

