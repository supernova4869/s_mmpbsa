/**
 * @defgroup Node Node class
 * @brief    Class Node: a node object.
 */

/**
 *  @file       node.h
 *  @brief      Class Node: a node object.
 *  @author     Michael Holst
 *  @note       None
 *  @version    $Id: node.h,v 1.11 2010/08/12 05:18:22 fetk Exp $ 
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

#ifndef _NODE_H_
#define _NODE_H_

#include <mc/mc_base.h>

#include <mc/bam.h>

/**
 * @ingroup global_mc
 * @author  Michael Holst
 * @brief Class Node: Parameters and datatypes 
 */
typedef struct aNode {
  /** @brief My id */
    int        myid; 
  /** @brief My vertex association (-1=none,else k>=0) */
    int        myvert;
  /** @brief My node type (0=interior,odd=diri,even=neu) */
    char       type;         
  /** @brief My chart number */
    char       chart;        
  /** @brief number of coordinates for my chart */
    char       numx;        
  /** @brief  My nodal coordinates (wrt my chart) */ 
    double     x[3]; 
  /** @brief My nodal value */       
    double     val;  
} aNode;

/** 
 * @ingroup Node
 * @brief   Class Node: Definition 
 * @author  Michael Holst
 */
struct sNode {

    Vmem   *vmem;        /**< @brief the memory manager                  */
    int    iMadeVmem;    /**< @brief did i make vmem or was it inherited */

    int    numR;         /**< @brief the number of nodes                 */
    aNode  *data;        /**< @brief vector data itself                  */

};

/**
 * @ingroup Node
 * @brief   Declaration of the Node class as the Node structure 
 * @author  Michael Holst
 */

typedef struct sNode Node;


#if !defined(VINLINE_APRX)

    /** 
     * @ingroup Node
     * @brief   Return the number of nodes. 
     * @author  Michael Holst
     * @note    Class Node: Inlineable methods (node.c) 
     * @return  the number of nodes
     * @param   thee Pointer to a Node allocated memory location
     */
    VEXTERNC int Node_numR(Node *thee);

    /**
     * @ingroup Node
     * @brief   Return a pointer to the node data.            
     * @author  Michael Holst
     * @note    Class Node: Inlineable methods (node.c) 
     * @return  a pointer to the node data.
     * @param   thee Pointer to a Node allocated memory location
     */
    VEXTERNC aNode* Node_data(Node *thee);
#else /* if defined(VINLINE_APRX) */
    /** 
     * @ingroup Node
     * @brief   Class Node: Inlineable methods (node.c) if defined(VINLINE_APRX) 
     * @author  Michael Holst
     * @note    Class Node: Inlineable methods (node.c) 
     * @return  the number of nodes
     * @param   thee Pointer to a Node allocated memory location
     */
#   define Node_numR(thee)            ((thee)->numR)
    /**
     * @ingroup Node
     * @brief   Class Node: Inlineable methods (node.c) if defined(VINLINE_APRX) 
     * @author  Michael Holst
     * @note    Class Node: Inlineable methods (node.c) 
     * @return  a pointer to the node data.
     * @param   thee Pointer to a Node allocated memory location
     */
#   define Node_data(thee)            ((thee)->data)
#endif /* if !defined(VINLINE_APRX) */

/** 
 * @ingroup Node
 * @brief   The node constructor.  
 * @author  Michael Holst
 * @note    Class Node: Non-inlineable methods (node.c) 
 * @return  Pointer to a Node allocated memory location
 * @param   vmem    Memory management object
 * @param   pnumR   num of rows in the matrix
 */
VEXTERNC Node* Node_ctor(Vmem *vmem, int pnumR);

/** 
 * @ingroup Node
 * @brief   The node destructor.  
 * @author  Michael Holst
 * @note    Class Node: Non-inlineable methods (node.c) 
 * @param   thee   Pointer to a Node allocated memory location
 */
VEXTERNC void Node_dtor(Node **thee);

/** 
 * @ingroup Node
 * @brief   Print the exact current malloc usage.
 * @author  Michael Holst
 * @note    Class Node: Non-inlineable methods (node.c) 
 * @return  none
 * @param   thee   Pointer to a Node allocated memory location
 */
VEXTERNC void Node_memChk(Node *thee);

#endif /* _NODE_H_ */

