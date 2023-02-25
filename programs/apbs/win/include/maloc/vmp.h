/**
 * @defgroup Vmp Vmp class
 * @brief    A Virtual MPI communication layer object.
 */

/**
 *  @file       vmp.h
 *  @ingroup    Vmp
 *  @brief      Class Vmp: a Virtual MPI communication layer object.
 *  @author     Michael Holst
 *  @note       None
 *  @version    $Id: vmp.h,v 1.22 2010/08/12 05:40:23 fetk Exp $
 *  
 *  @attention
 *  @verbatim
 *
 * MALOC = < Minimal Abstraction Layer for Object-oriented C >
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

#ifndef _VMP_H_
#define _VMP_H_

#include <maloc/maloc_base.h>

#include <maloc/vsys.h>
#include <maloc/vmpi.h>
#include <maloc/vcom.h>

/*
 * ***************************************************************************
 * Class Vmp: Parameters and datatypes
 * ***************************************************************************
 */

/**
 * @ingroup Vmp
 * @brief   Contains public data members for Vmp class
 * @author  Michael Holst
 * @note    None
 */
struct sVmp {
    int  mpi_rank;     /**< my process ID                                   */
    int  mpi_size;     /**< number of processess in this execution          */
};

/**
 * @ingroup Vmp
 * @brief   Declaration of the Vmp class as teh Vmp structure 
 * @author  Michael Holst
 */
typedef struct sVmp Vmp;

/*
 * ***************************************************************************
 * Class Vmp: Inlineable methods (vmp.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_MALOC)
#else /* if defined(VINLINE_MALOC) */
#endif /* if !defined(VINLINE_MALOC) */


/**
 * @ingroup Vmp
 * @brief   The Vmp initializer. 
 * @author  Michael Holst
 * @note    Class Vmp: Non-inlineable method (vmp.c) 
 * @return  Success enumeration
 * @param   argc number of the command line arguments
 * @param   argv the command line arguments
 */
VEXTERNC int Vmp_init(int *argc, char ***argv);

/**
 * @ingroup Vmp
 * @brief   The Vmp finalizer. 
 * @author  Michael Holst
 * @note    Class Vmp: Non-inlineable method (vmp.c) 
 * @return  Success enumeration
 */
VEXTERNC int Vmp_finalize(void);

/**
 * @ingroup Vmp
 * @brief   The Vmp constructor. 
 * @author  Michael Holst
 * @note    Class Vmp: Non-inlineable method (vmp.c) 
 * @return  Success enumeration
 */
VEXTERNC Vmp* Vmp_ctor(void);

/**
 * @ingroup Vmp
 * @brief   The Vmp destructor. 
 * @author  Michael Holst
 * @note    Class Vmp: Non-inlineable method (vmp.c) 
 * @return  None
 * @param   thee Pointer to pointer of Vmp object
 */
VEXTERNC void Vmp_dtor(Vmp **thee);

/**
 * @ingroup Vmp
 * @brief   Return my processor ID.
 * @author  Michael Holst
 * @note    Class Vmp: Non-inlineable method (vmp.c) 
 * @return  Success enumeration
 * @param   thee Pointer to the Vmp object
 */
VEXTERNC int Vmp_rank(Vmp *thee);

/**
 * @ingroup Vmp
 * @brief   Return the number of processors involved.     
 * @author  Michael Holst
 * @note    Class Vmp: Non-inlineable method (vmp.c) 
 * @return  Success enumeration
 * @param   thee Pointer to the Vmp object
 */
VEXTERNC int Vmp_size(Vmp *thee);

/**
 * @ingroup Vmp
 * @brief   An MPI barrier.     
 * @author  Michael Holst
 * @note    Class Vmp: Non-inlineable method (vmp.c) 
 * @return  Success enumeration
 * @param   thee Pointer to the Vmp object
 */
VEXTERNC int Vmp_barr(Vmp *thee);

/**
 * @ingroup Vmp
 * @brief   An MPI blocking send.     
 * @author  Michael Holst
 * @note    Class Vmp: Non-inlineable method (vmp.c) 
 * @return  Success enumeration
 * @param   thee    Pointer to the Vmp object
 * @param   des     rank of receiving processor 
 * @param   buf     buffer containing message 
 * @param   bufsize number of items (of declared type) in buffer 
 */
VEXTERNC int Vmp_send(Vmp *thee, int des, char *buf, int bufsize);

/**
 * @ingroup Vmp
 * @brief   An MPI blocking receive.    
 * @author  Michael Holst
 * @note    Class Vmp: Non-inlineable method (vmp.c) 
 * @return  Success enumeration
 * @param   thee    Pointer to the Vmp object
 * @param   src     rank of receiving processor 
 * @param   buf     buffer containing message 
 * @param   bufsize number of items (of declared type) in buffer 
 */
VEXTERNC int Vmp_recv(Vmp *thee, int src, char *buf, int bufsize);

#endif /* _VMP_H_ */

