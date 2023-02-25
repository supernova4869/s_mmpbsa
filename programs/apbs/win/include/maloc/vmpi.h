/**
 * @defgroup Vmpi Vmpi class
 * @brief    A Virtual MPI communication lawyer object.
 */

/**
 *  @file       vmpi.h
 *  @ingroup    Vmpi
 *  @brief      Class Vmpi: a Virtual MPI communication layer object.
 *  @author     Michael Holst
 *  @note       None
 *  @version    $Id: vmpi.h,v 1.29 2010/08/12 05:40:23 fetk Exp $
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


#ifndef _VMPI_H_
#define _VMPI_H_

#include <maloc/maloc_base.h>

#include <maloc/vsys.h>

/*
 * ***************************************************************************
 * Class Vmpi: Parameters and datatypes
 * ***************************************************************************
 */


/** 
 * @ingroup Vmpi
 * @author  Michael Holst
 * @brief   Class Vmpi: Definition 
 */
struct sVmpi {
    int  mpi_rank;     /**< my process ID                                   */
    int  mpi_size;     /**< number of processess in this execution           */
};

/**
 * @brief   Declaration of the Vmpi class as the Vmpi structure
 * @ingroup Vmpi
 * @author  Michael Holst
 */
typedef struct sVmpi Vmpi;

/*
 * ***************************************************************************
 * Class Vmpi: Inlineable methods (vmpi.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_MALOC)
#else /* if defined(VINLINE_MALOC) */
#endif /* if !defined(VINLINE_MALOC) */

/**
 * @ingroup Vmpi
 * @brief   The Vmp initializer. 
 * @author  Michael Holst
 * @note    Class Vmpi: Non-inlineable method (vmpi.c) 
 * @return  Success enumeration
 * @param   argc number of the command line arguments
 * @param   argv the command line arguments
 */
VEXTERNC int Vmpi_init(int *argc, char ***argv);

/**
 * @ingroup Vmpi
 * @brief   The Vmp finalizer. 
 * @author  Michael Holst
 * @note    Class Vmpi: Non-inlineable method (vmpi.c) 
 * @return  Success enumeration
 */
VEXTERNC int Vmpi_finalize(void);

/**
 * @ingroup Vmpi
 * @brief   The Vmpi constructor. 
 * @author  Michael Holst
 * @note    Class Vmpi: Non-inlineable method (vmpi.c) 
 * @return  Success enumeration
 */
VEXTERNC Vmpi* Vmpi_ctor(void);

/**
 * @ingroup Vmpi
 * @brief   The Vmpi destructor.
 * @author  Michael Holst
 * @note    Class Vmpi: Non-inlineable method (vmpi.c) 
 * @return  None
 * @param   thee Pointer to pointer of the Vmpi object
 */
VEXTERNC void Vmpi_dtor(Vmpi **thee);

/**
 * @ingroup Vmpi
 * @brief   Return my processor ID.    
 * @author  Michael Holst
 * @note    Class Vmpi: Non-inlineable method (vmpi.c) 
 * @return  Success enumeration
 * @param   thee Pointer to the Vmpi object
 */
VEXTERNC int Vmpi_rank(Vmpi *thee);

/**
 * @ingroup Vmpi
 * @brief   Return the number of processors involved.
 * @author  Michael Holst
 * @note    Class Vmpi: Non-inlineable method (vmpi.c) 
 * @return  Success enumeration
 * @param   thee Pointer to the Vmpi object
 */
VEXTERNC int Vmpi_size(Vmpi *thee);

/**
 * @ingroup Vmpi
 * @brief   An MPI barrier.   
 * @author  Michael Holst
 * @note    Class Vmpi: Non-inlineable method (vmpi.c) 
 * @return  Success enumeration
 * @param   thee Pointer to the Vmpi object
 */
VEXTERNC int Vmpi_barr(Vmpi *thee);

/**
 * @ingroup Vmpi
 * @brief   An MPI blocking send. 
 * @author  Michael Holst
 * @note    Class Vmpi: Non-inlineable method (vmpi.c) 
 * @return  Success enumeration
 * @param   thee    Pointer to the Vmpi object
 * @param   des     rank of receiving processor 
 * @param   buf     buffer containing message 
 * @param   bufsize number of items (of declared type) in buffer 
 */
VEXTERNC int Vmpi_send(Vmpi *thee, int des, char *buf, int bufsize);

/**
 * @ingroup Vmpi
 * @brief   An MPI blocking receive.
 * @author  Michael Holst
 * @note    Class Vmpi: Non-inlineable method (vmpi.c) 
 * @return  Success enumeration
 * @param   thee    Pointer to the Vmpi object
 * @param   src     rank of receiving processor 
 * @param   buf     buffer containing message 
 * @param   bufsize number of items (of declared type) in buffer 
 */
VEXTERNC int Vmpi_recv(Vmpi *thee, int src, char *buf, int bufsize);

/**
 * @ingroup Vmpi
 * @brief   An MPI broadcast. 
 * @author  Michael Holst
 * @note    Class Vmpi: Non-inlineable method (vmpi.c) 
 * @return  Success enumeration
 * @param   thee    Pointer to the Vmpi object
 * @param   buf     buffer containing message 
 * @param   bufsize number of items (of declared type) in buffer 
 */
VEXTERNC int Vmpi_bcast(Vmpi *thee, char *buf, int bufsize);

/**
 * @ingroup Vmpi
 * @brief   An MPI reduce.       
 * @author  Michael Holst
 * @note    Class Vmpi: Non-inlineable method (vmpi.c) 
 * @return  Success enumeration
 * @param   thee     Pointer to the Vmpi object
 * @param   sbuf     address of send buffer (choice)
 * @param   rbuf     address of receiving buffer (choice)
 * @param   bufsize  number of items (of declared type) in buffer
 */
VEXTERNC int Vmpi_reduce(Vmpi *thee, char *sbuf, char *rbuf, int bufsize);

/**
 * @ingroup Vmpi
 * @brief   An MPI non-blocking send.
 * @author  Michael Holst
 * @note    Class Vmpi: Non-inlineable method (vmpi.c) 
 * @return  Success enumeration
 * @param   thee    Pointer to the Vmpi object
 * @param   des     rank of receiving processor 
 * @param   buf     buffer containing message 
 * @param   bufsize number of items (of declared type) in buffer 
 */
VEXTERNC int Vmpi_isend(Vmpi *thee, int des, char *buf, int bufsize);

#endif /* _VMPI_H_ */

