/**
 * @defgroup MCsh MCsh class 
 * @brief    the foundation MC context class
 */

/**
 *  @file       mcsh.h
 *  @ingroup    MCsh
 *  @brief      Class MCsh: the foundation MC context class.
 *  @author     Michael Holst
 *  @note       None
 *  @version    $Id: mcsh.h,v 1.16 2010/08/12 05:19:11 fetk Exp $ 
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

#ifndef _MCSH_H_
#define _MCSH_H_

#include <mc/mc_base.h>

#include <mc/nam.h>

/*
 * ***************************************************************************
 * Class MCsh: Parameters and datatypes
 * ***************************************************************************
 */

/** 
 * @ingroup MCsh
 * @brief   Class MCsh: Definition
 * @author  Michael Holst
 */
struct sMCsh {

    Vmem *vmem;     /**< @brief the memory manager                           */
                                                                              
    Gem  *gm;       /**< @brief the geometry machine object                  */
    PDE  *pde;      /**< @brief the differential equation object             */
    Aprx *aprx;     /**< @brief the approximation object                     */
    AM   *am;       /**< @brief the multilevel linear algebra manager object */

    Vsh  *vsh;      /**< @brief the environment and shell object             */

    char PR[80];    /**< @brief shell prompt                                 */

    /** @brief user defined shell object */
    int (*USER_shell)(void *thee, int argc, char **argv);

};

/**
 * @ingroup MCsh
 * @brief   Declaration of the MCsh class as the MCsh structure
 * @author  Michael Holst
 */
typedef struct sMCsh MCsh;

/*
 * ***************************************************************************
 * Class MCsh: Inlineable methods (mcsh.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_MCSH)
#else /* if defined(VINLINE_MCSH) */
#endif /* if !defined(VINLINE_MCSH) */

/*
 * ***************************************************************************
 * Class MCsh: Non-inlineable methods (mcsh.c)
 * ***************************************************************************
 */

/**
 * @ingroup MCsh
 * @brief   The MCsh constructor.
 * @author  Michael Holst
 * @note    Class MCsh: Non-inlineable methods (mcsh.c) 
 * @return  Pointer to a new allocated MCsh shell
 * @param   tpde Pointer to the PDE object
 * @param   argc number of the command line arguments
 * @param   argv the command line arguments
 */
VEXTERNC MCsh* MCsh_ctor(PDE *tpde, int argc, char **argv);

/**
 * @ingroup MCsh
 * @brief   The MCsh destructor.
 * @author  Michael Holst
 * @note    Class MCsh: Non-inlineable methods (mcsh.c) 
 * @return  None
 * @param   thee Pointer to a MCsh shell
 */
VEXTERNC void MCsh_dtor(MCsh **thee);

/**
 * @ingroup MCsh
 * @brief   The actuall shell
 * @author  Michael Holst
 * @note    Class MCsh: Non-inlineable methods (mcsh.c) 
 * @return  Success enumeration
 * @param   thee       Pointer to a MCsh shell
 * @param   USER_shell Pointer to a user-defined shell
 */
VEXTERNC int MCsh_shell(MCsh *thee,
    int (*USER_shell)(void *thee, int argc, char **argv));

/**
 * @ingroup MCsh
 * @brief   The actuall shell (version with parallel extensions).
 * @author  Michael Holst
 * @note    Class MCsh: Non-inlineable methods (mcsh.c) 
 * @return  Success enumeration
 * @param   thee       Pointer to a MCsh shell
 * @param   USER_shell Pointer to a user-defined shell
 */
VEXTERNC int MCsh_pshell(MCsh *thee,
    int (*USER_shell)(void *thee, int argc, char **argv));

/**
 * @ingroup MCsh
 * @brief   Print the exact current malloc usage.   
 * @author  Michael Holst
 * @note    Class MCsh: Non-inlineable methods (mcsh.c) 
 * @return  None
 * @param   thee Pointer to a MCsh shell
 */
VEXTERNC void MCsh_memChk(MCsh *thee);

#endif /* _MCSH_H_ */
