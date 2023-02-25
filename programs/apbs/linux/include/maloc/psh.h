/**
 *  @file       psh.h
 *  @ingroup    Vsh
 *  @brief      Header file for a simple parallel extension of ALOC's VSH.
 *  @author     Michael Holst
 *  @note       None
 *  @version    $Id: psh.h,v 1.28 2010/08/12 05:40:23 fetk Exp $
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

#ifndef _PSH_H_
#define _PSH_H_

#include <maloc/maloc_base.h>

#include <maloc/vsys.h>
#include <maloc/vsh.h>
#include <maloc/vmp.h>

/**
 * @ingroup Vsh
 * @brief   Drop-in replacement for Vsh_shell giving parallel extensions.
 * @author  Michael Holst
 * @note    None
 * @return  Success enumeration
 * @param   thee    Pointer to the Vsh object
 * @param   pPR     minimal prompt
 * @param   pthee   the externally supplied builtin object pointer
 * @param   builtin external supershell builtin function
 */
VEXTERNC int Vsh_pshell(Vsh *thee, char *pPR, void *pthee,
    int (*builtin)(void *thee, int argc, char **argv));

#endif /* _PSH_H_ */

