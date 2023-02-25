/**
 * @defgroup Vsh Vsh class
 * @brief A bourne-compatible shell.
 */

/**
 *  @file       vsh.h
 *  @ingroup    Vsh
 *  @brief      Header file for vsh, a bourne-compatible shell.
 *  @author     Michael Holst
 *  @note       None
 *  @version    $Id: vsh.h,v 1.30 2010/08/12 05:40:29 fetk Exp $
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


#ifndef _VSH_H_
#define _VSH_H_

#include <maloc/maloc_base.h>

#include <maloc/vsys.h>

/*
 * ***************************************************************************
 * Class Vsh: Parameters and datatypes
 * ***************************************************************************
 */

/**
 * @brief   Contains public data members for Vsh class
 * @author  Michael Holst
 * @ingroup Vsh
 */
struct sVsh {

    /** @brief the memory manager                            */
    Vmem   *vmem;        
    /** @brief did i make vmem or was it inherited           */
    int    iMadeVmem;    

    /** @brief whether the shell should process (argc,argv)  */
    char processArgs;  

    /** @brief number of environment variables               */
    int envValuLen;
    /** @brief number of environment variable help strings   */  
    int envInfoLen;
    /** @brief the environment variables                     */
    char **envValu; 
    /** @brief the environment variable help strings         */ 
    char **envInfo; 

    /** @brief input unit                                    */
    FILE *inUnit;  
    /** @brief script input unit                             */  
    FILE *scUnit;
    /** @brief input unit                                    */  
    FILE *clUnit;
    /** @brief input unit                                    */
    FILE *cinUnit; 
    /** @brief input unit                                    */
    char cinName[VMAX_ARGLEN];

    /** @brief minimal prompt (just the binary name)         */
    char PR[VMAX_ARGLEN];      
    /** @brief full prompt (user,hostname,path,etc)          */ 
    char PR_PATH[VMAX_ARGLEN]; 
    /** @brief the exit print string                         */
    char PR_EXIT[VMAX_ARGLEN];  

    /** @brief external supershell command key               */
    int cmdKey;       
    /** @brief external supershell object                    */   
    void *Ext_thee;   

    /** @brief internal buffer                               */
    char *buf;
    /** @brief internal buffer size                          */
    int bufsize;        

    /** @brief external supershell builtin function */
    int (*Ext_builtin)(void *thee, int argc, char **argv);

};

/**
 * @brief   Declaration of the Vsh class as the Vsh structure
 * @author  Michael Holst
 * @ingroup Vsh
 */
typedef struct sVsh Vsh;


/*
 * ***************************************************************************
 * Class Vsh: Inlineable methods (vsh.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_MALOC)
#else /* if defined(VINLINE_MALOC) */
#endif /* if !defined(VINLINE_MALOC) */


/**
 * @ingroup Vsh
 * @brief   Create the shell. 
 * @author  Michael Holst
 * @note    Class Vsh: Non-inlineable method (vsh.c) 
 * @return  Pointer to the new allocated Vsh object
 * @param   vmem Memory management object
 * @param   argc number of the command line arguments
 * @param   argv the command line arguments
 */ 
VEXTERNC Vsh* Vsh_ctor(Vmem *vmem, int argc, char **argv);

/**
 * @ingroup Vsh
 * @brief   Destroy the shell. 
 * @author  Michael Holst
 * @note    Class Vsh: Non-inlineable method (vsh.c) 
 * @return  None
 * @param   thee Pointer to the Vsh object
 */ 
VEXTERNC void Vsh_dtor(Vsh **thee);

/**
 * @ingroup Vsh
 * @brief   A bash-like shell with user-definable extensions.           
 * @author  Michael Holst
 * @note    Class Vsh: Non-inlineable method (vsh.c) 
 * @return  Success enumeration
 * @param   thee    Pointer to the Vsh object
 * @param   pPR     minimal prompt
 * @param   pthee   the externally supplied builtin object pointer
 * @param   builtin external supershell builtin function
 */ 
VEXTERNC int Vsh_shell(Vsh *thee, char *pPR, void *pthee,
    int (*builtin)(void *thee, int argc, char **argv));

/**
 * @ingroup Vsh
 * @brief   Place a variable with a value in the environment. 
 * @author  Michael Holst
 * @note    Class Vsh: Non-inlineable method (vsh.c) 
 * @return  Success enumeration
 * @param   thee Pointer to the Vsh object
 * @param   envi environment string
 * @param   valu value string
 */ 
VEXTERNC int Vsh_putenv(Vsh *thee, const char *envi, const char *valu);

/**
 * @ingroup Vsh
 * @brief   Place a variable with an info string in the environment. 
 * @author  Michael Holst
 * @note    Class Vsh: Non-inlineable method (vsh.c) 
 * @return  Success enumeration
 * @param   thee Pointer to the Vsh object
 * @param   envi environment string
 * @param   valu value string
 */ 
VEXTERNC int Vsh_putenvInfo(Vsh *thee, const char *envi, const char *valu);

/**
 * @ingroup Vsh
 * @brief   Place a variable with a value (integer) in the environment.  
 * @author  Michael Holst
 * @note    Class Vsh: Non-inlineable method (vsh.c) 
 * @return  Success enumeration
 * @param   thee Pointer to the Vsh object
 * @param   envi environment string
 * @param   valu value string
 */ 
VEXTERNC int Vsh_putenvInt(Vsh *thee, const char *envi, const int valu);

/**
 * @ingroup Vsh
 * @brief   Place a variable with a value (real) in the environment. 
 * @author  Michael Holst
 * @note    Class Vsh: Non-inlineable method (vsh.c) 
 * @return  Success enumeration
 * @param   thee Pointer to the Vsh object
 * @param   envi environment string
 * @param   valu value string
 */ 
VEXTERNC int Vsh_putenvReal(Vsh *thee, const char *envi, const double valu);

/**
 * @ingroup Vsh
 * @brief   Get a value of variable in the environment.   
 * @author  Michael Holst
 * @note    Class Vsh: Non-inlineable method (vsh.c) 
 * @return  value of variable in the environment.
 * @param   thee Pointer to the Vsh object
 * @param   envi environment string
 */ 
VEXTERNC char *Vsh_getenv(Vsh *thee, const char *envi);

/**
 * @ingroup Vsh
 * @brief   Get info associated with a variable in the environment.
 * @author  Michael Holst
 * @note    Class Vsh: Non-inlineable method (vsh.c) 
 * @return  info associated with a variable in the environment.
 * @param   thee Pointer to the Vsh object
 * @param   envi environment string
 */ 
VEXTERNC char *Vsh_getenvInfo(Vsh *thee, const char *envi);

/**
 * @ingroup Vsh
 * @brief   Get a value of variable in the environment as an integer.
 * @author  Michael Holst
 * @note    Class Vsh: Non-inlineable method (vsh.c) 
 * @return  value of variable in the environment as an integer.
 * @param   thee Pointer to the Vsh object
 * @param   envi environment string
 */ 
VEXTERNC int Vsh_getenvInt(Vsh *thee, const char *envi);

/**
 * @ingroup Vsh
 * @brief   Get a value of variable in the environment as a real. 
 * @author  Michael Holst
 * @note    Class Vsh: Non-inlineable method (vsh.c)
 * @return  value of variable in the environment as a real.
 * @param   thee Pointer to the Vsh object
 * @param   envi environment string
 */ 
VEXTERNC double Vsh_getenvReal(Vsh *thee, const char *envi);

/**
 * @ingroup Vsh
 * @brief   Remove a variable from the environment.  
 * @author  Michael Holst
 * @note    Class Vsh: Non-inlineable method (vsh.c) 
 * @return  None
 * @param   thee Pointer to the Vsh object
 * @param   envi environment string
 */ 
VEXTERNC void Vsh_remove(Vsh *thee, const char *envi);

/**
 * @ingroup Vsh
 * @brief   Wipe the environment.
 * @author  Michael Holst
 * @note    Class Vsh: Non-inlineable method (vsh.c) 
 * @return  None
 * @param   thee Pointer to the Vsh object
 */ 
VEXTERNC void Vsh_wipe(Vsh *thee);

/**
 * @ingroup Vsh
 * @brief   Print the exact current malloc usage.         
 * @author  Michael Holst
 * @note    Class Vsh: Non-inlineable method (vsh.c) 
 * @return  None
 * @param   thee Pointer to the Vsh object
 */ 
VEXTERNC void Vsh_memChk(Vsh *thee);

/**
 * @ingroup Vsh
 * @brief   Setup for an I/O command.   
 * @author  Michael Holst
 * @note    Class Vsh: Non-inlineable method (vsh.c) 
 * @return  socket for reading/writing the external data
 * @param   thee Pointer to the Vsh object
 * @param   key  Pointer to different read/write option
 */ 
VEXTERNC Vio *Vsh_ioSetup(Vsh *thee, char *key);

/**
 * @ingroup Vsh
 * @brief   Cleanup an I/O command.      
 * @author  Michael Holst
 * @note    Class Vsh: Non-inlineable method (vsh.c) 
 * @return  None
 * @param   thee Pointer to the Vsh object
 * @param   sock socket for reading/writing the external data
 */ 
VEXTERNC void Vsh_ioCleanup(Vsh *thee, Vio **sock);

#endif /* _VSH_H_ */

