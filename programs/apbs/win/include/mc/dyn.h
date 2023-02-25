/** 
 * @defgroup SVio SVio class
 * @brief    Vio socket container object
 */

/**
 *  @file       dyn.h
 *  @ingroup    SVio global_mc
 *  @brief      Class Dyn: dynamics library
 *  @author     Stephen Bond and Michael Holst
 *  @note       None
 *  @version    $Id: dyn.h,v 1.9 2010/08/12 05:18:44 fetk Exp $ 
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

#ifndef _DYN_H
#define _DYN_H

#include <mc/mc_base.h>

#include <mc/aprx.h>
#include <mc/nam.h>

/** 
 * @ingroup global_mc
 * @brief   Class Dyn: Parameters and datatypes 
 * @author  Michael Holst 
 */

typedef enum DYNtype {
    SIMP_TYPE,
    TDEP_TYPE,
    NLIN_TYPE
} DYNtype;

/** 
 * @ingroup SVio
 * @brief   Contains public data memebers for the SVio class
 * @author  Michael Holst
 */
struct sSVio {

    /** @brief read/write key.\n
     *  "r" = read, "w" = write                  */
    char rwkey[1];   
    /** @brief device type.                      \n                       
     *  "SDIO" = standard I/O.                   \n
     *  "FILE" = file I/O.                       \n
     *  "BUFF" = buffer I/O.                     \n
     *  "UNIX" = UNIX (domain) socket I/O.       \n
     *  "INET" = INET (network) socket I/O       */
    char iodev[VMAX_BUFSIZE];
    /** @brief data format.                      \n
     *  "ASC" = ASCII (FILE,BUFF,UNIX,INET).     \n 
     *  "XDR" = BINARY (FILE,BUFF,UNIX,INET)     */ 
    char iofmt[VMAX_BUFSIZE];  
    /** @brief local hostname (me) (UNIX,INET)           */
    char iohost[VMAX_BUFSIZE];
    /** @brief file or device name (FILE,BUFF,UNIX,INET) */
    char iofile[VMAX_BUFSIZE]; 

    /** @brief GV colKey   (see Gem_writeGV)             */
    int colKeyGV;    
    /** @brief GV chartKey (see Gem_writeGV)             */
    int chartKeyGV; 
    /** @brief GV fkey     (see Gem_writeGV)             */
    int fkeyGV;  
    /** @brief GV gluVal   (see Gem_writeGV)             */
    double gluValGV; 

    /** @brief print format.                     \n          
     *      0 = Geomview (GV)                    \n
     *      1 = GMV                              \n
     *      2 = OpenDX (DX)                      \n
     *      3 = Matlab                           \n
     *      4 = Raw Data                         */
    int printtype;   

    /** @brief Vio socket pointer */
    Vio * sock;

    /** @brief (BUFF) */
    char *buf;
    /** @brief (BUFF) */
    int bufsize;   

};

/** 
 * @ingroup SVio
 * @brief   Delcaration of the SVio class as the SVio structure
 * @author  Michael Holst
 */
typedef struct sSVio SVio;

/**
 * @ingroup global_mc
 * @brief   Class DynPDE definition
 * @author  Michael Holst
 */
typedef struct DynPDE {

    /** @brief Current Time */
    double time;     
    /** @brief Current Energy Key */
    int ekey; 
    /** @brief NAB/YHS:  Function called at each time step to evaluate
     *         user-defined observables (can be VNULL) */ 
    void (*userStepHook)(PDE *thee, AM *am, int ekey);  
    /** @brief NAB/YHS:  \n
     *         Set to 1 if userStepHook function defined; \n
     *         set to 0 otherwise */
    int haveUserStepHook;  

} DynPDE;

/*
 * ***************************************************************************
 * Class Dyn: Inlineable methods (dyn.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_DYN)
#else /* if defined(VINLINE_DYN) */
#endif /* if !defined(VINLINE_DYN) */

/*
 * ***************************************************************************
 * Class Dyn: Non-Inlineable methods (dyn.c)
 * ***************************************************************************
 */

/**
 * @ingroup AM
 * @brief   Solution of Time-Dependent problems by Method of Lines.
 * @authors Stephen Bond and Kaihsu Tai 
 * @note    Class Dyn: Non-Inlineable methods (dyn.c) 
 * @return  None
 * @param   thee      Pointer to class AM
 * @param   meth      method choice (0=Forward Euler;1=Backward Euler;2=Trapezoidal Rule)
 * @param   dt        the time step
 * @param   t0        the initial time
 * @param   numstep   number of the time steps
 * @param   pfreq     frequency of the printing contour
 * @param   efreq     frequency of energy printing
 * @param   ekeytotal total number of printing energies
 * @param   pdetype   index for different types of PDEs:SIMP_TYPE/TDEP_TYPE/NLIN_TYPE
 * @param   ltol      error tolerance
 * @param   lmax      number of iterations to do (the maximum allowed)
 * @param   vsock     socket for writing a finite element mesh or mesh function
 */
VEXTERNC void AM_tSolve(AM *thee, int meth, double dt, double t0, int numstep, 
    int pfreq, int efreq, int ekeytotal, DYNtype pdetype, 
    double ltol, int lmax, SVio *vsock);

/**
 * @ingroup PDE
 * @brief   Construct the Dyn PDE datastructure inside the PDE structure.   
 * @author  Stephen Bond 
 * @note    Class Dyn: Non-Inlineable methods (dyn.c) 
 * @return  None
 * @param   thee  Pointer to the differential equation object
 */
VEXTERNC void PDE_initDyn(PDE *thee);

/**
 * @ingroup PDE
 * @brief   Destruct the Dyn PDE datastructure inside the PDE structure.  
 * @author  Stephen Bond 
 * @note    Class Dyn: Non-Inlineable methods (dyn.c) 
 * @return  None
 * @param   thee  Pointer to the differential equation object
 */
VEXTERNC void PDE_killDyn(PDE *thee);

/**
 * @ingroup PDE
 * @brief   Checks to see if the Dyn PDE structure has been created. 
 * @author  Stephen Bond 
 * @note    Class Dyn: Non-Inlineable methods (dyn.c) 
 * @return  Success enumeration
 * @param   thee  Pointer to the differential equation object
 */
VEXTERNC int PDE_checkDyn(PDE *thee);

/**
 * @ingroup PDE
 * @brief   Set the time variable in the Dyn PDE structure  
 * @author  Stephen Bond 
 * @note    Class Dyn: Non-Inlineable methods (dyn.c) 
 * @return  None
 * @param   thee    Pointer to the differential equation object
 * @param   mytime  the time variable in the Dyn PDE structure
 */
VEXTERNC void PDE_setTime(PDE *thee, double mytime);

/**
 * @ingroup PDE
 * @brief   Get the time variable from the Dyn PDE structure  
 * @author  Stephen Bond 
 * @note    Class Dyn: Non-Inlineable methods (dyn.c) 
 * @return  the time variable from the Dyn PDE structure
 * @param   thee    Pointer to the differential equation object
 */
VEXTERNC double PDE_getTime(PDE *thee);

/**
 * @ingroup PDE
 * @brief   Set the energy key variable in the Dyn PDE structure  
 * @author  Stephen Bond 
 * @note    Class Dyn: Non-Inlineable methods (dyn.c) 
 * @return  None
 * @param   thee          Pointer to the differential equation object
 * @param   UserStepHook  Pointer to user defined external function
 */
VEXTERNC void PDE_setUserStepHook(PDE *thee, 
    void (*UserStepHook)(PDE *thee, AM *am, int ekey));

/**
 * @ingroup PDE
 * @brief   Set the UserStepHook function pointer to VNULL     
 * @author  Stephen Bond 
 * @note    Class Dyn: Non-Inlineable methods (dyn.c) 
 * @return  None
 * @param   thee  Pointer to the differential equation object
 */
VEXTERNC void PDE_nullUserStepHook(PDE *thee);

/**
 * @ingroup PDE
 * @brief   Call the userStepHook function  
 * @author  Stephen Bond 
 * @note    Class Dyn: Non-Inlineable methods (dyn.c) 
 * @return  None
 * @param   thee  Pointer to the differential equation object
 * @param   am    Pointer to class AM
 * @param   ekey  the energy key variable in the Dyn PDE structure
 */
VEXTERNC void PDE_userStepHook(PDE *thee, AM *am, int ekey);  

/**
 * @ingroup PDE
 * @brief   Set the energy key variable in the Dyn PDE structure 
 * @author  Stephen Bond 
 * @note    Class Dyn: Non-Inlineable methods (dyn.c) 
 * @return  None
 * @param   thee  Pointer to the differential equation object
 * @param   ekey  the energy key variable in the Dyn PDE structure
 */
VEXTERNC void PDE_setEnergyKey(PDE *thee, int ekey);

/**
 * @ingroup PDE
 * @brief   Get the energy key variable from the Dyn PDE structure        
 * @author  Stephen Bond 
 * @note    Class Dyn: Non-Inlineable methods (dyn.c) 
 * @return  None
 * @param   thee  Pointer to the differential equation object
 */
VEXTERNC int PDE_getEnergyKey(PDE *thee);

/*
 * ***************************************************************************
 * Class Dyn: Non-Inlineable methods (svio.c)
 * ***************************************************************************
 */

/**
 * @ingroup SVio
 * @brief   Construct the Vio socket container object. 
 * @author  Stephen Bond 
 * @note    Class Dyn: Non-Inlineable methods (svio.c) 
 * @return  None
 */
VEXTERNC SVio* SVio_ctor(void);

/**
 * @ingroup SVio
 * @brief   Destroy the Vio socket container object.
 * @author  Stephen Bond 
 * @note    Class Dyn: Non-Inlineable methods (svio.c) 
 * @return  None
 * @param   vsock  socket for reading/writing a finite element mesh or mesh function
 */
VEXTERNC void SVio_dtor(SVio **vsock);

/**
 * @ingroup Vsh
 * @brief   Setup for I/O commands using SVio. 
 * @author  Stephen Bond 
 * @note    Class Dyn: Non-Inlineable methods (svio.c) 
 * @return  socket for reading/writing a finite element mesh or mesh function
 * @param   thee  Pointer to the shell with environment variables
 * @param   key   Pointer to I/O commands using SVio
 */
VEXTERNC SVio* Vsh_SVioSetup(Vsh *thee, char *key);

/**
 * @ingroup Vsh
 * @brief   Cleanup after I/O commands using SVio.
 * @author  Stephen Bond 
 * @note    Class Dyn: Non-Inlineable methods (svio.c) 
 * @return  None
 * @param   thee   Pointer to the shell with environment variables
 * @param   vsock  socket for reading/writing a finite element mesh or mesh function
 */ 
VEXTERNC void Vsh_SVioCleanup(Vsh *thee, SVio **vsock);

/**
 * @ingroup SVio
 * @brief   Initialize the Vio socket container object.         
 * @author  Stephen Bond 
 * @note    Class Dyn: Non-Inlineable methods (svio.c) 
 * @return  None
 * @param   thee     Pointer to the shell with environment variables
 * @param   rwkey    read/write key "r" = read, "w" = write
 * @param   iodev    device type:
 *                   "SDIO" = standard I/O
 *                   "FILE" = file I/O
 *                   "BUFF" = buffer I/O
 *                   "UNIX" = UNIX (domain) socket I/O
 *                   "INET" = INET (network) socket I/O
 * @param   iofmt    data format:"ASC" = ASCII (FILE,BUFF,UNIX,IN
 * @param   iohost   local hostname (me) (UNIX,INET)
 * @param   iofile   file or device name (FILE,BUFF,UNIX,INET)
 * @param   buf      Pointer to BUFF
 * @param   bufsize  size of BUFF
 * @param   ptype    print format(0 = Geomview (GV),1 = GMV,2 = OpenDX (DX),
 *                   3 = Matlab,4 = Raw Data)
 */
VEXTERNC void SVio_initStructure(SVio *thee, const char *rwkey,
    const char *iodev, const char *iofmt, const char *iohost,
    const char *iofile, char *buf, int bufsize, int ptype);

/**
 * @ingroup SVio
 * @brief   Return the length of the internal buffer.
 * @author  Stephen Bond 
 * @note    Class Dyn: Non-Inlineable methods (svio.c) 
 * @return  the length of the internal buffer.
 * @param   thee  Pointer to the shell with environment variables
 */
VEXTERNC int SVio_bufSize(SVio *thee);

/**
 * @ingroup SVio
 * @brief   Return the pointer to the internal buffer.   
 * @author  Stephen Bond 
 * @note    Class Dyn: Non-Inlineable methods (svio.c) 
 * @return  the pointer to the internal buffer.   
 * @param   thee  Pointer to the shell with environment variables
 */
VEXTERNC char* SVio_bufGive(SVio *thee);

/**
 * @ingroup SVio
 * @brief   Set the pointer to the internal buffer.           
 * @author  Stephen Bond 
 * @note    Class Dyn: Non-Inlineable methods (svio.c) 
 * @return  None
 * @param   thee     Pointer to the shell with environment variables
 * @param   buf      Pointer to BUFF
 * @param   bufsize  size of BUFF
 */
VEXTERNC void SVio_bufTake(SVio *thee, char *buf, int bufsize);

/**
 * @ingroup SVio
 * @brief   Write a finite element mesh or mesh function to a socket.
 * @author  Stephen Bond 
 * @note    Class Dyn: Non-Inlineable methods (svio.c) 
 * @return  None
 * @param   thee    Pointer to an Aprx allocated memory location
 * @param   plabel  index for printing label
 * @param   w0      Pointer to the block vector
 * @param   vsock   socket for reading/writing a finite element mesh or mesh function
 */
VEXTERNC void Aprx_writeSVio(Aprx *thee, int plabel, Bvec *w0, SVio *vsock);

#endif /* _DYN_H_ */

