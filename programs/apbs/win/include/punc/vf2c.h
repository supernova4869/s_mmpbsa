/**
 *  @file       vf2c.h
 *  @brief      Standard Fortran to C header file 
 *  @note       barf  [ba:rf]  
 *              2.  "He suggested using FORTRAN, and everybody barfed."
 *    	- From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition) 
 *  @version    $Id: vf2c.h,v 1.7 2010/08/12 05:52:33 fetk Exp $
 *  @author     Michael Holst
 *  
 *  @attention
 *  @verbatim
 *
 * PUNC = < Portable Understructure for Numerical Computing >
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

#ifndef F2C_INCLUDE
#define F2C_INCLUDE

typedef long int integer;
typedef unsigned long int uinteger;
typedef char *address;
typedef short int shortint;
typedef float real;
typedef double doublereal;
/** @brief datastructure for single-precision complex numbers */	
typedef struct { real r, i; } realcomplex;
/** @brief datastructure for double-precision complex numbers */	
typedef struct { doublereal r, i; } doublecomplex;
typedef long int logical;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;
#ifdef INTEGER_STAR_8	
/** @brief system-dependent. Adjust for integer*8. */
typedef long long longint;	
/** @brief system-dependent */	
typedef unsigned long long ulongint;   
#define qbit_clear(a,b)	((a) & ~((ulongint)1 << (b)))
#define qbit_set(a,b)	((a) |  ((ulongint)1 << (b)))
#endif

#define TRUE_ (1)
#define FALSE_ (0)

#ifndef Extern
/** @brief Extern is for use with -E */
#define Extern extern
#endif

/* I/O stuff */

#ifdef f2c_i2
/** @brief for -i2 */
typedef short flag;
typedef short ftnlen;
typedef short ftnint;
#else
typedef long int flag;
typedef long int ftnlen;
typedef long int ftnint;
#endif

/** @brief external read, write*/
typedef struct
{	flag cierr;
	ftnint ciunit;
	flag ciend;
	char *cifmt;
	ftnint cirec;
} cilist;

/** @brief internal read, write*/
typedef struct
{	flag icierr;
	char *iciunit;
	flag iciend;
	char *icifmt;
	ftnint icirlen;
	ftnint icirnum;
} icilist;

/** @brief open*/
typedef struct
{	flag oerr;
	ftnint ounit;
	char *ofnm;
	ftnlen ofnmlen;
	char *osta;
	char *oacc;
	char *ofm;
	ftnint orl;
	char *oblnk;
} olist;

/** @brief close*/
typedef struct
{	flag cerr;
	ftnint cunit;
	char *csta;
} cllist;

/** @brief rewind, backspace, endfile*/
typedef struct
{	flag aerr;
	ftnint aunit;
} alist;

/** @brief inquire */
typedef struct
{	flag inerr;
	ftnint inunit;
	char *infile;
	ftnlen infilen;
	ftnint	*inex;	/**< @brief parameters in standard's order*/
	ftnint	*inopen;
	ftnint	*innum;
	ftnint	*innamed;
	char	*inname;
	ftnlen	innamlen;
	char	*inacc;
	ftnlen	inacclen;
	char	*inseq;
	ftnlen	inseqlen;
	char 	*indir;
	ftnlen	indirlen;
	char	*infmt;
	ftnlen	infmtlen;
	char	*inform;
	ftnint	informlen;
	char	*inunf;
	ftnlen	inunflen;
	ftnint	*inrecl;
	ftnint	*innrec;
	char	*inblank;
	ftnlen	inblanklen;
} inlist;

#define VOID void

/** @brief for multiple entry points */
union Multitype {	
	integer1 g;
	shortint h;
	integer i;
	/* longint j; */
	real r;
	doublereal d;
	realcomplex c;
	doublecomplex z;
	};

typedef union Multitype Multitype;

/*typedef long int Long;*/	/* No longer used; formerly in Namelist */

/** @brief for Vardesc */
struct Vardesc {	
	char *name;
	char *addr;
	ftnlen *dims;
	int  type;
	};
typedef struct Vardesc Vardesc;

/** @brief for Namelist */
struct Namelist {
	char *name;
	Vardesc **vars;
	int nvars;
	};
typedef struct Namelist Namelist;

/*
 * SOME CHANGES FROM STANDARD F2C.H:
 *
 *    Some of the macros below are pre-defined in low-level
 *    ANSI-C headers on some (broken) ANSI-C platforms.
 *    A particular example is the IBM Visual Age C compiler,
 *    as pointed out by N. Baker.
 *
 *    Our solution is to undef these symbols before re-defining them.
 *    This is safe on non-broken ANSI-C platforms which do not
 *    pre-define them, and is the safest possible course on
 *    platforms which do pre-define them.  If we just avoid
 *    re-defining them by checking whether they are already defined,
 *    we will have no idea what implementation is actually used.
 *    Therefore, we brutally un-define and then re-define them.
 *
 * -michael holst (Thu Oct 24 09:54:46 PDT 2002)
 */

#undef abs
#undef dabs
#undef min
#undef max
#undef dmin
#undef dmax

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doublereal)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (doublereal)min(a,b)
#define dmax(a,b) (doublereal)max(a,b)

#define bit_test(a,b)	((a) >> (b) & 1)
#define bit_clear(a,b)	((a) & ~((uinteger)1 << (b)))
#define bit_set(a,b)	((a) |  ((uinteger)1 << (b)))

/* procedure parameter types for -A and -C++ */

#define F2C_proc_par_types 1
#ifdef __cplusplus
typedef int /* Unknown procedure type */ (*U_fp)(...);
typedef shortint (*J_fp)(...);
typedef integer (*I_fp)(...);
typedef real (*R_fp)(...);
typedef doublereal (*D_fp)(...), (*E_fp)(...);
typedef /* Complex */ VOID (*C_fp)(...);
typedef /* Double Complex */ VOID (*Z_fp)(...);
typedef logical (*L_fp)(...);
typedef shortlogical (*K_fp)(...);
typedef /* Character */ VOID (*H_fp)(...);
typedef /* Subroutine */ int (*S_fp)(...);
#else
typedef int /* Unknown procedure type */ (*U_fp)();
typedef shortint (*J_fp)();
typedef integer (*I_fp)();
typedef real (*R_fp)();
typedef doublereal (*D_fp)(), (*E_fp)();
typedef /* Complex */ VOID (*C_fp)();
typedef /* Double Complex */ VOID (*Z_fp)();
typedef logical (*L_fp)();
typedef shortlogical (*K_fp)();
typedef /* Character */ VOID (*H_fp)();
typedef /* Subroutine */ int (*S_fp)();
#endif
/* E_fp is for real functions when -R is not specified */
typedef VOID C_f;	/**< @brief complex function */
typedef VOID H_f;	/**< @brief character function */
typedef VOID Z_f;	/**< @brief double complex function */
typedef doublereal E_f;	/**< @brief real function with -R not specified */

/* undef any lower-case symbols that your C compiler predefines, e.g.: */

#ifndef Skip_f2c_Undefs
#undef cray
#undef gcos
#undef mc68010
#undef mc68020
#undef mips
#undef pdp11
#undef sgi
#undef sparc
#undef sun
#undef sun2
#undef sun3
#undef sun4
#undef u370
#undef u3b
#undef u3b2
#undef u3b5
#undef unix
#undef vax
#endif

/* If you are using a C++ compiler, append the following to f2c.h
   for compiling libF77 and libI77. */

#ifdef __cplusplus
extern "C" {
extern int abort_(void);
extern double c_abs(realcomplex *);
extern void c_cos(realcomplex *, realcomplex *);
extern void c_div(realcomplex *, realcomplex *, realcomplex *);
extern void c_exp(realcomplex *, realcomplex *);
extern void c_log(realcomplex *, realcomplex *);
extern void c_sin(realcomplex *, realcomplex *);
extern void c_sqrt(realcomplex *, realcomplex *);
extern double d_abs(double *);
extern double d_acos(double *);
extern double d_asin(double *);
extern double d_atan(double *);
extern double d_atn2(double *, double *);
extern void d_cnjg(doublecomplex *, doublecomplex *);
extern double d_cos(double *);
extern double d_cosh(double *);
extern double d_dim(double *, double *);
extern double d_exp(double *);
extern double d_imag(doublecomplex *);
extern double d_int(double *);
extern double d_lg10(double *);
extern double d_log(double *);
extern double d_mod(double *, double *);
extern double d_nint(double *);
extern double d_prod(float *, float *);
extern double d_sign(double *, double *);
extern double d_sin(double *);
extern double d_sinh(double *);
extern double d_sqrt(double *);
extern double d_tan(double *);
extern double d_tanh(double *);
extern double derf_(double *);
extern double derfc_(double *);
extern integer do_fio(ftnint *, char *, ftnlen);
extern integer do_lio(ftnint *, ftnint *, char *, ftnlen);
extern integer do_uio(ftnint *, char *, ftnlen);
extern integer e_rdfe(void);
extern integer e_rdue(void);
extern integer e_rsfe(void);
extern integer e_rsfi(void);
extern integer e_rsle(void);
extern integer e_rsli(void);
extern integer e_rsue(void);
extern integer e_wdfe(void);
extern integer e_wdue(void);
extern integer e_wsfe(void);
extern integer e_wsfi(void);
extern integer e_wsle(void);
extern integer e_wsli(void);
extern integer e_wsue(void);
extern int ef1asc_(ftnint *, ftnlen *, ftnint *, ftnlen *);
extern integer ef1cmc_(ftnint *, ftnlen *, ftnint *, ftnlen *);
extern double erf(double);
extern double erf_(float *);
extern double erfc(double);
extern double erfc_(float *);
extern integer f_back(alist *);
extern integer f_clos(cllist *);
extern integer f_end(alist *);
extern void f_exit(void);
extern integer f_inqu(inlist *);
extern integer f_open(olist *);
extern integer f_rew(alist *);
extern int flush_(void);
extern void getarg_(integer *, char *, ftnlen);
extern void getenv_(char *, char *, ftnlen, ftnlen);
extern short h_abs(short *);
extern short h_dim(short *, short *);
extern short h_dnnt(double *);
extern short h_indx(char *, char *, ftnlen, ftnlen);
extern short h_len(char *, ftnlen);
extern short h_mod(short *, short *);
extern short h_nint(float *);
extern short h_sign(short *, short *);
extern short hl_ge(char *, char *, ftnlen, ftnlen);
extern short hl_gt(char *, char *, ftnlen, ftnlen);
extern short hl_le(char *, char *, ftnlen, ftnlen);
extern short hl_lt(char *, char *, ftnlen, ftnlen);
extern integer i_abs(integer *);
extern integer i_dim(integer *, integer *);
extern integer i_dnnt(double *);
extern integer i_indx(char *, char *, ftnlen, ftnlen);
extern integer i_len(char *, ftnlen);
extern integer i_mod(integer *, integer *);
extern integer i_nint(float *);
extern integer i_sign(integer *, integer *);
extern integer iargc_(void);
extern ftnlen l_ge(char *, char *, ftnlen, ftnlen);
extern ftnlen l_gt(char *, char *, ftnlen, ftnlen);
extern ftnlen l_le(char *, char *, ftnlen, ftnlen);
extern ftnlen l_lt(char *, char *, ftnlen, ftnlen);
extern void pow_ci(realcomplex *, realcomplex *, integer *);
extern double pow_dd(double *, double *);
extern double pow_di(double *, integer *);
extern short pow_hh(short *, shortint *);
extern integer pow_ii(integer *, integer *);
extern double pow_ri(float *, integer *);
extern void pow_zi(doublecomplex *, doublecomplex *, integer *);
extern void pow_zz(doublecomplex *, doublecomplex *, doublecomplex *);
extern double r_abs(float *);
extern double r_acos(float *);
extern double r_asin(float *);
extern double r_atan(float *);
extern double r_atn2(float *, float *);
extern void r_cnjg(realcomplex *, realcomplex *);
extern double r_cos(float *);
extern double r_cosh(float *);
extern double r_dim(float *, float *);
extern double r_exp(float *);
extern double r_imag(realcomplex *);
extern double r_int(float *);
extern double r_lg10(float *);
extern double r_log(float *);
extern double r_mod(float *, float *);
extern double r_nint(float *);
extern double r_sign(float *, float *);
extern double r_sin(float *);
extern double r_sinh(float *);
extern double r_sqrt(float *);
extern double r_tan(float *);
extern double r_tanh(float *);
extern void s_cat(char *, char **, integer *, integer *, ftnlen);
extern integer s_cmp(char *, char *, ftnlen, ftnlen);
extern void s_copy(char *, char *, ftnlen, ftnlen);
extern int s_paus(char *, ftnlen);
extern integer s_rdfe(cilist *);
extern integer s_rdue(cilist *);
extern integer s_rnge(char *, integer, char *, integer);
extern integer s_rsfe(cilist *);
extern integer s_rsfi(icilist *);
extern integer s_rsle(cilist *);
extern integer s_rsli(icilist *);
extern integer s_rsne(cilist *);
extern integer s_rsni(icilist *);
extern integer s_rsue(cilist *);
extern int s_stop(char *, ftnlen);
extern integer s_wdfe(cilist *);
extern integer s_wdue(cilist *);
extern integer s_wsfe(cilist *);
extern integer s_wsfi(icilist *);
extern integer s_wsle(cilist *);
extern integer s_wsli(icilist *);
extern integer s_wsne(cilist *);
extern integer s_wsni(icilist *);
extern integer s_wsue(cilist *);
extern void sig_die(char *, int);
extern integer signal_(integer *, void (*)(int));
extern integer system_(char *, ftnlen);
extern double z_abs(doublecomplex *);
extern void z_cos(doublecomplex *, doublecomplex *);
extern void z_div(doublecomplex *, doublecomplex *, doublecomplex *);
extern void z_exp(doublecomplex *, doublecomplex *);
extern void z_log(doublecomplex *, doublecomplex *);
extern void z_sin(doublecomplex *, doublecomplex *);
extern void z_sqrt(doublecomplex *, doublecomplex *);
	}
#endif /* ifdef __cplusplus */

#endif /* ifndef F2C_INCLUDE */

