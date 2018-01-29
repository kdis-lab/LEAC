/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* All rights reserved.                                              */
/*                                                                   */
/* Redistribution and use in source and binary forms, with or        */
/* without modification, are permitted provided that the following   */
/* conditions are met:                                               */
/*                                                                   */
/*   1. Redistributions of source code must retain the above         */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer.                                                  */
/*                                                                   */
/*   2. Redistributions in binary form must reproduce the above      */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer in the documentation and/or other materials       */
/*      provided with the distribution.                              */
/*                                                                   */
/*    THIS  SOFTWARE IS PROVIDED  BY THE  UNIVERSITY OF  TEXAS AT    */
/*    AUSTIN  ``AS IS''  AND ANY  EXPRESS OR  IMPLIED WARRANTIES,    */
/*    INCLUDING, BUT  NOT LIMITED  TO, THE IMPLIED  WARRANTIES OF    */
/*    MERCHANTABILITY  AND FITNESS FOR  A PARTICULAR  PURPOSE ARE    */
/*    DISCLAIMED.  IN  NO EVENT SHALL THE UNIVERSITY  OF TEXAS AT    */
/*    AUSTIN OR CONTRIBUTORS BE  LIABLE FOR ANY DIRECT, INDIRECT,    */
/*    INCIDENTAL,  SPECIAL, EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES    */
/*    (INCLUDING, BUT  NOT LIMITED TO,  PROCUREMENT OF SUBSTITUTE    */
/*    GOODS  OR  SERVICES; LOSS  OF  USE,  DATA,  OR PROFITS;  OR    */
/*    BUSINESS INTERRUPTION) HOWEVER CAUSED  AND ON ANY THEORY OF    */
/*    LIABILITY, WHETHER  IN CONTRACT, STRICT  LIABILITY, OR TORT    */
/*    (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY WAY OUT    */
/*    OF  THE  USE OF  THIS  SOFTWARE,  EVEN  IF ADVISED  OF  THE    */
/*    POSSIBILITY OF SUCH DAMAGE.                                    */
/*                                                                   */
/* The views and conclusions contained in the software and           */
/* documentation are those of the authors and should not be          */
/* interpreted as representing official policies, either expressed   */
/* or implied, of The University of Texas at Austin.                 */
/*********************************************************************/

#if defined(__PGI) || defined(__PGIC__)
COMPILER_PGI
#endif

#if defined(__PATHSCALE__) || defined(__PATHCC__)
COMPILER_PATHSCALE
#endif

#if defined(__INTEL_COMPILER) || defined(__ICC) || defined(__ECC)
COMPILER_INTEL
#endif

#if defined(__OPENCC__)
COMPILER_OPEN64
#endif

#if defined(__SUNPRO_C)
COMPILER_SUN
#endif

#if defined(__IBMC__) || defined(__xlc__)
COMPILER_IBM
#endif

#if defined(__DECCC__)
COMPILER_DEC
#endif

#if defined(__GNUC__)
COMPILER_GNU
#endif

#if defined(__linux__)
OS_LINUX
#endif

#if defined(__FreeBSD__)
OS_FreeBSD
#endif

#if defined(__NetBSD__)
OS_NetBSD
#endif

#if defined(__sun)
OS_SunOS
#endif

#if defined(__APPLE__)
OS_Darwin
#endif

#if defined(_AIX)
OS_AIX
#endif

#if defined(__OSF)
OS_OSF
#endif

#if defined(__WIN32) || defined(__WIN64) || defined(__WINNT)
OS_WINNT
#endif

#if defined(__CYGWIN__)
OS_CYGWIN
#endif

#if defined(__INTERIX)
OS_INTERIX
#endif

#if defined(__i386) || defined(_X86)
ARCH_X86
#endif

#if defined(__x86_64__) || defined(__amd64__)
ARCH_X86_64
#endif

#if defined(__powerpc___) || defined(__PPC__) || defined(_POWER)
ARCH_POWER
#endif

#ifdef __mips64
ARCH_MIPS64
#endif

#if defined(__mips32) || defined(__mips)
ARCH_MIPS32
#endif

#ifdef __alpha
ARCH_ALPHA
#endif

#if defined(__sparc) || defined(__sparc__)
ARCH_SPARC
#endif

#if defined(__ia64__) || defined(__ia64)
ARCH_IA64
#endif

#if defined(__LP64) || defined(__LP64__) || defined(__ptr64) || defined(__x86_64__) || defined(__amd64__) || defined(__64BIT__)
BINARY_64
#endif
