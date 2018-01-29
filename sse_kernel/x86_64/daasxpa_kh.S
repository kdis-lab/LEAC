/*! \file daasxpa_kh.S	
 * This file is part of the LEAC.
 *
 * Implementation of the 
 *	
 *	a_ij = a_ij + alpha(a_ij - xj),
 *
 *  function for double
 *
 *  Based on Optimized BLAS libraries
 *
 * (c)  Hermes Robles Berumen <hermes@uaz.edu.mx>
 *
 * For the full copyright and license information, please view the LICENSE
 * file that was distributed with this source code.
 */
	
/*********************************************************************/
/*                                                                   */
/*             Optimized BLAS libraries                              */
/*                     By Kazushige Goto <kgoto@tacc.utexas.edu>     */
/*                                                                   */
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


#define ASSEMBLER
#include "common_arc.h"	
#include "l2param.h"
#ifndef WINDOWS_ABI

#define STACKSIZE	64
	
#define OLD_M	  %rdi
#define OLD_N	  %rsi
#define OLD_A	  %rcx
#define OLD_LDA	  %r8
#define STACK_INCX	 8 + STACKSIZE(%rsp)
#define STACK_Y		16 + STACKSIZE(%rsp)
#define STACK_INCY	24 + STACKSIZE(%rsp)
#define STACK_BUFFER	32 + STACKSIZE(%rsp)

#else

#define STACKSIZE	256
	
#define OLD_M	  %rcx
#define OLD_N	  %rdx
#define OLD_A		 40 + STACKSIZE(%rsp)
#define OLD_LDA		 48 + STACKSIZE(%rsp)
#define OLD_X		 56 + STACKSIZE(%rsp)
#define STACK_INCX	 64 + STACKSIZE(%rsp)
#define STACK_Y		 72 + STACKSIZE(%rsp)
#define STACK_INCY	 80 + STACKSIZE(%rsp)
#define STACK_BUFFER	 88 + STACKSIZE(%rsp)

#endif

#define LDA	%r8
#define X	%r9

#define INCX	%rsi
#define INCY	%rdi

#define M	%r10
#define N	%r11
#define A	%r12
#define Y	%r14
#define BUFFER	%r13

#define I	%rax
#define A1	%rbx
#define A2	%rcx
#define LDA3	%rdx
#define Y1	%rbp
#define X1	%r15
	
	
#define MM	M
	
#define SIZE    8
	
#define ALPHA	%xmm15

/*	PROLOGUE
	PROFCODE
*/
.file	"daasxpa_kh.S"
	.text
	.p2align 4,,15
.globl daasxpa_kh
	.type	daasxpa_kh, @function
daasxpa_kh:
.LFB2:
	
	subq	$STACKSIZE, %rsp
	movq	%rbx,  0(%rsp)
	movq	%rbp,  8(%rsp)
	movq	%r12, 16(%rsp)
	movq	%r13, 24(%rsp)
	movq	%r14, 32(%rsp)
	movq	%r15, 40(%rsp)

#ifdef WINDOWS_ABI
	movq	%rdi,    48(%rsp)
	movq	%rsi,    56(%rsp)
	movups	%xmm6,   64(%rsp)
	movups	%xmm7,   80(%rsp)
	movups	%xmm8,   96(%rsp)
	movups	%xmm9,  112(%rsp)
	movups	%xmm10, 128(%rsp)
	movups	%xmm11, 144(%rsp)
	movups	%xmm12, 160(%rsp)
	movups	%xmm13, 176(%rsp)
	movups	%xmm14, 192(%rsp)
	movups	%xmm15, 208(%rsp)

	movq	OLD_M,	      M
	movq	OLD_N,        N
	movq	OLD_A,        A
	movq	OLD_LDA,      LDA
	movq	OLD_X,        X
#else
	movq	OLD_M,	      M
	movq	OLD_N,        N
	movq	OLD_A,        A
	movq	OLD_LDA,      LDA
#endif

	leaq	(,LDA,  SIZE), LDA
	leaq	(LDA, LDA, 2), LDA3

	subq	$-16 * SIZE, A
	subq	$-16 * SIZE, X

	/*CHANGE SIGNE*/
	
#ifdef HAVE_SSE3
#ifndef WINDOWS_ABI
	movddup	%xmm0, ALPHA
#else
	movddup	%xmm3, ALPHA
#endif /*WINDOWS_ABI*/
#else
#ifndef WINDOWS_ABI
	movapd	%xmm0, ALPHA
#else
	movapd	%xmm3, ALPHA
#endif /*WINDOWS_ABI*/
	unpcklpd ALPHA, ALPHA
#endif /*HAVE_SSE3*/

	
	testq	M, M
	jle	.L999
	testq	N, N
	jle	.L999

	movq	X, X1  
	movq	M,  I
	
	cmpq	$4, N
	jl	.L30
	ALIGN_3
	
.L21:

	subq	$4, N 

	movq	X, X1  
	
	movq	A, A1
	leaq	(A1, LDA, 2), A2
	leaq	(A1, LDA, 4), A

#ifdef PREFETCHW
	PREFETCHW	3 * SIZE(Y)   
#endif

	movq	M,  I  
	sarq	$3,  I
	jle	.L25
	
	MOVUPS_A1(-16 * SIZE, A1, %xmm8)
	MOVUPS_A2(-16 * SIZE, A1, LDA, 1, %xmm9)
	MOVUPS_A1(-16 * SIZE, A2, %xmm10)
	MOVUPS_A2(-16 * SIZE, A2, LDA, 1, %xmm11)

	
	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)
	MOVUPS_XL1(-14 * SIZE, X1, %xmm13)

	decq	I
	jle	.L23
	ALIGN_4

.L22:
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A1)
#endif

	movapd  %xmm8,  %xmm0
	subpd   %xmm12, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm8

	movdqu  %xmm8,-16 * SIZE(A1)
	
	MOVUPS_A1(-14 * SIZE, A1, %xmm8)
	movapd  %xmm9,  %xmm1
	subpd   %xmm12, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm9
	
	movdqu  %xmm9,-16 * SIZE(A1,LDA,1)	
	
	MOVUPS_A2(-14 * SIZE, A1, LDA, 1, %xmm9)
	movapd  %xmm10,  %xmm0
	subpd   %xmm12, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm10
	movdqu  %xmm10, -16 * SIZE(A2)
	
	MOVUPS_A1(-14 * SIZE, A2, %xmm10)
	movapd  %xmm11,  %xmm1
	subpd   %xmm12, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm11
	movdqu  %xmm11,-16 * SIZE(A2,LDA,1)
	
	MOVUPS_XL1(-12 * SIZE, X1, %xmm12)
	MOVUPS_A2(-14 * SIZE, A2, LDA, 1, %xmm11)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A1, LDA)
#endif

	movapd  %xmm8,  %xmm0
	subpd   %xmm13, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm8
	movdqu	%xmm8, -14 * SIZE(A1)
	MOVUPS_A1(-12 * SIZE, A1, %xmm8)

	movapd  %xmm9,  %xmm1
	subpd   %xmm13, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm9
	movdqu  %xmm9, -14 * SIZE(A1,LDA,1)
	
	MOVUPS_A2(-12 * SIZE, A1, LDA, 1, %xmm9)
	movapd  %xmm10,  %xmm0
	subpd   %xmm13, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm10
	movdqu  %xmm10,-14 * SIZE(A2)
	MOVUPS_A1(-12 * SIZE, A2, %xmm10)

	movapd  %xmm11,  %xmm1
	subpd   %xmm13, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm11
	MOVUPS_XL1(-10 * SIZE, X1, %xmm13)
	movdqu  %xmm11,-14 * SIZE(A2,LDA,1)
	MOVUPS_A2(-12 * SIZE, A2, LDA, 1, %xmm11)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A2)
#endif

	movapd  %xmm8,  %xmm0
	subpd   %xmm12, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm8
	movdqu  %xmm8, -12 * SIZE(A1)
	
	MOVUPS_A1(-10 * SIZE, A1, %xmm8)
	movapd  %xmm9,  %xmm1
	subpd   %xmm12, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm9
	movdqu  %xmm9, -12 * SIZE(A1,LDA,1)
	
	MOVUPS_A2(-10 * SIZE, A1, LDA, 1, %xmm9)
	movapd  %xmm10,  %xmm0
	subpd   %xmm12, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm10
	movdqu  %xmm10, -12 * SIZE(A2)
	
	MOVUPS_A1(-10 * SIZE, A2, %xmm10)
	movapd  %xmm11,  %xmm1
	subpd   %xmm12, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm11
	MOVUPS_XL1( -8 * SIZE, X1, %xmm12)
	movdqu 	%xmm11, -12 * SIZE(A2,LDA,1)
	MOVUPS_A2(-10 * SIZE, A2, LDA, 1, %xmm11)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A2, LDA)
#endif

	movapd  %xmm8,  %xmm0
	subpd   %xmm13, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm8
	movdqu  %xmm8, -10 * SIZE(A1)
	
	MOVUPS_A1( -8 * SIZE, A1, %xmm8)
	movapd  %xmm9,  %xmm1
	subpd   %xmm13, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm9
	movdqu  %xmm9,-10 * SIZE(A1,LDA,1)
	MOVUPS_A2( -8 * SIZE, A1, LDA, 1, %xmm9)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(X1)
#endif

	movapd  %xmm10,  %xmm0
	subpd   %xmm13, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm10
	movdqu  %xmm10, -10 * SIZE(A2)
	MOVUPS_A1( -8 * SIZE, A2, %xmm10)

	movapd  %xmm11,  %xmm1
	subpd   %xmm13, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm11
	MOVUPS_XL1( -6 * SIZE, X1, %xmm13)
	movdqu  %xmm11,-10 * SIZE(A2,LDA,1)
	MOVUPS_A2( -8 * SIZE, A2, LDA, 1, %xmm11)

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1

	decq	I
	jg	.L22
	ALIGN_4

.L23:
	movapd  %xmm8,  %xmm0
	subpd   %xmm12, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm8
	movdqu  %xmm8, -16 * SIZE(A1)
	
	MOVUPS_A1(-14 * SIZE, A1, %xmm8)
	movapd  %xmm9,  %xmm1
	subpd   %xmm12, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm9
	movdqu  %xmm9, -16 * SIZE(A1,LDA,1)
	
	MOVUPS_A2(-14 * SIZE, A1, LDA, 1, %xmm9)
	movapd  %xmm10,  %xmm0
	subpd   %xmm12, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm10
	movdqu  %xmm10, -16 * SIZE(A2)
	
	MOVUPS_A1(-14 * SIZE, A2, %xmm10)
	movapd  %xmm11,  %xmm1
	subpd   %xmm12, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm11
	MOVUPS_XL1(-12 * SIZE, X1, %xmm12)
	movdqu  %xmm11, -16 * SIZE(A2,LDA,1) 
	MOVUPS_A2(-14 * SIZE, A2, LDA, 1, %xmm11)

	movapd  %xmm8,  %xmm0
	subpd   %xmm13, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm8
	movdqu  %xmm8, -14 * SIZE(A1)
	
	MOVUPS_A1(-12 * SIZE, A1, %xmm8)
	movapd  %xmm9,  %xmm1
	subpd   %xmm13, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm9
	movdqu  %xmm9, -14 * SIZE(A1,LDA,1)

	MOVUPS_A2(-12 * SIZE, A1, LDA, 1, %xmm9)
	movapd  %xmm10,  %xmm0
	subpd   %xmm13, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm10
	movdqu  %xmm10,-14 * SIZE(A2)
	
	MOVUPS_A1(-12 * SIZE, A2, %xmm10)
	movapd  %xmm11,  %xmm1
	subpd   %xmm13, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm11
	MOVUPS_XL1(-10 * SIZE, X1, %xmm13)
	movdqu  %xmm11, -14 * SIZE(A2,LDA,1)  
	MOVUPS_A2(-12 * SIZE, A2, LDA, 1, %xmm11)

	movapd  %xmm8,  %xmm0
	subpd   %xmm12, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm8
	movdqu  %xmm8, -12 * SIZE(A1)
	
	MOVUPS_A1(-10 * SIZE, A1, %xmm8)
	movapd  %xmm9,  %xmm1
	subpd   %xmm12, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm9
	movdqu  %xmm9, -12 * SIZE(A1,LDA,1)
	
	MOVUPS_A2(-10 * SIZE, A1, LDA, 1, %xmm9)
	movapd  %xmm10,  %xmm0
	subpd   %xmm12, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm10
	movdqu  %xmm10, -12 * SIZE(A2)
	MOVUPS_A1(-10 * SIZE, A2, %xmm10)
	movapd  %xmm11,  %xmm1
	subpd   %xmm12, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm11
	movdqu  %xmm11, -12 * SIZE(A2,LDA,1)
	MOVUPS_A2(-10 * SIZE, A2, LDA, 1, %xmm11)

	movapd  %xmm8,  %xmm0
	subpd   %xmm13, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm8
	movdqu  %xmm8, -10 * SIZE(A1)
	
	movapd  %xmm9,  %xmm1
	subpd   %xmm13, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm9
	movdqu  %xmm9, -10 * SIZE(A1,LDA,1)
	
	movapd  %xmm10,  %xmm0
	subpd   %xmm13, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm10
	movdqu  %xmm10, -10 * SIZE(A2)
	
	movapd  %xmm11,  %xmm1
	subpd   %xmm13, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm11
	movdqu  %xmm11, -10 * SIZE(A2,LDA,1)

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1
	ALIGN_4

.L25:
	testq	$4, M   
	jle	.L26

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)  

	MOVUPS_A1(-16 * SIZE, A1, %xmm8)
	movapd  %xmm8,  %xmm0
	subpd   %xmm12, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm8
	movdqu  %xmm8,-16 * SIZE(A1)
	
	MOVUPS_A2(-16 * SIZE, A1, LDA, 1, %xmm9)
	movapd  %xmm9,  %xmm1
	subpd   %xmm12, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm9
	movdqu  %xmm9,-16 * SIZE(A1,LDA,1)
	
	MOVUPS_A1(-16 * SIZE, A2, %xmm10)
	movapd  %xmm10,  %xmm0
	subpd   %xmm12, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm10
	movdqu  %xmm10,-16 * SIZE(A2)
	
	MOVUPS_A2(-16 * SIZE, A2, LDA, 1, %xmm11)
	movapd  %xmm11,  %xmm1
	subpd   %xmm12, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm11
	movdqu  %xmm11,-16 * SIZE(A2,LDA,1)

	MOVUPS_XL1(-14 * SIZE, X1, %xmm13)

	MOVUPS_A1(-14 * SIZE, A1, %xmm8)
	movapd  %xmm8,  %xmm0
	subpd   %xmm13, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm8
	movdqu  %xmm8,-14 * SIZE(A1)
	
	MOVUPS_A2(-14 * SIZE, A1, LDA, 1, %xmm9)
	movapd  %xmm9,  %xmm1
	subpd   %xmm13, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm9
	movdqu  %xmm9,-14 * SIZE(A1,LDA,1)
	
	MOVUPS_A1(-14 * SIZE, A2, %xmm10)
        movapd  %xmm10,  %xmm0
	subpd   %xmm13, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm10
	movdqu  %xmm10,-14 * SIZE(A2)
	
	MOVUPS_A2(-14 * SIZE, A2, LDA, 1, %xmm11)
	movapd  %xmm11,  %xmm1
	subpd   %xmm13, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm11
	movdqu  %xmm11,-14 * SIZE(A2,LDA,1)

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, A2
	addq	$4 * SIZE, X1
	ALIGN_4

.L26:
	testq	$2, M 
	jle	.L27

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)

	MOVUPS_A1(-16 * SIZE, A1, %xmm8)
	movapd  %xmm8,  %xmm0
	subpd   %xmm12, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm8
	movdqu  %xmm8,-16 * SIZE(A1)
	
	MOVUPS_A2(-16 * SIZE, A1, LDA, 1, %xmm9)
	movapd  %xmm9,  %xmm1
	subpd   %xmm12, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm9
	movdqu  %xmm9,-16 * SIZE(A1,LDA,1)
	
	MOVUPS_A1(-16 * SIZE, A2, %xmm10)
	movapd  %xmm10,  %xmm0
	subpd   %xmm12, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm10
	movdqu  %xmm10,-16 * SIZE(A2)
	MOVUPS_A2(-16 * SIZE, A2, LDA, 1, %xmm11)
	movapd  %xmm11,  %xmm1
	subpd   %xmm12, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm11
	movdqu  %xmm11,-16 * SIZE(A2,LDA,1)

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, A2
	addq	$2 * SIZE, X1
	ALIGN_4

.L27:
	testq	$1, M  /*ONE COL*/
	je	.L28

	movsd	-16 * SIZE(X1), %xmm12

	movsd	-16 * SIZE(A1), %xmm8
	movsd   %xmm8,  %xmm0
	subsd   %xmm12, %xmm0
	mulsd   ALPHA, %xmm0
	addsd   %xmm0, %xmm8
	movsd  %xmm8,-16 * SIZE(A1)
	
	movsd	-16 * SIZE(A1, LDA), %xmm9
	movsd   %xmm9,  %xmm1
	subsd   %xmm12, %xmm1
	mulsd   ALPHA, %xmm1
	addsd   %xmm1, %xmm9
	movsd  %xmm9,-16 * SIZE(A1,LDA)
	movsd	-16 * SIZE(A2), %xmm10
	movsd   %xmm10,  %xmm0
	subsd   %xmm12, %xmm0
	mulsd   ALPHA, %xmm0
	addsd   %xmm0, %xmm10
	movsd  %xmm10,-16 * SIZE(A2)
	movsd	-16 * SIZE(A2, LDA), %xmm11
	movsd   %xmm11,  %xmm1
	subsd   %xmm12, %xmm1
	mulsd   ALPHA, %xmm1
	addsd   %xmm1, %xmm11
	movsd  %xmm11,-16 * SIZE(A2,LDA)
	ALIGN_4

.L28:
	cmpq	$4, N    /*LOOP FOR ROWS*/
	jge	.L21
	
.L30:
	cmpq	$2, N   
	jl	.L40
	
	subq	$2, N  

	movq	X, X1  
	movq	A, A1
	leaq	(A1, LDA), A2
	leaq	(A1, LDA, 2), A
	
	movq	M,  I 
	sarq	$3,  I
	jle	.L35

	MOVUPS_A1(-16 * SIZE, A1, %xmm8)    /*2x32*/
	MOVUPS_A1(-16 * SIZE, A2, %xmm9)
	MOVUPS_A1(-14 * SIZE, A1, %xmm10)
	MOVUPS_A1(-14 * SIZE, A2, %xmm11)

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12) 
	MOVUPS_XL1(-14 * SIZE, X1, %xmm13)

	decq	I
	jle	.L33
	ALIGN_4

	
.L32:
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(A1)
#endif
	
	movapd  %xmm8,  %xmm0
	subpd   %xmm12, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm8
	movdqu  %xmm8,-16 * SIZE(A1)
	
	MOVUPS_A1(-12 * SIZE, A1, %xmm8)
	movapd  %xmm9,  %xmm1
	subpd   %xmm12, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm9
	movdqu  %xmm9,-16 * SIZE(A2)
	
	MOVUPS_XL1(-12 * SIZE, X1, %xmm12)
	MOVUPS_A1(-12 * SIZE, A2, %xmm9)
	movapd  %xmm10,  %xmm0
	subpd   %xmm13, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm10
	movdqu  %xmm10,-14 * SIZE(A1)
	MOVUPS_A1(-10 * SIZE, A1, %xmm10)
	movapd  %xmm11,  %xmm1
	subpd   %xmm13, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm11
	movdqu  %xmm11,-14 * SIZE(A2)     
	MOVUPS_XL1(-10 * SIZE, X1, %xmm13)
	MOVUPS_A1(-10 * SIZE, A2, %xmm11)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(A2)
#endif

	movapd  %xmm8,  %xmm0
	subpd   %xmm12, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm8
	movdqu  %xmm8, -12 * SIZE(A1)
	MOVUPS_A1( -8 * SIZE, A1, %xmm8)
	movapd  %xmm9,  %xmm1
	subpd   %xmm12, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm9
	movdqu  %xmm9, -12 * SIZE(A2)
	
	MOVUPS_XL1( -8 * SIZE, X1, %xmm12)
	MOVUPS_A1( -8 * SIZE, A2, %xmm9)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(X1)
#endif

	movapd  %xmm10,  %xmm0
	subpd   %xmm13, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm10
	movdqu  %xmm10, -10 * SIZE(A1)
	MOVUPS_A1( -6 * SIZE, A1, %xmm10)
	movapd  %xmm11, %xmm1
	subpd   %xmm13, %xmm1
	mulpd   ALPHA,  %xmm1
	addpd   %xmm1,  %xmm11
	movdqu  %xmm11, -10 * SIZE(A2)
	MOVUPS_XL1( -6 * SIZE, X1, %xmm13)
	MOVUPS_A1( -6 * SIZE, A2, %xmm11)

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1

	decq	I
	jg	.L32
	ALIGN_4

.L33:

	movapd  %xmm8,  %xmm0
	subpd   %xmm12, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm8
	movdqu  %xmm8,-16 * SIZE(A1)
	MOVUPS_A1(-12 * SIZE, A1, %xmm8)
	movapd  %xmm9,  %xmm1
	subpd   %xmm12, %xmm1
	mulpd   ALPHA,  %xmm1
	addpd   %xmm1,  %xmm9
	movdqu  %xmm9, -16 * SIZE(A2)
	MOVUPS_XL1(-12 * SIZE, X1, %xmm12)
	MOVUPS_A1(-12 * SIZE, A2, %xmm9)

	movapd  %xmm10,  %xmm0
	subpd   %xmm13, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm10
	movdqu  %xmm10, -14 * SIZE(A1)
	MOVUPS_A1(-10 * SIZE, A1, %xmm10)
	movapd  %xmm11,  %xmm1
	subpd   %xmm13, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm11
	movdqu  %xmm11, -14 * SIZE(A2)
	MOVUPS_XL1(-10 * SIZE, X1, %xmm13)
	MOVUPS_A1(-10 * SIZE, A2, %xmm11)

	movapd  %xmm8,  %xmm0
	subpd   %xmm12, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm8
	movdqu  %xmm8, -12 * SIZE(A1)

	movapd  %xmm9,  %xmm1
	subpd   %xmm12, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm9
	movdqu %xmm9,   -12 * SIZE(A2)

	movapd  %xmm10,  %xmm0
	subpd   %xmm13, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm10
	movdqu  %xmm10, -10 * SIZE(A1)

	movapd  %xmm11,  %xmm1
	subpd   %xmm13, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm11
	movdqu  %xmm11, -10 * SIZE(A2)

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1
	ALIGN_4

.L35:	
	testq	$4, M
	jle	.L36

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)

	MOVUPS_A1(-16 * SIZE, A1, %xmm8)
	movapd  %xmm8,  %xmm0
	subpd   %xmm12, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm8
	movdqu  %xmm8,-16 * SIZE(A1)
	MOVUPS_A1(-16 * SIZE, A2, %xmm9)
	movapd  %xmm9,  %xmm1
	subpd   %xmm12, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm9
	movdqu  %xmm9,-16 * SIZE(A2)

	MOVUPS_XL1(-14 * SIZE, X1, %xmm13)

	MOVUPS_A1(-14 * SIZE, A1, %xmm10)
	movapd  %xmm10,  %xmm0
	subpd   %xmm13, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm10
	movdqu  %xmm10,-14 * SIZE(A1)
	MOVUPS_A1(-14 * SIZE, A2, %xmm11)
	movapd  %xmm11,  %xmm1
	subpd   %xmm13, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm11
	movdqu  %xmm11,-14 * SIZE(A2)

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, A2
	addq	$4 * SIZE, X1
	ALIGN_4

.L36:
	
	testq	$2, M
	jle	.L37

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12) 

	MOVUPS_A1(-16 * SIZE, A1, %xmm8)   
	movapd  %xmm8,  %xmm0
	subpd   %xmm12, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm8
	movdqu  %xmm8,-16 * SIZE(A1)
	MOVUPS_A1(-16 * SIZE, A2, %xmm9)   
	movapd  %xmm9,  %xmm1
	subpd   %xmm12, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm9
	movdqu  %xmm9,-16 * SIZE(A2)

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, A2
	addq	$2 * SIZE, X1
	ALIGN_4

.L37:
	
	testq	$1, M
	je	.L40

	movsd	-16 * SIZE(X1), %xmm12

	movsd	-16 * SIZE(A1), %xmm8
	movsd   %xmm8,  %xmm0
	subsd   %xmm12, %xmm0
	mulsd   ALPHA, %xmm0
	addsd   %xmm0, %xmm8
	movsd   %xmm8, -16 * SIZE(A1)
	movsd	-16 * SIZE(A2), %xmm9
	movsd   %xmm9,  %xmm1
	subsd   %xmm12, %xmm1
	mulsd   ALPHA, %xmm1
	addsd   %xmm1, %xmm9
	movsd   %xmm9, -16 * SIZE(A2)
	ALIGN_4
	
.L40:
	cmpq	$1, N  /*ROW 1*/
	jl	.L999

	movq	X, X1  
	movq	A, A1

	movq	M,  I
	sarq	$3,  I
	jle	.L45

	MOVUPS_A1(-16 * SIZE, A1, %xmm8)
	MOVUPS_A1(-14 * SIZE, A1, %xmm9)
	MOVUPS_A1(-12 * SIZE, A1, %xmm10)
	MOVUPS_A1(-10 * SIZE, A1, %xmm11)

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)
	MOVUPS_XL1(-14 * SIZE, X1, %xmm13)

	decq	I
	jle	.L43
	ALIGN_4

.L42:
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 8 - 128 + PREOFFSET(A1)
#endif

	movapd  %xmm8,  %xmm0
	subpd   %xmm12, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm8
	movdqu  %xmm8,-16 * SIZE(A1)
	
	MOVUPS_XL1(-12 * SIZE, X1, %xmm12)
	MOVUPS_A1( -8 * SIZE, A1, %xmm8)

	movapd  %xmm9,  %xmm1
	subpd   %xmm13, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm9
	movdqu  %xmm9,-14 * SIZE(A1)
	
	MOVUPS_XL1(-10 * SIZE, X1, %xmm13)
	MOVUPS_A1( -6 * SIZE, A1, %xmm9)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 8 - 128 + PREOFFSET(X1)
#endif

	movapd  %xmm10,  %xmm0
	subpd   %xmm12, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm10
	movdqu  %xmm10,-12 * SIZE(A1)
	
	MOVUPS_XL1( -8 * SIZE, X1, %xmm12)
	MOVUPS_A1( -4 * SIZE, A1, %xmm10)
	movapd  %xmm11,  %xmm1
	subpd   %xmm13, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm11
	movdqu  %xmm11,-10 * SIZE(A1)

	MOVUPS_XL1( -6 * SIZE, X1, %xmm13)
	MOVUPS_A1( -2 * SIZE, A1, %xmm11)

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, X1

	decq	I
	jg	.L42
	ALIGN_4

.L43:
	movapd  %xmm8,  %xmm0
	subpd   %xmm12, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm8
	movdqu  %xmm8,-16 * SIZE(A1)
	MOVUPS_XL1(-12 * SIZE, X1, %xmm12)

	movapd  %xmm9,  %xmm1
	subpd   %xmm13, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm9
	movdqu  %xmm9,-14 * SIZE(A1)
	MOVUPS_XL1(-10 * SIZE, X1, %xmm13)

	movapd  %xmm10,  %xmm0
	subpd   %xmm12, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm10
	movdqu  %xmm10,-12 * SIZE(A1)
	movapd  %xmm11,  %xmm1
	subpd   %xmm13, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm11
	movdqu  %xmm11,-10 * SIZE(A1)

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, X1
	ALIGN_4

.L45:
	testq	$4, M     /*CONTINUE ROW = 4*/
	jle	.L46

	MOVUPS_A1(-16 * SIZE, A1, %xmm8)
	MOVUPS_A1(-14 * SIZE, A1, %xmm9)

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)
	MOVUPS_XL1(-14 * SIZE, X1, %xmm13)

	movapd  %xmm8,  %xmm0
	subpd   %xmm12, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm8
 	
	movapd  %xmm9,  %xmm1
	subpd   %xmm13, %xmm1
	mulpd   ALPHA, %xmm1
	addpd   %xmm1, %xmm9

	movdqu  %xmm8,-16 * SIZE(A1)
	movdqu  %xmm9,-14 * SIZE(A1)
	
	addq	$4 * SIZE, A1
	addq	$4 * SIZE, X1
	ALIGN_4

.L46:
	testq	$2, M  /*LAST ROW*/
	jle	.L47

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12) /*movups -0x80(%r15),%xmm12*/

	MOVUPS_A1(-16 * SIZE, A1, %xmm8)   /*movups -0x80(%rbx),%xmm8*/
	movapd  %xmm8,  %xmm0
	subpd   %xmm12, %xmm0
	mulpd   ALPHA, %xmm0
	addpd   %xmm0, %xmm8
	movdqu  %xmm8,-16 * SIZE(A1)

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, X1
	ALIGN_4

.L47:
	testq	$1, M
	je	.L999

	movsd	-16 * SIZE(X1), %xmm12   /*DATA LAST*/

	movsd	-16 * SIZE(A1), %xmm8
	movsd   %xmm8,  %xmm0
	subsd   %xmm12, %xmm0
	mulsd   ALPHA, %xmm0
	addsd   %xmm0, %xmm8
	movsd   %xmm8, -16 * SIZE(A1)
	ALIGN_4
	
.L999:
		
	movq	  0(%rsp), %rbx
	movq	  8(%rsp), %rbp
	movq	 16(%rsp), %r12
	movq	 24(%rsp), %r13
	movq	 32(%rsp), %r14
	movq	 40(%rsp), %r15

#ifdef WINDOWS_ABI
	movq	 48(%rsp), %rdi
	movq	 56(%rsp), %rsi
	movups	 64(%rsp), %xmm6
	movups	 80(%rsp), %xmm7
	movups	 96(%rsp), %xmm8
	movups	112(%rsp), %xmm9
	movups	128(%rsp), %xmm10
	movups	144(%rsp), %xmm11
	movups	160(%rsp), %xmm12
	movups	176(%rsp), %xmm13
	movups	192(%rsp), %xmm14
	movups	208(%rsp), %xmm15
#endif

	addq	$STACKSIZE, %rsp
	ret
	ALIGN_4

.LFE2:
	.size	daasxpa_kh, .-daasxpa_kh
	.section	.eh_frame,"a",@progbits
.Lframe1:
	.long	.LECIE1-.LSCIE1
.LSCIE1:
	.long	0x0
	.byte	0x1
	.string	"zR"
	.uleb128 0x1
	.sleb128 -8
	.byte	0x10
	.uleb128 0x1
	.byte	0x3
	.byte	0xc
	.uleb128 0x7
	.uleb128 0x8
	.byte	0x90
	.uleb128 0x1
	.align 8
.LECIE1:
.LSFDE1:
	.long	.LEFDE1-.LASFDE1
.LASFDE1:
	.long	.LASFDE1-.Lframe1
	.long	.LFB2
	.long	.LFE2-.LFB2
	.uleb128 0x0
	.align 8
.LEFDE1:
	.ident	"GCC: (Debian 4.3.2-1.1) 4.3.2"
	.section	.note.GNU-stack,"",@progbits
	
/*	EPILOGUE
*/
