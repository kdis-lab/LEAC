/*! \file iaysxpy_kh.S	
 * This file is part of the LEAC.
 *
 * Implementation of the 
 *	
 *	y = y + a (y -x),
 *
 *  function for integers
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

#ifndef WINDOWS_ABI
#define M	ARG1 /*%rdi*/
#define X	ARG4
#define INCX	ARG5
#define Y	ARG6
#define INCY	ARG2
#else
#define M	ARG1
#define X	ARG2
#define INCX	ARG3
#define Y	ARG4
#define INCY	%r10
#endif

#define	YY	%r11
#define ALPHA	%xmm15

#define SIZE    4
	
#include "l1param.h"

	.file	 "iaysxpy_kh.S"
	.text
	.p2align 4,,15
	.globl   iaysxpy_kh
	.type	 iaysxpy_kh, @function
iaysxpy_kh:
.LFB2:

#ifndef WINDOWS_ABI
	movq	 8(%rsp), INCY
	movaps	%xmm0,  ALPHA
#else
	movaps	%xmm3,  ALPHA
	movq	40(%rsp), X
	movq	48(%rsp), INCX
	movq	56(%rsp), Y
	movq	64(%rsp), INCY
#endif

	SAVEREGISTERS

/*CHANGE SIGNED*/
	movsd	.LC0(%rip), %xmm1
	xorpd	%xmm1, ALPHA
	unpcklpd ALPHA, ALPHA
	
	leaq	(, INCX, SIZE), INCX
	leaq	(, INCY, SIZE), INCY
	
	cmpq	$SIZE, INCX
	jne	.L50
	cmpq	$SIZE, INCY
	jne	.L50

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y

	cmpq	$3, M
	jle	.L16

	testq	$SIZE, Y
	je	.L05

	movl   -32 * SIZE(X), %r13d
	movl   -32 * SIZE(Y), %r14d
	subl   %r14d, %r13d
	cvtsi2sd %r13d, %xmm1
	mulsd	ALPHA, %xmm1
	cvttsd2si %xmm1, %r13d
	addl   %r13d, %r14d
	movl   %r14d, -32 * SIZE(Y)
	
	addq	$1 * SIZE, X
	addq	$1 * SIZE, Y
	decq	M
	jle	.L19
	ALIGN_3

.L05:
	testq	$SIZE * 2, Y
	je	.L10

	movq	-32 * SIZE(X), %mm0
	movq	-32 * SIZE(Y), %mm1
	psubd   %mm1, %mm0
	cvtpi2pd %mm0, %xmm2
	mulpd   ALPHA, %xmm2
	cvtpd2pi %xmm2, %mm2
	paddd	%mm2, %mm1
	movq	%mm1, -32 * SIZE(Y)
	
	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	subq	$2, M
	jle	.L19
	ALIGN_3

.L10:
	cmpq    $3, M
	jle     .L16

	testq   $2 * SIZE, X
	jne     .L30

	testq   $1 * SIZE, X
	jne     .L20
	
	movq    M,  %rax
	sarq    $5, %rax
	jle     .L13
	ALIGN_3
	
.L11:
#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif
	
	movdqa	-32 * SIZE(X), %xmm0
	movdqa	-32 * SIZE(Y), %xmm1

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm2
	cvtpi2pd %mm1, %xmm3
	
	mulpd	ALPHA, %xmm2
	mulpd	ALPHA, %xmm3
	cvtpd2pi %xmm2, %mm0
	cvtpd2pi %xmm3, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm3
	punpcklqdq %xmm3, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)

	movdqa	-28 * SIZE(X), %xmm0
	movdqa	-28 * SIZE(Y), %xmm1

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm2
	cvtpi2pd %mm1, %xmm3
	
	mulpd	ALPHA, %xmm2
	mulpd	ALPHA, %xmm3
	cvtpd2pi %xmm2, %mm0
	cvtpd2pi %xmm3, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm3
	punpcklqdq %xmm3, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -28 * SIZE(Y)

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif
	
	movdqa	-24 * SIZE(X), %xmm0
	movdqa	-24 * SIZE(Y), %xmm1

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm2
	cvtpi2pd %mm1, %xmm3
	
	mulpd	ALPHA, %xmm2
	mulpd	ALPHA, %xmm3
	cvtpd2pi %xmm2, %mm0
	cvtpd2pi %xmm3, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm3
	punpcklqdq %xmm3, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -24 * SIZE(Y)
	
	movdqa	-20 * SIZE(X), %xmm0
	movdqa	-20 * SIZE(Y), %xmm1
	
	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm2
	cvtpi2pd %mm1, %xmm3
	
	mulpd	ALPHA, %xmm2
	mulpd	ALPHA, %xmm3
	cvtpd2pi %xmm2, %mm0
	cvtpd2pi %xmm3, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm3
	punpcklqdq %xmm3, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -20 * SIZE(Y)

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE + 64) - PREOFFSET(X)
#endif
	
	movdqa	-16 * SIZE(X), %xmm0
	movdqa	-16 * SIZE(Y), %xmm1
	
	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm2
	cvtpi2pd %mm1, %xmm3
	
	mulpd	ALPHA, %xmm2
	mulpd	ALPHA, %xmm3
	cvtpd2pi %xmm2, %mm0
	cvtpd2pi %xmm3, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm3
	punpcklqdq %xmm3, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -16 * SIZE(Y)

	movdqa	-12 * SIZE(X), %xmm0
	movdqa	-12 * SIZE(Y), %xmm1
	
	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm2
	cvtpi2pd %mm1, %xmm3
	
	mulpd	ALPHA, %xmm2
	mulpd	ALPHA, %xmm3
	cvtpd2pi %xmm2, %mm0
	cvtpd2pi %xmm3, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm3
	punpcklqdq %xmm3, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -12 * SIZE(Y)

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE + 64) - PREOFFSET(Y)
#endif

	movdqa	-8 * SIZE(X), %xmm0
	movdqa	-8 * SIZE(Y), %xmm1
	
	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm2
	cvtpi2pd %mm1, %xmm3
	
	mulpd	ALPHA, %xmm2
	mulpd	ALPHA, %xmm3
	cvtpd2pi %xmm2, %mm0
	cvtpd2pi %xmm3, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm3
	punpcklqdq %xmm3, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -8 * SIZE(Y)

	movdqa	-4 * SIZE(X), %xmm0
	movdqa	-4 * SIZE(Y), %xmm1
	
	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm2
	cvtpi2pd %mm1, %xmm3
	
	mulpd	ALPHA, %xmm2
	mulpd	ALPHA, %xmm3
	cvtpd2pi %xmm2, %mm0
	cvtpd2pi %xmm3, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm3
	punpcklqdq %xmm3, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -4 * SIZE(Y)

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	
	decq	%rax
	jg	.L11
	ALIGN_3
	
.L13:
	testq   $16, M
	jle     .L14

	movdqa	-32 * SIZE(X), %xmm0
	movdqa	-32 * SIZE(Y), %xmm1

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm2
	cvtpi2pd %mm1, %xmm3
	
	mulpd	ALPHA, %xmm2
	mulpd	ALPHA, %xmm3
	cvtpd2pi %xmm2, %mm0
	cvtpd2pi %xmm3, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm3
	punpcklqdq %xmm3, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)

	movdqa	-28 * SIZE(X), %xmm0
	movdqa	-28 * SIZE(Y), %xmm1

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm2
	cvtpi2pd %mm1, %xmm3
	
	mulpd	ALPHA, %xmm2
	mulpd	ALPHA, %xmm3
	cvtpd2pi %xmm2, %mm0
	cvtpd2pi %xmm3, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm3
	punpcklqdq %xmm3, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -28 * SIZE(Y)
	
	movdqa	-24 * SIZE(X), %xmm0
	movdqa	-24 * SIZE(Y), %xmm1

        psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm2
	cvtpi2pd %mm1, %xmm3
	
	mulpd	ALPHA, %xmm2
	mulpd	ALPHA, %xmm3
	cvtpd2pi %xmm2, %mm0
	cvtpd2pi %xmm3, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm3
	punpcklqdq %xmm3, %xmm0
	
	paddd	%xmm1, %xmm0
	
	movdqa	%xmm0, -24 * SIZE(Y)

	movdqa	-20 * SIZE(X), %xmm0
	movdqa	-20 * SIZE(Y), %xmm1
	
	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm2
	cvtpi2pd %mm1, %xmm3
	
	mulpd	ALPHA, %xmm2
	mulpd	ALPHA, %xmm3
	cvtpd2pi %xmm2, %mm0
	cvtpd2pi %xmm3, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm3
	punpcklqdq %xmm3, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -20 * SIZE(Y)

	addq    $16 * SIZE, X
	addq    $16 * SIZE, Y
	ALIGN_3
	
.L14:
	movq	M,  %rax
	andq	$8, %rax
	jle	.L15
	ALIGN_3

	movdqa	-32 * SIZE(X), %xmm0
	movdqa	-32 * SIZE(Y), %xmm1

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm2
	cvtpi2pd %mm1, %xmm3
	
	mulpd	ALPHA, %xmm2
	mulpd	ALPHA, %xmm3
	cvtpd2pi %xmm2, %mm0
	cvtpd2pi %xmm3, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm3
	punpcklqdq %xmm3, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)
	
	movdqa	-28 * SIZE(X), %xmm0
	movdqa	-28 * SIZE(Y), %xmm1

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm2
	cvtpi2pd %mm1, %xmm3
	
	mulpd	ALPHA, %xmm2
	mulpd	ALPHA, %xmm3
	cvtpd2pi %xmm2, %mm0
	cvtpd2pi %xmm3, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm3
	punpcklqdq %xmm3, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -28 * SIZE(Y)
	
	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L15:
	movq	M,  %rax
	andq	$4, %rax
	jle	.L16
	ALIGN_3

	movdqa	-32 * SIZE(X), %xmm0
	movdqa	-32 * SIZE(Y), %xmm1

	psubd   %xmm1, %xmm0
	
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm2
	cvtpi2pd %mm1, %xmm3
	
	mulpd	ALPHA, %xmm2
	mulpd	ALPHA, %xmm3
	cvtpd2pi %xmm2, %mm0
	cvtpd2pi %xmm3,	%mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm3
	punpcklqdq %xmm3, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)
	
	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L16:
	movq	M,  %rax
	andq	$2, %rax
	jle	.L17
	ALIGN_3

	movq	-32 * SIZE(X), %mm0
	movq	-32 * SIZE(Y), %mm1

	psubd   %mm1, %mm0

	cvtpi2pd %mm0, %xmm2
	mulpd   ALPHA, %xmm2
	cvtpd2pi %xmm2, %mm2
	
	paddd	%mm2, %mm1

	movq	%mm1, -32 * SIZE(Y)

	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

.L17:
	movq	M,  %rax
	andq	$1, %rax
	jle	.L19
	ALIGN_3

	movl   -32 * SIZE(X), %r13d
	movl   -32 * SIZE(Y), %r14d

	subl   %r14d, %r13d
	cvtsi2sd %r13d, %xmm1
	
	mulsd	ALPHA, %xmm1
	cvttsd2si %xmm1, %r13d
	addl %r13d, %r14d
	movl   %r14d, -32 * SIZE(Y)

	ALIGN_3

.L19:
	xorq	%rax,%rax

	RESTOREREGISTERS

#ifdef WINDOWS_ABI
	popq	%rbx
#endif
	
	ret
	ALIGN_3

.L20:

#ifdef ALIGNED_ACCESS


	movdqa   -33 * SIZE(X), %xmm0
	movdqa   -32 * SIZE(Y), %xmm1
	
	movq    M,  %rax
	sarq    $5, %rax
	jle     .L23
	ALIGN_4

.L21:
#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movdqa  -29 * SIZE(X), %xmm2
	movdqa  -28 * SIZE(Y), %xmm3
	
	psrldq  $4, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $12, %xmm8
	por     %xmm8, %xmm0

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)

	movdqa  -25 * SIZE(X), %xmm0
	movdqa  -24 * SIZE(Y), %xmm1

	psrldq  $4, %xmm2  /* X */
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $12, %xmm8
	por     %xmm8, %xmm2

	psubd   %xmm3, %xmm2
	movdq2q  %xmm2, %mm0
	psrldq   $8,    %xmm2
	movdq2q  %xmm2, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm2
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm2
	
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, -28 * SIZE(Y)

#ifdef PREFETCH
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	movdqa  -21 * SIZE(X), %xmm2
	movdqa  -20 * SIZE(Y), %xmm3
	
	psrldq  $4, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $12, %xmm8
	por     %xmm8, %xmm0

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -24 * SIZE(Y)

	movdqa  -17 * SIZE(X), %xmm0
	movdqa  -16 * SIZE(Y), %xmm1

	psrldq  $4, %xmm2
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $12, %xmm8
	por     %xmm8, %xmm2

	psubd   %xmm3, %xmm2
	movdq2q  %xmm2, %mm0
	psrldq   $8,    %xmm2
	movdq2q  %xmm2, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm2
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm2
	
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, -20 * SIZE(Y)

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE + 64) - PREOFFSET(X)
#endif
	
	movdqa  -13 * SIZE(X), %xmm2
	movdqa  -12 * SIZE(Y), %xmm3
	
	psrldq  $4, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $12, %xmm8
	por     %xmm8, %xmm0

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -16 * SIZE(Y)

	movdqa   -9 * SIZE(X), %xmm0
	movdqa   -8 * SIZE(Y), %xmm1

	psrldq  $4, %xmm2
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $12, %xmm8
	por     %xmm8, %xmm2

	psubd   %xmm3, %xmm2
	movdq2q  %xmm2, %mm0
	psrldq   $8,    %xmm2
	movdq2q  %xmm2, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm2
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm2
	
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, -12 * SIZE(Y)

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE + 64) - PREOFFSET(Y)
#endif

	movdqa   -5 * SIZE(X), %xmm2
	movdqa   -4 * SIZE(Y), %xmm3
	
	psrldq  $4, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $12, %xmm8
	por     %xmm8, %xmm0

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -8 * SIZE(Y)

	movdqa   -1 * SIZE(X), %xmm0
	movdqa    0 * SIZE(Y), %xmm1
	
	psrldq  $4, %xmm2
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $12, %xmm8
	por     %xmm8, %xmm2

	psubd   %xmm3, %xmm2
	movdq2q  %xmm2, %mm0
	psrldq   $8,    %xmm2
	movdq2q  %xmm2, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm2
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm2
	
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, -4  * SIZE(Y)
	
	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	decq	%rax
	jg	.L21
	ALIGN_3

.L23:

	testq   $16, M
        jle     .L24

        movdqa  -29 * SIZE(X), %xmm2
	movdqa  -28 * SIZE(Y), %xmm3

	psrldq  $4, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $12, %xmm8
	por     %xmm8, %xmm0

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)
	
	movdqa  -25 * SIZE(X), %xmm0
	movdqa  -24 * SIZE(Y), %xmm1

	psrldq  $4, %xmm2  /* X */
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $12, %xmm8
	por     %xmm8, %xmm2

	psubd   %xmm3, %xmm2
	movdq2q  %xmm2, %mm0
	psrldq   $8,    %xmm2
	movdq2q  %xmm2, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm2
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm2
	
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, -28 * SIZE(Y)

	movdqa  -21 * SIZE(X), %xmm2
	movdqa  -20 * SIZE(Y), %xmm3

	psrldq  $4, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $12, %xmm8
	por     %xmm8, %xmm0

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -24 * SIZE(Y)

	movdqa  -17 * SIZE(X), %xmm0
	movdqa  -16 * SIZE(Y), %xmm1

	psrldq  $4, %xmm2
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $12, %xmm8
	por     %xmm8, %xmm2

	psubd   %xmm3, %xmm2
	movdq2q  %xmm2, %mm0
	psrldq   $8,    %xmm2
	movdq2q  %xmm2, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm2
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm2
	
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, -20 * SIZE(Y)

	addq    $16 * SIZE, X
        addq    $16 * SIZE, Y
	ALIGN_3

.L24:
	testq   $8, M
	jle     .L25

	movdqa  -29 * SIZE(X), %xmm2
	movdqa  -28 * SIZE(Y), %xmm3

	psrldq  $4, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $12, %xmm8
	por     %xmm8, %xmm0

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)
	
	movdqa  -25 * SIZE(X), %xmm0
	movdqa  -24 * SIZE(Y), %xmm1

	psrldq  $4, %xmm2  /* X */
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $12, %xmm8
	por     %xmm8, %xmm2

	psubd   %xmm3, %xmm2
	movdq2q  %xmm2, %mm0
	psrldq   $8,    %xmm2
	movdq2q  %xmm2, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm2
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm2
	
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, -28 * SIZE(Y)

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L25:
	testq   $4, M
	jle     .L26

	movdqa  -29 * SIZE(X), %xmm2

	psrldq  $4, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $12, %xmm8
	por     %xmm8, %xmm0

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)
	
	movdqa  %xmm2, %xmm0
	movdqa  -28 * SIZE(Y), %xmm1
	
	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3
	
.L26:
	movq	M,  %rax
	andq	$2, %rax
	jle	.L27
	ALIGN_3

	psrldq  $4, %xmm0
	movdq2q	%xmm0, %mm0
	movdq2q	%xmm1, %mm1
	
	psubd   %mm1, %mm0

	cvtpi2pd %mm0, %xmm2
	mulpd   ALPHA, %xmm2
	cvtpd2pi %xmm2, %mm2
	
	paddd	%mm2, %mm1

	movq	%mm1, -32 * SIZE(Y)
	
	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

.L27:
	movq	M,  %rax
	andq	$1, %rax
	jle	.L29
	ALIGN_3

	movl   -32 * SIZE(X),%r13d
	movl   -32 * SIZE(Y), %r14d

	subl   %r14d, %r13d
	cvtsi2sd %r13d, %xmm1
	
	mulsd	ALPHA, %xmm1
	cvttsd2si %xmm1, %r13d
	addl %r13d, %r14d
	movl   %r14d, -32 * SIZE(Y)

	ALIGN_3

.L29:
	xorq	%rax,%rax

	RESTOREREGISTERS

#ifdef WINDOWS_ABI
	popq    %rbx
#endif
	ret
	ALIGN_3
	
.L30:
	testq	$1 * SIZE, X
	jne	.L40

	movdqa   -34 * SIZE(X), %xmm0
	movdqa   -32 * SIZE(Y), %xmm1
	
	movq    M,  %rax
	sarq    $5, %rax
	jle     .L33
	ALIGN_4

.L31:
#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movdqa  -30 * SIZE(X), %xmm2
	movdqa  -28 * SIZE(Y), %xmm3

	psrldq   $8, %xmm0
	punpcklqdq %xmm2, %xmm0
	
	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)
	
	movaps	-26 * SIZE(X), %xmm0
	movaps	-24 * SIZE(Y), %xmm1

	psrldq   $8, %xmm2
	punpcklqdq %xmm0, %xmm2
	
	psubd   %xmm3, %xmm2
	movdq2q  %xmm2, %mm0
	psrldq   $8,    %xmm2
	movdq2q  %xmm2, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm2
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm2
	
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, -28 * SIZE(Y)
	
	
#ifdef PREFETCH
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif
	movdqa  -22 * SIZE(X), %xmm2
	movdqa  -20 * SIZE(Y), %xmm3

	psrldq   $8, %xmm0
	punpcklqdq %xmm2, %xmm0

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -24 * SIZE(Y)
	
	movdqa  -18 * SIZE(X), %xmm0
	movdqa  -16 * SIZE(Y), %xmm1

	psrldq   $8, %xmm2
	punpcklqdq %xmm0, %xmm2

	psubd   %xmm3, %xmm2
	movdq2q  %xmm2, %mm0
	psrldq   $8,    %xmm2
	movdq2q  %xmm2, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm2
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm2
	
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, -20 * SIZE(Y)
	
	
#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE + 64) - PREOFFSET(X)
#endif

	movdqa  -14 * SIZE(X), %xmm2
	movdqa  -12 * SIZE(Y), %xmm3

	psrldq   $8, %xmm0
	punpcklqdq %xmm2, %xmm0

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -16 * SIZE(Y)
	
	movdqa  -10 * SIZE(X), %xmm0
	movdqa   -8 * SIZE(Y), %xmm1

	psrldq   $8, %xmm2
	punpcklqdq %xmm0, %xmm2

	psubd   %xmm3, %xmm2
	movdq2q  %xmm2, %mm0
	psrldq   $8,    %xmm2
	movdq2q  %xmm2, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm2
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm2
	
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, -12 * SIZE(Y)
	

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  64) - PREOFFSET(Y)
#endif
	
	movdqa   -6 * SIZE(X), %xmm2
	movdqa   -4 * SIZE(Y), %xmm3

	psrldq   $8, %xmm0
	punpcklqdq %xmm2, %xmm0

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -8 * SIZE(Y)
	
	movdqa   -2 * SIZE(X), %xmm0
	movdqa    0 * SIZE(Y), %xmm1

	psrldq   $8, %xmm2
	punpcklqdq %xmm0, %xmm2

	psubd   %xmm3, %xmm2
	movdq2q  %xmm2, %mm0
	psrldq   $8,    %xmm2
	movdq2q  %xmm2, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm2
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm2
	
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, -4 * SIZE(Y)
		
	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	
	decq	%rax
	jg	.L31
	ALIGN_3

.L33:
	testq   $16, M
	jle     .L34
	
	movdqa  -30 * SIZE(X), %xmm2
	movdqa  -28 * SIZE(Y), %xmm3

	psrldq   $8, %xmm0
	punpcklqdq %xmm2, %xmm0
	
	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)

	movaps	-26 * SIZE(X), %xmm0
	movaps	-24 * SIZE(Y), %xmm1
	
	psrldq   $8, %xmm2
	punpcklqdq %xmm0, %xmm2
	
	psubd   %xmm3, %xmm2
	movdq2q  %xmm2, %mm0
	psrldq   $8,    %xmm2
	movdq2q  %xmm2, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm2
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm2
	
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, -28 * SIZE(Y)
	
	movdqa  -22 * SIZE(X), %xmm2
	movdqa  -20 * SIZE(Y), %xmm3

	psrldq   $8, %xmm0
	punpcklqdq %xmm2, %xmm0

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -24 * SIZE(Y)
	
	movdqa  -18 * SIZE(X), %xmm0
	movdqa  -16 * SIZE(Y), %xmm1

	psrldq   $8, %xmm2
	punpcklqdq %xmm0, %xmm2

	psubd   %xmm3, %xmm2
	movdq2q  %xmm2, %mm0
	psrldq   $8,    %xmm2
	movdq2q  %xmm2, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm2
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm2
	
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, -20 * SIZE(Y)
	

	addq    $16 * SIZE, X
	addq    $16 * SIZE, Y
	ALIGN_3

.L34:
	testq   $8, M
	jle     .L35
	ALIGN_3

	movdqa  -30 * SIZE(X), %xmm2
	movdqa  -28 * SIZE(Y), %xmm3

	psrldq   $8, %xmm0
	punpcklqdq %xmm2, %xmm0
	
	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)
	
	movaps	-26 * SIZE(X), %xmm0
	movaps	-24 * SIZE(Y), %xmm1

	psrldq   $8, %xmm2
	punpcklqdq %xmm0, %xmm2
	
	psubd   %xmm3, %xmm2
	movdq2q  %xmm2, %mm0
	psrldq   $8,    %xmm2
	movdq2q  %xmm2, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm2
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm2
	
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, -28 * SIZE(Y)
	
	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L35:
	testq   $4, M
	jle	.L36
	ALIGN_3

	movdqa  -30 * SIZE(X), %xmm2

	psrldq   $8, %xmm0
	punpcklqdq %xmm2, %xmm0
	
	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)

	pshufd  $0xe4, %xmm2, %xmm0
	movdqa  -28 * SIZE(Y), %xmm1

	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L36:
	movq	M,  %rax
	andq	$2, %rax
	jle	.L37
	ALIGN_3

	psrldq   $8, %xmm0
	movdq2q  %xmm0, %mm0
	movdq2q  %xmm1, %mm1
	
	psubd    %mm1, %mm0
	cvtpi2pd %mm0, %xmm10
	mulpd	ALPHA, %xmm10
	cvtpd2pi %xmm10, %mm0
	
	paddd	%mm1, %mm0
	movq	%mm0, -32 * SIZE(Y)
	
	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

.L37:
	movq	M,  %rax
	andq	$1, %rax
	jle	.L39
	ALIGN_3

	movl   -32 * SIZE(X),%r13d
	movl   -32 * SIZE(Y), %r14d

	subl   %r14d, %r13d
	cvtsi2sd %r13d, %xmm1
	
	mulsd	ALPHA, %xmm1
	cvttsd2si %xmm1, %r13d
	addl %r13d, %r14d
	movl   %r14d, -32 * SIZE(Y)
	
	ALIGN_3

.L39:
	xorq	%rax,%rax

	RESTOREREGISTERS

#ifdef WINDOWS_ABI
	popq    %rbx
#endif
	
	ret
	ALIGN_3

.L40:
	movdqa  -35 * SIZE(X), %xmm0
	movdqa  -32 * SIZE(Y), %xmm1

	movq	M,  %rax
	sarq	$5, %rax
	jle	.L43
	ALIGN_4
	
.L41:
#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movdqa  -31 * SIZE(X), %xmm2
	movdqa  -28 * SIZE(Y), %xmm3

	psrldq  $12, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $4, %xmm8
	por     %xmm8, %xmm0

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)

	movdqa  -27 * SIZE(X), %xmm0
	movdqa  -24 * SIZE(Y), %xmm1

	psrldq  $12, %xmm2  /* X */
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $4, %xmm8
	por     %xmm8, %xmm2

	psubd   %xmm3, %xmm2
	movdq2q  %xmm2, %mm0
	psrldq   $8,    %xmm2
	movdq2q  %xmm2, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm2
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm2
	
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, -28 * SIZE(Y)

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	movdqa  -23 * SIZE(X), %xmm2
	movdqa  -20 * SIZE(Y), %xmm3
	
	psrldq  $12, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $4, %xmm8
	por     %xmm8, %xmm0

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -24 * SIZE(Y)

	movdqa  -19 * SIZE(X), %xmm0
	movdqa  -16 * SIZE(Y), %xmm1

	psrldq  $12, %xmm2
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $4, %xmm8
	por     %xmm8, %xmm2

	psubd   %xmm3, %xmm2
	movdq2q  %xmm2, %mm0
	psrldq   $8,    %xmm2
	movdq2q  %xmm2, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm2
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm2
	
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, -20 * SIZE(Y)

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE + 64) - PREOFFSET(X)
#endif
	
	movdqa  -15 * SIZE(X), %xmm2
	movdqa  -12 * SIZE(Y), %xmm3
	
	psrldq  $12, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $4, %xmm8
	por     %xmm8, %xmm0

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -16 * SIZE(Y)

	movdqa   -11 * SIZE(X), %xmm0
	movdqa   -8 * SIZE(Y), %xmm1

	psrldq  $12, %xmm2
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $4, %xmm8
	por     %xmm8, %xmm2

	psubd   %xmm3, %xmm2
	movdq2q  %xmm2, %mm0
	psrldq   $8,    %xmm2
	movdq2q  %xmm2, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm2
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm2
	
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, -12 * SIZE(Y)

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE + 64) - PREOFFSET(Y)
#endif

	movdqa   -7 * SIZE(X), %xmm2
	movdqa   -4 * SIZE(Y), %xmm3
	
	psrldq  $12, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $4, %xmm8
	por     %xmm8, %xmm0

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -8 * SIZE(Y)

	movdqa   -3 * SIZE(X), %xmm0
	movdqa    0 * SIZE(Y), %xmm1
	
	psrldq  $12, %xmm2
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $4, %xmm8
	por     %xmm8, %xmm2

	psubd   %xmm3, %xmm2
	movdq2q  %xmm2, %mm0
	psrldq   $8,    %xmm2
	movdq2q  %xmm2, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm2
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm2
	
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, -4  * SIZE(Y)
	
	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	
	decq	%rax
	jg	.L41
	ALIGN_3

.L43:
	testq   $16, M
        jle     .L44

        movdqa  -31 * SIZE(X), %xmm2
	movdqa  -28 * SIZE(Y), %xmm3

	psrldq  $12, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $4, %xmm8
	por     %xmm8, %xmm0

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)
	
	movdqa  -27 * SIZE(X), %xmm0
	movdqa  -24 * SIZE(Y), %xmm1

	psrldq  $12, %xmm2  /* X */
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $4, %xmm8
	por     %xmm8, %xmm2

	psubd   %xmm3, %xmm2
	movdq2q  %xmm2, %mm0
	psrldq   $8,    %xmm2
	movdq2q  %xmm2, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm2
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm2
	
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, -28 * SIZE(Y)

	movdqa  -23 * SIZE(X), %xmm2
	movdqa  -20 * SIZE(Y), %xmm3

	psrldq  $12, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $4, %xmm8
	por     %xmm8, %xmm0

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -24 * SIZE(Y)

	movdqa  -19 * SIZE(X), %xmm0
	movdqa  -16 * SIZE(Y), %xmm1

	psrldq  $12, %xmm2
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $4, %xmm8
	por     %xmm8, %xmm2

	psubd   %xmm3, %xmm2
	movdq2q  %xmm2, %mm0
	psrldq   $8,    %xmm2
	movdq2q  %xmm2, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm2
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm2
	
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, -20 * SIZE(Y)

	addq    $16 * SIZE, X
        addq    $16 * SIZE, Y
	ALIGN_3
	
.L44:
	testq   $8, M
	jle     .L45

	movdqa  -31 * SIZE(X), %xmm2
	movdqa  -28 * SIZE(Y), %xmm3

	psrldq  $12, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $4, %xmm8
	por     %xmm8, %xmm0

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)
	
	movdqa  -27 * SIZE(X), %xmm0
	movdqa  -24 * SIZE(Y), %xmm1

	psrldq  $12, %xmm2  /* X */
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $4, %xmm8
	por     %xmm8, %xmm2

	psubd   %xmm3, %xmm2
	movdq2q  %xmm2, %mm0
	psrldq   $8,    %xmm2
	movdq2q  %xmm2, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm2
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm2
	
	paddd	%xmm3, %xmm2
	movdqa	%xmm2, -28 * SIZE(Y)

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L45:
	testq   $4, M
	jle     .L46

	movdqa  -31 * SIZE(X), %xmm2

	psrldq  $12, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $4, %xmm8
	por     %xmm8, %xmm0

	psubd   %xmm1, %xmm0
	movdq2q  %xmm0, %mm0
	psrldq   $8,    %xmm0
	movdq2q  %xmm0, %mm1
	cvtpi2pd %mm0, %xmm10
	cvtpi2pd %mm1, %xmm11
	
	mulpd	ALPHA, %xmm10
	mulpd	ALPHA, %xmm11
	cvtpd2pi %xmm10, %mm0
	cvtpd2pi %xmm11, %mm1

	movq2dq %mm0, %xmm0
	movq2dq %mm1, %xmm11
	punpcklqdq %xmm11, %xmm0
	
	paddd	%xmm1, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)
	
	movdqa  %xmm2, %xmm0
	movdqa  -28 * SIZE(Y), %xmm1
	
	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L46:
	movq	M,  %rax
	andq	$2, %rax
	jle	.L47
	ALIGN_3

	movq    -32 * SIZE(X), %mm0
	movdq2q	%xmm1, %mm1

	psubd   %mm1, %mm0

	cvtpi2pd %mm0, %xmm2
	mulpd   ALPHA, %xmm2
	cvtpd2pi %xmm2, %mm2
	
	paddd	%mm2, %mm1

	movq	%mm1, -32 * SIZE(Y)

	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3
	

.L47:
	movq	M,  %rax
	andq	$1, %rax
	jle	.L49
	ALIGN_3

	movl   -32 * SIZE(X),%r13d
	movl   -32 * SIZE(Y), %r14d

	subl   %r14d, %r13d
	cvtsi2sd %r13d, %xmm1
	
	mulsd	ALPHA, %xmm1
	cvttsd2si %xmm1, %r13d
	addl %r13d, %r14d
	movl   %r14d, -32 * SIZE(Y)
	
	ALIGN_3

.L49:
	xorq	%rax,%rax

	RESTOREREGISTERS

#ifdef WINDOWS_ABI
	popq    %rbx
#endif
	
	ret

#else

	movq	M,  %rax
	sarq	$5, %rax
	jle	.L23

	movsd	-32 * SIZE(X), %xmm0
	movhps	-30 * SIZE(X), %xmm0
	movsd	-28 * SIZE(X), %xmm1
	movhps	-26 * SIZE(X), %xmm1
	movsd	-24 * SIZE(X), %xmm2
	movhps	-22 * SIZE(X), %xmm2
	movsd	-20 * SIZE(X), %xmm3
	movhps	-18 * SIZE(X), %xmm3

	decq	%rax
	jle .L22
	ALIGN_4

.L21:
	movsd	-16 * SIZE(X), %xmm4
	movhps	-14 * SIZE(X), %xmm4
	movsd	-12 * SIZE(X), %xmm5
	movhps	-10 * SIZE(X), %xmm5

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	mulps	ALPHA, %xmm0
	addps	-32 * SIZE(Y), %xmm0
	movaps	%xmm0, -32 * SIZE(Y)

	mulps	ALPHA, %xmm1
	addps	-28 * SIZE(Y), %xmm1
	movaps	%xmm1, -28 * SIZE(Y)

	movsd	 -8 * SIZE(X), %xmm6
	movhps	 -6 * SIZE(X), %xmm6
	movsd	 -4 * SIZE(X), %xmm7
	movhps	 -2 * SIZE(X), %xmm7

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	mulps	ALPHA, %xmm2
	addps	-24 * SIZE(Y), %xmm2
	movaps	%xmm2, -24 * SIZE(Y)

	mulps	ALPHA, %xmm3
	addps	-20 * SIZE(Y), %xmm3
	movaps	%xmm3, -20 * SIZE(Y)

	movsd	  0 * SIZE(X), %xmm0
	movhps	  2 * SIZE(X), %xmm0
	movsd	  4 * SIZE(X), %xmm1
	movhps	  6 * SIZE(X), %xmm1

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  64) - PREOFFSET(Y)
#endif

	mulps	ALPHA, %xmm4
	addps	-16 * SIZE(Y), %xmm4
	movaps	%xmm4, -16 * SIZE(Y)

	mulps	ALPHA, %xmm5
	addps	-12 * SIZE(Y), %xmm5
	movaps	%xmm5, -12 * SIZE(Y)

	movsd	  8 * SIZE(X), %xmm2
	movhps	 10 * SIZE(X), %xmm2
	movsd	 12 * SIZE(X), %xmm3
	movhps	 14 * SIZE(X), %xmm3

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(X)
#endif

	mulps	ALPHA, %xmm6
	addps	 -8 * SIZE(Y), %xmm6
	movaps	%xmm6,  -8 * SIZE(Y)

	mulps	ALPHA, %xmm7
	addps	 -4 * SIZE(Y), %xmm7
	movaps	%xmm7,  -4 * SIZE(Y)

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	decq	%rax
	jg	.L21
	ALIGN_3

.L22:
	movsd	-16 * SIZE(X), %xmm4
	movhps	-14 * SIZE(X), %xmm4
	movsd	-12 * SIZE(X), %xmm5
	movhps	-10 * SIZE(X), %xmm5

	mulps	ALPHA, %xmm0
	addps	-32 * SIZE(Y), %xmm0
	movaps	%xmm0, -32 * SIZE(Y)

	mulps	ALPHA, %xmm1
	addps	-28 * SIZE(Y), %xmm1
	movaps	%xmm1, -28 * SIZE(Y)

	movsd	 -8 * SIZE(X), %xmm6
	movhps	 -6 * SIZE(X), %xmm6
	movsd	 -4 * SIZE(X), %xmm7
	movhps	 -2 * SIZE(X), %xmm7

	mulps	ALPHA, %xmm2
	addps	-24 * SIZE(Y), %xmm2
	movaps	%xmm2, -24 * SIZE(Y)

	mulps	ALPHA, %xmm3
	addps	-20 * SIZE(Y), %xmm3
	movaps	%xmm3, -20 * SIZE(Y)

	mulps	ALPHA, %xmm4
	addps	-16 * SIZE(Y), %xmm4
	movaps	%xmm4, -16 * SIZE(Y)

	mulps	ALPHA, %xmm5
	addps	-12 * SIZE(Y), %xmm5
	movaps	%xmm5, -12 * SIZE(Y)

	mulps	ALPHA, %xmm6
	addps	 -8 * SIZE(Y), %xmm6
	movaps	%xmm6,  -8 * SIZE(Y)

	mulps	ALPHA, %xmm7
	addps	 -4 * SIZE(Y), %xmm7
	movaps	%xmm7,  -4 * SIZE(Y)

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	ALIGN_3

.L23:
	movq	M,  %rax
	andq	$16, %rax
	jle	.L24
	ALIGN_3

	movsd	-32 * SIZE(X), %xmm0
	movhps	-30 * SIZE(X), %xmm0
	movsd	-28 * SIZE(X), %xmm1
	movhps	-26 * SIZE(X), %xmm1

	mulps	ALPHA, %xmm0
	addps	-32 * SIZE(Y), %xmm0
	movaps	%xmm0, -32 * SIZE(Y)
	mulps	ALPHA, %xmm1
	addps	-28 * SIZE(Y), %xmm1
	movaps	%xmm1, -28 * SIZE(Y)

	movsd	-24 * SIZE(X), %xmm2
	movhps	-22 * SIZE(X), %xmm2
	movsd	-20 * SIZE(X), %xmm3
	movhps	-18 * SIZE(X), %xmm3

	mulps	ALPHA, %xmm2
	addps	-24 * SIZE(Y), %xmm2
	movaps	%xmm2, -24 * SIZE(Y)
	mulps	ALPHA, %xmm3
	addps	-20 * SIZE(Y), %xmm3
	movaps	%xmm3, -20 * SIZE(Y)

	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
	ALIGN_3

.L24:
	movq	M,  %rax
	andq	$8, %rax
	jle	.L25
	ALIGN_3

	movsd	-32 * SIZE(X), %xmm0
	movhps	-30 * SIZE(X), %xmm0
	movsd	-28 * SIZE(X), %xmm1
	movhps	-26 * SIZE(X), %xmm1

	mulps	ALPHA, %xmm0
	addps	-32 * SIZE(Y), %xmm0
	mulps	ALPHA, %xmm1
	addps	-28 * SIZE(Y), %xmm1

	movaps	%xmm0, -32 * SIZE(Y)
	movaps	%xmm1, -28 * SIZE(Y)

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L25:
	movq	M,  %rax
	andq	$4, %rax
	jle	.L26
	ALIGN_3

	movsd	-32 * SIZE(X), %xmm0
	movhps	-30 * SIZE(X), %xmm0

	mulps	ALPHA, %xmm0

	addps	-32 * SIZE(Y), %xmm0

	movaps	%xmm0, -32 * SIZE(Y)

	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L26:
	movq	M,  %rax
	andq	$2, %rax
	jle	.L27
	ALIGN_3

	movsd	-32 * SIZE(X), %xmm0
	movsd	-32 * SIZE(Y), %xmm4

	mulps	ALPHA, %xmm0
	addps	%xmm4, %xmm0

	movsd	%xmm0, -32 * SIZE(Y)

	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

.L27:
	movq	M,  %rax
	andq	$1, %rax
	jle	.L29
	ALIGN_3

	movss	-32 * SIZE(X), %xmm0
	mulss	ALPHA, %xmm0
	addss	-32 * SIZE(Y), %xmm0

	movss	%xmm0, -32 * SIZE(Y)
	addq	$SIZE, Y
	ALIGN_3

.L29:
	xorq	%rax,%rax

	RESTOREREGISTERS

#ifdef WINDOWS_ABI
	popq    %rbx
#endif
	ret
#endif
	ALIGN_3


.L50:
	movq	M, %rax
	movq	Y, YY
	sarq	$3,   %rax
	jle	.L55
	ALIGN_3

.L51:
	movl    (X), %r13d
	addq	INCX, X
	movl    (YY), %r14d
	addq	INCY, YY
	subl   %r14d, %r13d
	cvtsi2sd %r13d, %xmm1
	mulsd   ALPHA, %xmm1
	cvttsd2si %xmm1, %r13d
	addl %r13d, %r14d
	movl   %r14d, (Y)
	addq	INCY, Y

	movl    (X), %r13d
	addq	INCX, X
	movl    (YY), %r14d
	addq	INCY, YY
	subl   %r14d, %r13d
	cvtsi2sd %r13d, %xmm1
	mulsd   ALPHA, %xmm1
	cvttsd2si %xmm1, %r13d
	addl %r13d, %r14d
	movl   %r14d, (Y)
	addq	INCY, Y
	
	movl    (X), %r13d
	addq	INCX, X
	movl    (YY), %r14d
	addq	INCY, YY
	subl   %r14d, %r13d
	cvtsi2sd %r13d, %xmm1
	mulsd   ALPHA, %xmm1
	cvttsd2si %xmm1, %r13d
	addl %r13d, %r14d
	movl   %r14d, (Y)
	addq	INCY, Y
	
	movl    (X), %r13d
	addq	INCX, X
	movl    (YY), %r14d
	addq	INCY, YY
	subl   %r14d, %r13d
	cvtsi2sd %r13d, %xmm1
	mulsd   ALPHA, %xmm1
	cvttsd2si %xmm1, %r13d
	addl %r13d, %r14d
	movl   %r14d, (Y)
	addq	INCY, Y
	
	movl    (X), %r13d
	addq	INCX, X
	movl    (YY), %r14d
	addq	INCY, YY
	subl   %r14d, %r13d
	cvtsi2sd %r13d, %xmm1
	mulsd   ALPHA, %xmm1
	cvttsd2si %xmm1, %r13d
	addl %r13d, %r14d
	movl   %r14d, (Y)
	addq	INCY, Y

	movl    (X), %r13d
	addq	INCX, X
	movl    (YY), %r14d
	addq	INCY, YY
	subl   %r14d, %r13d
	cvtsi2sd %r13d, %xmm1
	mulsd   ALPHA, %xmm1
	cvttsd2si %xmm1, %r13d
	addl %r13d, %r14d
	movl   %r14d, (Y)
	addq	INCY, Y

	movl    (X), %r13d
	addq	INCX, X
	movl    (YY), %r14d
	addq	INCY, YY
	subl   %r14d, %r13d
	cvtsi2sd %r13d, %xmm1
	mulsd   ALPHA, %xmm1
	cvttsd2si %xmm1, %r13d
	addl %r13d, %r14d
	movl   %r14d, (Y)
	addq	INCY, Y
	
	movl    (X), %r13d
	addq	INCX, X
	movl    (YY), %r14d
	addq	INCY, YY
	subl   %r14d, %r13d
	cvtsi2sd %r13d, %xmm1
	mulsd   ALPHA, %xmm1
	cvttsd2si %xmm1, %r13d
	addl %r13d, %r14d
	movl   %r14d, (Y)
	addq	INCY, Y

	decq	%rax
	jg	.L51
	ALIGN_3

.L55:
	movq	M, %rax
	andq	$7,   %rax
	jle	.L59
	ALIGN_3

.L56:
	movl    (X), %r13d
	movl    (Y), %r14d
	subl   %r14d, %r13d
	cvtsi2sd %r13d, %xmm1
	mulsd   ALPHA, %xmm1
	cvttsd2si %xmm1, %r13d
	addl %r13d, %r14d
	movl   %r14d, (Y)
	addq	INCX, X
	addq	INCY, Y
	
	decq	%rax
	jg	.L56
	ALIGN_3

.L59:
	xorq	%rax,%rax

	RESTOREREGISTERS

#ifdef WINDOWS_ABI
	popq    %rbx
#endif
	
	ret
	ALIGN_3


/*	EPILOGUE*/
.LFE2:
	.size	iaysxpy_kh, .-iaysxpy_kh
	.section	.eh_frame,"a",@progbits
.LC0:
	.long	0
	.long	-2147483648
	.long	0
	.long	0
	.ident	"GCC: (Debian 4.7.2-5) 4.7.2"
	.section	.note.GNU-stack,"",@progbits
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

