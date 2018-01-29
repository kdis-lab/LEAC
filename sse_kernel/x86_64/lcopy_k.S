/*! \file lcopy_k.S	
 * This file is part of the LEAC.
 *
 * Implementation of the 
 *	
 *	y <-- x, (long int <-- long int)
 *
 *  function for long int
 *
 *
 * (c)  Hermes Robles Berumen <hermes@uaz.edu.mx>
 *
 * For the full copyright and license information, please view the LICENSE
 * file that was distributed with this source code.
 */

#define ASSEMBLER
#include "common_arc.h"

#define M	ARG1	/* rdi */
#define X	ARG2	/* rsi */
#define INCX	ARG3	/* rdx */
#define Y	ARG4	/* rcx */
#ifndef WINDOWS_ABI
#define INCY	ARG5	/* r8  */
#else
#define INCY	%r10
#endif
#define SIZE    8
	
#include "l1param.h"

#ifdef OPTERON
#define LOAD(OFFSET, ADDR, REG)		xorps	REG, REG; addpd	OFFSET(ADDR), REG
#else
#define LOAD(OFFSET, ADDR, REG)		movdqa	OFFSET(ADDR), REG
#endif


        .file	"lcopy_k.S"
	.text
	.p2align 4,,15
.globl lcopy_k
	.type	lcopy_k, @function
lcopy_k:
.LFB2:


#ifdef WINDOWS_ABI
	movq	40(%rsp), INCY
#endif

	SAVEREGISTERS

	leaq	(, INCX, SIZE), INCX
	leaq	(, INCY, SIZE), INCY

	cmpq	$SIZE, INCX
	jne	.L40
	cmpq	$SIZE, INCY
	jne	.L40

#ifdef ALIGNED_ACCESS
	testq	$SIZE, Y
#else
	testq	$SIZE, X
#endif
	je	.L10

	movq	(X), %xmm0
	movq	%xmm0, (Y)
	addq	$1 * SIZE, X
	addq	$1 * SIZE, Y
	decq	M
	jle	.L19
	ALIGN_4

.L10:
	subq	$-16 * SIZE, X
	subq	$-16 * SIZE, Y

#ifdef ALIGNED_ACCESS
	testq	$SIZE, X
#else
	testq	$SIZE, Y
#endif
	jne	.L20

	movq	M,  %rax
	sarq	$4, %rax
	jle	.L13

	movdqa	-16 * SIZE(X), %xmm0
	movdqa	-14 * SIZE(X), %xmm1
	movdqa	-12 * SIZE(X), %xmm2
	movdqa	-10 * SIZE(X), %xmm3
	movdqa	 -8 * SIZE(X), %xmm4
	movdqa	 -6 * SIZE(X), %xmm5
	movdqa	 -4 * SIZE(X), %xmm6
	movdqa	 -2 * SIZE(X), %xmm7

	decq	%rax
	jle .L12
	ALIGN_3

.L11:
#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	movdqa	%xmm0, -16 * SIZE(Y)
	LOAD( 0 * SIZE, X, %xmm0)
	movdqa	%xmm1, -14 * SIZE(Y)
	LOAD( 2 * SIZE, X, %xmm1)

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movdqa	%xmm2, -12 * SIZE(Y)
	LOAD( 4 * SIZE, X, %xmm2)
	movdqa	%xmm3, -10 * SIZE(Y)
	LOAD( 6 * SIZE, X, %xmm3)

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  64) - PREOFFSET(Y)
#endif

	movdqa	%xmm4, -8 * SIZE(Y)
	LOAD( 8 * SIZE, X, %xmm4)
	movdqa	%xmm5, -6 * SIZE(Y)
	LOAD(10 * SIZE, X, %xmm5)

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(X)
#endif

	movdqa	%xmm6, -4 * SIZE(Y)
	LOAD(12 * SIZE, X, %xmm6)
	movdqa	%xmm7, -2 * SIZE(Y)
	LOAD(14 * SIZE, X, %xmm7)

	subq	$-16 * SIZE, Y
	subq	$-16 * SIZE, X
	decq	%rax
	jg	.L11
	ALIGN_3

.L12:
	movdqa	%xmm0, -16 * SIZE(Y)
	movdqa	%xmm1, -14 * SIZE(Y)
	movdqa	%xmm2, -12 * SIZE(Y)
	movdqa	%xmm3, -10 * SIZE(Y)
	movdqa	%xmm4,  -8 * SIZE(Y)
	movdqa	%xmm5,  -6 * SIZE(Y)
	movdqa	%xmm6,  -4 * SIZE(Y)
	movdqa	%xmm7,  -2 * SIZE(Y)

	subq	$-16 * SIZE, Y
	subq	$-16 * SIZE, X
	ALIGN_3

.L13:
	testq	$8, M
	jle	.L14
	ALIGN_3

	movdqa	-16 * SIZE(X), %xmm0
	movdqa	-14 * SIZE(X), %xmm1
	movdqa	-12 * SIZE(X), %xmm2
	movdqa	-10 * SIZE(X), %xmm3

	movdqa	%xmm0, -16 * SIZE(Y)
	movdqa	%xmm1, -14 * SIZE(Y)
	movdqa	%xmm2, -12 * SIZE(Y)
	movdqa	%xmm3, -10 * SIZE(Y)

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L14:
	testq	$4, M
	jle	.L15
	ALIGN_3

	movdqa	-16 * SIZE(X), %xmm0
	movdqa	-14 * SIZE(X), %xmm1

	movdqa	%xmm0, -16 * SIZE(Y)
	movdqa	%xmm1, -14 * SIZE(Y)

	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L15:
	testq	$2, M
	jle	.L16
	ALIGN_3

	movdqa	-16 * SIZE(X), %xmm0
	movdqa	%xmm0, -16 * SIZE(Y)

	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

.L16:
	testq	$1, M
	jle	.L19
	ALIGN_3

	movq	-16 * SIZE(X), %xmm0
	movq	%xmm0, 	-16 * SIZE(Y)
	ALIGN_3

.L19:
	xorq	%rax,%rax

	RESTOREREGISTERS

	ret
	ALIGN_3

.L20:
#ifdef ALIGNED_ACCESS

	//movhps	-16 * SIZE(X), %xmm0
	movdqa	-17 * SIZE(X), %xmm0

	movq	M,  %rax
	sarq	$4, %rax
	jle	.L23

	movdqa	-15 * SIZE(X), %xmm1
	movdqa	-13 * SIZE(X), %xmm2
	movdqa	-11 * SIZE(X), %xmm3
	movdqa	 -9 * SIZE(X), %xmm4
	movdqa	 -7 * SIZE(X), %xmm5
	movdqa	 -5 * SIZE(X), %xmm6
	movdqa	 -3 * SIZE(X), %xmm7

	decq	%rax
	jle .L22
	ALIGN_4

.L21:
#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif
	psrldq     $8, %xmm0
	punpcklqdq %xmm1, %xmm0
	//SHUFPD_1 %xmm1, %xmm0
	movdqa	%xmm0, -16 * SIZE(Y)
	LOAD(-1 * SIZE, X, %xmm0)

	psrldq     $8, %xmm1
	punpcklqdq %xmm2, %xmm1
	//SHUFPD_1 %xmm2, %xmm1
	movdqa	%xmm1, -14 * SIZE(Y)
	LOAD( 1 * SIZE, X, %xmm1)

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	psrldq     $8, %xmm2
	punpcklqdq %xmm3, %xmm2
	//SHUFPD_1 %xmm3, %xmm2
	movdqa	%xmm2, -12 * SIZE(Y)
	LOAD( 3 * SIZE, X, %xmm2)

	psrldq     $8, %xmm3
	punpcklqdq %xmm4, %xmm3
	//SHUFPD_1 %xmm4, %xmm3
	movdqa	%xmm3, -10 * SIZE(Y)
	LOAD( 5 * SIZE, X, %xmm3)

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  64) - PREOFFSET(Y)
#endif

	psrldq     $8, %xmm4
	punpcklqdq %xmm5, %xmm4
	//SHUFPD_1 %xmm5, %xmm4
	movdqa	%xmm4,  -8 * SIZE(Y)
	LOAD( 7 * SIZE, X, %xmm4)

	psrldq     $8, %xmm5
	punpcklqdq %xmm6, %xmm5
	//SHUFPD_1 %xmm6, %xmm5
	movdqa	%xmm5, -6 * SIZE(Y)
	LOAD( 9 * SIZE, X, %xmm5)

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(X)
#endif

	psrldq     $8, %xmm6
	punpcklqdq %xmm7, %xmm6
	//SHUFPD_1 %xmm7, %xmm6
	movdqa	%xmm6, -4 * SIZE(Y)
	LOAD(11 * SIZE, X, %xmm6)

	psrldq     $8, %xmm7
	punpcklqdq %xmm0, %xmm7
	//SHUFPD_1 %xmm0, %xmm7
	movdqa	%xmm7, -2 * SIZE(Y)
	LOAD(13 * SIZE, X, %xmm7)

	subq	$-16 * SIZE, X
	subq	$-16 * SIZE, Y
	decq	%rax
	jg	.L21
	ALIGN_3

.L22:
	psrldq     $8, %xmm0
	punpcklqdq %xmm1, %xmm0
	//SHUFPD_1 %xmm1, %xmm0
	movdqa	%xmm0, -16 * SIZE(Y)
	LOAD(-1 * SIZE, X, %xmm0)

	psrldq     $8, %xmm1
	punpcklqdq %xmm2, %xmm1
	//SHUFPD_1 %xmm2, %xmm1
	movdqa	%xmm1, -14 * SIZE(Y)

	psrldq     $8, %xmm2
	punpcklqdq %xmm3, %xmm2
	//SHUFPD_1 %xmm3, %xmm2
	movdqa	%xmm2, -12 * SIZE(Y)
	psrldq     $8, %xmm3
	punpcklqdq %xmm4, %xmm3
	//SHUFPD_1 %xmm4, %xmm3
	movdqa	%xmm3, -10 * SIZE(Y)

	psrldq     $8, %xmm4
	punpcklqdq %xmm5, %xmm4
	//SHUFPD_1 %xmm5, %xmm4
	movdqa	%xmm4,  -8 * SIZE(Y)

	psrldq     $8, %xmm5
	punpcklqdq %xmm6, %xmm5
	//SHUFPD_1 %xmm6, %xmm5
	movdqa	%xmm5,  -6 * SIZE(Y)

	psrldq     $8, %xmm6
	punpcklqdq %xmm7, %xmm6
	//SHUFPD_1 %xmm7, %xmm6
	movdqa	%xmm6,  -4 * SIZE(Y)
	psrldq     $8, %xmm7
	punpcklqdq %xmm0, %xmm7
	//SHUFPD_1 %xmm0, %xmm7
	movdqa	%xmm7,  -2 * SIZE(Y)

	subq	$-16 * SIZE, X
	subq	$-16 * SIZE, Y
	ALIGN_3

.L23:
	testq	$8, M
	jle	.L24
	ALIGN_3

	movdqa	-15 * SIZE(X), %xmm1
	movdqa	-13 * SIZE(X), %xmm2
	movdqa	-11 * SIZE(X), %xmm3
	movdqa	 -9 * SIZE(X), %xmm8

	psrldq     $8, %xmm0
	punpcklqdq %xmm1, %xmm0
	//SHUFPD_1 %xmm1, %xmm0
	movdqa	%xmm0, -16 * SIZE(Y)

	psrldq     $8, %xmm1
	punpcklqdq %xmm2, %xmm1
	//SHUFPD_1 %xmm2, %xmm1
	movdqa	%xmm1, -14 * SIZE(Y)

	psrldq     $8, %xmm2
	punpcklqdq %xmm3, %xmm2
	//SHUFPD_1 %xmm3, %xmm2
	movdqa	%xmm2, -12 * SIZE(Y)

	psrldq     $8, %xmm3
	punpcklqdq %xmm8, %xmm3
	//SHUFPD_1 %xmm8, %xmm3
	movdqa	%xmm3, -10 * SIZE(Y)

	movdqa	%xmm8, %xmm0

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L24:
	testq	$4, M
	jle	.L25
	ALIGN_3

	movdqa	-15 * SIZE(X), %xmm1
	movdqa	-13 * SIZE(X), %xmm2

	//SHUFPD_1 %xmm1, %xmm0
	//SHUFPD_1 %xmm2, %xmm1
	psrldq     $8, %xmm0
	punpcklqdq %xmm1, %xmm0
	movdqa	%xmm0, -16 * SIZE(Y)

	psrldq     $8, %xmm1
	punpcklqdq %xmm2, %xmm1
	movdqa	%xmm1, -14 * SIZE(Y)
	
	movdqa	%xmm2, %xmm0

	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L25:
	testq	$2, M
	jle	.L26
	ALIGN_3

	movdqa	-15 * SIZE(X), %xmm1
	psrldq     $8, %xmm0
	punpcklqdq %xmm1, %xmm0
	//SHUFPD_1 %xmm1, %xmm0

	movdqa	%xmm0, -16 * SIZE(Y)

	movdqa	%xmm1, %xmm0
	
	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

.L26:
	testq	$1, M
	jle	.L29
	ALIGN_3

	psrldq   $8, %xmm0
	//movq	-16 * SIZE(X), %xmm0
	movq	%xmm0, 	-16 * SIZE(Y)
	ALIGN_3

.L29:
	xorq	%rax,%rax

	RESTOREREGISTERS

	ret
	ALIGN_3

#else

	movq	M,  %rax
	sarq	$4, %rax
	jle	.L23

	movdqa	-16 * SIZE(X), %xmm0
	movdqa	-14 * SIZE(X), %xmm1
	movdqa	-12 * SIZE(X), %xmm2
	movdqa	-10 * SIZE(X), %xmm3
	movdqa	 -8 * SIZE(X), %xmm4
	movdqa	 -6 * SIZE(X), %xmm5
	movdqa	 -4 * SIZE(X), %xmm6
	movdqa	 -2 * SIZE(X), %xmm7

	decq	%rax
	jle .L22
	ALIGN_3

.L21:
#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif
	
	movq	%xmm0, 	-16 * SIZE(Y)
	psrldq   $8, %xmm0
	movq	%xmm0, 	-15 * SIZE(Y)
	//movlps	%xmm0, -16 * SIZE(Y)
	//movhps	%xmm0, -15 * SIZE(Y)
	LOAD( 0 * SIZE, X, %xmm0)
	movq	%xmm1, 	-14 * SIZE(Y)
	psrldq   $8, %xmm1
	movq	%xmm1, 	-13 * SIZE(Y)
	//movlps	%xmm1, -14 * SIZE(Y)
	//movhps	%xmm1, -13 * SIZE(Y)
	LOAD( 2 * SIZE, X, %xmm1)

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movq	%xmm2, 	-12 * SIZE(Y)
	psrldq   $8, %xmm2
	movq	%xmm2, 	-11 * SIZE(Y)
	//movlps	%xmm2, -12 * SIZE(Y)
	//movhps	%xmm2, -11 * SIZE(Y)
	LOAD( 4 * SIZE, X, %xmm2)
	movq	%xmm3, 	-10 * SIZE(Y)
	psrldq   $8, %xmm3
	movq	%xmm3, 	-9 * SIZE(Y)
	//movlps	%xmm3, -10 * SIZE(Y)
	//movhps	%xmm3,  -9 * SIZE(Y)
	LOAD( 6 * SIZE, X, %xmm3)

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  64) - PREOFFSET(Y)
#endif

	movq	%xmm4, 	-8 * SIZE(Y)
	psrldq   $8, %xmm4
	movq	%xmm4, 	-7 * SIZE(Y)
	//movlps	%xmm4,  -8 * SIZE(Y)
	//movhps	%xmm4,  -7 * SIZE(Y)
	LOAD( 8 * SIZE, X, %xmm4)
	movq	%xmm5, 	-6 * SIZE(Y)
	psrldq   $8, %xmm5
	movq	%xmm5, 	-5 * SIZE(Y)
	//movlps	%xmm5,  -6 * SIZE(Y)
	//movhps	%xmm5,  -5 * SIZE(Y)
	LOAD(10 * SIZE, X, %xmm5)

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(X)
#endif

	movq	%xmm6, 	-4 * SIZE(Y)
	psrldq   $8, %xmm6
	movq	%xmm6, 	-3 * SIZE(Y)
	//movlps	%xmm6,  -4 * SIZE(Y)
	//movhps	%xmm6,  -3 * SIZE(Y)
	LOAD(12 * SIZE, X, %xmm6)
	movq	%xmm7, 	-2 * SIZE(Y)
	psrldq   $8, %xmm7
	movq	%xmm7, 	-1 * SIZE(Y)
	//movlps	%xmm7,  -2 * SIZE(Y)
	//movhps	%xmm7,  -1 * SIZE(Y)
	LOAD(14 * SIZE, X, %xmm7)

	subq	$-16 * SIZE, Y
	subq	$-16 * SIZE, X
	decq	%rax
	jg	.L21
	ALIGN_3

.L22:
	movq	%xmm0, 	-16 * SIZE(Y)
	psrldq   $8, %xmm0
	movq	%xmm0, 	-15 * SIZE(Y)
	//movlps	%xmm0, -16 * SIZE(Y)
	//movhps	%xmm0, -15 * SIZE(Y)
	movq	%xmm1, 	-14 * SIZE(Y)
	psrldq   $8, %xmm1
	movq	%xmm1, 	-13 * SIZE(Y)
	//movlps	%xmm1, -14 * SIZE(Y)
	//movhps	%xmm1, -13 * SIZE(Y)
	movq	%xmm2, 	-12 * SIZE(Y)
	psrldq   $8, %xmm2
	movq	%xmm2, 	-11 * SIZE(Y)
	//movlps	%xmm2, -12 * SIZE(Y)
	//movhps	%xmm2, -11 * SIZE(Y)
	movq	%xmm3, 	-10 * SIZE(Y)
	psrldq   $8, %xmm3
	movq	%xmm3, 	-9 * SIZE(Y)
	//movlps	%xmm3, -10 * SIZE(Y)
	//movhps	%xmm3,  -9 * SIZE(Y)
	movq	%xmm4, 	-8 * SIZE(Y)
	psrldq   $8, %xmm4
	movq	%xmm4, 	-7 * SIZE(Y)
	//movlps	%xmm4,  -8 * SIZE(Y)
	//movhps	%xmm4,  -7 * SIZE(Y)
	movq	%xmm5, 	-6 * SIZE(Y)
	psrldq   $8, %xmm5
	movq	%xmm5, 	-5 * SIZE(Y)
	//movlps	%xmm5,  -6 * SIZE(Y)
	//movhps	%xmm5,  -5 * SIZE(Y)
	movq	%xmm6, 	-4 * SIZE(Y)
	psrldq   $8, %xmm6
	movq	%xmm6, 	-3 * SIZE(Y)
	//movlps	%xmm6,  -4 * SIZE(Y)
	//movhps	%xmm6,  -3 * SIZE(Y)
	movq	%xmm7, 	-2 * SIZE(Y)
	psrldq   $8, %xmm7
	movq	%xmm7, 	-1 * SIZE(Y)
	//movlps	%xmm7,  -2 * SIZE(Y)
	//movhps	%xmm7,  -1 * SIZE(Y)

	subq	$-16 * SIZE, Y
	subq	$-16 * SIZE, X
	ALIGN_3

.L23:
	testq	$8, M
	jle	.L24
	ALIGN_3

	movdqa	-16 * SIZE(X), %xmm0
	movq	%xmm0, 	-16 * SIZE(Y)
	psrldq   $8, %xmm0
	movq	%xmm0, 	-15 * SIZE(Y)
	//movlps	%xmm0, -16 * SIZE(Y)
	//movhps	%xmm0, -15 * SIZE(Y)
	movdqa	-14 * SIZE(X), %xmm1
	movq	%xmm1, 	-14 * SIZE(Y)
	psrldq   $8, %xmm1
	movq	%xmm1, 	-13 * SIZE(Y)
	//movlps	%xmm1, -14 * SIZE(Y)
	//movhps	%xmm1, -13 * SIZE(Y)
	movdqa	-12 * SIZE(X), %xmm2
	movq	%xmm2, 	-12 * SIZE(Y)
	psrldq   $8, %xmm2
	movq	%xmm2, 	-11 * SIZE(Y)
	//movlps	%xmm2, -12 * SIZE(Y)
	//movhps	%xmm2, -11 * SIZE(Y)
	movdqa	-10 * SIZE(X), %xmm3
	movq	%xmm3, 	-10 * SIZE(Y)
	psrldq   $8, %xmm3
	movq	%xmm3, 	-9 * SIZE(Y)
	//movlps	%xmm3, -10 * SIZE(Y)
	//movhps	%xmm3,  -9 * SIZE(Y)

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L24:
	testq	$4, M
	jle	.L25
	ALIGN_3

	movdqa	-16 * SIZE(X), %xmm0
	movq	%xmm0, 	-16 * SIZE(Y)
	psrldq   $8, %xmm0
	movq	%xmm0, 	-15 * SIZE(Y)
	//movlps	%xmm0, -16 * SIZE(Y)
	//movhps	%xmm0, -15 * SIZE(Y)
	movdqa	-14 * SIZE(X), %xmm1
	movq	%xmm1, 	-14 * SIZE(Y)
	psrldq   $8, %xmm1
	movq	%xmm1, 	-13 * SIZE(Y)
	//movlps	%xmm1, -14 * SIZE(Y)
	//movhps	%xmm1, -13 * SIZE(Y)

	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L25:
	testq	$2, M
	jle	.L26
	ALIGN_3

	movdqa	-16 * SIZE(X), %xmm0
	movq	%xmm0, 	-16 * SIZE(Y)
	psrldq   $8, %xmm0
	movq	%xmm0, 	-15 * SIZE(Y)
	//movlps	%xmm0, -16 * SIZE(Y)
	//movhps	%xmm0, -15 * SIZE(Y)
	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

.L26:
	testq	$1, M
	jle	.L29
	ALIGN_3

	movq	-16 * SIZE(X), %xmm0
	movq	%xmm0, 	-16 * SIZE(Y)
	ALIGN_3

.L29:
	xorq	%rax,%rax

	RESTOREREGISTERS

	ret
	ALIGN_3

#endif

.L40:
	movq	M,  %rax
	sarq	$3, %rax
	jle	.L45
	ALIGN_3

.L41:
	movq	(X), %xmm0
	addq	INCX, X
	movq	(X), %xmm1
	addq	INCX, X
	movq	(X), %xmm2
	addq	INCX, X
	movq	(X), %xmm3
	addq	INCX, X
	movq	(X), %xmm4
	addq	INCX, X
	movq	(X), %xmm5
	addq	INCX, X
	movq	(X), %xmm6
	addq	INCX, X
	movq	(X), %xmm7
	addq	INCX, X

	movlps	%xmm0, (Y)
	addq	INCY, Y
	movq	%xmm1, (Y)
	addq	INCY, Y
	movlps	%xmm2, (Y)
	addq	INCY, Y
	movq	%xmm3, (Y)
	addq	INCY, Y
	movlps	%xmm4, (Y)
	addq	INCY, Y
	movq	%xmm5, (Y)
	addq	INCY, Y
	movlps	%xmm6, (Y)
	addq	INCY, Y
	movq	%xmm7, (Y)
	addq	INCY, Y

	decq	%rax
	jg	.L41
	ALIGN_3

.L45:
	movq	M,  %rax
	andq	$7, %rax
	jle	.L47
	ALIGN_3

.L46:
	movq	(X), %xmm0
	addq	INCX, X
	movq	%xmm0, (Y)
	addq	INCY, Y
	decq	%rax
	jg	.L46
	ALIGN_3

.L47:
	xorq	%rax, %rax

	RESTOREREGISTERS

	ret

.LFE2:
	.size	lcopy_k, .-lcopy_k
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


