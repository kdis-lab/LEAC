/*! \file icopy_k.S	
 * This file is part of the LEAC.
 *
 * Implementation of the 
 *	
 *	y <-- x, (int <-- int)
 *
 *  function for int
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
#define SIZE    4
	
#include "l1param.h"

#ifdef OPTERON
#define LOAD(OFFSET, ADDR, REG)		xorps	REG, REG; paddq	OFFSET(ADDR), REG
#else
#define LOAD(OFFSET, ADDR, REG)		movdqa	OFFSET(ADDR), REG
#endif

        .file	"icopy_k.S"
	.text
	.p2align 4,,15
.globl icopy_k
	.type	icopy_k, @function
icopy_k:
.LFB2:

	
#ifdef WINDOWS_ABI
	movq	40(%rsp), INCY
#endif

	SAVEREGISTERS

	leaq	(, INCX, SIZE), INCX
	leaq	(, INCY, SIZE), INCY

	cmpq	$SIZE, INCX
	jne	.L50
	cmpq	$SIZE, INCY
	jne	.L50

	cmpq	$3, M
	jle	.L55

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y

	testq	$SIZE, Y
	je	.L05

	movd	-32 * SIZE(X), %xmm0
	movd	%xmm0, -32 * SIZE(Y)
	addq	$1 * SIZE, X
	addq	$1 * SIZE, Y
	decq	M
	ALIGN_4

.L05:
	testq	$2 * SIZE, Y
	je	.L10

	movq	-32 * SIZE(X), %xmm0
	movlps	%xmm0, -32 * SIZE(Y)
	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	subq	$2, M
	jle	.L19
	ALIGN_4

.L10:
	testq	$3 * SIZE, X
	jne	.L20

	movq	M,  %rax
	sarq	$5, %rax
	jle	.L13

	movdqa	-32 * SIZE(X), %xmm0
	movdqa	-28 * SIZE(X), %xmm1
	movdqa	-24 * SIZE(X), %xmm2
	movdqa	-20 * SIZE(X), %xmm3
	movdqa	-16 * SIZE(X), %xmm4
	movdqa	-12 * SIZE(X), %xmm5
	movdqa	 -8 * SIZE(X), %xmm6
	movdqa	 -4 * SIZE(X), %xmm7

	decq	%rax
	jle .L12
	ALIGN_3

.L11:
#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	movdqa	%xmm0, -32 * SIZE(Y)
	LOAD( 0 * SIZE, X, %xmm0)
	movdqa	%xmm1, -28 * SIZE(Y)
	LOAD( 4 * SIZE, X, %xmm1)

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movdqa	%xmm2, -24 * SIZE(Y)
	LOAD( 8 * SIZE, X, %xmm2)
	movdqa	%xmm3, -20 * SIZE(Y)
	LOAD(12 * SIZE, X, %xmm3)

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  64) - PREOFFSET(Y)
#endif

	movdqa	%xmm4,-16 * SIZE(Y)
	LOAD(16 * SIZE, X, %xmm4)
	movdqa	%xmm5,-12 * SIZE(Y)
	LOAD(20 * SIZE, X, %xmm5)

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(X)
#endif

	movdqa	%xmm6, -8 * SIZE(Y)
	LOAD(24 * SIZE, X, %xmm6)
	movdqa	%xmm7, -4 * SIZE(Y)
	LOAD(28 * SIZE, X, %xmm7)

	subq	$-32 * SIZE, Y
	subq	$-32 * SIZE, X
	decq	%rax
	jg	.L11
	ALIGN_3

.L12:
	movdqa	%xmm0, -32 * SIZE(Y)
	movdqa	%xmm1, -28 * SIZE(Y)
	movdqa	%xmm2, -24 * SIZE(Y)
	movdqa	%xmm3, -20 * SIZE(Y)
	movdqa	%xmm4, -16 * SIZE(Y)
	movdqa	%xmm5, -12 * SIZE(Y)
	movdqa	%xmm6,  -8 * SIZE(Y)
	movdqa	%xmm7,  -4 * SIZE(Y)

	subq	$-32 * SIZE, Y
	subq	$-32 * SIZE, X
	ALIGN_3

.L13:
	testq	$16, M
	jle	.L14

	movdqa	-32 * SIZE(X), %xmm0
	movdqa	-28 * SIZE(X), %xmm1
	movdqa	-24 * SIZE(X), %xmm2
	movdqa	-20 * SIZE(X), %xmm3

	movdqa	%xmm0, -32 * SIZE(Y)
	movdqa	%xmm1, -28 * SIZE(Y)
	movdqa	%xmm2, -24 * SIZE(Y)
	movdqa	%xmm3, -20 * SIZE(Y)

	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
	ALIGN_3

.L14:
	testq	$8, M
	jle	.L15

	movdqa	-32 * SIZE(X), %xmm0
	movdqa	-28 * SIZE(X), %xmm1

	movdqa	%xmm0, -32 * SIZE(Y)
	movdqa	%xmm1, -28 * SIZE(Y)

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L15:
	testq	$4, M
	jle	.L16

	movdqa	-32 * SIZE(X), %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)

	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L16:
	testq	$2, M
	jle	.L17

	movq	-32 * SIZE(X), %xmm0
	movlps	%xmm0, -32 * SIZE(Y)

	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

.L17:
	testq	$1, M
	jle	.L19

	movd	-32 * SIZE(X), %xmm0
	movd	%xmm0, 	-32 * SIZE(Y)
	ALIGN_3

.L19:
	xorq	%rax,%rax

	RESTOREREGISTERS

	ret
	ALIGN_3


.L20:
	testq	$SIZE, X
	jne	.L30

	movdqa	-34 * SIZE(X), %xmm0

	movq	M,  %rax
	sarq	$5, %rax
	jle	.L23

	movdqa	-30 * SIZE(X), %xmm1
	movdqa	-26 * SIZE(X), %xmm2
	movdqa	-22 * SIZE(X), %xmm3
	movdqa	-18 * SIZE(X), %xmm4
	movdqa	-14 * SIZE(X), %xmm5
	movdqa	-10 * SIZE(X), %xmm6
	movdqa	 -6 * SIZE(X), %xmm7

	decq	%rax
	jle .L22
	ALIGN_4

.L21:

#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

        psrldq     $8, %xmm0
	punpcklqdq %xmm1, %xmm0	
	movdqa	%xmm0, -32 * SIZE(Y)
	movdqa	-2 * SIZE(X), %xmm0

	psrldq     $8, %xmm1
	punpcklqdq %xmm2, %xmm1
	movdqa	%xmm1, -28 * SIZE(Y)
	movdqa	  2 * SIZE(X), %xmm1

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	psrldq     $8, %xmm2
	punpcklqdq %xmm3, %xmm2		
	movdqa	%xmm2, -24 * SIZE(Y)
	movdqa	  6 * SIZE(X), %xmm2

	psrldq     $8, %xmm3
	punpcklqdq %xmm4, %xmm3
	movdqa	%xmm3, -20 * SIZE(Y)
	movdqa	 10 * SIZE(X), %xmm3

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  64) - PREOFFSET(Y)
#endif

	psrldq     $8, %xmm4
	punpcklqdq %xmm5, %xmm4
	movdqa	%xmm4, -16 * SIZE(Y)
	movdqa	 14 * SIZE(X), %xmm4

	psrldq     $8, %xmm5
	punpcklqdq %xmm6, %xmm5
	movdqa	%xmm5, -12 * SIZE(Y)
	movdqa	 18 * SIZE(X), %xmm5

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(X)
#endif

	psrldq     $8, %xmm6
	punpcklqdq %xmm7, %xmm6
	movdqa	%xmm6,  -8 * SIZE(Y)
	movdqa	 22 * SIZE(X), %xmm6

	psrldq     $8, %xmm7
	punpcklqdq %xmm0, %xmm7
	movdqa	%xmm7,  -4 * SIZE(Y)
	movdqa	 26 * SIZE(X), %xmm7

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	decq	%rax
	jg	.L21
	ALIGN_3

.L22:
	psrldq     $8, %xmm0
	punpcklqdq %xmm1, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)
	movdqa	 -2 * SIZE(X), %xmm0

	psrldq     $8, %xmm1
	punpcklqdq %xmm2, %xmm1
	movdqa	%xmm1, -28 * SIZE(Y)

	psrldq     $8, %xmm2
	punpcklqdq %xmm3, %xmm2
	movdqa	%xmm2, -24 * SIZE(Y)

	psrldq     $8, %xmm3
	punpcklqdq %xmm4, %xmm3
	movdqa	%xmm3, -20 * SIZE(Y)

	psrldq     $8, %xmm4
	punpcklqdq %xmm5, %xmm4
	movdqa	%xmm4, -16 * SIZE(Y)

	psrldq     $8, %xmm5
	punpcklqdq %xmm6, %xmm5
	movdqa	%xmm5, -12 * SIZE(Y)

	psrldq     $8, %xmm6
	punpcklqdq %xmm7, %xmm6
	movdqa	%xmm6,  -8 * SIZE(Y)

	psrldq     $8, %xmm7
	punpcklqdq %xmm0, %xmm7
	movdqa	%xmm7,  -4 * SIZE(Y)

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	ALIGN_3

.L23:
	testq	$16, M
	jle	.L24
	ALIGN_3

	movdqa	-30 * SIZE(X), %xmm1
	movdqa	-26 * SIZE(X), %xmm2
	movdqa	-22 * SIZE(X), %xmm3
	movdqa	-18 * SIZE(X), %xmm4

	psrldq     $8, %xmm0
	punpcklqdq %xmm1, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)
	psrldq     $8, %xmm1
	punpcklqdq %xmm2, %xmm1
	movdqa	%xmm1, -28 * SIZE(Y)
	psrldq     $8, %xmm2
	punpcklqdq %xmm3, %xmm2
	movdqa	%xmm2, -24 * SIZE(Y)
	psrldq     $8, %xmm3
	punpcklqdq %xmm4, %xmm3
	movdqa	%xmm3, -20 * SIZE(Y)

	movdqa	%xmm4, %xmm0

	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
	ALIGN_3

.L24:
	testq	$8, M
	jle	.L25
	ALIGN_3

	movdqa	-30 * SIZE(X), %xmm1
	movdqa	-26 * SIZE(X), %xmm2

	psrldq     $8, %xmm0
	punpcklqdq %xmm1, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)
	
	psrldq     $8, %xmm1
	punpcklqdq %xmm2, %xmm1
	movdqa	%xmm1, -28 * SIZE(Y)

	movdqa	%xmm2, %xmm0

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L25:
	testq	$4, M
	jle	.L26
	ALIGN_3

	movdqa	-30 * SIZE(X), %xmm1
	psrldq     $8, %xmm0
	punpcklqdq %xmm1, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)

	movdqa	%xmm1, %xmm0
	
	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L26:
	testq	$2, M
	jle	.L27
	ALIGN_3

	psrldq   $8, %xmm0
	movq	 %xmm0, -32 * SIZE(Y)
	
	movdqa	-30 * SIZE(X), %xmm0
	pslldq     $8, %xmm0
	
	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

.L27:
	testq	$1, M
	jle	.L29
	ALIGN_3
	
        psrldq     $8, %xmm0
	movd	%xmm0, 	-32 * SIZE(Y)
	addq	$SIZE, Y
	ALIGN_3

.L29:
	xorq	%rax,%rax

	RESTOREREGISTERS

	ret
	ALIGN_3

.L30:
	testq	$2 * SIZE, X
	jne	.L40

	movdqa	-33 * SIZE(X), %xmm0

	movq	M,  %rax
	sarq	$5, %rax
	jle	.L33

	movdqa	-29 * SIZE(X), %xmm1
	movdqa	-25 * SIZE(X), %xmm2
	movdqa	-21 * SIZE(X), %xmm3
	movdqa	-17 * SIZE(X), %xmm4
	movdqa	-13 * SIZE(X), %xmm5
	movdqa	 -9 * SIZE(X), %xmm6
	movdqa	 -5 * SIZE(X), %xmm7

	decq	%rax
	jle .L32
	ALIGN_4

.L31:
#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	psrldq  $4, %xmm0
	pshufd  $0xe4,%xmm1, %xmm8
	pslldq  $12, %xmm8
        por     %xmm8, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)
	movdqa	 -1 * SIZE(X), %xmm0

	psrldq  $4, %xmm1
	pshufd  $0xe4,%xmm2, %xmm9
	pslldq  $12, %xmm9
        por     %xmm9, %xmm1
	movdqa	%xmm1, -28 * SIZE(Y)
	movdqa	  3 * SIZE(X), %xmm1

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	psrldq  $4, %xmm2
	pshufd  $0xe4,%xmm3, %xmm10
	pslldq  $12, %xmm10
        por     %xmm10, %xmm2
	movdqa	%xmm2, -24 * SIZE(Y)
	movdqa	  7 * SIZE(X), %xmm2

	psrldq  $4, %xmm3
	pshufd  $0xe4,%xmm4, %xmm11
	pslldq  $12, %xmm11
        por     %xmm11, %xmm3
	movdqa	%xmm3, -20 * SIZE(Y)
	movdqa	 11 * SIZE(X), %xmm3

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  64) - PREOFFSET(Y)
#endif

	psrldq  $4, %xmm4
	pshufd  $0xe4,%xmm5, %xmm12
	pslldq  $12, %xmm12
        por     %xmm12, %xmm4
	movdqa	%xmm4, -16 * SIZE(Y)
	movdqa	 15 * SIZE(X), %xmm4

	psrldq  $4, %xmm5
	pshufd  $0xe4,%xmm6, %xmm13
	pslldq  $12, %xmm13
        por     %xmm13, %xmm5
	movdqa	%xmm5, -12 * SIZE(Y)
	movdqa	 19 * SIZE(X), %xmm5

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(X)
#endif
	psrldq  $4, %xmm6
	pshufd  $0xe4,%xmm7, %xmm14
	pslldq  $12, %xmm14
        por     %xmm14, %xmm6
	movdqa	%xmm6,  -8 * SIZE(Y)
	movdqa	 23 * SIZE(X), %xmm6

	psrldq  $4, %xmm7
	pshufd  $0xe4,%xmm0, %xmm15
	pslldq  $12, %xmm15
        por     %xmm15, %xmm7
	movdqa	%xmm7,  -4 * SIZE(Y)
	movdqa	 27 * SIZE(X), %xmm7

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	decq	%rax
	jg	.L31
	ALIGN_3

.L32:
	pshufd  $0xe4,%xmm1, %xmm8
	psrldq  $4, %xmm0
	pslldq  $12, %xmm8
        por     %xmm8, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)
	movdqa	 -1 * SIZE(X), %xmm0

	pshufd  $0xe4,%xmm2, %xmm9
	psrldq  $4, %xmm1
	pslldq  $12, %xmm9
        por     %xmm9, %xmm1
	movdqa	%xmm1, -28 * SIZE(Y)

	psrldq  $4, %xmm2
	pshufd  $0xe4,%xmm3, %xmm10
	pslldq  $12, %xmm10
        por     %xmm10, %xmm2
	movdqa	%xmm2, -24 * SIZE(Y)

	psrldq  $4, %xmm3
	pshufd  $0xe4,%xmm4, %xmm11
	pslldq  $12, %xmm11
        por     %xmm11, %xmm3
	movdqa	%xmm3, -20 * SIZE(Y)

	psrldq  $4, %xmm4
	pshufd  $0xe4,%xmm5, %xmm12
	pslldq  $12, %xmm12
        por     %xmm12, %xmm4
	movdqa	%xmm4, -16 * SIZE(Y)

	psrldq  $4, %xmm5
	pshufd  $0xe4,%xmm6, %xmm13
	pslldq  $12, %xmm13
        por     %xmm13, %xmm5
	movdqa	%xmm5, -12 * SIZE(Y)

	pshufd  $0xe4,%xmm7, %xmm14
	psrldq  $4, %xmm6
	pslldq  $12, %xmm14
        por     %xmm14, %xmm6
	movdqa	%xmm6,  -8 * SIZE(Y)

	pshufd  $0xe4,%xmm0, %xmm15
	psrldq  $4, %xmm7
	pslldq  $12, %xmm15
        por     %xmm15, %xmm7
	movdqa	%xmm7,  -4 * SIZE(Y)

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	ALIGN_3

.L33:
	testq	$16, M
	jle	.L34
	ALIGN_3

	movdqa	-29 * SIZE(X), %xmm1
	movdqa	-25 * SIZE(X), %xmm2
	movdqa	-21 * SIZE(X), %xmm3
	movdqa	-17 * SIZE(X), %xmm4

	pshufd  $0xe4,%xmm1, %xmm5
	psrldq  $4, %xmm0
	pslldq  $12,%xmm5
        por     %xmm5,%xmm0
 	movdqa	%xmm0, -32 * SIZE(Y)

	psrldq  $4, %xmm1
	pshufd  $0xe4,%xmm2, %xmm6
	pslldq  $12, %xmm6
        por     %xmm6, %xmm1
	movdqa	%xmm1, -28 * SIZE(Y)

	psrldq  $4, %xmm2
	pshufd  $0xe4,%xmm3, %xmm7
	pslldq  $12, %xmm7
        por     %xmm7, %xmm2
 	movdqa	%xmm2, -24 * SIZE(Y)

	psrldq  $4, %xmm3
	pshufd  $0xe4,%xmm4, %xmm8
	pslldq  $12, %xmm8
        por     %xmm8, %xmm3
	movdqa	%xmm3, -20 * SIZE(Y)

	pshufd  $0xe4,%xmm4, %xmm0
	

	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
	ALIGN_3

.L34:
	testq	$8, M
	jle	.L35
	ALIGN_3

	movdqa	-29 * SIZE(X), %xmm1
	movdqa	-25 * SIZE(X), %xmm2

	psrldq  $4, %xmm0
	pshufd  $0xe4,%xmm1, %xmm3
	pslldq  $12, %xmm3
        por     %xmm3, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)

	psrldq  $4, %xmm1
	pshufd  $0xe4,%xmm2, %xmm0
	pslldq  $12, %xmm2
        por     %xmm2, %xmm1

	movdqa	%xmm1, -28 * SIZE(Y)

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L35:
	testq	$4, M
	jle	.L36
	ALIGN_3

	movdqa	-29 * SIZE(X), %xmm1

	pshufd  $0xe4,%xmm1, %xmm2
	psrldq  $4, %xmm0
	pslldq  $12, %xmm2
        por     %xmm2, %xmm0

	movdqa	%xmm0, -32 * SIZE(Y)
	pshufd  $0xe4,%xmm1, %xmm0

	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L36:
	testq	$2, M
	jle	.L37
	ALIGN_3

	psrldq  $4, %xmm0
	movq	%xmm0, -32 * SIZE(Y)
	psrldq  $4, %xmm0
	
	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

.L37:
	testq	$1, M
	jle	.L39
	ALIGN_3

	psrldq  $4, %xmm0
	movd	%xmm0, 	-32 * SIZE(Y)
	addq	$SIZE, Y
	ALIGN_3

.L39:
	xorq	%rax,%rax

	RESTOREREGISTERS

	ret
	ALIGN_3

.L40:
	movdqa	-35 * SIZE(X), %xmm0

	movq	M,  %rax
	sarq	$5, %rax
	jle	.L43

	movdqa	-31 * SIZE(X), %xmm1
	movdqa	-27 * SIZE(X), %xmm2
	movdqa	-23 * SIZE(X), %xmm3
	movdqa	-19 * SIZE(X), %xmm4
	movdqa	-15 * SIZE(X), %xmm5
	movdqa	-11 * SIZE(X), %xmm6
	movdqa	 -7 * SIZE(X), %xmm7

	decq	%rax
	jle .L42
	ALIGN_4

.L41:
#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif
	
	pshufd  $0xe4,%xmm1,%xmm8
	psrldq  $12,%xmm0
	pslldq  $4, %xmm8
	por     %xmm8, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)
	movdqa	 -3 * SIZE(X), %xmm0

	pshufd  $0xe4,%xmm2, %xmm9
	psrldq  $12,%xmm1
	pslldq  $4, %xmm9
	por     %xmm9, %xmm1
	movdqa	%xmm1, -28 * SIZE(Y)
	movdqa	  1 * SIZE(X), %xmm1

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	pshufd  $0xe4,%xmm3, %xmm10
	psrldq  $12,  %xmm2
	pslldq  $4, %xmm10
	por     %xmm10, %xmm2
	movdqa	%xmm2, -24 * SIZE(Y)
	movdqa	  5 * SIZE(X), %xmm2

	pshufd  $0xe4,%xmm4, %xmm11
	psrldq  $12,  %xmm3
	pslldq  $4, %xmm11
	por     %xmm11, %xmm3
	movdqa	%xmm3, -20 * SIZE(Y)
	movdqa	  9 * SIZE(X), %xmm3

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  64) - PREOFFSET(Y)
#endif

	pshufd  $0xe4,%xmm5, %xmm12
	psrldq  $12,  %xmm4
	pslldq  $4, %xmm12
	por     %xmm12, %xmm4
	movdqa	%xmm4, -16 * SIZE(Y)
	movdqa	 13 * SIZE(X), %xmm4

	pshufd  $0xe4,%xmm6, %xmm13
	psrldq  $12,  %xmm5
	pslldq  $4, %xmm13
	por     %xmm13, %xmm5
	movdqa	%xmm5, -12 * SIZE(Y)
	movdqa	 17 * SIZE(X), %xmm5

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(X)
#endif

	pshufd  $0xe4,%xmm7, %xmm14
	psrldq  $12,  %xmm6
	pslldq  $4, %xmm14
	por     %xmm14, %xmm6
	movdqa	%xmm6, -8 * SIZE(Y)
	movdqa	 21 * SIZE(X), %xmm6

	pshufd  $0xe4,%xmm0, %xmm15
	psrldq  $12,  %xmm7
	pslldq  $4, %xmm15
	por     %xmm15, %xmm7
	movdqa	%xmm7, -4 * SIZE(Y)
	movdqa	 25 * SIZE(X), %xmm7

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	decq	%rax
	jg	.L41
	ALIGN_3

.L42:
	pshufd  $0xe4,%xmm1, %xmm8
	psrldq  $12, %xmm0
	pslldq  $4, %xmm8
        por     %xmm8, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)
	movdqa	 -3 * SIZE(X), %xmm0

	pshufd  $0xe4,%xmm2, %xmm9
	psrldq  $12, %xmm1
	pslldq  $4, %xmm9
        por     %xmm9, %xmm1
	movdqa	%xmm1, -28 * SIZE(Y)

	pshufd  $0xe4,%xmm3, %xmm10
	psrldq  $12, %xmm2
	pslldq  $4, %xmm10
        por     %xmm10, %xmm2
	movdqa	%xmm2, -24 * SIZE(Y)

	pshufd  $0xe4,%xmm4, %xmm11
	psrldq  $12, %xmm3
	pslldq  $4, %xmm11
        por     %xmm11, %xmm3
	movdqa	%xmm3, -20 * SIZE(Y)

	pshufd  $0xe4,%xmm5, %xmm12
	psrldq  $12, %xmm4
	pslldq  $4, %xmm12
        por     %xmm12, %xmm4
	movdqa	%xmm4, -16 * SIZE(Y)

	pshufd  $0xe4,%xmm6, %xmm13
	psrldq  $12, %xmm5
	pslldq  $4, %xmm13
        por     %xmm13, %xmm5
	movdqa	%xmm5, -12 * SIZE(Y)

	pshufd  $0xe4,%xmm7, %xmm14
	psrldq  $12, %xmm6
	pslldq  $4, %xmm14
        por     %xmm14, %xmm6
	movdqa	%xmm6, -8 * SIZE(Y)

	pshufd  $0xe4,%xmm0, %xmm15
	psrldq  $12, %xmm7
	pslldq  $4, %xmm15
        por     %xmm15, %xmm7
	movdqa	%xmm7, -4 * SIZE(Y)

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	ALIGN_3

.L43:
	testq	$16, M
	jle	.L44
	ALIGN_3

	movdqa	-31 * SIZE(X), %xmm1
	movdqa	-27 * SIZE(X), %xmm2
	movdqa	-23 * SIZE(X), %xmm3
	movdqa	-19 * SIZE(X), %xmm4

	pshufd  $0xe4,%xmm1, %xmm5
	psrldq  $12, %xmm0
	pslldq  $4, %xmm5
        por     %xmm5, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)

	pshufd  $0xe4,%xmm2, %xmm6
	psrldq  $12, %xmm1
	pslldq  $4, %xmm6
        por     %xmm6, %xmm1
	movdqa	%xmm1, -28 * SIZE(Y)

	pshufd  $0xe4,%xmm3, %xmm7
	psrldq  $12, %xmm2
	pslldq  $4, %xmm7
        por     %xmm7, %xmm2
	movdqa	%xmm2, -24 * SIZE(Y)

	pshufd  $0xe4,%xmm4, %xmm8
	psrldq  $12, %xmm3
	pslldq  $4, %xmm8
        por     %xmm8, %xmm3
	movdqa	%xmm3, -20 * SIZE(Y)

	pshufd  $0xe4,%xmm4, %xmm0

	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
	ALIGN_3

.L44:
	testq	$8, M
	jle	.L45
	ALIGN_3

	movdqa	-31 * SIZE(X), %xmm1
	movdqa	-27 * SIZE(X), %xmm2

	pshufd  $0xe4,%xmm1, %xmm3
	psrldq  $12, %xmm0
	pslldq  $4, %xmm3
        por     %xmm3, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)

	pshufd  $0xe4,%xmm2, %xmm4
	psrldq  $12, %xmm1
	pslldq  $4, %xmm4
        por     %xmm4, %xmm1
	movdqa	%xmm1, -28 * SIZE(Y)

	pshufd  $0xe4,%xmm2, %xmm0

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L45:
	testq	$4, M
	jle	.L46
	ALIGN_3

	movdqa	-31 * SIZE(X), %xmm1
	
	pshufd  $0xe4,%xmm1, %xmm2
	psrldq  $12, %xmm0
	pslldq  $4, %xmm2
        por     %xmm2, %xmm0

	movdqa	%xmm0, -32 * SIZE(Y)
	pshufd  $0xe4,%xmm1, %xmm0

	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L46:
	testq	$2, M
	jle	.L47
	ALIGN_3

	movdqa	-31 * SIZE(X), %xmm1

	pshufd  $0xe4,%xmm1, %xmm2
	psrldq  $12, %xmm0
	pslldq  $4, %xmm2
        por     %xmm2, %xmm0
	
	movq	%xmm0, -32 * SIZE(Y)
	pslldq  $4, %xmm0

	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

.L47:
	testq	$1, M
	jle	.L49
	ALIGN_3

	psrldq  $12, %xmm0
	movd	%xmm0, 	-32 * SIZE(Y)
	addq	$SIZE, Y
	ALIGN_3

.L49:
	xorq	%rax,%rax

	RESTOREREGISTERS

	ret
	ALIGN_4

.L50:
	movq	M,  %rax
	sarq	$3, %rax
	jle	.L55
	ALIGN_3

.L51:
	movd	(X), %xmm0
	addq	INCX, X
	movd	(X), %xmm1
	addq	INCX, X
	movd	(X), %xmm2
	addq	INCX, X
	movd	(X), %xmm3
	addq	INCX, X
	movd	(X), %xmm4
	addq	INCX, X
	movd	(X), %xmm5
	addq	INCX, X
	movd	(X), %xmm6
	addq	INCX, X
	movd	(X), %xmm7
	addq	INCX, X

	movd	%xmm0, (Y)
	addq	INCY, Y
	movd	%xmm1, (Y)
	addq	INCY, Y
	movd	%xmm2, (Y)
	addq	INCY, Y
	movd	%xmm3, (Y)
	addq	INCY, Y
	movd	%xmm4, (Y)
	addq	INCY, Y
	movd	%xmm5, (Y)
	addq	INCY, Y
	movd	%xmm6, (Y)
	addq	INCY, Y
	movd	%xmm7, (Y)
	addq	INCY, Y

	decq	%rax
	jg	.L51
	ALIGN_3

.L55:
	movq	M,  %rax
	andq	$7, %rax
	jle	.L57
	ALIGN_3

.L56:
	movd	(X), %xmm0
	addq	INCX, X
	movd	%xmm0, (Y)
	addq	INCY, Y
	decq	%rax
	jg	.L56
	ALIGN_3

.L57:
	xorq	%rax, %rax

	RESTOREREGISTERS

	ret

.LFE2:
	.size	icopy_k, .-icopy_k
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


