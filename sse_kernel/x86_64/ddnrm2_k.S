/*! \file ddnrm2_k.S	
 * This file is part of the LEAC.
 *
 * Implementation of the 
 *	
 *  nom2 <-- ||x||2, (double <-- ||double-double||2)
 *
 *  function for double
 *
 *
 * (c)  Hermes Robles Berumen <hermes@uaz.edu.mx>
 *
 * For the full copyright and license information, please view the LICENSE
 * file that was distributed with this source code.
 */
	
#define ASSEMBLER
#include "common_arc.h"

#define ARG1	%rdi
#define ARG2	%rsi
#define ARG3	%rdx
#define ARG4	%rcx
#define ARG5	%r8
#define ARG6	%r9

#define N	ARG1	/* rdi */
#define X	ARG2	/* rsi */
#define INCX	ARG3    /* rdx */
#define Y	ARG4	/* rcx */
#define INCY	ARG5	/* r8  */
#define SIZE    8

#include "l1param.h"

        .file	"ddnrm2_k.S"
	.text
	.p2align 4,,15
        .globl  ddnrm2_k
	.type	ddnrm2_k, @function
ddnrm2_k:
.LFB2:
	leaq	(, INCX, SIZE), INCX	
	leaq	(, INCY, SIZE), INCY	

	xorps	%xmm0, %xmm0
	xorps	%xmm1, %xmm1
	xorps	%xmm2, %xmm2
	xorps	%xmm3, %xmm3

	cmpq	$0, N
	jle	.L999

	cmpq	$SIZE, INCX
	jne	.L50
	cmpq	$SIZE, INCY
	jne	.L50

	subq	$-16 * SIZE, X
	subq	$-16 * SIZE, Y
  
        testq	$SIZE, Y
	je	.L10

	movsd	-16 * SIZE(X), %xmm0
	subsd	-16 * SIZE(Y), %xmm0
        mulsd   %xmm0, %xmm0
	addq	$1 * SIZE, X
	addq	$1 * SIZE, Y
	decq	N
	ALIGN_2
.L10:
	testq	$SIZE, X
	jne	.L20

        movq	N,  %rax
	sarq	$4, %rax
	jle	.L14

	movapd	-16 * SIZE(X), %xmm4
	movapd	-14 * SIZE(X), %xmm5
	movapd	-12 * SIZE(X), %xmm6
	movapd	-10 * SIZE(X), %xmm7

	movapd	 -8 * SIZE(X), %xmm8
	movapd	 -6 * SIZE(X), %xmm9
	movapd	 -4 * SIZE(X), %xmm10
	movapd	 -2 * SIZE(X), %xmm11

	decq	%rax
	jle	.L12

	ALIGN_3

.L11:
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	subpd	-16 * SIZE(Y), %xmm4
        mulpd    %xmm4, %xmm4
	addpd	 %xmm4, %xmm0
	movapd	  0 * SIZE(X), %xmm4

	subpd	-14 * SIZE(Y), %xmm5
        mulpd    %xmm5, %xmm5
	addpd	 %xmm5,  %xmm1
	movapd	  2 * SIZE(X), %xmm5

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	subpd	-12 * SIZE(Y), %xmm6
        mulpd    %xmm6, %xmm6
        addpd	 %xmm6,  %xmm2
	movapd	  4 * SIZE(X), %xmm6

	subpd	-10 * SIZE(Y), %xmm7
        mulpd    %xmm7, %xmm7
	addpd	 %xmm7,  %xmm3
	movapd	  6 * SIZE(X), %xmm7

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(X)
#endif

	subpd	 -8 * SIZE(Y), %xmm8
        mulpd    %xmm8, %xmm8
	addpd	 %xmm8,  %xmm0
	movapd	  8 * SIZE(X), %xmm8

	subpd	 -6 * SIZE(Y), %xmm9
        mulpd    %xmm9, %xmm9
	addpd	 %xmm9,  %xmm1
	movapd	 10 * SIZE(X), %xmm9

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(Y)
#endif

	subpd	 -4 * SIZE(Y), %xmm10
        mulpd    %xmm10, %xmm10
	addpd	 %xmm10, %xmm2
	movapd	 12 * SIZE(X), %xmm10

	subpd	 -2 * SIZE(Y), %xmm11
        mulpd    %xmm11, %xmm11
	addpd	 %xmm11, %xmm3
	movapd	 14 * SIZE(X), %xmm11

	subq	$-16 * SIZE, X
	subq	$-16 * SIZE, Y

	decq	%rax
	jg	.L11
	ALIGN_3

.L12:
	subpd	-16 * SIZE(Y), %xmm4
        mulpd    %xmm4, %xmm4
	addpd	 %xmm4,  %xmm0
	subpd	-14 * SIZE(Y), %xmm5
        mulpd    %xmm5, %xmm5
	addpd	 %xmm5,  %xmm1
	subpd	-12 * SIZE(Y), %xmm6
        mulpd    %xmm6, %xmm6
	addpd	 %xmm6,  %xmm2
	subpd	-10 * SIZE(Y), %xmm7
        mulpd    %xmm7, %xmm7
	addpd	 %xmm7,  %xmm3

	subpd	 -8 * SIZE(Y), %xmm8
        mulpd    %xmm8, %xmm8
	addpd	 %xmm8,  %xmm0
	subpd	 -6 * SIZE(Y), %xmm9
        mulpd    %xmm9, %xmm9
	addpd	 %xmm9,  %xmm1
	subpd	 -4 * SIZE(Y), %xmm10
        mulpd    %xmm10, %xmm10
	addpd	 %xmm10, %xmm2
	subpd	 -2 * SIZE(Y), %xmm11
        mulpd    %xmm11, %xmm11
	addpd	 %xmm11, %xmm3

	subq	$-16 * SIZE, X
	subq	$-16 * SIZE, Y
	ALIGN_3

.L14:
	testq	$15, N
	jle	.L999

	testq	$8, N
	jle	.L15

	movapd	-16 * SIZE(X), %xmm4
	movapd	-14 * SIZE(X), %xmm5
	movapd	-12 * SIZE(X), %xmm6
	movapd	-10 * SIZE(X), %xmm7

	subpd	-16 * SIZE(Y), %xmm4
        mulpd    %xmm4, %xmm4
	addpd	 %xmm4,  %xmm0
	subpd	-14 * SIZE(Y), %xmm5
        mulpd    %xmm5, %xmm5
	addpd	 %xmm5,  %xmm1
	subpd	-12 * SIZE(Y), %xmm6
        mulpd    %xmm6, %xmm6
	addpd	 %xmm6,  %xmm2
	subpd	-10 * SIZE(Y), %xmm7
        mulpd    %xmm7, %xmm7
	addpd	 %xmm7,  %xmm3

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L15:
	testq	$4, N
	jle	.L16

	movapd	-16 * SIZE(X), %xmm4
	movapd	-14 * SIZE(X), %xmm5

        subpd	-16 * SIZE(Y), %xmm4
        mulpd    %xmm4, %xmm4
      	addpd	 %xmm4, %xmm0
	subpd	-14 * SIZE(Y), %xmm5
        mulpd    %xmm5, %xmm5
        addpd	 %xmm5, %xmm1

	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L16:
	testq	$2, N
	jle	.L17

	movapd	-16 * SIZE(X), %xmm4
	subpd	-16 * SIZE(Y), %xmm4
        mulpd	 %xmm4, %xmm4
	addpd	 %xmm4, %xmm0

	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

.L17:
	testq	$1, N
	jle	.L999

	movsd	-16 * SIZE(X), %xmm4
        subsd	-16 * SIZE(Y), %xmm4
	mulsd	%xmm4, %xmm4
	addsd	%xmm4,  %xmm0
	jmp	.L999
	ALIGN_3

.L20:

#ifdef ALIGNED_ACCESS

	movhpd	-16 * SIZE(X), %xmm4
	addq	$SIZE, X

	movq	N,  %rax
	sarq	$4, %rax
	jle	.L24

	movapd	-16 * SIZE(X), %xmm5
	movapd	-14 * SIZE(X), %xmm6
	movapd	-12 * SIZE(X), %xmm7
	movapd	-10 * SIZE(X), %xmm8

	movapd	 -8 * SIZE(X), %xmm9
	movapd	 -6 * SIZE(X), %xmm10
	movapd	 -4 * SIZE(X), %xmm11

	decq	%rax
	jle	.L22

	ALIGN_3

.L21:
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	SHUFPD_1 %xmm5, %xmm4
	subpd	-16 * SIZE(Y), %xmm4
        mulpd    %xmm4, %xmm4
	addpd	 %xmm4,  %xmm0
	movaps	 -2 * SIZE(X), %xmm4

	SHUFPD_1 %xmm6, %xmm5
	subpd	-14 * SIZE(Y), %xmm5
        mulpd    %xmm5, %xmm5
	addpd	 %xmm5, %xmm1
	movaps	  0 * SIZE(X), %xmm5

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	SHUFPD_1 %xmm7, %xmm6
	subpd	-12 * SIZE(Y), %xmm6
        mulpd    %xmm6, %xmm6
	addpd	 %xmm6,  %xmm2
	movaps	  2 * SIZE(X), %xmm6

	SHUFPD_1 %xmm8, %xmm7
	subpd	-10 * SIZE(Y), %xmm7
        mulpd    %xmm7, %xmm7
	addpd	 %xmm7,  %xmm3
	movaps	  4 * SIZE(X), %xmm7

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(X)
#endif

	SHUFPD_1 %xmm9, %xmm8
	subpd	 -8 * SIZE(Y), %xmm8
        mulpd    %xmm8, %xmm8
	addpd	 %xmm8, %xmm0
	movaps	  6 * SIZE(X), %xmm8

	SHUFPD_1 %xmm10, %xmm9
	subpd	 -6 * SIZE(Y), %xmm9
        mulpd    %xmm9, %xmm9
	addpd	 %xmm9,  %xmm1
	movaps	  8 * SIZE(X), %xmm9

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(Y)
#endif

	SHUFPD_1 %xmm11, %xmm10
	subpd	 -4 * SIZE(Y), %xmm10
        mulpd    %xmm10, %xmm10
	addpd	 %xmm10, %xmm2
	movaps	 10 * SIZE(X), %xmm10

	SHUFPD_1 %xmm4, %xmm11
	subpd	 -2 * SIZE(Y), %xmm11
        mulpd    %xmm11, %xmm11
	addpd	 %xmm11, %xmm3
	movaps	 12 * SIZE(X), %xmm11

	subq	$-16 * SIZE, X
	subq	$-16 * SIZE, Y

	decq	%rax
	jg	.L21
	ALIGN_3

.L22:
	SHUFPD_1 %xmm5, %xmm4
	subpd	-16 * SIZE(Y), %xmm4
        mulpd    %xmm4, %xmm4
	addpd	 %xmm4,  %xmm0
	movaps	 -2 * SIZE(X), %xmm4

	SHUFPD_1 %xmm6, %xmm5
	subpd	-14 * SIZE(Y), %xmm5
        mulpd    %xmm5, %xmm5
	addpd	 %xmm5,  %xmm1

	SHUFPD_1 %xmm7, %xmm6
	subpd	-12 * SIZE(Y), %xmm6
        mulpd    %xmm6, %xmm6
	addpd	 %xmm6,  %xmm2

	SHUFPD_1 %xmm8, %xmm7
	subpd	-10 * SIZE(Y), %xmm7
        mulpd    %xmm7, %xmm7
	addpd	 %xmm7,  %xmm3

	SHUFPD_1 %xmm9, %xmm8
	subpd	 -8 * SIZE(Y), %xmm8
        mulpd    %xmm8, %xmm8
	addpd	 %xmm8,  %xmm0

	SHUFPD_1 %xmm10, %xmm9
	subpd	 -6 * SIZE(Y), %xmm9
        mulpd    %xmm9, %xmm9
	addpd	 %xmm9,  %xmm1

	SHUFPD_1 %xmm11, %xmm10
	subpd	 -4 * SIZE(Y), %xmm10
        mulpd    %xmm10, %xmm10
	addpd	 %xmm10, %xmm2

	SHUFPD_1 %xmm4, %xmm11
	subpd	 -2 * SIZE(Y), %xmm11
        mulpd    %xmm11, %xmm11
	addpd	 %xmm11, %xmm3

	subq	$-16 * SIZE, X
	subq	$-16 * SIZE, Y
	ALIGN_3

.L24:
	testq	$15, N
	jle	.L999

	testq	$8, N
	jle	.L25

	movaps	-16 * SIZE(X), %xmm5
	movaps	-14 * SIZE(X), %xmm6
	movaps	-12 * SIZE(X), %xmm7

	SHUFPD_1 %xmm5, %xmm4
	subpd	-16 * SIZE(Y), %xmm4
        mulpd    %xmm4, %xmm4
	addpd	 %xmm4,  %xmm0
	movaps	-10 * SIZE(X), %xmm4

	SHUFPD_1 %xmm6, %xmm5
	subpd	-14 * SIZE(Y), %xmm5
        mulpd    %xmm5, %xmm5
	addpd	 %xmm5,  %xmm1

	SHUFPD_1 %xmm7, %xmm6
	subpd	-12 * SIZE(Y), %xmm6
        mulpd    %xmm6, %xmm6
	addpd	 %xmm6,  %xmm2

	SHUFPD_1 %xmm4, %xmm7
	subpd	-10 * SIZE(Y), %xmm7
        mulpd    %xmm7, %xmm7
	addpd	 %xmm7,  %xmm3

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L25:
	testq	$4, N
	jle	.L26

	movaps	-16 * SIZE(X), %xmm5
	movaps	-14 * SIZE(X), %xmm6

	SHUFPD_1 %xmm5, %xmm4
	subpd	-16 * SIZE(Y), %xmm4
        mulpd    %xmm4, %xmm4
	addpd	 %xmm4, %xmm0
	SHUFPD_1 %xmm6, %xmm5
	subpd	-14 * SIZE(Y), %xmm5
        mulpd    %xmm5, %xmm5
	addpd	 %xmm5, %xmm1
	movapd	 %xmm6, %xmm4

	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L26:
	testq	$2, N
	jle	.L27

	movaps	-16 * SIZE(X), %xmm5

	SHUFPD_1 %xmm5, %xmm4
	subpd	-16 * SIZE(Y), %xmm4
        mulpd    %xmm4, %xmm4
	addpd	 %xmm4,  %xmm0
	movapd	 %xmm5, %xmm4

	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

.L27:
	testq	$1, N
	jle	.L999

	SHUFPD_1 %xmm4, %xmm4
	subsd	-16 * SIZE(Y), %xmm4
        mulsd   %xmm4, %xmm4
	addsd	%xmm4, %xmm0
	jmp	.L999
	ALIGN_3

#else

	movq	N,  %rax
	sarq	$4, %rax
	jle	.L24

	movlps	-16 * SIZE(X), %xmm4
	movhps	-15 * SIZE(X), %xmm4
	movlps	-14 * SIZE(X), %xmm5
	movhps	-13 * SIZE(X), %xmm5
	movlps	-12 * SIZE(X), %xmm6
	movhps	-11 * SIZE(X), %xmm6
	movlps	-10 * SIZE(X), %xmm7
	movhps	 -9 * SIZE(X), %xmm7

	movlps	 -8 * SIZE(X), %xmm8
	movhps	 -7 * SIZE(X), %xmm8
	movlps	 -6 * SIZE(X), %xmm9
	movhps	 -5 * SIZE(X), %xmm9
	movlps	 -4 * SIZE(X), %xmm10
	movhps	 -3 * SIZE(X), %xmm10
	movlps	 -2 * SIZE(X), %xmm11
	movhps	 -1 * SIZE(X), %xmm11

	decq	%rax
	jle	.L22

	ALIGN_3

.L21:
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	subpd	-16 * SIZE(Y), %xmm4
        mulpd    %xmm4, %xmm4
	addpd	 %xmm4, %xmm0
	movlps	  0 * SIZE(X), %xmm4
	movhps	  1 * SIZE(X), %xmm4

	subpd	-14 * SIZE(Y), %xmm5
        mulpd    %xmm5, %xmm5
	addpd	 %xmm5,  %xmm1
	movlps	  2 * SIZE(X), %xmm5
	movhps	  3 * SIZE(X), %xmm5

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	subpd	-12 * SIZE(Y), %xmm6
        mulpd    %xmm6, %xmm6
	addpd	 %xmm6,  %xmm2
	movlps	  4 * SIZE(X), %xmm6
	movhps	  5 * SIZE(X), %xmm6

	subpd	-10 * SIZE(Y), %xmm7
        mulpd    %xmm7, %xmm7
	addpd	 %xmm7,  %xmm3
	movlps	  6 * SIZE(X), %xmm7
	movhps	  7 * SIZE(X), %xmm7

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(X)
#endif

	subpd	 -8 * SIZE(Y), %xmm8
        mulpd    %xmm8, %xmm8
	addpd	 %xmm8,  %xmm0
	movlps	  8 * SIZE(X), %xmm8
	movhps	  9 * SIZE(X), %xmm8

	subpd	 -6 * SIZE(Y), %xmm9
        mulpd    %xmm9, %xmm9
	addpd	 %xmm9,  %xmm1
	movlps	 10 * SIZE(X), %xmm9
	movhps	 11 * SIZE(X), %xmm9

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(Y)
#endif

	subpd	 -4 * SIZE(Y), %xmm10
        mulpd    %xmm10, %xmm10
	addpd	 %xmm10, %xmm2
	movlps	 12 * SIZE(X), %xmm10
	movhps	 13 * SIZE(X), %xmm10

	subpd	 -2 * SIZE(Y), %xmm11
        mulpd    %xmm11, %xmm11
	addpd	 %xmm11, %xmm3
	movlps	 14 * SIZE(X), %xmm11
	movhps	 15 * SIZE(X), %xmm11

	subq	$-16 * SIZE, X
	subq	$-16 * SIZE, Y

	decq	%rax
	jg	.L21
	ALIGN_3

.L22:
	subpd	-16 * SIZE(Y), %xmm4
        mulpd    %xmm4, %xmm4
	addpd	 %xmm4,  %xmm0
	subpd	-14 * SIZE(Y), %xmm5
        mulpd    %xmm5, %xmm5
	addpd	 %xmm5,  %xmm1
	subpd	-12 * SIZE(Y), %xmm6
        mulpd    %xmm6, %xmm6
	addpd	 %xmm6,  %xmm2
	subpd	-10 * SIZE(Y), %xmm7
        mulpd    %xmm7, %xmm7
	addpd	 %xmm7,  %xmm3

	subpd	 -8 * SIZE(Y), %xmm8
        mulpd    %xmm8, %xmm8
	addpd	 %xmm8,  %xmm0
	subpd	 -6 * SIZE(Y), %xmm9
        mulpd    %xmm9, %xmm9
	addpd	 %xmm9,  %xmm1
	subpd	 -4 * SIZE(Y), %xmm10
        mulpd    %xmm10, %xmm10
	addpd	 %xmm10, %xmm2
	subpd	 -2 * SIZE(Y), %xmm11
        mulpd    %xmm11, %xmm11
	addpd	 %xmm11, %xmm3

	subq	$-16 * SIZE, X
	subq	$-16 * SIZE, Y
	ALIGN_3

.L24:
	testq	$15, N
	jle	.L999

	testq	$8, N
	jle	.L25

	movlps	-16 * SIZE(X), %xmm4
	movhps	-15 * SIZE(X), %xmm4
	movlps	-14 * SIZE(X), %xmm5
	movhps	-13 * SIZE(X), %xmm5
	movlps	-12 * SIZE(X), %xmm6
	movhps	-11 * SIZE(X), %xmm6
	movlps	-10 * SIZE(X), %xmm7
	movhps	 -9 * SIZE(X), %xmm7

	subpd	-16 * SIZE(Y), %xmm4
        mulpd    %xmm4, %xmm4
	addpd	 %xmm4,  %xmm0
	subpd	-14 * SIZE(Y), %xmm5
        mulpd    %xmm5, %xmm5
	addpd	 %xmm5,  %xmm1
	subpd	-12 * SIZE(Y), %xmm6
        mulpd    %xmm6, %xmm6
	addpd	 %xmm6, %xmm2
	subpd	-10 * SIZE(Y), %xmm7
        mulpd    %xmm7, %xmm7
	addpd	 %xmm7, %xmm3

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L25:
	testq	$4, N
	jle	.L26

	movlps	-16 * SIZE(X), %xmm4
	movhps	-15 * SIZE(X), %xmm4
	movlps	-14 * SIZE(X), %xmm5
	movhps	-13 * SIZE(X), %xmm5

	subpd	-16 * SIZE(Y), %xmm4
        mulpd    %xmm4, %xmm4
	addpd	 %xmm4,  %xmm0
	subpd	-14 * SIZE(Y), %xmm5
        mulpd    %xmm5, %xmm5
	addpd	 %xmm5, %xmm1

	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L26:
	testq	$2, N
	jle	.L27

	movlps	-16 * SIZE(X), %xmm4
	movhps	-15 * SIZE(X), %xmm4

	subpd	-16 * SIZE(Y), %xmm4
        mulpd    %xmm4, %xmm4
	addpd	 %xmm4, %xmm0

	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

.L27:
	testq	$1, N
	jle	.L999

	movsd	-16 * SIZE(X), %xmm4
	subsd	-16 * SIZE(Y), %xmm4
        mulsd	%xmm4, %xmm4
	addsd	%xmm4, %xmm0
	jmp	.L999
	ALIGN_3
#endif


.L50:
	movq	N,  %rax
	sarq	$2, %rax
	jle	.L55
	ALIGN_3

.L53:
	movsd	0 * SIZE(X), %xmm4
	addq	INCX, X
	subsd	0 * SIZE(Y), %xmm4
        mulsd	%xmm4, %xmm4
	addq	INCY, Y
	movsd	0 * SIZE(X), %xmm5
	addq	INCX, X
	subsd	0 * SIZE(Y), %xmm5
        mulsd	%xmm5, %xmm5
	addq	INCY, Y
	movsd	0 * SIZE(X), %xmm6
	addq	INCX, X
	subsd	0 * SIZE(Y), %xmm6
        mulsd	%xmm6, %xmm6
	addq	INCY, Y
	movsd	0 * SIZE(X), %xmm7
	addq	INCX, X
	subsd	0 * SIZE(Y), %xmm7
        mulsd	%xmm7, %xmm7
	addq	INCY, Y

	addsd	%xmm4, %xmm0
	addsd	%xmm5, %xmm1
	addsd	%xmm6, %xmm2
	addsd	%xmm7, %xmm3

	decq	%rax
	jg	.L53
	ALIGN_3

.L55:
	movq	N, %rax
	andq	$3,   %rax
	jle	.L999
	ALIGN_3

.L56:
	movsd	0 * SIZE(X), %xmm4
	addq	INCX, X
	subsd	0 * SIZE(Y), %xmm4
        mulsd	%xmm4, %xmm4
	addq	INCY, Y
	addsd	%xmm4, %xmm0
	decq	%rax
	jg	.L56
	ALIGN_3

.L999:
	addpd	%xmm1, %xmm0
	addpd	%xmm3, %xmm2
	addpd	%xmm2, %xmm0

#ifndef HAVE_SSE3
	pshufd	$0xe, %xmm0, %xmm1
	addsd	%xmm1, %xmm0
#else
	haddpd	%xmm0, %xmm0
#endif

	ret
.LFE2:
	.size	ddnrm2_k, .-ddnrm2_k
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
