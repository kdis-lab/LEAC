/*! \file iswap_kh.S	
 * This file is part of the LEAC.
 *
 * Implementation of the 
 *	
 *	y <--> x,
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

#ifndef WINDOWS_ABI
#define M	ARG1  /*%rdi*/
#define X	ARG5  /*%r8*/
#define INCX	%r10
#define Y	ARG3  /*%rdx*/
#define INCY	ARG6  /*%r9*/
#else
#define M	ARG1   /*%rcx*/
#define X	ARG2   /*%rdx*/
#define INCX	ARG3   /*%r8*/
#define Y	ARG4   /*%r9*/
#define INCY	%rbx  
#endif

#define SIZE    4
	
#include "l1param.h"

.file	"iswap_kh.S"
	.text
	.p2align 4,,15
.globl iswap_kh
	.type	iswap_kh, @function
iswap_kh:
.LFB2:	

/*	PROLOGUE
	PROFCODE
*/

#ifndef WINDOWS_ABI
	movq	16(%rsp), INCX
	movq	8(%rsp),  Y    /*INCY*/
	
#else
	pushq	%rbx

	movq	48(%rsp), X
	movq	56(%rsp), INCX
	movq	64(%rsp), Y
	movq	72(%rsp), INCY
#endif

	SAVEREGISTERS

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

	movd	-32 * SIZE(X), %xmm0
	movd	-32 * SIZE(Y), %xmm1

	movd	%xmm1, -32 * SIZE(X)
	movd	%xmm0, -32 * SIZE(Y)

	addq	$1 * SIZE, X
	addq	$1 * SIZE, Y
	decq	M
	ALIGN_3

.L05:
	testq	$2 * SIZE, Y
	je	.L10

	movq	-32 * SIZE(X), %xmm0
	movq	-32 * SIZE(Y), %xmm1

	movq	%xmm1, -32 * SIZE(X)
	movq	%xmm0, -32 * SIZE(Y)

	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	subq	$2, M
	jle	.L19
	ALIGN_3

.L10:
	cmpq	$3, M
	jle	.L16

	testq	$2 * SIZE, X
	jne	.L30

	testq	$1 * SIZE, X
	jne	.L20

	movq	M,  %rax
	sarq	$5, %rax
	jle	.L13
	ALIGN_3

.L11:
#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movdqa	-32 * SIZE(X), %xmm0
	movdqa	-32 * SIZE(Y), %xmm1

	movdqa	%xmm0, 	-32 * SIZE(Y)
	movdqa	%xmm1, 	-32 * SIZE(X)

	movdqa	-28 * SIZE(X), %xmm0
	movdqa	-28 * SIZE(Y), %xmm1

	movdqa	%xmm0, 	-28 * SIZE(Y)
	movdqa	%xmm1, 	-28 * SIZE(X)

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	movdqa	-24 * SIZE(X), %xmm0
	movdqa	-24 * SIZE(Y), %xmm1

	movdqa	%xmm0, 	-24 * SIZE(Y)
	movdqa	%xmm1, 	-24 * SIZE(X)

	movdqa	-20 * SIZE(X), %xmm0
	movdqa	-20 * SIZE(Y), %xmm1

	movdqa	%xmm0, 	-20 * SIZE(Y)
	movdqa	%xmm1, 	-20 * SIZE(X)

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE + 64) - PREOFFSET(X)
#endif

	movdqa	-16 * SIZE(X), %xmm0
	movdqa	-16 * SIZE(Y), %xmm1

	movdqa	%xmm0, 	-16 * SIZE(Y)
	movdqa	%xmm1, 	-16 * SIZE(X)

	movdqa	-12 * SIZE(X), %xmm0
	movdqa	-12 * SIZE(Y), %xmm1

	movdqa	%xmm0, 	-12 * SIZE(Y)
	movdqa	%xmm1, 	-12 * SIZE(X)

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE + 64) - PREOFFSET(Y)
#endif

	movdqa	 -8 * SIZE(X), %xmm0
	movdqa	 -8 * SIZE(Y), %xmm1

	movdqa	%xmm0, 	 -8 * SIZE(Y)
	movdqa	%xmm1, 	 -8 * SIZE(X)

	movdqa	 -4 * SIZE(X), %xmm0
	movdqa	 -4 * SIZE(Y), %xmm1

	movdqa	%xmm0, 	 -4 * SIZE(Y)
	movdqa	%xmm1, 	 -4 * SIZE(X)

	subq	$-32 * SIZE, Y
	subq	$-32 * SIZE, X

	decq	%rax
	jg	.L11
	ALIGN_3

.L13:
	testq	$16, M
	jle	.L14

	movdqa	-32 * SIZE(X), %xmm0
	movdqa	-32 * SIZE(Y), %xmm1

	movdqa	%xmm0, 	-32 * SIZE(Y)
	movdqa	%xmm1, 	-32 * SIZE(X)

	movdqa	-28 * SIZE(X), %xmm0
	movdqa	-28 * SIZE(Y), %xmm1

	movdqa	%xmm0, 	-28 * SIZE(Y)
	movdqa	%xmm1, 	-28 * SIZE(X)

	movdqa	-24 * SIZE(X), %xmm0
	movdqa	-24 * SIZE(Y), %xmm1

	movdqa	%xmm0, 	-24 * SIZE(Y)
	movdqa	%xmm1, 	-24 * SIZE(X)

	movdqa	-20 * SIZE(X), %xmm0
	movdqa	-20 * SIZE(Y), %xmm1

	movdqa	%xmm0, 	-20 * SIZE(Y)
 	movdqa	%xmm1, 	-20 * SIZE(X)

	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
	ALIGN_3

.L14:
	testq	$8, M
	jle	.L15

	movdqa	-32 * SIZE(X), %xmm0
	movdqa	-32 * SIZE(Y), %xmm1

	movdqa	%xmm0, 	-32 * SIZE(Y)
	movdqa	%xmm1, 	-32 * SIZE(X)

	movdqa	-28 * SIZE(X), %xmm0
	movdqa	-28 * SIZE(Y), %xmm1

	movdqa	%xmm0, 	-28 * SIZE(Y)
	movdqa	%xmm1, 	-28 * SIZE(X)

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L15:
	testq	$4, M
	jle	.L16

	movdqa	-32 * SIZE(X), %xmm0
	movdqa	-32 * SIZE(Y), %xmm1

	movdqa	%xmm0, 	-32 * SIZE(Y)
	movdqa	%xmm1, 	-32 * SIZE(X)

	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L16:
	testq	$2, M
	jle	.L17

	movq	-32 * SIZE(X), %xmm0
	movq	-32 * SIZE(Y), %xmm1

	movq	%xmm1, 	-32 * SIZE(X)
	addq	$2 * SIZE, X
	movq	%xmm0, 	-32 * SIZE(Y)
	addq	$2 * SIZE, Y
	ALIGN_3

.L17:
	testq	$1, M
	jle	.L19

	movd	-32 * SIZE(X), %xmm0
	movd	-32 * SIZE(Y), %xmm1

	movd	%xmm1, 	-32 * SIZE(X)
	movd	%xmm0, 	-32 * SIZE(Y)
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
	movdqa	 -33 * SIZE(X), %xmm0
	movdqa	 -32 * SIZE(Y), %xmm1
	
	movd	%xmm1, -32 * SIZE(X)
	pshufd	$0x39, %xmm1, %xmm3
	movq	%xmm3, -31 * SIZE(X)

	subq	$3, M

	movq	M,  %rax
	sarq	$5, %rax
	jle	.L23
	ALIGN_4

.L21:
#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movdqa	-29 * SIZE(X), %xmm2
	movdqa	-28 * SIZE(Y), %xmm3

	psrldq  $4, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $12, %xmm8
        por     %xmm8, %xmm0
	//movss	%xmm2, %xmm0
	//shufps $0x39, %xmm0, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)

	psrldq  $12, %xmm1
	pshufd  $0xe4,%xmm3, %xmm9
	pslldq  $4, %xmm9
        por     %xmm9, %xmm1
	//movss	%xmm3, %xmm1
	//shufps	$0x93, %xmm3, %xmm1
	movdqa	%xmm1, -29 * SIZE(X)

	movdqa	-25 * SIZE(X), %xmm0
	movdqa	-24 * SIZE(Y), %xmm1

	psrldq  $4, %xmm2
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $12, %xmm8
        por     %xmm8, %xmm2
	//movss	%xmm0, %xmm2
	//shufps	$0x39, %xmm2, %xmm2
	movdqa	%xmm2, -28 * SIZE(Y)

	psrldq  $12, %xmm3
	pshufd  $0xe4,%xmm1, %xmm9
	pslldq  $4, %xmm9
        por     %xmm9, %xmm3
	//movss	%xmm1, %xmm3
	//shufps	$0x93, %xmm1, %xmm3
	movdqa	%xmm3, -25 * SIZE(X)

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	movdqa	-21 * SIZE(X), %xmm2
	movdqa	-20 * SIZE(Y), %xmm3

	psrldq  $4, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $12, %xmm8
        por     %xmm8, %xmm0
	//movss	%xmm2, %xmm0
	//shufps	$0x39, %xmm0, %xmm0
	movdqa	%xmm0, -24 * SIZE(Y)

	psrldq  $12, %xmm1
	pshufd  $0xe4,%xmm3, %xmm9
	pslldq  $4, %xmm9
        por     %xmm9, %xmm1
	//movss	%xmm3, %xmm1
	//shufps	$0x93, %xmm3, %xmm1
	movdqa	%xmm1, -21 * SIZE(X)

	movdqa	-17 * SIZE(X), %xmm0
	movdqa	-16 * SIZE(Y), %xmm1

	psrldq  $4, %xmm2
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $12, %xmm8
        por     %xmm8, %xmm2
	//movss	%xmm0, %xmm2
	//shufps	$0x39, %xmm2, %xmm2
	movdqa	%xmm2, -20 * SIZE(Y)

	psrldq  $12, %xmm3
	pshufd  $0xe4,%xmm1, %xmm9
	pslldq  $4, %xmm9
        por     %xmm9, %xmm3
	//movss	%xmm1, %xmm3
	//shufps	$0x93, %xmm1, %xmm3
	movdqa	%xmm3, -17 * SIZE(X)

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE + 64) - PREOFFSET(X)
#endif

	movdqa	-13 * SIZE(X), %xmm2
	movdqa	-12 * SIZE(Y), %xmm3

	psrldq  $4, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $12, %xmm8
        por     %xmm8, %xmm0
	//movss	%xmm2, %xmm0
	//shufps	$0x39, %xmm0, %xmm0
	movdqa	%xmm0, -16 * SIZE(Y)

	psrldq  $12, %xmm1
	pshufd  $0xe4,%xmm3, %xmm9
	pslldq  $4, %xmm9
        por     %xmm9, %xmm1
	//movss	%xmm3, %xmm1
	//shufps	$0x93, %xmm3, %xmm1
	movdqa	%xmm1, -13 * SIZE(X)

	movdqa	 -9 * SIZE(X), %xmm0
	movdqa	 -8 * SIZE(Y), %xmm1

	psrldq  $4, %xmm2
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $12, %xmm8
        por     %xmm8, %xmm2
	//movss	%xmm0, %xmm2
	//shufps	$0x39, %xmm2, %xmm2
	movdqa	%xmm2, -12 * SIZE(Y)

	psrldq  $12, %xmm3
	pshufd  $0xe4,%xmm1, %xmm9
	pslldq  $4, %xmm9
        por     %xmm9, %xmm3
	//movss	%xmm1, %xmm3
	//shufps	$0x93, %xmm1, %xmm3
	movdqa	%xmm3,  -9 * SIZE(X)

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE + 64) - PREOFFSET(Y)
#endif

	movdqa	 -5 * SIZE(X), %xmm2
	movdqa	 -4 * SIZE(Y), %xmm3

	psrldq  $4, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $12, %xmm8
        por     %xmm8, %xmm0
	//movss	%xmm2, %xmm0
	//shufps	$0x39, %xmm0, %xmm0
	movdqa	%xmm0,  -8 * SIZE(Y)

	psrldq  $12, %xmm1
	pshufd  $0xe4,%xmm3, %xmm9
	pslldq  $4, %xmm9
        por     %xmm9, %xmm1
	//movss	%xmm3, %xmm1
	//shufps	$0x93, %xmm3, %xmm1
	movdqa	%xmm1,  -5 * SIZE(X)

	movdqa	 -1 * SIZE(X), %xmm0
	movdqa	  0 * SIZE(Y), %xmm1

	psrldq  $4, %xmm2
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $12, %xmm8
        por     %xmm8, %xmm2
	//movss	%xmm0, %xmm2
	//shufps	$0x39, %xmm2, %xmm2
	movdqa	%xmm2,  -4 * SIZE(Y)

	psrldq  $12, %xmm3
	pshufd  $0xe4,%xmm1, %xmm9
	pslldq  $4, %xmm9
        por     %xmm9, %xmm3
	//movss	%xmm1, %xmm3
	//shufps	$0x93, %xmm1, %xmm3
	movdqa	%xmm3,  -1 * SIZE(X)

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y

	decq	%rax
	jg	.L21
	ALIGN_3

.L23:
	testq	$16, M
	jle	.L24

	movdqa	-29 * SIZE(X), %xmm2
	movdqa	-28 * SIZE(Y), %xmm3

	psrldq  $4,   %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $12, %xmm8
        por     %xmm8, %xmm0
	//movss	%xmm2, %xmm0
	//shufps	$0x39, %xmm0, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)

	psrldq  $12, %xmm1
	pshufd  $0xe4,%xmm3, %xmm9
	pslldq  $4,%xmm9
        por     %xmm9, %xmm1
	//movss	%xmm3, %xmm1
	//shufps	$0x93, %xmm3, %xmm1
	movdqa	%xmm1, -29 * SIZE(X)

	movdqa	-25 * SIZE(X), %xmm0
	movdqa	-24 * SIZE(Y), %xmm1

	psrldq  $4, %xmm2
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $12, %xmm8
        por     %xmm8, %xmm2
	//movss	%xmm0, %xmm2
	//shufps	$0x39, %xmm2, %xmm2
	movdqa	%xmm2, -28 * SIZE(Y)

	psrldq  $12, %xmm3
	pshufd  $0xe4,%xmm1, %xmm9
	pslldq  $4, %xmm9
        por     %xmm9, %xmm3
	//movss	%xmm1, %xmm3
	//shufps	$0x93, %xmm1, %xmm3
	movdqa	%xmm3, -25 * SIZE(X)

	movdqa	-21 * SIZE(X), %xmm2
	movdqa	-20 * SIZE(Y), %xmm3

	psrldq  $4, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $12, %xmm8
        por     %xmm8, %xmm0
	//movss	%xmm2, %xmm0
	//shufps	$0x39, %xmm0, %xmm0
	movdqa	%xmm0, -24 * SIZE(Y)

	psrldq  $12, %xmm1
	pshufd  $0xe4,%xmm3, %xmm9
	pslldq  $4, %xmm9
        por     %xmm9, %xmm1
	//movss	%xmm3, %xmm1
	//shufps	$0x93, %xmm3, %xmm1
	movdqa	%xmm1, -21 * SIZE(X)

	movdqa	-17 * SIZE(X), %xmm0
	movdqa	-16 * SIZE(Y), %xmm1

	psrldq  $4, %xmm2
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $12, %xmm8
        por     %xmm8, %xmm2
	//movss	%xmm0, %xmm2
	//shufps	$0x39, %xmm2, %xmm2
	movdqa	%xmm2, -20 * SIZE(Y)

	psrldq  $12, %xmm3
	pshufd  $0xe4,%xmm1, %xmm9
	pslldq  $4, %xmm9
        por     %xmm9, %xmm3
	//movss	%xmm1, %xmm3
	//shufps	$0x93, %xmm1, %xmm3
	movdqa	%xmm3, -17 * SIZE(X)

	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
	ALIGN_3

.L24:
	testq	$8, M
	jle	.L25

	movdqa	-29 * SIZE(X), %xmm2
	movdqa	-28 * SIZE(Y), %xmm3

	psrldq  $4,   %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $12, %xmm8
        por     %xmm8, %xmm0
	//movss	%xmm2, %xmm0
	//shufps	$0x39, %xmm0, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)

	psrldq  $12, %xmm1
	pshufd  $0xe4,%xmm3, %xmm9
	pslldq  $4,%xmm9
        por     %xmm9, %xmm1
	//movss	%xmm3, %xmm1
	//shufps	$0x93, %xmm3, %xmm1
	movdqa	%xmm1, -29 * SIZE(X)

	movdqa	-25 * SIZE(X), %xmm0
	movdqa	-24 * SIZE(Y), %xmm1

	psrldq  $4, %xmm2
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $12, %xmm8
        por     %xmm8, %xmm2
	//movss	%xmm0, %xmm2
	//shufps	$0x39, %xmm2, %xmm2
	movdqa	%xmm2, -28 * SIZE(Y)

	psrldq  $12, %xmm3
	pshufd  $0xe4,%xmm1, %xmm9
	pslldq  $4, %xmm9
        por     %xmm9, %xmm3
	//movss	%xmm1, %xmm3
	//shufps	$0x93, %xmm1, %xmm3
	movdqa	%xmm3, -25 * SIZE(X)

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L25:
	testq	$4, M
	jle	.L26

	movdqa	-29 * SIZE(X), %xmm2
	movdqa	-28 * SIZE(Y), %xmm3

	psrldq  $4, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $12, %xmm8
        por     %xmm8, %xmm0
	//movss	%xmm2, %xmm0
	//shufps	$0x39, %xmm0, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)

	psrldq  $12, %xmm1
	pshufd  $0xe4,%xmm3, %xmm9
	pslldq  $4, %xmm9
        por     %xmm9, %xmm1
	//movss	%xmm3, %xmm1
	//shufps	$0x93, %xmm3, %xmm1
	movdqa	%xmm1, -29 * SIZE(X)

	movdqa	%xmm2, %xmm0
	movdqa	%xmm3, %xmm1

	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L26:
	pshufd	$0x39, %xmm0, %xmm2
	pshufd	$0xff, %xmm0, %xmm0

	movq	%xmm2, -32 * SIZE(Y)
	movd	%xmm0, -30 * SIZE(Y)

	testq	$2, M
	jle	.L27

	movq	-29 * SIZE(X), %xmm0
	movq	-29 * SIZE(Y), %xmm1

	movq	%xmm0, -29 * SIZE(Y)
	movq	%xmm1, -29 * SIZE(X)

	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

.L27:
	testq	$1, M
	jle	.L29

	movd	-29 * SIZE(X), %xmm0
	movd	-29 * SIZE(Y), %xmm1

	movd	%xmm0, -29 * SIZE(Y)
	movd	%xmm1, -29 * SIZE(X)
	ALIGN_3

.L29:
	xorq	%rax,%rax

	RESTOREREGISTERS

#ifdef WINDOWS_ABI
	popq	%rbx
#endif

	ret
	ALIGN_3

.L30:
	testq	$1 * SIZE, X
	jne	.L40

	movdqa	 -34 * SIZE(X), %xmm0
	psrldq    $8, %xmm0
	//movq	 -32 * SIZE(X), %xmm0
	//movhps	 -32 * SIZE(X), %xmm0
	movdqa	 -32 * SIZE(Y), %xmm1

	movq	%xmm1, -32 * SIZE(X)
	//movlps	%xmm1, -32 * SIZE(X)
	subq	$2, M

	movq	M,  %rax
	sarq	$5, %rax
	jle	.L33
	ALIGN_4

.L31:
#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movdqa	-30 * SIZE(X), %xmm2
	movdqa	-28 * SIZE(Y), %xmm3

	punpcklqdq %xmm2, %xmm0
	//SHUFPD_1  %xmm2, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)
	psrldq     $8, %xmm2

	psrldq     $8, %xmm1
	punpcklqdq %xmm3, %xmm1
	//SHUFPD_1  %xmm3, %xmm1
	movdqa	%xmm1, -30 * SIZE(X)

	movdqa	-26 * SIZE(X), %xmm0
	movdqa	-24 * SIZE(Y), %xmm1

	punpcklqdq %xmm0, %xmm2
	//SHUFPD_1  %xmm0, %xmm2
	movdqa	%xmm2, -28 * SIZE(Y)
	psrldq     $8, %xmm0

	psrldq     $8, %xmm3
	punpcklqdq %xmm1, %xmm3
	//SHUFPD_1  %xmm1, %xmm3
	movdqa	%xmm3, -26 * SIZE(X)

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	movdqa	-22 * SIZE(X), %xmm2
	movdqa	-20 * SIZE(Y), %xmm3

	punpcklqdq %xmm2, %xmm0
	//SHUFPD_1  %xmm2, %xmm0
	movdqa	%xmm0, -24 * SIZE(Y)
	psrldq     $8, %xmm2

	psrldq     $8, %xmm1
	punpcklqdq %xmm3, %xmm1
	//SHUFPD_1  %xmm3, %xmm1
	movdqa	%xmm1, -22 * SIZE(X)

	movdqa	-18 * SIZE(X), %xmm0
	movdqa	-16 * SIZE(Y), %xmm1

	punpcklqdq %xmm0, %xmm2
	//SHUFPD_1  %xmm0, %xmm2
	movdqa	%xmm2, -20 * SIZE(Y)
	psrldq     $8, %xmm0

	psrldq     $8, %xmm3
	punpcklqdq %xmm1, %xmm3
	//SHUFPD_1  %xmm1, %xmm3
	movdqa	%xmm3, -18 * SIZE(X)

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE + 64) - PREOFFSET(X)
#endif

	movdqa	-14 * SIZE(X), %xmm2
	movdqa	-12 * SIZE(Y), %xmm3

	punpcklqdq %xmm2, %xmm0
	//SHUFPD_1  %xmm2, %xmm0
	movdqa	%xmm0, -16 * SIZE(Y)
	psrldq     $8, %xmm2

	psrldq     $8, %xmm1
	punpcklqdq %xmm3, %xmm1
	//SHUFPD_1  %xmm3, %xmm1
	movdqa	%xmm1, -14 * SIZE(X)

	movdqa	-10 * SIZE(X), %xmm0
	movdqa	 -8 * SIZE(Y), %xmm1

	punpcklqdq %xmm0, %xmm2
	//SHUFPD_1  %xmm0, %xmm2
	movdqa	%xmm2, -12 * SIZE(Y)
	psrldq     $8, %xmm0

	psrldq     $8, %xmm3
	punpcklqdq %xmm1, %xmm3
	//SHUFPD_1  %xmm1, %xmm3
	movdqa	%xmm3, -10 * SIZE(X)

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE + 64) - PREOFFSET(Y)
#endif

	movdqa	 -6 * SIZE(X), %xmm2
	movdqa	 -4 * SIZE(Y), %xmm3

	punpcklqdq %xmm2, %xmm0
	//SHUFPD_1  %xmm2, %xmm0
	movdqa	%xmm0,  -8 * SIZE(Y)
	psrldq     $8, %xmm2
	
	psrldq     $8, %xmm1
	punpcklqdq %xmm3, %xmm1
	//SHUFPD_1  %xmm3, %xmm1
	movdqa	%xmm1,  -6 * SIZE(X)

	movdqa	 -2 * SIZE(X), %xmm0
	movdqa	  0 * SIZE(Y), %xmm1

	punpcklqdq %xmm0, %xmm2
	//SHUFPD_1  %xmm0, %xmm2
	movdqa	%xmm2,  -4 * SIZE(Y)
	psrldq     $8, %xmm0

	psrldq     $8, %xmm3
	punpcklqdq %xmm1, %xmm3
	//SHUFPD_1  %xmm1, %xmm3
	movdqa	%xmm3,  -2 * SIZE(X)

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y

	decq	%rax
	jg	.L31
	ALIGN_3

.L33:
	testq	$16, M
	jle	.L34

	movdqa	-30 * SIZE(X), %xmm2
	movdqa	-28 * SIZE(Y), %xmm3

	punpcklqdq %xmm2, %xmm0
	//SHUFPD_1  %xmm2, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)
	psrldq     $8, %xmm2

	psrldq     $8, %xmm1
	punpcklqdq %xmm3, %xmm1
	//SHUFPD_1  %xmm3, %xmm1
	movdqa	%xmm1, -30 * SIZE(X)

	movdqa	-26 * SIZE(X), %xmm0
	movdqa	-24 * SIZE(Y), %xmm1

	punpcklqdq %xmm0, %xmm2
	//SHUFPD_1  %xmm0, %xmm2
	movdqa	%xmm2, -28 * SIZE(Y)
	psrldq     $8, %xmm0
	
	psrldq     $8, %xmm3
	punpcklqdq %xmm1, %xmm3
	//SHUFPD_1  %xmm1, %xmm3
	movdqa	%xmm3, -26 * SIZE(X)

	movdqa	-22 * SIZE(X), %xmm2
	movdqa	-20 * SIZE(Y), %xmm3

	punpcklqdq %xmm2, %xmm0
	//SHUFPD_1  %xmm2, %xmm0
	movdqa	%xmm0, -24 * SIZE(Y)
	psrldq     $8, %xmm2

	psrldq     $8, %xmm1
	punpcklqdq %xmm3, %xmm1
	//SHUFPD_1  %xmm3, %xmm1
	movdqa	%xmm1, -22 * SIZE(X)

	movdqa	-18 * SIZE(X), %xmm0
	movdqa	-16 * SIZE(Y), %xmm1

	punpcklqdq %xmm0, %xmm2
	//SHUFPD_1  %xmm0, %xmm2
	movdqa	%xmm2, -20 * SIZE(Y)
	psrldq     $8, %xmm0
	
	psrldq     $8, %xmm3
	punpcklqdq %xmm1, %xmm3
	//SHUFPD_1  %xmm1, %xmm3
	movdqa	%xmm3, -18 * SIZE(X)

	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
	ALIGN_3

.L34:
	testq	$8, M
	jle	.L35

	movdqa	-30 * SIZE(X), %xmm2
	movdqa	-28 * SIZE(Y), %xmm3

	punpcklqdq %xmm2, %xmm0
	//SHUFPD_1  %xmm2, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)
	psrldq     $8, %xmm2

	psrldq     $8, %xmm1
	punpcklqdq %xmm3, %xmm1
	//SHUFPD_1  %xmm3, %xmm1
	movdqa	%xmm1, -30 * SIZE(X)

	movdqa	-26 * SIZE(X), %xmm0
	movdqa	-24 * SIZE(Y), %xmm1

	punpcklqdq %xmm0, %xmm2
	//SHUFPD_1  %xmm0, %xmm2
	movdqa	%xmm2, -28 * SIZE(Y)
	psrldq     $8, %xmm0
	
	psrldq     $8, %xmm3
	punpcklqdq %xmm1, %xmm3
	//SHUFPD_1  %xmm1, %xmm3
	movdqa	%xmm3, -26 * SIZE(X)

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L35:
	testq	$4, M
	jle	.L36

	movdqa	-30 * SIZE(X), %xmm2
	movdqa	-28 * SIZE(Y), %xmm3

	// psrldq     $8, %xmm0
	punpcklqdq %xmm2, %xmm0
	//SHUFPD_1  %xmm2, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)

	psrldq     $8, %xmm1
	punpcklqdq %xmm3, %xmm1
	//SHUFPD_1  %xmm3, %xmm1
	movdqa	%xmm1, -30 * SIZE(X)

	pshufd  $0xee, %xmm2, %xmm0
	//movdqa	%xmm2, %xmm0
	pshufd  $0x4e, %xmm3, %xmm1
	//movdqa	%xmm3, %xmm1

	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L36:
	movq	%xmm0, -32 * SIZE(Y)
	//movhps	%xmm0, -32 * SIZE(Y)

	testq	$2, M
	jle	.L37

	movq	-30 * SIZE(X), %xmm0
	//movsd	-30 * SIZE(X), %xm0m
	movq	-30 * SIZE(Y), %xmm1
	//movsd	-30 * SIZE(Y), %xmm1

	movq	%xmm0, -30 * SIZE(Y)
	//movlps	%xmm0, -30 * SIZE(Y)
	movq	%xmm1, -30 * SIZE(X)
	//movlps	%xmm1, -30 * SIZE(X)

	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

.L37:
	testq	$1, M
	jle	.L39

	movd	-30 * SIZE(X), %xmm0
	movd	-30 * SIZE(Y), %xmm1

	movd	%xmm0, -30 * SIZE(Y)
	movd	%xmm1, -30 * SIZE(X)
	ALIGN_3

.L39:
	xorq	%rax,%rax

	RESTOREREGISTERS

#ifdef WINDOWS_ABI
	popq	%rbx
#endif

	ret
	ALIGN_3

.L40:
	movdqa	 -35 * SIZE(X), %xmm0
	movdqa	 -32 * SIZE(Y), %xmm1
	
	movd	%xmm1, -32 * SIZE(X)

	subq	$3, M

	movq	M,  %rax
	sarq	$5, %rax
	jle	.L43
	ALIGN_4

.L41:
#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movdqa	-31 * SIZE(X), %xmm2
	movdqa	-28 * SIZE(Y), %xmm3

	psrldq  $12, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $4, %xmm8
        por     %xmm8, %xmm0
	//movss	%xmm2, %xmm0
	//shufps	$0x93, %xmm2, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)


	psrldq  $4, %xmm1
	pshufd  $0xe4,%xmm3, %xmm9
	pslldq  $12, %xmm9
        por     %xmm9, %xmm1
	//movss	%xmm3, %xmm1
	//shufps	$0x39, %xmm1, %xmm1
	movdqa	%xmm1, -31 * SIZE(X)

	movdqa	-27 * SIZE(X), %xmm0
	movdqa	-24 * SIZE(Y), %xmm1

	psrldq  $12, %xmm2
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $4, %xmm8
        por     %xmm8, %xmm2
	//movss	%xmm0, %xmm2
	//shufps	$0x93, %xmm0, %xmm2
	movdqa	%xmm2, -28 * SIZE(Y)

	psrldq  $4, %xmm3
	pshufd  $0xe4,%xmm1, %xmm9
	pslldq  $12, %xmm9
        por     %xmm9, %xmm3
	//movss	%xmm1, %xmm3
	//shufps	$0x39, %xmm3, %xmm3
	movdqa	%xmm3, -27 * SIZE(X)

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	movdqa	-23 * SIZE(X), %xmm2
	movdqa	-20 * SIZE(Y), %xmm3

	psrldq  $12, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $4, %xmm8
        por     %xmm8, %xmm0
	//movss	%xmm2, %xmm0
	//shufps	$0x93, %xmm2, %xmm0
	movdqa	%xmm0, -24 * SIZE(Y)

	psrldq  $4, %xmm1
	pshufd  $0xe4,%xmm3, %xmm9
	pslldq  $12, %xmm9
        por     %xmm9, %xmm1
	//movss	%xmm3, %xmm1
	//shufps	$0x39, %xmm1, %xmm1
	movdqa	%xmm1, -23 * SIZE(X)

	movdqa	-19 * SIZE(X), %xmm0
	movdqa	-16 * SIZE(Y), %xmm1

	psrldq  $12, %xmm2
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $4, %xmm8
        por     %xmm8, %xmm2
	//movss	%xmm0, %xmm2
	//shufps	$0x93, %xmm0, %xmm2
	movdqa	%xmm2, -20 * SIZE(Y)

	psrldq  $4, %xmm3
	pshufd  $0xe4,%xmm1, %xmm9
	pslldq  $12, %xmm9
        por     %xmm9, %xmm3
	//movss	%xmm1, %xmm3
	//shufps	$0x39, %xmm3, %xmm3
	movdqa	%xmm3, -19 * SIZE(X)

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE + 64) - PREOFFSET(X)
#endif

	movdqa	-15 * SIZE(X), %xmm2
	movdqa	-12 * SIZE(Y), %xmm3

	psrldq  $12, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $4, %xmm8
        por     %xmm8, %xmm0
	//movss	%xmm2, %xmm0
	//shufps	$0x93, %xmm2, %xmm0
	movdqa	%xmm0, -16 * SIZE(Y)

	psrldq  $4, %xmm1
	pshufd  $0xe4,%xmm3, %xmm9
	pslldq  $12, %xmm9
        por     %xmm9, %xmm1
	//movss	%xmm3, %xmm1
	//shufps	$0x39, %xmm1, %xmm1
	movdqa	%xmm1, -15 * SIZE(X)

	movdqa	-11 * SIZE(X), %xmm0
	movdqa	 -8 * SIZE(Y), %xmm1

	psrldq  $12, %xmm2
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $4, %xmm8
        por     %xmm8, %xmm2
	//movss	%xmm0, %xmm2
	//shufps	$0x93, %xmm0, %xmm2
	movdqa	%xmm2, -12 * SIZE(Y)

	psrldq  $4, %xmm3
	pshufd  $0xe4,%xmm1, %xmm9
	pslldq  $12, %xmm9
        por     %xmm9, %xmm3
	//movss	%xmm1, %xmm3
	//shufps	$0x39, %xmm3, %xmm3
	movdqa	%xmm3, -11 * SIZE(X)

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE + 64) - PREOFFSET(Y)
#endif

	movdqa	 -7 * SIZE(X), %xmm2
	movdqa	 -4 * SIZE(Y), %xmm3

	psrldq  $12, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $4, %xmm8
        por     %xmm8, %xmm0
	//movss	%xmm2, %xmm0
	//shufps	$0x93, %xmm2, %xmm0
	movdqa	%xmm0,  -8 * SIZE(Y)

	psrldq  $4, %xmm1
	pshufd  $0xe4,%xmm3, %xmm9
	pslldq  $12, %xmm9
        por     %xmm9, %xmm1
	//movss	%xmm3, %xmm1
	//shufps	$0x39, %xmm1, %xmm1
	movdqa	%xmm1,  -7 * SIZE(X)

	movdqa	 -3 * SIZE(X), %xmm0
	movdqa	  0 * SIZE(Y), %xmm1

	psrldq  $12, %xmm2
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $4, %xmm8
        por     %xmm8, %xmm2
	//movss	%xmm0, %xmm2
	//shufps	$0x93, %xmm0, %xmm2
	movdqa	%xmm2,  -4 * SIZE(Y)

	psrldq  $4, %xmm3
	pshufd  $0xe4,%xmm1, %xmm9
	pslldq  $12, %xmm9
        por     %xmm9, %xmm3
	//movss	%xmm1, %xmm3
	//shufps	$0x39, %xmm3, %xmm3
	movdqa	%xmm3,  -3 * SIZE(X)

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y

	decq	%rax
	jg	.L41
	ALIGN_3

.L43:
	testq	$16, M
	jle	.L44

	movdqa	-31 * SIZE(X), %xmm2
	movdqa	-28 * SIZE(Y), %xmm3

	psrldq  $12, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $4, %xmm8
        por     %xmm8, %xmm0
	//movss	%xmm2, %xmm0
	//shufps	$0x93, %xmm2, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)

	psrldq  $4, %xmm1
	pshufd  $0xe4,%xmm3, %xmm9
	pslldq  $12, %xmm9
        por     %xmm9, %xmm1
	//movss	%xmm3, %xmm1
	//shufps	$0x39, %xmm1, %xmm1
	movdqa	%xmm1, -31 * SIZE(X)

	movdqa	-27 * SIZE(X), %xmm0
	movdqa	-24 * SIZE(Y), %xmm1

	psrldq  $12, %xmm2
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $4, %xmm8
        por     %xmm8, %xmm2
	//movss	%xmm0, %xmm2
	//shufps	$0x93, %xmm0, %xmm2
	movdqa	%xmm2, -28 * SIZE(Y)

	psrldq  $4, %xmm3
	pshufd  $0xe4,%xmm1, %xmm9
	pslldq  $12, %xmm9
        por     %xmm9, %xmm3
	//movss	%xmm1, %xmm3
	//shufps	$0x39, %xmm3, %xmm3
	movdqa	%xmm3, -27 * SIZE(X)

	movdqa	-23 * SIZE(X), %xmm2
	movdqa	-20 * SIZE(Y), %xmm3

	psrldq  $12, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $4, %xmm8
        por     %xmm8, %xmm0
	//movss	%xmm2, %xmm0
	//shufps	$0x93, %xmm2, %xmm0
	movdqa	%xmm0, -24 * SIZE(Y)

	psrldq  $4, %xmm1
	pshufd  $0xe4,%xmm3, %xmm9
	pslldq  $12, %xmm9
        por     %xmm9, %xmm1
	//movss	%xmm3, %xmm1
	//shufps	$0x39, %xmm1, %xmm1
	movdqa	%xmm1, -23 * SIZE(X)

	movdqa	-19 * SIZE(X), %xmm0
	movdqa	-16 * SIZE(Y), %xmm1

	psrldq  $12, %xmm2
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $4, %xmm8
        por     %xmm8, %xmm2
	//movss	%xmm0, %xmm2
	//shufps	$0x93, %xmm0, %xmm2
	movdqa	%xmm2, -20 * SIZE(Y)

	psrldq  $4, %xmm3
	pshufd  $0xe4,%xmm1, %xmm9
	pslldq  $12, %xmm9
        por     %xmm9, %xmm3
	//movss	%xmm1, %xmm3
	//shufps	$0x39, %xmm3, %xmm3
	movdqa	%xmm3, -19 * SIZE(X)

	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
	ALIGN_3

.L44:
	testq	$8, M
	jle	.L45

	movdqa	-31 * SIZE(X), %xmm2
	movdqa	-28 * SIZE(Y), %xmm3


	psrldq  $12, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $4, %xmm8
        por     %xmm8, %xmm0
	//movss	%xmm2, %xmm0
	//shufps	$0x93, %xmm2, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)

	psrldq  $4, %xmm1
	pshufd  $0xe4,%xmm3, %xmm9
	pslldq  $12, %xmm9
        por     %xmm9, %xmm1
	//movss	%xmm3, %xmm1
	//shufps	$0x39, %xmm1, %xmm1
	movdqa	%xmm1, -31 * SIZE(X)

	movdqa	-27 * SIZE(X), %xmm0
	movdqa	-24 * SIZE(Y), %xmm1

	psrldq  $12, %xmm2
	pshufd  $0xe4,%xmm0, %xmm8
	pslldq  $4, %xmm8
        por     %xmm8, %xmm2
	//movss	%xmm0, %xmm2
	//shufps	$0x93, %xmm0, %xmm2
	movdqa	%xmm2, -28 * SIZE(Y)

	psrldq  $4, %xmm3
	pshufd  $0xe4,%xmm1, %xmm9
	pslldq  $12, %xmm9
        por     %xmm9, %xmm3
	//movss	%xmm1, %xmm3
	//shufps	$0x39, %xmm3, %xmm3
	movdqa	%xmm3, -27 * SIZE(X)

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L45:
	testq	$4, M
	jle	.L46

	movdqa	-31 * SIZE(X), %xmm2
	movdqa	-28 * SIZE(Y), %xmm3

	psrldq  $12, %xmm0
	pshufd  $0xe4,%xmm2, %xmm8
	pslldq  $4, %xmm8
        por     %xmm8, %xmm0
	//movss	%xmm2, %xmm0
	//shufps	$0x93, %xmm2, %xmm0
	movdqa	%xmm0, -32 * SIZE(Y)

	psrldq  $4, %xmm1
	pshufd  $0xe4,%xmm3, %xmm9
	pslldq  $12, %xmm9
        por     %xmm9, %xmm1
	//movss	%xmm3, %xmm1
	//shufps	$0x39, %xmm1, %xmm1
	movdqa	%xmm1, -31 * SIZE(X)

	movdqa	%xmm2, %xmm0
	movdqa	%xmm3, %xmm1

	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L46:
	movsd	-31 * SIZE(X), %xmm2

	pshufd	$0x39, %xmm1, %xmm1
	movlps	%xmm1, -31 * SIZE(X)

	pshufd	$0xff, %xmm0, %xmm0

	movd	%xmm0, -32 * SIZE(Y)
	movlps	%xmm2, -31 * SIZE(Y)

	addq	$3 * SIZE, X
	addq	$3 * SIZE, Y

	testq	$2, M
	jle	.L47

	movsd	-32 * SIZE(X), %xmm0
	movsd	-32 * SIZE(Y), %xmm1

	movlps	%xmm0, -32 * SIZE(Y)
	movlps	%xmm1, -32 * SIZE(X)

	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

.L47:
	testq	$1, M
	jle	.L49

	movd	-32 * SIZE(X), %xmm0
	movd	-32 * SIZE(Y), %xmm1

	movd	%xmm0, -32 * SIZE(Y)
	movd	%xmm1, -32 * SIZE(X)
	ALIGN_3

.L49:
	xorq	%rax,%rax

	RESTOREREGISTERS

#ifdef WINDOWS_ABI
	popq	%rbx
#endif

	ret
	ALIGN_3

.L50:
	movq	M,  %rax
	sarq	$3, %rax
	jle	.L55
	ALIGN_3

.L51:
	movd	(X), %xmm0
	movd	(Y), %xmm1

	movd	%xmm1, (X)
	addq	INCX, X
	movd	%xmm0, (Y)
	addq	INCY, Y

	movd	(X), %xmm0
	movd	(Y), %xmm1

	movd	%xmm1, (X)
	addq	INCX, X
	movd	%xmm0, (Y)
	addq	INCY, Y

	movd	(X), %xmm0
	movd	(Y), %xmm1

	movd	%xmm1, (X)
	addq	INCX, X
	movd	%xmm0, (Y)
	addq	INCY, Y

	movd	(X), %xmm0
	movd	(Y), %xmm1

	movd	%xmm1, (X)
	addq	INCX, X
	movd	%xmm0, (Y)
	addq	INCY, Y

	movd	(X), %xmm0
	movd	(Y), %xmm1

	movd	%xmm1, (X)
	addq	INCX, X
	movd	%xmm0, (Y)
	addq	INCY, Y

	movd	(X), %xmm0
	movd	(Y), %xmm1

	movd	%xmm1, (X)
	addq	INCX, X
	movd	%xmm0, (Y)
	addq	INCY, Y

	movd	(X), %xmm0
	movd	(Y), %xmm1

	movd	%xmm1, (X)
	addq	INCX, X
	movd	%xmm0, (Y)
	addq	INCY, Y

	movd	(X), %xmm0
	movd	(Y), %xmm1

	movd	%xmm1, (X)
	addq	INCX, X
	movd	%xmm0, (Y)
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
	movd	(Y), %xmm1

	movd	%xmm1, (X)
	movd	%xmm0, (Y)

	addq	INCX, X
	addq	INCY, Y
	decq	%rax
	jg	.L56
	ALIGN_3

.L57:
	xorq	%rax, %rax

	RESTOREREGISTERS

#ifdef WINDOWS_ABI
	popq	%rbx
#endif

	ret

/*	EPILOGUE
*/

.LFE2:
	.size	iswap_kh, .-iswap_kh
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
