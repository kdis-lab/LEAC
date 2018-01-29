/*! \file dsum_kh.S	
 * This file is part of the LEAC.
 *
 * Implementation of the 
 *	
 *	sum x,
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

#define N	ARG1	/* rdi */
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
#define LOAD(OFFSET, ADDR, REG)		movaps	OFFSET(ADDR), REG
#endif
	
/*	PROLOGUE
	PROFCODE
*/

	.file	"dsum_kh.S"
	.text
	.p2align 4,,15
.globl dsum_kh
	.type	dsum_kh, @function
dsum_kh:
.LFB2:


	SAVEREGISTERS

	leaq	(, INCX, SIZE), INCX	

	xorps	%xmm0, %xmm0
	xorps	%xmm1, %xmm1
	xorps	%xmm2, %xmm2
	xorps	%xmm3, %xmm3
	
	cmpq	$0, N
	jle	.L999

	cmpq	$SIZE, INCX
	jne	.L50
	
	subq	$-16 * SIZE, X

	testq	$SIZE, X
	je	.L10

	addsd	-16 * SIZE(X), %xmm0
	addq	$1 * SIZE, X
	decq	N
	ALIGN_2

.L10:

	movq	N,  %rax
	sarq	$4, %rax
	jle	.L14

	addpd	-16 * SIZE(X), %xmm0
	addpd	-14 * SIZE(X), %xmm1
	addpd	-12 * SIZE(X), %xmm2
	addpd	-10 * SIZE(X), %xmm3

	decq	%rax
	jle	.L12
	ALIGN_3

.L11:
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	addpd	 -8 * SIZE(X), %xmm0
	addpd	 -6 * SIZE(X), %xmm1

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	addpd	 -4 * SIZE(X), %xmm2
	addpd	 -2 * SIZE(X), %xmm3

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(X)
#endif

	addpd	  0 * SIZE(X), %xmm0
	addpd	  2 * SIZE(X), %xmm1

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(X)
#endif


	addpd	  4 * SIZE(X), %xmm2
	addpd	  6 * SIZE(X), %xmm3

	subq	$-16 * SIZE, X

	decq	%rax
	jg	.L11
	ALIGN_3

.L12:

	addpd	 -8 * SIZE(X), %xmm0
	addpd	 -6 * SIZE(X), %xmm1

	addpd	 -4 * SIZE(X), %xmm2
	addpd	 -2 * SIZE(X), %xmm3

	subq	$-16 * SIZE, X
	ALIGN_3

.L14:
	testq	$15, N
	jle	.L999

	testq	$8, N
	jle	.L15

	addpd	-16 * SIZE(X), %xmm0
	addpd	-14 * SIZE(X), %xmm1
	addpd	-12 * SIZE(X), %xmm2
	addpd	-10 * SIZE(X), %xmm3

	addq	$8 * SIZE, X
	ALIGN_3

.L15:
	testq	$4, N
	jle	.L16

	addpd	-16 * SIZE(X), %xmm0
	addpd	-14 * SIZE(X), %xmm1


	addq	$4 * SIZE, X
	ALIGN_3

.L16:
	testq	$2, N
	jle	.L17

	addpd	-16 * SIZE(X), %xmm0

	addq	$2 * SIZE, X
	
	ALIGN_3

.L17:
	testq	$1, N
	jle	.L999

	addsd	-16 * SIZE(X), %xmm0
	jmp	.L999
	ALIGN_3

.L50:
	movq	N,  %rax
	sarq	$2, %rax
	jle	.L55
	ALIGN_3

.L53:
	addsd	(X), %xmm0
	addq	INCX, X
	addsd	(X), %xmm1
	addq	INCX, X
	addsd	(X), %xmm2
	addq	INCX, X
	addsd	(X), %xmm3
	addq	INCX, X
	
	decq	%rax
	jg	.L53
	ALIGN_3

.L55:
	movq	N,  %rax
	andq	$3, %rax
	jle	.L999
	ALIGN_3

.L56:
	addsd	(X), %xmm0
	addq	INCX, X
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
	
/*	EPILOGUE
*/
.LFE2:
	.size	dsum_kh, .-dsum_kh
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
