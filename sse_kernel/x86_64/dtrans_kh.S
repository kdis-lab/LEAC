/*! \file iaysxpy_kh.S	
 * This file is part of the LEAC.
 *
 * Implementation of the 
 *	
 *	y = y + a,
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

#ifndef WINDOWS_ABI
#define M	ARG1
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
#define I	%rax
	
#define SIZE    8
	
#include "l1param.h"

        .file	"dtrans_kh.S"
	.text
	.p2align 4,,15
.globl dtrans_kh
	.type	dtrans_kh, @function
dtrans_kh:
.LFB2:

#ifndef WINDOWS_ABI
#ifndef XDOUBLE
	movq	 8(%rsp), INCY
#else
	movq	24(%rsp), INCY
#endif
	movaps	%xmm0,  ALPHA
#else
	movaps	%xmm3,  ALPHA

	movq	40(%rsp), X
	movq	48(%rsp), INCX
	movq	56(%rsp), Y
	movq	64(%rsp), INCY
#endif
	
	SAVEREGISTERS

	leaq	(, INCX, SIZE), INCX
	
	testq	M, M 
	jle	.L999

	movlhps ALPHA, ALPHA
	
	cmpq	$SIZE, INCX
	jne	.L50

/* INCX == 1 */
	testq	$15, X		# aligned for quad word?
	je	.L05

	movsd   0 * SIZE(X), %xmm0
	addsd   ALPHA,%xmm0
	movsd   %xmm0, 0 * SIZE(X)
	addq	$SIZE, X
	decq	M
	jle	.L999
	ALIGN_3
.L05:

/* Aligned Mode */
	movq	M,  I	# rcx = n
	sarq	$4, I
	jle	.L12
	ALIGN_4

.L11:
#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif
	movaps	0 * SIZE(X), %xmm0
	movaps	2 * SIZE(X), %xmm1
	movaps	4 * SIZE(X), %xmm2
	movaps	6 * SIZE(X), %xmm3
	
	addpd	ALPHA, %xmm0
	addpd	ALPHA, %xmm1
	addpd	ALPHA, %xmm2
	addpd	ALPHA, %xmm3
	
	movaps	%xmm0,  0 * SIZE(X)
	movaps	%xmm1,  2 * SIZE(X)
	movaps	%xmm2,  4 * SIZE(X)
	movaps	%xmm3,  6 * SIZE(X)

#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE + 64) - PREOFFSET(X)
#endif

	movaps	8 * SIZE(X), %xmm0
	movaps	10 * SIZE(X), %xmm1
	movaps	12 * SIZE(X), %xmm2
	movaps	14 * SIZE(X), %xmm3
	
	addpd	ALPHA, %xmm0
	addpd	ALPHA, %xmm1
	addpd	ALPHA, %xmm2
	addpd	ALPHA, %xmm3
	
	movaps	%xmm0,  8 * SIZE(X)
	movaps	%xmm1,  10 * SIZE(X)
	movaps	%xmm2,  12 * SIZE(X)
	movaps	%xmm3,  14 * SIZE(X)
	
	addq	$16 * SIZE, X
	decq	I
	jg	.L11
	ALIGN_4

.L12:
	testq	$15, M
	je	.L999
	testq	$8, M
	je	.L13

	movaps	0 * SIZE(X), %xmm0
	movaps	2 * SIZE(X), %xmm1
	movaps	4 * SIZE(X), %xmm2
	movaps	6 * SIZE(X), %xmm3
	
	addpd	ALPHA, %xmm0
	addpd	ALPHA, %xmm1
	addpd	ALPHA, %xmm2
	addpd	ALPHA, %xmm3
	
	movaps	%xmm0,  0 * SIZE(X)
	movaps	%xmm1,  2 * SIZE(X)
	movaps	%xmm2,  4 * SIZE(X)
	movaps	%xmm3,  6 * SIZE(X)
	
	addq	$8 * SIZE, X
	ALIGN_3

.L13:
	testq  $4, M
	je    .L14

	movaps	0 * SIZE(X), %xmm0
	movaps	2 * SIZE(X), %xmm1
	
	addpd	ALPHA, %xmm0
	addpd	ALPHA, %xmm1
	
	movaps	%xmm0,  0 * SIZE(X)
	movaps	%xmm1,  2 * SIZE(X)
	
	addq	$4 * SIZE, X
	ALIGN_3

.L14:
	testq  $2, M
	je    .L15

	movaps	0 * SIZE(X), %xmm0
	addpd	ALPHA, %xmm0
	movaps	%xmm0,  0 * SIZE(X)
	
	addq	$2 * SIZE, X
	ALIGN_3

.L15:
	testq  $1, M
	je    .L999

	movsd   0 * SIZE(X), %xmm0
	addsd   ALPHA,%xmm0
	movsd   %xmm0, 0 * SIZE(X)
	
	jmp	.L999
	ALIGN_4

.L50:
	movq	M,  I
	sarq	$3, I
	jle	.L52
	ALIGN_4

.L51:
#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movsd   (X), %xmm0
	addsd	ALPHA, %xmm0
	movsd   %xmm0, (X)
	addq	INCX, X

	movsd   (X), %xmm0
	addsd	ALPHA, %xmm0
	movsd   %xmm0, (X)
	addq	INCX, X

	movsd   (X), %xmm0
	addsd	ALPHA, %xmm0
	movsd   %xmm0, (X)
	addq	INCX, X

	movsd   (X), %xmm0
	addsd	ALPHA, %xmm0
	movsd   %xmm0, (X)
	addq	INCX, X

	movsd   (X), %xmm0
	addsd	ALPHA, %xmm0
	movsd   %xmm0, (X)
	addq	INCX, X

	movsd   (X), %xmm0
	addsd	ALPHA, %xmm0
	movsd   %xmm0, (X)
	addq	INCX, X

	movsd   (X), %xmm0
	addsd	ALPHA, %xmm0
	movsd   %xmm0, (X)
	addq	INCX, X

	movsd   (X), %xmm0
	addsd	ALPHA, %xmm0
	movsd   %xmm0, (X)
	addq	INCX, X

	decq	I
	jg	.L51
	ALIGN_4

.L52:
	testq	$7, M
	je	.L999

	testq	$4, M
	je	.L53

	movsd   (X), %xmm0
	addsd	ALPHA, %xmm0
	movsd   %xmm0, (X)
	addq	INCX, X
	
	movsd   (X), %xmm0
	addsd	ALPHA, %xmm0
	movsd   %xmm0, (X)
	addq	INCX, X
	
	movsd   (X), %xmm0
	addsd	ALPHA, %xmm0
	movsd   %xmm0, (X)
	addq	INCX, X
	
	movsd   (X), %xmm0
	addsd	ALPHA, %xmm0
	movsd   %xmm0, (X)
	addq	INCX, X
	ALIGN_3

.L53:
	testq	$2, M
	je	.L54

	movsd   (X), %xmm0
	addsd	ALPHA, %xmm0
	movsd   %xmm0, (X)
	addq	INCX, X
	
	movsd   (X), %xmm0
	addsd	ALPHA, %xmm0
	movsd   %xmm0, (X)
	addq	INCX, X
	ALIGN_3

.L54:
	testq	$1, M
	je	.L999

	movsd   (X), %xmm0
	addsd	ALPHA, %xmm0
	movsd   %xmm0, (X)
	
	jmp	.L999
	ALIGN_4

.L999:
	xorq	%rax, %rax

	RESTOREREGISTERS

	ret

.LFE2:
	.size	dtrans_kh, .-dtrans_kh
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

