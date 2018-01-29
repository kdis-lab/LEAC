/*! \file dcopya_k.S	
 * This file is part of the LEAC.
 *
 * Implementation of the 
 *	
 *	y <-- a, (double <--double) 
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

#define ARG1	%rdi  /*N*/
#define ARG2	%rsi  /*X*/
#define ARG3	%rdx  /*INCX*/
#define ARG4	%rcx
#define ARG5	%r8
#define ARG6	%r9
	
#define M	ARG1
#define X	ARG2
#define INCX	ARG3

#define XX	%r10
#define I	%rax

#define SIZE    8
	
#include "l1param.h"

        .file	"dcopya_k.S"
	.text
	.p2align 4,,15
.globl dcopya_k
	.type	dcopya_k, @function
dcopya_k:
.LFB2:
	
	SAVEREGISTERS

	leaq	(, INCX, SIZE), INCX
	
	testq	M, M 
	jle	.L999

	movlhps %xmm0, %xmm0
	
	cmpq	$SIZE, INCX
	jne	.L50

/* INCX == 1 */
	testq	$15, X		# aligned for quad word?
	je	.L05

	movsd	%xmm0, 0 * SIZE(X)
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

	movaps	%xmm0,  0 * SIZE(X)
	movaps	%xmm0,  2 * SIZE(X)
	movaps	%xmm0,  4 * SIZE(X)
	movaps	%xmm0,  6 * SIZE(X)

#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE + 64) - PREOFFSET(X)
#endif

	movaps	%xmm0,  8 * SIZE(X)
	movaps	%xmm0, 10 * SIZE(X)
	movaps	%xmm0, 12 * SIZE(X)
	movaps	%xmm0, 14 * SIZE(X)

	addq	$16 * SIZE, X
	decq	I
	jg	.L11
	ALIGN_4

.L12:
	testq	$15, M
	je	.L999
	testq	$8, M
	je	.L13

	movaps	%xmm0,  0 * SIZE(X)
	movaps	%xmm0,  2 * SIZE(X)
	movaps	%xmm0,  4 * SIZE(X)
	movaps	%xmm0,  6 * SIZE(X)
	addq	$8 * SIZE, X
	ALIGN_3

.L13:
	testq  $4, M
	je    .L14

	movaps	%xmm0,  0 * SIZE(X)
	movaps	%xmm0,  2 * SIZE(X)
	addq	$4 * SIZE, X
	ALIGN_3

.L14:
	testq  $2, M
	je    .L15

	movaps	%xmm0,  0 * SIZE(X)
	addq	$2 * SIZE, X
	ALIGN_3

.L15:
	testq  $1, M
	je    .L999

	movsd	%xmm0,  0 * SIZE(X)
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

	movsd	%xmm0, (X)
	addq	INCX, X
	movsd	%xmm0, (X)
	addq	INCX, X
	movsd	%xmm0, (X)
	addq	INCX, X
	movsd	%xmm0, (X)
	addq	INCX, X
	movsd	%xmm0, (X)
	addq	INCX, X
	movsd	%xmm0, (X)
	addq	INCX, X
	movsd	%xmm0, (X)
	addq	INCX, X
	movsd	%xmm0, (X)
	addq	INCX, X

	decq	I
	jg	.L51
	ALIGN_4

.L52:
	testq	$7, M
	je	.L999

	testq	$4, M
	je	.L53

	movsd	%xmm0, (X)
	addq	INCX, X
	movsd	%xmm0, (X)
	addq	INCX, X
	movsd	%xmm0, (X)
	addq	INCX, X
	movsd	%xmm0, (X)
	addq	INCX, X
	ALIGN_3

.L53:
	testq	$2, M
	je	.L54

	movsd	%xmm0, (X)
	addq	INCX, X
	movsd	%xmm0, (X)
	addq	INCX, X
	ALIGN_3

.L54:
	testq	$1, M
	je	.L999

	movsd	%xmm0, (X)
	jmp	.L999
	ALIGN_4

.L999:
	xorq	%rax, %rax

	RESTOREREGISTERS

	ret

.LFE2:
	.size	dcopya_k, .-dcopya_k
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

