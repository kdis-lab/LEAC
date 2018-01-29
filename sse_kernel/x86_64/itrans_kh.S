/*! \file iaysxpy_kh.S	
 * This file is part of the LEAC.
 *
 * Implementation of the 
 *	
 *	y = y + a,
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


#define M	ARG1  /* rdi */
#define ALPHA32 %ecx  
#define X	ARG5  /* r8  */
#define INCX	ARG6  /* r9  */
#define Y	ARG3  /* rdx */
#define INCY	%r10  /* r10 */
	
#define XX	%r10		
#define ALPHA	%xmm15	
#define I	%rax

#define SIZE    4
	
#include "l1param.h"

        .file	"itrans_kh.S"
	.text
	.p2align 4,,15
.globl itrans_kh
	.type	itrans_kh, @function
itrans_kh:
.LFB2:	

	SAVEREGISTERS
	
	lea	(, INCX, SIZE), INCX
			
	testq	M, M
	jle	.L999

	movd    ALPHA32, ALPHA
	pshufd  $0x0,ALPHA,ALPHA

	cmpq	$3, M
	jle	.L14
	
	cmpq	$SIZE, INCX
	jne	.L50

/* INCX == 1 */
	
	testq	$4, X		# aligned for double word?
	je	.L05

	movl    0 * SIZE(X), %esi
	addl    ALPHA32, %esi
	movl    %esi, 0 * SIZE(X)
	addq	$SIZE, X
	decq	M
	jle	.L999
	ALIGN_3

.L05:
	testq	$8, X		# aligned for quad word?
	je	.L06

	movq   0 * SIZE(X), %xmm0 
	paddd  ALPHA, %xmm0
	movq   %xmm0, 0 * SIZE(X) 
	
	addq	$2 * SIZE, X
	subq	$2, M
	jle	.L999
	ALIGN_3

.L06:
	movq	M,  I
	sarq	$4, I
	jle	.L12
	ALIGN_4

.L11:
#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif


	movdqa	0 * SIZE(X),  %xmm0
	movdqa	4 * SIZE(X),  %xmm1
	movdqa	8 * SIZE(X),  %xmm2
	movdqa	12 * SIZE(X), %xmm3

	paddd  ALPHA, %xmm0
	paddd  ALPHA, %xmm1
	paddd  ALPHA, %xmm2
	paddd  ALPHA, %xmm3

	movdqa	%xmm0,  0 * SIZE(X)
	movdqa	%xmm1,  4 * SIZE(X)
	movdqa	%xmm2,  8 * SIZE(X)
	movdqa	%xmm3, 12 * SIZE(X)

	
	addq	$16 * SIZE, X
	decq	I
	jg	.L11
	ALIGN_4

.L12:
	testq	$15, M
	je	.L999
	
	testq	$8, M
	je	.L13

	movdqa	0 * SIZE(X),  %xmm0
	movdqa	4 * SIZE(X),  %xmm1
	paddd  ALPHA, %xmm0
	paddd  ALPHA, %xmm1
	movdqa	%xmm0,  0 * SIZE(X)
	movdqa	%xmm1,  4 * SIZE(X)
	
	addq	$8 * SIZE, X
	ALIGN_3

.L13:
	testq  $4, M
	je    .L14

	movdqa  0 * SIZE(X),  %xmm0
	paddd  ALPHA, %xmm0
	movdqa	%xmm0,  0 * SIZE(X)
	
	addq	$4 * SIZE, X
	ALIGN_3

.L14:
	testq  $2, M
	je    .L15

	movq	0 * SIZE(X), %xmm0
	paddd  ALPHA, %xmm0
	movq   %xmm0, 0 * SIZE(X)
	
	addq	$2 * SIZE, X
	ALIGN_3

.L15:
	testq  $1, M
	je    .L999

	movl    0 * SIZE(X), %esi
	addl    ALPHA32, %esi
	movl    %esi, 0 * SIZE(X)
	jmp	.L999
	ALIGN_4

/* incx != 1 */
.L50:
	movq	M,  I		# rcx = n
	sarq	$3, I		# (n >> 3)
	jle	.L52
	ALIGN_4

.L51:
#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movd	%xmm0, (X)
	addq	INCX, X
	movd	%xmm0, (X)
	addq	INCX, X
	movd	%xmm0, (X)
	addq	INCX, X
	movd	%xmm0, (X)
	addq	INCX, X
	movd	%xmm0, (X)
	addq	INCX, X
	movd	%xmm0, (X)
	addq	INCX, X
	movd	%xmm0, (X)
	addq	INCX, X
	movd	%xmm0, (X)
	addq	INCX, X

	decq	I
	jg	.L51
	ALIGN_4

.L52:
	testq	$7, M
	je	.L999

	testq	$4, M
	je	.L53

	movd	%xmm0, (X)
	addq	INCX, X
	movd	%xmm0, (X)
	addq	INCX, X
	movd	%xmm0, (X)
	addq	INCX, X
	movd	%xmm0, (X)
	addq	INCX, X
	ALIGN_3

.L53:
	testq	$2, M
	je	.L54

	movd	%xmm0, (X)
	addq	INCX, X
	movd	%xmm0, (X)
	addq	INCX, X
	ALIGN_3

.L54:
	testq	$1, M
	je	.L999

	movd	%xmm0, (X)
	jmp	.L999
	ALIGN_4

.L999:
	xorq	%rax, %rax

	RESTOREREGISTERS

	ret

.LFE2:
	.size	itrans_kh, .-itrans_kh
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

	
