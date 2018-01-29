/*! \file isum_kh.S	
 * This file is part of the LEAC.
 *
 * Implementation of the 
 *	
 *	sum x, (long int <-- sum int)
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

#define SIZE    4

#define ZERO1   %xmm14
#define ZERO2   %xmm15
	
#include "l1param.h"

        .file	"isum_kh.S"
	.text
	.p2align 4,,15
.globl isum_kh
	.type	isum_kh, @function
isum_kh:
.LFB2:

	SAVEREGISTERS
	
	leaq	(, INCX, SIZE), INCX	

	pxor	%xmm0, %xmm0
	pxor	%xmm1, %xmm1
	pxor	%xmm2, %xmm2
	pxor	%xmm3, %xmm3

	cmpq	$0, N
	jle	.L999


	cmpq	$SIZE, INCX
	jne	.L50

	subq	$-32 * SIZE, X

	cmpq	$3, N
	jle	.L17

        testq	$SIZE, X
	je	.L05
        
        movd	-32 * SIZE(X), %xmm4
	
	pxor      %xmm8, %xmm8
	pcmpgtd	  %xmm4, %xmm8
	punpckldq %xmm8, %xmm4
	
	paddq	%xmm4, %xmm0
     
     	addq	$1 * SIZE, X
	decq	N
	ALIGN_2

.L05:
	testq	$2 * SIZE, X
	je	.L10

	movq	-32 * SIZE(X), %xmm4
	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm4, %xmm8
	punpckldq %xmm8, %xmm4
	
	paddq	%xmm4, %xmm0

	addq	$2 * SIZE, X
	subq	$2, N
	jle	.L999
	ALIGN_2

.L10:		
 
	movq	N,  %rax
	sarq	$5, %rax
	jle	.L14

	movdqa	-32 * SIZE(X), %xmm4
	movdqa	-28 * SIZE(X), %xmm5
	movdqa	-24 * SIZE(X), %xmm6
	movdqa	-20 * SIZE(X), %xmm7

	decq	%rax
	jle	.L12

	ALIGN_3

.L11:
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif
	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm4, %xmm8
	movdqa	  %xmm4, %xmm9
	punpckldq %xmm8, %xmm4
	punpckhdq %xmm8, %xmm9
	
        paddq	%xmm4, %xmm0
        paddq   %xmm9, %xmm1
	
	movdqa	-16 * SIZE(X), %xmm4	

	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm5, %xmm8
	movdqa	  %xmm5, %xmm9
	punpckldq %xmm8, %xmm5
	punpckhdq %xmm8, %xmm9
	
        paddq	%xmm5, %xmm2
        paddq   %xmm9, %xmm3
	
	movdqa	-12 * SIZE(X), %xmm5

	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm6, %xmm8
	movdqa	  %xmm6, %xmm9
	punpckldq %xmm8, %xmm6
	punpckhdq %xmm8, %xmm9
	
        paddq	%xmm6, %xmm0
        paddq   %xmm9, %xmm1
	
	movdqa	-8 * SIZE(X), %xmm6

	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm7, %xmm8
	movdqa	  %xmm7, %xmm9
	punpckldq %xmm8, %xmm7
	punpckhdq %xmm8, %xmm9
	
        paddq	%xmm7, %xmm2
        paddq   %xmm9, %xmm3
	
	movdqa	-4 * SIZE(X), %xmm7
	
#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(X)
#endif

	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm4, %xmm8
	movdqa	  %xmm4, %xmm9
	punpckldq %xmm8, %xmm4
	punpckhdq %xmm8, %xmm9
	
        paddq	%xmm4, %xmm0
        paddq   %xmm9, %xmm1

	movdqa	 0 * SIZE(X), %xmm4

	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm5, %xmm8
	movdqa	  %xmm5, %xmm9
	punpckldq %xmm8, %xmm5
	punpckhdq %xmm8, %xmm9
	
        paddq	%xmm5, %xmm2
        paddq   %xmm9, %xmm3
	
	movdqa	4 * SIZE(X), %xmm5

	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm6, %xmm8
	movdqa	  %xmm6, %xmm9
	punpckldq %xmm8, %xmm6
	punpckhdq %xmm8, %xmm9
	
        paddq	%xmm6, %xmm0
        paddq   %xmm9, %xmm1
	
    	movdqa	 8 * SIZE(X), %xmm6

	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm7, %xmm8
	movdqa	  %xmm7, %xmm9
	punpckldq %xmm8, %xmm7
	punpckhdq %xmm8, %xmm9
	
        paddq	%xmm7, %xmm2
        paddq   %xmm9, %xmm3
	
        movdqa	12 * SIZE(X), %xmm7
	
	subq	$-32 * SIZE, X
	
	decq	%rax
	jg	.L11
	ALIGN_3

.L12:
	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm4, %xmm8
	movdqa	  %xmm4, %xmm9
	punpckldq %xmm8, %xmm4
	punpckhdq %xmm8, %xmm9
	
        paddq	%xmm4, %xmm0
        paddq   %xmm9, %xmm1
	
	movdqa	-16 * SIZE(X), %xmm4

	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm5, %xmm8
	movdqa	  %xmm5, %xmm9
	punpckldq %xmm8, %xmm5
	punpckhdq %xmm8, %xmm9
	
        paddq	%xmm5, %xmm2
        paddq   %xmm9, %xmm3
	
	movdqa	-12 * SIZE(X), %xmm5

	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm6, %xmm8
	movdqa	  %xmm6, %xmm9
	punpckldq %xmm8, %xmm6
	punpckhdq %xmm8, %xmm9
	
        paddq	%xmm6, %xmm0
        paddq   %xmm9, %xmm1
	
	movdqa	 -8 * SIZE(X), %xmm6

	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm7, %xmm8
	movdqa	  %xmm7, %xmm9
	punpckldq %xmm8, %xmm7
	punpckhdq %xmm8, %xmm9
	
        paddq	%xmm7, %xmm2
        paddq   %xmm9, %xmm3
	
	movdqa	 -4 * SIZE(X), %xmm7
	
	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm4, %xmm8
	movdqa	  %xmm4, %xmm9
	punpckldq %xmm8, %xmm4
	punpckhdq %xmm8, %xmm9
	
        paddq	%xmm4, %xmm0
        paddq   %xmm9, %xmm1
	
	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm5, %xmm8
	movdqa	  %xmm5, %xmm9
	punpckldq %xmm8, %xmm5
	punpckhdq %xmm8, %xmm9
	
        paddq	%xmm5, %xmm2
        paddq   %xmm9, %xmm3

	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm6, %xmm8
	movdqa	  %xmm6, %xmm9
	punpckldq %xmm8, %xmm6
	punpckhdq %xmm8, %xmm9
	
        paddq	%xmm6, %xmm0
        paddq   %xmm9, %xmm1

	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm7, %xmm8
	movdqa	  %xmm7, %xmm9
	punpckldq %xmm8, %xmm7
	punpckhdq %xmm8, %xmm9
	
        paddq	%xmm7, %xmm2
        paddq   %xmm9, %xmm3
	
	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y

	ALIGN_3

.L14:
	testq	$31, N
	jle	.L999

	testq	$16, N
	jle	.L15

	movdqa	-32 * SIZE(X), %xmm4
	movdqa	-28 * SIZE(X), %xmm5

	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm4, %xmm8
	movdqa	  %xmm4, %xmm9
	punpckldq %xmm8, %xmm4
	punpckhdq %xmm8, %xmm9
	
        paddq	%xmm4, %xmm0
        paddq   %xmm9, %xmm1
	
	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm5, %xmm8
	movdqa	  %xmm5, %xmm9
	punpckldq %xmm8, %xmm5
	punpckhdq %xmm8, %xmm9
	
        paddq	%xmm5, %xmm2
        paddq   %xmm9, %xmm3
	
	movdqa	-24 * SIZE(X), %xmm4
	movdqa	-20 * SIZE(X), %xmm5

	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm4, %xmm8
	movdqa	  %xmm4, %xmm9
	punpckldq %xmm8, %xmm4
	punpckhdq %xmm8, %xmm9
	
        paddq	%xmm4, %xmm0
        paddq   %xmm9, %xmm1

	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm5, %xmm8
	movdqa	  %xmm5, %xmm9
	punpckldq %xmm8, %xmm5
	punpckhdq %xmm8, %xmm9
	
        paddq	%xmm5, %xmm2
        paddq   %xmm9, %xmm3
		
	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
	ALIGN_3

.L15:
	testq	$8, N
	jle	.L16

	movdqa	-32 * SIZE(X), %xmm4
	movdqa	-28 * SIZE(X), %xmm5
	
	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm4, %xmm8
	movdqa	  %xmm4, %xmm9
	punpckldq %xmm8, %xmm4
	punpckhdq %xmm8, %xmm9
	
        paddq	%xmm4, %xmm0
        paddq   %xmm9, %xmm1

	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm5, %xmm8
	movdqa	  %xmm5, %xmm9
	punpckldq %xmm8, %xmm5
	punpckhdq %xmm8, %xmm9
	
        paddq	%xmm5, %xmm2
        paddq   %xmm9, %xmm3
	
	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L16:
	testq	$4, N
	jle	.L17

	movdqa	-32 * SIZE(X), %xmm4

	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm4, %xmm8
	movdqa	  %xmm4, %xmm9
	punpckldq %xmm8, %xmm4
	punpckhdq %xmm8, %xmm9
	
        paddq	%xmm4, %xmm0
        paddq   %xmm9, %xmm1
	
	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L17:
	testq	$2, N
	jle	.L18

        movq	-32 * SIZE(X), %xmm4

	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm4, %xmm8
	punpckldq %xmm8, %xmm4
	
	paddq	%xmm4, %xmm0
	  
	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

.L18:
	testq	$1, N
	jle	.L999

        movd	-32 * SIZE(X), %xmm4
	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm4, %xmm8
	punpckldq %xmm8, %xmm4
	
	paddq	%xmm4, %xmm0

	jmp	.L999
	ALIGN_3


.L50:
	movq	N,  %rax
	sarq	$2, %rax
	jle	.L55
	ALIGN_3

.L53:
        movd	0 * SIZE(X), %xmm4
        addq	INCX, X
        movd	0 * SIZE(X), %xmm5
        addq	INCX, X
	punpckldq %xmm5, %xmm4

	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm4, %xmm8
	punpckldq %xmm8, %xmm4
	
	paddq	%xmm4, %xmm0
	
        movd	0 * SIZE(X), %xmm6
        addq	INCX, X
        movd	0 * SIZE(X), %xmm7
        addq	INCX, X
	punpckldq %xmm7, %xmm6

	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm6, %xmm8
	punpckldq %xmm8, %xmm6
	
	paddq	%xmm6, %xmm1
          
	decq	%rax
	jg	.L53
	ALIGN_3

.L55:
	movq	N, %rax
	andq	$3,   %rax
	jle	.L999
	ALIGN_3

.L56:
        movd	0 * SIZE(X), %xmm4
	addq	INCX, X

	pxor      %xmm8,%xmm8
	pcmpgtd	  %xmm4, %xmm8
	punpckldq %xmm8, %xmm4
	
	paddq	%xmm4, %xmm0
	
	decq	%rax
	jg	.L56

	ALIGN_3

.L999:
	paddq	%xmm1, %xmm0
	paddq	%xmm3, %xmm2
	paddq	%xmm2, %xmm0
   
        movq    %xmm0,%rbx
        psrldq  $8, %xmm0
        movq    %xmm0,%rax
        addq    %rbx,%rax

	RESTOREREGISTERS
     
	ret
.LFE2:
	.size	isum_kh, .-isum_kh
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
