/*! \file shdot_k.S	
 * This file is part of the LEAC.
 *
 * Implementation of the 
 *	
 *	dot = x^Ty, (int <-- short ^T short)
 *
 *  function for integers
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
#define SIZE    2

#include "l1param.h"

/* NO SE NOS OLVIDE QUITAR ESTA INSTRUCCION
#undef  ALIGNED_ACCESS
*/
        .file	"shdot_k.S"
	.text
	.p2align 4,,15
.globl shdot_k
	.type	shdot_k, @function
shdot_k:
.LFB2:
	leaq	(, INCX, SIZE), INCX	
	leaq	(, INCY, SIZE), INCY	

	pxor	%xmm0, %xmm0
	pxor	%xmm1, %xmm1
	pxor	%xmm2, %xmm2
        pxor	%xmm3, %xmm3

	cmpq	$0, N
	jle	.L999

	cmpq	$SIZE, INCX
	jne	.L90
	cmpq	$SIZE, INCY
	jne	.L90

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y

	cmpq	$7, N
	jle	.L15

        testq	$SIZE, Y
	je	.L05

        movw    -32 * SIZE(X), %bx
        subw    -32 * SIZE(Y), %bx
        movswl  %bx, %eax
        imull   %eax, %eax
        movd    %eax, %xmm0
        
        addq    $1 * SIZE, X
        addq    $1 * SIZE, Y
        decq    N
        ALIGN_2

.L05:
	testq	$2 * SIZE, Y
	je	.L06

        movd	-32 * SIZE(X), %xmm1
        movd	-32 * SIZE(Y), %xmm4
        psubw    %xmm4, %xmm1
        pmaddwd  %xmm1, %xmm1
        
        addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
        subq	$2, N
	jle	.L999
	ALIGN_2

.L06:
	testq	$4 * SIZE, Y
	je	.L10

        movq    -32 * SIZE(X), %xmm2
        movq    -32 * SIZE(Y), %xmm4
        psubw    %xmm4, %xmm2
        pmaddwd  %xmm2, %xmm2
        
        addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
        subq	$4, N
	jle	.L999
	ALIGN_2

.L10:
#ifdef ALIGNED_ACCESS
        /* 8  A C E */
        testq	$4 * SIZE, X
	jne	.L50

        /* 4  6 */
	testq	$2 * SIZE, X
	jne	.L30

        /* 2 */
        testq	$1 * SIZE, X
	jne	.L20
  
#else
	testq	$7 * SIZE, X
	jne	.L20
#endif

	movq	N,  %rax
	sarq	$5, %rax
	jle	.L13

	movdqa	-32 * SIZE(X), %xmm4
	movdqa	-24 * SIZE(X), %xmm5
	movdqa	-16 * SIZE(X), %xmm6
	movdqa	 -8 * SIZE(X), %xmm7

	decq	%rax
	jle	.L12

	ALIGN_3

.L11:
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0
	movdqa	 0 * SIZE(X), %xmm4

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	psubw	 -24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5 
        paddd    %xmm5, %xmm1
	movdqa	 8 * SIZE(X), %xmm5

	psubw	-16 * SIZE(Y), %xmm6
        pmaddwd  %xmm6, %xmm6
        paddd    %xmm6, %xmm2
       	movdqa	 16 * SIZE(X), %xmm6

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  32) - PREOFFSET(X)
#endif

	psubw	-8 * SIZE(Y), %xmm7
        pmaddwd  %xmm7, %xmm7
        paddd    %xmm7,%xmm3
        movdqa	 24 * SIZE(X), %xmm7

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  32) - PREOFFSET(Y)
#endif

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y

	decq	%rax
	jg	.L11
	ALIGN_3

.L12:
	psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0

	psubw	-24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5 
        paddd    %xmm5, %xmm1

	psubw	-16 * SIZE(Y), %xmm6
        pmaddwd  %xmm6, %xmm6 
        paddd    %xmm6, %xmm2

        psubw	-8 * SIZE(Y), %xmm7
        pmaddwd  %xmm7, %xmm7 
        paddd    %xmm7, %xmm3
	
	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	ALIGN_3

.L13:
        testq	$31, N
	jle	.L999

	testq	$16, N
	jle	.L14
  
	movdqa	-32 * SIZE(X), %xmm4
	movdqa	-24 * SIZE(X), %xmm5

	psubw    -32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4 
        paddd    %xmm4, %xmm0
	
	psubw	 -24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5
        paddd    %xmm5, %xmm1
	
	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
     	ALIGN_3

.L14:
	testq	$8, N
	jle	.L15
  
	movdqa	-32 * SIZE(X), %xmm4
	psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0
        
	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y

	ALIGN_3

.L15:
        testq	$4, N
	jle	.L16        

        movq    -32 * SIZE(X), %xmm4
        movq    -32 * SIZE(Y), %xmm5
        psubw    %xmm5, %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0

        addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
        
	ALIGN_3

.L16:
        testq	$2, N
	jle	.L17

        movd	-32 * SIZE(X), %xmm4
        movd	-32 * SIZE(Y), %xmm5
        psubw    %xmm5, %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0
        
        addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
        
	ALIGN_3

.L17:
        testq	$1, N
	jle	.L999

        movw    -32 * SIZE(X), %bx
        subw    -32 * SIZE(Y), %bx
        movswl  %bx, %eax
        imull   %eax, %eax
        movd    %eax, %xmm4
        paddd    %xmm4, %xmm0
        jmp	.L999
        ALIGN_3
        
.L20:

#ifdef ALIGNED_ACCESS
  
        /* 2 */
        movdqa  -33 * SIZE(X), %xmm4
        psrldq  $2, %xmm4
       	addq	$7 * SIZE, X
      
        movq	N,  %rax
	sarq	$5, %rax
	jle	.L23

	movdqa	-32 * SIZE(X), %xmm5
	movdqa	-24 * SIZE(X), %xmm6
	movdqa	-16 * SIZE(X), %xmm7
	movdqa	 -8 * SIZE(X), %xmm8

	decq	%rax
	jle	.L22

	ALIGN_3

.L21:
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

        pshufd  $0xe4,%xmm5, %xmm9
        pslldq  $14, %xmm9
        por     %xmm9, %xmm4
	psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0

        psrldq  $2, %xmm5
        
        
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

        pshufd  $0xe4,%xmm6, %xmm10
        pslldq  $14, %xmm10
        por     %xmm10, %xmm5
       	psubw	 -24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5 
        paddd    %xmm5, %xmm1
	movdqa	 0 * SIZE(X), %xmm5

        psrldq   $2,  %xmm6

        pshufd  $0xe4,%xmm7, %xmm11
        pslldq  $14, %xmm11
        por     %xmm11, %xmm6
        psubw	-16 * SIZE(Y), %xmm6
        pmaddwd  %xmm6, %xmm6
        paddd    %xmm6, %xmm2
        movdqa	 8 * SIZE(X), %xmm6

        psrldq   $2,  %xmm7

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  32) - PREOFFSET(X)
#endif
        pshufd  $0xe4,%xmm8, %xmm12
        pslldq  $14, %xmm12
        por     %xmm12, %xmm7
        psubw	-8 * SIZE(Y), %xmm7
        pmaddwd  %xmm7, %xmm7
        paddd    %xmm7,%xmm3

        psrldq   $2, %xmm8
        pshufd   $0xe4, %xmm8, %xmm4
        
        movdqa	 16 * SIZE(X), %xmm7
        movdqa	 24 * SIZE(X), %xmm8

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  32) - PREOFFSET(Y)
#endif

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y

	decq	%rax
	jg	.L21
	ALIGN_3

.L22:
        pshufd  $0xe4,%xmm5, %xmm9
        pslldq  $14, %xmm9
        por     %xmm9, %xmm4
      	psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0

        psrldq   $2, %xmm5

        pshufd  $0xe4,%xmm6, %xmm10
        pslldq  $14, %xmm10
        por     %xmm10, %xmm5
	psubw	-24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5 
        paddd    %xmm5, %xmm1

        psrldq   $2,  %xmm6

        pshufd  $0xe4,%xmm7, %xmm11
        pslldq  $14, %xmm11
        por     %xmm11, %xmm6
	psubw	-16 * SIZE(Y), %xmm6
        pmaddwd  %xmm6, %xmm6 
        paddd    %xmm6, %xmm2

        psrldq   $2,  %xmm7

        pshufd  $0xe4,%xmm8, %xmm12
        pslldq  $14, %xmm12
        por     %xmm12, %xmm7
        psubw	-8 * SIZE(Y), %xmm7
        pmaddwd  %xmm7, %xmm7 
        paddd    %xmm7, %xmm3

        psrldq   $2, %xmm8
        pshufd   $0xe4, %xmm8, %xmm4
       
	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	ALIGN_3

.L23:
        testq	$31, N
	jle	.L999

	testq	$16, N
	jle	.L24
  
	movdqa	-32 * SIZE(X), %xmm5
	movdqa	-24 * SIZE(X), %xmm6
     
        pshufd  $0xe4,%xmm5, %xmm9
        pslldq  $14, %xmm9
        por     %xmm9, %xmm4
	psubw    -32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4 
        paddd    %xmm4, %xmm0

        psrldq  $2, %xmm5
        pshufd  $0xe4,%xmm6, %xmm10

        pslldq  $14, %xmm10
        por     %xmm10, %xmm5
	psubw	 -24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5
        paddd    %xmm5, %xmm1

        psrldq   $2,    %xmm6
        pshufd   $0xe4, %xmm6, %xmm4
        
    	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
     	ALIGN_3

.L24:
	testq	$8, N
	jle	.L25
        
        movdqa	-32 * SIZE(X), %xmm5

        pshufd  $0xe4,%xmm5, %xmm9
        pslldq  $14, %xmm9
        por     %xmm9, %xmm4     
       	psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0

        psrldq   $2, %xmm5
        pshufd  $0xe4,%xmm5, %xmm4
               
	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y

	ALIGN_3

.L25:
        testq	$4, N
	jle	.L26        
         
        movq     %xmm4,%xmm5
        movq     -32 * SIZE(Y), %xmm7
        psubw    %xmm7, %xmm5
        pmaddwd  %xmm5, %xmm5
        paddd    %xmm5, %xmm0

        psrldq   $8, %xmm4
        
        addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3
       
.L26:
        testq	$2, N
	jle	.L27

        pslldq  $4, %xmm4
        movq    %xmm4, %xmm5
        psrldq  $4, %xmm5
        
        movd	-32 * SIZE(Y), %xmm6
        psubw    %xmm6, %xmm5
        pmaddwd  %xmm5, %xmm5
        paddd    %xmm5, %xmm1

        psrldq  $8, %xmm4
        
        addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
        
	ALIGN_3

.L27:
        testq	$1, N
	jle	.L999

        pslldq  $6, %xmm4
        movq    %xmm4, %xmm5
        psrldq  $6, %xmm5

        movd	-33 * SIZE(Y), %xmm6
        psrldq  $2, %xmm6

        psubw    %xmm6, %xmm5
        pmaddwd  %xmm5, %xmm5
        paddd    %xmm5, %xmm3
        jmp	.L999
        ALIGN_3

.L30:
      /* 4  6 */
         / * 6 */   
        testq	$1 * SIZE, X
	jne	.L40

         /*  4  */   
        movdqa  -34 * SIZE(X), %xmm4
        psrldq  $4, %xmm4
       	addq	$6 * SIZE, X
      
        movq	N,  %rax
	sarq	$5, %rax
	jle	.L33

	movdqa	-32 * SIZE(X), %xmm5
	movdqa	-24 * SIZE(X), %xmm6
	movdqa	-16 * SIZE(X), %xmm7
	movdqa	 -8 * SIZE(X), %xmm8

	decq	%rax
	jle	.L32

	ALIGN_3

.L31:
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

        pshufd  $0xe4,%xmm5, %xmm9
        pslldq  $12, %xmm9
        por     %xmm9, %xmm4
	psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0

        psrldq  $4, %xmm5
        
        
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

        pshufd  $0xe4,%xmm6, %xmm10
        pslldq  $12, %xmm10
        por     %xmm10, %xmm5
       	psubw	 -24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5 
        paddd    %xmm5, %xmm1
	movdqa	 0 * SIZE(X), %xmm5

        psrldq   $4,  %xmm6

        pshufd  $0xe4,%xmm7, %xmm11
        pslldq  $12, %xmm11
        por     %xmm11, %xmm6
        psubw	-16 * SIZE(Y), %xmm6
        pmaddwd  %xmm6, %xmm6
        paddd    %xmm6, %xmm2
        movdqa	 8 * SIZE(X), %xmm6

        psrldq   $4,  %xmm7

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  32) - PREOFFSET(X)
#endif
        pshufd  $0xe4,%xmm8, %xmm12
        pslldq  $12, %xmm12
        por     %xmm12, %xmm7
        psubw	-8 * SIZE(Y), %xmm7
        pmaddwd  %xmm7, %xmm7
        paddd    %xmm7,%xmm3

        psrldq   $4, %xmm8
        pshufd   $0xe4, %xmm8, %xmm4
        
        movdqa	 16 * SIZE(X), %xmm7
        movdqa	 24 * SIZE(X), %xmm8

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  32) - PREOFFSET(Y)
#endif

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y

	decq	%rax
	jg	.L31
	ALIGN_3

.L32:
        pshufd  $0xe4,%xmm5, %xmm9
        pslldq  $12, %xmm9
        por     %xmm9, %xmm4
      	psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0

        psrldq   $4, %xmm5

        pshufd  $0xe4,%xmm6, %xmm10
        pslldq  $12, %xmm10
        por     %xmm10, %xmm5
	psubw	-24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5 
        paddd    %xmm5, %xmm1

        psrldq   $4,  %xmm6

        pshufd  $0xe4,%xmm7, %xmm11
        pslldq  $12, %xmm11
        por     %xmm11, %xmm6
	psubw	-16 * SIZE(Y), %xmm6
        pmaddwd  %xmm6, %xmm6 
        paddd    %xmm6, %xmm2

        psrldq   $4,  %xmm7

        pshufd  $0xe4,%xmm8, %xmm12
        pslldq  $12, %xmm12
        por     %xmm12, %xmm7
        psubw	-8 * SIZE(Y), %xmm7
        pmaddwd  %xmm7, %xmm7 
        paddd    %xmm7, %xmm3

        psrldq   $4, %xmm8
        pshufd   $0xe4, %xmm8, %xmm4
       
	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	ALIGN_3

.L33:
        testq	$31, N
	jle	.L999

	testq	$16, N
	jle	.L34
  
	movdqa	-32 * SIZE(X), %xmm5
	movdqa	-24 * SIZE(X), %xmm6
     
        pshufd  $0xe4,%xmm5, %xmm9
        pslldq  $12, %xmm9
        por     %xmm9, %xmm4
	psubw    -32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4 
        paddd    %xmm4, %xmm0

        psrldq  $4, %xmm5

        pshufd  $0xe4,%xmm6, %xmm10
        pslldq  $12, %xmm10
        por     %xmm10, %xmm5
	psubw	 -24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5
        paddd    %xmm5, %xmm1

        psrldq   $4,    %xmm6
        pshufd   $0xe4, %xmm6, %xmm4
        
    	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
     	ALIGN_3

.L34:
	testq	$8, N
	jle	.L35
        
        movdqa	-32 * SIZE(X), %xmm5

        pshufd  $0xe4,%xmm5, %xmm9
        pslldq  $12, %xmm9
        por     %xmm9, %xmm4     
       	psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0

        psrldq   $4, %xmm5
        pshufd  $0xe4,%xmm5, %xmm4
               
	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y

	ALIGN_3

.L35:
        testq	$4, N
	jle	.L36        

        movdqa  -32 * SIZE(X), %xmm5

        movq     %xmm4, %xmm6
        movq     -32 * SIZE(Y), %xmm7
        psubw    %xmm7, %xmm6
        pmaddwd  %xmm6, %xmm6
        paddd    %xmm6, %xmm0

        /*punpckhqdq %xmm6, %xmm4*/

        psrldq   $8, %xmm4
        pslldq   $4, %xmm5
        por     %xmm5, %xmm4     

        addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3
       
.L36:
        testq	$2, N
	jle	.L37

        pslldq  $4, %xmm4
        movq    %xmm4, %xmm5
        psrldq  $4, %xmm5
        
        movd	-32 * SIZE(Y), %xmm6
        psubw    %xmm6, %xmm5
        pmaddwd  %xmm5, %xmm5
        paddd    %xmm5, %xmm1

        psrldq  $8, %xmm4
        
        addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
        
	ALIGN_3

.L37:
        testq	$1, N
	jle	.L999

        pslldq  $6, %xmm4
        movq    %xmm4, %xmm5
        psrldq  $6, %xmm5

        movd	-33 * SIZE(Y), %xmm6
        psrldq  $2, %xmm6

        psubw    %xmm6, %xmm5
        pmaddwd  %xmm5, %xmm5
        paddd    %xmm5, %xmm3
        jmp	.L999
        ALIGN_3

.L40:
        / * 6 */   
        movdqa  -35 * SIZE(X), %xmm4
        psrldq  $6, %xmm4
       	addq	$5 * SIZE, X
      
        movq	N,  %rax
	sarq	$5, %rax
	jle	.L43

	movdqa	-32 * SIZE(X), %xmm5
	movdqa	-24 * SIZE(X), %xmm6
	movdqa	-16 * SIZE(X), %xmm7
	movdqa	 -8 * SIZE(X), %xmm8

	decq	%rax
	jle	.L42

	ALIGN_3

.L41:
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

        pshufd  $0xe4,%xmm5, %xmm9
        pslldq  $10, %xmm9
        por     %xmm9, %xmm4
	psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0

        psrldq  $6, %xmm5
        
        
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

        pshufd  $0xe4,%xmm6, %xmm10
        pslldq  $10, %xmm10
        por     %xmm10, %xmm5
       	psubw	 -24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5 
        paddd    %xmm5, %xmm1
	movdqa	 0 * SIZE(X), %xmm5

        psrldq   $6,  %xmm6

        pshufd  $0xe4,%xmm7, %xmm11
        pslldq  $10, %xmm11
        por     %xmm11, %xmm6
        psubw	-16 * SIZE(Y), %xmm6
        pmaddwd  %xmm6, %xmm6
        paddd    %xmm6, %xmm2
        movdqa	 8 * SIZE(X), %xmm6

        psrldq   $6,  %xmm7

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  32) - PREOFFSET(X)
#endif
        pshufd  $0xe4,%xmm8, %xmm12
        pslldq  $10, %xmm12
        por     %xmm12, %xmm7
        psubw	-8 * SIZE(Y), %xmm7
        pmaddwd  %xmm7, %xmm7
        paddd    %xmm7,%xmm3

        psrldq   $6, %xmm8
        pshufd   $0xe4, %xmm8, %xmm4
        
        movdqa	 16 * SIZE(X), %xmm7
        movdqa	 24 * SIZE(X), %xmm8

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  32) - PREOFFSET(Y)
#endif

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y

	decq	%rax
	jg	.L41
	ALIGN_3

.L42:
        pshufd  $0xe4,%xmm5, %xmm9
        pslldq  $10, %xmm9
        por     %xmm9, %xmm4
      	psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0

        psrldq   $6, %xmm5

        pshufd  $0xe4,%xmm6, %xmm10
        pslldq  $10, %xmm10
        por     %xmm10, %xmm5
	psubw	-24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5 
        paddd    %xmm5, %xmm1

        psrldq   $6,  %xmm6

        pshufd  $0xe4,%xmm7, %xmm11
        pslldq  $10, %xmm11
        por     %xmm11, %xmm6
	psubw	-16 * SIZE(Y), %xmm6
        pmaddwd  %xmm6, %xmm6 
        paddd    %xmm6, %xmm2

        psrldq   $6,  %xmm7

        pshufd  $0xe4,%xmm8, %xmm12
        pslldq  $10, %xmm12
        por     %xmm12, %xmm7
        psubw	-8 * SIZE(Y), %xmm7
        pmaddwd  %xmm7, %xmm7 
        paddd    %xmm7, %xmm3

        psrldq   $6, %xmm8
        pshufd   $0xe4, %xmm8, %xmm4
       
	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	ALIGN_3

.L43:
        testq	$31, N
	jle	.L999

	testq	$16, N
	jle	.L44
  
	movdqa	-32 * SIZE(X), %xmm5
	movdqa	-24 * SIZE(X), %xmm6
     
        pshufd  $0xe4,%xmm5, %xmm9
        pslldq  $10, %xmm9
        por     %xmm9, %xmm4
	psubw    -32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4 
        paddd    %xmm4, %xmm0

        psrldq  $6, %xmm5

        pshufd  $0xe4,%xmm6, %xmm10
        pslldq  $10, %xmm10
        por     %xmm10, %xmm5
	psubw	 -24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5
        paddd    %xmm5, %xmm1

        psrldq   $6,    %xmm6
        pshufd   $0xe4, %xmm6, %xmm4
        
    	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
     	ALIGN_3

.L44:
	testq	$8, N
	jle	.L45
        
        movdqa	-32 * SIZE(X), %xmm5

        pshufd  $0xe4,%xmm5, %xmm9
        pslldq  $10, %xmm9
        por     %xmm9, %xmm4     
       	psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0

        psrldq   $6, %xmm5
        pshufd  $0xe4,%xmm5, %xmm4
               
	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y

	ALIGN_3

.L45:
        testq	$4, N
	jle	.L46        

        movdqa  -32 * SIZE(X), %xmm5

        movq     %xmm4, %xmm6
        movq     -32 * SIZE(Y), %xmm7
        psubw    %xmm7, %xmm6
        pmaddwd  %xmm6, %xmm6
        paddd    %xmm6, %xmm0

        psrldq   $8, %xmm4
        pslldq   $2, %xmm5
        por     %xmm5, %xmm4     

        addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3
       
.L46:
        testq	$2, N
	jle	.L47

        pslldq  $4, %xmm4
        movq    %xmm4, %xmm5
        psrldq  $4, %xmm5
        
        movd	-32 * SIZE(Y), %xmm6
        psubw    %xmm6, %xmm5
        pmaddwd  %xmm5, %xmm5
        paddd    %xmm5, %xmm1

        psrldq  $8, %xmm4
        
        addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
        
	ALIGN_3

.L47:
        testq	$1, N
	jle	.L999

        pslldq  $6, %xmm4
        movq    %xmm4, %xmm5
        psrldq  $6, %xmm5

        movd	-33 * SIZE(Y), %xmm6
        psrldq  $2, %xmm6

        psubw    %xmm6, %xmm5
        pmaddwd  %xmm5, %xmm5
        paddd    %xmm5, %xmm3
        jmp	.L999
        ALIGN_3

.L50:
        /* 8  A C E */

         /* C E */
	testq	$2 * SIZE, X
	jne	.L70

        /* 8  A  */
  
           /*  A  */
        testq	$1 * SIZE, X
	jne	.L60

           /* 8  */
        movq    -32 * SIZE(X), %xmm4
	addq	$4 * SIZE, X
      
        movq	N,  %rax
	sarq	$5, %rax
	jle	.L53

	movdqa	-32 * SIZE(X), %xmm5
	movdqa	-24 * SIZE(X), %xmm6
	movdqa	-16 * SIZE(X), %xmm7
	movdqa	 -8 * SIZE(X), %xmm8

	decq	%rax
	jle	.L52

	ALIGN_3

.L51:
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

        punpcklqdq %xmm5, %xmm4
	psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0
        pshufd  $0x4e, %xmm5, %xmm5
        
	/*movdqa	 0 * SIZE(X), %xmm4*/

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

        punpcklqdq %xmm6, %xmm5
	psubw	 -24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5 
        paddd    %xmm5, %xmm1
        pshufd  $0x4e, %xmm6, %xmm6
	movdqa	 0 * SIZE(X), %xmm5

        punpcklqdq %xmm7, %xmm6
	psubw	-16 * SIZE(Y), %xmm6
        pmaddwd  %xmm6, %xmm6
        paddd    %xmm6, %xmm2
        pshufd  $0x4e, %xmm7, %xmm7
       	movdqa	 8 * SIZE(X), %xmm6

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  32) - PREOFFSET(X)
#endif

        punpcklqdq %xmm8, %xmm7
      	psubw	-8 * SIZE(Y), %xmm7
        pmaddwd  %xmm7, %xmm7
        paddd    %xmm7,%xmm3

        pshufd  $0x4e, %xmm8, %xmm4
        movdqa	 16 * SIZE(X), %xmm7
        movdqa	 24 * SIZE(X), %xmm8

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  32) - PREOFFSET(Y)
#endif

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y

	decq	%rax
	jg	.L51
	ALIGN_3

.L52:
        punpcklqdq %xmm5, %xmm4
	psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0
        pshufd  $0x4e, %xmm5, %xmm5
        
        punpcklqdq %xmm6, %xmm5
	psubw	-24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5 
        paddd    %xmm5, %xmm1
        pshufd  $0x4e, %xmm6, %xmm6

        punpcklqdq %xmm7, %xmm6
	psubw	-16 * SIZE(Y), %xmm6
        pmaddwd  %xmm6, %xmm6 
        paddd    %xmm6, %xmm2
        pshufd  $0x4e, %xmm7, %xmm7

        punpcklqdq %xmm8, %xmm7
        psubw	-8 * SIZE(Y), %xmm7
        pmaddwd  %xmm7, %xmm7 
        paddd    %xmm7, %xmm3
        pshufd  $0x4e, %xmm8, %xmm4
	
	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	ALIGN_3

.L53:
        testq	$31, N
	jle	.L999

	testq	$16, N
	jle	.L54
  
	movdqa	-32 * SIZE(X), %xmm5
	movdqa	-24 * SIZE(X), %xmm6

        punpcklqdq %xmm5, %xmm4 
	psubw    -32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4 
        paddd    %xmm4, %xmm0
        pshufd  $0x4e, %xmm5, %xmm5

        punpcklqdq %xmm6, %xmm5
	psubw	 -24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5
        paddd    %xmm5, %xmm1
        pshufd  $0x4e, %xmm6, %xmm4
	
	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
     	ALIGN_3

.L54:
	testq	$8, N
	jle	.L55

        movdqa	-32 * SIZE(X), %xmm5
        punpcklqdq %xmm5, %xmm4
       	psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0
        pshufd  $0x4e, %xmm5, %xmm4

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y

	ALIGN_3

.L55:
        testq	$4, N
	jle	.L56        

        movq      %xmm4, %xmm5
        movdqa     -32 * SIZE(Y), %xmm6
        movq      %xmm6, %xmm7
        
        psubw    %xmm7, %xmm5
        pmaddwd  %xmm5, %xmm5
        paddd    %xmm5, %xmm0
        movdqa    -32 * SIZE(X), %xmm4

        addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3
       
.L56:
        testq	$2, N
	jle	.L57

        movq      %xmm4, %xmm5
        pshufd  $0xa8, %xmm5, %xmm5        
        movd	-32 * SIZE(Y), %xmm6
        psubw    %xmm6, %xmm5
        pmaddwd  %xmm5, %xmm5
        paddd    %xmm5, %xmm1
        psrldq   $4, %xmm4
        
        addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
        
	ALIGN_3

.L57:
        testq	$1, N
	jle	.L999

        pslldq   $14, %xmm4
        psrldq   $14, %xmm4

        movd	-33 * SIZE(Y), %xmm5
        psrldq  $2, %xmm5

        psubw    %xmm5, %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm3
        jmp	.L999

        ALIGN_3

.L60:
          /*  A  */
        movq    -33 * SIZE(X), %xmm4
        psrldq  $2, %xmm4
	addq	$3 * SIZE, X
      
        movq	N,  %rax
	sarq	$5, %rax
	jle	.L63

	movdqa	-32 * SIZE(X), %xmm5
	movdqa	-24 * SIZE(X), %xmm6
	movdqa	-16 * SIZE(X), %xmm7
	movdqa	 -8 * SIZE(X), %xmm8

	decq	%rax
	jle	.L62

	ALIGN_3

.L61:
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

        pshufd  $0xe4,%xmm5, %xmm9
        pslldq  $6, %xmm9
        por     %xmm9, %xmm4
	psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0
        psrldq   $10, %xmm5

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif
        pshufd  $0xe4,%xmm6, %xmm10
        pslldq  $6, %xmm10
        por     %xmm10, %xmm5
       	psubw	 -24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5 
        paddd    %xmm5, %xmm1
        psrldq   $10,  %xmm6
	movdqa	 0 * SIZE(X), %xmm5

        pshufd  $0xe4,%xmm7, %xmm11
        pslldq  $6, %xmm11
        por     %xmm11, %xmm6
        psubw	-16 * SIZE(Y), %xmm6
        pmaddwd  %xmm6, %xmm6
        paddd    %xmm6, %xmm2
        psrldq   $10,  %xmm7
        movdqa	 8 * SIZE(X), %xmm6

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  32) - PREOFFSET(X)
#endif

        pshufd  $0xe4,%xmm8, %xmm12
        pslldq  $6, %xmm12
        por     %xmm12, %xmm7
        psubw	-8 * SIZE(Y), %xmm7
        pmaddwd  %xmm7, %xmm7
        paddd    %xmm7,%xmm3

        psrldq   $10,  %xmm8
        pshufd   $0xe4,%xmm8, %xmm4

        movdqa	 16 * SIZE(X), %xmm7
        movdqa	 24 * SIZE(X), %xmm8

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  32) - PREOFFSET(Y)
#endif

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y

	decq	%rax
	jg	.L61
	ALIGN_3

.L62:
        pshufd  $0xe4,%xmm5, %xmm9
        pslldq  $6, %xmm9
        por     %xmm9, %xmm4
      	psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0
        psrldq   $10, %xmm5
              
        pshufd  $0xe4,%xmm6, %xmm10
        pslldq  $6, %xmm10
        por     %xmm10, %xmm5
	psubw	-24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5 
        paddd    %xmm5, %xmm1
        psrldq   $10,  %xmm6

        pshufd  $0xe4,%xmm7, %xmm11
        pslldq  $6, %xmm11
        por     %xmm11, %xmm6
	psubw	-16 * SIZE(Y), %xmm6
        pmaddwd  %xmm6, %xmm6 
        paddd    %xmm6, %xmm2
        psrldq   $10,  %xmm7

        pshufd  $0xe4,%xmm8, %xmm12
        pslldq  $6, %xmm12
        por     %xmm12, %xmm7
        psubw	-8 * SIZE(Y), %xmm7
        pmaddwd  %xmm7, %xmm7 
        paddd    %xmm7, %xmm3

        psrldq   $10,  %xmm8
        pshufd   $0xe4,%xmm8, %xmm4

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	ALIGN_3

.L63:
        testq	$31, N
	jle	.L999

	testq	$16, N
	jle	.L64
  
	movdqa	-32 * SIZE(X), %xmm5
	movdqa	-24 * SIZE(X), %xmm6
     
        pshufd  $0xe4,%xmm5, %xmm9
        pslldq  $6, %xmm9
        por     %xmm9, %xmm4
 
	psubw    -32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4 
        paddd    %xmm4, %xmm0
        psrldq   $10, %xmm5
        
        pshufd  $0xe4,%xmm6, %xmm10
        pslldq  $6, %xmm10
        por     %xmm10, %xmm5

	psubw	 -24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5
        paddd    %xmm5, %xmm1

        psrldq   $10,  %xmm6
        pshufd   $0xe4,%xmm6, %xmm4
        
    	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
     	ALIGN_3

.L64:
	testq	$8, N
	jle	.L65
        
        movdqa	-32 * SIZE(X), %xmm5
        pshufd  $0xe4,%xmm5, %xmm9
        pslldq  $6, %xmm9
        por     %xmm9, %xmm4
     
       	psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0
        psrldq   $10, %xmm5
        pshufd   $0xe4, %xmm5, %xmm4
        

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y

	ALIGN_3

.L65:
        testq	$4, N
	jle	.L66        

        movdqa	-32 * SIZE(X), %xmm5
        pslldq  $6, %xmm5
        movq	%xmm5, %xmm6
        
        por     %xmm6, %xmm4
   
        movdqa  -32 * SIZE(Y), %xmm7
        movq    %xmm7, %xmm8

        psubw    %xmm8, %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0

        psrldq  $2, %xmm5
        pslldq  $2, %xmm5
        pshufd  $0xe, %xmm5, %xmm4
        

        addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3
       
.L66:
        testq	$2, N
	jle	.L67

        pslldq  $4, %xmm4
        movq    %xmm4, %xmm5
        psrldq  $4, %xmm5

        movd	-32 * SIZE(Y), %xmm6
        psubw    %xmm6, %xmm5
        pmaddwd  %xmm5, %xmm5
        paddd    %xmm5, %xmm1
        pshufd  $0x2, %xmm4, %xmm4
        
        addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
        
	ALIGN_3

.L67:
        testq	$1, N
	jle	.L999

        pslldq  $6, %xmm4
        movq    %xmm4, %xmm5
        psrldq  $6, %xmm5

        movd	-33 * SIZE(Y), %xmm6
        psrldq  $2, %xmm6

        psubw    %xmm6, %xmm5
        pmaddwd  %xmm5, %xmm5
        paddd    %xmm5, %xmm3
        jmp	.L999
        ALIGN_3
.L70:
        /* C E */
  
           /*  E  */
        testq	$1 * SIZE, X
	jne	.L80

          /*  C  */
        movd    -32 * SIZE(X), %xmm4
       	addq	$2 * SIZE, X
      
        movq	N,  %rax
	sarq	$5, %rax
	jle	.L73

	movdqa	-32 * SIZE(X), %xmm5
	movdqa	-24 * SIZE(X), %xmm6
	movdqa	-16 * SIZE(X), %xmm7
	movdqa	 -8 * SIZE(X), %xmm8

	decq	%rax
	jle	.L72

	ALIGN_3

.L71:
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

        pshufd  $0xe4,%xmm5, %xmm9
        pslldq  $4, %xmm9
        por     %xmm9, %xmm4
	psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0
        psrldq   $12, %xmm5

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif
        pshufd  $0xe4,%xmm6, %xmm10
        pslldq  $4, %xmm10
        por     %xmm10, %xmm5
       	psubw	 -24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5 
        paddd    %xmm5, %xmm1
        psrldq   $12,  %xmm6
	movdqa	 0 * SIZE(X), %xmm5

        pshufd  $0xe4,%xmm7, %xmm11
        pslldq  $4, %xmm11
        por     %xmm11, %xmm6
        psubw	-16 * SIZE(Y), %xmm6
        pmaddwd  %xmm6, %xmm6
        paddd    %xmm6, %xmm2
        psrldq   $12,  %xmm7
        movdqa	 8 * SIZE(X), %xmm6

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  32) - PREOFFSET(X)
#endif

        pshufd  $0xe4,%xmm8, %xmm12
        pslldq  $4, %xmm12
        por     %xmm12, %xmm7
        psubw	-8 * SIZE(Y), %xmm7
        pmaddwd  %xmm7, %xmm7
        paddd    %xmm7,%xmm3

        psrldq   $12,  %xmm8
        pshufd   $0xe4,%xmm8, %xmm4

        movdqa	 16 * SIZE(X), %xmm7
        movdqa	 24 * SIZE(X), %xmm8

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  32) - PREOFFSET(Y)
#endif

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y

	decq	%rax
	jg	.L71
	ALIGN_3

.L72:
        pshufd  $0xe4,%xmm5, %xmm9
        pslldq  $4, %xmm9
        por     %xmm9, %xmm4
      	psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0
        psrldq   $12, %xmm5
              
        pshufd  $0xe4,%xmm6, %xmm10
        pslldq  $4, %xmm10
        por     %xmm10, %xmm5
	psubw	-24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5 
        paddd    %xmm5, %xmm1
        psrldq   $12,  %xmm6

        pshufd  $0xe4,%xmm7, %xmm11
        pslldq  $4, %xmm11
        por     %xmm11, %xmm6
	psubw	-16 * SIZE(Y), %xmm6
        pmaddwd  %xmm6, %xmm6 
        paddd    %xmm6, %xmm2
        psrldq   $12,  %xmm7

        pshufd  $0xe4,%xmm8, %xmm12
        pslldq  $4, %xmm12
        por     %xmm12, %xmm7
        psubw	-8 * SIZE(Y), %xmm7
        pmaddwd  %xmm7, %xmm7 
        paddd    %xmm7, %xmm3

        psrldq   $12,  %xmm8
        pshufd   $0xe4,%xmm8, %xmm4

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	ALIGN_3

.L73:
        testq	$31, N
	jle	.L999

	testq	$16, N
	jle	.L74
  
	movdqa	-32 * SIZE(X), %xmm5
	movdqa	-24 * SIZE(X), %xmm6
     
        pshufd  $0xe4,%xmm5, %xmm9
        pslldq  $4, %xmm9
        por     %xmm9, %xmm4
	psubw    -32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4 
        paddd    %xmm4, %xmm0
        psrldq   $12, %xmm5
  
        pshufd  $0xe4,%xmm6, %xmm10
        pslldq  $4, %xmm10
        por     %xmm10, %xmm5
	psubw	 -24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5
        paddd    %xmm5, %xmm1

        psrldq   $12,  %xmm6
        pshufd   $0xe4,%xmm6, %xmm4
        
    	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
     	ALIGN_3

.L74:
	testq	$8, N
	jle	.L75
        
        movdqa	-32 * SIZE(X), %xmm5
        pshufd  $0xe4,%xmm5, %xmm9
        pslldq  $4, %xmm9
        por     %xmm9, %xmm4
     
       	psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0

        psrldq   $12, %xmm5
        pshufd   $0xe4, %xmm5, %xmm4
        
	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y

	ALIGN_3

.L75:
        testq	$4, N
	jle	.L76        

        movdqa	-32 * SIZE(X), %xmm5
        pslldq  $4, %xmm5
        movq	%xmm5, %xmm6        
        por     %xmm6, %xmm4
   
        movdqa  -32 * SIZE(Y), %xmm7
        movq    %xmm7, %xmm8

        psubw    %xmm8, %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0

        pshufd  $0x2, %xmm5, %xmm4
        
        addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3
       
.L76:
        testq	$2, N
	jle	.L77

        /*pslldq  $4, %xmm4
        movq    %xmm4, %xmm5
        psrldq  $4, %xmm5*/

        movd	-32 * SIZE(Y), %xmm6
        psubw    %xmm6, %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm1

        movd	-33 * SIZE(X), %xmm4
        psrldq  $2, %xmm4

        /*pshufd  $0x2, %xmm4, %xmm4*/
        
        addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
        
	ALIGN_3

.L77:
        testq	$1, N
	jle	.L999

        pslldq  $6, %xmm4
        movq    %xmm4, %xmm5
        psrldq  $6, %xmm5

        movd	-33 * SIZE(Y), %xmm6
        psrldq  $2, %xmm6

        psubw    %xmm6, %xmm5
        pmaddwd  %xmm5, %xmm5
        paddd    %xmm5, %xmm3
        jmp	.L999
        ALIGN_3

.L80:
        /*  E  */
        movd    -33 * SIZE(X), %xmm4
        psrldq  $2, %xmm4
       	addq	$1 * SIZE, X
      
        movq	N,  %rax
	sarq	$5, %rax
	jle	.L83

	movdqa	-32 * SIZE(X), %xmm5
	movdqa	-24 * SIZE(X), %xmm6
	movdqa	-16 * SIZE(X), %xmm7
	movdqa	 -8 * SIZE(X), %xmm8

	decq	%rax
	jle	.L82

	ALIGN_3

.L81:
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

        pshufd  $0xe4,%xmm5, %xmm9
        pslldq  $2, %xmm9
        por     %xmm9, %xmm4
	psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0

        psrldq   $14, %xmm5

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif
        pshufd  $0xe4,%xmm6, %xmm10
        pslldq  $2, %xmm10
        por     %xmm10, %xmm5
       	psubw	 -24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5 
        paddd    %xmm5, %xmm1
        psrldq   $14,  %xmm6
	movdqa	 0 * SIZE(X), %xmm5

        pshufd  $0xe4,%xmm7, %xmm11
        pslldq  $2, %xmm11
        por     %xmm11, %xmm6
        psubw	-16 * SIZE(Y), %xmm6
        pmaddwd  %xmm6, %xmm6
        paddd    %xmm6, %xmm2
        psrldq   $14,  %xmm7
        movdqa	 8 * SIZE(X), %xmm6

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  32) - PREOFFSET(X)
#endif

        pshufd  $0xe4,%xmm8, %xmm12
        pslldq  $2, %xmm12
        por     %xmm12, %xmm7
        psubw	-8 * SIZE(Y), %xmm7
        pmaddwd  %xmm7, %xmm7
        paddd    %xmm7,%xmm3

        psrldq   $6, %xmm8
        pshufd   $0xfe, %xmm8, %xmm4
        
        movdqa	 16 * SIZE(X), %xmm7
        movdqa	 24 * SIZE(X), %xmm8

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  32) - PREOFFSET(Y)
#endif

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y

	decq	%rax
	jg	.L81
	ALIGN_3

.L82:
        pshufd  $0xe4,%xmm5, %xmm9
        pslldq  $2, %xmm9
        por     %xmm9, %xmm4
      	psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0

        psrldq   $14, %xmm5
              
        pshufd  $0xe4,%xmm6, %xmm10
        pslldq  $2, %xmm10
        por     %xmm10, %xmm5
	psubw	-24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5 
        paddd    %xmm5, %xmm1

        psrldq   $14,  %xmm6

        pshufd  $0xe4,%xmm7, %xmm11
        pslldq  $2, %xmm11
        por     %xmm11, %xmm6
	psubw	-16 * SIZE(Y), %xmm6
        pmaddwd  %xmm6, %xmm6 
        paddd    %xmm6, %xmm2

        psrldq   $14,  %xmm7

        pshufd  $0xe4,%xmm8, %xmm12
        pslldq  $2, %xmm12
        por     %xmm12, %xmm7
        psubw	-8 * SIZE(Y), %xmm7
        pmaddwd  %xmm7, %xmm7 
        paddd    %xmm7, %xmm3

        psrldq   $6, %xmm8
        pshufd   $0xfe, %xmm8, %xmm4
       
	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	ALIGN_3

.L83:
        testq	$31, N
	jle	.L999

	testq	$16, N
	jle	.L84
  
	movdqa	-32 * SIZE(X), %xmm5
	movdqa	-24 * SIZE(X), %xmm6
     
        pshufd  $0xe4,%xmm5, %xmm9
        pslldq  $2, %xmm9
        por     %xmm9, %xmm4
	psubw    -32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4 
        paddd    %xmm4, %xmm0

        psrldq   $14, %xmm5
  
        pshufd  $0xe4,%xmm6, %xmm10
        pslldq  $2, %xmm10
        por     %xmm10, %xmm5
	psubw	 -24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5
        paddd    %xmm5, %xmm1

        psrldq   $6, %xmm6
        pshufd   $0xfe, %xmm6, %xmm4
        
    	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
     	ALIGN_3

.L84:
	testq	$8, N
	jle	.L85
        
        movdqa	-32 * SIZE(X), %xmm5
        pshufd  $0xe4,%xmm5, %xmm9
        pslldq  $2, %xmm9
        por     %xmm9, %xmm4
     
       	psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0

        psrldq   $6, %xmm5
        pshufd   $0xfe, %xmm5, %xmm4
               
	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y

	ALIGN_3

.L85:
        testq	$4, N
	jle	.L86        

        /*movdqa	-32 * SIZE(X), %xmm5*/
        movq	-32 * SIZE(X), %xmm5
        pslldq  $2, %xmm5
        movq	%xmm5, %xmm6        
        por     %xmm6, %xmm4

        movq   -32 * SIZE(Y), %xmm7
        /*movdqa  -32 * SIZE(Y), %xmm7
        movq    %xmm7, %xmm8*/

        psubw    %xmm7, %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0

        pshufd  $0xfe, %xmm5, %xmm4
        
        addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3
       
.L86:
        testq	$2, N
	jle	.L87

        movd    -32 * SIZE(X), %xmm5
        pslldq  $6, %xmm5
        movq    %xmm5, %xmm6
        psrldq  $4, %xmm6
        por     %xmm6, %xmm4

        movd	-32 * SIZE(Y), %xmm7
        psubw    %xmm7, %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm1

        pshufd  $0x2, %xmm5, %xmm4
        
        addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
        
	ALIGN_3

.L87:
        testq	$1, N
	jle	.L999

        pslldq  $6, %xmm4
        movq    %xmm4, %xmm5
        psrldq  $6, %xmm5

        movd	-33 * SIZE(Y), %xmm6
        psrldq  $2, %xmm6

        psubw    %xmm6, %xmm5
        pmaddwd  %xmm5, %xmm5
        paddd    %xmm5, %xmm3
        jmp	.L999
        ALIGN_3

#else /* L20 NO ALIGNED_ACCESS*/

       	movq	N,  %rax
	sarq	$5, %rax
	jle	.L24

	movq	-32 * SIZE(X), %xmm4
	movq	-28 * SIZE(X), %xmm8
        punpcklqdq %xmm8, %xmm4
	movq	-24 * SIZE(X), %xmm5
	movq	-20 * SIZE(X), %xmm9
        punpcklqdq %xmm9, %xmm5
	movq	-16 * SIZE(X), %xmm6
	movq	-12 * SIZE(X), %xmm10
        punpcklqdq %xmm10, %xmm6
	movq	-8 * SIZE(X), %xmm7
	movq	-4 * SIZE(X), %xmm11
        punpcklqdq %xmm11, %xmm7

	decq	%rax
	jle	.L22

	ALIGN_3

.L21:

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

        psubw	-32 * SIZE(Y), %xmm4
        pmaddwd %xmm4, %xmm4
        paddd   %xmm4, %xmm0
        movq	0 * SIZE(X), %xmm4
	movq	4 * SIZE(X), %xmm8
        punpcklqdq %xmm8, %xmm4

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

        psubw	 -24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5 
        paddd    %xmm5, %xmm1
        movq	8  * SIZE(X), %xmm5
	movq	12 * SIZE(X), %xmm9
        punpcklqdq %xmm9, %xmm5


        psubw	-16 * SIZE(Y), %xmm6
        pmaddwd  %xmm6, %xmm6
        paddd    %xmm6, %xmm2
        movq	 16 * SIZE(X), %xmm6
	movq	 20 * SIZE(X), %xmm10
        punpcklqdq %xmm10, %xmm6

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  32) - PREOFFSET(X)
#endif

        psubw	-8 * SIZE(Y), %xmm7
        pmaddwd  %xmm7, %xmm7
        paddd    %xmm7,%xmm3
        movq	 24 * SIZE(X), %xmm7
	movq	 28 * SIZE(X), %xmm11
        punpcklqdq %xmm11, %xmm7
        
#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  32) - PREOFFSET(Y)
#endif
       
	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y

	decq	%rax
	jg	.L21
	ALIGN_3


.L22:
        psubw	-32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0

	psubw	-24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5 
        paddd    %xmm5, %xmm1

	psubw	-16 * SIZE(Y), %xmm6
        pmaddwd  %xmm6, %xmm6 
        paddd    %xmm6, %xmm2

        psubw	-8 * SIZE(Y), %xmm7
        pmaddwd  %xmm7, %xmm7 
        paddd    %xmm7, %xmm3
	
	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	ALIGN_3

.L24:
        testq	$31, N
	jle	.L999

	testq	$16, N
	jle	.L25
 
        movq	-32 * SIZE(X), %xmm4
	movq	-28 * SIZE(X), %xmm8
        punpcklqdq %xmm8, %xmm4
	psubw    -32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4 
        paddd    %xmm4, %xmm0

        movq	-24  * SIZE(X), %xmm5
	movq	-20 * SIZE(X), %xmm9
        punpcklqdq %xmm9, %xmm5
	psubw	 -24 * SIZE(Y), %xmm5
        pmaddwd  %xmm5, %xmm5
        paddd    %xmm5, %xmm1
	
	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y

     	ALIGN_3	

.L25:
	testq	$8, N
	jle	.L26

	movq	-32 * SIZE(X), %xmm4
	movq	-28 * SIZE(X), %xmm8
        punpcklqdq %xmm8, %xmm4

        psubw    -32 * SIZE(Y), %xmm4
        pmaddwd  %xmm4, %xmm4 
        paddd    %xmm4, %xmm0
        
	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L26:
	testq	$4, N
	jle	.L27

        movq    -32 * SIZE(X), %xmm4
        movq    -32 * SIZE(Y), %xmm5
        psubw    %xmm5, %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0
        
	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L27:
	testq	$2, N
	jle	.L28

        movd	-32 * SIZE(X), %xmm4
        movd	-32 * SIZE(Y), %xmm5
        psubw    %xmm5, %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0

	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y

	ALIGN_3

.L28:
	testq	$1, N
	jle	.L999

        movw    -32 * SIZE(X), %bx
        subw    -32 * SIZE(Y), %bx
        movswl  %bx, %eax
        imull   %eax, %eax
        movd    %eax, %xmm4
        paddd    %xmm4, %xmm0
        jmp	.L999
        ALIGN_3
               
#endif

.L90:
	movq	N,  %rax
	sarq	$3, %rax
	jle	.L94
	ALIGN_3

.L93:
        pinsrw  $0, (X), %xmm4
        addq	INCX, X
        pinsrw  $1, (X), %xmm4
        addq	INCX, X
        pinsrw  $2, (X), %xmm4
        addq	INCX, X
        pinsrw  $3, (X), %xmm4
        addq	INCX, X
        pinsrw  $4, (X), %xmm4
        addq	INCX, X
        pinsrw  $5, (X), %xmm4
        addq	INCX, X
        pinsrw  $6, (X), %xmm4
        addq	INCX, X
        pinsrw  $7, (X), %xmm4
        addq	INCX, X

        pinsrw  $0, (Y), %xmm5
        addq	INCY, Y
        pinsrw  $1, (Y), %xmm5
        addq	INCY, Y
        pinsrw  $2, (Y), %xmm5
        addq	INCY, Y
        pinsrw  $3, (Y), %xmm5
        addq	INCY, Y
        pinsrw  $4, (Y), %xmm5
        addq	INCY, Y
        pinsrw  $5, (Y), %xmm5
        addq	INCY, Y
        pinsrw  $6, (Y), %xmm5
        addq	INCY, Y
        pinsrw  $7, (Y), %xmm5
        addq	INCY, Y

        psubw    %xmm5, %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm0
        
	decq	%rax
	jg	.L93

	ALIGN_3

.L94:        
        testq	$4, N
	jle	.L95

        pxor	%xmm4, %xmm4
        pxor	%xmm5, %xmm5

        pinsrw  $0, (X), %xmm4
        addq	INCX, X
        pinsrw  $1, (X), %xmm4
        addq	INCX, X
        pinsrw  $2, (X), %xmm4
        addq	INCX, X
        pinsrw  $3, (X), %xmm4
        addq	INCX, X
        
        pinsrw  $0, (Y), %xmm5
        addq	INCY, Y
        pinsrw  $1, (Y), %xmm5
        addq	INCY, Y
        pinsrw  $2, (Y), %xmm5
        addq	INCY, Y
        pinsrw  $3, (Y), %xmm5
        addq	INCY, Y

        psubw    %xmm5, %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm3

      	ALIGN_3


.L95:
        testq	$2, N
	jle	.L96

        pxor	%xmm4, %xmm4
        pxor	%xmm5, %xmm5

        pinsrw  $4, (X), %xmm4
        addq	INCX, X
        pinsrw  $5, (X), %xmm4
        addq	INCX, X
        
        pinsrw  $4, (Y), %xmm5
        addq	INCY, Y
        pinsrw  $5, (Y), %xmm5
        addq	INCY, Y
        psubw    %xmm5, %xmm4
        pmaddwd  %xmm4, %xmm4
        paddd    %xmm4, %xmm1

      	ALIGN_3

.L96:
	testq	$1, N
	jle	.L999

        movw    0 * SIZE(X), %bx
        subw    0 * SIZE(Y), %bx
        movswl  %bx, %eax
        imull   %eax, %eax
        movd    %eax, %xmm4
        paddd   %xmm4, %xmm0       
        ALIGN_3
.L999:
        
        /*pxor     %xmm5, %xmm5
        pxor     %xmm6, %xmm6
        pshufd    $0xe4,%xmm0, %xmm4
        punpcklwd %xmm5, %xmm0
        punpckhwd %xmm5, %xmm4

        paddd    %xmm4, %xmm0

        punpcklwd %xmm2,%xmm5
        punpckhwd %xmm2,%xmm6

        paddd    %xmm6, %xmm5

        paddd    %xmm5, %xmm0*/
          
	paddd	%xmm1, %xmm0
	paddd	%xmm3, %xmm2

        paddd	%xmm2, %xmm0

        phaddd  %xmm0, %xmm0
        phaddd  %xmm0, %xmm0

        movd    %xmm0, %eax
        
        /*movq    %xmm0,%rbx
        psrldq  $8, %xmm0
        movq    %xmm0,%rax
        addq    %rbx,%rax*/
     
	ret
.LFE2:
	.size	shdot_k, .-shdot_k
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
