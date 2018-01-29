/*! \file idot_kh.S	
 * This file is part of the LEAC.
 *
 * Implementation of the 
 *	
 *	dot = x^Ty, (long int <-- int ^T int)
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

#include "l1param.h"

/* NO SE NOS OLVIDE QUITAR ESTA INSTRUCCION
#undef  ALIGNED_ACCESS
*/  
        .file	"idot_kh.S"
	.text
	.p2align 4,,15
.globl idot_kh
	.type	idot_kh, @function
idot_kh:
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
	jne	.L50
	cmpq	$SIZE, INCY
	jne	.L50

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y

	cmpq	$3, N
	jle	.L17

        testq	$SIZE, Y
	je	.L05
        
        movd	-32 * SIZE(X), %xmm4
        movd	-32 * SIZE(Y), %xmm5

	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm5, %xmm4
#else
	
#endif
	paddq	%xmm4, %xmm0

     	addq	$1 * SIZE, X
	addq	$1 * SIZE, Y
	decq	N
	ALIGN_2

.L05:
	testq	$2 * SIZE, Y
	je	.L10

	movq	-32 * SIZE(X), %xmm6
        movq	-32 * SIZE(Y), %xmm7
	pshufd  $0x50,%xmm6,%xmm6
	pshufd  $0x50,%xmm7,%xmm7
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm7, %xmm6
#else
	
#endif
        paddq   %xmm6,%xmm0

	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	subq	$2, N
	jle	.L999
	ALIGN_2

.L10:
#ifdef ALIGNED_ACCESS	
	testq	$2 * SIZE, X
	jne	.L30

	testq	$SIZE, X
	jne	.L20
#else
	testq	$3 * SIZE, X
	jne	.L20
#endif

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
	
	movdqa	-32 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm4,%xmm13
	pshufd  $0x50,%xmm4,%xmm4
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm4
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm4,  %xmm0
        paddq   %xmm13, %xmm1
	movdqa	-16 * SIZE(X), %xmm4
		
	movdqa	-28 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm5,%xmm13
	pshufd  $0x50,%xmm5,%xmm5
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm5
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm5,  %xmm2
        paddq   %xmm13, %xmm3
	movdqa	-12 * SIZE(X), %xmm5

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	movdqa	-24 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm6,%xmm13
	pshufd  $0x50,%xmm6,%xmm6
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm6
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm6,  %xmm0
        paddq   %xmm13, %xmm1
	movdqa	-8 * SIZE(X), %xmm6
	
	movdqa	-20 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm7,%xmm13
	pshufd  $0x50,%xmm7,%xmm7
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm7
	pmuldq   %xmm14, %xmm13
#else
	
#endif	
        paddq	%xmm7,  %xmm2
        paddq   %xmm13, %xmm3
	movdqa	-4 * SIZE(X), %xmm7
	
#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(X)
#endif
	
	movdqa	-16 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm4,%xmm13
	pshufd  $0x50,%xmm4,%xmm4
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12

#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm4
	pmuldq   %xmm14, %xmm13
#else
	
#endif	
        paddq	%xmm4,  %xmm0
        paddq   %xmm13, %xmm1
	movdqa	 0 * SIZE(X), %xmm4
	
	movdqa	-12 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm5,%xmm13
	pshufd  $0x50,%xmm5,%xmm5
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12

#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm5
	pmuldq   %xmm14, %xmm13
#else
	
#endif	
        paddq	%xmm5,  %xmm2
        paddq   %xmm13, %xmm3
	movdqa	4 * SIZE(X), %xmm5

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(Y)
#endif

	movdqa	-8 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm6,%xmm13
	pshufd  $0x50,%xmm6,%xmm6
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12

#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm6
	pmuldq   %xmm14, %xmm13
#else
	
#endif	
        paddq	%xmm6,  %xmm0
        paddq   %xmm13,  %xmm1
    	movdqa	 8 * SIZE(X), %xmm6

	movdqa	-4 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm7,%xmm13
	pshufd  $0x50,%xmm7,%xmm7
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12

#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm7
	pmuldq   %xmm14, %xmm13
#else
	
#endif	
        paddq	%xmm7,  %xmm2
        paddq   %xmm13,  %xmm3
        movdqa	12 * SIZE(X), %xmm7
	
	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y

	decq	%rax
	jg	.L11
	ALIGN_3

.L12:
	movdqa	-32 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm4,%xmm13
	pshufd  $0x50,%xmm4,%xmm4
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12

#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm4
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm4,  %xmm0
        paddq   %xmm13, %xmm1
	
	movdqa	-16 * SIZE(X), %xmm4
	movdqa	-28 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm5,%xmm13
	pshufd  $0x50,%xmm5,%xmm5
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12

#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm5
	pmuldq   %xmm14, %xmm13
#else
	
#endif	
        paddq	%xmm5,  %xmm2
        paddq   %xmm13, %xmm3
	movdqa	-12 * SIZE(X), %xmm5

	movdqa	-24 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm6,%xmm13
	pshufd  $0x50,%xmm6,%xmm6
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12

#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm6
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm6,  %xmm0
        paddq   %xmm13, %xmm1
	movdqa	 -8 * SIZE(X), %xmm6

	movdqa	-20 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm7,%xmm13
	pshufd  $0x50,%xmm7,%xmm7
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12

#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm7
	pmuldq   %xmm14, %xmm13
#else
	
#endif	
        paddq	%xmm7,  %xmm2
        paddq   %xmm13, %xmm3
	
	movdqa	 -4 * SIZE(X), %xmm7
	movdqa	-16 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm4,%xmm13
	pshufd  $0x50,%xmm4,%xmm4
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12

#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm4
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm4,  %xmm0
        paddq   %xmm13, %xmm1

	movdqa	-12 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm5,%xmm13
	pshufd  $0x50,%xmm5,%xmm5
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12

#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm5
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm5,  %xmm2
        paddq   %xmm13, %xmm3

	movdqa	 -8 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm6,%xmm13
	pshufd  $0x50,%xmm6,%xmm6
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12

#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm6
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm6,  %xmm0
        paddq   %xmm13, %xmm1

	movdqa	 -4 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm7,%xmm13
	pshufd  $0x50,%xmm7,%xmm7
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12

#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm7
	pmuldq   %xmm14, %xmm13
#else
	
#endif	
        paddq	%xmm7,  %xmm2
        paddq   %xmm13, %xmm3
	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y

	ALIGN_3

.L14:
	testq	$31, N
	jle	.L999

	testq	$16, N
	jle	.L15

	movdqa	-32 * SIZE(X), %xmm4
	movdqa	-28 * SIZE(X), %xmm6
	movdqa	-32 * SIZE(Y), %xmm5
	movdqa	-28 * SIZE(Y), %xmm7
		
	pshufd  $0xfa,%xmm4,%xmm8
	pshufd  $0x50,%xmm4,%xmm4

	pshufd  $0xfa,%xmm6,%xmm10
	pshufd  $0x50,%xmm6,%xmm6

	pshufd  $0xfa,%xmm5,%xmm9
	pshufd  $0x50,%xmm5,%xmm5

	pshufd  $0xfa,%xmm7,%xmm11
	pshufd  $0x50,%xmm7,%xmm7
	
#ifdef HAVE_SSE4_1
	
	pmuldq   %xmm5, %xmm4
	pmuldq   %xmm7, %xmm6
	pmuldq   %xmm9, %xmm8
	pmuldq   %xmm11, %xmm10
#else
	
#endif
	paddq	%xmm4,%xmm0
	paddq	%xmm6,%xmm1
	paddq	%xmm8,%xmm2
	paddq	%xmm10,%xmm3

	
	movdqa	-24 * SIZE(X), %xmm4
	movdqa	-20 * SIZE(X), %xmm6
	movdqa	-24 * SIZE(Y), %xmm5
	movdqa	-20 * SIZE(Y), %xmm7
	
	pshufd  $0xfa,%xmm4,%xmm8
	pshufd  $0x50,%xmm4,%xmm4

	pshufd  $0xfa,%xmm6,%xmm10
	pshufd  $0x50,%xmm6,%xmm6

	pshufd  $0xfa,%xmm5,%xmm9
	pshufd  $0x50,%xmm5,%xmm5

	pshufd  $0xfa,%xmm7,%xmm11
	pshufd  $0x50,%xmm7,%xmm7
	
#ifdef HAVE_SSE4_1
	
	pmuldq   %xmm5, %xmm4
	pmuldq   %xmm7, %xmm6
	pmuldq   %xmm9, %xmm8
	pmuldq   %xmm11, %xmm10
#else
	
#endif
	paddq	%xmm4,%xmm0
	paddq	%xmm6,%xmm1
	paddq	%xmm8,%xmm2
	paddq	%xmm10,%xmm3
        
	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
	ALIGN_3

.L15:
	testq	$8, N
	jle	.L16

	movdqa	-32 * SIZE(X), %xmm4
	movdqa	-28 * SIZE(X), %xmm6
	movdqa	-32 * SIZE(Y), %xmm5
	movdqa	-28 * SIZE(Y), %xmm7
	
	pshufd  $0xfa,%xmm4,%xmm8
	pshufd  $0x50,%xmm4,%xmm4

	pshufd  $0xfa,%xmm6,%xmm10
	pshufd  $0x50,%xmm6,%xmm6

	pshufd  $0xfa,%xmm5,%xmm9
	pshufd  $0x50,%xmm5,%xmm5

	pshufd  $0xfa,%xmm7,%xmm11
	pshufd  $0x50,%xmm7,%xmm7
	
#ifdef HAVE_SSE4_1
	
	pmuldq   %xmm5, %xmm4
	pmuldq   %xmm7, %xmm6
	pmuldq   %xmm9, %xmm8
	pmuldq   %xmm11, %xmm10
#else
	
#endif
	paddq	%xmm4,%xmm0
	paddq	%xmm6,%xmm1
	paddq	%xmm8,%xmm2
	paddq	%xmm10,%xmm3
	
	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L16:
	testq	$4, N
	jle	.L17

	movdqa	-32 * SIZE(X), %xmm4
	movdqa	-32 * SIZE(Y), %xmm5

	pshufd  $0xfa,%xmm4,%xmm6
	pshufd  $0x50,%xmm4,%xmm4

	pshufd  $0xfa,%xmm5,%xmm7
	pshufd  $0x50,%xmm5,%xmm5
		
#ifdef HAVE_SSE4_1
	pmuldq   %xmm5, %xmm4
	pmuldq   %xmm7, %xmm6
#else
	
#endif
	
        paddq	%xmm4, %xmm0
        paddq   %xmm6, %xmm1
	
	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L17:
	testq	$2, N
	jle	.L18

        movq	-32 * SIZE(X), %xmm6
        movq	-32 * SIZE(Y), %xmm7
	pshufd  $0x50,%xmm6,%xmm6
	pshufd  $0x50,%xmm7,%xmm7
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm7, %xmm6
#else
	
#endif
        paddq   %xmm6,%xmm0
        
	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

.L18:
	testq	$1, N
	jle	.L999

        movd	-32 * SIZE(X), %xmm4
        movd	-32 * SIZE(Y), %xmm5
#ifdef HAVE_SSE4_1
	pmuldq   %xmm5, %xmm4
#else
       /*
	failure to deploy before SSE4_1
       */	
#endif
	paddq	%xmm4, %xmm0

	jmp	.L999
	ALIGN_3

.L20:
#ifdef ALIGNED_ACCESS

	movdqa	-33 * SIZE(X), %xmm4
	addq	$3 * SIZE, X

	movq	N,  %rax
	sarq	$5, %rax
	jle	.L24

	movdqa	-32 * SIZE(X), %xmm5
	movdqa	-28 * SIZE(X), %xmm6
	movdqa	-24 * SIZE(X), %xmm7

	decq	%rax
	jle	.L22

	ALIGN_3

.L21:
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movdqa	-32 * SIZE(Y), %xmm12 // 4/16 
	pshufd  $0xe4,%xmm5, %xmm15
	psrldq  $4, %xmm4
        pslldq  $12, %xmm15
        por     %xmm15, %xmm4

	pshufd  $0xfa,%xmm4,%xmm8
	pshufd  $0x50,%xmm4,%xmm4
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12

#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm4
	pmuldq   %xmm14, %xmm8
#else

#endif
	paddq	%xmm4,%xmm0
	paddq	%xmm8,%xmm1
	movdqa	-20 * SIZE(X), %xmm4 

	movdqa	-28 * SIZE(Y), %xmm12 
	pshufd  $0xe4,%xmm6, %xmm15
	psrldq  $4, %xmm5
        pslldq  $12, %xmm15
        por     %xmm15, %xmm5
	
	pshufd  $0xfa,%xmm5,%xmm8
	pshufd  $0x50,%xmm5,%xmm5
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm5 
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm5,%xmm2
	paddq	%xmm8,%xmm3
	movdqa	-16 * SIZE(X), %xmm5 
	
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	movdqa	-24 * SIZE(Y), %xmm12 //12/16
	pshufd  $0xe4,%xmm7, %xmm15
	psrldq  $4, %xmm6
        pslldq  $12, %xmm15
        por     %xmm15, %xmm6
	
	pshufd  $0xfa,%xmm6,%xmm8
	pshufd  $0x50,%xmm6,%xmm6
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm6
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm6,%xmm0
	paddq	%xmm8,%xmm1
	movdqa	-12 * SIZE(X), %xmm6 
	
	movdqa	-20 * SIZE(Y), %xmm12 
	pshufd  $0xe4,%xmm4, %xmm15
	psrldq  $4, %xmm7
        pslldq  $12, %xmm15
        por     %xmm15, %xmm7
	
	pshufd  $0xfa,%xmm7,%xmm8
	pshufd  $0x50,%xmm7,%xmm7
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm7
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm7,%xmm2
	paddq	%xmm8,%xmm3
	movdqa	-8 * SIZE(X), %xmm7 

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(X)
#endif
	
	movdqa	-16 * SIZE(Y), %xmm12 
	pshufd  $0xe4,%xmm5, %xmm15
	psrldq  $4, %xmm4
        pslldq  $12, %xmm15
        por     %xmm15, %xmm4

	pshufd  $0xfa,%xmm4,%xmm8
	pshufd  $0x50,%xmm4,%xmm4
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12

#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm4
	pmuldq   %xmm14, %xmm8
#else

#endif
	paddq	%xmm4,%xmm0
	paddq	%xmm8,%xmm1
	movdqa	-4 * SIZE(X), %xmm4 
	
	movdqa	-12 * SIZE(Y), %xmm12 
	pshufd  $0xe4,%xmm6, %xmm15
	psrldq  $4, %xmm5
        pslldq  $12, %xmm15
        por     %xmm15, %xmm5
	
	pshufd  $0xfa,%xmm5,%xmm8
	pshufd  $0x50,%xmm5,%xmm5
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm5 
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm5,%xmm2
	paddq	%xmm8,%xmm3
	movdqa	0 * SIZE(X), %xmm5

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(Y)
#endif

	movdqa	-8 * SIZE(Y), %xmm12 
	pshufd  $0xe4,%xmm7, %xmm15
	psrldq  $4, %xmm6
        pslldq  $12, %xmm15
        por     %xmm15, %xmm6
	
	pshufd  $0xfa,%xmm6,%xmm8
	pshufd  $0x50,%xmm6,%xmm6
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm6
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm6,%xmm0
	paddq	%xmm8,%xmm1
	movdqa	4 * SIZE(X), %xmm6 
	
	movdqa	-4 * SIZE(Y), %xmm12 
	pshufd  $0xe4,%xmm4, %xmm15
	psrldq  $4, %xmm7
        pslldq  $12, %xmm15
        por     %xmm15, %xmm7
	
	pshufd  $0xfa,%xmm7,%xmm8
	pshufd  $0x50,%xmm7,%xmm7
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm7
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm7,%xmm2
	paddq	%xmm8,%xmm3
	movdqa	8 * SIZE(X), %xmm7 
	
	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y

	decq	%rax
	jg	.L21
	ALIGN_3

.L22:

	movdqa	-32 * SIZE(Y), %xmm12 
	pshufd  $0xe4,%xmm5, %xmm15
	psrldq  $4, %xmm4
        pslldq  $12, %xmm15
        por     %xmm15, %xmm4

	pshufd  $0xfa,%xmm4,%xmm8
	pshufd  $0x50,%xmm4,%xmm4
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12

#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm4
	pmuldq   %xmm14, %xmm8
#else

#endif
	paddq	%xmm4,%xmm0
	paddq	%xmm8,%xmm1
	movdqa	-20 * SIZE(X), %xmm4 

	movdqa	-28 * SIZE(Y), %xmm12
	pshufd  $0xe4,%xmm6, %xmm15
	psrldq  $4, %xmm5
        pslldq  $12, %xmm15
        por     %xmm15, %xmm5
	
	pshufd  $0xfa,%xmm5,%xmm8
	pshufd  $0x50,%xmm5,%xmm5
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm5 
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm5,%xmm2
	paddq	%xmm8,%xmm3
	movdqa	-16 * SIZE(X), %xmm5 
	
	movdqa	-24 * SIZE(Y), %xmm12 
	pshufd  $0xe4,%xmm7, %xmm15
	psrldq  $4, %xmm6
        pslldq  $12, %xmm15
        por     %xmm15, %xmm6
	
	pshufd  $0xfa,%xmm6,%xmm8
	pshufd  $0x50,%xmm6,%xmm6
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm6
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm6,%xmm0
	paddq	%xmm8,%xmm1
	movdqa	-12 * SIZE(X), %xmm6 
	
	movdqa	-20 * SIZE(Y), %xmm12 
	pshufd  $0xe4,%xmm4, %xmm15
	psrldq  $4, %xmm7
        pslldq  $12, %xmm15
        por     %xmm15, %xmm7
	
	pshufd  $0xfa,%xmm7,%xmm8
	pshufd  $0x50,%xmm7,%xmm7
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm7
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm7,%xmm2
	paddq	%xmm8,%xmm3
	movdqa	-8 * SIZE(X), %xmm7 
		
	movdqa	-16 * SIZE(Y), %xmm12
	pshufd  $0xe4,%xmm5, %xmm15
	psrldq  $4, %xmm4
        pslldq  $12, %xmm15
        por     %xmm15, %xmm4

	pshufd  $0xfa,%xmm4,%xmm8
	pshufd  $0x50,%xmm4,%xmm4
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12

#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm4
	pmuldq   %xmm14, %xmm8
#else

#endif
	paddq	%xmm4,%xmm0
	paddq	%xmm8,%xmm1
	movdqa	-4 * SIZE(X), %xmm4 
	
	movdqa	-12 * SIZE(Y), %xmm12 
	pshufd  $0xe4,%xmm6, %xmm15
	psrldq  $4, %xmm5
        pslldq  $12, %xmm15
        por     %xmm15, %xmm5
	
	pshufd  $0xfa,%xmm5,%xmm8
	pshufd  $0x50,%xmm5,%xmm5
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm5 
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm5,%xmm2
	paddq	%xmm8,%xmm3
	
	movdqa	-8 * SIZE(Y), %xmm12 
	pshufd  $0xe4,%xmm7, %xmm15
	psrldq  $4, %xmm6
        pslldq  $12, %xmm15
        por     %xmm15, %xmm6
	
	pshufd  $0xfa,%xmm6,%xmm8
	pshufd  $0x50,%xmm6,%xmm6
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm6
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm6,%xmm0
	paddq	%xmm8,%xmm1
	
	movdqa	-4 * SIZE(Y), %xmm12 
	pshufd  $0xe4,%xmm4, %xmm15
	psrldq  $4, %xmm7
        pslldq  $12, %xmm15
        por     %xmm15, %xmm7
	
	pshufd  $0xfa,%xmm7,%xmm8
	pshufd  $0x50,%xmm7,%xmm7
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm7
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm7,%xmm2
	paddq	%xmm8,%xmm3
	
	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	ALIGN_3

.L24:
	testq	$31, N
	jle	.L999

	testq	$16, N
	jle	.L25

	movdqa	-32 * SIZE(X), %xmm5 
	movdqa	-28 * SIZE(X), %xmm6 
	movdqa	-24 * SIZE(X), %xmm7 

	movdqa	-32 * SIZE(Y), %xmm12 
	pshufd  $0xe4,%xmm5, %xmm15
	psrldq  $4, %xmm4
        pslldq  $12, %xmm15
        por     %xmm15, %xmm4

	pshufd  $0xfa,%xmm4,%xmm8
	pshufd  $0x50,%xmm4,%xmm4
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12

#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm4
	pmuldq   %xmm14, %xmm8
#else

#endif
	paddq	%xmm4,%xmm0
	paddq	%xmm8,%xmm1
	movdqa	-20 * SIZE(X), %xmm4 
	
	movdqa	-28 * SIZE(Y), %xmm12 
	pshufd  $0xe4,%xmm6, %xmm15
	psrldq  $4, %xmm5
        pslldq  $12, %xmm15
        por     %xmm15, %xmm5
	
	pshufd  $0xfa,%xmm5,%xmm8
	pshufd  $0x50,%xmm5,%xmm5
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm5 
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm5,%xmm2
	paddq	%xmm8,%xmm3

	movdqa	-24 * SIZE(Y), %xmm12 
	pshufd  $0xe4,%xmm7, %xmm15
	psrldq  $4, %xmm6
        pslldq  $12, %xmm15
        por     %xmm15, %xmm6
	
	pshufd  $0xfa,%xmm6,%xmm8
	pshufd  $0x50,%xmm6,%xmm6
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm6
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm6,%xmm0
	paddq	%xmm8,%xmm1

	movdqa	-20 * SIZE(Y), %xmm12 
	pshufd  $0xe4,%xmm4, %xmm15
	psrldq  $4, %xmm7
        pslldq  $12, %xmm15
        por     %xmm15, %xmm7
	
	pshufd  $0xfa,%xmm7,%xmm8
	pshufd  $0x50,%xmm7,%xmm7
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm7
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm7,%xmm2
	paddq	%xmm8,%xmm3
	
	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
	ALIGN_3

.L25:
	testq	$8, N
	jle	.L26

	movdqa	-32 * SIZE(X), %xmm5 
	movdqa	-32 * SIZE(Y), %xmm6
	
	pshufd  $0xe4,%xmm5, %xmm15
	psrldq  $4, %xmm4
        pslldq  $12, %xmm15
        por     %xmm15, %xmm4

	pshufd  $0xfa,%xmm4,%xmm8
	pshufd  $0x50,%xmm4,%xmm4
	pshufd  $0xfa,%xmm6,%xmm10
	pshufd  $0x50,%xmm6,%xmm6

#ifdef HAVE_SSE4_1	
	pmuldq   %xmm6, %xmm4
	pmuldq   %xmm10, %xmm8
#else

#endif
	paddq	%xmm4,%xmm0
	paddq	%xmm8,%xmm1
	
	movdqa	-28 * SIZE(X), %xmm4 
	movdqa	-28 * SIZE(Y), %xmm7
	
	pshufd  $0xe4,%xmm4, %xmm15
	psrldq  $4, %xmm5
        pslldq  $12, %xmm15
        por     %xmm15, %xmm5
	
	pshufd  $0xfa,%xmm5,%xmm9
	pshufd  $0x50,%xmm5,%xmm5
	pshufd  $0xfa,%xmm7,%xmm11
	pshufd  $0x50,%xmm7,%xmm7
	
	
#ifdef HAVE_SSE4_1
	
	pmuldq   %xmm7, %xmm5 
	pmuldq   %xmm11, %xmm9
#else
	
#endif

	paddq	%xmm5,%xmm2
	paddq	%xmm9,%xmm3
	       
	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L26:
	testq	$4, N
	jle	.L27

	movdqa	-32 * SIZE(X), %xmm5
	pshufd  $0xe4,%xmm5, %xmm15
	psrldq  $4, %xmm4
        pslldq  $12, %xmm15
        por     %xmm15, %xmm4

	movdqa	-32 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm4,%xmm13
	pshufd  $0x50,%xmm4,%xmm4
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm4
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm4,  %xmm0
        paddq   %xmm13, %xmm1
	
        pshufd  $0xe4,%xmm5, %xmm4
	
	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L27:
	testq	$2, N
	jle	.L28

	movdqa	-32 * SIZE(Y), %xmm6
	psrldq  $4, %xmm4

	pshufd  $0x50,%xmm4,%xmm5
	pshufd  $0x50,%xmm6,%xmm6
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm6, %xmm5
#else
	
#endif
        paddq   %xmm5,%xmm0
	psrldq  $4, %xmm4
        
	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3
	
.L28:
	testq	$1, N
	jle	.L999

        psrldq  $4, %xmm4
        movd	-32 * SIZE(Y), %xmm5
#ifdef HAVE_SSE4_1
	pmuldq   %xmm5, %xmm4
#else
	
#endif
	paddq	%xmm4, %xmm0

	jmp	.L999
	ALIGN_3

	
.L30:
	testq	$SIZE, X
	jne	.L40

	movdqa  -34 * SIZE(X), %xmm4
	addq	$2 * SIZE, X

	movq	N,  %rax
	sarq	$5, %rax
	jle	.L34

	movdqa	-32 * SIZE(X), %xmm5
	movdqa	-28 * SIZE(X), %xmm6
	movdqa	-24 * SIZE(X), %xmm7
	
	decq	%rax   
	jle	.L32

	ALIGN_3

.L31:
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movdqa	-32 * SIZE(Y), %xmm12
	
	pshufd  $0xfa,%xmm4,%xmm4 //x0,x1
	pshufd  $0x50,%xmm5,%xmm13 //x2,x3
	pshufd  $0xfa,%xmm12,%xmm14 //y2,y3
	pshufd  $0x50,%xmm12,%xmm12 //y1,y2
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm4
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm4,  %xmm0
        paddq   %xmm13, %xmm1
	movdqa	-20 * SIZE(X), %xmm4
		
	movdqa	-28 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm5,%xmm5 //x0,x1
	pshufd  $0x50,%xmm6,%xmm13 //x2,x3
	pshufd  $0xfa,%xmm12,%xmm14 //y2,y3
	pshufd  $0x50,%xmm12,%xmm12 //y1,y2
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm5
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm5,  %xmm2
        paddq   %xmm13, %xmm3
	movdqa	-16 * SIZE(X), %xmm5

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	movdqa	-24 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm6,%xmm6 //x0,x1
	pshufd  $0x50,%xmm7,%xmm13 //x2,x3
	pshufd  $0xfa,%xmm12,%xmm14 //y2,y3
	pshufd  $0x50,%xmm12,%xmm12 //y1,y2
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm6
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm6,  %xmm0
        paddq   %xmm13, %xmm1
	movdqa	-12 * SIZE(X), %xmm6


	movdqa	-20 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm7,%xmm7 //x0,x1
	pshufd  $0x50,%xmm4,%xmm13 //x2,x3
	pshufd  $0xfa,%xmm12,%xmm14 //y2,y3
	pshufd  $0x50,%xmm12,%xmm12 //y1,y2
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm7
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm7,  %xmm2
        paddq   %xmm13, %xmm3
	movdqa	-8 * SIZE(X), %xmm7

	
#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(X)
#endif	

	movdqa	-16 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm4,%xmm4 //x0,x1
	pshufd  $0x50,%xmm5,%xmm13 //x2,x3
	pshufd  $0xfa,%xmm12,%xmm14 //y2,y3
	pshufd  $0x50,%xmm12,%xmm12 //y1,y2
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm4
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm4,  %xmm0
        paddq   %xmm13, %xmm1
	movdqa	-4 * SIZE(X), %xmm4

	movdqa	-12 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm5,%xmm5 //x0,x1
	pshufd  $0x50,%xmm6,%xmm13 //x2,x3
	pshufd  $0xfa,%xmm12,%xmm14 //y2,y3
	pshufd  $0x50,%xmm12,%xmm12 //y1,y2
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm5
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm5,  %xmm2
        paddq   %xmm13, %xmm3
	movdqa	0 * SIZE(X), %xmm5

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(Y)
#endif

	movdqa	-8 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm6,%xmm6 //x0,x1
	pshufd  $0x50,%xmm7,%xmm13 //x2,x3
	pshufd  $0xfa,%xmm12,%xmm14 //y2,y3
	pshufd  $0x50,%xmm12,%xmm12 //y1,y2
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm6
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm6,  %xmm0
        paddq   %xmm13, %xmm1
	movdqa	4 * SIZE(X), %xmm6

	movdqa	-4 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm7,%xmm7 //x0,x1
	pshufd  $0x50,%xmm4,%xmm13 //x2,x3
	pshufd  $0xfa,%xmm12,%xmm14 //y2,y3
	pshufd  $0x50,%xmm12,%xmm12 //y1,y2
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm7
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm7,  %xmm2
        paddq   %xmm13, %xmm3
	movdqa	8 * SIZE(X), %xmm7
	
	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	ALIGN_3

	decq	%rax
	jg	.L31
	ALIGN_3

.L32:
	movdqa	-32 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm4,%xmm4 //x0,x1
	pshufd  $0x50,%xmm5,%xmm13 //x2,x3
	pshufd  $0xfa,%xmm12,%xmm14 //y2,y3
	pshufd  $0x50,%xmm12,%xmm12 //y1,y2
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm4
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm4,  %xmm0
        paddq   %xmm13, %xmm1
	movdqa	-20 * SIZE(X), %xmm4
		
	movdqa	-28 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm5,%xmm5 //x0,x1
	pshufd  $0x50,%xmm6,%xmm13 //x2,x3
	pshufd  $0xfa,%xmm12,%xmm14 //y2,y3
	pshufd  $0x50,%xmm12,%xmm12 //y1,y2
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm5
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm5,  %xmm2
        paddq   %xmm13, %xmm3
	movdqa	-16 * SIZE(X), %xmm5

	movdqa	-24 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm6,%xmm6 //x0,x1
	pshufd  $0x50,%xmm7,%xmm13 //x2,x3
	pshufd  $0xfa,%xmm12,%xmm14 //y2,y3
	pshufd  $0x50,%xmm12,%xmm12 //y1,y2
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm6
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm6,  %xmm0
        paddq   %xmm13, %xmm1
	movdqa	-12 * SIZE(X), %xmm6

	movdqa	-20 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm7,%xmm7 //x0,x1
	pshufd  $0x50,%xmm4,%xmm13 //x2,x3
	pshufd  $0xfa,%xmm12,%xmm14 //y2,y3
	pshufd  $0x50,%xmm12,%xmm12 //y1,y2
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm7
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm7,  %xmm2
        paddq   %xmm13, %xmm3
	movdqa	-8 * SIZE(X), %xmm7

	movdqa	-16 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm4,%xmm4 //x0,x1
	pshufd  $0x50,%xmm5,%xmm13 //x2,x3
	pshufd  $0xfa,%xmm12,%xmm14 //y2,y3
	pshufd  $0x50,%xmm12,%xmm12 //y1,y2
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm4
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm4,  %xmm0
        paddq   %xmm13, %xmm1
	movdqa	-4 * SIZE(X), %xmm4
		
	movdqa	-12 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm5,%xmm5 //x0,x1
	pshufd  $0x50,%xmm6,%xmm13 //x2,x3
	pshufd  $0xfa,%xmm12,%xmm14 //y2,y3
	pshufd  $0x50,%xmm12,%xmm12 //y1,y2
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm5
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm5,  %xmm2
        paddq   %xmm13, %xmm3

	movdqa	-8 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm6,%xmm6 //x0,x1
	pshufd  $0x50,%xmm7,%xmm13 //x2,x3
	pshufd  $0xfa,%xmm12,%xmm14 //y2,y3
	pshufd  $0x50,%xmm12,%xmm12 //y1,y2
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm6
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm6,  %xmm0
        paddq   %xmm13, %xmm1
	
	movdqa	-4 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm7,%xmm7 //x0,x1
	pshufd  $0x50,%xmm4,%xmm13 //x2,x3
	pshufd  $0xfa,%xmm12,%xmm14 //y2,y3
	pshufd  $0x50,%xmm12,%xmm12 //y1,y2
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm7
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm7,  %xmm2
        paddq   %xmm13, %xmm3

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	ALIGN_3

.L34:
	testq	$31, N
	jle	.L999

	testq	$16, N
	jle	.L35

	movdqa	-32 * SIZE(X), %xmm5
	movdqa	-28 * SIZE(X), %xmm6
	movdqa	-24 * SIZE(X), %xmm7

	movdqa	-32 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm4,%xmm4 //x0,x1
	pshufd  $0x50,%xmm5,%xmm13 //x2,x3
	pshufd  $0xfa,%xmm12,%xmm14 //y2,y3
	pshufd  $0x50,%xmm12,%xmm12 //y1,y2
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm4
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm4,  %xmm0
        paddq   %xmm13, %xmm1
	movdqa	-20 * SIZE(X), %xmm4
		
	movdqa	-28 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm5,%xmm5 //x0,x1
	pshufd  $0x50,%xmm6,%xmm13 //x2,x3
	pshufd  $0xfa,%xmm12,%xmm14 //y2,y3
	pshufd  $0x50,%xmm12,%xmm12 //y1,y2
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm5
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm5,  %xmm2
        paddq   %xmm13, %xmm3

	movdqa	-24 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm6,%xmm6 //x0,x1
	pshufd  $0x50,%xmm7,%xmm13 //x2,x3
	pshufd  $0xfa,%xmm12,%xmm14 //y2,y3
	pshufd  $0x50,%xmm12,%xmm12 //y1,y2
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm6
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm6,  %xmm0
        paddq   %xmm13, %xmm1

	movdqa	-20 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm7,%xmm7 //x0,x1
	pshufd  $0x50,%xmm4,%xmm13 //x2,x3
	pshufd  $0xfa,%xmm12,%xmm14 //y2,y3
	pshufd  $0x50,%xmm12,%xmm12 //y1,y2
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm7
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm7,  %xmm2
        paddq   %xmm13, %xmm3
	
	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
	ALIGN_3

.L35:
	testq	$8, N
	jle	.L36

	movdqa	-32 * SIZE(X), %xmm5
	movdqa	-28 * SIZE(X), %xmm6

	movdqa	-32 * SIZE(Y), %xmm12
	
	pshufd  $0xfa,%xmm4,%xmm4 //x0,x1
	pshufd  $0x50,%xmm5,%xmm13 //x2,x3
	pshufd  $0xfa,%xmm12,%xmm14 //y2,y3
	pshufd  $0x50,%xmm12,%xmm12 //y1,y2
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm4
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm4,  %xmm0
        paddq   %xmm13, %xmm1
		
	movdqa	-28 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm5,%xmm5 //x0,x1
	pshufd  $0x50,%xmm6,%xmm13 //x2,x3
	pshufd  $0xfa,%xmm12,%xmm14 //y2,y3
	pshufd  $0x50,%xmm12,%xmm12 //y1,y2
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm5
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm5,  %xmm2
        paddq   %xmm13, %xmm3

	pshufd  $0xe4,%xmm6,%xmm4 //x0,x1
	
	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L36:
	testq	$4, N
	jle	.L37

	movdqa	-32 * SIZE(X), %xmm5
	movdqa	-32 * SIZE(Y), %xmm12

	pshufd  $0xfa,%xmm4,%xmm4 //x0,x1
	pshufd  $0x50,%xmm5,%xmm13 //x2,x3
	pshufd  $0xfa,%xmm12,%xmm14 //y2,y3
	pshufd  $0x50,%xmm12,%xmm12 //y1,y2
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm4
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm4,  %xmm0
        paddq   %xmm13, %xmm1
	
	pshufd  $0xe4,%xmm5,%xmm4 //x0,x1
	
	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L37:
	testq	$2, N
	jle	.L38

        movdqa	-32 * SIZE(Y), %xmm5
	pshufd  $0xfa,%xmm4,%xmm4
	pshufd  $0x50,%xmm5,%xmm5
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm5, %xmm4
#else
	
#endif
        paddq	%xmm4, %xmm0

	movdqa	-32 * SIZE(X), %xmm4
	pslldq  $8, %xmm4
       
	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

.L38:
	testq	$1, N
	jle	.L999

        //movd	-34 * SIZE(X), %xmm4
        movd	-32 * SIZE(Y), %xmm5
	pslldq  $4, %xmm4
	psrldq  $12, %xmm4
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm5, %xmm4
#else
	
#endif
	paddq	%xmm4, %xmm0

	jmp	.L999
	ALIGN_3

.L40:
	movdqa	-35 * SIZE(X), %xmm4
	addq	$SIZE, X

	movq	N,  %rax
	sarq	$5, %rax
	jle	.L44
 
	movdqa	-32 * SIZE(X), %xmm5
	movdqa	-28 * SIZE(X), %xmm6
	movdqa	-24 * SIZE(X), %xmm7

	decq	%rax
	jle	.L42

	ALIGN_3

.L41:
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movdqa	-32 * SIZE(Y), %xmm12 // 4/16 
	pshufd  $0xe4,%xmm5, %xmm15
	psrldq  $12, %xmm4
        pslldq  $4, %xmm15
        por     %xmm15, %xmm4

	pshufd  $0xfa,%xmm4,%xmm8
	pshufd  $0x50,%xmm4,%xmm4
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12

#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm4
	pmuldq   %xmm14, %xmm8
#else

#endif
	paddq	%xmm4,%xmm0
	paddq	%xmm8,%xmm1
	movdqa	-20 * SIZE(X), %xmm4 //4

	movdqa	-28 * SIZE(Y), %xmm12 //8/16
	pshufd  $0xe4,%xmm6, %xmm15
	psrldq  $12, %xmm5
        pslldq  $4, %xmm15
        por     %xmm15, %xmm5
	
	pshufd  $0xfa,%xmm5,%xmm8
	pshufd  $0x50,%xmm5,%xmm5
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm5 
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm5,%xmm2
	paddq	%xmm8,%xmm3
	movdqa	-16 * SIZE(X), %xmm5 
	
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	movdqa	-24 * SIZE(Y), %xmm12 //12/16
	pshufd  $0xe4,%xmm7, %xmm15
	psrldq  $12, %xmm6
        pslldq  $4, %xmm15
        por     %xmm15, %xmm6
	
	pshufd  $0xfa,%xmm6,%xmm8
	pshufd  $0x50,%xmm6,%xmm6
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm6
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm6,%xmm0
	paddq	%xmm8,%xmm1
	movdqa	-12 * SIZE(X), %xmm6 //4
	
	movdqa	-20 * SIZE(Y), %xmm12 //16/16
	pshufd  $0xe4,%xmm4, %xmm15
	psrldq  $12, %xmm7
        pslldq  $4, %xmm15
        por     %xmm15, %xmm7
	
	pshufd  $0xfa,%xmm7,%xmm8
	pshufd  $0x50,%xmm7,%xmm7
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm7
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm7,%xmm2
	paddq	%xmm8,%xmm3
	movdqa	-8 * SIZE(X), %xmm7 //4

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(X)
#endif
	
	movdqa	-16 * SIZE(Y), %xmm12 // 4/16 
	pshufd  $0xe4,%xmm5, %xmm15
	psrldq  $12, %xmm4
        pslldq  $4, %xmm15
        por     %xmm15, %xmm4

	pshufd  $0xfa,%xmm4,%xmm8
	pshufd  $0x50,%xmm4,%xmm4
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12

#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm4
	pmuldq   %xmm14, %xmm8
#else

#endif
	paddq	%xmm4,%xmm0
	paddq	%xmm8,%xmm1
	movdqa	-4 * SIZE(X), %xmm4 //4
	
	movdqa	-12 * SIZE(Y), %xmm12 //8/16
	pshufd  $0xe4,%xmm6, %xmm15
	psrldq  $12, %xmm5
        pslldq  $4, %xmm15
        por     %xmm15, %xmm5
	
	pshufd  $0xfa,%xmm5,%xmm8
	pshufd  $0x50,%xmm5,%xmm5
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm5 
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm5,%xmm2
	paddq	%xmm8,%xmm3
	movdqa	0 * SIZE(X), %xmm5 //4

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(Y)
#endif

	movdqa	-8 * SIZE(Y), %xmm12 //12/16
	pshufd  $0xe4,%xmm7, %xmm15
	psrldq  $12, %xmm6
        pslldq  $4, %xmm15
        por     %xmm15, %xmm6
	
	pshufd  $0xfa,%xmm6,%xmm8
	pshufd  $0x50,%xmm6,%xmm6
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm6
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm6,%xmm0
	paddq	%xmm8,%xmm1
	movdqa	4 * SIZE(X), %xmm6 //4
	
	movdqa	-4 * SIZE(Y), %xmm12 //16/16
	pshufd  $0xe4,%xmm4, %xmm15
	psrldq  $12, %xmm7
        pslldq  $4, %xmm15
        por     %xmm15, %xmm7
	
	pshufd  $0xfa,%xmm7,%xmm8
	pshufd  $0x50,%xmm7,%xmm7
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm7
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm7,%xmm2
	paddq	%xmm8,%xmm3
	movdqa	8 * SIZE(X), %xmm7 //4

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y

	decq	%rax
	jg	.L41
	ALIGN_3

.L42:
	movdqa	-32 * SIZE(Y), %xmm12 // 4/16 
	pshufd  $0xe4,%xmm5, %xmm15
	psrldq  $12, %xmm4
        pslldq  $4, %xmm15
        por     %xmm15, %xmm4

	pshufd  $0xfa,%xmm4,%xmm8
	pshufd  $0x50,%xmm4,%xmm4
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12

#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm4
	pmuldq   %xmm14, %xmm8
#else

#endif
	paddq	%xmm4,%xmm0
	paddq	%xmm8,%xmm1
	movdqa	-20 * SIZE(X), %xmm4 //4

	movdqa	-28 * SIZE(Y), %xmm12 //8/16
	pshufd  $0xe4,%xmm6, %xmm15
	psrldq  $12, %xmm5
        pslldq  $4, %xmm15
        por     %xmm15, %xmm5
	
	pshufd  $0xfa,%xmm5,%xmm8
	pshufd  $0x50,%xmm5,%xmm5
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm5 
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm5,%xmm2
	paddq	%xmm8,%xmm3
	movdqa	-16 * SIZE(X), %xmm5 //4
	
	movdqa	-24 * SIZE(Y), %xmm12 //12/16
	pshufd  $0xe4,%xmm7, %xmm15
	psrldq  $12, %xmm6
        pslldq  $4, %xmm15
        por     %xmm15, %xmm6
	
	pshufd  $0xfa,%xmm6,%xmm8
	pshufd  $0x50,%xmm6,%xmm6
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm6
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm6,%xmm0
	paddq	%xmm8,%xmm1
	movdqa	-12 * SIZE(X), %xmm6 //4
	
	movdqa	-20 * SIZE(Y), %xmm12 //16/16
	pshufd  $0xe4,%xmm4, %xmm15
	psrldq  $12, %xmm7
        pslldq  $4, %xmm15
        por     %xmm15, %xmm7
	
	pshufd  $0xfa,%xmm7,%xmm8
	pshufd  $0x50,%xmm7,%xmm7
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm7
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm7,%xmm2
	paddq	%xmm8,%xmm3
	movdqa	-8 * SIZE(X), %xmm7 //4
		
	movdqa	-16 * SIZE(Y), %xmm12 // 4/16 
	pshufd  $0xe4,%xmm5, %xmm15
	psrldq  $12, %xmm4
        pslldq  $4, %xmm15
        por     %xmm15, %xmm4

	pshufd  $0xfa,%xmm4,%xmm8
	pshufd  $0x50,%xmm4,%xmm4
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12

#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm4
	pmuldq   %xmm14, %xmm8
#else

#endif
	paddq	%xmm4,%xmm0
	paddq	%xmm8,%xmm1
	movdqa	-4 * SIZE(X), %xmm4 //4
	
	movdqa	-12 * SIZE(Y), %xmm12 //8/16
	pshufd  $0xe4,%xmm6, %xmm15
	psrldq  $12, %xmm5
        pslldq  $4, %xmm15
        por     %xmm15, %xmm5
	
	pshufd  $0xfa,%xmm5,%xmm8
	pshufd  $0x50,%xmm5,%xmm5
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm5 
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm5,%xmm2
	paddq	%xmm8,%xmm3
	
	movdqa	-8 * SIZE(Y), %xmm12 //12/16
	pshufd  $0xe4,%xmm7, %xmm15
	psrldq  $12, %xmm6
        pslldq  $4, %xmm15
        por     %xmm15, %xmm6
	
	pshufd  $0xfa,%xmm6,%xmm8
	pshufd  $0x50,%xmm6,%xmm6
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm6
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm6,%xmm0
	paddq	%xmm8,%xmm1
	
	movdqa	-4 * SIZE(Y), %xmm12 //16/16
	pshufd  $0xe4,%xmm4, %xmm15
	psrldq  $12, %xmm7
        pslldq  $4, %xmm15
        por     %xmm15, %xmm7
	
	pshufd  $0xfa,%xmm7,%xmm8
	pshufd  $0x50,%xmm7,%xmm7
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm7
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm7,%xmm2
	paddq	%xmm8,%xmm3
	
	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	ALIGN_3

.L44:
	testq	$31, N
	jle	.L999

	testq	$16, N
	jle	.L45

	movdqa	-32 * SIZE(X), %xmm5 //4
	movdqa	-28 * SIZE(X), %xmm6 //4
	movdqa	-24 * SIZE(X), %xmm7 //4

	movdqa	-32 * SIZE(Y), %xmm12 // 4/16 
	pshufd  $0xe4,%xmm5, %xmm15
	psrldq  $12, %xmm4
        pslldq  $4, %xmm15
        por     %xmm15, %xmm4

	pshufd  $0xfa,%xmm4,%xmm8
	pshufd  $0x50,%xmm4,%xmm4
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12

#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm4
	pmuldq   %xmm14, %xmm8
#else

#endif
	paddq	%xmm4,%xmm0
	paddq	%xmm8,%xmm1
	movdqa	-20 * SIZE(X), %xmm4 //4o
	
	movdqa	-28 * SIZE(Y), %xmm12 //8/16
	pshufd  $0xe4,%xmm6, %xmm15
	psrldq  $12, %xmm5
        pslldq  $4, %xmm15
        por     %xmm15, %xmm5
	
	pshufd  $0xfa,%xmm5,%xmm8
	pshufd  $0x50,%xmm5,%xmm5
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm5 
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm5,%xmm2
	paddq	%xmm8,%xmm3

	movdqa	-24 * SIZE(Y), %xmm12 //12/16
	pshufd  $0xe4,%xmm7, %xmm15
	psrldq  $12, %xmm6
        pslldq  $4, %xmm15
        por     %xmm15, %xmm6
	
	pshufd  $0xfa,%xmm6,%xmm8
	pshufd  $0x50,%xmm6,%xmm6
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm6
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm6,%xmm0
	paddq	%xmm8,%xmm1

	movdqa	-20 * SIZE(Y), %xmm12 //16/16
	pshufd  $0xe4,%xmm4, %xmm15
	psrldq  $12, %xmm7
        pslldq  $4, %xmm15
        por     %xmm15, %xmm7
	
	pshufd  $0xfa,%xmm7,%xmm8
	pshufd  $0x50,%xmm7,%xmm7
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
	
#ifdef HAVE_SSE4_1	
	pmuldq   %xmm12, %xmm7
	pmuldq   %xmm14, %xmm8
#else
	
#endif
	paddq	%xmm7,%xmm2
	paddq	%xmm8,%xmm3
	
	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
	ALIGN_3

.L45:
	testq	$8, N
	jle	.L46

	movdqa	-32 * SIZE(X), %xmm5 //4
	movdqa	-32 * SIZE(Y), %xmm6
	
	pshufd  $0xe4,%xmm5, %xmm15
	psrldq  $12, %xmm4
        pslldq  $4, %xmm15
        por     %xmm15, %xmm4

	pshufd  $0xfa,%xmm4,%xmm8
	pshufd  $0x50,%xmm4,%xmm4
	pshufd  $0xfa,%xmm6,%xmm10
	pshufd  $0x50,%xmm6,%xmm6

#ifdef HAVE_SSE4_1	
	pmuldq   %xmm6, %xmm4
	pmuldq   %xmm10, %xmm8
#else

#endif
	paddq	%xmm4,%xmm0
	paddq	%xmm8,%xmm1
	
	movdqa	-28 * SIZE(X), %xmm4 //4
	movdqa	-28 * SIZE(Y), %xmm7
	
	pshufd  $0xe4,%xmm4, %xmm15
	psrldq  $12, %xmm5
        pslldq  $4, %xmm15
        por     %xmm15, %xmm5
	
	pshufd  $0xfa,%xmm5,%xmm9
	pshufd  $0x50,%xmm5,%xmm5
	pshufd  $0xfa,%xmm7,%xmm11
	pshufd  $0x50,%xmm7,%xmm7
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm7, %xmm5 
	pmuldq   %xmm11, %xmm9
#else
	
#endif
	paddq	%xmm5,%xmm2
	paddq	%xmm9,%xmm3

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L46:
	testq	$4, N
	jle	.L47

	movdqa	-32 * SIZE(X), %xmm5
	pshufd  $0xe4,%xmm5, %xmm15
	psrldq  $12, %xmm4
        pslldq  $4, %xmm15
        por     %xmm15, %xmm4

	movdqa	-32 * SIZE(Y), %xmm12
	pshufd  $0xfa,%xmm4,%xmm13
	pshufd  $0x50,%xmm4,%xmm4
	pshufd  $0xfa,%xmm12,%xmm14
	pshufd  $0x50,%xmm12,%xmm12
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm4
	pmuldq   %xmm14, %xmm13
#else
	
#endif
        paddq	%xmm4,  %xmm0
        paddq   %xmm13, %xmm1
	
        pshufd  $0xe4,%xmm5, %xmm4


	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L47:
	testq	$2, N
	jle	.L48

	movdqa	-32 * SIZE(X), %xmm5
	pshufd  $0xe4,%xmm5, %xmm15
	psrldq  $12, %xmm4
        pslldq  $4, %xmm15
        por     %xmm15, %xmm4

	movdqa	-32 * SIZE(Y), %xmm12
	pshufd  $0x50,%xmm4,%xmm4
	pshufd  $0x50,%xmm12,%xmm12
#ifdef HAVE_SSE4_1
	pmuldq   %xmm12, %xmm4
#else
	
#endif
        paddq	%xmm4,  %xmm0
        
        pshufd  $0x55,%xmm5, %xmm4
	
	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

.L48:
	testq	$1, N
	jle	.L999

        movd	-32 * SIZE(Y), %xmm5
        psrldq  $12, %xmm4
#ifdef HAVE_SSE4_1
	pmuldq   %xmm5, %xmm4
#else
	
#endif
	paddq	%xmm4, %xmm0
	
	jmp	.L999
	ALIGN_4

#else /* L20 NO ALIGNED_ACCESS*/

	movq	N,  %rax
	sarq	$5, %rax
	jle	.L24

	movq	-32 * SIZE(X), %xmm4
	movq	-30 * SIZE(X), %xmm12
        punpcklqdq %xmm12, %xmm4
	movq	-28 * SIZE(X), %xmm5
	movq	-26 * SIZE(X), %xmm13
        punpcklqdq %xmm13, %xmm5
	movq	-24 * SIZE(X), %xmm6
	movq	-22 * SIZE(X), %xmm14
        punpcklqdq %xmm14, %xmm6
	movq	-20 * SIZE(X), %xmm7
	movq	-18 * SIZE(X), %xmm15
        punpcklqdq %xmm15, %xmm7

	movq	-16 * SIZE(X), %xmm8
	movq	-14 * SIZE(X), %xmm12
        punpcklqdq %xmm12, %xmm8
	movq	-12 * SIZE(X), %xmm9
	movq	-10 * SIZE(X), %xmm13
        punpcklqdq %xmm13, %xmm9
	movq	 -8 * SIZE(X), %xmm10
	movq	 -6 * SIZE(X), %xmm14
        punpcklqdq %xmm14, %xmm10
	movq	 -4 * SIZE(X), %xmm11
	movq	 -2 * SIZE(X), %xmm15
        punpcklqdq %xmm15, %xmm11

	decq	%rax
	jle	.L22

	ALIGN_3

.L21:
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

        movdqa	-32 * SIZE(Y), %xmm4
#ifdef HAVE_SSE3
	pabsd   %xmm4, %xmm4
#else
	
#endif
        pshufd  $0xf5, %xmm4, %xmm12
        pmuludq %xmm4, %xmm4
        pmuludq %xmm12, %xmm12
        paddq   %xmm4,  %xmm0
	paddq	%xmm12, %xmm1
        movq	0 * SIZE(X), %xmm4
	movq	2 * SIZE(X), %xmm12
        punpcklqdq %xmm12, %xmm4

        movdqa	-28 * SIZE(Y), %xmm5
#ifdef HAVE_SSE3
	pabsd   %xmm5, %xmm5
#else
	
#endif
        pshufd  $0xf5, %xmm5, %xmm13
        pmuludq %xmm5, %xmm5
        pmuludq %xmm13, %xmm13
        paddq   %xmm5,  %xmm2
	paddq	%xmm13, %xmm3
        movq	4 * SIZE(X), %xmm5
	movq	6 * SIZE(X), %xmm13
        punpcklqdq %xmm13, %xmm5
	
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

        movdqa	-24 * SIZE(Y), %xmm6
#ifdef HAVE_SSE3
	pabsd   %xmm6, %xmm6
#else
	
#endif
        pshufd $0xf5,  %xmm6, %xmm14
        pmuludq %xmm6, %xmm6
        pmuludq %xmm14,%xmm14
        paddq	%xmm6, %xmm0
        paddq   %xmm14,%xmm1
        movq	8 * SIZE(X), %xmm6
	movq	10 * SIZE(X), %xmm14
        punpcklqdq %xmm14, %xmm6

        movdqa	-20 * SIZE(Y), %xmm7
#ifdef HAVE_SSE3
	pabsd   %xmm7, %xmm7
#else
	
#endif
        pshufd $0xf5,  %xmm7, %xmm15
        pmuludq %xmm7, %xmm7
        pmuludq %xmm15,%xmm15
        paddq	%xmm7, %xmm2
        paddq   %xmm15,%xmm3
        movq	12 * SIZE(X), %xmm7
	movq	14 * SIZE(X), %xmm15
        punpcklqdq %xmm15, %xmm7
       
#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(X)
#endif

        movdqa	-16 * SIZE(Y), %xmm8
#ifdef HAVE_SSE3
	pabsd   %xmm8, %xmm8
#else
	
#endif
        pshufd  $0xf5, %xmm8, %xmm12
        pmuludq %xmm8, %xmm8
        pmuludq %xmm12, %xmm12
        paddq	%xmm8,  %xmm0
        paddq   %xmm12, %xmm1
        movq	16 * SIZE(X), %xmm8
	movq	18 * SIZE(X), %xmm12
        punpcklqdq %xmm12, %xmm8

        movdqa	-12 * SIZE(Y), %xmm9
#ifdef HAVE_SSE3
	pabsd   %xmm9, %xmm9
#else
	
#endif
        pshufd $0xf5,  %xmm9, %xmm13
        pmuludq %xmm9, %xmm9
        pmuludq %xmm13,%xmm13
        paddq	%xmm9, %xmm2
        paddq   %xmm13,%xmm3
        movq	20 * SIZE(X), %xmm9
	movq	22 * SIZE(X), %xmm13
        punpcklqdq %xmm13, %xmm9

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(Y)
#endif

        movdqa	-8 * SIZE(Y), %xmm10
#ifdef HAVE_SSE3
	pabsd   %xmm10, %xmm10
#else
	
#endif
        pshufd  $0xf5, %xmm10, %xmm14
        pmuludq %xmm10,%xmm10
        pmuludq %xmm14,%xmm14
        paddq	%xmm10,%xmm0
        paddq   %xmm14,%xmm1
        movq	24 * SIZE(X), %xmm10
	movq	26 * SIZE(X), %xmm14
        punpcklqdq %xmm14, %xmm10

        movdqa	-4 * SIZE(Y), %xmm11
#ifdef HAVE_SSE3
	pabsd   %xmm11, %xmm11
#else
	
#endif
        pshufd  $0xf5, %xmm11, %xmm15
        pmuludq %xmm11,%xmm11
        pmuludq %xmm15,%xmm15
        paddq	%xmm11,%xmm2
        paddq   %xmm15,%xmm3
        movq	28 * SIZE(X), %xmm11
	movq	30 * SIZE(X), %xmm15
        punpcklqdq %xmm15, %xmm11
       
	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y

	decq	%rax
	jg	.L21
	ALIGN_3

.L22:
        movdqa	-32 * SIZE(Y), %xmm4
#ifdef HAVE_SSE3
	pabsd   %xmm4, %xmm4
#else
	
#endif
        pshufd  $0xf5, %xmm4, %xmm12
        pmuludq %xmm4, %xmm4
        pmuludq %xmm12, %xmm12
        paddq	%xmm4,  %xmm0
        paddq   %xmm12, %xmm1
	
	movdqa	-28 * SIZE(Y), %xmm5
#ifdef HAVE_SSE3
	pabsd   %xmm5, %xmm5
#else
	
#endif
        pshufd  $0xf5, %xmm5, %xmm13
        pmuludq %xmm5, %xmm5
        pmuludq %xmm13, %xmm13
        paddq	%xmm5, %xmm2
        paddq   %xmm13, %xmm3
	
	movdqa	-24 * SIZE(Y), %xmm6
#ifdef HAVE_SSE3
	pabsd   %xmm6, %xmm6
#else
	
#endif
        pshufd $0xf5, %xmm6, %xmm14
        pmuludq %xmm6,%xmm6
        pmuludq %xmm14,%xmm14
        paddq	%xmm6,%xmm0
        paddq   %xmm14,%xmm1
	
        movdqa	-20 * SIZE(Y), %xmm7
#ifdef HAVE_SSE3
	pabsd   %xmm7, %xmm7
#else
	
#endif
        pshufd $0xf5, %xmm7, %xmm15
        pmuludq %xmm7,%xmm7
        pmuludq %xmm15,%xmm15
        paddq	%xmm7,%xmm2
        paddq   %xmm15,%xmm3
	
	movdqa	-16 * SIZE(Y), %xmm8
#ifdef HAVE_SSE3
	pabsd   %xmm8, %xmm8
#else
	
#endif
        pshufd  $0xf5, %xmm8, %xmm12
        pmuludq %xmm8, %xmm8
        pmuludq %xmm12, %xmm12
        paddq	%xmm8,  %xmm0
        paddq   %xmm12, %xmm1
	
	movdqa	-12 * SIZE(Y), %xmm9
#ifdef HAVE_SSE3
	pabsd   %xmm9, %xmm9
#else
	
#endif
        pshufd $0xf5, %xmm9, %xmm13
        pmuludq %xmm9, %xmm9
        pmuludq %xmm13, %xmm13
        paddq	%xmm9, %xmm2
        paddq   %xmm13, %xmm3
	
	movdqa	 -8 * SIZE(Y), %xmm10
#ifdef HAVE_SSE3
	pabsd   %xmm10, %xmm10
#else
	
#endif
        pshufd $0xf5, %xmm10, %xmm14
        pmuludq %xmm10,%xmm10
        pmuludq %xmm14,%xmm14
        paddq	%xmm10,%xmm0
        paddq   %xmm14,%xmm1
	
	movdqa	-4 * SIZE(Y), %xmm11
#ifdef HAVE_SSE3
	pabsd   %xmm11, %xmm11
#else
	
#endif
        pshufd  $0xf5, %xmm11, %xmm15
        pmuludq %xmm11,%xmm11
        pmuludq %xmm15,%xmm15
        paddq	%xmm11,%xmm2
        paddq   %xmm15,%xmm3

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	ALIGN_3

.L24:
	testq	$31, N
	jle	.L999

	testq	$16, N
	jle	.L25

	movq	-32 * SIZE(X), %xmm4
	movq	-30 * SIZE(X), %xmm12
        punpcklqdq %xmm12, %xmm4
	movq	-28 * SIZE(X), %xmm5
	movq	-26 * SIZE(X), %xmm13
        punpcklqdq %xmm13, %xmm5
	movq	-24 * SIZE(X), %xmm6
	movq	-22 * SIZE(X), %xmm14
        punpcklqdq %xmm14, %xmm6
	movq	-20 * SIZE(X), %xmm7
	movq	-18 * SIZE(X), %xmm15
        punpcklqdq %xmm15, %xmm7

        movdqa	-32 * SIZE(Y), %xmm4
#ifdef HAVE_SSE3
	pabsd   %xmm4, %xmm4
#else
	
#endif
        pshufd  $0xf5, %xmm4, %xmm8
        pmuludq %xmm4,%xmm4
        pmuludq %xmm8,%xmm8
        paddq	%xmm4,%xmm0
        paddq   %xmm8,%xmm1
        
	movdqa  -28 * SIZE(Y), %xmm5
#ifdef HAVE_SSE3
	pabsd   %xmm5, %xmm5
#else
	
#endif
        pshufd $0xf5, %xmm5, %xmm9
        pmuludq %xmm5,%xmm5
        pmuludq %xmm9,%xmm9
        paddq	%xmm5,%xmm2
        paddq   %xmm9,%xmm3
	
	movdqa	-24 * SIZE(Y), %xmm6
#ifdef HAVE_SSE3
	pabsd   %xmm6, %xmm6
#else
	
#endif
        pshufd $0xf5, %xmm6, %xmm10
        pmuludq %xmm6,%xmm6
        pmuludq %xmm10,%xmm10
        paddq	%xmm6, %xmm0
        paddq   %xmm10,%xmm1
	
	movdqa	-20 * SIZE(Y), %xmm7
#ifdef HAVE_SSE3
	pabsd   %xmm7, %xmm7
#else
	
#endif
        pshufd $0xf5, %xmm7, %xmm11
        pmuludq %xmm7,%xmm7
        pmuludq %xmm11,%xmm11
        paddq	%xmm7,%xmm2
        paddq   %xmm11,%xmm3

	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
	ALIGN_3

.L25:
	testq	$8, N
	jle	.L26

	movq	-32 * SIZE(X), %xmm4
	movq	-30 * SIZE(X), %xmm12
        punpcklqdq %xmm12, %xmm4
	movq	-28 * SIZE(X), %xmm5
	movq	-26 * SIZE(X), %xmm13
        punpcklqdq %xmm13, %xmm5

        movdqa  -32 * SIZE(Y), %xmm4
#ifdef HAVE_SSE3
	pabsd   %xmm4, %xmm4
#else
	
#endif
        pshufd $0xf5, %xmm4, %xmm6
        pmuludq %xmm4,%xmm4
        pmuludq %xmm6,%xmm6
        paddq	%xmm4,%xmm0
        paddq   %xmm6,%xmm1
	
	movdqa	-28 * SIZE(Y), %xmm5
#ifdef HAVE_SSE3
	pabsd   %xmm5, %xmm5
#else
	
#endif
        pshufd $0xf5, %xmm5,%xmm7
        pmuludq %xmm5,%xmm5
        pmuludq %xmm7,%xmm7
        paddq	%xmm5,%xmm2
        paddq   %xmm7,%xmm3

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

.L26:
	testq	$4, N
	jle	.L27

        movq	-32 * SIZE(X), %xmm4
        movq	-30 * SIZE(X), %xmm5

	pshufd  $0x50,%xmm4,%xmm4
	pshufd  $0x50,%xmm5,%xmm5
       
        movdqa	-32 * SIZE(Y), %xmm6
	pshufd  $0xfa,%xmm6,%xmm7
	pshufd  $0x50,%xmm6,%xmm6

#ifdef HAVE_SSE4_1
	pmuldq   %xmm6, %xmm4
	pmuldq   %xmm7, %xmm5
#else
	
#endif
	
        paddq	%xmm4, %xmm0
        paddq   %xmm5, %xmm1
	        
	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

.L27:
	testq	$2, N
	jle	.L28

        movq	-32 * SIZE(X), %xmm6
        movq	-32 * SIZE(Y), %xmm7

	pshufd  $0x50,%xmm6,%xmm6
	pshufd  $0x50,%xmm7,%xmm7
#ifdef HAVE_SSE4_1
	pmuldq   %xmm7, %xmm6
#else
	
#endif
        paddq   %xmm6,%xmm0

	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

.L28:
	testq	$1, N
	jle	.L999

        movd	-32 * SIZE(X), %xmm4
        movd	-32 * SIZE(Y), %xmm5

#ifdef HAVE_SSE4_1
	pmuldq   %xmm5, %xmm4
#else
	
#endif
	paddq	%xmm4, %xmm0

	jmp	.L999
	ALIGN_3
	
#endif /*END ALINEADO*/

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
	punpcklqdq %xmm5, %xmm4
        
        movd	0 * SIZE(X), %xmm6
        addq	INCX, X
        movd	0 * SIZE(X), %xmm7
        addq	INCX, X
	punpcklqdq %xmm7, %xmm6
        
        movd	0 * SIZE(Y), %xmm8
        addq	INCY, Y
        movd	0 * SIZE(Y), %xmm9
        addq	INCY, Y
	punpcklqdq %xmm9, %xmm8
	
        movd	0 * SIZE(Y), %xmm10
        addq	INCY, Y
        movd	0 * SIZE(Y), %xmm11
        addq	INCY, Y
	punpcklqdq %xmm11, %xmm10
        
#ifdef HAVE_SSE4_1
	pmuldq   %xmm8, %xmm4
	pmuldq   %xmm10, %xmm6
#else
	
#endif
        paddq   %xmm4,%xmm0
	paddq   %xmm6,%xmm1
	  
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
        movd	0 * SIZE(Y), %xmm5
	addq	INCY, Y
	
#ifdef HAVE_SSE4_1
	pmuldq   %xmm5, %xmm4
#else
	
#endif
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
     
	ret
.LFE2:
	.size	idot_kh, .-idot_kh
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
