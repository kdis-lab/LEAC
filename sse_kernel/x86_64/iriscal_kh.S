/*! \file iriscal_kh.S	
 * This file is part of the LEAC.
 *
 * Implementation of the 
 *	
 *	x <-- a x, ( int <-- a int )
 *
 *  function for unsigned long int
 *
 *
 * (c)  Hermes Robles Berumen <hermes@uaz.edu.mx>
 *
 * For the full copyright and license information, please view the LICENSE
 * file that was distributed with this source code.
 */
	
	.file	"iriscal_kh.S"
	.section	.text.unlikely,"ax",@progbits
.LCOLDB0:
	.text
.LHOTB0:
	.p2align 4,,15
	.globl	iriscal_kh
	.type	iriscal_kh, @function
iriscal_kh:
.LFB0:
	.cfi_startproc
	cmpq	$1, %r8
	pushq	%rbx
	.cfi_def_cfa_offset 16
	.cfi_offset 3, -16
	je	.L2
	salq	$2, %r8
	testq	%rdi, %rdi
	jle	.L20
	pxor	%xmm1, %xmm1
	leaq	-1(%rdi), %rdx
	movl	$1, %ebx
	andl	$7, %edx
	cvtsi2sd	(%rcx), %xmm1
	mulsd	%xmm0, %xmm1
	cvttsd2si	%xmm1, %eax
	movl	%eax, (%rcx)
	addq	%r8, %rcx
	cmpq	%rdi, %rbx
	je	.L20
	testq	%rdx, %rdx
	je	.L21
	cmpq	$1, %rdx
	je	.L61
	cmpq	$2, %rdx
	je	.L62
	cmpq	$3, %rdx
	je	.L63
	cmpq	$4, %rdx
	je	.L64
	cmpq	$5, %rdx
	je	.L65
	cmpq	$6, %rdx
	je	.L66
	pxor	%xmm2, %xmm2
	addq	$1, %rbx
	cvtsi2sd	(%rcx), %xmm2
	mulsd	%xmm0, %xmm2
	cvttsd2si	%xmm2, %esi
	movl	%esi, (%rcx)
	addq	%r8, %rcx
.L66:
	pxor	%xmm3, %xmm3
	addq	$1, %rbx
	cvtsi2sd	(%rcx), %xmm3
	mulsd	%xmm0, %xmm3
	cvttsd2si	%xmm3, %r9d
	movl	%r9d, (%rcx)
	addq	%r8, %rcx
.L65:
	pxor	%xmm4, %xmm4
	addq	$1, %rbx
	cvtsi2sd	(%rcx), %xmm4
	mulsd	%xmm0, %xmm4
	cvttsd2si	%xmm4, %r10d
	movl	%r10d, (%rcx)
	addq	%r8, %rcx
.L64:
	pxor	%xmm5, %xmm5
	addq	$1, %rbx
	cvtsi2sd	(%rcx), %xmm5
	mulsd	%xmm0, %xmm5
	cvttsd2si	%xmm5, %r11d
	movl	%r11d, (%rcx)
	addq	%r8, %rcx
.L63:
	pxor	%xmm6, %xmm6
	addq	$1, %rbx
	cvtsi2sd	(%rcx), %xmm6
	mulsd	%xmm0, %xmm6
	cvttsd2si	%xmm6, %edx
	movl	%edx, (%rcx)
	addq	%r8, %rcx
.L62:
	pxor	%xmm7, %xmm7
	addq	$1, %rbx
	cvtsi2sd	(%rcx), %xmm7
	mulsd	%xmm0, %xmm7
	cvttsd2si	%xmm7, %eax
	movl	%eax, (%rcx)
	addq	%r8, %rcx
.L61:
	pxor	%xmm8, %xmm8
	addq	$1, %rbx
	cvtsi2sd	(%rcx), %xmm8
	mulsd	%xmm0, %xmm8
	cvttsd2si	%xmm8, %esi
	movl	%esi, (%rcx)
	addq	%r8, %rcx
	cmpq	%rdi, %rbx
	je	.L20
.L21:
	pxor	%xmm9, %xmm9
	addq	$8, %rbx
	pxor	%xmm10, %xmm10
	pxor	%xmm11, %xmm11
	cvtsi2sd	(%rcx), %xmm9
	mulsd	%xmm0, %xmm9
	pxor	%xmm12, %xmm12
	pxor	%xmm13, %xmm13
	pxor	%xmm14, %xmm14
	cvttsd2si	%xmm9, %r9d
	pxor	%xmm15, %xmm15
	pxor	%xmm1, %xmm1
	movl	%r9d, (%rcx)
	addq	%r8, %rcx
	cvtsi2sd	(%rcx), %xmm10
	mulsd	%xmm0, %xmm10
	cvttsd2si	%xmm10, %r10d
	movl	%r10d, (%rcx)
	addq	%r8, %rcx
	cvtsi2sd	(%rcx), %xmm11
	mulsd	%xmm0, %xmm11
	cvttsd2si	%xmm11, %r11d
	movl	%r11d, (%rcx)
	addq	%r8, %rcx
	cvtsi2sd	(%rcx), %xmm12
	mulsd	%xmm0, %xmm12
	cvttsd2si	%xmm12, %edx
	movl	%edx, (%rcx)
	addq	%r8, %rcx
	cvtsi2sd	(%rcx), %xmm13
	mulsd	%xmm0, %xmm13
	cvttsd2si	%xmm13, %eax
	movl	%eax, (%rcx)
	addq	%r8, %rcx
	cvtsi2sd	(%rcx), %xmm14
	mulsd	%xmm0, %xmm14
	cvttsd2si	%xmm14, %esi
	movl	%esi, (%rcx)
	addq	%r8, %rcx
	cvtsi2sd	(%rcx), %xmm15
	mulsd	%xmm0, %xmm15
	cvttsd2si	%xmm15, %r9d
	movl	%r9d, (%rcx)
	addq	%r8, %rcx
	cvtsi2sd	(%rcx), %xmm1
	mulsd	%xmm0, %xmm1
	cvttsd2si	%xmm1, %r10d
	movl	%r10d, (%rcx)
	addq	%r8, %rcx
	cmpq	%rdi, %rbx
	jne	.L21
	.p2align 4,,10
	.p2align 3
.L20:
	movl	$1, %eax
	popq	%rbx
	.cfi_remember_state
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L2:
	.cfi_restore_state
	testq	%rdi, %rdi
	jle	.L20
	movq	%rcx, %rax
	andl	$15, %eax
	shrq	$2, %rax
	negq	%rax
	andl	$3, %eax
	cmpq	%rdi, %rax
	cmova	%rdi, %rax
	cmpq	$4, %rdi
	jg	.L82
	movq	%rdi, %rax
.L15:
	pxor	%xmm2, %xmm2
	cmpq	$1, %rax
	movl	$1, %edx
	cvtsi2sd	(%rcx), %xmm2
	mulsd	%xmm0, %xmm2
	cvttsd2si	%xmm2, %r8d
	movl	%r8d, (%rcx)
	je	.L7
	pxor	%xmm3, %xmm3
	cmpq	$2, %rax
	movl	$2, %edx
	cvtsi2sd	4(%rcx), %xmm3
	mulsd	%xmm0, %xmm3
	cvttsd2si	%xmm3, %ebx
	movl	%ebx, 4(%rcx)
	je	.L7
	pxor	%xmm4, %xmm4
	cmpq	$4, %rax
	movl	$3, %edx
	cvtsi2sd	8(%rcx), %xmm4
	mulsd	%xmm0, %xmm4
	cvttsd2si	%xmm4, %r11d
	movl	%r11d, 8(%rcx)
	jne	.L7
	pxor	%xmm5, %xmm5
	cvtsi2sd	12(%rcx), %xmm5
	mulsd	%xmm0, %xmm5
	cvttsd2si	%xmm5, %edx
	movl	%edx, 12(%rcx)
	movl	$4, %edx
.L7:
	cmpq	%rdi, %rax
	je	.L20
.L6:
	movq	%rdi, %r11
	leaq	-1(%rdi), %r9
	subq	%rax, %r11
	leaq	-4(%r11), %rsi
	subq	%rax, %r9
	shrq	$2, %rsi
	cmpq	$2, %r9
	leaq	1(%rsi), %r10
	leaq	0(,%r10,4), %rbx
	jbe	.L9
	leaq	(%rcx,%rax,4), %rax
	movapd	%xmm0, %xmm8
	movl	$1, %r9d
	andl	$3, %esi
	movl	$16, %r8d
	movdqa	(%rax), %xmm6
	cmpq	%r10, %r9
	unpcklpd	%xmm8, %xmm8
	pshufd	$238, %xmm6, %xmm9
	cvtdq2pd	%xmm6, %xmm7
	mulpd	%xmm8, %xmm7
	cvtdq2pd	%xmm9, %xmm10
	mulpd	%xmm8, %xmm10
	cvttpd2dq	%xmm7, %xmm11
	cvttpd2dq	%xmm10, %xmm12
	movdqa	%xmm11, %xmm13
	punpcklqdq	%xmm12, %xmm13
	movaps	%xmm13, (%rax)
	jnb	.L79
	testq	%rsi, %rsi
	je	.L10
	cmpq	$1, %rsi
	je	.L59
	cmpq	$2, %rsi
	je	.L60
	movdqa	16(%rax), %xmm15
	movl	$2, %r9d
	movl	$32, %r8d
	pshufd	$238, %xmm15, %xmm1
	cvtdq2pd	%xmm15, %xmm2
	mulpd	%xmm8, %xmm2
	cvtdq2pd	%xmm1, %xmm3
	mulpd	%xmm8, %xmm3
	cvttpd2dq	%xmm2, %xmm4
	cvttpd2dq	%xmm3, %xmm6
	movdqa	%xmm4, %xmm5
	punpcklqdq	%xmm6, %xmm5
	movaps	%xmm5, 16(%rax)
.L60:
	addq	$1, %r9
	movdqa	(%rax,%r8), %xmm9
	pshufd	$238, %xmm9, %xmm11
	cvtdq2pd	%xmm9, %xmm10
	mulpd	%xmm8, %xmm10
	cvtdq2pd	%xmm11, %xmm12
	mulpd	%xmm8, %xmm12
	cvttpd2dq	%xmm10, %xmm13
	cvttpd2dq	%xmm12, %xmm14
	movdqa	%xmm13, %xmm15
	punpcklqdq	%xmm14, %xmm15
	movaps	%xmm15, (%rax,%r8)
	addq	$16, %r8
.L59:
	addq	$1, %r9
	movdqa	(%rax,%r8), %xmm1
	pshufd	$238, %xmm1, %xmm4
	cvtdq2pd	%xmm1, %xmm3
	mulpd	%xmm8, %xmm3
	cvtdq2pd	%xmm4, %xmm6
	mulpd	%xmm8, %xmm6
	cvttpd2dq	%xmm3, %xmm5
	cvttpd2dq	%xmm6, %xmm9
	movdqa	%xmm5, %xmm7
	punpcklqdq	%xmm9, %xmm7
	movaps	%xmm7, (%rax,%r8)
	addq	$16, %r8
	cmpq	%r10, %r9
	jnb	.L79
.L10:
	addq	$4, %r9
	movdqa	(%rax,%r8), %xmm11
	pshufd	$238, %xmm11, %xmm13
	cvtdq2pd	%xmm11, %xmm12
	movdqa	16(%rax,%r8), %xmm4
	mulpd	%xmm8, %xmm12
	pshufd	$238, %xmm4, %xmm5
	cvtdq2pd	%xmm13, %xmm14
	cvtdq2pd	%xmm4, %xmm6
	mulpd	%xmm8, %xmm14
	movdqa	32(%rax,%r8), %xmm13
	mulpd	%xmm8, %xmm6
	cvtdq2pd	%xmm5, %xmm9
	movdqa	48(%rax,%r8), %xmm5
	mulpd	%xmm8, %xmm9
	cvttpd2dq	%xmm12, %xmm15
	cvttpd2dq	%xmm14, %xmm2
	cvtdq2pd	%xmm13, %xmm14
	cvttpd2dq	%xmm6, %xmm7
	cvttpd2dq	%xmm9, %xmm10
	cvtdq2pd	%xmm5, %xmm9
	movdqa	%xmm15, %xmm1
	pshufd	$238, %xmm13, %xmm15
	punpcklqdq	%xmm2, %xmm1
	mulpd	%xmm8, %xmm14
	movdqa	%xmm7, %xmm11
	pshufd	$238, %xmm5, %xmm7
	mulpd	%xmm8, %xmm9
	movaps	%xmm1, (%rax,%r8)
	cvtdq2pd	%xmm15, %xmm1
	punpcklqdq	%xmm10, %xmm11
	cvtdq2pd	%xmm7, %xmm10
	mulpd	%xmm8, %xmm1
	mulpd	%xmm8, %xmm10
	movaps	%xmm11, 16(%rax,%r8)
	cvttpd2dq	%xmm14, %xmm2
	cvttpd2dq	%xmm9, %xmm11
	cvttpd2dq	%xmm1, %xmm3
	cvttpd2dq	%xmm10, %xmm12
	movdqa	%xmm2, %xmm4
	movdqa	%xmm11, %xmm13
	punpcklqdq	%xmm3, %xmm4
	punpcklqdq	%xmm12, %xmm13
	movaps	%xmm4, 32(%rax,%r8)
	movaps	%xmm13, 48(%rax,%r8)
	addq	$64, %r8
	cmpq	%r10, %r9
	jb	.L10
	.p2align 4,,10
	.p2align 3
.L79:
	addq	%rbx, %rdx
	cmpq	%rbx, %r11
	je	.L20
.L9:
	leaq	0(,%rdx,4), %rsi
	pxor	%xmm8, %xmm8
	leaq	1(%rdx), %r9
	leaq	(%rcx,%rsi), %r10
	cmpq	%r9, %rdi
	cvtsi2sd	(%r10), %xmm8
	mulsd	%xmm0, %xmm8
	cvttsd2si	%xmm8, %eax
	movl	%eax, (%r10)
	jle	.L20
	pxor	%xmm15, %xmm15
	leaq	4(%rcx,%rsi), %r8
	addq	$2, %rdx
	cmpq	%rdx, %rdi
	cvtsi2sd	(%r8), %xmm15
	mulsd	%xmm0, %xmm15
	cvttsd2si	%xmm15, %ebx
	movl	%ebx, (%r8)
	jle	.L20
	pxor	%xmm1, %xmm1
	leaq	8(%rcx,%rsi), %rcx
	movl	$1, %eax
	cvtsi2sd	(%rcx), %xmm1
	mulsd	%xmm1, %xmm0
	cvttsd2si	%xmm0, %edi
	movl	%edi, (%rcx)
	popq	%rbx
	.cfi_remember_state
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L82:
	.cfi_restore_state
	testq	%rax, %rax
	jne	.L15
	xorl	%edx, %edx
	jmp	.L6
	.cfi_endproc
.LFE0:
	.size	iriscal_kh, .-iriscal_kh
	.section	.text.unlikely
.LCOLDE0:
	.text
.LHOTE0:
	.ident	"GCC: (Debian 4.9.2-10) 4.9.2"
	.section	.note.GNU-stack,"",@progbits
