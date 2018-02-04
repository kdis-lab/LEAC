/*! \file saasxpa_kh.S	
 * This file is part of the LEAC.
 *
 * Implementation of the 
 *	
 *	a_ij = a_ij + alpha(a_ij - xj),
 *
 *  function for float
 *
 *  Based on Optimized BLAS libraries
 *
 * (c)  Hermes Robles Berumen <hermes@uaz.edu.mx>
 *
 * For the full copyright and license information, please view the LICENSE
 * file that was distributed with this source code.
 */

	.file	"saasxpa_kh.S"
	.section	.text.unlikely,"ax",@progbits
.LCOLDB0:
	.text
.LHOTB0:
	.p2align 4,,15
	.globl	saasxpa_kh
	.type	saasxpa_kh, @function
saasxpa_kh:
.LFB0:
	.cfi_startproc
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	testq	%rsi, %rsi
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	leaq	(%rcx,%rdx,4), %rbp
	movq	%rsi, -48(%rsp)
	jle	.L1
	leaq	4(%rcx), %rdx
	movq	%rbp, %r12
	leaq	16(%rcx), %rsi
	movaps	%xmm0, %xmm3
	subq	%rdx, %r12
	movq	%rdx, -16(%rsp)
	movq	%rsi, -32(%rsp)
	shrq	$2, %r12
	shufps	$0, %xmm3, %xmm3
	addq	$1, %r12
	leaq	0(,%r12,4), %rbx
	cmpq	$6, %r12
	seta	-17(%rsp)
	movq	%rbx, -40(%rsp)
	xorl	%ebx, %ebx
	.p2align 4,,10
	.p2align 3
.L3:
	cmpq	%rbp, %rcx
	je	.L12
	leaq	16(%rdi), %r8
	cmpq	%r8, %rcx
	setnb	%r9b
	cmpq	-32(%rsp), %rdi
	setnb	%r10b
	orb	%r10b, %r9b
	je	.L18
	cmpb	$0, -17(%rsp)
	je	.L18
	movq	%rdi, %rax
	andl	$15, %eax
	shrq	$2, %rax
	negq	%rax
	andl	$3, %eax
	cmpq	%rax, %r12
	cmovbe	%r12, %rax
	testq	%rax, %rax
	je	.L16
	movss	(%rdi), %xmm1
	cmpq	$1, %rax
	leaq	4(%rdi), %rsi
	movaps	%xmm1, %xmm2
	movq	-16(%rsp), %r11
	subss	(%rcx), %xmm2
	mulss	%xmm0, %xmm2
	addss	%xmm2, %xmm1
	movss	%xmm1, (%rdi)
	je	.L6
	movss	4(%rdi), %xmm4
	cmpq	$3, %rax
	leaq	8(%rdi), %rsi
	movaps	%xmm4, %xmm5
	leaq	8(%rcx), %r11
	subss	4(%rcx), %xmm5
	mulss	%xmm0, %xmm5
	addss	%xmm5, %xmm4
	movss	%xmm4, 4(%rdi)
	jne	.L6
	movss	8(%rdi), %xmm6
	leaq	12(%rdi), %rsi
	leaq	12(%rcx), %r11
	movaps	%xmm6, %xmm7
	subss	8(%rcx), %xmm7
	mulss	%xmm0, %xmm7
	addss	%xmm7, %xmm6
	movss	%xmm6, 8(%rdi)
.L6:
	movq	%r12, %r14
	movl	$1, %r10d
	movl	$16, %edx
	subq	%rax, %r14
	salq	$2, %rax
	leaq	(%rdi,%rax), %r8
	addq	%rcx, %rax
	leaq	-4(%r14), %r9
	movups	(%rax), %xmm9
	shrq	$2, %r9
	movaps	(%r8), %xmm8
	leaq	1(%r9), %r13
	andl	$7, %r9d
	movaps	%xmm8, %xmm10
	cmpq	$1, %r13
	leaq	0(,%r13,4), %r15
	subps	%xmm9, %xmm10
	movaps	%xmm10, %xmm11
	mulps	%xmm3, %xmm11
	addps	%xmm11, %xmm8
	movaps	%xmm8, (%r8)
	jbe	.L97
	testq	%r9, %r9
	je	.L8
	cmpq	$1, %r9
	je	.L75
	cmpq	$2, %r9
	je	.L76
	cmpq	$3, %r9
	je	.L77
	cmpq	$4, %r9
	je	.L78
	cmpq	$5, %r9
	je	.L79
	cmpq	$6, %r9
	je	.L80
	movaps	16(%r8), %xmm12
	movl	$2, %r10d
	movl	$32, %edx
	movups	16(%rax), %xmm13
	movaps	%xmm12, %xmm14
	subps	%xmm13, %xmm14
	movaps	%xmm14, %xmm15
	mulps	%xmm3, %xmm15
	addps	%xmm15, %xmm12
	movaps	%xmm12, 16(%r8)
.L80:
	addq	$1, %r10
	movaps	(%r8,%rdx), %xmm1
	movups	(%rax,%rdx), %xmm2
	movaps	%xmm1, %xmm4
	subps	%xmm2, %xmm4
	movaps	%xmm4, %xmm5
	mulps	%xmm3, %xmm5
	addps	%xmm5, %xmm1
	movaps	%xmm1, (%r8,%rdx)
	addq	$16, %rdx
.L79:
	addq	$1, %r10
	movaps	(%r8,%rdx), %xmm6
	movups	(%rax,%rdx), %xmm8
	movaps	%xmm6, %xmm7
	subps	%xmm8, %xmm7
	movaps	%xmm7, %xmm9
	mulps	%xmm3, %xmm9
	addps	%xmm9, %xmm6
	movaps	%xmm6, (%r8,%rdx)
	addq	$16, %rdx
.L78:
	addq	$1, %r10
	movaps	(%r8,%rdx), %xmm10
	movups	(%rax,%rdx), %xmm11
	movaps	%xmm10, %xmm12
	subps	%xmm11, %xmm12
	movaps	%xmm12, %xmm13
	mulps	%xmm3, %xmm13
	addps	%xmm13, %xmm10
	movaps	%xmm10, (%r8,%rdx)
	addq	$16, %rdx
.L77:
	addq	$1, %r10
	movaps	(%r8,%rdx), %xmm14
	movups	(%rax,%rdx), %xmm15
	movaps	%xmm14, %xmm1
	subps	%xmm15, %xmm1
	movaps	%xmm1, %xmm2
	mulps	%xmm3, %xmm2
	addps	%xmm2, %xmm14
	movaps	%xmm14, (%r8,%rdx)
	addq	$16, %rdx
.L76:
	addq	$1, %r10
	movaps	(%r8,%rdx), %xmm4
	movups	(%rax,%rdx), %xmm5
	movaps	%xmm4, %xmm6
	subps	%xmm5, %xmm6
	movaps	%xmm6, %xmm8
	mulps	%xmm3, %xmm8
	addps	%xmm8, %xmm4
	movaps	%xmm4, (%r8,%rdx)
	addq	$16, %rdx
.L75:
	addq	$1, %r10
	movaps	(%r8,%rdx), %xmm7
	movups	(%rax,%rdx), %xmm9
	movaps	%xmm7, %xmm10
	subps	%xmm9, %xmm10
	movaps	%xmm10, %xmm11
	mulps	%xmm3, %xmm11
	addps	%xmm11, %xmm7
	movaps	%xmm7, (%r8,%rdx)
	addq	$16, %rdx
	cmpq	%r10, %r13
	jbe	.L97
.L8:
	addq	$8, %r10
	movaps	(%r8,%rdx), %xmm12
	movups	(%rax,%rdx), %xmm13
	movaps	%xmm12, %xmm14
	movaps	16(%r8,%rdx), %xmm1
	subps	%xmm13, %xmm14
	movaps	%xmm1, %xmm4
	movaps	32(%r8,%rdx), %xmm8
	movaps	%xmm8, %xmm6
	movaps	48(%r8,%rdx), %xmm10
	movaps	%xmm14, %xmm15
	movaps	64(%r8,%rdx), %xmm14
	mulps	%xmm3, %xmm15
	addps	%xmm15, %xmm12
	movaps	%xmm12, (%r8,%rdx)
	movaps	%xmm10, %xmm12
	movups	16(%rax,%rdx), %xmm2
	subps	%xmm2, %xmm4
	movaps	%xmm4, %xmm5
	movaps	80(%r8,%rdx), %xmm4
	mulps	%xmm3, %xmm5
	addps	%xmm5, %xmm1
	movaps	%xmm1, 16(%r8,%rdx)
	movaps	%xmm14, %xmm1
	movups	32(%rax,%rdx), %xmm7
	subps	%xmm7, %xmm6
	movaps	%xmm6, %xmm9
	mulps	%xmm3, %xmm9
	addps	%xmm9, %xmm8
	movaps	96(%r8,%rdx), %xmm9
	movaps	%xmm9, %xmm6
	movaps	%xmm8, 32(%r8,%rdx)
	movaps	%xmm4, %xmm8
	movups	48(%rax,%rdx), %xmm11
	subps	%xmm11, %xmm12
	movaps	%xmm12, %xmm13
	movaps	112(%r8,%rdx), %xmm12
	mulps	%xmm3, %xmm13
	addps	%xmm13, %xmm10
	movaps	%xmm10, 48(%r8,%rdx)
	movups	64(%rax,%rdx), %xmm15
	subps	%xmm15, %xmm1
	movaps	%xmm1, %xmm2
	mulps	%xmm3, %xmm2
	addps	%xmm2, %xmm14
	movaps	%xmm14, 64(%r8,%rdx)
	movaps	%xmm12, %xmm14
	movups	80(%rax,%rdx), %xmm5
	subps	%xmm5, %xmm8
	movaps	%xmm8, %xmm7
	mulps	%xmm3, %xmm7
	addps	%xmm7, %xmm4
	movaps	%xmm4, 80(%r8,%rdx)
	movups	96(%rax,%rdx), %xmm10
	subps	%xmm10, %xmm6
	movaps	%xmm6, %xmm11
	mulps	%xmm3, %xmm11
	addps	%xmm11, %xmm9
	movaps	%xmm9, 96(%r8,%rdx)
	movups	112(%rax,%rdx), %xmm13
	subps	%xmm13, %xmm14
	movaps	%xmm14, %xmm15
	mulps	%xmm3, %xmm15
	addps	%xmm15, %xmm12
	movaps	%xmm12, 112(%r8,%rdx)
	subq	$-128, %rdx
	cmpq	%r10, %r13
	ja	.L8
	.p2align 4,,10
	.p2align 3
.L97:
	leaq	0(,%r15,4), %rax
	addq	%rax, %rsi
	addq	%rax, %r11
	cmpq	%r15, %r14
	je	.L5
	movss	(%rsi), %xmm1
	leaq	4(%r11), %r9
	movaps	%xmm1, %xmm2
	cmpq	%r9, %rbp
	subss	(%r11), %xmm2
	mulss	%xmm0, %xmm2
	addss	%xmm2, %xmm1
	movss	%xmm1, (%rsi)
	je	.L5
	movss	4(%rsi), %xmm4
	leaq	8(%r11), %r13
	movaps	%xmm4, %xmm5
	cmpq	%r13, %rbp
	subss	4(%r11), %xmm5
	mulss	%xmm0, %xmm5
	addss	%xmm5, %xmm4
	movss	%xmm4, 4(%rsi)
	je	.L5
	movss	8(%rsi), %xmm8
	movaps	%xmm8, %xmm7
	subss	8(%r11), %xmm7
	mulss	%xmm0, %xmm7
	addss	%xmm7, %xmm8
	movss	%xmm8, 8(%rsi)
.L5:
	addq	-40(%rsp), %rdi
.L12:
	addq	$1, %rbx
	cmpq	-48(%rsp), %rbx
	jne	.L3
.L1:
	popq	%rbx
	.cfi_remember_state
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L16:
	.cfi_restore_state
	movq	%rcx, %r11
	movq	%rdi, %rsi
	jmp	.L6
	.p2align 4,,10
	.p2align 3
.L18:
	movss	(%rdi), %xmm9
	movq	%rbp, %r8
	movaps	%xmm9, %xmm10
	subq	%rcx, %r8
	subq	$4, %r8
	leaq	4(%rcx), %rdx
	leaq	4(%rdi), %r10
	subss	(%rcx), %xmm10
	shrq	$2, %r8
	andl	$7, %r8d
	cmpq	%rdx, %rbp
	mulss	%xmm0, %xmm10
	addss	%xmm10, %xmm9
	movss	%xmm9, (%rdi)
	je	.L5
	testq	%r8, %r8
	je	.L4
	cmpq	$1, %r8
	je	.L69
	cmpq	$2, %r8
	je	.L70
	cmpq	$3, %r8
	je	.L71
	cmpq	$4, %r8
	je	.L72
	cmpq	$5, %r8
	je	.L73
	cmpq	$6, %r8
	je	.L74
	movss	(%r10), %xmm6
	leaq	8(%rcx), %rdx
	movaps	%xmm6, %xmm11
	subss	4(%rcx), %xmm11
	mulss	%xmm0, %xmm11
	addss	%xmm11, %xmm6
	movss	%xmm6, (%r10)
	leaq	8(%rdi), %r10
.L74:
	movss	(%r10), %xmm12
	addq	$4, %rdx
	addq	$4, %r10
	movaps	%xmm12, %xmm13
	subss	-4(%rdx), %xmm13
	mulss	%xmm0, %xmm13
	addss	%xmm13, %xmm12
	movss	%xmm12, -4(%r10)
.L73:
	movss	(%r10), %xmm14
	addq	$4, %rdx
	addq	$4, %r10
	movaps	%xmm14, %xmm15
	subss	-4(%rdx), %xmm15
	mulss	%xmm0, %xmm15
	addss	%xmm15, %xmm14
	movss	%xmm14, -4(%r10)
.L72:
	movss	(%r10), %xmm1
	addq	$4, %rdx
	addq	$4, %r10
	movaps	%xmm1, %xmm2
	subss	-4(%rdx), %xmm2
	mulss	%xmm0, %xmm2
	addss	%xmm2, %xmm1
	movss	%xmm1, -4(%r10)
.L71:
	movss	(%r10), %xmm4
	addq	$4, %rdx
	addq	$4, %r10
	movaps	%xmm4, %xmm5
	subss	-4(%rdx), %xmm5
	mulss	%xmm0, %xmm5
	addss	%xmm5, %xmm4
	movss	%xmm4, -4(%r10)
.L70:
	movss	(%r10), %xmm8
	addq	$4, %rdx
	addq	$4, %r10
	movaps	%xmm8, %xmm7
	subss	-4(%rdx), %xmm7
	mulss	%xmm0, %xmm7
	addss	%xmm7, %xmm8
	movss	%xmm8, -4(%r10)
.L69:
	movss	(%r10), %xmm9
	addq	$4, %rdx
	addq	$4, %r10
	movaps	%xmm9, %xmm10
	subss	-4(%rdx), %xmm10
	mulss	%xmm0, %xmm10
	addss	%xmm10, %xmm9
	movss	%xmm9, -4(%r10)
	cmpq	%rdx, %rbp
	je	.L5
.L4:
	movss	(%r10), %xmm6
	addq	$32, %rdx
	addq	$32, %r10
	movaps	%xmm6, %xmm11
	movss	-28(%r10), %xmm12
	subss	-32(%rdx), %xmm11
	movaps	%xmm12, %xmm13
	movss	-24(%r10), %xmm14
	movaps	%xmm14, %xmm15
	movss	-20(%r10), %xmm1
	movaps	%xmm1, %xmm2
	movss	-16(%r10), %xmm4
	movaps	%xmm4, %xmm5
	movss	-12(%r10), %xmm8
	mulss	%xmm0, %xmm11
	movaps	%xmm8, %xmm7
	movss	-8(%r10), %xmm9
	movaps	%xmm9, %xmm10
	addss	%xmm11, %xmm6
	movss	%xmm6, -32(%r10)
	movss	-4(%r10), %xmm6
	subss	-28(%rdx), %xmm13
	movaps	%xmm6, %xmm11
	mulss	%xmm0, %xmm13
	addss	%xmm13, %xmm12
	movss	%xmm12, -28(%r10)
	subss	-24(%rdx), %xmm15
	mulss	%xmm0, %xmm15
	addss	%xmm15, %xmm14
	movss	%xmm14, -24(%r10)
	subss	-20(%rdx), %xmm2
	mulss	%xmm0, %xmm2
	addss	%xmm2, %xmm1
	movss	%xmm1, -20(%r10)
	subss	-16(%rdx), %xmm5
	mulss	%xmm0, %xmm5
	addss	%xmm5, %xmm4
	movss	%xmm4, -16(%r10)
	subss	-12(%rdx), %xmm7
	mulss	%xmm0, %xmm7
	addss	%xmm7, %xmm8
	movss	%xmm8, -12(%r10)
	subss	-8(%rdx), %xmm10
	mulss	%xmm0, %xmm10
	addss	%xmm10, %xmm9
	movss	%xmm9, -8(%r10)
	subss	-4(%rdx), %xmm11
	mulss	%xmm0, %xmm11
	addss	%xmm11, %xmm6
	movss	%xmm6, -4(%r10)
	cmpq	%rdx, %rbp
	jne	.L4
	jmp	.L5
	.cfi_endproc
.LFE0:
	.size	saasxpa_kh, .-saasxpa_kh
	.section	.text.unlikely
.LCOLDE0:
	.text
.LHOTE0:
	.ident	"GCC: (Debian 4.9.2-10) 4.9.2"
	.section	.note.GNU-stack,"",@progbits
