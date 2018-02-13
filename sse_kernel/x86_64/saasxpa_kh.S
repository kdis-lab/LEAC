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
	.p2align 5,,31
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
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	leaq	(%r9,%rdi,4), %rbp
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	movq	%rsi, -48(%rsp)
	testq	%rsi, %rsi
	jle	.L2
	leaq	4(%r9), %rax
	movq	%rbp, %r12
	movaps	%xmm0, %xmm3
	leaq	16(%r9), %rsi
	subq	%rax, %r12
	movq	%rax, -16(%rsp)
	shufps	$0, %xmm3, %xmm3
	shrq	$2, %r12
	movq	%rsi, -32(%rsp)
	addq	$1, %r12
	leaq	0(,%r12,4), %rdx
	cmpq	$6, %r12
	seta	-17(%rsp)
	movq	%rdx, -40(%rsp)
	xorl	%ebx, %ebx
	.p2align 5
.L3:
	cmpq	%rbp, %r9
	je	.L12
	leaq	16(%rcx), %rdi
	cmpq	%rdi, %r9
	setnb	%r8b
	cmpq	-32(%rsp), %rcx
	setnb	%r10b
	orb	%r10b, %r8b
	je	.L18
	cmpb	$0, -17(%rsp)
	je	.L18
	movq	%rcx, %rax
	andl	$15, %eax
	shrq	$2, %rax
	negq	%rax
	andl	$3, %eax
	cmpq	%r12, %rax
	cmova	%r12, %rax
	testq	%rax, %rax
	je	.L16
	movss	(%rcx), %xmm1
	leaq	4(%rcx), %rsi
	movq	-16(%rsp), %r11
	movaps	%xmm1, %xmm2
	subss	(%r9), %xmm2
	mulss	%xmm0, %xmm2
	addss	%xmm2, %xmm1
	movss	%xmm1, (%rcx)
	cmpq	$1, %rax
	je	.L6
	movss	4(%rcx), %xmm4
	leaq	8(%rcx), %rsi
	leaq	8(%r9), %r11
	movaps	%xmm4, %xmm5
	subss	4(%r9), %xmm5
	mulss	%xmm0, %xmm5
	addss	%xmm5, %xmm4
	movss	%xmm4, 4(%rcx)
	cmpq	$3, %rax
	jne	.L6
	movss	8(%rcx), %xmm6
	leaq	12(%rcx), %rsi
	leaq	12(%r9), %r11
	movaps	%xmm6, %xmm7
	subss	8(%r9), %xmm7
	mulss	%xmm0, %xmm7
	addss	%xmm7, %xmm6
	movss	%xmm6, 8(%rcx)
.L6:
	movq	%r12, %r14
	movl	$1, %r10d
	movl	$16, %edx
	subq	%rax, %r14
	salq	$2, %rax
	leaq	(%rcx,%rax), %rdi
	addq	%r9, %rax
	movaps	(%rdi), %xmm8
	leaq	-4(%r14), %r8
	movups	(%rax), %xmm9
	shrq	$2, %r8
	leaq	1(%r8), %r13
	movaps	%xmm8, %xmm10
	andl	$7, %r8d
	subps	%xmm9, %xmm10
	leaq	0(,%r13,4), %r15
	movaps	%xmm10, %xmm11
	mulps	%xmm3, %xmm11
	addps	%xmm11, %xmm8
	movaps	%xmm8, (%rdi)
	cmpq	$1, %r13
	jbe	.L97
	testq	%r8, %r8
	je	.L8
	cmpq	$1, %r8
	je	.L75
	cmpq	$2, %r8
	je	.L76
	cmpq	$3, %r8
	je	.L77
	cmpq	$4, %r8
	je	.L78
	cmpq	$5, %r8
	je	.L79
	cmpq	$6, %r8
	je	.L80
	movaps	16(%rdi), %xmm12
	movl	$2, %r10d
	movl	$32, %edx
	movups	16(%rax), %xmm13
	movaps	%xmm12, %xmm14
	subps	%xmm13, %xmm14
	movaps	%xmm14, %xmm15
	mulps	%xmm3, %xmm15
	addps	%xmm15, %xmm12
	movaps	%xmm12, 16(%rdi)
.L80:
	movaps	(%rdi,%rdx), %xmm1
	addq	$1, %r10
	movups	(%rax,%rdx), %xmm2
	movaps	%xmm1, %xmm4
	subps	%xmm2, %xmm4
	movaps	%xmm4, %xmm5
	mulps	%xmm3, %xmm5
	addps	%xmm5, %xmm1
	movaps	%xmm1, (%rdi,%rdx)
	addq	$16, %rdx
.L79:
	movaps	(%rdi,%rdx), %xmm6
	addq	$1, %r10
	movups	(%rax,%rdx), %xmm8
	movaps	%xmm6, %xmm7
	subps	%xmm8, %xmm7
	movaps	%xmm7, %xmm9
	mulps	%xmm3, %xmm9
	addps	%xmm9, %xmm6
	movaps	%xmm6, (%rdi,%rdx)
	addq	$16, %rdx
.L78:
	movaps	(%rdi,%rdx), %xmm10
	addq	$1, %r10
	movups	(%rax,%rdx), %xmm11
	movaps	%xmm10, %xmm12
	subps	%xmm11, %xmm12
	movaps	%xmm12, %xmm13
	mulps	%xmm3, %xmm13
	addps	%xmm13, %xmm10
	movaps	%xmm10, (%rdi,%rdx)
	addq	$16, %rdx
.L77:
	movaps	(%rdi,%rdx), %xmm14
	addq	$1, %r10
	movups	(%rax,%rdx), %xmm15
	movaps	%xmm14, %xmm1
	subps	%xmm15, %xmm1
	movaps	%xmm1, %xmm2
	mulps	%xmm3, %xmm2
	addps	%xmm2, %xmm14
	movaps	%xmm14, (%rdi,%rdx)
	addq	$16, %rdx
.L76:
	movaps	(%rdi,%rdx), %xmm4
	addq	$1, %r10
	movups	(%rax,%rdx), %xmm5
	movaps	%xmm4, %xmm6
	subps	%xmm5, %xmm6
	movaps	%xmm6, %xmm8
	mulps	%xmm3, %xmm8
	addps	%xmm8, %xmm4
	movaps	%xmm4, (%rdi,%rdx)
	addq	$16, %rdx
.L75:
	movaps	(%rdi,%rdx), %xmm7
	addq	$1, %r10
	movups	(%rax,%rdx), %xmm9
	movaps	%xmm7, %xmm10
	subps	%xmm9, %xmm10
	movaps	%xmm10, %xmm11
	mulps	%xmm3, %xmm11
	addps	%xmm11, %xmm7
	movaps	%xmm7, (%rdi,%rdx)
	addq	$16, %rdx
	cmpq	%r10, %r13
	jbe	.L97
.L8:
	movaps	(%rdi,%rdx), %xmm12
	addq	$8, %r10
	movups	(%rax,%rdx), %xmm13
	movaps	%xmm12, %xmm14
	movaps	16(%rdi,%rdx), %xmm1
	subps	%xmm13, %xmm14
	movaps	32(%rdi,%rdx), %xmm8
	movaps	%xmm1, %xmm4
	movaps	48(%rdi,%rdx), %xmm10
	movaps	%xmm8, %xmm6
	movaps	%xmm14, %xmm15
	movaps	64(%rdi,%rdx), %xmm14
	mulps	%xmm3, %xmm15
	addps	%xmm15, %xmm12
	movaps	%xmm12, (%rdi,%rdx)
	movups	16(%rax,%rdx), %xmm2
	movaps	%xmm10, %xmm12
	subps	%xmm2, %xmm4
	movaps	%xmm4, %xmm5
	movaps	80(%rdi,%rdx), %xmm4
	mulps	%xmm3, %xmm5
	addps	%xmm5, %xmm1
	movaps	%xmm1, 16(%rdi,%rdx)
	movups	32(%rax,%rdx), %xmm7
	movaps	%xmm14, %xmm1
	subps	%xmm7, %xmm6
	movaps	%xmm6, %xmm9
	mulps	%xmm3, %xmm9
	addps	%xmm9, %xmm8
	movaps	96(%rdi,%rdx), %xmm9
	movaps	%xmm9, %xmm6
	movaps	%xmm8, 32(%rdi,%rdx)
	movups	48(%rax,%rdx), %xmm11
	movaps	%xmm4, %xmm8
	subps	%xmm11, %xmm12
	movaps	%xmm12, %xmm13
	movaps	112(%rdi,%rdx), %xmm12
	mulps	%xmm3, %xmm13
	addps	%xmm13, %xmm10
	movaps	%xmm10, 48(%rdi,%rdx)
	movups	64(%rax,%rdx), %xmm15
	subps	%xmm15, %xmm1
	movaps	%xmm1, %xmm2
	mulps	%xmm3, %xmm2
	addps	%xmm2, %xmm14
	movaps	%xmm14, 64(%rdi,%rdx)
	movups	80(%rax,%rdx), %xmm5
	movaps	%xmm12, %xmm14
	subps	%xmm5, %xmm8
	movaps	%xmm8, %xmm7
	mulps	%xmm3, %xmm7
	addps	%xmm7, %xmm4
	movaps	%xmm4, 80(%rdi,%rdx)
	movups	96(%rax,%rdx), %xmm10
	subps	%xmm10, %xmm6
	movaps	%xmm6, %xmm11
	mulps	%xmm3, %xmm11
	addps	%xmm11, %xmm9
	movaps	%xmm9, 96(%rdi,%rdx)
	movups	112(%rax,%rdx), %xmm13
	subps	%xmm13, %xmm14
	movaps	%xmm14, %xmm15
	mulps	%xmm3, %xmm15
	addps	%xmm15, %xmm12
	movaps	%xmm12, 112(%rdi,%rdx)
	subq	$-128, %rdx
	cmpq	%r10, %r13
	ja	.L8
	.p2align 5
.L97:
	leaq	0(,%r15,4), %rax
	addq	%rax, %rsi
	addq	%rax, %r11
	cmpq	%r15, %r14
	je	.L5
	movss	(%rsi), %xmm1
	leaq	4(%r11), %r8
	movaps	%xmm1, %xmm2
	subss	(%r11), %xmm2
	mulss	%xmm0, %xmm2
	addss	%xmm2, %xmm1
	movss	%xmm1, (%rsi)
	cmpq	%r8, %rbp
	je	.L5
	movss	4(%rsi), %xmm4
	leaq	8(%r11), %r13
	movaps	%xmm4, %xmm5
	subss	4(%r11), %xmm5
	mulss	%xmm0, %xmm5
	addss	%xmm5, %xmm4
	movss	%xmm4, 4(%rsi)
	cmpq	%r13, %rbp
	je	.L5
	movss	8(%rsi), %xmm8
	movaps	%xmm8, %xmm7
	subss	8(%r11), %xmm7
	mulss	%xmm0, %xmm7
	addss	%xmm7, %xmm8
	movss	%xmm8, 8(%rsi)
.L5:
	addq	-40(%rsp), %rcx
.L12:
	addq	$1, %rbx
	cmpq	%rbx, -48(%rsp)
	jne	.L3
.L2:
	popq	%rbx
	.cfi_remember_state
	.cfi_def_cfa_offset 48
	movl	$1, %eax
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
	movq	%r9, %r11
	movq	%rcx, %rsi
	jmp	.L6
	.p2align 4,,10
	.p2align 3
.L18:
	movss	(%rcx), %xmm9
	leaq	4(%r9), %rdx
	movq	%rbp, %rdi
	leaq	4(%rcx), %r10
	subq	%r9, %rdi
	movaps	%xmm9, %xmm10
	subss	(%r9), %xmm10
	subq	$4, %rdi
	shrq	$2, %rdi
	andl	$7, %edi
	mulss	%xmm0, %xmm10
	addss	%xmm10, %xmm9
	movss	%xmm9, (%rcx)
	cmpq	%rdx, %rbp
	je	.L5
	testq	%rdi, %rdi
	je	.L4
	cmpq	$1, %rdi
	je	.L69
	cmpq	$2, %rdi
	je	.L70
	cmpq	$3, %rdi
	je	.L71
	cmpq	$4, %rdi
	je	.L72
	cmpq	$5, %rdi
	je	.L73
	cmpq	$6, %rdi
	je	.L74
	movss	(%r10), %xmm6
	leaq	8(%r9), %rdx
	movaps	%xmm6, %xmm11
	subss	4(%r9), %xmm11
	mulss	%xmm0, %xmm11
	addss	%xmm11, %xmm6
	movss	%xmm6, (%r10)
	leaq	8(%rcx), %r10
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
	movss	-28(%r10), %xmm12
	movaps	%xmm6, %xmm11
	subss	-32(%rdx), %xmm11
	movaps	%xmm12, %xmm13
	movss	-24(%r10), %xmm14
	movss	-20(%r10), %xmm1
	mulss	%xmm0, %xmm11
	movaps	%xmm14, %xmm15
	movss	-16(%r10), %xmm4
	movaps	%xmm1, %xmm2
	movss	-12(%r10), %xmm8
	movaps	%xmm4, %xmm5
	movss	-8(%r10), %xmm9
	movaps	%xmm8, %xmm7
	addss	%xmm11, %xmm6
	movaps	%xmm9, %xmm10
	movss	%xmm6, -32(%r10)
	subss	-28(%rdx), %xmm13
	movss	-4(%r10), %xmm6
	mulss	%xmm0, %xmm13
	movaps	%xmm6, %xmm11
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
	.ident	"GCC: (GNU) 5.5.0"
	.section	.note.GNU-stack,"",@progbits
