/*! \file laxpy_kh.S	
 * This file is part of the LEAC.
 *
 * Implementation of the 
 *	
 *	y <-- ax + y, (long <-- int * int + long)
 *
 *  function for int
 *
 *
 * (c)  Hermes Robles Berumen <hermes@uaz.edu.mx>
 *
 * For the full copyright and license information, please view the LICENSE
 * file that was distributed with this source code.
*/
	
	.file	"laxpy_kh.S"
	.section	.text.unlikely,"ax",@progbits
.LCOLDB0:
	.text
.LHOTB0:
	.p2align 4,,15
	.globl	laxpy_kh
	.type	laxpy_kh, @function
laxpy_kh:
.LFB0:
	.cfi_startproc
	cmpq	$1, %r9
	pushq	%r13
	.cfi_def_cfa_offset 16
	.cfi_offset 13, -16
	pushq	%r12
	.cfi_def_cfa_offset 24
	.cfi_offset 12, -24
	pushq	%rbp
	.cfi_def_cfa_offset 32
	.cfi_offset 6, -32
	pushq	%rbx
	.cfi_def_cfa_offset 40
	.cfi_offset 3, -40
	movq	40(%rsp), %rax
	jne	.L21
	cmpq	$1, 48(%rsp)
	jne	.L21
	testq	%rdi, %rdi
	jle	.L23
	movq	%rax, %rdx
	salq	$60, %rdx
	shrq	$63, %rdx
	cmpq	%rdi, %rdx
	cmova	%rdi, %rdx
	cmpq	$4, %rdi
	jg	.L85
	movq	%rdi, %rdx
.L16:
	movl	(%r8), %esi
	imull	%ecx, %esi
	movslq	%esi, %rbx
	addq	%rbx, (%rax)
	cmpq	$1, %rdx
	movl	$1, %esi
	je	.L8
	movl	4(%r8), %ebp
	movl	$2, %esi
	imull	%ecx, %ebp
	movslq	%ebp, %r9
	addq	%r9, 8(%rax)
	cmpq	$2, %rdx
	je	.L8
	movl	8(%r8), %r10d
	movl	$3, %esi
	imull	%ecx, %r10d
	movslq	%r10d, %r11
	addq	%r11, 16(%rax)
	cmpq	$4, %rdx
	jne	.L8
	movl	12(%r8), %r12d
	movl	$4, %esi
	imull	%ecx, %r12d
	movslq	%r12d, %r13
	addq	%r13, 24(%rax)
.L8:
	cmpq	%rdi, %rdx
	je	.L23
.L7:
	movq	%rdi, %rbx
	leaq	-1(%rdi), %r9
	subq	%rdx, %rbx
	leaq	-4(%rbx), %r10
	subq	%rdx, %r9
	shrq	$2, %r10
	cmpq	$2, %r9
	leaq	1(%r10), %r11
	leaq	0(,%r11,4), %rbp
	jbe	.L10
	movl	%ecx, -12(%rsp)
	leaq	(%r8,%rdx,4), %r13
	leaq	(%rax,%rdx,8), %r9
	movd	-12(%rsp), %xmm6
	pxor	%xmm3, %xmm3
	movdqu	0(%r13), %xmm1
	andl	$3, %r10d
	pshufd	$0, %xmm6, %xmm2
	cmpq	$1, %r11
	movdqa	%xmm1, %xmm0
	psrlq	$32, %xmm1
	movl	$1, %r12d
	movdqa	%xmm3, %xmm8
	movl	$16, %edx
	movdqa	%xmm2, %xmm4
	pmuludq	%xmm2, %xmm0
	pshufd	$8, %xmm0, %xmm7
	psrlq	$32, %xmm4
	pmuludq	%xmm4, %xmm1
	pshufd	$8, %xmm1, %xmm5
	punpckldq	%xmm5, %xmm7
	pcmpgtd	%xmm7, %xmm8
	movdqa	%xmm7, %xmm9
	punpckldq	%xmm8, %xmm7
	punpckhdq	%xmm8, %xmm9
	paddq	(%r9), %xmm7
	movaps	%xmm7, (%r9)
	paddq	16(%r9), %xmm9
	movaps	%xmm9, 16(%r9)
	jbe	.L82
	testq	%r10, %r10
	je	.L11
	cmpq	$1, %r10
	je	.L67
	cmpq	$2, %r10
	je	.L68
	movdqu	16(%r13), %xmm10
	movdqa	%xmm3, %xmm14
	movl	$2, %r12d
	movl	$32, %edx
	movdqa	%xmm10, %xmm11
	psrlq	$32, %xmm10
	pmuludq	%xmm4, %xmm10
	pshufd	$8, %xmm10, %xmm13
	pmuludq	%xmm2, %xmm11
	pshufd	$8, %xmm11, %xmm12
	punpckldq	%xmm13, %xmm12
	pcmpgtd	%xmm12, %xmm14
	movdqa	%xmm12, %xmm15
	punpckldq	%xmm14, %xmm12
	punpckhdq	%xmm14, %xmm15
	paddq	32(%r9), %xmm12
	movaps	%xmm12, 32(%r9)
	paddq	48(%r9), %xmm15
	movaps	%xmm15, 48(%r9)
.L68:
	movdqa	%xmm3, %xmm5
	addq	$1, %r12
	movdqu	0(%r13,%rdx), %xmm6
	movdqa	%xmm6, %xmm1
	psrlq	$32, %xmm6
	pmuludq	%xmm4, %xmm6
	pshufd	$8, %xmm6, %xmm7
	pmuludq	%xmm2, %xmm1
	pshufd	$8, %xmm1, %xmm0
	punpckldq	%xmm7, %xmm0
	pcmpgtd	%xmm0, %xmm5
	movdqa	%xmm0, %xmm8
	punpckldq	%xmm5, %xmm0
	punpckhdq	%xmm5, %xmm8
	paddq	(%r9,%rdx,2), %xmm0
	movaps	%xmm0, (%r9,%rdx,2)
	paddq	16(%r9,%rdx,2), %xmm8
	movaps	%xmm8, 16(%r9,%rdx,2)
	addq	$16, %rdx
.L67:
	movdqa	%xmm3, %xmm13
	addq	$1, %r12
	movdqu	0(%r13,%rdx), %xmm9
	movdqa	%xmm9, %xmm10
	psrlq	$32, %xmm9
	pmuludq	%xmm4, %xmm9
	pshufd	$8, %xmm9, %xmm12
	pmuludq	%xmm2, %xmm10
	pshufd	$8, %xmm10, %xmm11
	punpckldq	%xmm12, %xmm11
	pcmpgtd	%xmm11, %xmm13
	movdqa	%xmm11, %xmm14
	punpckldq	%xmm13, %xmm11
	punpckhdq	%xmm13, %xmm14
	paddq	(%r9,%rdx,2), %xmm11
	movaps	%xmm11, (%r9,%rdx,2)
	paddq	16(%r9,%rdx,2), %xmm14
	movaps	%xmm14, 16(%r9,%rdx,2)
	addq	$16, %rdx
	cmpq	%r12, %r11
	jbe	.L82
.L11:
	movdqa	%xmm3, %xmm7
	movdqa	%xmm3, %xmm12
	movdqu	0(%r13,%rdx), %xmm15
	leaq	16(%rdx), %r10
	addq	$4, %r12
	movdqa	%xmm15, %xmm6
	movdqu	16(%r13,%rdx), %xmm8
	psrlq	$32, %xmm15
	pmuludq	%xmm4, %xmm15
	pshufd	$8, %xmm15, %xmm1
	pmuludq	%xmm2, %xmm6
	pshufd	$8, %xmm6, %xmm0
	movdqa	%xmm8, %xmm9
	psrlq	$32, %xmm8
	pmuludq	%xmm4, %xmm8
	pshufd	$8, %xmm8, %xmm11
	pmuludq	%xmm2, %xmm9
	pshufd	$8, %xmm9, %xmm10
	movdqu	32(%r13,%rdx), %xmm14
	punpckldq	%xmm1, %xmm0
	movdqa	%xmm14, %xmm15
	psrlq	$32, %xmm14
	pmuludq	%xmm4, %xmm14
	punpckldq	%xmm11, %xmm10
	movdqa	%xmm3, %xmm11
	pcmpgtd	%xmm0, %xmm7
	movdqa	%xmm0, %xmm5
	pmuludq	%xmm2, %xmm15
	pshufd	$8, %xmm15, %xmm6
	punpckhdq	%xmm7, %xmm5
	pcmpgtd	%xmm10, %xmm12
	punpckldq	%xmm7, %xmm0
	movdqa	%xmm10, %xmm13
	punpckldq	%xmm12, %xmm10
	paddq	16(%r9,%rdx,2), %xmm5
	movaps	%xmm5, 16(%r9,%rdx,2)
	paddq	(%r9,%rdx,2), %xmm0
	movdqa	%xmm3, %xmm7
	movaps	%xmm0, (%r9,%rdx,2)
	pshufd	$8, %xmm14, %xmm0
	punpckhdq	%xmm12, %xmm13
	paddq	(%r9,%r10,2), %xmm10
	movaps	%xmm10, (%r9,%r10,2)
	punpckldq	%xmm0, %xmm6
	paddq	16(%r9,%r10,2), %xmm13
	movaps	%xmm13, 16(%r9,%r10,2)
	leaq	32(%rdx), %r10
	pcmpgtd	%xmm6, %xmm7
	movdqa	%xmm6, %xmm1
	movdqu	48(%r13,%rdx), %xmm5
	punpckhdq	%xmm7, %xmm1
	movdqa	%xmm5, %xmm8
	psrlq	$32, %xmm5
	pmuludq	%xmm4, %xmm5
	pshufd	$8, %xmm5, %xmm10
	pmuludq	%xmm2, %xmm8
	pshufd	$8, %xmm8, %xmm9
	paddq	16(%r9,%r10,2), %xmm1
	punpckldq	%xmm7, %xmm6
	movaps	%xmm1, 16(%r9,%r10,2)
	punpckldq	%xmm10, %xmm9
	paddq	(%r9,%r10,2), %xmm6
	movaps	%xmm6, (%r9,%r10,2)
	leaq	48(%rdx), %r10
	addq	$64, %rdx
	cmpq	%r12, %r11
	pcmpgtd	%xmm9, %xmm11
	movdqa	%xmm9, %xmm12
	punpckldq	%xmm11, %xmm9
	punpckhdq	%xmm11, %xmm12
	paddq	(%r9,%r10,2), %xmm9
	movaps	%xmm9, (%r9,%r10,2)
	paddq	16(%r9,%r10,2), %xmm12
	movaps	%xmm12, 16(%r9,%r10,2)
	ja	.L11
	.p2align 4,,10
	.p2align 3
.L82:
	addq	%rbp, %rsi
	cmpq	%rbx, %rbp
	je	.L23
.L10:
	movl	(%r8,%rsi,4), %r11d
	leaq	1(%rsi), %r13
	imull	%ecx, %r11d
	movslq	%r11d, %r9
	addq	%r9, (%rax,%rsi,8)
	cmpq	%r13, %rdi
	jle	.L23
	movl	(%r8,%r13,4), %r12d
	addq	$2, %rsi
	imull	%ecx, %r12d
	movslq	%r12d, %rdx
	addq	%rdx, (%rax,%r13,8)
	cmpq	%rsi, %rdi
	jle	.L23
	imull	(%r8,%rsi,4), %ecx
	movslq	%ecx, %rcx
	addq	%rcx, (%rax,%rsi,8)
.L23:
	popq	%rbx
	.cfi_remember_state
	.cfi_def_cfa_offset 32
	movl	$1, %eax
	popq	%rbp
	.cfi_def_cfa_offset 24
	popq	%r12
	.cfi_def_cfa_offset 16
	popq	%r13
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L21:
	.cfi_restore_state
	testq	%rdi, %rdi
	jle	.L23
	movl	(%r8), %esi
	movq	48(%rsp), %rbx
	leaq	-1(%rdi), %r10
	salq	$2, %r9
	movl	$1, %edx
	andl	$7, %r10d
	addq	%r9, %r8
	imull	%ecx, %esi
	leaq	0(,%rbx,8), %rbp
	movslq	%esi, %r11
	addq	%r11, (%rax)
	addq	%rbp, %rax
	cmpq	%rdi, %rdx
	je	.L23
	testq	%r10, %r10
	je	.L15
	cmpq	$1, %r10
	je	.L61
	cmpq	$2, %r10
	je	.L62
	cmpq	$3, %r10
	je	.L63
	cmpq	$4, %r10
	je	.L64
	cmpq	$5, %r10
	je	.L65
	cmpq	$6, %r10
	je	.L66
	movl	(%r8), %r13d
	movl	$2, %edx
	addq	%r9, %r8
	imull	%ecx, %r13d
	movslq	%r13d, %r12
	addq	%r12, (%rax)
	addq	%rbp, %rax
.L66:
	movl	(%r8), %ebx
	addq	$1, %rdx
	addq	%r9, %r8
	imull	%ecx, %ebx
	movslq	%ebx, %r10
	addq	%r10, (%rax)
	addq	%rbp, %rax
.L65:
	movl	(%r8), %esi
	addq	$1, %rdx
	addq	%r9, %r8
	imull	%ecx, %esi
	movslq	%esi, %r11
	addq	%r11, (%rax)
	addq	%rbp, %rax
.L64:
	movl	(%r8), %r13d
	addq	$1, %rdx
	addq	%r9, %r8
	imull	%ecx, %r13d
	movslq	%r13d, %r12
	addq	%r12, (%rax)
	addq	%rbp, %rax
.L63:
	movl	(%r8), %ebx
	addq	$1, %rdx
	addq	%r9, %r8
	imull	%ecx, %ebx
	movslq	%ebx, %r10
	addq	%r10, (%rax)
	addq	%rbp, %rax
.L62:
	movl	(%r8), %esi
	addq	$1, %rdx
	addq	%r9, %r8
	imull	%ecx, %esi
	movslq	%esi, %r11
	addq	%r11, (%rax)
	addq	%rbp, %rax
.L61:
	movl	(%r8), %r13d
	addq	$1, %rdx
	addq	%r9, %r8
	imull	%ecx, %r13d
	movslq	%r13d, %r12
	addq	%r12, (%rax)
	addq	%rbp, %rax
	cmpq	%rdi, %rdx
	je	.L23
.L15:
	movl	(%r8), %ebx
	addq	%r9, %r8
	addq	$8, %rdx
	movl	(%r8), %esi
	addq	%r9, %r8
	movl	(%r8), %r13d
	addq	%r9, %r8
	imull	%ecx, %ebx
	imull	%ecx, %esi
	movslq	%ebx, %r10
	movl	(%r8), %ebx
	addq	%r9, %r8
	movslq	%esi, %r11
	movl	(%r8), %esi
	addq	%r9, %r8
	imull	%ecx, %r13d
	addq	%r10, (%rax)
	addq	%rbp, %rax
	imull	%ecx, %ebx
	addq	%r11, (%rax)
	addq	%rbp, %rax
	movslq	%r13d, %r12
	movl	(%r8), %r13d
	addq	%r9, %r8
	imull	%ecx, %esi
	movslq	%ebx, %r10
	movl	(%r8), %ebx
	addq	%r9, %r8
	addq	%r12, (%rax)
	addq	%rbp, %rax
	movslq	%esi, %r11
	movl	(%r8), %esi
	addq	%r10, (%rax)
	imull	%ecx, %r13d
	addq	%rbp, %rax
	addq	%r9, %r8
	imull	%ecx, %ebx
	addq	%r11, (%rax)
	addq	%rbp, %rax
	movslq	%r13d, %r12
	addq	%r12, (%rax)
	addq	%rbp, %rax
	imull	%ecx, %esi
	movslq	%ebx, %r10
	addq	%r10, (%rax)
	addq	%rbp, %rax
	movslq	%esi, %r11
	addq	%r11, (%rax)
	addq	%rbp, %rax
	cmpq	%rdi, %rdx
	jne	.L15
	jmp	.L23
	.p2align 4,,10
	.p2align 3
.L85:
	testq	%rdx, %rdx
	jne	.L16
	xorl	%esi, %esi
	jmp	.L7
	.cfi_endproc
.LFE0:
	.size	laxpy_kh, .-laxpy_kh
	.section	.text.unlikely
.LCOLDE0:
	.text
.LHOTE0:
	.ident	"GCC: (Debian 4.9.2-10) 4.9.2"
	.section	.note.GNU-stack,"",@progbits
