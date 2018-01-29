/*! \file inrm2_k.S	
 * This file is part of the LEAC.
 *
 * Implementation of the 
 *	
 *	nom2 <-- ||x-y||2, (unsigned long int <-- ||int -int||2)
 *
 *  function for int
 *
 *
 * (c)  Hermes Robles Berumen <hermes@uaz.edu.mx>
 *
 * For the full copyright and license information, please view the LICENSE
 * file that was distributed with this source code.
 */

	.file	"inrm2_k.S"
	.section	.text.unlikely,"ax",@progbits
.LCOLDB0:
	.text
.LHOTB0:
	.p2align 4,,15
	.globl	inrm2_k
	.type	inrm2_k, @function
inrm2_k:
.LFB0:
	.cfi_startproc
	cmpq	$1, %rdx
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
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	jne	.L30
	cmpq	$1, %r8
	jne	.L30
	testq	%rdi, %rdi
	jle	.L22
	movq	%rsi, %rdx
	andl	$15, %edx
	shrq	$2, %rdx
	negq	%rdx
	andl	$3, %edx
	cmpq	%rdi, %rdx
	cmova	%rdi, %rdx
	cmpq	$5, %rdi
	jg	.L56
	movq	%rdi, %rdx
.L6:
	movl	(%rsi), %eax
	subl	(%rcx), %eax
	cltq
	imulq	%rax, %rax
	cmpq	$1, %rdx
	je	.L25
	movl	4(%rsi), %r8d
	subl	4(%rcx), %r8d
	movslq	%r8d, %rbx
	imulq	%rbx, %rbx
	addq	%rbx, %rax
	cmpq	$2, %rdx
	je	.L26
	movl	8(%rsi), %ebp
	subl	8(%rcx), %ebp
	movslq	%ebp, %r9
	imulq	%r9, %r9
	addq	%r9, %rax
	cmpq	$3, %rdx
	je	.L27
	movl	12(%rsi), %r10d
	subl	12(%rcx), %r10d
	movslq	%r10d, %r11
	imulq	%r11, %r11
	addq	%r11, %rax
	cmpq	$5, %rdx
	jne	.L28
	movl	16(%rsi), %r12d
	subl	16(%rcx), %r12d
	movl	$5, %r8d
	movslq	%r12d, %r13
	imulq	%r13, %r13
	addq	%r13, %rax
.L8:
	cmpq	%rdx, %rdi
	je	.L5
.L7:
	movq	%rdi, %r11
	leaq	-1(%rdi), %rbp
	subq	%rdx, %r11
	leaq	-4(%r11), %r14
	subq	%rdx, %rbp
	shrq	$2, %r14
	cmpq	$2, %rbp
	leaq	1(%r14), %r15
	leaq	0(,%r15,4), %rbx
	jbe	.L10
	salq	$2, %rdx
	pxor	%xmm3, %xmm3
	leaq	(%rsi,%rdx), %r9
	addq	%rcx, %rdx
	andl	$1, %r14d
	movdqu	(%rdx), %xmm1
	cmpq	$1, %r15
	movdqa	%xmm3, %xmm2
	movl	$1, %r10d
	movdqa	(%r9), %xmm0
	movl	$16, %r12d
	psubd	%xmm1, %xmm0
	pcmpgtd	%xmm0, %xmm2
	movdqa	%xmm0, %xmm4
	punpckhdq	%xmm2, %xmm0
	punpckldq	%xmm2, %xmm4
	movdqa	%xmm0, %xmm7
	movdqa	%xmm4, %xmm5
	movdqa	%xmm4, %xmm6
	movdqa	%xmm0, %xmm8
	pmuludq	%xmm0, %xmm7
	pmuludq	%xmm4, %xmm5
	psrlq	$32, %xmm6
	pmuludq	%xmm6, %xmm4
	psrlq	$32, %xmm8
	paddq	%xmm4, %xmm4
	pmuludq	%xmm8, %xmm0
	psllq	$32, %xmm4
	paddq	%xmm0, %xmm0
	paddq	%xmm4, %xmm5
	psllq	$32, %xmm0
	paddq	%xmm7, %xmm0
	paddq	%xmm5, %xmm0
	jbe	.L48
	testq	%r14, %r14
	je	.L11
	movdqu	16(%rdx), %xmm9
	movdqa	%xmm3, %xmm11
	cmpq	$2, %r15
	movl	$2, %r10d
	movl	$32, %r12d
	movdqa	16(%r9), %xmm10
	psubd	%xmm9, %xmm10
	pcmpgtd	%xmm10, %xmm11
	movdqa	%xmm10, %xmm12
	punpckhdq	%xmm11, %xmm10
	punpckldq	%xmm11, %xmm12
	movdqa	%xmm10, %xmm15
	movdqa	%xmm12, %xmm13
	movdqa	%xmm12, %xmm14
	movdqa	%xmm10, %xmm1
	pmuludq	%xmm10, %xmm15
	pmuludq	%xmm12, %xmm13
	psrlq	$32, %xmm14
	pmuludq	%xmm14, %xmm12
	psrlq	$32, %xmm1
	paddq	%xmm12, %xmm12
	pmuludq	%xmm1, %xmm10
	psllq	$32, %xmm12
	paddq	%xmm10, %xmm10
	paddq	%xmm12, %xmm13
	psllq	$32, %xmm10
	paddq	%xmm13, %xmm0
	paddq	%xmm15, %xmm10
	paddq	%xmm10, %xmm0
	jbe	.L48
.L11:
	movdqu	(%rdx,%r12), %xmm2
	movdqa	%xmm3, %xmm4
	movdqa	%xmm3, %xmm14
	addq	$2, %r10
	movdqa	(%r9,%r12), %xmm6
	psubd	%xmm2, %xmm6
	pcmpgtd	%xmm6, %xmm4
	movdqa	%xmm6, %xmm7
	movdqa	16(%r9,%r12), %xmm13
	punpckldq	%xmm4, %xmm7
	punpckhdq	%xmm4, %xmm6
	movdqa	%xmm7, %xmm5
	movdqa	%xmm7, %xmm8
	movdqa	%xmm6, %xmm10
	movdqa	%xmm6, %xmm11
	pmuludq	%xmm7, %xmm5
	pmuludq	%xmm6, %xmm10
	psrlq	$32, %xmm8
	pmuludq	%xmm8, %xmm7
	psrlq	$32, %xmm11
	paddq	%xmm7, %xmm7
	pmuludq	%xmm11, %xmm6
	movdqa	%xmm7, %xmm9
	paddq	%xmm6, %xmm6
	psllq	$32, %xmm6
	paddq	%xmm10, %xmm6
	movdqa	%xmm6, %xmm12
	psllq	$32, %xmm9
	paddq	%xmm5, %xmm9
	paddq	%xmm9, %xmm0
	paddq	%xmm0, %xmm12
	movdqu	16(%rdx,%r12), %xmm0
	addq	$32, %r12
	cmpq	%r10, %r15
	psubd	%xmm0, %xmm13
	pcmpgtd	%xmm13, %xmm14
	movdqa	%xmm13, %xmm15
	punpckhdq	%xmm14, %xmm13
	punpckldq	%xmm14, %xmm15
	movdqa	%xmm13, %xmm0
	movdqa	%xmm15, %xmm1
	movdqa	%xmm15, %xmm2
	movdqa	%xmm13, %xmm6
	pmuludq	%xmm13, %xmm0
	pmuludq	%xmm15, %xmm1
	psrlq	$32, %xmm2
	pmuludq	%xmm2, %xmm15
	psrlq	$32, %xmm6
	paddq	%xmm15, %xmm15
	pmuludq	%xmm6, %xmm13
	psllq	$32, %xmm15
	paddq	%xmm13, %xmm13
	paddq	%xmm1, %xmm15
	psllq	$32, %xmm13
	paddq	%xmm12, %xmm15
	paddq	%xmm13, %xmm0
	paddq	%xmm15, %xmm0
	ja	.L11
	.p2align 4,,10
	.p2align 3
.L48:
	movdqa	%xmm0, %xmm3
	addq	%rbx, %r8
	psrldq	$8, %xmm3
	paddq	%xmm3, %xmm0
	movq	%xmm0, %rdx
	addq	%rdx, %rax
	cmpq	%r11, %rbx
	je	.L5
.L10:
	movl	(%rsi,%r8,4), %r14d
	subl	(%rcx,%r8,4), %r14d
	leaq	1(%r8), %rbp
	leaq	0(,%r8,4), %r13
	movslq	%r14d, %r15
	imulq	%r15, %r15
	addq	%r15, %rax
	cmpq	%rbp, %rdi
	jle	.L5
	movl	4(%rsi,%r13), %r9d
	subl	4(%rcx,%r13), %r9d
	addq	$2, %r8
	movslq	%r9d, %r10
	imulq	%r10, %r10
	addq	%r10, %rax
	cmpq	%r8, %rdi
	jle	.L5
	movl	8(%rsi,%r13), %esi
	subl	8(%rcx,%r13), %esi
	movslq	%esi, %rcx
	imulq	%rcx, %rcx
	addq	%rcx, %rax
.L5:
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
.L30:
	.cfi_restore_state
	testq	%rdi, %rdi
	jle	.L22
	leaq	-4(%rdi), %r15
	leaq	-1(%rdi), %rbx
	shrq	$2, %r15
	cmpq	$6, %rbx
	leaq	1(%r15), %r12
	leaq	0(,%r12,4), %rax
	movq	%rax, -24(%rsp)
	jbe	.L29
	leaq	0(,%r8,4), %r10
	leaq	0(,%rdx,4), %rbx
	movq	%r8, %r9
	salq	$4, %r9
	movd	(%rsi), %xmm0
	movq	%rdx, %rbp
	leaq	(%rcx,%r10), %rax
	movq	%r9, -16(%rsp)
	movq	%r10, %r14
	leaq	(%rsi,%rbx), %r9
	negq	%r14
	movd	(%rsi,%rbx,2), %xmm4
	movd	(%rax,%r10), %xmm9
	pxor	%xmm8, %xmm8
	movd	(%r9,%rbx,2), %xmm7
	salq	$4, %rbp
	movl	$1, %r13d
	movd	(%r9), %xmm5
	andl	$1, %r15d
	addq	%rbp, %r9
	movd	(%rax,%r10,2), %xmm10
	punpckldq	%xmm7, %xmm4
	movd	(%rax,%r14), %xmm11
	punpckldq	%xmm5, %xmm0
	movd	(%rax), %xmm12
	punpckldq	%xmm10, %xmm9
	movdqa	%xmm8, %xmm13
	addq	-16(%rsp), %rax
	punpckldq	%xmm12, %xmm11
	cmpq	%r12, %r13
	punpcklqdq	%xmm4, %xmm0
	leaq	(%rsi,%rbp), %r11
	punpcklqdq	%xmm9, %xmm11
	psubd	%xmm11, %xmm0
	pcmpgtd	%xmm0, %xmm13
	movdqa	%xmm0, %xmm14
	punpckhdq	%xmm13, %xmm0
	punpckldq	%xmm13, %xmm14
	movdqa	%xmm0, %xmm2
	movdqa	%xmm14, %xmm15
	movdqa	%xmm14, %xmm1
	movdqa	%xmm0, %xmm6
	pmuludq	%xmm0, %xmm2
	pmuludq	%xmm14, %xmm15
	psrlq	$32, %xmm1
	pmuludq	%xmm1, %xmm14
	psrlq	$32, %xmm6
	paddq	%xmm14, %xmm14
	pmuludq	%xmm6, %xmm0
	psllq	$32, %xmm14
	paddq	%xmm0, %xmm0
	paddq	%xmm14, %xmm15
	psllq	$32, %xmm0
	paddq	%xmm2, %xmm0
	paddq	%xmm15, %xmm0
	jnb	.L49
	testq	%r15, %r15
	je	.L51
	movd	(%r11,%rbx,2), %xmm3
	movdqa	%xmm8, %xmm13
	movd	(%r9,%rbx,2), %xmm4
	movl	$2, %r13d
	movd	(%r11), %xmm9
	addq	%rbp, %r11
	movd	(%r9), %xmm7
	punpckldq	%xmm4, %xmm3
	movd	(%rax,%r10), %xmm5
	addq	%rbp, %r9
	movd	(%rax,%r10,2), %xmm10
	punpckldq	%xmm7, %xmm9
	movd	(%rax,%r14), %xmm11
	movd	(%rax), %xmm12
	punpckldq	%xmm10, %xmm5
	punpcklqdq	%xmm3, %xmm9
	addq	-16(%rsp), %rax
	punpckldq	%xmm12, %xmm11
	cmpq	%r12, %r13
	punpcklqdq	%xmm5, %xmm11
	psubd	%xmm11, %xmm9
	pcmpgtd	%xmm9, %xmm13
	movdqa	%xmm9, %xmm14
	punpckhdq	%xmm13, %xmm9
	punpckldq	%xmm13, %xmm14
	movdqa	%xmm9, %xmm6
	movdqa	%xmm14, %xmm15
	movdqa	%xmm14, %xmm1
	movdqa	%xmm9, %xmm2
	pmuludq	%xmm9, %xmm6
	pmuludq	%xmm14, %xmm15
	psrlq	$32, %xmm1
	pmuludq	%xmm1, %xmm14
	psrlq	$32, %xmm2
	paddq	%xmm14, %xmm14
	pmuludq	%xmm2, %xmm9
	psllq	$32, %xmm14
	paddq	%xmm9, %xmm9
	paddq	%xmm14, %xmm15
	psllq	$32, %xmm9
	paddq	%xmm15, %xmm0
	paddq	%xmm6, %xmm9
	paddq	%xmm9, %xmm0
	jnb	.L49
.L51:
	movq	-16(%rsp), %r15
.L17:
	movd	(%r9,%rbx,2), %xmm4
	movdqa	%xmm8, %xmm13
	movd	(%r9), %xmm7
	addq	%rbp, %r9
	addq	$2, %r13
	movd	(%rax,%r10,2), %xmm10
	movd	(%rax), %xmm12
	movd	(%r11,%rbx,2), %xmm3
	movd	(%r11), %xmm9
	addq	%rbp, %r11
	movd	(%rax,%r10), %xmm5
	punpckldq	%xmm4, %xmm3
	movd	(%rax,%r14), %xmm11
	punpckldq	%xmm7, %xmm9
	punpckldq	%xmm10, %xmm5
	addq	%r15, %rax
	punpckldq	%xmm12, %xmm11
	movd	(%r9,%rbx,2), %xmm4
	punpcklqdq	%xmm3, %xmm9
	movd	(%r9), %xmm7
	movd	(%r11,%rbx,2), %xmm3
	addq	%rbp, %r9
	punpcklqdq	%xmm5, %xmm11
	movd	(%rax,%r10,2), %xmm10
	movd	(%rax,%r10), %xmm5
	punpckldq	%xmm4, %xmm3
	movd	(%rax), %xmm12
	psubd	%xmm11, %xmm9
	pcmpgtd	%xmm9, %xmm13
	movdqa	%xmm9, %xmm14
	movd	(%rax,%r14), %xmm11
	punpckldq	%xmm10, %xmm5
	punpckldq	%xmm13, %xmm14
	addq	%r15, %rax
	punpckhdq	%xmm13, %xmm9
	punpckldq	%xmm12, %xmm11
	movdqa	%xmm8, %xmm13
	movdqa	%xmm14, %xmm15
	movdqa	%xmm14, %xmm1
	movdqa	%xmm9, %xmm6
	movdqa	%xmm9, %xmm2
	pmuludq	%xmm14, %xmm15
	pmuludq	%xmm9, %xmm6
	psrlq	$32, %xmm1
	pmuludq	%xmm1, %xmm14
	psrlq	$32, %xmm2
	paddq	%xmm14, %xmm14
	pmuludq	%xmm2, %xmm9
	psllq	$32, %xmm14
	paddq	%xmm9, %xmm9
	paddq	%xmm15, %xmm14
	psllq	$32, %xmm9
	paddq	%xmm14, %xmm0
	paddq	%xmm6, %xmm9
	paddq	%xmm9, %xmm0
	movd	(%r11), %xmm9
	punpcklqdq	%xmm5, %xmm11
	addq	%rbp, %r11
	cmpq	%r12, %r13
	punpckldq	%xmm7, %xmm9
	punpcklqdq	%xmm3, %xmm9
	psubd	%xmm11, %xmm9
	pcmpgtd	%xmm9, %xmm13
	movdqa	%xmm9, %xmm14
	punpckhdq	%xmm13, %xmm9
	punpckldq	%xmm13, %xmm14
	movdqa	%xmm9, %xmm6
	movdqa	%xmm14, %xmm15
	movdqa	%xmm14, %xmm1
	movdqa	%xmm9, %xmm2
	pmuludq	%xmm9, %xmm6
	pmuludq	%xmm14, %xmm15
	psrlq	$32, %xmm1
	pmuludq	%xmm1, %xmm14
	psrlq	$32, %xmm2
	paddq	%xmm14, %xmm14
	pmuludq	%xmm2, %xmm9
	psllq	$32, %xmm14
	paddq	%xmm9, %xmm9
	paddq	%xmm15, %xmm14
	psllq	$32, %xmm9
	paddq	%xmm14, %xmm0
	paddq	%xmm6, %xmm9
	paddq	%xmm9, %xmm0
	jb	.L17
.L49:
	movq	-24(%rsp), %r11
	movdqa	%xmm0, %xmm8
	movq	%rdx, %r14
	movq	%r8, %r13
	psrldq	$8, %xmm8
	paddq	%xmm8, %xmm0
	movq	%xmm0, %rax
	imulq	%r11, %r14
	imulq	%r11, %r13
	cmpq	%r11, %rdi
	je	.L5
.L16:
	movl	(%rsi,%r14,4), %r12d
	subl	(%rcx,%r13,4), %r12d
	leaq	1(%r11), %rbx
	addq	%rdx, %r14
	addq	%r8, %r13
	movslq	%r12d, %rbp
	imulq	%rbp, %rbp
	addq	%rbp, %rax
	cmpq	%rbx, %rdi
	jle	.L5
	movl	(%rsi,%r14,4), %r10d
	subl	(%rcx,%r13,4), %r10d
	leaq	2(%r11), %r15
	addq	%rdx, %r14
	addq	%r8, %r13
	movslq	%r10d, %r9
	imulq	%r9, %r9
	addq	%r9, %rax
	cmpq	%r15, %rdi
	jle	.L5
	movl	(%rsi,%r14,4), %r12d
	subl	(%rcx,%r13,4), %r12d
	leaq	3(%r11), %rbx
	addq	%rdx, %r14
	addq	%r8, %r13
	movslq	%r12d, %rbp
	imulq	%rbp, %rbp
	addq	%rbp, %rax
	cmpq	%rbx, %rdi
	jle	.L5
	movl	(%rsi,%r14,4), %r10d
	subl	(%rcx,%r13,4), %r10d
	leaq	4(%r11), %r15
	addq	%rdx, %r14
	addq	%r8, %r13
	movslq	%r10d, %r9
	imulq	%r9, %r9
	addq	%r9, %rax
	cmpq	%r15, %rdi
	jle	.L5
	movl	(%rsi,%r14,4), %r12d
	subl	(%rcx,%r13,4), %r12d
	leaq	5(%r11), %rbx
	addq	%rdx, %r14
	addq	%r8, %r13
	movslq	%r12d, %rbp
	imulq	%rbp, %rbp
	addq	%rbp, %rax
	cmpq	%rbx, %rdi
	jle	.L5
	movl	(%rsi,%r14,4), %r10d
	subl	(%rcx,%r13,4), %r10d
	addq	$6, %r11
	addq	%r14, %rdx
	addq	%r13, %r8
	movslq	%r10d, %r9
	imulq	%r9, %r9
	addq	%r9, %rax
	cmpq	%r11, %rdi
	jle	.L5
	movl	(%rsi,%rdx,4), %edi
	subl	(%rcx,%r8,4), %edi
	popq	%rbx
	.cfi_remember_state
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	movslq	%edi, %r8
	imulq	%r8, %r8
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	addq	%r8, %rax
	ret
	.p2align 4,,10
	.p2align 3
.L22:
	.cfi_restore_state
	xorl	%eax, %eax
	jmp	.L5
	.p2align 4,,10
	.p2align 3
.L29:
	xorl	%r13d, %r13d
	xorl	%r14d, %r14d
	xorl	%r11d, %r11d
	xorl	%eax, %eax
	jmp	.L16
	.p2align 4,,10
	.p2align 3
.L56:
	testq	%rdx, %rdx
	jne	.L6
	xorl	%r8d, %r8d
	xorl	%eax, %eax
	jmp	.L7
	.p2align 4,,10
	.p2align 3
.L28:
	movl	$4, %r8d
	jmp	.L8
	.p2align 4,,10
	.p2align 3
.L25:
	movl	$1, %r8d
	jmp	.L8
	.p2align 4,,10
	.p2align 3
.L27:
	movl	$3, %r8d
	jmp	.L8
	.p2align 4,,10
	.p2align 3
.L26:
	movl	$2, %r8d
	jmp	.L8
	.cfi_endproc
.LFE0:
	.size	inrm2_k, .-inrm2_k
	.section	.text.unlikely
.LCOLDE0:
	.text
.LHOTE0:
	.ident	"GCC: (GNU) 5.5.0"
	.section	.note.GNU-stack,"",@progbits
