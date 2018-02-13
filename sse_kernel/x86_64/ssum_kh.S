/*! \file ssum_kh.S	
 * This file is part of the LEAC.
 *
 * Implementation of the 
 *	
 *	sum x,
 *
 *  function for float
 *
 *
 * (c)  Hermes Robles Berumen <hermes@uaz.edu.mx>
 *
 * For the full copyright and license information, please view the LICENSE
 * file that was distributed with this source code.
 */	
	.file	"ssum_kh.S"
	.section	.text.unlikely,"ax",@progbits
.LCOLDB1:
	.text
.LHOTB1:
	.p2align 5,,31
	.globl	ssum_kh
	.type	ssum_kh, @function
ssum_kh:
.LFB0:
	.cfi_startproc
	cmpq	$1, %rdx
	je	.L2
	testq	%rdi, %rdi
	jle	.L7
	salq	$2, %rdx
	pxor	%xmm0, %xmm0
	movl	$1, %eax
	leaq	-1(%rdi), %rcx
	addss	(%rsi), %xmm0
	addq	%rdx, %rsi
	andl	$7, %ecx
	cmpq	$1, %rdi
	je	.L83
	testq	%rcx, %rcx
	je	.L6
	cmpq	$1, %rcx
	je	.L62
	cmpq	$2, %rcx
	je	.L63
	cmpq	$3, %rcx
	je	.L64
	cmpq	$4, %rcx
	je	.L65
	cmpq	$5, %rcx
	je	.L66
	cmpq	$6, %rcx
	je	.L67
	addss	(%rsi), %xmm0
	movl	$2, %eax
	addq	%rdx, %rsi
.L67:
	addss	(%rsi), %xmm0
	addq	$1, %rax
	addq	%rdx, %rsi
.L66:
	addss	(%rsi), %xmm0
	addq	$1, %rax
	addq	%rdx, %rsi
.L65:
	addss	(%rsi), %xmm0
	addq	$1, %rax
	addq	%rdx, %rsi
.L64:
	addss	(%rsi), %xmm0
	addq	$1, %rax
	addq	%rdx, %rsi
.L63:
	addss	(%rsi), %xmm0
	addq	$1, %rax
	addq	%rdx, %rsi
.L62:
	addq	$1, %rax
	addss	(%rsi), %xmm0
	addq	%rdx, %rsi
	cmpq	%rax, %rdi
	je	.L84
.L6:
	addss	(%rsi), %xmm0
	addq	%rdx, %rsi
	addq	$8, %rax
	addss	(%rsi), %xmm0
	addq	%rdx, %rsi
	addss	(%rsi), %xmm0
	addq	%rdx, %rsi
	addss	(%rsi), %xmm0
	addq	%rdx, %rsi
	addss	(%rsi), %xmm0
	addq	%rdx, %rsi
	addss	(%rsi), %xmm0
	addq	%rdx, %rsi
	addss	(%rsi), %xmm0
	addq	%rdx, %rsi
	addss	(%rsi), %xmm0
	addq	%rdx, %rsi
	cmpq	%rax, %rdi
	jne	.L6
	ret
	.p2align 4,,10
	.p2align 3
.L2:
	testq	%rdi, %rdi
	jle	.L7
	leaq	(%rsi,%rdi,4), %rdx
	pxor	%xmm0, %xmm0
	addss	(%rsi), %xmm0
	addq	$4, %rsi
	leaq	-4(,%rdi,4), %rdi
	shrq	$2, %rdi
	andl	$7, %edi
	cmpq	%rsi, %rdx
	je	.L85
	testq	%rdi, %rdi
	je	.L5
	cmpq	$1, %rdi
	je	.L56
	cmpq	$2, %rdi
	je	.L57
	cmpq	$3, %rdi
	je	.L58
	cmpq	$4, %rdi
	je	.L59
	cmpq	$5, %rdi
	je	.L60
	cmpq	$6, %rdi
	je	.L61
	addss	(%rsi), %xmm0
	addq	$4, %rsi
.L61:
	addss	(%rsi), %xmm0
	addq	$4, %rsi
.L60:
	addss	(%rsi), %xmm0
	addq	$4, %rsi
.L59:
	addss	(%rsi), %xmm0
	addq	$4, %rsi
.L58:
	addss	(%rsi), %xmm0
	addq	$4, %rsi
.L57:
	addss	(%rsi), %xmm0
	addq	$4, %rsi
.L56:
	addss	(%rsi), %xmm0
	addq	$4, %rsi
	cmpq	%rsi, %rdx
	je	.L86
.L5:
	addss	(%rsi), %xmm0
	addq	$32, %rsi
	addss	-28(%rsi), %xmm0
	addss	-24(%rsi), %xmm0
	addss	-20(%rsi), %xmm0
	addss	-16(%rsi), %xmm0
	addss	-12(%rsi), %xmm0
	addss	-8(%rsi), %xmm0
	addss	-4(%rsi), %xmm0
	cmpq	%rsi, %rdx
	jne	.L5
	ret
	.p2align 4,,10
	.p2align 3
.L84:
	ret
	.p2align 4,,10
	.p2align 3
.L7:
	pxor	%xmm0, %xmm0
	ret
	.p2align 4,,10
	.p2align 3
.L83:
	ret
	.p2align 4,,10
	.p2align 3
.L86:
	ret
	.p2align 4,,10
	.p2align 3
.L85:
	ret
	.cfi_endproc
.LFE0:
	.size	ssum_kh, .-ssum_kh
	.section	.text.unlikely
.LCOLDE1:
	.text
.LHOTE1:
	.ident	"GCC: (GNU) 5.5.0"
	.section	.note.GNU-stack,"",@progbits
