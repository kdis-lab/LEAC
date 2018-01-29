/*! \file dsum_kh.S	
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
	.p2align 4,,15
	.globl	ssum_kh
	.type	ssum_kh, @function
ssum_kh:
.LFB0:
	.cfi_startproc
	testq	%rsi, %rsi
	jle	.L4
	leaq	-4(,%rsi,4), %rdx
	leaq	(%rdi,%rsi,4), %rcx
	leaq	4(%rdi), %rax
	pxor	%xmm0, %xmm0
	shrq	$2, %rdx
	andl	$7, %edx
	cmpq	%rcx, %rax
	addss	(%rdi), %xmm0
	je	.L43
	testq	%rdx, %rdx
	je	.L3
	cmpq	$1, %rdx
	je	.L29
	cmpq	$2, %rdx
	je	.L30
	cmpq	$3, %rdx
	je	.L31
	cmpq	$4, %rdx
	je	.L32
	cmpq	$5, %rdx
	je	.L33
	cmpq	$6, %rdx
	je	.L34
	addss	(%rax), %xmm0
	leaq	8(%rdi), %rax
.L34:
	addss	(%rax), %xmm0
	addq	$4, %rax
.L33:
	addss	(%rax), %xmm0
	addq	$4, %rax
.L32:
	addss	(%rax), %xmm0
	addq	$4, %rax
.L31:
	addss	(%rax), %xmm0
	addq	$4, %rax
.L30:
	addss	(%rax), %xmm0
	addq	$4, %rax
.L29:
	addss	(%rax), %xmm0
	addq	$4, %rax
	cmpq	%rcx, %rax
	je	.L44
.L3:
	addss	(%rax), %xmm0
	addq	$32, %rax
	addss	-28(%rax), %xmm0
	addss	-24(%rax), %xmm0
	addss	-20(%rax), %xmm0
	addss	-16(%rax), %xmm0
	addss	-12(%rax), %xmm0
	addss	-8(%rax), %xmm0
	addss	-4(%rax), %xmm0
	cmpq	%rcx, %rax
	jne	.L3
	rep ret
	.p2align 4,,10
	.p2align 3
.L44:
	rep ret
	.p2align 4,,10
	.p2align 3
.L4:
	pxor	%xmm0, %xmm0
	ret
	.p2align 4,,10
	.p2align 3
.L43:
	rep ret
	.cfi_endproc
.LFE0:
	.size	ssum_kh, .-ssum_kh
	.section	.text.unlikely
.LCOLDE1:
	.text
.LHOTE1:
	.ident	"GCC: (Debian 4.9.2-10) 4.9.2"
	.section	.note.GNU-stack,"",@progbits
