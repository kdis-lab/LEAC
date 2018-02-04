/*! \file daasxpa_kh.S	
 * This file is part of the LEAC.
 *
 * Implementation of the 
 *	
 *	y = y + a(y-x),
 *
 *  function for double
 *
 *  Based on Optimized BLAS libraries
 *
 * (c)  Hermes Robles Berumen <hermes@uaz.edu.mx>
 *
 * For the full copyright and license information, please view the LICENSE
 * file that was distributed with this source code.
 */
	.file	"daysxpy_kh.c"
	.section	.text.unlikely,"ax",@progbits
.LCOLDB0:
	.text
.LHOTB0:
	.p2align 4,,15
	.globl	daysxpy_kh
	.type	daysxpy_kh, @function
daysxpy_kh:
.LFB0:
	.cfi_startproc
	cmpq	$1, %r8
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	pushq	%rbx
	.cfi_def_cfa_offset 24
	.cfi_offset 3, -24
	jne	.L14
	cmpq	$1, 24(%rsp)
	jne	.L14
	testq	%rdi, %rdi
	jle	.L18
	leaq	16(%r9), %rax
	leaq	16(%rcx), %rbx
	cmpq	%rax, %rcx
	setnb	%dl
	cmpq	%rbx, %r9
	setnb	%sil
	orb	%sil, %dl
	je	.L12
	cmpq	$6, %rdi
	jbe	.L12
	movq	%r9, %rdx
	salq	$60, %rdx
	shrq	$63, %rdx
	cmpq	%rdi, %rdx
	cmova	%rdi, %rdx
	xorl	%r11d, %r11d
	testq	%rdx, %rdx
	je	.L7
	movsd	(%r9), %xmm2
	movb	$1, %r11b
	movapd	%xmm2, %xmm1
	subsd	(%rcx), %xmm1
	mulsd	%xmm0, %xmm1
	addsd	%xmm2, %xmm1
	movsd	%xmm1, (%r9)
.L7:
	subq	%rdx, %rdi
	salq	$3, %rdx
	movl	$1, %r8d
	leaq	(%r9,%rdx), %rsi
	addq	%rcx, %rdx
	leaq	-2(%rdi), %r10
	movupd	(%rdx), %xmm5
	movl	$16, %eax
	movddup	%xmm0, %xmm8
	shrq	%r10
	movapd	(%rsi), %xmm3
	leaq	1(%r10), %rbx
	andl	$7, %r10d
	movapd	%xmm3, %xmm4
	cmpq	%rbx, %r8
	leaq	(%rbx,%rbx), %rbp
	subpd	%xmm5, %xmm4
	movapd	%xmm4, %xmm6
	mulpd	%xmm8, %xmm6
	addpd	%xmm3, %xmm6
	movaps	%xmm6, (%rsi)
	jnb	.L137
	testq	%r10, %r10
	je	.L8
	cmpq	$1, %r10
	je	.L106
	cmpq	$2, %r10
	je	.L107
	cmpq	$3, %r10
	je	.L108
	cmpq	$4, %r10
	je	.L109
	cmpq	$5, %r10
	je	.L110
	cmpq	$6, %r10
	je	.L111
	movapd	16(%rsi), %xmm7
	movl	$2, %r8d
	movl	$32, %eax
	movupd	16(%rdx), %xmm9
	movapd	%xmm7, %xmm10
	subpd	%xmm9, %xmm10
	movapd	%xmm10, %xmm11
	mulpd	%xmm8, %xmm11
	addpd	%xmm7, %xmm11
	movaps	%xmm11, 16(%rsi)
.L111:
	addq	$1, %r8
	movapd	(%rsi,%rax), %xmm12
	movupd	(%rdx,%rax), %xmm13
	movapd	%xmm12, %xmm14
	subpd	%xmm13, %xmm14
	movapd	%xmm14, %xmm15
	mulpd	%xmm8, %xmm15
	addpd	%xmm12, %xmm15
	movaps	%xmm15, (%rsi,%rax)
	addq	$16, %rax
.L110:
	addq	$1, %r8
	movapd	(%rsi,%rax), %xmm1
	movupd	(%rdx,%rax), %xmm2
	movapd	%xmm1, %xmm3
	subpd	%xmm2, %xmm3
	movapd	%xmm3, %xmm5
	mulpd	%xmm8, %xmm5
	addpd	%xmm1, %xmm5
	movaps	%xmm5, (%rsi,%rax)
	addq	$16, %rax
.L109:
	addq	$1, %r8
	movapd	(%rsi,%rax), %xmm4
	movupd	(%rdx,%rax), %xmm6
	movapd	%xmm4, %xmm7
	subpd	%xmm6, %xmm7
	movapd	%xmm7, %xmm9
	mulpd	%xmm8, %xmm9
	addpd	%xmm4, %xmm9
	movaps	%xmm9, (%rsi,%rax)
	addq	$16, %rax
.L108:
	addq	$1, %r8
	movapd	(%rsi,%rax), %xmm10
	movupd	(%rdx,%rax), %xmm11
	movapd	%xmm10, %xmm12
	subpd	%xmm11, %xmm12
	movapd	%xmm12, %xmm13
	mulpd	%xmm8, %xmm13
	addpd	%xmm10, %xmm13
	movaps	%xmm13, (%rsi,%rax)
	addq	$16, %rax
.L107:
	addq	$1, %r8
	movapd	(%rsi,%rax), %xmm14
	movupd	(%rdx,%rax), %xmm15
	movapd	%xmm14, %xmm1
	subpd	%xmm15, %xmm1
	movapd	%xmm1, %xmm2
	mulpd	%xmm8, %xmm2
	addpd	%xmm14, %xmm2
	movaps	%xmm2, (%rsi,%rax)
	addq	$16, %rax
.L106:
	addq	$1, %r8
	movapd	(%rsi,%rax), %xmm3
	movupd	(%rdx,%rax), %xmm5
	movapd	%xmm3, %xmm4
	subpd	%xmm5, %xmm4
	movapd	%xmm4, %xmm6
	mulpd	%xmm8, %xmm6
	addpd	%xmm3, %xmm6
	movaps	%xmm6, (%rsi,%rax)
	addq	$16, %rax
	cmpq	%rbx, %r8
	jnb	.L137
.L8:
	addq	$8, %r8
	movapd	(%rsi,%rax), %xmm7
	movupd	(%rdx,%rax), %xmm9
	movapd	%xmm7, %xmm10
	movapd	16(%rsi,%rax), %xmm12
	subpd	%xmm9, %xmm10
	movapd	%xmm12, %xmm14
	movapd	32(%rsi,%rax), %xmm1
	movapd	%xmm1, %xmm3
	movapd	48(%rsi,%rax), %xmm4
	movapd	%xmm10, %xmm11
	movapd	64(%rsi,%rax), %xmm10
	mulpd	%xmm8, %xmm11
	addpd	%xmm7, %xmm11
	movapd	%xmm4, %xmm7
	movaps	%xmm11, (%rsi,%rax)
	movupd	16(%rdx,%rax), %xmm13
	subpd	%xmm13, %xmm14
	movapd	%xmm14, %xmm15
	movapd	80(%rsi,%rax), %xmm14
	mulpd	%xmm8, %xmm15
	addpd	%xmm12, %xmm15
	movapd	%xmm10, %xmm12
	movaps	%xmm15, 16(%rsi,%rax)
	movupd	32(%rdx,%rax), %xmm2
	subpd	%xmm2, %xmm3
	movapd	%xmm3, %xmm5
	movapd	96(%rsi,%rax), %xmm3
	mulpd	%xmm8, %xmm5
	addpd	%xmm1, %xmm5
	movapd	%xmm14, %xmm1
	movaps	%xmm5, 32(%rsi,%rax)
	movupd	48(%rdx,%rax), %xmm6
	subpd	%xmm6, %xmm7
	movapd	%xmm7, %xmm9
	movapd	112(%rsi,%rax), %xmm7
	mulpd	%xmm8, %xmm9
	addpd	%xmm4, %xmm9
	movapd	%xmm3, %xmm4
	movaps	%xmm9, 48(%rsi,%rax)
	movupd	64(%rdx,%rax), %xmm11
	subpd	%xmm11, %xmm12
	movapd	%xmm12, %xmm13
	mulpd	%xmm8, %xmm13
	addpd	%xmm10, %xmm13
	movapd	%xmm7, %xmm10
	movaps	%xmm13, 64(%rsi,%rax)
	movupd	80(%rdx,%rax), %xmm15
	subpd	%xmm15, %xmm1
	movapd	%xmm1, %xmm2
	mulpd	%xmm8, %xmm2
	addpd	%xmm14, %xmm2
	movaps	%xmm2, 80(%rsi,%rax)
	movupd	96(%rdx,%rax), %xmm5
	subpd	%xmm5, %xmm4
	movapd	%xmm4, %xmm6
	mulpd	%xmm8, %xmm6
	addpd	%xmm3, %xmm6
	movaps	%xmm6, 96(%rsi,%rax)
	movupd	112(%rdx,%rax), %xmm9
	subpd	%xmm9, %xmm10
	movapd	%xmm10, %xmm11
	mulpd	%xmm8, %xmm11
	addpd	%xmm7, %xmm11
	movaps	%xmm11, 112(%rsi,%rax)
	subq	$-128, %rax
	cmpq	%rbx, %r8
	jb	.L8
	.p2align 4,,10
	.p2align 3
.L137:
	addq	%rbp, %r11
	cmpq	%rbp, %rdi
	je	.L18
	leaq	(%r9,%r11,8), %rdi
	movsd	(%rdi), %xmm8
	movapd	%xmm8, %xmm12
	subsd	(%rcx,%r11,8), %xmm12
	mulsd	%xmm12, %xmm0
	addsd	%xmm8, %xmm0
	movsd	%xmm0, (%rdi)
.L18:
	movl	$1, %eax
	popq	%rbx
	.cfi_remember_state
	.cfi_def_cfa_offset 16
	popq	%rbp
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L14:
	.cfi_restore_state
	testq	%rdi, %rdi
	jle	.L18
	movsd	(%r9), %xmm15
	movq	24(%rsp), %r10
	leaq	-1(%rdi), %rbx
	salq	$3, %r8
	movl	$1, %esi
	movapd	%xmm15, %xmm2
	andl	$7, %ebx
	salq	$3, %r10
	subsd	(%rcx), %xmm2
	addq	%r8, %rcx
	mulsd	%xmm0, %xmm2
	addsd	%xmm15, %xmm2
	movsd	%xmm2, (%r9)
	addq	%r10, %r9
	cmpq	%rdi, %rsi
	je	.L18
	testq	%rbx, %rbx
	je	.L11
	cmpq	$1, %rbx
	je	.L94
	cmpq	$2, %rbx
	je	.L95
	cmpq	$3, %rbx
	je	.L96
	cmpq	$4, %rbx
	je	.L97
	cmpq	$5, %rbx
	je	.L98
	cmpq	$6, %rbx
	je	.L99
	movsd	(%r9), %xmm3
	movl	$2, %esi
	movapd	%xmm3, %xmm1
	subsd	(%rcx), %xmm1
	addq	%r8, %rcx
	mulsd	%xmm0, %xmm1
	addsd	%xmm3, %xmm1
	movsd	%xmm1, (%r9)
	addq	%r10, %r9
.L99:
	movsd	(%r9), %xmm5
	addq	$1, %rsi
	movapd	%xmm5, %xmm4
	subsd	(%rcx), %xmm4
	addq	%r8, %rcx
	mulsd	%xmm0, %xmm4
	addsd	%xmm5, %xmm4
	movsd	%xmm4, (%r9)
	addq	%r10, %r9
.L98:
	movsd	(%r9), %xmm6
	addq	$1, %rsi
	movapd	%xmm6, %xmm7
	subsd	(%rcx), %xmm7
	addq	%r8, %rcx
	mulsd	%xmm0, %xmm7
	addsd	%xmm6, %xmm7
	movsd	%xmm7, (%r9)
	addq	%r10, %r9
.L97:
	movsd	(%r9), %xmm9
	addq	$1, %rsi
	movapd	%xmm9, %xmm10
	subsd	(%rcx), %xmm10
	addq	%r8, %rcx
	mulsd	%xmm0, %xmm10
	addsd	%xmm9, %xmm10
	movsd	%xmm10, (%r9)
	addq	%r10, %r9
.L96:
	movsd	(%r9), %xmm11
	addq	$1, %rsi
	movapd	%xmm11, %xmm8
	subsd	(%rcx), %xmm8
	addq	%r8, %rcx
	mulsd	%xmm0, %xmm8
	addsd	%xmm11, %xmm8
	movsd	%xmm8, (%r9)
	addq	%r10, %r9
.L95:
	movsd	(%r9), %xmm12
	addq	$1, %rsi
	movapd	%xmm12, %xmm13
	subsd	(%rcx), %xmm13
	addq	%r8, %rcx
	mulsd	%xmm0, %xmm13
	addsd	%xmm12, %xmm13
	movsd	%xmm13, (%r9)
	addq	%r10, %r9
.L94:
	movsd	(%r9), %xmm14
	addq	$1, %rsi
	movapd	%xmm14, %xmm15
	subsd	(%rcx), %xmm15
	addq	%r8, %rcx
	mulsd	%xmm0, %xmm15
	addsd	%xmm14, %xmm15
	movsd	%xmm15, (%r9)
	addq	%r10, %r9
	cmpq	%rdi, %rsi
	je	.L18
.L11:
	movsd	(%r9), %xmm2
	addq	$8, %rsi
	movapd	%xmm2, %xmm3
	subsd	(%rcx), %xmm3
	addq	%r8, %rcx
	mulsd	%xmm0, %xmm3
	addsd	%xmm2, %xmm3
	movsd	%xmm3, (%r9)
	addq	%r10, %r9
	movsd	(%r9), %xmm5
	movapd	%xmm5, %xmm1
	subsd	(%rcx), %xmm1
	addq	%r8, %rcx
	mulsd	%xmm0, %xmm1
	addsd	%xmm5, %xmm1
	movsd	%xmm1, (%r9)
	addq	%r10, %r9
	movsd	(%r9), %xmm4
	movapd	%xmm4, %xmm6
	subsd	(%rcx), %xmm6
	addq	%r8, %rcx
	mulsd	%xmm0, %xmm6
	addsd	%xmm4, %xmm6
	movsd	%xmm6, (%r9)
	addq	%r10, %r9
	movsd	(%r9), %xmm7
	movapd	%xmm7, %xmm9
	subsd	(%rcx), %xmm9
	addq	%r8, %rcx
	mulsd	%xmm0, %xmm9
	addsd	%xmm7, %xmm9
	movsd	%xmm9, (%r9)
	addq	%r10, %r9
	movsd	(%r9), %xmm10
	movapd	%xmm10, %xmm11
	subsd	(%rcx), %xmm11
	addq	%r8, %rcx
	mulsd	%xmm0, %xmm11
	addsd	%xmm10, %xmm11
	movsd	%xmm11, (%r9)
	addq	%r10, %r9
	movsd	(%r9), %xmm8
	movapd	%xmm8, %xmm12
	subsd	(%rcx), %xmm12
	addq	%r8, %rcx
	mulsd	%xmm0, %xmm12
	addsd	%xmm8, %xmm12
	movsd	%xmm12, (%r9)
	addq	%r10, %r9
	movsd	(%r9), %xmm13
	movapd	%xmm13, %xmm14
	subsd	(%rcx), %xmm14
	addq	%r8, %rcx
	mulsd	%xmm0, %xmm14
	addsd	%xmm13, %xmm14
	movsd	%xmm14, (%r9)
	addq	%r10, %r9
	movsd	(%r9), %xmm15
	movapd	%xmm15, %xmm2
	subsd	(%rcx), %xmm2
	addq	%r8, %rcx
	mulsd	%xmm0, %xmm2
	addsd	%xmm15, %xmm2
	movsd	%xmm2, (%r9)
	addq	%r10, %r9
	cmpq	%rdi, %rsi
	jne	.L11
	jmp	.L18
	.p2align 4,,10
	.p2align 3
.L12:
	movsd	(%r9), %xmm13
	leaq	-1(%rdi), %rdx
	movl	$1, %r10d
	movapd	%xmm13, %xmm14
	andl	$7, %edx
	cmpq	%rdi, %r10
	subsd	(%rcx), %xmm14
	mulsd	%xmm0, %xmm14
	addsd	%xmm13, %xmm14
	movsd	%xmm14, (%r9)
	je	.L18
	testq	%rdx, %rdx
	je	.L6
	cmpq	$1, %rdx
	je	.L100
	cmpq	$2, %rdx
	je	.L101
	cmpq	$3, %rdx
	je	.L102
	cmpq	$4, %rdx
	je	.L103
	cmpq	$5, %rdx
	je	.L104
	cmpq	$6, %rdx
	je	.L105
	movsd	8(%r9), %xmm15
	movl	$2, %r10d
	movapd	%xmm15, %xmm1
	subsd	8(%rcx), %xmm1
	mulsd	%xmm0, %xmm1
	addsd	%xmm15, %xmm1
	movsd	%xmm1, 8(%r9)
.L105:
	movsd	(%r9,%r10,8), %xmm2
	movapd	%xmm2, %xmm3
	subsd	(%rcx,%r10,8), %xmm3
	mulsd	%xmm0, %xmm3
	addsd	%xmm2, %xmm3
	movsd	%xmm3, (%r9,%r10,8)
	addq	$1, %r10
.L104:
	movsd	(%r9,%r10,8), %xmm5
	movapd	%xmm5, %xmm4
	subsd	(%rcx,%r10,8), %xmm4
	mulsd	%xmm0, %xmm4
	addsd	%xmm5, %xmm4
	movsd	%xmm4, (%r9,%r10,8)
	addq	$1, %r10
.L103:
	movsd	(%r9,%r10,8), %xmm6
	movapd	%xmm6, %xmm7
	subsd	(%rcx,%r10,8), %xmm7
	mulsd	%xmm0, %xmm7
	addsd	%xmm6, %xmm7
	movsd	%xmm7, (%r9,%r10,8)
	addq	$1, %r10
.L102:
	movsd	(%r9,%r10,8), %xmm9
	movapd	%xmm9, %xmm10
	subsd	(%rcx,%r10,8), %xmm10
	mulsd	%xmm0, %xmm10
	addsd	%xmm9, %xmm10
	movsd	%xmm10, (%r9,%r10,8)
	addq	$1, %r10
.L101:
	movsd	(%r9,%r10,8), %xmm11
	movapd	%xmm11, %xmm8
	subsd	(%rcx,%r10,8), %xmm8
	mulsd	%xmm0, %xmm8
	addsd	%xmm11, %xmm8
	movsd	%xmm8, (%r9,%r10,8)
	addq	$1, %r10
.L100:
	movsd	(%r9,%r10,8), %xmm12
	movapd	%xmm12, %xmm13
	subsd	(%rcx,%r10,8), %xmm13
	mulsd	%xmm0, %xmm13
	addsd	%xmm12, %xmm13
	movsd	%xmm13, (%r9,%r10,8)
	addq	$1, %r10
	cmpq	%rdi, %r10
	je	.L18
.L6:
	movsd	(%r9,%r10,8), %xmm14
	leaq	1(%r10), %rbx
	leaq	2(%r10), %rsi
	leaq	3(%r10), %r8
	leaq	4(%r10), %rax
	leaq	5(%r10), %rbp
	movapd	%xmm14, %xmm15
	leaq	6(%r10), %r11
	leaq	7(%r10), %rdx
	subsd	(%rcx,%r10,8), %xmm15
	mulsd	%xmm0, %xmm15
	addsd	%xmm14, %xmm15
	movsd	%xmm15, (%r9,%r10,8)
	addq	$8, %r10
	cmpq	%rdi, %r10
	movsd	(%r9,%rbx,8), %xmm2
	movapd	%xmm2, %xmm1
	subsd	(%rcx,%rbx,8), %xmm1
	mulsd	%xmm0, %xmm1
	addsd	%xmm2, %xmm1
	movsd	%xmm1, (%r9,%rbx,8)
	movsd	(%r9,%rsi,8), %xmm3
	movapd	%xmm3, %xmm5
	subsd	(%rcx,%rsi,8), %xmm5
	mulsd	%xmm0, %xmm5
	addsd	%xmm3, %xmm5
	movsd	%xmm5, (%r9,%rsi,8)
	movsd	(%r9,%r8,8), %xmm4
	movapd	%xmm4, %xmm6
	subsd	(%rcx,%r8,8), %xmm6
	mulsd	%xmm0, %xmm6
	addsd	%xmm4, %xmm6
	movsd	%xmm6, (%r9,%r8,8)
	movsd	(%r9,%rax,8), %xmm7
	movapd	%xmm7, %xmm9
	subsd	(%rcx,%rax,8), %xmm9
	mulsd	%xmm0, %xmm9
	addsd	%xmm7, %xmm9
	movsd	%xmm9, (%r9,%rax,8)
	movsd	(%r9,%rbp,8), %xmm10
	movapd	%xmm10, %xmm11
	subsd	(%rcx,%rbp,8), %xmm11
	mulsd	%xmm0, %xmm11
	addsd	%xmm10, %xmm11
	movsd	%xmm11, (%r9,%rbp,8)
	movsd	(%r9,%r11,8), %xmm8
	movapd	%xmm8, %xmm12
	subsd	(%rcx,%r11,8), %xmm12
	mulsd	%xmm0, %xmm12
	addsd	%xmm8, %xmm12
	movsd	%xmm12, (%r9,%r11,8)
	movsd	(%r9,%rdx,8), %xmm13
	movapd	%xmm13, %xmm14
	subsd	(%rcx,%rdx,8), %xmm14
	mulsd	%xmm0, %xmm14
	addsd	%xmm13, %xmm14
	movsd	%xmm14, (%r9,%rdx,8)
	jne	.L6
	jmp	.L18
	.cfi_endproc
.LFE0:
	.size	daysxpy_kh, .-daysxpy_kh
	.section	.text.unlikely
.LCOLDE0:
	.text
.LHOTE0:
	.ident	"GCC: (Debian 4.9.2-10) 4.9.2"
	.section	.note.GNU-stack,"",@progbits
