/*! \file iiaxpy_kh.S	
 * This file is part of the LEAC.
 *
 * Implementation of the 
 *	
 *	y <-- ax + y, (int <-- int * int + int)
 *
 *  function for int
 *
 *
 * (c)  Hermes Robles Berumen <hermes@uaz.edu.mx>
 *
 * For the full copyright and license information, please view the LICENSE
 * file that was distributed with this source code.
 */
	
	
	.file	"iiaxpy_kh.S"
	.section	.text.unlikely,"ax",@progbits
.LCOLDB0:
	.text
.LHOTB0:
	.p2align 4,,15
	.globl	iiaxpy_kh
	.type	iiaxpy_kh, @function
iiaxpy_kh:
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
	jne	.L19
	cmpq	$1, 48(%rsp)
	jne	.L19
	testq	%rdi, %rdi
	jle	.L23
	leaq	16(%rax), %rdx
	leaq	16(%r8), %rbx
	cmpq	%rdx, %r8
	setnb	%sil
	cmpq	%rbx, %rax
	setnb	%bpl
	orb	%bpl, %sil
	je	.L15
	cmpq	$7, %rdi
	jbe	.L15
	movq	%rax, %rdx
	andl	$15, %edx
	shrq	$2, %rdx
	negq	%rdx
	andl	$3, %edx
	cmpq	%rdi, %rdx
	cmova	%rdi, %rdx
	xorl	%esi, %esi
	testq	%rdx, %rdx
	je	.L7
	movl	(%r8), %r9d
	movl	$1, %esi
	imull	%ecx, %r9d
	addl	%r9d, (%rax)
	cmpq	$1, %rdx
	je	.L7
	movl	4(%r8), %r10d
	movl	$2, %esi
	imull	%ecx, %r10d
	addl	%r10d, 4(%rax)
	cmpq	$3, %rdx
	jne	.L7
	movl	8(%r8), %r11d
	movl	$3, %esi
	imull	%ecx, %r11d
	addl	%r11d, 8(%rax)
.L7:
	movl	%ecx, -12(%rsp)
	movq	%rdi, %rbx
	movl	$1, %r13d
	movd	-12(%rsp), %xmm4
	subq	%rdx, %rbx
	salq	$2, %rdx
	leaq	(%rax,%rdx), %r11
	addq	%r8, %rdx
	leaq	-4(%rbx), %r9
	pshufd	$0, %xmm4, %xmm2
	movl	$16, %r10d
	movdqu	(%rdx), %xmm0
	shrq	$2, %r9
	leaq	1(%r9), %r12
	andl	$3, %r9d
	movdqa	%xmm0, %xmm1
	movdqa	%xmm2, %xmm3
	psrlq	$32, %xmm0
	movdqa	(%r11), %xmm7
	cmpq	%r12, %r13
	psrlq	$32, %xmm3
	pmuludq	%xmm2, %xmm1
	pmuludq	%xmm3, %xmm0
	pshufd	$8, %xmm1, %xmm5
	leaq	0(,%r12,4), %rbp
	pshufd	$8, %xmm0, %xmm6
	punpckldq	%xmm6, %xmm5
	paddd	%xmm5, %xmm7
	movaps	%xmm7, (%r11)
	jb	.L84
.L123:
	addq	%rbp, %rsi
	cmpq	%rbp, %rbx
	je	.L23
	movl	(%r8,%rsi,4), %r9d
	leaq	0(,%rsi,4), %rdx
	leaq	1(%rsi), %r12
	imull	%ecx, %r9d
	addl	%r9d, (%rax,%rdx)
	cmpq	%r12, %rdi
	jle	.L23
	movl	4(%r8,%rdx), %r11d
	addq	$2, %rsi
	imull	%ecx, %r11d
	addl	%r11d, 4(%rax,%rdx)
	cmpq	%rsi, %rdi
	jle	.L23
	imull	8(%r8,%rdx), %ecx
	addl	%ecx, 8(%rax,%rdx)
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
.L19:
	.cfi_restore_state
	testq	%rdi, %rdi
	jle	.L23
	movq	48(%rsp), %rbp
	movl	(%r8), %ebx
	leaq	-1(%rdi), %r10
	salq	$2, %r9
	movl	$1, %edx
	andl	$7, %r10d
	addq	%r9, %r8
	leaq	0(,%rbp,4), %r13
	imull	%ecx, %ebx
	addl	%ebx, (%rax)
	addq	%r13, %rax
	cmpq	%rdi, %rdx
	je	.L23
	testq	%r10, %r10
	je	.L14
	cmpq	$1, %r10
	je	.L87
	cmpq	$2, %r10
	je	.L88
	cmpq	$3, %r10
	je	.L89
	cmpq	$4, %r10
	je	.L90
	cmpq	$5, %r10
	je	.L91
	cmpq	$6, %r10
	je	.L92
	movl	(%r8), %esi
	movl	$2, %edx
	addq	%r9, %r8
	imull	%ecx, %esi
	addl	%esi, (%rax)
	addq	%r13, %rax
.L92:
	movl	(%r8), %r12d
	addq	$1, %rdx
	addq	%r9, %r8
	imull	%ecx, %r12d
	addl	%r12d, (%rax)
	addq	%r13, %rax
.L91:
	movl	(%r8), %r11d
	addq	$1, %rdx
	addq	%r9, %r8
	imull	%ecx, %r11d
	addl	%r11d, (%rax)
	addq	%r13, %rax
.L90:
	movl	(%r8), %ebp
	addq	$1, %rdx
	addq	%r9, %r8
	imull	%ecx, %ebp
	addl	%ebp, (%rax)
	addq	%r13, %rax
.L89:
	movl	(%r8), %r10d
	addq	$1, %rdx
	addq	%r9, %r8
	imull	%ecx, %r10d
	addl	%r10d, (%rax)
	addq	%r13, %rax
.L88:
	movl	(%r8), %ebx
	addq	$1, %rdx
	addq	%r9, %r8
	imull	%ecx, %ebx
	addl	%ebx, (%rax)
	addq	%r13, %rax
.L87:
	movl	(%r8), %esi
	addq	$1, %rdx
	addq	%r9, %r8
	imull	%ecx, %esi
	addl	%esi, (%rax)
	addq	%r13, %rax
	cmpq	%rdi, %rdx
	je	.L23
.L14:
	movl	(%r8), %r12d
	addq	%r9, %r8
	addq	$8, %rdx
	imull	%ecx, %r12d
	addl	%r12d, (%rax)
	addq	%r13, %rax
	movl	(%r8), %r11d
	addq	%r9, %r8
	imull	%ecx, %r11d
	addl	%r11d, (%rax)
	addq	%r13, %rax
	movl	(%r8), %ebp
	addq	%r9, %r8
	imull	%ecx, %ebp
	addl	%ebp, (%rax)
	addq	%r13, %rax
	movl	(%r8), %r10d
	addq	%r9, %r8
	imull	%ecx, %r10d
	addl	%r10d, (%rax)
	addq	%r13, %rax
	movl	(%r8), %ebx
	addq	%r9, %r8
	imull	%ecx, %ebx
	addl	%ebx, (%rax)
	addq	%r13, %rax
	movl	(%r8), %esi
	addq	%r9, %r8
	imull	%ecx, %esi
	addl	%esi, (%rax)
	addq	%r13, %rax
	movl	(%r8), %r12d
	addq	%r9, %r8
	imull	%ecx, %r12d
	addl	%r12d, (%rax)
	addq	%r13, %rax
	movl	(%r8), %r11d
	addq	%r9, %r8
	imull	%ecx, %r11d
	addl	%r11d, (%rax)
	addq	%r13, %rax
	cmpq	%rdi, %rdx
	jne	.L14
	jmp	.L23
	.p2align 4,,10
	.p2align 3
.L84:
	testq	%r9, %r9
	je	.L9
	cmpq	$1, %r9
	je	.L99
	cmpq	$2, %r9
	je	.L100
	movdqu	16(%rdx), %xmm8
	movl	$2, %r13d
	movl	$32, %r10d
	movdqa	%xmm8, %xmm9
	psrlq	$32, %xmm8
	pmuludq	%xmm3, %xmm8
	pshufd	$8, %xmm8, %xmm11
	movdqa	16(%r11), %xmm12
	pmuludq	%xmm2, %xmm9
	pshufd	$8, %xmm9, %xmm10
	punpckldq	%xmm11, %xmm10
	paddd	%xmm10, %xmm12
	movaps	%xmm12, 16(%r11)
.L100:
	addq	$1, %r13
	movdqu	(%rdx,%r10), %xmm13
	movdqa	%xmm13, %xmm14
	psrlq	$32, %xmm13
	pmuludq	%xmm3, %xmm13
	pshufd	$8, %xmm13, %xmm4
	movdqa	(%r11,%r10), %xmm0
	pmuludq	%xmm2, %xmm14
	pshufd	$8, %xmm14, %xmm15
	punpckldq	%xmm4, %xmm15
	paddd	%xmm15, %xmm0
	movaps	%xmm0, (%r11,%r10)
	addq	$16, %r10
.L99:
	addq	$1, %r13
	movdqu	(%rdx,%r10), %xmm5
	movdqa	%xmm5, %xmm1
	psrlq	$32, %xmm5
	pmuludq	%xmm3, %xmm5
	pshufd	$8, %xmm5, %xmm7
	movdqa	(%r11,%r10), %xmm8
	pmuludq	%xmm2, %xmm1
	pshufd	$8, %xmm1, %xmm6
	punpckldq	%xmm7, %xmm6
	paddd	%xmm6, %xmm8
	movaps	%xmm8, (%r11,%r10)
	addq	$16, %r10
	cmpq	%r12, %r13
	jnb	.L123
.L9:
	addq	$4, %r13
	movdqu	(%rdx,%r10), %xmm9
	movdqa	%xmm9, %xmm10
	psrlq	$32, %xmm9
	pmuludq	%xmm3, %xmm9
	pshufd	$8, %xmm9, %xmm12
	movdqa	(%r11,%r10), %xmm13
	pmuludq	%xmm2, %xmm10
	pshufd	$8, %xmm10, %xmm11
	movdqa	16(%r11,%r10), %xmm5
	punpckldq	%xmm12, %xmm11
	movdqa	32(%r11,%r10), %xmm9
	paddd	%xmm11, %xmm13
	movaps	%xmm13, (%r11,%r10)
	movdqu	16(%rdx,%r10), %xmm14
	movdqa	%xmm14, %xmm15
	psrlq	$32, %xmm14
	pmuludq	%xmm3, %xmm14
	pshufd	$8, %xmm14, %xmm0
	movdqa	48(%r11,%r10), %xmm14
	pmuludq	%xmm2, %xmm15
	pshufd	$8, %xmm15, %xmm4
	punpckldq	%xmm0, %xmm4
	paddd	%xmm4, %xmm5
	movaps	%xmm5, 16(%r11,%r10)
	movdqu	32(%rdx,%r10), %xmm6
	movdqa	%xmm6, %xmm1
	psrlq	$32, %xmm6
	pmuludq	%xmm3, %xmm6
	pshufd	$8, %xmm6, %xmm8
	pmuludq	%xmm2, %xmm1
	pshufd	$8, %xmm1, %xmm7
	punpckldq	%xmm8, %xmm7
	paddd	%xmm7, %xmm9
	movaps	%xmm9, 32(%r11,%r10)
	movdqu	48(%rdx,%r10), %xmm10
	movdqa	%xmm10, %xmm11
	psrlq	$32, %xmm10
	pmuludq	%xmm3, %xmm10
	pshufd	$8, %xmm10, %xmm13
	pmuludq	%xmm2, %xmm11
	pshufd	$8, %xmm11, %xmm12
	punpckldq	%xmm13, %xmm12
	paddd	%xmm12, %xmm14
	movaps	%xmm14, 48(%r11,%r10)
	addq	$64, %r10
	cmpq	%r12, %r13
	jb	.L9
	jmp	.L123
	.p2align 4,,10
	.p2align 3
.L15:
	movl	(%r8), %r10d
	leaq	-1(%rdi), %r13
	movl	$1, %ebp
	andl	$7, %r13d
	imull	%ecx, %r10d
	addl	%r10d, (%rax)
	cmpq	%rdi, %rbp
	je	.L23
	testq	%r13, %r13
	je	.L6
	cmpq	$1, %r13
	je	.L93
	cmpq	$2, %r13
	je	.L94
	cmpq	$3, %r13
	je	.L95
	cmpq	$4, %r13
	je	.L96
	cmpq	$5, %r13
	je	.L97
	cmpq	$6, %r13
	je	.L98
	movl	4(%r8), %ebx
	movl	$2, %ebp
	imull	%ecx, %ebx
	addl	%ebx, 4(%rax)
.L98:
	movl	(%r8,%rbp,4), %esi
	imull	%ecx, %esi
	addl	%esi, (%rax,%rbp,4)
	addq	$1, %rbp
.L97:
	movl	(%r8,%rbp,4), %edx
	imull	%ecx, %edx
	addl	%edx, (%rax,%rbp,4)
	addq	$1, %rbp
.L96:
	movl	(%r8,%rbp,4), %r9d
	imull	%ecx, %r9d
	addl	%r9d, (%rax,%rbp,4)
	addq	$1, %rbp
.L95:
	movl	(%r8,%rbp,4), %r12d
	imull	%ecx, %r12d
	addl	%r12d, (%rax,%rbp,4)
	addq	$1, %rbp
.L94:
	movl	(%r8,%rbp,4), %r11d
	imull	%ecx, %r11d
	addl	%r11d, (%rax,%rbp,4)
	addq	$1, %rbp
.L93:
	movl	(%r8,%rbp,4), %r13d
	imull	%ecx, %r13d
	addl	%r13d, (%rax,%rbp,4)
	addq	$1, %rbp
	cmpq	%rdi, %rbp
	je	.L23
.L6:
	movl	(%r8,%rbp,4), %r10d
	leaq	1(%rbp), %rbx
	leaq	2(%rbp), %rdx
	leaq	3(%rbp), %r12
	leaq	4(%rbp), %r13
	imull	%ecx, %r10d
	addl	%r10d, (%rax,%rbp,4)
	movl	(%r8,%rbx,4), %esi
	imull	%ecx, %esi
	addl	%esi, (%rax,%rbx,4)
	leaq	5(%rbp), %rbx
	movl	(%r8,%rdx,4), %r9d
	imull	%ecx, %r9d
	addl	%r9d, (%rax,%rdx,4)
	leaq	6(%rbp), %rdx
	movl	(%r8,%r12,4), %r11d
	imull	%ecx, %r11d
	addl	%r11d, (%rax,%r12,4)
	leaq	7(%rbp), %r12
	movl	(%r8,%r13,4), %r10d
	addq	$8, %rbp
	imull	%ecx, %r10d
	addl	%r10d, (%rax,%r13,4)
	movl	(%r8,%rbx,4), %esi
	imull	%ecx, %esi
	addl	%esi, (%rax,%rbx,4)
	movl	(%r8,%rdx,4), %r9d
	imull	%ecx, %r9d
	addl	%r9d, (%rax,%rdx,4)
	movl	(%r8,%r12,4), %r11d
	imull	%ecx, %r11d
	addl	%r11d, (%rax,%r12,4)
	cmpq	%rdi, %rbp
	jne	.L6
	jmp	.L23
	.cfi_endproc
.LFE0:
	.size	iiaxpy_kh, .-iiaxpy_kh
	.section	.text.unlikely
.LCOLDE0:
	.text
.LHOTE0:
	.ident	"GCC: (Debian 4.9.2-10) 4.9.2"
	.section	.note.GNU-stack,"",@progbits
