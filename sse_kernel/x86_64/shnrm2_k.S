/*! \file inrm2_k.S	
 * This file is part of the LEAC.
 *
 * Implementation of the 
 *	
 *	nom2 <-- ||x-y||2, (unsigned int int <-- ||short -short||2)
 *
 *  function for int
 *
 *
 * (c)  Hermes Robles Berumen <hermes@uaz.edu.mx>
 *
 * For the full copyright and license information, please view the LICENSE
 * file that was distributed with this source code.
 */
	
	.file	"shnrm2_k.S"
	.section	.text.unlikely,"ax",@progbits
.LCOLDB0:
	.text
.LHOTB0:
	.p2align 4,,15
	.globl	shnrm2_k
	.type	shnrm2_k, @function
shnrm2_k:
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
	jne	.L33
	cmpq	$1, %r8
	jne	.L33
	testq	%rdi, %rdi
	jle	.L23
	movq	%rsi, %rdx
	andl	$15, %edx
	shrq	%rdx
	negq	%rdx
	andl	$7, %edx
	cmpq	%rdi, %rdx
	cmova	%rdi, %rdx
	cmpq	$8, %rdi
	jg	.L53
	movq	%rdi, %rdx
.L22:
	movzwl	(%rsi), %eax
	subw	(%rcx), %ax
	cwtl
	imull	%eax, %eax
	cmpq	$1, %rdx
	je	.L25
	movzwl	2(%rsi), %r8d
	subw	2(%rcx), %r8w
	movswl	%r8w, %ebx
	imull	%ebx, %ebx
	addl	%ebx, %eax
	cmpq	$2, %rdx
	je	.L26
	movzwl	4(%rsi), %ebp
	subw	4(%rcx), %bp
	movswl	%bp, %r9d
	imull	%r9d, %r9d
	addl	%r9d, %eax
	cmpq	$3, %rdx
	je	.L27
	movzwl	6(%rsi), %r10d
	subw	6(%rcx), %r10w
	movswl	%r10w, %r11d
	imull	%r11d, %r11d
	addl	%r11d, %eax
	cmpq	$4, %rdx
	je	.L28
	movzwl	8(%rsi), %r12d
	subw	8(%rcx), %r12w
	movswl	%r12w, %r13d
	imull	%r13d, %r13d
	addl	%r13d, %eax
	cmpq	$5, %rdx
	je	.L29
	movzwl	10(%rsi), %r14d
	subw	10(%rcx), %r14w
	movswl	%r14w, %r15d
	imull	%r15d, %r15d
	addl	%r15d, %eax
	cmpq	$6, %rdx
	je	.L30
	movzwl	12(%rsi), %r8d
	subw	12(%rcx), %r8w
	movswl	%r8w, %ebx
	imull	%ebx, %ebx
	addl	%ebx, %eax
	cmpq	$8, %rdx
	jne	.L31
	movzwl	14(%rsi), %ebp
	movl	$8, %r8d
	subw	14(%rcx), %bp
	movswl	%bp, %r9d
	imull	%r9d, %r9d
	addl	%r9d, %eax
.L8:
	cmpq	%rdi, %rdx
	je	.L35
.L7:
	movq	%rdi, %r11
	leaq	-1(%rdi), %r10
	subq	%rdx, %r11
	leaq	-8(%r11), %r14
	subq	%rdx, %r10
	shrq	$3, %r14
	cmpq	$6, %r10
	leaq	1(%r14), %r15
	leaq	0(,%r15,8), %rbx
	jbe	.L10
	addq	%rdx, %rdx
	pxor	%xmm3, %xmm3
	leaq	(%rsi,%rdx), %r12
	addq	%rcx, %rdx
	andl	$1, %r14d
	movdqu	(%rdx), %xmm0
	cmpq	$1, %r15
	movdqa	%xmm3, %xmm2
	movl	$1, %r13d
	movdqa	(%r12), %xmm1
	movl	$16, %ebp
	psubw	%xmm0, %xmm1
	pcmpgtw	%xmm1, %xmm2
	movdqa	%xmm1, %xmm5
	punpckhwd	%xmm2, %xmm1
	punpcklwd	%xmm2, %xmm5
	movdqa	%xmm1, %xmm9
	movdqa	%xmm5, %xmm4
	pmuludq	%xmm1, %xmm9
	psrlq	$32, %xmm1
	pshufd	$8, %xmm9, %xmm0
	pmuludq	%xmm5, %xmm4
	pmuludq	%xmm1, %xmm1
	psrlq	$32, %xmm5
	pshufd	$8, %xmm4, %xmm6
	pmuludq	%xmm5, %xmm5
	pshufd	$8, %xmm1, %xmm10
	pshufd	$8, %xmm5, %xmm7
	punpckldq	%xmm10, %xmm0
	punpckldq	%xmm7, %xmm6
	paddd	%xmm6, %xmm0
	jbe	.L46
	testq	%r14, %r14
	je	.L11
	movdqu	16(%rdx), %xmm11
	movdqa	%xmm3, %xmm13
	cmpq	$2, %r15
	movl	$2, %r13d
	movl	$32, %ebp
	movdqa	16(%r12), %xmm12
	psubw	%xmm11, %xmm12
	pcmpgtw	%xmm12, %xmm13
	movdqa	%xmm12, %xmm14
	punpckhwd	%xmm13, %xmm12
	punpcklwd	%xmm13, %xmm14
	movdqa	%xmm12, %xmm4
	movdqa	%xmm14, %xmm15
	pmuludq	%xmm12, %xmm4
	psrlq	$32, %xmm12
	pshufd	$8, %xmm4, %xmm6
	pmuludq	%xmm14, %xmm15
	pmuludq	%xmm12, %xmm12
	psrlq	$32, %xmm14
	pshufd	$8, %xmm15, %xmm1
	pmuludq	%xmm14, %xmm14
	pshufd	$8, %xmm12, %xmm7
	pshufd	$8, %xmm14, %xmm2
	punpckldq	%xmm7, %xmm6
	punpckldq	%xmm2, %xmm1
	paddd	%xmm1, %xmm0
	paddd	%xmm6, %xmm0
	jbe	.L46
	.p2align 4,,10
	.p2align 3
.L11:
	movdqu	(%rdx,%rbp), %xmm9
	movdqa	%xmm3, %xmm11
	movdqa	%xmm3, %xmm7
	addq	$2, %r13
	movdqa	(%r12,%rbp), %xmm10
	psubw	%xmm9, %xmm10
	pcmpgtw	%xmm10, %xmm11
	movdqa	%xmm10, %xmm12
	movdqa	16(%r12,%rbp), %xmm4
	punpckhwd	%xmm11, %xmm10
	punpcklwd	%xmm11, %xmm12
	movdqa	%xmm10, %xmm1
	movdqa	%xmm12, %xmm13
	pmuludq	%xmm10, %xmm1
	psrlq	$32, %xmm10
	pshufd	$8, %xmm1, %xmm2
	pmuludq	%xmm10, %xmm10
	pshufd	$8, %xmm10, %xmm5
	pmuludq	%xmm12, %xmm13
	psrlq	$32, %xmm12
	pshufd	$8, %xmm13, %xmm14
	pmuludq	%xmm12, %xmm12
	pshufd	$8, %xmm12, %xmm15
	punpckldq	%xmm5, %xmm2
	punpckldq	%xmm15, %xmm14
	movdqa	%xmm2, %xmm6
	paddd	%xmm14, %xmm0
	paddd	%xmm0, %xmm6
	movdqu	16(%rdx,%rbp), %xmm0
	addq	$32, %rbp
	cmpq	%r13, %r15
	psubw	%xmm0, %xmm4
	pcmpgtw	%xmm4, %xmm7
	movdqa	%xmm4, %xmm8
	punpckhwd	%xmm7, %xmm4
	punpcklwd	%xmm7, %xmm8
	movdqa	%xmm4, %xmm13
	movdqa	%xmm8, %xmm9
	pmuludq	%xmm4, %xmm13
	psrlq	$32, %xmm4
	pshufd	$8, %xmm13, %xmm0
	pmuludq	%xmm8, %xmm9
	pmuludq	%xmm4, %xmm4
	psrlq	$32, %xmm8
	pshufd	$8, %xmm9, %xmm10
	pmuludq	%xmm8, %xmm8
	pshufd	$8, %xmm4, %xmm14
	pshufd	$8, %xmm8, %xmm11
	punpckldq	%xmm14, %xmm0
	punpckldq	%xmm11, %xmm10
	paddd	%xmm6, %xmm10
	paddd	%xmm10, %xmm0
	ja	.L11
.L46:
	movdqa	%xmm0, %xmm3
	addq	%rbx, %r8
	psrldq	$8, %xmm3
	paddd	%xmm3, %xmm0
	movdqa	%xmm0, %xmm15
	psrldq	$4, %xmm15
	paddd	%xmm15, %xmm0
	movd	%xmm0, %edx
	addl	%edx, %eax
	cmpq	%r11, %rbx
	je	.L35
.L10:
	movzwl	(%rsi,%r8,2), %r9d
	leaq	1(%r8), %r10
	leaq	(%r8,%r8), %r14
	subw	(%rcx,%r8,2), %r9w
	movswl	%r9w, %r15d
	imull	%r15d, %r15d
	addl	%r15d, %eax
	cmpq	%r10, %rdi
	jle	.L35
	movzwl	2(%rsi,%r14), %r12d
	leaq	2(%r8), %rbp
	subw	2(%rcx,%r14), %r12w
	movswl	%r12w, %r13d
	imull	%r13d, %r13d
	addl	%r13d, %eax
	cmpq	%rbp, %rdi
	jle	.L35
	movzwl	4(%rsi,%r14), %ebx
	leaq	3(%r8), %rdx
	subw	4(%rcx,%r14), %bx
	movswl	%bx, %r11d
	imull	%r11d, %r11d
	addl	%r11d, %eax
	cmpq	%rdx, %rdi
	jle	.L35
	movzwl	6(%rsi,%r14), %r9d
	leaq	4(%r8), %r10
	subw	6(%rcx,%r14), %r9w
	movswl	%r9w, %r15d
	imull	%r15d, %r15d
	addl	%r15d, %eax
	cmpq	%r10, %rdi
	jle	.L35
	movzwl	8(%rsi,%r14), %r12d
	leaq	5(%r8), %rbp
	subw	8(%rcx,%r14), %r12w
	movswl	%r12w, %r13d
	imull	%r13d, %r13d
	addl	%r13d, %eax
	cmpq	%rbp, %rdi
	jle	.L35
	movzwl	10(%rsi,%r14), %ebx
	addq	$6, %r8
	subw	10(%rcx,%r14), %bx
	movswl	%bx, %r11d
	imull	%r11d, %r11d
	addl	%r11d, %eax
	cmpq	%r8, %rdi
	jle	.L35
	movzwl	12(%rsi,%r14), %esi
	subw	12(%rcx,%r14), %si
	popq	%rbx
	.cfi_remember_state
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	movswl	%si, %ecx
	imull	%ecx, %ecx
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	addl	%ecx, %eax
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L33:
	.cfi_restore_state
	testq	%rdi, %rdi
	jle	.L23
	leaq	-8(%rdi), %rax
	leaq	-1(%rdi), %r14
	shrq	$3, %rax
	addq	$1, %rax
	movq	%rax, -64(%rsp)
	salq	$3, %rax
	cmpq	$6, %r14
	movq	%rax, -72(%rsp)
	jbe	.L32
	movq	%rdx, %r12
	leaq	(%rdx,%rdx), %r11
	leaq	(%r8,%r8), %r10
	salq	$4, %r12
	movq	%r8, %r13
	movq	%rsi, %rbp
	movq	%r12, %r9
	salq	$4, %r13
	leaq	(%r11,%rdx,4), %rax
	movq	%r9, %r15
	leaq	(%r10,%r8,4), %r9
	movq	%r13, %r14
	pxor	%xmm5, %xmm5
	movq	%r12, -56(%rsp)
	pxor	%xmm6, %xmm6
	movq	%r13, -48(%rsp)
	addq	%rax, %rax
	subq	%r11, %r15
	addq	%r9, %r9
	subq	%r10, %r14
	leaq	(%rsi,%r11), %r13
	leaq	(%rcx,%r10), %r12
	movq	%rcx, %rbx
	movq	%rax, -40(%rsp)
	movq	%r15, -32(%rsp)
	movq	%r9, -24(%rsp)
	movq	%r14, -16(%rsp)
	movq	$0, -112(%rsp)
.L17:
	movzwl	0(%r13,%r11,2), %r9d
	movzwl	0(%rbp,%r11,2), %eax
	movq	-32(%rsp), %r14
	movzwl	0(%r13), %r15d
	movdqa	%xmm6, %xmm1
	addq	$1, -112(%rsp)
	salq	$16, %r9
	orq	%rax, %r9
	movzwl	0(%rbp,%r14), %eax
	movq	-40(%rsp), %r14
	salq	$16, %r9
	orq	%r15, %r9
	movzwl	0(%rbp), %r15d
	salq	$16, %r9
	movzwl	0(%rbp,%r14), %r14d
	salq	$16, %rax
	orq	%r15, %r9
	movzwl	(%rbx,%r10,2), %r15d
	movq	%r9, -104(%rsp)
	movzwl	(%r12,%r10,2), %r9d
	orq	%r14, %rax
	movzwl	0(%r13,%r11,4), %r14d
	salq	$16, %rax
	salq	$16, %r9
	orq	%r14, %rax
	movzwl	0(%rbp,%r11,4), %r14d
	orq	%r15, %r9
	salq	$16, %rax
	salq	$16, %r9
	movzwl	(%rbx), %r15d
	orq	%r14, %rax
	movq	-16(%rsp), %r14
	movq	%rax, -96(%rsp)
	movzwl	(%r12), %eax
	movdqa	-104(%rsp), %xmm2
	orq	%rax, %r9
	movzwl	(%rbx,%r14), %eax
	movq	-24(%rsp), %r14
	salq	$16, %r9
	orq	%r15, %r9
	movq	-56(%rsp), %r15
	movzwl	(%rbx,%r14), %r14d
	movq	%r9, -88(%rsp)
	salq	$16, %rax
	movq	-112(%rsp), %r9
	addq	%r15, %rbp
	addq	%r15, %r13
	orq	%r14, %rax
	movzwl	(%r12,%r10,4), %r14d
	salq	$16, %rax
	orq	%r14, %rax
	movzwl	(%rbx,%r10,4), %r14d
	salq	$16, %rax
	orq	%r14, %rax
	movq	%rax, -80(%rsp)
	movq	-48(%rsp), %rax
	psubw	-88(%rsp), %xmm2
	pcmpgtw	%xmm2, %xmm1
	movdqa	%xmm2, %xmm4
	punpckhwd	%xmm1, %xmm2
	addq	%rax, %rbx
	punpcklwd	%xmm1, %xmm4
	addq	%rax, %r12
	cmpq	-64(%rsp), %r9
	movdqa	%xmm2, %xmm10
	movdqa	%xmm4, %xmm7
	pmuludq	%xmm2, %xmm10
	psrlq	$32, %xmm2
	pshufd	$8, %xmm10, %xmm11
	pmuludq	%xmm4, %xmm7
	pmuludq	%xmm2, %xmm2
	psrlq	$32, %xmm4
	pshufd	$8, %xmm7, %xmm8
	pmuludq	%xmm4, %xmm4
	pshufd	$8, %xmm2, %xmm12
	pshufd	$8, %xmm4, %xmm9
	punpckldq	%xmm12, %xmm11
	punpckldq	%xmm9, %xmm8
	paddd	%xmm8, %xmm5
	paddd	%xmm11, %xmm5
	jb	.L17
	movq	-72(%rsp), %r10
	movdqa	%xmm5, %xmm13
	movq	%rdx, %r15
	movq	%r8, %r9
	psrldq	$8, %xmm13
	paddd	%xmm13, %xmm5
	movdqa	%xmm5, %xmm0
	imulq	%r10, %r15
	imulq	%r10, %r9
	cmpq	%r10, %rdi
	psrldq	$4, %xmm0
	paddd	%xmm0, %xmm5
	movd	%xmm5, %eax
	je	.L35
.L16:
	movzwl	(%rsi,%r15,2), %r11d
	leaq	1(%r10), %rbx
	addq	%rdx, %r15
	subw	(%rcx,%r9,2), %r11w
	addq	%r8, %r9
	movswl	%r11w, %ebp
	imull	%ebp, %ebp
	addl	%ebp, %eax
	cmpq	%rbx, %rdi
	jle	.L35
	movzwl	(%rsi,%r15,2), %r13d
	leaq	2(%r10), %r14
	addq	%rdx, %r15
	subw	(%rcx,%r9,2), %r13w
	addq	%r8, %r9
	movswl	%r13w, %r12d
	imull	%r12d, %r12d
	addl	%r12d, %eax
	cmpq	%r14, %rdi
	jle	.L35
	movzwl	(%rsi,%r15,2), %r11d
	leaq	3(%r10), %rbx
	addq	%rdx, %r15
	subw	(%rcx,%r9,2), %r11w
	addq	%r8, %r9
	movswl	%r11w, %ebp
	imull	%ebp, %ebp
	addl	%ebp, %eax
	cmpq	%rbx, %rdi
	jle	.L35
	movzwl	(%rsi,%r15,2), %r13d
	leaq	4(%r10), %r14
	addq	%rdx, %r15
	subw	(%rcx,%r9,2), %r13w
	addq	%r8, %r9
	movswl	%r13w, %r12d
	imull	%r12d, %r12d
	addl	%r12d, %eax
	cmpq	%r14, %rdi
	jle	.L35
	movzwl	(%rsi,%r15,2), %r11d
	leaq	5(%r10), %rbx
	addq	%rdx, %r15
	subw	(%rcx,%r9,2), %r11w
	addq	%r8, %r9
	movswl	%r11w, %ebp
	imull	%ebp, %ebp
	addl	%ebp, %eax
	cmpq	%rbx, %rdi
	jle	.L35
	movzwl	(%rsi,%r15,2), %r13d
	addq	$6, %r10
	addq	%r15, %rdx
	subw	(%rcx,%r9,2), %r13w
	addq	%r9, %r8
	movswl	%r13w, %r12d
	imull	%r12d, %r12d
	addl	%r12d, %eax
	cmpq	%r10, %rdi
	jle	.L35
	movzwl	(%rsi,%rdx,2), %edi
	subw	(%rcx,%r8,2), %di
	movswl	%di, %r8d
	imull	%r8d, %r8d
	addl	%r8d, %eax
.L35:
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
.L23:
	.cfi_restore_state
	xorl	%eax, %eax
	jmp	.L35
	.p2align 4,,10
	.p2align 3
.L32:
	xorl	%r9d, %r9d
	xorl	%r15d, %r15d
	xorl	%r10d, %r10d
	xorl	%eax, %eax
	jmp	.L16
	.p2align 4,,10
	.p2align 3
.L53:
	testq	%rdx, %rdx
	jne	.L22
	xorl	%r8d, %r8d
	xorl	%eax, %eax
	jmp	.L7
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
.L28:
	movl	$4, %r8d
	jmp	.L8
	.p2align 4,,10
	.p2align 3
.L26:
	movl	$2, %r8d
	jmp	.L8
	.p2align 4,,10
	.p2align 3
.L31:
	movl	$7, %r8d
	jmp	.L8
	.p2align 4,,10
	.p2align 3
.L29:
	movl	$5, %r8d
	jmp	.L8
	.p2align 4,,10
	.p2align 3
.L30:
	movl	$6, %r8d
	jmp	.L8
	.cfi_endproc
.LFE0:
	.size	shnrm2_k, .-shnrm2_k
	.section	.text.unlikely
.LCOLDE0:
	.text
.LHOTE0:
	.ident	"GCC: (Debian 4.9.2-10) 4.9.2"
	.section	.note.GNU-stack,"",@progbits
