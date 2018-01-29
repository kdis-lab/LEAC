	.file	"ilaxpy_kh.S"
	.section	.text.unlikely,"ax",@progbits
.LCOLDB0:
	.text
.LHOTB0:
	.p2align 4,,15
	.globl	ilaxpy_kh
	.type	ilaxpy_kh, @function
ilaxpy_kh:
.LFB0:
	.cfi_startproc
	cmpq	$1, %r8
	jne	.L8
	cmpq	$1, 8(%rsp)
	jne	.L8
	testq	%rdi, %rdi
	jle	.L10
	pxor	%xmm1, %xmm1
	leaq	-1(%rdi), %rdx
	movl	$1, %r8d
	andl	$7, %edx
	cvtsi2sdq	(%rcx), %xmm1
	mulsd	%xmm0, %xmm1
	cvttsd2si	%xmm1, %eax
	addl	%eax, (%r9)
	cmpq	%rdi, %r8
	je	.L10
	testq	%rdx, %rdx
	je	.L6
	cmpq	$1, %rdx
	je	.L66
	cmpq	$2, %rdx
	je	.L67
	cmpq	$3, %rdx
	je	.L68
	cmpq	$4, %rdx
	je	.L69
	cmpq	$5, %rdx
	je	.L70
	cmpq	$6, %rdx
	je	.L71
	pxor	%xmm2, %xmm2
	movl	$2, %r8d
	cvtsi2sdq	8(%rcx), %xmm2
	mulsd	%xmm0, %xmm2
	cvttsd2si	%xmm2, %esi
	addl	%esi, 4(%r9)
.L71:
	pxor	%xmm3, %xmm3
	cvtsi2sdq	(%rcx,%r8,8), %xmm3
	mulsd	%xmm0, %xmm3
	cvttsd2si	%xmm3, %r10d
	addl	%r10d, (%r9,%r8,4)
	addq	$1, %r8
.L70:
	pxor	%xmm4, %xmm4
	cvtsi2sdq	(%rcx,%r8,8), %xmm4
	mulsd	%xmm0, %xmm4
	cvttsd2si	%xmm4, %r11d
	addl	%r11d, (%r9,%r8,4)
	addq	$1, %r8
.L69:
	pxor	%xmm5, %xmm5
	cvtsi2sdq	(%rcx,%r8,8), %xmm5
	mulsd	%xmm0, %xmm5
	cvttsd2si	%xmm5, %edx
	addl	%edx, (%r9,%r8,4)
	addq	$1, %r8
.L68:
	pxor	%xmm6, %xmm6
	cvtsi2sdq	(%rcx,%r8,8), %xmm6
	mulsd	%xmm0, %xmm6
	cvttsd2si	%xmm6, %eax
	addl	%eax, (%r9,%r8,4)
	addq	$1, %r8
.L67:
	pxor	%xmm7, %xmm7
	cvtsi2sdq	(%rcx,%r8,8), %xmm7
	mulsd	%xmm0, %xmm7
	cvttsd2si	%xmm7, %esi
	addl	%esi, (%r9,%r8,4)
	addq	$1, %r8
.L66:
	pxor	%xmm8, %xmm8
	cvtsi2sdq	(%rcx,%r8,8), %xmm8
	mulsd	%xmm0, %xmm8
	cvttsd2si	%xmm8, %r10d
	addl	%r10d, (%r9,%r8,4)
	addq	$1, %r8
	cmpq	%rdi, %r8
	je	.L10
.L6:
	pxor	%xmm9, %xmm9
	leaq	1(%r8), %rax
	leaq	2(%r8), %rsi
	pxor	%xmm10, %xmm10
	pxor	%xmm11, %xmm11
	cvtsi2sdq	(%rcx,%r8,8), %xmm9
	mulsd	%xmm0, %xmm9
	cvtsi2sdq	(%rcx,%rax,8), %xmm10
	mulsd	%xmm0, %xmm10
	cvtsi2sdq	(%rcx,%rsi,8), %xmm11
	mulsd	%xmm0, %xmm11
	pxor	%xmm12, %xmm12
	cvttsd2si	%xmm9, %r11d
	pxor	%xmm13, %xmm13
	cvttsd2si	%xmm10, %edx
	cvttsd2si	%xmm11, %r10d
	pxor	%xmm14, %xmm14
	pxor	%xmm15, %xmm15
	pxor	%xmm1, %xmm1
	addl	%r11d, (%r9,%r8,4)
	leaq	3(%r8), %r11
	addl	%edx, (%r9,%rax,4)
	cvtsi2sdq	(%rcx,%r11,8), %xmm12
	mulsd	%xmm0, %xmm12
	addl	%r10d, (%r9,%rsi,4)
	leaq	4(%r8), %rsi
	leaq	5(%r8), %r10
	cvtsi2sdq	(%rcx,%rsi,8), %xmm13
	mulsd	%xmm0, %xmm13
	cvtsi2sdq	(%rcx,%r10,8), %xmm14
	cvttsd2si	%xmm12, %eax
	mulsd	%xmm0, %xmm14
	cvttsd2si	%xmm13, %edx
	addl	%eax, (%r9,%r11,4)
	cvttsd2si	%xmm14, %r11d
	leaq	6(%r8), %rax
	addl	%edx, (%r9,%rsi,4)
	cvtsi2sdq	(%rcx,%rax,8), %xmm15
	mulsd	%xmm0, %xmm15
	cvttsd2si	%xmm15, %esi
	addl	%r11d, (%r9,%r10,4)
	leaq	7(%r8), %r10
	addq	$8, %r8
	cvtsi2sdq	(%rcx,%r10,8), %xmm1
	mulsd	%xmm0, %xmm1
	cvttsd2si	%xmm1, %edx
	addl	%esi, (%r9,%rax,4)
	addl	%edx, (%r9,%r10,4)
	cmpq	%rdi, %r8
	jne	.L6
	.p2align 4,,10
	.p2align 3
.L10:
	movl	$1, %eax
	ret
	.p2align 4,,10
	.p2align 3
.L8:
	testq	%rdi, %rdi
	jle	.L10
	pxor	%xmm2, %xmm2
	movq	8(%rsp), %r11
	leaq	-1(%rdi), %rax
	salq	$3, %r8
	movl	$1, %edx
	andl	$7, %eax
	cvtsi2sdq	(%rcx), %xmm2
	mulsd	%xmm0, %xmm2
	salq	$2, %r11
	addq	%r8, %rcx
	cvttsd2si	%xmm2, %esi
	addl	%esi, (%r9)
	addq	%r11, %r9
	cmpq	%rdi, %rdx
	je	.L10
	testq	%rax, %rax
	je	.L7
	cmpq	$1, %rax
	je	.L60
	cmpq	$2, %rax
	je	.L61
	cmpq	$3, %rax
	je	.L62
	cmpq	$4, %rax
	je	.L63
	cmpq	$5, %rax
	je	.L64
	cmpq	$6, %rax
	je	.L65
	pxor	%xmm3, %xmm3
	movl	$2, %edx
	cvtsi2sdq	(%rcx), %xmm3
	mulsd	%xmm0, %xmm3
	addq	%r8, %rcx
	cvttsd2si	%xmm3, %r10d
	addl	%r10d, (%r9)
	addq	%r11, %r9
.L65:
	pxor	%xmm4, %xmm4
	addq	$1, %rdx
	cvtsi2sdq	(%rcx), %xmm4
	mulsd	%xmm0, %xmm4
	addq	%r8, %rcx
	cvttsd2si	%xmm4, %eax
	addl	%eax, (%r9)
	addq	%r11, %r9
.L64:
	pxor	%xmm5, %xmm5
	addq	$1, %rdx
	cvtsi2sdq	(%rcx), %xmm5
	mulsd	%xmm0, %xmm5
	addq	%r8, %rcx
	cvttsd2si	%xmm5, %esi
	addl	%esi, (%r9)
	addq	%r11, %r9
.L63:
	pxor	%xmm6, %xmm6
	addq	$1, %rdx
	cvtsi2sdq	(%rcx), %xmm6
	mulsd	%xmm0, %xmm6
	addq	%r8, %rcx
	cvttsd2si	%xmm6, %r10d
	addl	%r10d, (%r9)
	addq	%r11, %r9
.L62:
	pxor	%xmm7, %xmm7
	addq	$1, %rdx
	cvtsi2sdq	(%rcx), %xmm7
	mulsd	%xmm0, %xmm7
	addq	%r8, %rcx
	cvttsd2si	%xmm7, %eax
	addl	%eax, (%r9)
	addq	%r11, %r9
.L61:
	pxor	%xmm8, %xmm8
	addq	$1, %rdx
	cvtsi2sdq	(%rcx), %xmm8
	mulsd	%xmm0, %xmm8
	addq	%r8, %rcx
	cvttsd2si	%xmm8, %esi
	addl	%esi, (%r9)
	addq	%r11, %r9
.L60:
	pxor	%xmm9, %xmm9
	addq	$1, %rdx
	cvtsi2sdq	(%rcx), %xmm9
	mulsd	%xmm0, %xmm9
	addq	%r8, %rcx
	cvttsd2si	%xmm9, %r10d
	addl	%r10d, (%r9)
	addq	%r11, %r9
	cmpq	%rdi, %rdx
	je	.L10
.L7:
	pxor	%xmm10, %xmm10
	addq	$8, %rdx
	pxor	%xmm11, %xmm11
	pxor	%xmm12, %xmm12
	cvtsi2sdq	(%rcx), %xmm10
	addq	%r8, %rcx
	mulsd	%xmm0, %xmm10
	cvtsi2sdq	(%rcx), %xmm11
	mulsd	%xmm0, %xmm11
	addq	%r8, %rcx
	pxor	%xmm13, %xmm13
	cvtsi2sdq	(%rcx), %xmm12
	addq	%r8, %rcx
	pxor	%xmm14, %xmm14
	cvttsd2si	%xmm10, %eax
	mulsd	%xmm0, %xmm12
	cvttsd2si	%xmm11, %esi
	cvtsi2sdq	(%rcx), %xmm13
	addq	%r8, %rcx
	mulsd	%xmm0, %xmm13
	cvtsi2sdq	(%rcx), %xmm14
	addq	%r8, %rcx
	mulsd	%xmm0, %xmm14
	cvttsd2si	%xmm12, %r10d
	pxor	%xmm15, %xmm15
	pxor	%xmm1, %xmm1
	pxor	%xmm2, %xmm2
	cvtsi2sdq	(%rcx), %xmm15
	addq	%r8, %rcx
	mulsd	%xmm0, %xmm15
	addl	%eax, (%r9)
	addq	%r11, %r9
	cvtsi2sdq	(%rcx), %xmm1
	cvttsd2si	%xmm13, %eax
	addl	%esi, (%r9)
	addq	%r8, %rcx
	cvttsd2si	%xmm14, %esi
	mulsd	%xmm0, %xmm1
	cvtsi2sdq	(%rcx), %xmm2
	mulsd	%xmm0, %xmm2
	addq	%r11, %r9
	addq	%r8, %rcx
	addl	%r10d, (%r9)
	addq	%r11, %r9
	cvttsd2si	%xmm15, %r10d
	addl	%eax, (%r9)
	addq	%r11, %r9
	cvttsd2si	%xmm1, %eax
	addl	%esi, (%r9)
	addq	%r11, %r9
	cvttsd2si	%xmm2, %esi
	addl	%r10d, (%r9)
	addq	%r11, %r9
	addl	%eax, (%r9)
	addq	%r11, %r9
	addl	%esi, (%r9)
	addq	%r11, %r9
	cmpq	%rdi, %rdx
	jne	.L7
	jmp	.L10
	.cfi_endproc
.LFE0:
	.size	ilaxpy_kh, .-ilaxpy_kh
	.section	.text.unlikely
.LCOLDE0:
	.text
.LHOTE0:
	.ident	"GCC: (Debian 4.9.2-10) 4.9.2"
	.section	.note.GNU-stack,"",@progbits
