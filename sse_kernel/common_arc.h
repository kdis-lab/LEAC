#ifndef COMMON_H
#define COMMON_H

#include "config.h"

#ifdef ARCH_ALPHA
#include "common_alpha.h"
#endif

#ifdef ARCH_X86
#include "common_x86.h"
#endif

#ifdef ARCH_X86_64
#include "common_x86_64.h"
#endif

#ifdef ARCH_IA64
#include "common_ia64.h"
#endif

#ifdef ARCH_POWER
#include "common_power.h"
#endif

#ifdef sparc
#include "common_sparc.h"
#endif

#ifdef ARCH_MIPS64
#include "common_mips64.h"
#endif

#ifdef OS_LINUX
#include "common_linux.h"
#endif

#include "param.h"

#endif
