#USE: 
#make -k [DEBUG=yes] 
#
#

TOPDIR = .
include ./Makefile.system

# # Generating Makefile.conf and config.h
# DUMMY := $(shell make -C $(TOPDIR) -f Makefile.getarch CC="$(CC)" FC="$(FC)" HOSTCC="$(HOSTCC)" CFLAGS=$(GETARCH_FLAGS) BINARY=$(BINARY) USE_OPENMP=$(USE_OPENMP) TARGET_CORE=$(TARGET_CORE) all)

include $(TOPDIR)/Makefile.conf

KERNELDIR = $(TOPDIR)/$(ARCH)

$(info   Kernel build)
$(info               )
$(info   OS               ... $(OSNAME))
$(info   Architecture     ... $(ARCH))
ifndef BINARY64
$(info   BINARY           ... 32bit)
else
$(info   BINARY           ... 64bit)
PFLAGS    += -I$(TOPDIR) -m128bit-long-double  -m64 -fPIC
endif
$(info   Number of cores  ... $(NUM_CORES))
$(info   Multi threaded; Max num-threads is $(NUM_THREADS))
$(info   Library name     ... $(LIBNAME))
$(info   Library so name  ... $(LIBSONAME))
$(info   C compiler       ... $(C_COMPILER)  (command line : $(CC)))
$(info   Fortran compiler ... $(F_COMPILER)  (command line : $(FC)))
$(info   CFLAGS           ... $(CFLAGS))

DNAME   =  $(LIBNAME).so
LDFLAGS =  
LDLIBS  = 

.SUFFIXES: .S
#
KERNELDIR	= $(TOPDIR)/$(ARCH)
#
SOURCE_S  = $(KERNELDIR)/uicopya_k.S \
	    $(KERNELDIR)/icopya_k.S \
	    $(KERNELDIR)/lcopya_k.S \
	    $(KERNELDIR)/scopya_k.S \
            $(KERNELDIR)/dcopya_k.S \
	    $(KERNELDIR)/iswap_kh.S \
	    $(KERNELDIR)/iriscal_kh.S \
	    $(KERNELDIR)/icopy_k.S \
	    $(KERNELDIR)/lcopy_k.S \
	    $(KERNELDIR)/idot_kh.S \
            $(KERNELDIR)/ilaxpy_kh.S \
	    $(KERNELDIR)/iiaxpy_kh.S \
	    $(KERNELDIR)/laxpy_kh.S \
	    $(KERNELDIR)/iaysxpy_kh.S \
	    $(KERNELDIR)/daysxpy_kh.S \
	    $(KERNELDIR)/itrans_kh.S \
	    $(KERNELDIR)/dtrans_kh.S \
	    $(KERNELDIR)/isum_kh.S \
            $(KERNELDIR)/ssum_kh.S \
            $(KERNELDIR)/dsum_kh.S \
	    $(KERNELDIR)/ddnrm2_k.S \
	    $(KERNELDIR)/ssnrm2_k.S \
	    $(KERNELDIR)/inrm2_k.S \
            $(KERNELDIR)/shnrm2_k.S \
            $(KERNELDIR)/saasxpa_kh.S \
	    $(KERNELDIR)/daasxpa_kh.S

OBJECTS_S = ${SOURCE_S:.S=.o}

all:  $(LIBNAME); ln -s $(LIBNAME) ./libssekernel.a

$(LIBNAME): ${OBJECTS_S}
	ar rcs $@ $^

.S.o:   
	${CC} -c ${CFLAGS}  $< -o $@

#
clean:
	rm libssekernel.a $(LIBNAME) getarch getarch_2nd Makefile.conf config.h  $(KERNELDIR)/*.o

#install:
#	cp -v  $(LIBNAME) ./libssekernel.a

