#-----------------------------------------------------------------------------
# File    : makefile
# Contents: build MEX-files
# Author  : Kristian Loewe
#
# Usage   : make
#           MEX_FLAGS='-v' make
#           make | grep -v 'Warning.*gcc version'
#           make -B | grep -v 'Warning.*gcc version'
#           DEBUG=1 make -B | grep -v 'Warning.*gcc version'
#-----------------------------------------------------------------------------
.SUFFIXES:
MAKEFLAGS   += -r

CC          ?= gcc
CFBASE       = -std=c99 -Wall -Wextra -Wno-unused-parameter -Wconversion \
               -Wshadow -pedantic
DEFS        ?=
MEX_FLAGS   ?=

DEBUG       ?= 0
ifeq ($(DEBUG), 1)
  CFBASE    += -g
  CFOPT     ?= -O0
else
  CFOPT     ?= -O2
  DEFS      += -DNDEBUG
endif
CFLAGS       = $(CFBASE) -fPIC $(DEFS)

MATLABROOT   = $(dir $(realpath $(shell which matlab)))
MEXCC        = $(realpath $(MATLABROOT))/mex -largeArrayDims $(MEX_FLAGS) \
               CFLAGS='$(CFLAGS)'
MEXLD        = $(realpath $(MATLABROOT))/mex -largeArrayDims $(MEX_FLAGS)
SUF          = mexa64

_OBJDIR      = obj/$(shell uname -m)/matlab
OBJDIR       = ../$(_OBJDIR)
OUTDIR       = ../private
_DUMMY      := $(shell mkdir -p $(OBJDIR); mkdir -p $(OUTDIR))

#-----------------------------------------------------------------------------

eq = $(if $(or $(1),$(2)),$(and $(findstring $(1),$(2)),\
                                $(findstring $(2),$(1))),1)

#-----------------------------------------------------------------------------

CORRDIR      = ../../corr
CPUINFODIR   = ../../cpuinfo

INCS         = -I$(CORRDIR)/src -I$(CPUINFODIR)/src

MEX_FILES    = mxPccXxFlt mxTetraccXxFlt \
               mxPccXxDbl mxTetraccXxDbl

#-----------------------------------------------------------------------------
# Link MEX files
#-----------------------------------------------------------------------------
all: $(MEX_FILES)

mxPccXxFlt:                      $(OUTDIR)/mxPccXxFlt.$(SUF)
$(OUTDIR)/mxPccXxFlt.$(SUF):     $(CORRDIR)/$(_OBJDIR)/pcc_flt.o \
                                 $(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o
$(OUTDIR)/mxPccXxFlt.$(SUF):     $(OBJDIR)/mxPccXxFlt.o makefile
	$(MEXLD)                       $(CORRDIR)/$(_OBJDIR)/pcc_flt.o \
                                 $(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o \
                                 $(OBJDIR)/mxPccXxFlt.o -output $@

mxPccXxDbl:                      $(OUTDIR)/mxPccXxDbl.$(SUF)
$(OUTDIR)/mxPccXxDbl.$(SUF):     $(CORRDIR)/$(_OBJDIR)/pcc_dbl.o \
                                 $(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o
$(OUTDIR)/mxPccXxDbl.$(SUF):     $(OBJDIR)/mxPccXxDbl.o makefile
	$(MEXLD)                       $(CORRDIR)/$(_OBJDIR)/pcc_dbl.o \
                                 $(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o \
                                 $(OBJDIR)/mxPccXxDbl.o -output $@

mxTetraccXxFlt:                  $(OUTDIR)/mxTetraccXxFlt.$(SUF)
$(OUTDIR)/mxTetraccXxFlt.$(SUF): $(CORRDIR)/$(_OBJDIR)/tetracc_flt.o \
                                 $(CORRDIR)/$(_OBJDIR)/binarize_flt.o \
                                 $(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o
$(OUTDIR)/mxTetraccXxFlt.$(SUF): $(OBJDIR)/mxTetraccXxFlt.o makefile
	$(MEXLD)                       $(CORRDIR)/$(_OBJDIR)/tetracc_flt.o \
                                 $(CORRDIR)/$(_OBJDIR)/binarize_flt.o \
                                 $(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o \
                                 $(OBJDIR)/mxTetraccXxFlt.o -output $@

mxTetraccXxDbl:                  $(OUTDIR)/mxTetraccXxDbl.$(SUF)
$(OUTDIR)/mxTetraccXxDbl.$(SUF): $(CORRDIR)/$(_OBJDIR)/tetracc_dbl.o \
                                 $(CORRDIR)/$(_OBJDIR)/binarize_dbl.o \
                                 $(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o
$(OUTDIR)/mxTetraccXxDbl.$(SUF): $(OBJDIR)/mxTetraccXxDbl.o makefile
	$(MEXLD)                       $(CORRDIR)/$(_OBJDIR)/tetracc_dbl.o \
                                 $(CORRDIR)/$(_OBJDIR)/binarize_dbl.o \
                                 $(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o \
                                 $(OBJDIR)/mxTetraccXxDbl.o -output $@

#-----------------------------------------------------------------------------
# corr-m objects
#-----------------------------------------------------------------------------
mxPccXxFlt.o:                    $(OBJDIR)/mxPccXxFlt.o
$(OBJDIR)/mxPccXxFlt.o:          $(CORRDIR)/src/pcc.* \
                                 $(CPUINFODIR)/src/cpuinfo.*
$(OBJDIR)/mxPccXxFlt.o:          mxPccXx.c makefile
	$(MEXCC) COPTIMFLAGS='$(CFOPT)' -DREAL=float $(INCS) \
    -c mxPccXx.c -outdir $(OBJDIR); \
  mv $(OBJDIR)/mxPccXx.o $(OBJDIR)/mxPccXxFlt.o

mxTetraccXxFlt.o:                $(OBJDIR)/mxTetraccXxFlt.o
$(OBJDIR)/mxTetraccXxFlt.o:      $(CORRDIR)/src/tetracc.* \
                                 $(CORRDIR)/src/binarize.* \
                                 $(CPUINFODIR)/src/cpuinfo.*
$(OBJDIR)/mxTetraccXxFlt.o:      mxTetraccXx.c makefile
	$(MEXCC) COPTIMFLAGS='$(CFOPT)' -DREAL=float $(INCS) \
    -c mxTetraccXx.c -outdir $(OBJDIR); \
  mv $(OBJDIR)/mxTetraccXx.o $(OBJDIR)/mxTetraccXxFlt.o
#-----------------------------------------------------------------------------
mxPccXxDbl.o:                    $(OBJDIR)/mxPccXxDbl.o
$(OBJDIR)/mxPccXxDbl.o:          $(CORRDIR)/src/pcc.* \
                                 $(CPUINFODIR)/src/cpuinfo.*
$(OBJDIR)/mxPccXxDbl.o:          mxPccXx.c makefile
	$(MEXCC) COPTIMFLAGS='$(CFOPT)' -DREAL=double $(INCS) \
    -c mxPccXx.c -outdir $(OBJDIR); \
  mv $(OBJDIR)/mxPccXx.o $(OBJDIR)/mxPccXxDbl.o

mxTetraccXxDbl.o:                $(OBJDIR)/mxTetraccXxDbl.o
$(OBJDIR)/mxTetraccXxDbl.o:      $(CORRDIR)/src/tetracc.* \
                                 $(CORRDIR)/src/binarize.* \
                                 $(CPUINFODIR)/src/cpuinfo.*
$(OBJDIR)/mxTetraccXxDbl.o:      mxTetraccXx.c makefile
	$(MEXCC) COPTIMFLAGS='$(CFOPT)' -DREAL=float $(INCS) \
    -c mxTetraccXx.c -outdir $(OBJDIR); \
  mv $(OBJDIR)/mxTetraccXx.o $(OBJDIR)/mxTetraccXxDbl.o

#-----------------------------------------------------------------------------
# corr objects
#-----------------------------------------------------------------------------
MKARGS_CORR1F = -C $(CORRDIR)/src -f makefile-mex pcc_flt.o
REMAKE_CORR1F = $(shell $(MAKE) -q $(MKARGS_CORR1F) > /dev/null 2>&1; \
               echo $$?\>0 | bc)
ifneq ($(or $(call eq,$(REMAKE_CORR1F),1),$(findstring B,$(MAKEFLAGS))),)
$(CORRDIR)/$(_OBJDIR)/pcc_flt.o: FORCE
	$(MAKE) $(MKARGS_CORR1F)
FORCE:
endif

MKARGS_CORR2F = -C $(CORRDIR)/src -f makefile-mex tetracc_flt.o
REMAKE_CORR2F = $(shell $(MAKE) -q $(MKARGS_CORR2F) > /dev/null 2>&1; \
               echo $$?\>0 | bc)
ifneq ($(or $(call eq,$(REMAKE_CORR2F),1),$(findstring B,$(MAKEFLAGS))),)
$(CORRDIR)/$(_OBJDIR)/tetracc_flt.o: FORCE
	$(MAKE) $(MKARGS_CORR2F)
FORCE:
endif

MKARGS_CORR3F = -C $(CORRDIR)/src -f makefile-mex binarize_flt.o
REMAKE_CORR3F = $(shell $(MAKE) -q $(MKARGS_CORR3F) > /dev/null 2>&1; \
               echo $$?\>0 | bc)
ifneq ($(or $(call eq,$(REMAKE_CORR3F),1),$(findstring B,$(MAKEFLAGS))),)
$(CORRDIR)/$(_OBJDIR)/binarize_flt.o: FORCE
	$(MAKE) $(MKARGS_CORR3F)
FORCE:
endif
#-----------------------------------------------------------------------------
MKARGS_CORR1D = -C $(CORRDIR)/src -f makefile-mex pcc_dbl.o
REMAKE_CORR1D = $(shell $(MAKE) -q $(MKARGS_CORR1D) > /dev/null 2>&1; \
               echo $$?\>0 | bc)
ifneq ($(or $(call eq,$(REMAKE_CORR1D),1),$(findstring B,$(MAKEFLAGS))),)
$(CORRDIR)/$(_OBJDIR)/pcc_dbl.o: FORCE
	$(MAKE) $(MKARGS_CORR1D)
FORCE:
endif

MKARGS_CORR2D = -C $(CORRDIR)/src -f makefile-mex tetracc_dbl.o
REMAKE_CORR2D = $(shell $(MAKE) -q $(MKARGS_CORR2D) > /dev/null 2>&1; \
               echo $$?\>0 | bc)
ifneq ($(or $(call eq,$(REMAKE_CORR2D),1),$(findstring B,$(MAKEFLAGS))),)
$(CORRDIR)/$(_OBJDIR)/tetracc_dbl.o: FORCE
	$(MAKE) $(MKARGS_CORR2D)
FORCE:
endif

MKARGS_CORR3D = -C $(CORRDIR)/src -f makefile-mex binarize_dbl.o
REMAKE_CORR3D = $(shell $(MAKE) -q $(MKARGS_CORR3D) > /dev/null 2>&1; \
               echo $$?\>0 | bc)
ifneq ($(or $(call eq,$(REMAKE_CORR3D),1),$(findstring B,$(MAKEFLAGS))),)
$(CORRDIR)/$(_OBJDIR)/binarize_dbl.o: FORCE
	$(MAKE) $(MKARGS_CORR3D)
FORCE:
endif

#-----------------------------------------------------------------------------
# cpuinfo objects
#-----------------------------------------------------------------------------
MKARGS_CPUINFO = -C $(CPUINFODIR)/src -f makefile-mex cpuinfo.o
REMAKE_CPUINFO = $(shell $(MAKE) -q $(MKARGS_CPUINFO) > /dev/null 2>&1; \
                 echo $$?\>0 | bc)
ifneq ($(or $(call eq,$(REMAKE_CPUINFO),1),$(findstring B,$(MAKEFLAGS))),)
$(CPUINFODIR)/$(_OBJDIR)/cpuinfo.o: FORCE
	$(MAKE) $(MKARGS_CPUINFO)
FORCE:
endif
