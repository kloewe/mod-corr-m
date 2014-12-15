#-----------------------------------------------------------------------------
# File    : makefile
# Contents: build MEX-files
#           NB: tested with Matlab R2011b
# Author  : Kristian Loewe
#-----------------------------------------------------------------------------
.SUFFIXES:                                   # remove built-in rules
MAKEFLAGS    += -r

SHELL           = /bin/bash
CORRDIR         = ../../corr/src
CPUINFODIR      = ../../cpuinfo/src

MATLABROOT      = $(dir $(realpath $(shell which matlab)))
MEX             = $(realpath $(MATLABROOT))/mex
OPTS            = -largeArrayDims -O
# OPTS          = -v -largeArrayDims -O      # verbose mex output
# OPTS          = -g -v -largeArrayDims -O   # verbose mex output & debugging
DEFS            = -DNDEBUG
CFLAGS          = CFLAGS='-Wall -Wextra -Wno-unused-parameter -Wconversion \
                  -pedantic -fPIC -std=c99 -march=native'
COPTIMFLAGS_PCC = COPTIMFLAGS='-O3 -funroll-loops'
COPTIMFLAGS_TCC = COPTIMFLAGS='-O3'
# COPTIMFLAGS_PCC = COPTIMFLAGS='-O2'        # debugging
# COPTIMFLAGS_TCC = $(COPTIMFLAGS_PCC)       # debugging
INCS            = -I$(CORRDIR) -I$(CPUINFODIR)
SUF             = mexa64
OUTDIR          = ../private

LIBS            = -lm -lpthread

PRGS            = $(addprefix $(OUTDIR)/,\
                    mxPccXxFlt.$(SUF)     mxPccXxDbl.$(SUF) \
                    mxTetraccXxFlt.$(SUF) mxTetraccXxDbl.$(SUF) )

#-----------------------------------------------------------------------------
# Build Programs
#-----------------------------------------------------------------------------
all: $(PRGS)

$(PRGS): | $(OUTDIR)

$(OUTDIR):
	mkdir -p $(OUTDIR)

# pcc (flt)
$(OUTDIR)/mxPccXxFlt.$(SUF): mxPccXx.c \
$(CORRDIR)/pcc.* $(CPUINFODIR)/cpuinfo.* makefile
	$(MEX) $(OPTS) $(DEFS) -DREAL=float \
    $(CFLAGS) $(COPTIMFLAGS_PCC) $(INCS) $(LIBS) \
    $< $(CORRDIR)/pcc.c $(CPUINFODIR)/cpuinfo.c \
    -outdir $(OUTDIR) -output $@

# pcc (dbl)
$(OUTDIR)/mxPccXxDbl.$(SUF): mxPccXx.c \
$(CORRDIR)/pcc.* $(CPUINFODIR)/cpuinfo.* makefile
	$(MEX) $(OPTS) $(DEFS) -DREAL=double \
    $(CFLAGS) $(COPTIMFLAGS_PCC) $(INCS) $(LIBS) \
    $< $(CORRDIR)/pcc.c $(CPUINFODIR)/cpuinfo.c \
    -outdir $(OUTDIR) -output $@

# tetracc (flt)
$(OUTDIR)/mxTetraccXxFlt.$(SUF): mxTetraccXx.c \
$(CORRDIR)/tetracc.* $(CORRDIR)/binarize.* $(CPUINFODIR)/cpuinfo.* makefile
	$(MEX) $(OPTS) $(DEFS) -DREAL=float \
    $(CFLAGS) $(COPTIMFLAGS_TCC) $(INCS) $(LIBS) \
    $< $(CORRDIR)/tetracc.c $(CORRDIR)/binarize.c $(CPUINFODIR)/cpuinfo.c \
    -outdir $(OUTDIR) -output $@

# tetracc (dbl)
$(OUTDIR)/mxTetraccXxDbl.$(SUF): mxTetraccXx.c \
$(CORRDIR)/tetracc.* $(CORRDIR)/binarize.* $(CPUINFODIR)/cpuinfo.* makefile
	$(MEX) $(OPTS) $(DEFS) -DREAL=double \
    $(CFLAGS) $(COPTIMFLAGS_TCC) $(INCS) $(LIBS) \
    $< $(CORRDIR)/tetracc.c $(CORRDIR)/binarize.c $(CPUINFODIR)/cpuinfo.c \
    -outdir $(OUTDIR) -output $@
