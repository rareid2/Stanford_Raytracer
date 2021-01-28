SUBDIRS = xform xform_double iri2007 tsyganenko tricubic-for gcpm fortran
# export G95 = g95
export G95 = gfortran

ifeq "${G95}" "g95"
export FLAGS = -Wall -fstatic -ffixed-line-length-132 -ffree-line-length-huge
endif
ifeq "${G95}" "gfortran"
# export FLAGS = -Wall -fno-automatic -ffixed-line-length-132 -ffree-line-length-132 -fd-lines-as-comments -finit-local-zero
export FLAGS = -Wall -ffixed-line-length-132 -ffree-line-length-132 -fd-lines-as-comments -finit-local-zero
endif
ifeq "${G95}" "gcc"
export FLAGS = -Wall
endif


.PHONY: all
all:
	$(MAKE) -C lapack-3.2.1
	$(MAKE) -C xform
	$(MAKE) -C xform_double
	$(MAKE) -C iri2007
	$(MAKE) -C tsyganenko
	$(MAKE) -C tricubic-for
	$(MAKE) -C gcpm
	$(MAKE) -C fortran

test:
	$(MAKE) -C lapack-3.2.1
	$(MAKE) -C xform
	$(MAKE) -C xform_double
	$(MAKE) -C iri2007
	$(MAKE) -C fortran
# 	$(MAKE) -C tsyganenko
# 	$(MAKE) -C tricubic-for
# 	$(MAKE) -C gcpm

clean:
	$(MAKE) -C lapack-3.2.1 clean
	$(MAKE) -C xform clean
	$(MAKE) -C xform_double clean
	$(MAKE) -C iri2007 clean
	$(MAKE) -C tsyganenko clean
	$(MAKE) -C tricubic-for clean
	$(MAKE) -C gcpm clean
	$(MAKE) -C fortran clean

tidy:
	# $(MAKE) -C lapack-3.2.1 clean
	# $(MAKE) -C xform clean
	# $(MAKE) -C xform_double clean
	$(MAKE) -C iri2007 clean
	# $(MAKE) -C tsyganenko clean
	# $(MAKE) -C tricubic-for clean
	# $(MAKE) -C gcpm clean
	$(MAKE) -C fortran clean
