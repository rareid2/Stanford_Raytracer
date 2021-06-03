.SUFFIXES : .f95 .mod .for .FOR

ifeq "$(OS)" "Windows_NT"
RM = del
else
RM = rm -f
endif

sources = \
	bulge.for \
	gcpm_v24.for \
	iri_ps_bridge.for \
	iri_ps_eq_bridge.for \
	iri_sm.for \
	ne_inner_ps_trough.for \
	ne_iri_cap.for \
	ne_iri_ps_trough.for \
	ne_iri_ps_trough_eq.for \
	pp_profile.for \
	switchon.for



OBJECTS = ${sources:.for=.o}

all: libgcpm.a

clean:
	${RM} *.o
	${RM} *.mod
	${RM} *.a

libgcpm.a: ${OBJECTS}
	ar rc libgcpm.a ${OBJECTS}

.for.o:
	f77 -Wall -fno-automatic -ffixed-line-length-132 -ffree-line-length-132 -fd-lines-as-comments -finit-local-zero -c $<

%.o : %.mod  