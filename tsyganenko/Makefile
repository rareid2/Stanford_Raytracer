.SUFFIXES : .f95 .mod .for .for

ifeq "$(OS)" "Windows_NT"
RM = del
else
RM = rm -f
endif

sources = TS05_aka_TS04.for geopack2008.for geopack0508_adapter.for #geopack2005.for
# Uncomment the above and comment out geopack2008.for and geopack0508_adapter.for to switch back to geopack 2005

# Use the -fstatic option with g95.  With gfortran, use -fno-automatic
ifeq "$(G95)" "g95"
FLAGS = -O3 -Wall -fstatic
endif
ifeq "$(G95)" "gfortran"
FLAGS = -O3 -Wall -fno-automatic
endif

OBJECTS = TS05_aka_TS04.o geopack2008.o geopack0508_adapter.o #geopack2005.o
# Uncomment the above and comment out geopack2008.for and geopack0508_adapter.for to switch back to geopack 2005


all: libtsy.a

clean:
	${RM} *.o
	${RM} *.mod
	${RM} *.a

libtsy.a: ${OBJECTS}
	ar rc libtsy.a ${OBJECTS}

.for.o:
	${G95} ${FLAGS} -c $<

.f95.o:
	${G95} ${FLAGS} -c $<

%.o : %.mod  
