.SUFFIXES : .f95 .mod

# FLAGS += -static-libgfortran

ifeq "$(OS)" "Windows_NT"
RM = del
EXT= .exe
# necessary to keep it from linking to a silly dll
FLAGS += -static-libgcc 
# ifeq "$(G95)" "gfortran"
# FLAGS += -static-libgfortran
# endif
else
RM = rm -f
EXT = 
endif

sources = \
	bmodel_dipole.f95 \
	pp_profile_d.f95 \
	switch_d.f95 \
	constants.f95 \
	types.f95 \
	gcpm_dens_model_adapter.f95 \
	gcpm_dens_model_buildgrid_random_helpermod.f95 \
	simple_3d_model_adapter.f95 \
	interp_dens_model_adapter.f95 \
	scattered_interp_dens_model_adapter.f95 \
	ngo_dens_model.f95 \
	ngo_3d_dens_model.f95 \
	randomsampling_mod.f95 \
	ngo_dens_model_adapter.f95 \
	ngo_3d_dens_model_adapter.f95 \
	kdtree_mod.f95 \
	util.f95 \
	raytracer.f95 \
	lsinterp_mod.f95 \
	blas.f95


FLAGS += -O3 -Wall 
#FLAGS += -fdefault-real-8 
#FLAGS += -g -fbounds-check

INCLUDES = -I../tricubic-for

LIBS = -L../xform_double -lxformd -L../xform -lxform -L../gcpm -lgcpm -L../iri2007 -liri -L../xform -lxform  -L../xform_double -lxformd -L../tricubic-for -ltricubic -L../tsyganenko -ltsy -L../lapack-3.2.1 -lblas -llapack -lblas

OBJECTS = ${sources:.f95=.o}

all: ../bin/raytracer${EXT} ../bin/gcpm_dens_model_buildgrid${EXT} ../bin/dumpmodel${EXT} ../bin/gcpm_dens_model_buildgrid_random${EXT}

clean:
	${RM} *.o
	${RM} *.mod
	${RM} *.a

../bin/dumpmodel${EXT}: dumpmodel.f95 ${OBJECTS}
	${G95} ${FLAGS} ${INCLUDES} -o ../bin/dumpmodel${EXT} dumpmodel.f95 ${OBJECTS} ${LIBS} 

../bin/gcpm_dens_model_buildgrid${EXT}: gcpm_dens_model_buildgrid.f95 ${OBJECTS}
	${G95} ${FLAGS} ${INCLUDES} -o ../bin/gcpm_dens_model_buildgrid${EXT} gcpm_dens_model_buildgrid.f95 ${OBJECTS} ${LIBS} 

../bin/gcpm_dens_model_buildgrid_random${EXT}: gcpm_dens_model_buildgrid_random.f95 ${OBJECTS}
	${G95} ${FLAGS} ${INCLUDES} -o ../bin/gcpm_dens_model_buildgrid_random${EXT} gcpm_dens_model_buildgrid_random.f95 ${OBJECTS} ${LIBS} 

../bin/raytracer${EXT}: raytracer_driver.f95 ${OBJECTS}
	${G95} ${FLAGS} ${INCLUDES} -o ../bin/raytracer${EXT} raytracer_driver.f95 ${OBJECTS} ${LIBS} 

bmodel_dipole.f95 : util.o constants.o types.o

gcpm_dens_model_adapter.f95 : util.o constants.o bmodel_dipole.o types.o

gcpm_dens_model_buildgrid_random_helpermod.f95 : util.o constants.o bmodel_dipole.o types.o

scattered_interp_dens_model_adapter.f95 : types.o kdtree_mod.o util.o constants.o bmodel_dipole.o lsinterp_mod.o

interp_dens_model_adapter.f95 : util.o constants.o bmodel_dipole.o types.o

ngo_dens_model.f95 : util.o constants.o types.o

ngo_dens_model_adapter.f95 : util.o constants.o ngo_dens_model.o bmodel_dipole.o types.o

raytracer.f95 : util.o constants.o types.o

raytracer_driver.f95 : types.o util.o constants.o ngo_dens_model_adapter.o gcpm_dens_model_adapter.o interp_dens_model_adapter.o scattered_interp_dens_model_adapter.o raytracer.o

dumpmodel.f95 : types.o util.o constants.o ngo_dens_model_adapter.o gcpm_dens_model_adapter.o interp_dens_model_adapter.o scattered_interp_dens_model_adapter.o

constants.f95 : types.o

util.f95 : types.o

randomsampling_mod.f95 : types.o kdtree_mod.o constants.o 

gcpm_dens_model_buildgrid_random.f95 : constants.o util.o kdtree_mod.o gcpm_dens_model_adapter.o randomsampling_mod.o gcpm_dens_model_buildgrid_random_helpermod.o

kdtree_mod.f95 : types.o

lsinterp_mod.f95 : kdtree_mod.o types.o util.o constants.o bmodel_dipole.o kdtree_mod.o blas.o

%.o : %.mod  

.f95.o:
	${G95} ${FLAGS} ${INCLUDES} -c $<
