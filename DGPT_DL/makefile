FFLAGS= -O0 -fcheck=all -ffree-line-length-none
#FFLAGS= -O3 -ffree-line-length-none

# Debugging
#FFLAGS = -g -ffree-line-length-none -ffpe-trap='invalid','zero','overflow' -fbacktrace -Wunused

dg:
	gfortran -c $(FFLAGS) f90_kind.f90
	gfortran -c $(FFLAGS) physconst.f90
	gfortran -c $(FFLAGS) necessary_funcs.f90
	gfortran -c $(FFLAGS) shape_fun.f90
	gfortran -c $(FFLAGS) lapack_ops.f90
	gfortran -c $(FFLAGS) quadrature.f90
	gfortran -c $(FFLAGS) heterogeneous.f90
	gfortran -c $(FFLAGS) fem.f90
	gfortran -c $(FFLAGS) CN.f90
	gfortran -c $(FFLAGS) main.f90
	gfortran -o dg *.o -lblas -llapack
	rm -f *.o *.mod

