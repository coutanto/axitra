################## PARAMETERS TO MODIFY *********************************

CC=icc		#	C Compiler
CFLAGS=-O3 	# 	and its options

#
# GFORTRAN compiler
# uncomment the following lines to compile with gcc
# fortran compiler
FC=gfortran
FFLAGS = -O3 -cpp -fopenmp
# end of GFORTRAN


#
# Intel ifort compiler (>V15)
# uncomment the following lines to compile with Intel
# compiler
# FC=ifort
# debug with openmp
#FFLAGS = -g -qopenmp -cpp -assume byterecl
# debug without openmp
#FFLAGS = -g -cpp -assume byterecl
#
# run without openmp
# FFLAGS = -O3 -cpp -assume byterecl
# run with openmp
# FFLAGS = -O3 -cpp -qopenmp -assume byterecl
# end of ifort

################## END OF MODIFICATIONS *************************
################## Nothing to be changed below #################
.SUFFIXES: .o .f90
.SUFFIXES: .o .f

.c.o:
	$(CC) $(CFLAGS) -c $<

.f90.o:	$(INC)
	$(FC) $(FFLAGS) -c $<

all:	axitra convms

OBJ1=dimension1.o dimension2.o parameter.o allocateArray.o initdata.o reflect0.o reflect1.o reflect2.o reflect3.o \
reflect4.o reflect5.o ff0ad.o dll2km.o axitra.o
OBJ2=dll2km.o dimension1.o dimension2.o parameter.o fft2cd.o fsource.o convm.o
#######################################
#			AXITRA
#######################################

axitra:	$(INC) $(OBJ1)
	$(FC) $(FFLAGS) -o $@ $(OBJ1) 

######################################
#			CONVM
######################################

convms: $(OBJ2) writesac.o
	$(FC) $(FFLAGS) -o $@ writesac.o $(OBJ2)

######################################
#			CLEAN, TAR, ...
######################################
clean: 
	- $(RM) *.o *.mod *.a axitra convm convms 2>/dev/null

tar:	
	tar cvf axitra.tar * ;gzip axitra.tar\

