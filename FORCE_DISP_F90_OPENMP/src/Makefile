#CC compiler
CC=cc
CFLAGS=

#
# GFORTRAN compiler
#
FC=gfortran
F2PY=-lgomp --opt="-Ofast -march=native" --f90flags="-fopenmp -cpp"
FFLAGS = -Ofast -cpp -fopenmp


#
# Intel ifort compiler (>V15)
#
#FC=ifort
#F2PY=-L/opt/intel/compilers_and_libraries_2017.4.181/mac/compiler/lib -liomp5 --f90flags="-qopenmp -fpp"
# debug with openmp
#FFLAGS = -g -qopenmp -cpp -assume byterecl
# debug without openmp
#FFLAGS = -g -cpp -assume byterecl -check
#
# run without openmp
#FFLAGS = -O3 -cpp -assume byterecl
# run with openmp
#FFLAGS = -O3 -cpp -qopenmp -assume byterecl

################## Nothing to be changed below #################
.SUFFIXES: .o .f90
.SUFFIXES: .o .f

.c.o:
	$(CC) $(CFLAGS) -c $<

.f90.o:	$(INC)
	$(FC) $(FFLAGS) -c $<
	
.f.o:	$(INC)
	$(FC) $(FFLAGS) -c $<
	
clean:
	/bin/rm *.o *.mod

veryclean:
	/bin/rm *.o axitra convms *.mod

all: axitra convms

axitra :	dimension1.o dimension2.o parameter.o allocateArray.o initdata.o reflect0.o reflect1.o reflect2.o reflect3.o reflect4.o reflect5.o ff0ad.o axitra.o dll2km.o
	$(FC) $(FFLAGS) -o axitra axitra.o initdata.o reflect0.o reflect1.o reflect2.o reflect3.o reflect4.o reflect5.o ff0ad.o dimension1.o dimension2.o parameter.o allocateArray.o dll2km.o

convms :	parameter.o fsource.o convm.o fft2cd.o writesac.o  dll2km.o
	$(FC) $(FFLAGS) -o convms parameter.o convm.o fft2cd.o writesac.o fsource.o dll2km.o

######################################
#                       PYTHON
######################################
python: parameter.o fft2cd.o fsource.o dll2km.o
	f2py -c -m convfPy  -DF2PY_REPORT_ON_ARRAY_COPY=1 $(F2PY) convfPy.f90 fft2cd.o fsource.o dll2km.o

