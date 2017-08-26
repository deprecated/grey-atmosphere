F90C=gfortran
F90FLAGS=-O3
EXECS=ftest fqtest

all: $(EXECS)

clean: 
	rm -fr $(EXECS) *.o *.mod

.SUFFIXES:
.SUFFIXES: .f90 .o
.f90.o:  
	$(F90C) $(F90FLAGS) -c $<

ftest: ftest.o fspecial.o special_functions.o
	$(F90C) $(F90FLAGS) -o $@ $^ 

fqtest: fqtest.o fquad.o quadpack.o
	$(F90C) $(F90FLAGS) -o $@ $^ 

ftest.o: fspecial.o special_functions.o
fspecial.o: special_functions.o
fqtest.o: fquad.o quadpack.o
fquad.o: quadpack.o

