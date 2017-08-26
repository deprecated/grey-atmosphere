F90C=gfortran
F90FLAGS=-O3
EXECS=ftest

all: $(EXECS)

clean: 
	rm -fr $(EXECS) *.o *.mod

.SUFFIXES: .f90
.f90.o:  
	$(F90C) $(F90FLAGS) -c $<

ftest.o: fspecial.o special_functions.o

ftest: ftest.o fspecial.o special_functions.o
	$(F90C) $(F90FLAGS) -o $@ $^ 


