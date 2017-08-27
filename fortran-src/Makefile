F90C=gfortran
M=mod
F90FLAGS=-O3
EXECS=ftest fqtest fgtest
LIBRARIES=libspecfunc.so libquad.so libfgrey.so
TESTS=test-ftest test-fqtest test-fgtest

.PHONY: test $(TESTS) all clean

all: $(LIBRARIES) $(EXECS) fg4py

clean: 
	rm -fr $(EXECS) $(LIBRARIES) *.o *.$(M)

.SUFFIXES:
.SUFFIXES: .f90 .o .$(M) .so
.f90.o:  
	$(F90C) $(F90FLAGS) -c $<
.f90.$(M):  
	$(F90C) $(F90FLAGS) -c $<

# Rule to make a python extension library with f2py
fg4py: libfgrey.so
	f2py -L. -lfgrey -m fg4py -c fg4py.f90 

# Rules for building shared libraries from the third-party code
libspecfunc.so: special_functions.o
	$(F90C) $(F90FLAGS) -shared -o $@ $^ -fPIC

libquad.so: quadpack.o
	$(F90C) $(F90FLAGS) -shared -o $@ $^ -fPIC

# And make a single shared library for all the Gray Atmosphere
# functionality.  It seems that we only need to link to this, and not
# to libquad or libspecfunc
libfgrey.so: fgrey.o fquad.o qp_interface.o fspecial.o
	$(F90C) $(F90FLAGS) -shared -o $@ $^ -fPIC -L. -lquad -lspecfunc

# Rules for linking executables
ftest: ftest.o fspecial.o
	$(F90C) $(F90FLAGS) -o $@ $^ -L. -lspecfunc

fqtest: fqtest.o fquad.o qp_interface.o
	$(F90C) $(F90FLAGS) -o $@ $^ -L. -lquad

fgtest: fgtest.o		# Everything now stuffed in shared library libfgrey
	$(F90C) $(F90FLAGS) -o $@ $^ -L. -lfgrey 

# Module dependencies
fquad.o fquad.$(M): qp_interface.$(M)
fgtest.o fgtest.$(M): fgrey.$(M) fquad.$(M) qp_interface.$(M) fspecial.$(M)
fgrey.o fgrey.$(M): fquad.$(M) qp_interface.$(M) fspecial.$(M)
ftest.o ftest.$(M): fspecial.$(M)
fqtest.o fqtest.$(M): fquad.$(M) qp_interface.$(M)

# Rules for running tests
test: $(LIBRARIES) $(EXECS) $(TESTS)

test-ftest: ftest
	@echo
	@echo "------------------------------------------------"
	@echo "               Running ./ftest"
	@echo "------------------------------------------------"
	@./ftest
	@echo

test-fqtest: fqtest
	@echo
	@echo "------------------------------------------------"
	@echo "               Running ./fqtest"
	@echo "------------------------------------------------"
	@./fqtest
	@echo

test-fgtest: fgtest
	@echo
	@echo "------------------------------------------------"
	@echo "               Running ./fgtest"
	@echo "------------------------------------------------"
	@./fgtest
	@echo
