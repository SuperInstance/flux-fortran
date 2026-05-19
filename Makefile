# Fortran 2008 build
FC = gfortran
FFLAGS = -Wall -Wextra -O2

SRCS = src/flux_fracture.f90 src/flux_sediment.f90
TEST = src/flux_test.f90

.PHONY: test clean

test: flux_test
	./flux_test

flux_test: $(SRCS) $(TEST)
	$(FC) $(FFLAGS) -o $@ $^

clean:
	rm -f flux_test *.mod
