FC = pgf90
FCFLAGS = -O3 

all: bw-s1p2-cpu

bw-s1p2-cpu: *.f90 Makefile
	$(FC) $(FCFLAGS) \
	parameters.f90 \
	bw-s1p2-functions.f90 \
	bw-s1p2-cpu.f90 \
	-o bw-s1p2-cpu
	
clean:
	rm -f *.o *.mod bw-s1p2-cpu
