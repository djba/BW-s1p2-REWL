FC = pgf90
FCFLAGS = -g -fast -Mcuda=cc60,ptxinfo,lineinfo

all: bw-s1p2-gpu

bw-s1p2-gpu: *.f90 *.cuf Makefile
	$(FC) $(FCFLAGS) \
	parameters.cuf \
	host-functions.f90 \
	device-functions.cuf \
	bw-s1p2-gpu.cuf \
	-o bw-s1p2-gpu

clean:
	rm -f *.o *.mod bw-s1p2-gpu
