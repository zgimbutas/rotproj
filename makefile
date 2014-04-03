FC=gfortran
FFLAGS=-O2

all: rotviarecur3 rotproj 

rotproj:
	$(FC) $(FFLAGS) rotproj_cmpl_dr.f rotproj_cmpl.f dfft.f yrecursion.f 
	./a.out

rotviarecur3:
	$(FC) $(FFLAGS) rotviarecur3_dr.f rotviarecur3.f 
	./a.out

clean:
	rm -f a.out fort.13 *.o

