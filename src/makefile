all: vbmode.o mapping.o main.o
	ifort -mkl -O2 -o a.out main.o mapping.o vbmode.o
	rm -rf *.o *.mod
	mv a.out ../

mapping.o: mapping.f90
	ifort -c -O2 mapping.f90

vbmode.o: vbmode.f90
	ifort -c -O2 vbmode.f90

main.o: main.f90
	ifort -c -O2 main.f90

clean:
	rm -rf *.o *.mod
