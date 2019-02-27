DBG = -g #Debugger Option
OPT = #-O2 # -pg 
CC  = g++-4.9 -std=c++11 
FF  = gfortran-4.9
#NOWARN = 2>&1 >/dev/null | grep -v -e '^/var/folders/*' -e '^[[:space:]]*\.section' -e '^[[:space:]]*\^[[:space:]]*~*'

hiv: hivShell.cpp Vector.o Matrix.o Problem.o Quadrature.o VirusShell.o libLBFGS MultiMin.o
	$(CC) $(DBG) $(OPT) -c hivShell.cpp
	$(CC) $(DBG) $(OPT) hivShell.o Vector.o Matrix.o Problem.o MultiMin.o  VirusShell.o Quadrature.o -L./ -lLBFGS -lgfortran -lgsl -lgslcblas -o hivShell.out

test: test.cpp Vector.o Matrix.o Problem.o libLBFGS MultiMin.o Group.o
	$(CC) $(DBG)  $(OPT) test.cpp Vector.o Matrix.o Problem.o MultiMin.o Group.o -L./ -lLBFGS -lgfortran -lgsl -lgslcblas -o test.out 

VirusShell.o: VirusShell.h VirusShell.cpp
	$(CC) $(DBG) $(OPT) -c VirusShell.h VirusShell.cpp 

MultiMin.o: MultiMin.h MultiMin.cpp libLBFGS
	$(CC) $(DBG) $(OPT) -L./ -lLBFGS -lgfortran -c MultiMin.h MultiMin.cpp 

Group.o: Group.h Group.cpp libLBFGS
	$(CC) $(DBG) $(OPT) -L./ -lLBFGS -lgfortran -c Group.h Group.cpp 

Problem.o: Problem.h Problem.cpp
	$(CC) $(DBG) $(OPT) -c Problem.h Problem.cpp 

Matrix.o: Matrix.h Matrix.cpp
	$(CC) $(DBG) $(OPT) -c Matrix.h Matrix.cpp 

Vector.o: Vector.h Vector.cpp
	$(CC) $(DBG) $(OPT) -c Vector.h Vector.cpp 

Quadrature.o: Quadrature.h Quadrature.cpp
	$(CC) $(DBG) $(OPT) -c Quadrature.h Quadrature.cpp 

blas.o: blas.f
	$(FF) -c blas.f -o blas.o

timer.o: timer.f
	$(FF) -c timer.f -o timer.o

routines.o: routines.f
	$(FF) -c routines.f -o routines.o

linpack.o: linpack.f
	$(FF) -c linpack.f -o linpack.o

lbfgsb-routines.o: lbfgsb-routines.f
	$(FF) -c lbfgsb-routines.f -o lbfgsb-routines.o

libLBFGS: lbfgsb-routines.o blas.o linpack.o timer.o
	ar cru libLBFGS.a lbfgsb-routines.o blas.o linpack.o timer.o
	ranlib libLBFGS.a
clean:
	rm -r *.o *.gch *.vtk *.out *.a *.txt *.dSYM
