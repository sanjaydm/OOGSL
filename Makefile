DBG = #-g #Debugger Option
OPT = -O0 # -pg 
CC  = g++-9 -std=c++11 
FF  = gfortran-9
#NOWARN = 2>&1 >/dev/null | grep -v -e '^/var/folders/*' -e '^[[:space:]]*\.section' -e '^[[:space:]]*\^[[:space:]]*~*'

hiv: hivShell.cpp Vector.o Matrix.o Problem.o Quadrature.o VirusShell.o VirusShellBC.o libLBFGS MultiMin.o Shape.o
	$(CC) $(DBG) $(OPT) -c hivShell.cpp
	$(CC) $(DBG) $(OPT) hivShell.o Vector.o Shape.o Matrix.o Problem.o MultiMin.o  VirusShell.o VirusShellBC.o Quadrature.o -L./ -lLBFGS -lgfortran -lgsl -lgslcblas -o hivShell.out

indented: indentedShell.cpp Vector.o Matrix.o Problem.o Quadrature.o Indented.o IndentedBC.o libLBFGS MultiMin.o Shape.o MultiRoot.o
	$(CC) $(DBG) $(OPT) -c indentedShell.cpp
	$(CC) $(DBG) $(OPT) indentedShell.o Vector.o Shape.o Matrix.o Problem.o MultiMin.o  MultiRoot.o Indented.o IndentedBC.o  Quadrature.o -L./ -lLBFGS -lgfortran -lgsl -lgslcblas -o indentedShell.out

jamming: jamming.cpp myfdf.o libLBFGS
	$(CC) $(DBG)  $(OPT) jamming.cpp myfdf.o -L./ -lLBFGS -lgfortran -lgsl -lgslcblas -o driver.out

Wigner: Wigner.cpp Vector.o Matrix.o libLBFGS
	$(CC) $(DBG) $(OPT) -c Matrix.o Vector.o Wigner.cpp -L./ -lLBFGS -lgfortran -lgsl -lgslcblas -o Wigner.o

testgroup: testgroup.cpp Vector.o Matrix.o Problem.o Group.o
	$(CC) $(DBG)  $(OPT) testgroup.cpp Vector.o Matrix.o Problem.o Group.o -L./ -lLBFGS -lgfortran -lgsl -lgslcblas -o testgroup.out 

test: test.cpp Vector.o Matrix.o Problem.o libLBFGS MultiMin.o MultiRoot.o Group.o 
	$(CC) $(DBG)  $(OPT) test.cpp Vector.o Matrix.o Problem.o MultiMin.o MultiRoot.o Group.o -L./ -lLBFGS -lgfortran -lgsl -lgslcblas -o test.out 

MembLJ.o: MembLJ.cpp
	$(CC) $(DBG) $(OPT) -c MembLJ.cpp 

VirusShell.o: VirusShell.cpp
	$(CC) $(DBG) $(OPT) -c VirusShell.cpp 

Indented.o: Indented.cpp
	$(CC) $(DBG) $(OPT) -c Indented.cpp 

VirusShellBC.o: VirusShellBC.cpp
	$(CC) $(DBG) $(OPT) -c  VirusShellBC.cpp 

IndentedBC.o: IndentedBC.cpp
	$(CC) $(DBG) $(OPT) -c IndentedBC.cpp 

MultiMin.o: MultiMin.cpp libLBFGS
	$(CC) $(DBG) $(OPT) -L./ -lLBFGS -lgfortran -c  MultiMin.cpp 

MultiRoot.o: MultiRoot.cpp libLBFGS
	$(CC) $(DBG) $(OPT) -L./ -lLBFGS -lgfortran -c MultiRoot.cpp 

Group.o: Group.cpp libLBFGS
	$(CC) $(DBG) $(OPT) -L./ -lLBFGS -lgfortran -c Group.cpp 

Problem.o: Problem.cpp
	$(CC) $(DBG) $(OPT) -c Problem.cpp 

Matrix.o: Matrix.cpp
	$(CC) $(DBG) $(OPT) -c Matrix.cpp 

Vector.o: Vector.cpp
	$(CC) $(DBG) $(OPT) -c Vector.cpp 

Shape.o: Shape.cpp
	$(CC) $(DBG) $(OPT) -c Shape.cpp 

Quadrature.o: Quadrature.cpp
	$(CC) $(DBG) $(OPT) -c Quadrature.cpp 

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
