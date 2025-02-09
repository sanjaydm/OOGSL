DBG = -g #Debugger Option
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

compute_D_matrix: compute_D_matrix.cpp Vector.o Matrix.o libLBFGS
	$(CC) $(DBG) $(OPT) -c Matrix.o Vector.o compute_D_matrix.cpp -L./ -lLBFGS -lgfortran -lgsl -lgslcblas -o compute_D_matrix.o

testgroup: testgroup.cpp Vector.o Matrix.o Problem.o Group.o SymmReduced.o
	$(CC) $(DBG)  $(OPT) testgroup.cpp Vector.o Matrix.o Problem.o Group.o SymmReduced.o -L./ -lLBFGS -lgfortran -lgsl -lgslcblas -o testgroup.out 

test: test.cpp Vector.o Matrix.o Problem.o libLBFGS MultiMin.o MultiRoot.o Group.o SymmReduced.o 
	$(CC) $(DBG)  $(OPT) test.cpp Vector.o Matrix.o Problem.o SymmReduced.o MultiMin.o MultiRoot.o Group.o -L./ -lLBFGS -lgfortran -lgsl -lgslcblas -o test.out 

symmHelfLJ: symmHelfLJ.cpp Vector.o Matrix.o Problem.o libLBFGS MultiMin.o MultiRoot.o Group.o SymmReduced.o MembLJ.o libCont
	$(CC) $(DBG)  $(OPT) symmHelfLJ.cpp Vector.o Matrix.o Problem.o SymmReduced.o MultiMin.o MultiRoot.o Group.o MembLJ.o -L./ -lCont -lLBFGS -lgfortran -lgsl -lgslcblas -o symmHelfLJ.out 

symmHelfLJSJ: symmHelfLJSJ.cpp Vector.o Matrix.o Problem.o libLBFGS MultiMin.o MultiRoot.o Group.o SymmReduced.o MembLJ.o libCont
	$(CC) $(DBG)  $(OPT) symmHelfLJSJ.cpp Vector.o Matrix.o Problem.o SymmReduced.o MultiMin.o MultiRoot.o Group.o MembLJ.o -L./ -lCont -lLBFGS -lgfortran -lgsl -lgslcblas -o symmHelfLJSJ.out

ReadData: ReadData.cpp Vector.o Matrix.o Problem.o MembLJ.o
	$(CC) $(DBG)  $(OPT) ReadData.cpp Vector.o Matrix.o Problem.o MembLJ.o -lgfortran -lgsl -lgslcblas -o ReadData.out

contHelfLJ: contHelfLJ.cpp libCont
	$(CC) contHelfLJ.cpp -L ./ -lCont -lgsl -lgslcblas -o test.out $(DBG) -O2

SymmReduced.o: SymmReduced.cpp
	$(CC) $(DBG) $(OPT) -c SymmReduced.cpp 

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

libCont: continuation.o kbhit.o
	ar cru libCont.a continuation.o kbhit.o
	ranlib libCont.a

continuation.o: continuation.h continuation.cpp
	$(CC) -c continuation.h continuation.cpp -lgsl -lgslcblas $(DBG) $(OPT)

kbhit.o: kbhit.h kbhit.c
	$(CC) -c kbhit.h kbhit.c $(DBG) $(OPT)

clean:
	rm -r *.o *.gch *.vtk *.out *.a *.txt *.dSYM
