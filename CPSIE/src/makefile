objects=SolverPipeline.o KmerHash.o EMAlgorithm.o SLSQP.o
boost=/cygdrive/c/Boost/boost_1_59_0/
Eigen=/cygdrive/c/Eigen/

SolverPipeline.exe: $(objects)
	g++ -O3 -lgfortran -std=c++11 $(objects) -o SolverPipeline.exe 
	
SolverPipeline.o: SolverPipeline.cpp EMAlgorithm.cpp EMAlgorithm.h \
    KmerHash.cpp KmerHash.h SLSQP.f
	g++ -O3 -std=c++11 -I $(boost) -I $(Eigen) -c SolverPipeline.cpp

EMAlgorithm.o: EMAlgorithm.cpp EMAlgorithm.h KmerHash.cpp KmerHash.h \
    SLSQP.f
	g++ -O3 -std=c++11 -I $(boost) -I $(Eigen) -c EMAlgorithm.cpp

KmerHash.o: KmerHash.cpp KmerHash.h
	g++ -O3 -std=c++11 -I $(boost) -c KmerHash.cpp

SLSQP.o: SLSQP.f
	gfortran -O3 -c -w SLSQP.f

clean:
	rm -f $(objects) SolverPipeline.exe
