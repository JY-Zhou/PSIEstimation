objects=SolverPipeline.o KmerHash.o EMAlgorithm.o SLSQP.o
#boost=/cygdrive/c/Boost/boost_1_59_0/
#eigen=/cygdrive/c/eigen/
boost=../../../libs/boost_1_59_0/
eigen=../../../libs/eigen/
OptLev=-O3

SolverPipeline.exe: $(objects)
	g++ $(OptLev) -std=c++11 -fopenmp -lgfortran $(objects) -o SolverPipeline.exe -g
	
SolverPipeline.o: SolverPipeline.cpp EMAlgorithm.cpp EMAlgorithm.h KmerHash.cpp KmerHash.h SLSQP.f
	g++ $(OptLev) -std=c++11 -I $(boost) -I $(eigen) -c SolverPipeline.cpp -g

EMAlgorithm.o: EMAlgorithm.cpp EMAlgorithm.h KmerHash.cpp KmerHash.h SLSQP.f
	g++ $(OptLev) -std=c++11 -fopenmp -I $(boost) -I $(eigen) -c EMAlgorithm.cpp -g

KmerHash.o: KmerHash.cpp KmerHash.h
	g++ $(OptLev) -std=c++11 -I $(boost) -c KmerHash.cpp -g

SLSQP.o: SLSQP.f
	gfortran $(OptLev) -c -w SLSQP.f

clean:
	rm -f $(objects) SolverPipeline.exe
