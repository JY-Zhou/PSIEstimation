objects=SolverPipeline.o KmerHash.o EMAlgorithm.o
boost=/cygdrive/c/Boost/boost_1_59_0/

SolverPipeline.exe: $(objects)
	g++ -O2 -std=c++11 $(objects) -o SolverPipeline.exe
	
SolverPipeline.o: SolverPipeline.cpp
	g++ -O2 -std=c++11 -I $(boost) -c SolverPipeline.cpp

KmerHash.o: KmerHash.cpp KmerHash.h
	g++ -O2 -std=c++11 -I $(boost) -c KmerHash.cpp

EMAlgorithm.o: EMAlgorithm.cpp EMAlgorithm.h
	g++ -O2 -std=c++11 -I $(boost) -c EMAlgorithm.cpp

clean:
	rm -f $(objects) SolverPipeline.exe
