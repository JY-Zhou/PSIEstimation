#!/bin/bash

OUTPUT=call.prof
CALLGRAPH=callGraph.png

valgrind --tool=callgrind --callgrind-out-file=${OUTPUT} ./SolverPipeline.exe
./gprof2dot.py -f callgrind ${OUTPUT} | dot -Tpng -o ${CALLGRAPH}
eog ${CALLGRAPH}
