#!/bin/bash

OUTPUT=call.prof
CALLGRAPH=callGraph.png

valgrind --tool=callgrind --callgrind-out-file=${OUTPUT} ./SolverPipeline.exe \
    -g ../input/genome.fa \
    -a ../input/exonBoundary.bed \
    -r ../input/reads.fa \
    -k 30 \
    -o ./PsiResult.json \
    2> profiling.out
./gprof2dot.py -f callgrind ${OUTPUT} | dot -Tpng -o ${CALLGRAPH}
eog ${CALLGRAPH}
