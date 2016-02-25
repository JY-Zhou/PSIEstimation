#!/bin/bash
set -e

make

echo ''

GENOME=../input/chrSim.fa
EXONBND=../input/exonBoundarySim.bed
READ=../kits/readsFlux.fasta
GROUNDTRUTH=../kits/PsiGroundTruthFlux.json
#READ=../input/readsSim.fa
#GROUNDTRUTH=../kits/PsiGroundTruthSim.json
OUTPUT=./PsiResult.json

time ./SolverPipeline.exe -g ${GENOME} -a ${EXONBND} -r ${READ} -o ${OUTPUT} -k 30

echo ''
echo ''
echo --- RESULT CHECKER ---
python3 ./checkPsi.py ${GROUNDTRUTH} ${OUTPUT}
