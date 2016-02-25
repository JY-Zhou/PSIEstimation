#!/bin/bash
rm -rf readsFlux.pro readsFlux.fasta readsFlux.bed readsFlux.lib isoformSim_sorted.gtf
./flux-simulator -p readsFlux.par

python3 ComputeFluxGroundTruth.py
echo 'Grounth Truth are calculated too..'
