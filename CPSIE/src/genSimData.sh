#!/bin/bash
python3 ../generator/GeneSetting.py $1
python3 ../generator/GenerationPipeline.py
cd ../kits/
./SimulateReads.sh
