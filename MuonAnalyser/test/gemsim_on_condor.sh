#!/bin/sh

# Steps from "runTheMatrix.py -w upgrade -l 10811.0"

cd /cms/scratch/iwatson/GEM/CMSSW_10_3_0_pre3_vanilla/src/GEMSimulation
eval `scramv1 runtime -sh`

cmsDriver.py TenMuExtendedE_0_200_pythia8_cfi --python_filename=python/step1_$1.py  --conditions auto:phase1_2018_realistic -n 500 --era Run2_2018 --eventcontent FEVTDEBUG --relval 10000,100 -s GEN,SIM --datatier GEN-SIM --beamspot Realistic25ns13TeVEarly2018Collision --geometry DB:Extended --fileout root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/iawatson/GemSimulation/step1/$1.root   > log/step1_$1.log  2>&1

cmsDriver.py step2  --python_filename=python/step2_$1.py  --conditions auto:phase1_2018_realistic -s DIGI:pdigi_valid,L1,DIGI2RAW,HLT:@relval2018 --datatier GEN-SIM-DIGI-RAW -n -1 --geometry DB:Extended --era Run2_2018 --eventcontent FEVTDEBUGHLT --filein root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/iawatson/GemSimulation/step1/$1.root  --fileout root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/iawatson/GemSimulation/step2/$1.root  > log/step2_$1.log  2>&1

