#!/bin/bash -l
#Surfactant simulations
# expects two command line variables: the prefix of the surfactant  i.e. DDC  and the gpu number
# Call: nohup ./processSurfactant.sh basename 1 > ../basename/basenameprocess.log &


source ~/anaconda3/etc/profile.d/conda.sh
conda activate OpenMM77

for var in 25 50 75 100 125
do
	python InterfaceSimulation.py ../$1/$1$var/$1$var 4 $2 > ../$1/$1$var/$1$var"_sim.log"
	python RunFullAnalysis.py ../$1/$1$var/$1$var  > ../$1/$1$var/$1$var"_anal.log"
done
