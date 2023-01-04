"""
Ashley Ung 
runPythonScript.sh

This shell program will speed up the process of running the python program for analysis. The shell script will 
loop through the different saturations for analysis.

To make this script an executable use the following command: chmod +x runPythonScript.sh

To run the script, use the following command: ./runPythonScript.sh
"""

#!/bin/bash
for i in 25 50 75 100 125 
do
	python Interface_Ionic_Analysis_Copy.py /Users/Shared/SurfactantProject/Sim2022/CCl4_Sims/DDC/DDC$i/DDC$i 1
done
