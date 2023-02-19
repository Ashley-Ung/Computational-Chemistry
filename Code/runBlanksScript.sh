#!/bin/bash

# Ashley Ung 
# runBlanksScript.sh
# This shell program will speed up the process of running the python program for analysis. The shell script will 
# loop through the different blanks for analysis.
# Locate the folder where this script is found first. In this case, "WorkingCode" folder. 
# To make this script an executable use the following command: chmod +x runPythonScript.sh
# To run the script, use the following command: ./runPythonScript.sh

#!/bin/bash

for j in {1..4} # This loop will run from 1 to 4 for each value of j
do
	# Water Cube Blank 
	python Interface_Ionic_Analysis_Copy.py /Users/Shared/SurfactantProject/Sim2022/Blanks/waterCube/waterCube $j
	
	# CCl4 NaCl Blank 
	python Interface_Ionic_Analysis_Copy.py /Users/Shared/SurfactantProject/Sim2022/Blanks/CCl4_NaCl/CCl4_NaCl $j
	
	# NaCl 0.5M Blank
	python Interface_Ionic_Analysis_Copy.py /Users/Shared/SurfactantProject/Sim2022/Blanks/NaCl_0.5M/NaCl_0.5M $j
	
	# Neat Interface Blank
	python Interface_Ionic_Analysis_Copy.py /Users/Shared/SurfactantProject/Sim2022/Blanks/neatInterface/neat $j
done
