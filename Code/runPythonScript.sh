#!/bin/bash

# Ashley Ung 
# runPythonScript.sh
# This shell program will speed up the process of running the python program for analysis. The shell script will 
# loop through the different saturations for analysis.
# Locate the folder where this script is found first. In this case, "WorkingCode" folder. 
# To make this script an executable use the following command: chmod +x runPythonScript.sh
# To run the script, use the following command: ./runPythonScript.sh

#!/bin/bash

for i in 25 50 75 100 125
do
	for j in {1..4} # The inner loop will run from 1 to 4 for each value of i
	do
		# CCl4_Sims for DDC
		python Interface_Ionic_Analysis_Copy.py /Users/Shared/SurfactantProject/Sim2022/CCl4_Sims/DDC/DDC$i/DDC$i $j

		# CCl4_Sims for DDS
		python Interface_Ionic_Analysis_Copy.py /Users/Shared/SurfactantProject/Sim2022/CCl4_Sims/DDS/DDS$i/DDS$i $j
		
		# CCl4_Sims for DBC
		python Interface_Ionic_Analysis_Copy.py /Users/Shared/SurfactantProject/Sim2022/CCl4_Sims/DBC/DBC$i/DBC$i $j
		
		# CCl4_Sims for DBS
		python Interface_Ionic_Analysis_Copy.py /Users/Shared/SurfactantProject/Sim2022/CCl4_Sims/DBS/DBS$i/DBS$i $j

		# CCl4_SaltSims for NaCl DDC
		python Interface_Ionic_Analysis_Copy.py /Users/Shared/SurfactantProject/Sim2022/CCl4_SaltSims/NaCl_DDC/NaCl_DDC$i/NaCl_DDC$i $j
		
		# CCl4_SaltSims for NaCl DDS
		python Interface_Ionic_Analysis_Copy.py /Users/Shared/SurfactantProject/Sim2022/CCl4_SaltSims/NaCl_DDS/NaCl_DDS$i/NaCl_DDS$i $j
		
		# CCl4_SaltSims for NaCl DBC
		python Interface_Ionic_Analysis_Copy.py /Users/Shared/SurfactantProject/Sim2022/CCl4_SaltSims/NaCl_DBC/NaCl_DBC$i/NaCl_DBC$i $j
		
		# CCl4_SaltSims for NaCl DBS
		python Interface_Ionic_Analysis_Copy.py /Users/Shared/SurfactantProject/Sim2022/CCl4_SaltSims/NaCl_DBS/NaCl_DBS$i/NaCl_DBS$i $j
	done
done
