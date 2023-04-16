#!/usr/bin/env python3

#!/usr/bin/env python3
"""

Purpose:	Reads interface analysis reqults for a single trajectory and averages. Write average file.
Types:		Mass Density, Charge Density, Force Profile

Interface analysis: For Mass density the water density is fit to a hyperbolic tangent function for each trajectory. fitting parameters and covarience are waved to a file for each tranjectory and for the average

input format: $python Ingerface_Traj_Average.py basefilename

Code will count the number of each type of file with the basefilename followed by file type.

"""
import sys


import os
import subprocess


def main():
	frameextention="_traj.framebin"
	
# Check the input line
	if len(sys.argv)!=2:
		print("One input argument are expected: The input base file name for the simulation. Exiting.")
		sys.exit(1)
		
	basefile = sys.argv[1]
	if basefile[len(basefile)-4:]==".pdb": basefile=basefile[0:len(basefile)-4]
	
# Count the  number of bibnary trajectory files.
	trajcount = 0
	while os.path.isfile("{}_{:02d}{}".format(basefile,trajcount+1,frameextention)):
		trajcount+=1
			
	if trajcount==0:
		print("No trajectory files found. Exiting.")
		sys.exit(1)
			
# Loop through density analysis
		
	for i in range(trajcount):
		if not os.path.isfile("{}_{:02d}_mass_density.txt".format(basefile,i+1)):
			print("Executing Density analysis for file: {}_{:02d}{}".format(basefile,i+1,frameextention))
			subprocess.call("python3 MassAndChargeDensityAnalysis.py {} {:d}".format(basefile,i+1), shell=True)
		
# Loop through water force analysis
		
	for i in range(trajcount):
		if not os.path.isfile("{}_{:02d}_force_profile.txt".format(basefile,i+1)):
			print("Executing Water Force analysis for file: {}_{:02d}{}".format(basefile,i+1,frameextention))
			subprocess.call("python3 WaterForceAnalysis.py {} {:d}".format(basefile,i+1), shell=True)
		
				
# Loop Through Interface Analysis
		
	for i in range(trajcount):
		if not os.path.isfile("{}_{:02d}_interface_comp.txt".format(basefile,i+1)):
			print("Executing Interface analysis for file: {}_{:02d}{}".format(basefile,i+1,frameextention))
			subprocess.call("python3 Interface_Ionic_Analysis.py {} {:d}".format(basefile,i+1), shell=True)

# Execute the twist-tilt analysis.  This internally averages all simulations in the local directory
		print(basefile)
		subprocess.call("python3 Tilt_and_Twist.py {}".format(basefile), shell=True)	
		subprocess.call("python3 Plot_Twist_Tilt.py {}_TwistTilt.txt".format(basefile), shell=True)
		
# Now average 
	print("Averaging all analysis files for {}".format(basefile))
	subprocess.call("python3 Traj_Averager.py {}".format(basefile), shell=True)
	print("All Done!")	
		
	return True


#########

if __name__ == '__main__':  
	main()						
	
	
