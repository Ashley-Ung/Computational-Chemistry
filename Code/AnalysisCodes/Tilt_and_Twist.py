"""
Twist_Tilt.py
Kevin Johnson Pacific University
Summer 2022



This Program will read in all of the atom data from a full set of repeat simulations simulations using FrameDataClass.py and will:
	• determine the surfactant present.  Assumes only one, DDC, DDS, DBC, or DBS
	• Creates an array of atom numbers for tilt and twist analysis. Tilt: cos(angle of axis to interface) Twist cos(angle of HG plane to the interface).
	• Histogram data collected in one 3-D array: first dimension is depth defined by C or S HG atom, second the twist bins, and third the tilt bins.
	• prints integer count for full frame analysis.  i.e. (# fraames)x the probability

Multiprocessing not implemented: There's a max of 30 surfactant molecules per frame.  Better use would be to run myltiple instances of this python code for multiple simulations.



"""
import sys

import os.path

# import numpy for extensive array and math functions
import numpy as np


import FrameDataClass as fd



import math
		
def main():
	print("Welcome To Twist and Tilt Analysis")
	
#This is the size the analysis bins
	global depthBins
	depthBins = 12
	global depthMin
	depthMin = -1.0 # nm
	global depthMax
	depthMax =  0.5 #nm

# Calculate the function variables to determine the depth bin number
# formula is int((z - depthMin)/(depthMax-depthMin)*depthBins
# note that int is the same as floor() and only rounds down

	global depthConst
	depthConst = depthBins/(depthMax-depthMin)
	
	global tiltBins
	tiltBins = 40
	global twistBins
	twistBins = 20
	
# Create an integer array of zeros for binning
	global TTBins
	TTBins=np.zeros([depthBins,twistBins,tiltBins],dtype=np.int16)
	

# Check for presence of two parameters in input line
	global basefile
	basefile=''
	for i in range(len(sys.argv)-1):
		basefile=basefile+sys.argv[i+1]+' '
	basefile=basefile[0:len(basefile)-1]
	if basefile[len(basefile)-4:]==".pdb" or basefile[len(basefile)-4:]==".chk": basefile=basefile[0:len(basefile)-4]

# First check for presence of epected data files
	if not os.path.isfile(basefile+".siminfo"):
		print("Simulation info file {} does not exist.  Exiting!".format(basefile+".siminfo"))
		sys.exit(1)
	if not os.path.isfile(basefile+".simtopo"):
		print("Simulation info file {} does not exist.  Exiting!".format(basefile+".simtopo"))
		sys.exit(1)

# Count the nmber of simulations present
		
	global simCount
	simCount=0
	while os.path.isfile("{}_{:02d}{}".format(basefile,simCount+1,"_traj.framebin")):
		simCount+=1
		
	if simCount==0:
		print("No trajectory binary files found")
		sys.exit(1)
		

	
	global sim_data
#Read the data for the simulation
	sim_data=fd.FrameData(basefile,1)
	
	global frameCount	
	frameCount=sim_data.framecount

#The number of atoms in the simulation
	global atomCount
	atomCount = sim_data.atomcount

	
#List of Surfactant Molecu;e labels
	
	global interface_molecules
	surfactants = [b'DDC',b'DDS',b'DBC',b'DBS']
	
# Atom names for analysis in order: depth key, tilt1, tilt2, twist1, twist2, and for DDS needed twist3
	
	global analysisAtoms
	analysisAtoms = {
	b'DDC': (b'C12',b'C11',b'C12',b'O1',b'O2'),
	b'DDS' : (b'S1',b'C12',b'S1',b'O1',b'O2',b'O3'),
	b'DBC' : (b'C19',b'C13',b'C18',b'C16',b'C17'),
	b'DBS' : (b'S1',b'C13',b'C18',b'C16',b'C17')
	}
	


#Fills List with atomic number from atom index
	print("Creating Atom Listing")
	global SurfactantAtoms
	global SurfacantCount
	SurfacantCount = 0
# Creates list of atomic numbers to analyze for each surfactant
# Parse all atoms. When a Surfactant is identified then the atoms are sequentially checked for the key atoms for the twost and tilt analysis, the atom numbers are then saved to an array
	atom=0
	while atom <atomCount:
		molecule= sim_data.atominfo[atom]['molecule']
		moleculeIndex=sim_data.atominfo[atom]['moleculeindex']
		if molecule not in surfactants:
			atom+=1
		else:
			if SurfacantCount == 0:
# Initialize array
				SurfactantAtoms = np.zeros([1,len(analysisAtoms[molecule])],dtype=np.int16)
			else:
				SurfactantAtoms=np.append(SurfactantAtoms,np.zeros([1,len(analysisAtoms[molecule])],dtype=np.int16),axis=0)
# Now parse the next atoms in the array to find the atoms used for tilt/twist alculations			
			while atom < atomCount and sim_data.atominfo[atom]['moleculeindex'] == moleculeIndex:
				if sim_data.atominfo[atom]['type'] in analysisAtoms[molecule]:
					for i in range(len(analysisAtoms[molecule])):
						if sim_data.atominfo[atom]['type']==analysisAtoms[molecule][i]:
							SurfactantAtoms[SurfacantCount,i ]= atom
				atom+=1
			SurfacantCount+=1
		
	print(SurfacantCount," Surfactants detected and parsed")
	sim_data.done()

	print("Starting Twist/Tilt Analysis.")


# Simulation Loop
	for sim in range(simCount):
		print("Starting analysis for sim # {:d}".format(sim+1))
		sim_data=fd.FrameData(basefile,sim+1)
	
#frame loop	
		for frame in range(frameCount):
			sim_data.next()
			for s in range(SurfacantCount):
				Assign_Bin(SurfactantAtoms[s,:])
				
		sim_data.done()
	
	
	print("Writing Data File")
	write_TwistTilt()
	
	
	print("Twist/Tilt Analysis Complete for"+basefile)	
	
	return
### End Main

#######	
# This function calculates the twist and tilt quantities and increments the corresponding bin based on the HG depth
def Assign_Bin(atomArray):

	z_bin = int((sim_data.positions[atomArray[0],2]-depthMin)*depthConst)
	if z_bin<0:
		print ("Warning: Surfactant HG is below of bin range")
		z_bin=0
	if z_bin>depthBins-1:
		print ("Warning: Surfactant HG is above of bin range")
		z_bin=depthBins-1

# define unit vector perpendicular to the interface
	z_ref=np.array([0.0,0.0,1.0])

	tiltvec=sim_data.positions[atomArray[2],:]-sim_data.positions[atomArray[1],:]
	tilt=np.dot(tiltvec,z_ref)/np.linalg.norm(tiltvec)

	tilt_bin=max(0,min(tiltBins-1,int((-tilt+1)/2*tiltBins)))

	
# the twist vector calculation depends on the reference

# For DBC and DBS the reference is the benzene ring. For DDC the reference is the two O atoms
	if len(atomArray) == 5:
		twistvec=sim_data.positions[atomArray[4],:]-sim_data.positions[atomArray[3],:]
	else:
# This is the case of DDS, where the tilt is determione by the O atom closest to the water layer -- lowest z value.
		if sim_data.positions[atomArray[5],2] <  sim_data.positions[atomArray[4],2]:
			if sim_data.positions[atomArray[5],2] <  sim_data.positions[atomArray[3],2]:
				twistvec=sim_data.positions[atomArray[5],:]-(sim_data.positions[atomArray[4],:]+sim_data.positions[atomArray[3],:])/2
			else:
				twistvec=sim_data.positions[atomArray[3],:]-(sim_data.positions[atomArray[4],:]+sim_data.positions[atomArray[5],:])/2
		else:
			twistvec=sim_data.positions[atomArray[4],:]-(sim_data.positions[atomArray[5],:]+sim_data.positions[atomArray[3],:])/2
	twist=abs(np.dot(twistvec,z_ref)/np.linalg.norm(twistvec))
	
	twist_bin=max(0,min(twistBins-1,int(twist*twistBins)))

	TTBins[z_bin,twist_bin,tilt_bin] +=1
						
	return True
#####



def write_TwistTilt():
	
	f=open(basefile+"_TwistTilt.txt",'w')
	
		
	f.write(basefile+"\n")
	
	zbin_increment = (depthMax-depthMin)/depthBins
	normbins=TTBins/frameCount/simCount
	f.write("{:.3f}\n".format(np.amax(normbins)))
	
	for counter in range(depthBins):
		f.write("{:.3f} – {:.3f} nm\n".format(depthMin+counter*zbin_increment,depthMin+(counter+1)*zbin_increment))
		for i in range(tiltBins):
			for j in range(twistBins):
				f.write("{:3f}\t".format(normbins[counter,j,i]))
			f.write("\n")
	f.close()
	
	return True
	
	


	
if __name__ == '__main__':  
		main()
	