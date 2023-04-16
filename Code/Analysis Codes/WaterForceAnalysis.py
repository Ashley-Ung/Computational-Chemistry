#!/usr/bin/python
"""

WaterForceAnalysis.py

	Purpose is to read frame data from a simulation and generate an output file with water force vs depth data readable by excel:

		1. Check for necessary data files
		2. parse and read simulation info and topology files
		3. setup for density analysis
			a. Parse atom list for HOH molecules
			b. Create a numpy array of those molecules iwht column 0  O atom index, Column 1 H1 atom index, and column 2 H2 atom index
			
		4. The data is kept in a numpy array of rows determined by z value of O atom position;
			Columns are sum of x force for the molecule, sum of Y force, sum of z force, and sum of magnitude of force, and sum of magnitude in x-y plane.
			
			
		5. Iterate through the binary frame data
			a. in paralell for iteate through all water molecules listed in the water array, and add to z bin 
		5. Divide by frames to get average	
		
		6. Save as labeled space delimited ascii file
		
		
		
"""


import sys

import os.path
# String regular expressions
import re
# import numpy for extensive array and math functions
import numpy as np
#f multiprocessing
import multiprocessing as mp
# necessary for data sharing in multiprocessing
import ctypes as c
# Need abs math function
import math

# Garbage Collector!
import gc

import FrameDataClass as fd

# First check for presence of expected data files		
		
def main():
	
# Check the imput line
	if len(sys.argv)!=3:
		print("Two input arguments are expected: The input file for the simulations, and the integer number of the simulation.")
		sys.exit(1)
		
	basefile = sys.argv[1]
	
	simnumber= int(sys.argv[2])
		
# First check for presence of epected data files
	if not os.path.isfile(basefile+".siminfo"):
		print("Simulation info file {} does not exist.  Exiting!".format(basefile+".siminfo"))
		sys.exit(1)
	if not os.path.isfile(basefile+".simtopo"):
		print("Simulation info file {} does not exist.  Exiting!".format(basefile+".simtopo"))
		sys.exit(1)

	if not os.path.isfile("{}_{:02d}_traj.framebin".format(basefile,simnumber)):
		print("Simulation info file {}_{:02d}_traj.framebin does not exist.  Exiting!".format(basefile,simnumber))
		sys.exit(1)

# This reads both simulation data and topology into structured data
	global sim_data
	sim_data=fd.FrameData(basefile,simnumber)

# Output relevant simulation data	
	print('Simulation dimensions {} x {} x {} {}'.format(sim_data.simdimensions.x, sim_data.simdimensions.y, sim_data.simdimensions.z, sim_data.simdimensions.units))
	
	print('Frame count ', sim_data.framecount)


#  First set up the z bins for count data

# Editable Constant window +/- zero

	global z_window
	global z_bins
	z_window = 2.0  # 0 +/- this distamce nm units
	z_bins = 201  # Odd value, centered at 0
	
#Calculated Constants  (Calculate once only here! and not in function)
	global bin_width
	bin_width = 2.*z_window / (z_bins-1.)
	global z_limit
	z_limit = float(z_window)+bin_width/2
	global z_multiplier
	z_multiplier =(z_bins-1.)/z_window/2.
###

# Now parse the topology for the atom types.  
	water_lookup = Parse_waters(sim_data.atominfo,sim_data.atomcount)
	water_count = np.shape(water_lookup)[0]
	print('Analyzing forces for {} water molecules'.format(water_count))
	
	
# Create counting bins for z_bins rows and atom_type_count columns.  This requires the array be created in the parallel processing memory
	global force_bins
	
	mp_arr = mp.Array('d',z_bins*5)  # initializes as zero
	np_arr = np.frombuffer(mp_arr.get_obj()) #,dtype="float64"
	force_bins=np_arr.reshape(z_bins,5)
	
# This loops through the frames	
	for framecounter in range(sim_data.framecount):
		print('Analyzing frame {} of {}'.format(framecounter+1,sim_data.framecount))

# Load the next froame of binary data
		sim_data.next()

# Set up parallel processing  
		pool = mp.Pool(mp.cpu_count())
		pool.map(assign_z_bin, [ row for row in water_lookup])
		pool.close()
		pool.join()
		
	print('Frame analysis complete')
		
		
# Next major step is to take count data to averaged density
	
# Save the force Output Files
	print('Writing Force Output File')
	write_Force_Results("{}_{:02d}_force_profile.txt".format(basefile,simnumber), force_bins)

	return

	
# END MAIN
###############################################################################

#######
# Function returns a dictionary with a tuple key (res type, atom type) and data tuple (index, atomic number)
def Parse_waters(atom_list,atom_count):

	water_dict = {}
	watercount = 0
	for i in range(atom_count):
		if atom_list[i]['molecule'] == b'HOH':
			if not (atom_list[i]['moleculeindex'] in water_dict):
				water_dict[atom_list[i]['moleculeindex']]=watercount
				watercount +=1
	
	water_list=np.zeros([watercount,3],int)
	
	for i in range(atom_count):
		if (atom_list[i]['moleculeindex'] in water_dict):
			j = water_dict[atom_list[i]['moleculeindex']]
			if atom_list[i]['type']==b'O':
				water_list[j,0]=i
			elif atom_list[i]['type']==b'H1':
				water_list[j,1]=i
			elif atom_list[i]['type']==b'H2':
				water_list[j,2]=i
			
	return water_list
#######



# Define the inner loop function to execute in parallel
# this takes the row of the water list array and adds all the forces to the z-slab in which the O atom is located
# Sums all the x, y, and z forces, plus the force magnitude in the x-y plane

def assign_z_bin(water):
	global force_bins
	z_pos=sim_data.positions[water[0],2]
# Checkif position is in analysis ragnge
	if abs(z_pos) <= z_limit:
# Increment count in zbin , atom type bin
		bin_index=int(round((z_pos+z_window)*z_multiplier))
		
		x_force=sim_data.forces[water[0],0]+sim_data.forces[water[1],0]+sim_data.forces[water[2],0]
		y_force=sim_data.forces[water[0],1]+sim_data.forces[water[1],1]+sim_data.forces[water[2],1]
		z_force=sim_data.forces[water[0],2]+sim_data.forces[water[1],2]+sim_data.forces[water[2],2]
		force_bins[bin_index, 1] += x_force
		force_bins[bin_index, 2] += y_force
		force_bins[bin_index, 3] += z_force
		force_bins[bin_index, 4] += math.sqrt(x_force**2 + y_force**2)
		force_bins[bin_index, 0] +=1
	return
######

#######

# generate z-labels for output
def z_label(index):
	return (float(index) - float(z_bins-1.)/2)* bin_width

## end label function


#######
# Format and write force data
# Data output format:
#	first row is labels
#	First Column is z-depth
#	Columns are molecules for ions and solvents, and fragments of surfactants
def write_Force_Results(filename,fbins):
	
# Average the per molecule force, and replace the counter in the 0 column with the z value

	for b in range(z_bins):
		if fbins[b,0]>0:
			fbins[b,1:5]=fbins[b,1:5]/fbins[b,0]
		fbins[b,0]=z_label(b)
		
	
		
# Create Header string
	hdr='nm\tx_force\ty_force\tz_force\tmagnitude_x-y_force'
		
	np.savetxt(filename, fbins, delimiter='\t', newline='\n', header=hdr)

	return
########		

if __name__ == '__main__':  
	mp.set_start_method('fork')
	main()						
	
