#!/usr/bin/python
"""

DensityAnalysis.py
	Purpose is to read frame data from a simulation and generate an output file with density / depth data readable by excel:

		1. Check for necessary data files
		2. parse and read simulation info and topology files
		3. setup for density analysis
			a. Define atom bins amd bin labels based on the topology atom types and residues.  THIS FUNCTION must be edited for new molecule types.  Also create the bin mass multiplier
			b. Create an atom bin assignment ARRAY that will be the second index of the bin
			c create a function that takes z value and converts to bin matrix.  Additionally a list of the bin z values/labels
			
		4. Iterate through the binary frame data
			a. in paralell for iteate through all atoms and increment appropriate mass bins
			
		5. take int counting array and in float array average the counts.  Multiply by bin mass
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

# First check for presence of epected data files		
		
def main():
	
# Check the imput line
	if len(sys.argv)!=3:
		print("Two input arguments are expected: The input base file name for the simulations, and the integer number of the simulation.")
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
	print('Atom count ', sim_data.atomcount)
#	print('atom info:\n',sim_data.atominfo)

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
	atom_lookup = Parse_Atoms(sim_data.atominfo)
	atom_type_count = len(atom_lookup) 
	print('Analyzing for {} atom types'.format(atom_type_count))
	
# Create the lookup array
	global atom_index
	atom_index = create_atom_bins(sim_data.atomcount,sim_data.atominfo,atom_lookup)
	
# Create counting bins for z_bins rows and atom_type_count columns.  This requires the array be created in the parallel processing memory
	global counting_bins
	
	mp_arr = mp.Array('i',z_bins*atom_type_count)  # initializes as zero
	np_arr = np.frombuffer(mp_arr.get_obj(),dtype="int32")
	counting_bins=np_arr.reshape(z_bins,atom_type_count)
	


# This loops through the frames	
	for framecounter in range(sim_data.framecount):
		print('Analyzing frame {} of {}'.format(framecounter+1,sim_data.framecount))

# Load the next froame of binary data
		sim_data.next()

# Set up parallel processing  
		pool = mp.Pool(mp.cpu_count())
		pool.map(assign_z_bin, [ a for a in range(sim_data.atomcount)])
		pool.close()
		pool.join()
		
	print('Frame analysis complete')
	sim_data.done()
	print('Count of assigned atoms = ',np.sum(counting_bins))
		
		
# Next major step is to take count data to averaged density

	print('Calculating density profiles')
	mass_density_bins = calculateMassDensity(counting_bins, atom_lookup, sim_data.framecount)
	print('Calculating charge density profiles')
	charge_density_bins = calculateChargeDensity(counting_bins, atom_lookup, sim_data.framecount)
	
	
# Save the Mass Output Files
	print('Writing density Output Files')
	write_Mass_Density_Results("{}_{:02d}_mass_density.txt".format(basefile,simnumber), mass_density_bins, atom_lookup, sim_data.atomcount)

	write_Charge_Density_Results("{}_{:02d}_charge_density.txt".format(basefile,simnumber), charge_density_bins, atom_lookup, sim_data.atomcount)

	return
	
# END MAIN
###############################################################################

#######
# Function returns a dictionary with a tuple key (res type, atom type) and data tuple (index, atomic number)
def Parse_Atoms(atom_list):
	mydict={}
	count = 0
	for atom in atom_list:
		if not((atom['molecule'],atom['type']) in mydict):
			mydict[(atom['molecule'],atom['type'])]=(count,atom['atomicno'])
			count +=1
	return mydict
#######

######
# Creates the numpy array that by index returns the index of atom type column of the data matrix
def create_atom_bins(count,atom_lookup,atom_dict):
	array= np.empty(count,int)
	index=0
	for atom in atom_lookup:
		array[index]=atom_dict[atom['molecule'],atom['type']][0]
		index +=1
	return array
#######


# Define the inner loop function to execute in parallel

def assign_z_bin(a_index):
	z_pos=sim_data.positions[a_index,2]
# Checkif position is in analysis ragnge
	if abs(z_pos) <= z_limit:
# Increment count in zbin , atom type bin
		counting_bins[int(round((z_pos+z_window)*z_multiplier)), atom_index[a_index]] +=1
	return
######

#######
# Function returns a discionary with a tuple key (res type, atom type) and data tuple (index, molar mass)

# generate z-labels for output
def z_label(index):
	return (float(index) - float(z_bins-1.)/2)* bin_width

## end label function

#######
# Function to convert count dta to density data
def calculateMassDensity(bins, atomlist, fcount):
	
# Data structure is float array same size as bins
	density = np.zeros(bins.shape,float)
	
#Calculate the multiplier for atom count to mole count per volume volume of a slab in the z direction converts to moles per cm^3
	slab_multiplier	= (1.6606e-3)/(bin_width * sim_data.simdimensions.x * sim_data.simdimensions.y )
# Define the dictionary with molar masses

	molarmasses = {
		0 : 0.0,
		1 : 1.00794,
		3 : 6.938,
		6 : 12.011,
		7 : 14.0067,
		8 : 15.9994,
		9 : 18.99840,
		11 : 22.98977,
		12 : 24.305,
		15 : 30.97376,
		16 : 32.066,
		17 : 35.4527,
		19 : 39.0983,
		20 : 40.078,
		30 : 65.38,
		35 : 79.904,
		53 : 126.90
	}


# create the molar mass multiplier array
	mm_mult = np.empty(len(atomlist),float)
	for a in atomlist:
		index=atomlist[a][0]
		mass=molarmasses[atomlist[a][1]]
		mm_mult[index]=mass
	
			
# Calculate the atom type density bins

	density = (bins*mm_mult)*slab_multiplier/fcount

	return density
#######

#######
# Format and write mass data
# Data output format:
#	first row is labels
#	First Column is z-depth
#	Columns are molecules for ions and solvents, and fragments of surfactants
def write_Mass_Density_Results(filename,dbins,atom_lookup,atomcount):
	
# Dictionary listing that will assign molecule,atom pairs to a column with label
	column_dict = {
		(b'HOH',b'H1') : 'H2O',
		(b'HOH',b'H2') : 'H2O',
		(b'HOH',b'O')  : 'H2O',
		
		(b'LI',b'Li')  : 'Li+',
		(b'NA',b'Na')  : 'Na+',
		(b'MG',b'Mg')  : 'Mg2+',
		(b'K',b'K')  : 'K+',
		(b'CA',b'Ca')  : 'Ca2+',
		(b'ZN',b'Zn')  : 'Zn2+',
		
		(b'F',b'F')  : 'F-',
		(b'Cl',b'Cl')  : 'Cl-',
		(b'Br',b'Br')  : 'Br-',
		(b'I',b'I')  : 'I-',
		
		(b'CTC',b'Cl1'): 'CCl4',
		(b'CTC',b'Cl2'): 'CCl4',
		(b'CTC',b'Cl3'): 'CCl4',
		(b'CTC',b'Cl4'): 'CCl4',
		(b'CTC',b'C')  : 'CCl4',
		
		(b'DDC',b'O2')  : 'DDC.HG',
		(b'DDC',b'O1')  : 'DDC.HG',
		(b'DDC',b'C12') : 'DDC.HG',
		(b'DDC',b'C11') : 'DDC.C2',
		(b'DDC',b'H22') : 'DDC.C2',
		(b'DDC',b'H23') : 'DDC.C2',
		(b'DDC',b'C10') : 'DDC.C3',
		(b'DDC',b'H20') : 'DDC.C3',
		(b'DDC',b'H21') : 'DDC.C3',
		(b'DDC',b'C9')  : 'DDC.C4',
		(b'DDC',b'H19') : 'DDC.C4',
		(b'DDC',b'H18') : 'DDC.C4',
		(b'DDC',b'C8')  : 'DDC.C5',
		(b'DDC',b'H16') : 'DDC.C5',
		(b'DDC',b'H17') : 'DDC.C5',
		(b'DDC',b'C7')  : 'DDC.C6',
		(b'DDC',b'H14') : 'DDC.C6',
		(b'DDC',b'H15') : 'DDC.C6',
		(b'DDC',b'C6')  : 'DDC.C7',
		(b'DDC',b'H12') : 'DDC.C7',
		(b'DDC',b'H13') : 'DDC.C7',
		(b'DDC',b'C5')  : 'DDC.C8',
		(b'DDC',b'H10') : 'DDC.C8',
		(b'DDC',b'H11') : 'DDC.C8',
		(b'DDC',b'C4')  : 'DDC.C9',
		(b'DDC',b'H8')  : 'DDC.C9',
		(b'DDC',b'H9')  : 'DDC.C9',
		(b'DDC',b'C3')  : 'DDC.C10',
		(b'DDC',b'H6')  : 'DDC.C10',
		(b'DDC',b'H7')  : 'DDC.C10',
		(b'DDC',b'C2')  : 'DDC.C11',
		(b'DDC',b'H4')  : 'DDC.C11',
		(b'DDC',b'H5')  : 'DDC.C11',
		(b'DDC',b'C1')  : 'DDC.C12',
		(b'DDC',b'H1')  : 'DDC.C12',
		(b'DDC',b'H2')  : 'DDC.C12',
		(b'DDC',b'H3')  : 'DDC.C12',
		
		(b'DDS',b'O3')  : 'DDS.HG',
		(b'DDS',b'O2')  : 'DDS.HG',
		(b'DDS',b'O1')  : 'DDS.HG',
		(b'DDS',b'S1')  : 'DDS.HG',
		(b'DDS',b'C12') : 'DDS.C1',
		(b'DDS',b'H25') : 'DDS.C1',
		(b'DDS',b'H24') : 'DDS.C1',
		(b'DDS',b'C11') : 'DDS.C2',
		(b'DDS',b'H22') : 'DDS.C2',
		(b'DDS',b'H23') : 'DDS.C2',
		(b'DDS',b'C10') : 'DDS.C3',
		(b'DDS',b'H20') : 'DDS.C3',
		(b'DDS',b'H21') : 'DDS.C3',
		(b'DDS',b'C9')  : 'DDS.C4',
		(b'DDS',b'H19') : 'DDS.C4',
		(b'DDS',b'H18') : 'DDS.C4',
		(b'DDS',b'C8')  : 'DDS.C5',
		(b'DDS',b'H16') : 'DDS.C5',
		(b'DDS',b'H17') : 'DDS.C5',
		(b'DDS',b'C7')  : 'DDS.C6',
		(b'DDS',b'H14') : 'DDS.C6',
		(b'DDS',b'H15') : 'DDS.C6',
		(b'DDS',b'C6')  : 'DDS.C7',
		(b'DDS',b'H12') : 'DDS.C7',
		(b'DDS',b'H13') : 'DDS.C7',
		(b'DDS',b'C5')  : 'DDS.C8',
		(b'DDS',b'H10') : 'DDS.C8',
		(b'DDS',b'H11') : 'DDS.C8',
		(b'DDS',b'C4')  : 'DDS.C9',
		(b'DDS',b'H8')  : 'DDS.C9',
		(b'DDS',b'H9')  : 'DDS.C9',
		(b'DDS',b'C3')  : 'DDS.C10',
		(b'DDS',b'H6')  : 'DDS.C10',
		(b'DDS',b'H7')  : 'DDS.C10',
		(b'DDS',b'C2')  : 'DDS.C11',
		(b'DDS',b'H4')  : 'DDS.C11',
		(b'DDS',b'H5')  : 'DDS.C11',
		(b'DDS',b'C1')  : 'DDS.C12',
		(b'DDS',b'H1')  : 'DDS.C12',
		(b'DDS',b'H2')  : 'DDS.C12',
		(b'DDS',b'H3')  : 'DDS.C12',
		
		(b'DBC',b'O2')  : 'DBC.HG',
		(b'DBC',b'O1')  : 'DBC.HG',
		(b'DBC',b'C19') : 'DBC.HG',
		(b'DBC',b'C18') : 'DBC.AR',
		(b'DBC',b'C17') : 'DBC.AR',
		(b'DBC',b'H29') : 'DBC.AR',
		(b'DBC',b'C16') : 'DBC.AR',
		(b'DBC',b'H28') : 'DBC.AR',
		(b'DBC',b'C15') : 'DBC.AR',
		(b'DBC',b'H27') : 'DBC.AR',
		(b'DBC',b'C14') : 'DBC.AR',
		(b'DBC',b'H26') : 'DBC.AR',
		(b'DBC',b'C13') : 'DBC.AR',
		(b'DBC',b'C12') : 'DBC.C1',
		(b'DBC',b'H24') : 'DBC.C1',
		(b'DBC',b'H25') : 'DBC.C1',
		(b'DBC',b'C11') : 'DBC.C2',
		(b'DBC',b'H22') : 'DBC.C2',
		(b'DBC',b'H23') : 'DBC.C2',
		(b'DBC',b'C10') : 'DBC.C3',
		(b'DBC',b'H20') : 'DBC.C3',
		(b'DBC',b'H21') : 'DBC.C3',
		(b'DBC',b'C9')  : 'DBC.C4',
		(b'DBC',b'H19') : 'DBC.C4',
		(b'DBC',b'H18') : 'DBC.C4',
		(b'DBC',b'C8')  : 'DBC.C5',
		(b'DBC',b'H16') : 'DBC.C5',
		(b'DBC',b'H17') : 'DBC.C5',
		(b'DBC',b'C7')  : 'DBC.C6',
		(b'DBC',b'H14') : 'DBC.C6',
		(b'DBC',b'H15') : 'DBC.C6',
		(b'DBC',b'C6')  : 'DBC.C7',
		(b'DBC',b'H12') : 'DBC.C7',
		(b'DBC',b'H13') : 'DBC.C7',
		(b'DBC',b'C5')  : 'DBC.C8',
		(b'DBC',b'H10') : 'DBC.C8',
		(b'DBC',b'H11') : 'DBC.C8',
		(b'DBC',b'C4')  : 'DBC.C9',
		(b'DBC',b'H8')  : 'DBC.C9',
		(b'DBC',b'H9')  : 'DBC.C9',
		(b'DBC',b'C3')  : 'DBC.C10',
		(b'DBC',b'H6')  : 'DBC.C10',
		(b'DBC',b'H7')  : 'DBC.C10',
		(b'DBC',b'C2')  : 'DBC.C11',
		(b'DBC',b'H4')  : 'DBC.C11',
		(b'DBC',b'H5')  : 'DBC.C11',
		(b'DBC',b'C1')  : 'DBC.C12',
		(b'DBC',b'H1')  : 'DBC.C12',
		(b'DBC',b'H2')  : 'DBC.C12',
		(b'DBC',b'H3')  : 'DBC.C12',
		
		(b'DBS',b'O3')  : 'DBS.HG',
		(b'DBS',b'O2')  : 'DBS.HG',
		(b'DBS',b'O1')  : 'DBS.HG',
		(b'DBS',b'S1')  : 'DBS.HG',
		(b'DBS',b'C18') : 'DBS.AR',
		(b'DBS',b'C17') : 'DBS.AR',
		(b'DBS',b'H29') : 'DBS.AR',
		(b'DBS',b'C16') : 'DBS.AR',
		(b'DBS',b'H28') : 'DBS.AR',
		(b'DBS',b'C15') : 'DBS.AR',
		(b'DBS',b'H27') : 'DBS.AR',
		(b'DBS',b'C14') : 'DBS.AR',
		(b'DBS',b'H26') : 'DBS.AR',
		(b'DBS',b'C13') : 'DBS.AR',
		(b'DBS',b'C12') : 'DBS.C1',
		(b'DBS',b'H24') : 'DBS.C1',
		(b'DBS',b'H25') : 'DBS.C1',
		(b'DBS',b'C11') : 'DBS.C2',
		(b'DBS',b'H22') : 'DBS.C2',
		(b'DBS',b'H23') : 'DBS.C2',
		(b'DBS',b'C10') : 'DBS.C3',
		(b'DBS',b'H20') : 'DBS.C3',
		(b'DBS',b'H21') : 'DBS.C3',
		(b'DBS',b'C9')  : 'DBS.C4',
		(b'DBS',b'H19') : 'DBS.C4',
		(b'DBS',b'H18') : 'DBS.C4',
		(b'DBS',b'C8')  : 'DBS.C5',
		(b'DBS',b'H16') : 'DBS.C5',
		(b'DBS',b'H17') : 'DBS.C5',
		(b'DBS',b'C7')  : 'DBS.C6',
		(b'DBS',b'H14') : 'DBS.C6',
		(b'DBS',b'H15') : 'DBS.C6',
		(b'DBS',b'C6')  : 'DBS.C7',
		(b'DBS',b'H12') : 'DBS.C7',
		(b'DBS',b'H13') : 'DBS.C7',
		(b'DBS',b'C5')  : 'DBS.C8',
		(b'DBS',b'H10') : 'DBS.C8',
		(b'DBS',b'H11') : 'DBS.C8',
		(b'DBS',b'C4')  : 'DBS.C9',
		(b'DBS',b'H8')  : 'DBS.C9',
		(b'DBS',b'H9')  : 'DBS.C9',
		(b'DBS',b'C3')  : 'DBS.C10',
		(b'DBS',b'H6')  : 'DBS.C10',
		(b'DBS',b'H7')  : 'DBS.C10',
		(b'DBS',b'C2')  : 'DBS.C11',
		(b'DBS',b'H4')  : 'DBS.C11',
		(b'DBS',b'H5')  : 'DBS.C11',
		(b'DBS',b'C1')  : 'DBS.C12',
		(b'DBS',b'H1')  : 'DBS.C12',
		(b'DBS',b'H2')  : 'DBS.C12',
		(b'DBS',b'H3')  : 'DBS.C12'
			
	}


# Now create a column listing with reference to the density bins to be summed for a given molecule or fragment	
# This dictionary had key values of column title and entry a list of the columns in the density bin matrix
	column_list = {}

	for k,v in atom_lookup.items():
		if column_dict[k] in column_list:
			column_list[column_dict[k]].append(v[0])
		else:
			column_list[column_dict[k]]=[v[0],]


# Now create the array of numeric data: column 0 is the z-position data, columns 1... the molecule/fragment density

	save_data=np.empty([z_bins,len(column_list)+1],float)	
	for b in range(z_bins):
		save_data[b,0]=z_label(b)
		
# Sum the atom type densities into fragments and save in outpit array
	c_index=1
	for label,lst in column_list.items():
		save_data[:,c_index]= np.sum(dbins[:,lst],axis=1)
		c_index +=1
		
# Create Header string
	hdr=sim_data.simdimensions.units
	
	for label in column_list.keys():
		hdr +='\t'+label

		
	np.savetxt(filename, save_data, delimiter='\t', newline='\n', header=hdr)
		
#
#######
# Function to convert count dta to density data
def calculateChargeDensity(bins, atomlist, fcount):
	
# Data structure is float array same size as bins
	density = np.zeros(bins.shape,float)
	
#Calculate the multiplier for atom count to mole count per volume volume of a slab in the z direction converts to coulombs per m^3
	slab_multiplier	= (1.6021766e8)/(bin_width * sim_data.simdimensions.x * sim_data.simdimensions.y )
# Define the dictionary with partial charges based on molecule, atomtype pairs
#Use AMOEBA13 for HOH

	partialcharge = {
		(b'HOH',b'H1')  : 0.25983,
		(b'HOH',b'H2')  : 0.25983,
		(b'HOH',b'O')   : -0.51966,
		
		(b'CTC',b'C')   : 0.29898,
		(b'CTC',b'Cl1') : -0.07475,
		(b'CTC',b'Cl2') : -0.07475,
		(b'CTC',b'Cl3') : -0.07475,
		(b'CTC',b'Cl4') : -0.07475,
		
		(b'LI',b'Li')  : 1.0,
		(b'NA',b'Na')  : 1.0,
		(b'MG',b'Mg')  : 2.0,
		(b'K',b'K')  : 1.0,
		(b'CA',b'Ca')  : 2.0,
		(b'ZN',b'Zn')  : 2.0,
		
		(b'F',b'F')  : -1.0,
		(b'Cl',b'Cl')  : -1.0,
		(b'Br',b'Br')  : -1.0,
		(b'I',b'I')  : -1.0,
		
		
		(b'DDC',b'O2')  : -0.903,
		(b'DDC',b'O1')  : -0.903,
		(b'DDC',b'C12') : 1.050,
		(b'DDC',b'C11') : -0.36,
		(b'DDC',b'H22') : 0.058,
		(b'DDC',b'H23') : 0.058,
		(b'DDC',b'C10') : -0.127,
		(b'DDC',b'H20') : 0.0635,
		(b'DDC',b'H21') : 0.0635,
		(b'DDC',b'C9')  : -0.127,
		(b'DDC',b'H19') : 0.0635,
		(b'DDC',b'H18') : 0.0635,
		(b'DDC',b'C8')  : -0.127,
		(b'DDC',b'H16') : 0.0635,
		(b'DDC',b'H17') : 0.0635,
		(b'DDC',b'C7')  : -0.127,
		(b'DDC',b'H14') : 0.0635,
		(b'DDC',b'H15') : 0.0635,
		(b'DDC',b'C6')  : -0.127,
		(b'DDC',b'H12') : 0.0635,
		(b'DDC',b'H13') : 0.0635,
		(b'DDC',b'C5')  : -0.127,
		(b'DDC',b'H10') : 0.0635,
		(b'DDC',b'H11') : 0.0635,
		(b'DDC',b'C4')  : -0.127,
		(b'DDC',b'H8')  : 0.0635,
		(b'DDC',b'H9')  : 0.0635,
		(b'DDC',b'C3')  : -0.127,
		(b'DDC',b'H6')  : 0.0635,
		(b'DDC',b'H7')  : 0.0635,
		(b'DDC',b'C2')  : -0.127,
		(b'DDC',b'H4')  : 0.0635,
		(b'DDC',b'H5')  : 0.0635,
		(b'DDC',b'C1')  : -0.171,
		(b'DDC',b'H1')  : 0.057,
		(b'DDC',b'H2')  : 0.057,
		(b'DDC',b'H3')  : 0.057,
		
		
		(b'DDS',b'O3')  : -0.80,
		(b'DDS',b'O2')  : -0.80,
		(b'DDS',b'O1')  : -0.80,
		(b'DDS',b'S1')  : 1.306,
		(b'DDS',b'C12') : -0.10,
		(b'DDS',b'H25') : 0.097,
		(b'DDS',b'H24') : 0.097,
		(b'DDS',b'C11') : -0.127,
		(b'DDS',b'H22') : 0.0635,
		(b'DDS',b'H23') : 0.0635,
		(b'DDS',b'C10') : -0.127,
		(b'DDS',b'H20') : 0.0635,
		(b'DDS',b'H21') : 0.0635,
		(b'DDS',b'C9')  : -0.127,
		(b'DDS',b'H19') : 0.0635,
		(b'DDS',b'H18') : 0.0635,
		(b'DDS',b'C8')  : -0.127,
		(b'DDS',b'H16') : 0.0635,
		(b'DDS',b'H17') : 0.0635,
		(b'DDS',b'C7')  : -0.127,
		(b'DDS',b'H14') : 0.0635,
		(b'DDS',b'H15') : 0.0635,
		(b'DDS',b'C6')  : -0.127,
		(b'DDS',b'H12') : 0.0635,
		(b'DDS',b'H13') : 0.0635,
		(b'DDS',b'C5')  : -0.127,
		(b'DDS',b'H10') : 0.0635,
		(b'DDS',b'H11') : 0.0635,
		(b'DDS',b'C4')  : -0.127,
		(b'DDS',b'H8')  : 0.0635,
		(b'DDS',b'H9')  : 0.0635,
		(b'DDS',b'C3')  : -0.127,
		(b'DDS',b'H6')  : 0.0635,
		(b'DDS',b'H7')  : 0.0635,
		(b'DDS',b'C2')  : -0.127,
		(b'DDS',b'H4')  : 0.0635,
		(b'DDS',b'H5')  : 0.0635,
		(b'DDS',b'C1')  : -0.171,
		(b'DDS',b'H1')  : 0.057,
		(b'DDS',b'H2')  : 0.057,
		(b'DDS',b'H3')  : 0.057,
		

		(b'DBC',b'O2')  : -0.945,
		(b'DBC',b'O1')  : -0.945,
		(b'DBC',b'C19') : 1.25,
		(b'DBC',b'C18') : -0.288,
		(b'DBC',b'C17') : 0.0486,
		(b'DBC',b'H29') : -0.0028,
		(b'DBC',b'C16') : 0.0486,
		(b'DBC',b'H28') : -0.0028,
		(b'DBC',b'C15') : -0.0174,
		(b'DBC',b'H27') : -0.0482,
		(b'DBC',b'C14') : -0.0174,
		(b'DBC',b'H26') : -0.0482,
		(b'DBC',b'C13') : -0.052,
		(b'DBC',b'C12') : -0.087,
		(b'DBC',b'H24') : 0.0533,
		(b'DBC',b'H25') : 0.0533,
		(b'DBC',b'C11') : -0.127,
		(b'DBC',b'H22') : 0.0635,
		(b'DBC',b'H23') : 0.0635,
		(b'DBC',b'C10') : -0.127,
		(b'DBC',b'H20') : 0.0635,
		(b'DBC',b'H21') : 0.0635,
		(b'DBC',b'C9')  : -0.127,
		(b'DBC',b'H19') : 0.0635,
		(b'DBC',b'H18') : 0.0635,
		(b'DBC',b'C8')  : -0.127,
		(b'DBC',b'H16') : 0.0635,
		(b'DBC',b'H17') : 0.0635,
		(b'DBC',b'C7')  : -0.127,
		(b'DBC',b'H14') : 0.0635,
		(b'DBC',b'H15') : 0.0635,
		(b'DBC',b'C6')  : -0.127,
		(b'DBC',b'H12') : 0.0635,
		(b'DBC',b'H13') : 0.0635,
		(b'DBC',b'C5')  : -0.127,
		(b'DBC',b'H10') : 0.0635,
		(b'DBC',b'H11') : 0.0635,
		(b'DBC',b'C4')  : -0.127,
		(b'DBC',b'H8')  : 0.0635,
		(b'DBC',b'H9')  : 0.0635,
		(b'DBC',b'C3')  : -0.127,
		(b'DBC',b'H6')  : 0.0635,
		(b'DBC',b'H7')  : 0.0635,
		(b'DBC',b'C2')  : -0.127,
		(b'DBC',b'H4')  : 0.0635,
		(b'DBC',b'H5')  : 0.0635,
		(b'DBC',b'C1')  : -0.171,
		(b'DBC',b'H1')  : 0.057,
		(b'DBC',b'H2')  : 0.057,
		(b'DBC',b'H3')  : 0.057,
		
		
		(b'DBS',b'O3')  : -0.80164,
		(b'DBS',b'O2')  : -0.80164,
		(b'DBS',b'O1')  : -0.80164,
		(b'DBS',b'S1')  : 1.57285,
		(b'DBS',b'C18') : -0.17459,
		(b'DBS',b'C17') : -0.04022,
		(b'DBS',b'H29') : 0.02719,
		(b'DBS',b'C16') : -0.04022,
		(b'DBS',b'H28') : 0.02719,
		(b'DBS',b'C15') : 0.06060,
		(b'DBS',b'H27') : -0.0348,
		(b'DBS',b'C14') : 0.06060,
		(b'DBS',b'H26') : -0.0348,
		(b'DBS',b'C13') : -0.04091,
		(b'DBS',b'C12') : -0.09415,
		(b'DBS',b'H24') : 0.05817,
		(b'DBS',b'H25') : 0.05817,
		(b'DBS',b'C11') : -0.127,
		(b'DBS',b'H22') : 0.0635,
		(b'DBS',b'H23') : 0.0635,
		(b'DBS',b'C10') : -0.127,
		(b'DBS',b'H20') : 0.0635,
		(b'DBS',b'H21') : 0.0635,
		(b'DBS',b'C9')  : -0.127,
		(b'DBS',b'H19') : 0.0635,
		(b'DBS',b'H18') : 0.0635,
		(b'DBS',b'C8')  : -0.127,
		(b'DBS',b'H16') : 0.0635,
		(b'DBS',b'H17') : 0.0635,
		(b'DBS',b'C7')  : -0.127,
		(b'DBS',b'H14') : 0.0635,
		(b'DBS',b'H15') : 0.0635,
		(b'DBS',b'C6')  : -0.127,
		(b'DBS',b'H12') : 0.0635,
		(b'DBS',b'H13') : 0.0635,
		(b'DBS',b'C5')  : -0.127,
		(b'DBS',b'H10') : 0.0635,
		(b'DBS',b'H11') : 0.0635,
		(b'DBS',b'C4')  : -0.127,
		(b'DBS',b'H8')  : 0.0635,
		(b'DBS',b'H9')  : 0.0635,
		(b'DBS',b'C3')  : -0.127,
		(b'DBS',b'H6')  : 0.0635,
		(b'DBS',b'H7')  : 0.0635,
		(b'DBS',b'C2')  : -0.127,
		(b'DBS',b'H4')  : 0.0635,
		(b'DBS',b'H5')  : 0.0635,
		(b'DBS',b'C1')  : -0.171,
		(b'DBS',b'H1')  : 0.057,
		(b'DBS',b'H2')  : 0.057,
		(b'DBS',b'H3')  : 0.057
				
			}


# create the partial charge multiplier array
	pc_mult = np.empty(len(atomlist),float)
	for a in atomlist:
		index=atomlist[a][0]
		pc=partialcharge[a]
		pc_mult[index]=pc
	
			
# Calculate the atom type density bins

	density = (bins*pc_mult)*slab_multiplier/fcount

	return density
#######

#######
# Format and write chage density data
# Data output format:
#	first row is labels
#	First Column is z-depth
#	Columns are molecules for ions and solvents, and fragments of surfactants
def write_Charge_Density_Results(filename,cbins,atom_lookup,atomcount):
	
# Dictionary listing that will assign molecule,atom pairs to a column with label
	column_dict = {
		(b'HOH',b'H1') : 'H2O.H',
		(b'HOH',b'H2') : 'H2O.H',
		(b'HOH',b'O')  : 'H2O.O',
		
		(b'LI',b'Li')  : 'Li+',
		(b'NA',b'Na')  : 'Na+',
		(b'MG',b'Mg')  : 'Mg2+',
		(b'K',b'K')  : 'K+',
		(b'CA',b'Ca')  : 'Ca2+',
		(b'ZN',b'Zn')  : 'Zn2+',
		
		(b'F',b'F')  : 'F-',
		(b'Cl',b'Cl')  : 'Cl-',
		(b'Br',b'Br')  : 'Br-',
		(b'I',b'I')  : 'I-',
				
		(b'CTC',b'Cl1'): 'CCl4.Cl',
		(b'CTC',b'Cl2'): 'CCl4.Cl',
		(b'CTC',b'Cl3'): 'CCl4.Cl',
		(b'CTC',b'Cl4'): 'CCl4.Cl',
		(b'CTC',b'C')  : 'CCl4.C',
		
		(b'DDC',b'O2')  : 'DDC.HG',
		(b'DDC',b'O1')  : 'DDC.HG',
		(b'DDC',b'C12') : 'DDC.HG',
		(b'DDC',b'C11') : 'DDC.T',
		(b'DDC',b'H22') : 'DDC.T',
		(b'DDC',b'H23') : 'DDC.T',
		(b'DDC',b'C10') : 'DDC.T',
		(b'DDC',b'H20') : 'DDC.T',
		(b'DDC',b'H21') : 'DDC.T',
		(b'DDC',b'C9')  : 'DDC.T',
		(b'DDC',b'H19') : 'DDC.T',
		(b'DDC',b'H18') : 'DDC.T',
		(b'DDC',b'C8')  : 'DDC.T',
		(b'DDC',b'H16') : 'DDC.T',
		(b'DDC',b'H17') : 'DDC.T',
		(b'DDC',b'C7')  : 'DDC.T',
		(b'DDC',b'H14') : 'DDC.T',
		(b'DDC',b'H15') : 'DDC.T',
		(b'DDC',b'C6')  : 'DDC.T',
		(b'DDC',b'H12') : 'DDC.T',
		(b'DDC',b'H13') : 'DDC.T',
		(b'DDC',b'C5')  : 'DDC.T',
		(b'DDC',b'H10') : 'DDC.T',
		(b'DDC',b'H11') : 'DDC.T',
		(b'DDC',b'C4')  : 'DDC.T',
		(b'DDC',b'H8')  : 'DDC.T',
		(b'DDC',b'H9')  : 'DDC.T',
		(b'DDC',b'C3')  : 'DDC.T',
		(b'DDC',b'H6')  : 'DDC.T',
		(b'DDC',b'H7')  : 'DDC.T',
		(b'DDC',b'C2')  : 'DDC.T',
		(b'DDC',b'H4')  : 'DDC.T',
		(b'DDC',b'H5')  : 'DDC.T',
		(b'DDC',b'C1')  : 'DDC.T',
		(b'DDC',b'H1')  : 'DDC.T',
		(b'DDC',b'H2')  : 'DDC.T',
		(b'DDC',b'H3')  : 'DDC.T',
		
		(b'DDS',b'O3')  : 'DDS.HG',
		(b'DDS',b'O2')  : 'DDS.HG',
		(b'DDS',b'O1')  : 'DDS.HG',
		(b'DDS',b'S1')  : 'DDS.HG',
		(b'DDS',b'C12') : 'DDS.T',
		(b'DDS',b'H25') : 'DDS.T',
		(b'DDS',b'H24') : 'DDS.T',
		(b'DDS',b'C11') : 'DDS.T',
		(b'DDS',b'H22') : 'DDS.T',
		(b'DDS',b'H23') : 'DDS.T',
		(b'DDS',b'C10') : 'DDS.T',
		(b'DDS',b'H20') : 'DDS.T',
		(b'DDS',b'H21') : 'DDS.T',
		(b'DDS',b'C9')  : 'DDS.T',
		(b'DDS',b'H19') : 'DDS.T',
		(b'DDS',b'H18') : 'DDS.T',
		(b'DDS',b'C8')  : 'DDS.T',
		(b'DDS',b'H16') : 'DDS.T',
		(b'DDS',b'H17') : 'DDS.T',
		(b'DDS',b'C7')  : 'DDS.T',
		(b'DDS',b'H14') : 'DDS.T',
		(b'DDS',b'H15') : 'DDS.T',
		(b'DDS',b'C6')  : 'DDS.T',
		(b'DDS',b'H12') : 'DDS.T',
		(b'DDS',b'H13') : 'DDS.T',
		(b'DDS',b'C5')  : 'DDS.T',
		(b'DDS',b'H10') : 'DDS.T',
		(b'DDS',b'H11') : 'DDS.T',
		(b'DDS',b'C4')  : 'DDS.T',
		(b'DDS',b'H8')  : 'DDS.T',
		(b'DDS',b'H9')  : 'DDS.T',
		(b'DDS',b'C3')  : 'DDS.T',
		(b'DDS',b'H6')  : 'DDS.T',
		(b'DDS',b'H7')  : 'DDS.T',
		(b'DDS',b'C2')  : 'DDS.T',
		(b'DDS',b'H4')  : 'DDS.T',
		(b'DDS',b'H5')  : 'DDS.T',
		(b'DDS',b'C1')  : 'DDS.T',
		(b'DDS',b'H1')  : 'DDS.T',
		(b'DDS',b'H2')  : 'DDS.T',
		(b'DDS',b'H3')  : 'DDS.T',
		

		(b'DBC',b'O2')  : 'DBC.HG',
		(b'DBC',b'O1')  : 'DBC.HG',
		(b'DBC',b'C19') : 'DBC.HG',
		(b'DBC',b'C18') : 'DBC.AR',
		(b'DBC',b'C17') : 'DBC.AR',
		(b'DBC',b'H29') : 'DBC.AR',
		(b'DBC',b'C16') : 'DBC.AR',
		(b'DBC',b'H28') : 'DBC.AR',
		(b'DBC',b'C15') : 'DBC.AR',
		(b'DBC',b'H27') : 'DBC.AR',
		(b'DBC',b'C14') : 'DBC.AR',
		(b'DBC',b'H26') : 'DBC.AR',
		(b'DBC',b'C13') : 'DBC.AR',
		(b'DBC',b'C12') : 'DBC.T',
		(b'DBC',b'H24') : 'DBC.T',
		(b'DBC',b'H25') : 'DBC.T',
		(b'DBC',b'C11') : 'DBC.T',
		(b'DBC',b'H22') : 'DBC.T',
		(b'DBC',b'H23') : 'DBC.T',
		(b'DBC',b'C10') : 'DBC.T',
		(b'DBC',b'H20') : 'DBC.T',
		(b'DBC',b'H21') : 'DBC.T',
		(b'DBC',b'C9')  : 'DBC.T',
		(b'DBC',b'H19') : 'DBC.T',
		(b'DBC',b'H18') : 'DBC.T',
		(b'DBC',b'C8')  : 'DBC.T',
		(b'DBC',b'H16') : 'DBC.T',
		(b'DBC',b'H17') : 'DBC.T',
		(b'DBC',b'C7')  : 'DBC.T',
		(b'DBC',b'H14') : 'DBC.T',
		(b'DBC',b'H15') : 'DBC.T',
		(b'DBC',b'C6')  : 'DBC.T',
		(b'DBC',b'H12') : 'DBC.T',
		(b'DBC',b'H13') : 'DBC.T',
		(b'DBC',b'C5')  : 'DBC.T',
		(b'DBC',b'H10') : 'DBC.T',
		(b'DBC',b'H11') : 'DBC.T',
		(b'DBC',b'C4')  : 'DBC.T',
		(b'DBC',b'H8')  : 'DBC.T',
		(b'DBC',b'H9')  : 'DBC.T',
		(b'DBC',b'C3')  : 'DBC.T',
		(b'DBC',b'H6')  : 'DBC.T',
		(b'DBC',b'H7')  : 'DBC.T',
		(b'DBC',b'C2')  : 'DBC.T',
		(b'DBC',b'H4')  : 'DBC.T',
		(b'DBC',b'H5')  : 'DBC.T',
		(b'DBC',b'C1')  : 'DBC.T',
		(b'DBC',b'H1')  : 'DBC.T',
		(b'DBC',b'H2')  : 'DBC.T',
		(b'DBC',b'H3')  : 'DBC.T',

		(b'DBS',b'O3')  : 'DBS.HG',
		(b'DBS',b'O2')  : 'DBS.HG',
		(b'DBS',b'O1')  : 'DBS.HG',
		(b'DBS',b'S1')  : 'DBS.HG',
		(b'DBS',b'C18') : 'DBS.AR',
		(b'DBS',b'C17') : 'DBS.AR',
		(b'DBS',b'H29') : 'DBS.AR',
		(b'DBS',b'C16') : 'DBS.AR',
		(b'DBS',b'H28') : 'DBS.AR',
		(b'DBS',b'C15') : 'DBS.AR',
		(b'DBS',b'H27') : 'DBS.AR',
		(b'DBS',b'C14') : 'DBS.AR',
		(b'DBS',b'H26') : 'DBS.AR',
		(b'DBS',b'C13') : 'DBS.AR',
		(b'DBS',b'C12') : 'DBS.T',
		(b'DBS',b'H24') : 'DBS.T',
		(b'DBS',b'H25') : 'DBS.T',
		(b'DBS',b'C11') : 'DBS.T',
		(b'DBS',b'H22') : 'DBS.T',
		(b'DBS',b'H23') : 'DBS.T',
		(b'DBS',b'C10') : 'DBS.T',
		(b'DBS',b'H20') : 'DBS.T',
		(b'DBS',b'H21') : 'DBS.T',
		(b'DBS',b'C9')  : 'DBS.T',
		(b'DBS',b'H19') : 'DBS.T',
		(b'DBS',b'H18') : 'DBS.T',
		(b'DBS',b'C8')  : 'DBS.T',
		(b'DBS',b'H16') : 'DBS.T',
		(b'DBS',b'H17') : 'DBS.T',
		(b'DBS',b'C7')  : 'DBS.T',
		(b'DBS',b'H14') : 'DBS.T',
		(b'DBS',b'H15') : 'DBS.T',
		(b'DBS',b'C6')  : 'DBS.T',
		(b'DBS',b'H12') : 'DBS.T',
		(b'DBS',b'H13') : 'DBS.T',
		(b'DBS',b'C5')  : 'DBS.T',
		(b'DBS',b'H10') : 'DBS.T',
		(b'DBS',b'H11') : 'DBS.T',
		(b'DBS',b'C4')  : 'DBS.T',
		(b'DBS',b'H8')  : 'DBS.T',
		(b'DBS',b'H9')  : 'DBS.T',
		(b'DBS',b'C3')  : 'DBS.T',
		(b'DBS',b'H6')  : 'DBS.T',
		(b'DBS',b'H7')  : 'DBS.T',
		(b'DBS',b'C2')  : 'DBS.T',
		(b'DBS',b'H4')  : 'DBS.T',
		(b'DBS',b'H5')  : 'DBS.T',
		(b'DBS',b'C1')  : 'DBS.T',
		(b'DBS',b'H1')  : 'DBS.T',
		(b'DBS',b'H2')  : 'DBS.T',
		(b'DBS',b'H3')  : 'DBS.T'
		
	}


# Now create a column listing with reference to the density bins to be summed for a given molecule or fragment	
# This dictionary had key values of column title and entry a list of the columns in the density bin matrix
	column_list = {}

	for k,v in atom_lookup.items():
		if column_dict[k] in column_list:
			column_list[column_dict[k]].append(v[0])
		else:
			column_list[column_dict[k]]=[v[0],]


# Now create the array of numeric data: column 0 is the z-position data, columns 1... the molecule/fragment density

	save_data=np.empty([z_bins,len(column_list)+1],float)	
	for b in range(z_bins):
		save_data[b,0]=z_label(b)
		
# Sum the atom type densities into fragments and save in outpit array
	c_index=1
	for label,lst in column_list.items():
		save_data[:,c_index]= np.sum(cbins[:,lst],axis=1)
		c_index += 1
		
		
# Create Header string
	hdr=sim_data.simdimensions.units
	
	for label in column_list.keys():
		hdr +='\t'+label

		
	np.savetxt(filename, save_data, delimiter='\t', newline='\n', header=hdr)
		
#

if __name__ == '__main__':  
	mp.set_start_method('fork')
	main()						
	
