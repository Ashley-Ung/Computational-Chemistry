#!/usr/bin/env python3
# Ashley Ung: the interface_ionic_analysis code that i will be working on 

"""
InterfaceIonicAnalysis.py
Ashley Ung, Pacific University
Fall 2022 and Spring 2023 

Hydrophobic and Hydrophillic interface analysis for DDC 

This Program will read in all of the atom data from the simulations using FrameDataClass.py and will create a 2-d of the uppermost water atoms and ions (i.e. the interface) for each frame and save the z height and the atom type count. The resulting interface location data is analyzed with 2d FFT and fractal analysis.  Output three files:
		nn_##.interface_comp.txt including percenrt atomic composition of interface of O, H, ions, and holes
		nn_##.interface_FFT.txt averaged trace of 2-d fft with x axis and stdev
		nn_##.interface_FD.txt  averaged fractal area ratio, x axis and stdev
"""
import sys
import os.path
import re # String regular expressions
import numpy as np # import numpy for extensive array and math functions 
import multiprocessing as mp
import FrameDataClass as fd
import time
import math
from decimal import Decimal

#This is the size of your 2-d array for a surface z position and element type. Must be 2^n, and bigger is more often better
global size
size = 512

# Analysis is too slow to analyze every frame in a simulation.  This is the target number of frames to analyze.
global Analysis_Target
Analysis_Target = 100

# First check for the presence of expected data files		
def main():
	print("\nWelcome To Surface Analysis")
	# Check for presence of two parameters in input line
	if len(sys.argv)!=3:
		print("Two input arguments are expected: The input pdb file for the simulations, and the integer number of the simulation.")
		sys.exit(1)
	basefile = sys.argv[1]
	simnumber= int(sys.argv[2])
	# First check for presence of expected data files
	if not os.path.isfile(basefile+".siminfo"):
		print("Simulation info file {} does not exist.  Exiting!".format(basefile+".siminfo"))
		sys.exit(1)
	if not os.path.isfile(basefile+".simtopo"):
		print("Simulation info file {} does not exist.  Exiting!".format(basefile+".simtopo"))
		sys.exit(1)
		
	global filename
	filename="{}_{:02d}".format(basefile,simnumber)
	if not os.path.isfile(filename+"_traj.framebin"):
		print("Simulation info file "+filename+"_traj.framebin does not exist.  Exiting!")
		sys.exit(1)

	global sim_data
	#Read the data for the simulation
	sim_data=fd.FrameData(basefile,simnumber)
	print('Simulation dimensions {} x {} x {} {}'.format(sim_data.simdimensions.x, sim_data.simdimensions.y, sim_data.simdimensions.z, sim_data.simdimensions.units))
	if (sim_data.simdimensions.x != sim_data.simdimensions.y):
		print("Analysis assumes square dimensional simulation. Exiting!")
		sys.exit(1)

	global space_dim
	space_dim = sim_data.simdimensions.x

	print('Frame count ', sim_data.framecount)
	print('Atom count ', sim_data.atomcount)

	# Minimum z height for simulation
	global min_z
	min_z = -4.0 #in nanometers(-40 angstroms)

	# Variable to convert real position to array index and inverse.  These constants reduce calculations in time-critical core function
	global cell_dimension
	cell_dimension = space_dim/size
	global inverse_cell_dimension
	inverse_cell_dimension = size / space_dim
	
	# The number of atoms in the simulation
	global atomCount
	atomCount = sim_data.atomcount

	# The number of frames to Analyze
	global Analysis_Count
	Analysis_Count = min(Analysis_Target,sim_data.framecount)
	Frame_Skip = sim_data.framecount // Analysis_Count
	
	# The specific molecule we are looking for on our surface. The b signifies that this is a binary string.
	# Assumes water and ions w/ van der Waals Radius as specified.
	global interface_molecules
	global aqueousMolecules
	global hydrophobicMolecules
	aqueousMolecules = [b'HOH',b'LI',b'NA',b'K',b'MG',b'CA',b'ZN',b'F',b'Cl',b'Br',b'I']
	hydrophobicMolecules = [b'DDC', b'DBC', b'DDS', b'DBS', b'CTC']
	
	# A list of aqueous and hydrophobic atom types
	global aqueousTypes
	aqueousTypes = ["Holes", "H", "Li+", "O", "F-", "Na+", "Mg2+", "Cl-", "K+", "Ca2+", "Zn2+", "Br-", "I-"]
	global hydrophobicTypes
	#hydrophobicTypes = ["Holes", "Solvent", "DDC_Tail", "DDC_HeadGroup"]
	hydrophobicTypes = ["Holes", "Solvent", "DDC_Tail", "DDC_HeadGroup", "DDS_Tail", "DDS_HeadGroup"]
	#hydrophobicTypes = ["Holes", "Solvent", "DDC_Tail", "DDC_HeadGroup", "DBC_Tail", "DBC_HeadGroup", "DDS_Tail", "DDS_HeadGroup", "DBS_Tail", "DBS_HeadGroup"] # syntax error because other surfactants are not declared yet in the dictionary 
	
	# A dictionary of the hydrophobic and aqueous interface atom names and vdw radii values-- values are from Amoeba2018 data adjust to 0.7 of diameter
	adjust = 0.70 #Empirical Value 
	global interfaceData 
	interfaceData = {
		# Aqueous 
		(b"HOH", b"H1"):(aqueousTypes.index("H"), 0.2655*adjust),
		(b"HOH", b"H2"):(aqueousTypes.index("H"), 0.2655*adjust), 
		(b"LI", b"Li"):(aqueousTypes.index("Li+"), 0.22*adjust),
		(b"HOH", b"O"):(aqueousTypes.index("O"), 0.3405*adjust),
		(b"F", b"F"):(aqueousTypes.index("F-"), 0.343*adjust),
		(b"NA", b"Na"):(aqueousTypes.index("Na+"), 0.2955*adjust),
		(b"MG", b"Mg"):(aqueousTypes.index("Mg2+"), 0.29*adjust),
		(b"Cl", b"Cl"):(aqueousTypes.index("Cl-"), 0.412*adjust),
		(b"K", b"K"):(aqueousTypes.index("K+"), 0.368*adjust),
		(b"CA", b"Ca"):(aqueousTypes.index("Ca2+"), 0.359*adjust),
		(b"ZN", b"Zn"):(aqueousTypes.index("Zn2+"), 0.268*adjust),
		(b"BR", b"Br"):(aqueousTypes.index("Br-"), 0.432*adjust),
		(b"I", b"I"):(aqueousTypes.index("I-"), 0.461*adjust),
		
		# Tailgroups for Hydrophobic DDC
		(b"DDC", b"C1"):(hydrophobicTypes.index ("DDC_Tail"), 3.820e-01 * adjust),
		(b"DDC", b"H1"):(hydrophobicTypes.index ("DDC_Tail"), 2.980e-01 * adjust),
		(b"DDC", b"H2"):(hydrophobicTypes.index ("DDC_Tail"), 2.980e-01 * adjust),
		(b"DDC", b"H3"):(hydrophobicTypes.index ("DDC_Tail"), 2.980e-01 * adjust),
		(b"DDC", b"C2"):(hydrophobicTypes.index ("DDC_Tail"), 3.820e-01 * adjust),
		(b"DDC", b"H4"):(hydrophobicTypes.index ("DDC_Tail"), 2.980e-01 * adjust),
		(b"DDC", b"H5"):(hydrophobicTypes.index ("DDC_Tail"), 2.980e-01 * adjust),
		(b"DDC", b"C3"):(hydrophobicTypes.index ("DDC_Tail"), 3.820e-01 * adjust),
		(b"DDC", b"H6"):(hydrophobicTypes.index ("DDC_Tail"), 2.980e-01 * adjust),
		(b"DDC", b"H7"):(hydrophobicTypes.index ("DDC_Tail"), 2.980e-01 * adjust),
		(b"DDC", b"C4"):(hydrophobicTypes.index ("DDC_Tail"), 3.820e-01 * adjust),
		(b"DDC", b"H8"):(hydrophobicTypes.index ("DDC_Tail"), 2.980e-01 * adjust),
		(b"DDC", b"H9"):(hydrophobicTypes.index ("DDC_Tail"), 2.980e-01* adjust),
		(b"DDC", b"C5"):(hydrophobicTypes.index ("DDC_Tail"), 3.820e-01 * adjust),
		(b"DDC", b"H10"):(hydrophobicTypes.index ("DDC_Tail"), 2.980e-01 * adjust),
		(b"DDC", b"H11"):(hydrophobicTypes.index ("DDC_Tail"), 2.980e-01 * adjust),
		(b"DDC", b"C6"):(hydrophobicTypes.index ("DDC_Tail"), 3.820e-01 * adjust),
		(b"DDC", b"H12"):(hydrophobicTypes.index ("DDC_Tail"), 2.980e-01 * adjust),
		(b"DDC", b"H13"):(hydrophobicTypes.index ("DDC_Tail"), 2.980e-01 * adjust),
		(b"DDC", b"C7"):(hydrophobicTypes.index ("DDC_Tail"), 3.820e-01 * adjust),
		(b"DDC", b"H14"):(hydrophobicTypes.index ("DDC_Tail"), 2.980e-01 * adjust),
		(b"DDC", b"H15"):(hydrophobicTypes.index ("DDC_Tail"), 2.980e-01 * adjust),
		(b"DDC", b"C8"):(hydrophobicTypes.index ("DDC_Tail"), 3.820e-01 * adjust),
		(b"DDC", b"H16"):(hydrophobicTypes.index ("DDC_Tail"), 2.980e-01 * adjust),
		(b"DDC", b"H17"):(hydrophobicTypes.index ("DDC_Tail"), 2.980e-01 * adjust),
		(b"DDC", b"C9"):(hydrophobicTypes.index ("DDC_Tail"), 3.820e-01 * adjust),
		(b"DDC", b"H18"):(hydrophobicTypes.index ("DDC_Tail"), 2.980e-01 * adjust),
		(b"DDC", b"H19"):(hydrophobicTypes.index ("DDC_Tail"), 2.980e-01 * adjust),
		(b"DDC", b"C10"):(hydrophobicTypes.index ("DDC_Tail"), 3.820e-01 * adjust),
		(b"DDC", b"H20"):(hydrophobicTypes.index ("DDC_Tail"), 2.980e-01 * adjust),
		(b"DDC", b"H21"):(hydrophobicTypes.index ("DDC_Tail"), 2.980e-01 * adjust),
		(b"DDC", b"C11"):(hydrophobicTypes.index ("DDC_Tail"), 3.820e-01 * adjust),
		(b"DDC", b"H22"):(hydrophobicTypes.index ("DDC_Tail"), 2.980e-01 * adjust),
		(b"DDC", b"H23"):(hydrophobicTypes.index ("DDC_Tail"), 2.980e-01 * adjust),
		# Headgroups for Hydrophobic DDC
		(b"DDC", b"O1"):(hydrophobicTypes.index ("DDC_HeadGroup"), 0.355*adjust),
		(b"DDC", b"O2"):(hydrophobicTypes.index ("DDC_HeadGroup"), 0.355*adjust),
		(b"DDC", b"C12"):(hydrophobicTypes.index ("DDC_HeadGroup"), 0.382*adjust),
		
		# Tailgroups for Hydrophobic DDS 
		(b"DDS", b"C1"):(hydrophobicTypes.index ("DDS_Tail"), 3.820e-01 * adjust),
		(b"DDS", b"H1"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01 * adjust),
		(b"DDS", b"H2"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01 * adjust),
		(b"DDS", b"H3"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01 * adjust),
		(b"DDS", b"C2"):(hydrophobicTypes.index ("DDS_Tail"), 3.820e-01 * adjust),
		(b"DDS", b"H4"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01 * adjust),
		(b"DDS", b"H5"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01 * adjust),
		(b"DDS", b"C3"):(hydrophobicTypes.index ("DDS_Tail"), 3.820e-01 * adjust),
		(b"DDS", b"H6"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01 * adjust),
		(b"DDS", b"H7"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01 * adjust),
		(b"DDS", b"C4"):(hydrophobicTypes.index ("DDS_Tail"), 3.820e-01 * adjust),
		(b"DDS", b"H8"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01 * adjust),
		(b"DDS", b"H9"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01* adjust),
		(b"DDS", b"C5"):(hydrophobicTypes.index ("DDS_Tail"), 3.820e-01 * adjust),
		(b"DDS", b"H10"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01 * adjust),
		(b"DDS", b"H11"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01 * adjust),
		(b"DDS", b"C6"):(hydrophobicTypes.index ("DDS_Tail"), 3.820e-01 * adjust),
		(b"DDS", b"H12"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01 * adjust),
		(b"DDS", b"H13"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01 * adjust),
		(b"DDS", b"C7"):(hydrophobicTypes.index ("DDS_Tail"), 3.820e-01 * adjust),
		(b"DDS", b"H14"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01 * adjust),
		(b"DDS", b"H15"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01 * adjust),
		(b"DDS", b"C8"):(hydrophobicTypes.index ("DDS_Tail"), 3.820e-01 * adjust),
		(b"DDS", b"H16"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01 * adjust),
		(b"DDS", b"H17"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01 * adjust),
		(b"DDS", b"C9"):(hydrophobicTypes.index ("DDS_Tail"), 3.820e-01 * adjust),
		(b"DDS", b"H18"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01 * adjust),
		(b"DDS", b"H19"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01 * adjust),
		(b"DDS", b"C10"):(hydrophobicTypes.index ("DDS_Tail"), 3.820e-01 * adjust),
		(b"DDS", b"H20"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01 * adjust),
		(b"DDS", b"H21"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01 * adjust),
		(b"DDS", b"C11"):(hydrophobicTypes.index ("DDS_Tail"), 0.382 * adjust),   # dont know the vdw value for C12 currentlying using 3031
		(b"DDS", b"H22"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01 * adjust),
		(b"DDS", b"H23"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01 * adjust),
		(b"DDS", b"H24"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01 * adjust),					# is this correct vdw value ??? 
		(b"DDS", b"H25"):(hydrophobicTypes.index ("DDS_Tail"), 2.980e-01 * adjust),					# is this correct vdw value ???
		# Headgroups for Hydrophobic DDS
		(b"DDS", b"C12"):(hydrophobicTypes.index ("DDS_HeadGroup"), 2.980e-01 * adjust), # dont know the vdw value for C12 currentlying using 3032
		(b"DDS", b"S1"):(hydrophobicTypes.index ("DDS_HeadGroup"), 0.391 * adjust),
		(b"DDS", b"O1"):(hydrophobicTypes.index ("DDS_HeadGroup"), 0.351 * adjust),
		(b"DDS", b"O2"):(hydrophobicTypes.index ("DDS_HeadGroup"), 0.351 * adjust),
		(b"DDS", b"O3"):(hydrophobicTypes.index ("DDS_HeadGroup"), 0.351 * adjust),
		
		# Solvent 
		(b"CTC", b"C"):(hydrophobicTypes.index ("Solvent"), 0.36 * adjust),
		(b"CTC", b"Cl1"):(hydrophobicTypes.index ("Solvent"), 0.3898 * adjust),
		(b"CTC", b"Cl2"):(hydrophobicTypes.index ("Solvent"), 0.3898 * adjust),
		(b"CTC", b"Cl3"):(hydrophobicTypes.index ("Solvent"), 0.3898 * adjust),
		(b"CTC", b"Cl4"):(hydrophobicTypes.index ("Solvent"), 0.3898 * adjust)
	}
		
	
	# Looping through the data file and sorting the atoms by aqueous or hydrophobic by adding into designated array 
	global atom_type   # An array that stores the atom type 
	global atom_radius # An array that stores the van der waals radius 
	atom_type = mp.Array ('l',atomCount)
	atom_radius = mp.Array ('d', atomCount)
	global aqueousArray
	global hydrophobicArray 
	aqueousArray = []
	hydrophobicArray = []
	for atom in range (atomCount): 
		if sim_data.atominfo[atom]['molecule'] in aqueousMolecules:
			aqueousArray.append (atom)
		elif sim_data.atominfo[atom]['molecule'] in hydrophobicMolecules: 
			hydrophobicArray.append (atom)
		else: 
			sys.exit (1) # print atom number 
		atom_type[atom] = interfaceData[(sim_data.atominfo[atom]['molecule'], sim_data.atominfo[atom]['type'])][0]
		atom_radius[atom] = interfaceData[(sim_data.atominfo[atom]['molecule'], sim_data.atominfo[atom]['type'])][1]
		
	# Defines the range above and below the interface to test for surface atoms; in nm
		#global interfaceBound
		#interfaceBound = 2.0
	global aqueousInterfaceBound 
	aqueousInterfaceBound = 2.0 
	global hydrophobicInterfaceBound
	hydrophobicInterfaceBound = 3.0 
	
	# Aqueous Interface 2-d arrays for surface atom type and z height. Defined as 1-d array for purposes of multiprocessing, which only uses 1-d indexing	
	global z_type_aqueous										              
	mpArray_type_aqueous = mp.Array('l',size*size) #l = long integer data type
	z_type_aqueous = np.frombuffer(mpArray_type_aqueous.get_obj())
	global z_height_aqueous
	mpArray_height_aqueous = mp.Array('d',size*size) #d = decimal float data type 
	z_height_aqueous = np.frombuffer(mpArray_height_aqueous.get_obj(),dtype='float64')
	
	# Hydrophobic Interface 2-d arrays for surface atom type and z height. Defined as 1-d array for purposes of multiprocessing, which only uses 1-d indexing
	global z_type_hydrophobic 										              
	mpArray_type_hydrophobic = mp.Array('l',size*size) #l = long integer data type
	z_type_hydrophobic = np.frombuffer(mpArray_type_hydrophobic.get_obj())
	global z_height_hydrophobic
	mpArray_height_hydrophobic = mp.Array('d',size*size) #d = decimal float data type 
	z_height_hydrophobic = np.frombuffer(mpArray_height_hydrophobic.get_obj(),dtype='float64')
	
	# Roughness analysis for the aqueous interface FFT and Fractal
	# Array and constants for 2-d fft
	global interface_fft
	interface_fft = np.empty([Analysis_Count,size//2],float)
	
	global hanning2d
	han= np.hanning(size)
	hanning2d = np.sqrt(np.outer(han,han))
	
	global hanning_correct
	hanning_correct=12.0
	
	# Fractal Analysis constants
	# Box increment
	global box_increment
	box_increment = 2  # MUST BE EVEN
	global fd_size
	fd_size = size//box_increment
	
	# Fractal dimension data array 
	global interface_fd
	interface_fd = np.empty((Analysis_Count,fd_size))

	# Constants used for time display			
	global intervals
	intervals = (
		('weeks', 604800),  # 60 * 60 * 24 * 7
		('days', 86400),    # 60 * 60 * 24
		('hours', 3600),    # 60 * 60
		('minutes', 60),
		('seconds', 1),
		)

	# Interface composition (aqueous or hydrophobic) 
	global aqueousComposition 
	aqueousComposition = np.empty ([Analysis_Count, len (z_type_aqueous)],int) #atom_type? 
	global hydrophobicComposition 
	hydrophobicComposition = np.empty ([Analysis_Count, len (z_type_hydrophobic)], int) #atom_type? 
			
	#This portion fills the z_height and z_type arrays with the max z heights at every (x,y) and what atom type they are	
	print("Starting Interface Analysis.")

	start_time = time.time()
	Frames_Remaining = Analysis_Count
	
	# Frame loop	
	for frame in range(0,sim_data.framecount+1-Frame_Skip,Frame_Skip):
		print('Analysis of frame {} of {}'.format(frame + 1,sim_data.framecount))
		frame_start_time = time.time()

		# Load Frame Data --  skip over
		for i in range(Frame_Skip):
			sim_data.next()
		
		# Clear Result Arrays for Aqueous Interface 
		z_height_aqueous.fill(-1.0*aqueousInterfaceBound)
		z_type_aqueous.fill(0)
		
		# Clear Result Arrays for Hydrophobic Interface 
		z_height_hydrophobic.fill(+1.0*hydrophobicInterfaceBound)                     					#+1.0 for hydrophobic interface? 
		z_type_hydrophobic.fill(0)
		
		# Aqueous -- Set up the multiprocessing to use 1/2 of available cores.
		pool = mp.Pool(mp.cpu_count()//2)
		# The multiprocessing is called here and closed when all threads are finished.
		aqueousResults = pool.map (Test_Atom_Aqueous, [atom for atom in aqueousArray])
		pool.close () 
		pool.join ()
		
		# Hydrophobic -- Set up the multiprocessing to use 1/2 of available cores.
		pool = mp.Pool(mp.cpu_count()//2)
		# The multiprocessing is called here and closed when all threads are finished.
		hydrophobicResults = pool.map (Test_Atom_Hydrophobic, [atom for atom in hydrophobicArray])
		pool.close ()
		pool.join ()

		# Record surface composition
		print ("\nDetermining Aqueous Interface Composition and Location")
		append_aqueous_composition(frame//Frame_Skip) # returns the floor value for both integer and floating point arguments after division
		print ("\nDetermining Hydrophobic Interface Composition and Location")
		append_hydrophobic_composition(frame//Frame_Skip)
		
		# FFT and Fractal Dimension analysis for Aqueous Interface 
		print ("\nBegin FFT analysis for Aqueous Interface")
		Analysis_FFT (frame//Frame_Skip)
		print ("\nBegin Fractal Dimension Analysis for Aqueous Interface")
		Analysis_FD (frame//Frame_Skip)
		
		# This portion estimates the remaining time left in the data retrieval
		Frames_Remaining -= 1
		remainingSeconds = (time.time() - frame_start_time) * (Frames_Remaining)
		
		print ("\nFrame analysis Complete\n\tEstimated Remaining Time = {} ".format(display_time(remainingSeconds)))	
	
	sim_data.done ()
	print ("Writing Data Files")
	write_aqueous_composition()
	write_hydrophobic_composition()
	write_FFT ()
	write_FD ()
	
	print("\nExecution time = {} ".format(display_time(time.time() - start_time)))
	
	print("Surface Analysis Complete for {}_{:02d}".format(basefile,simnumber))	
	
	return
### End Main

#######	
# This function assumes atom is index of atom in water (including ions), and tests if the sphere of the atom is above the current interface. 
# Does not need to be called sequentially, and therefore is the core multiprocess function.
def Test_Atom_Aqueous(atom):
	global atomLookup
	global sim_data
	# First test if atom is within +/- interfacebound (8 angstrom defined above) of interface 0	
	atomType = int(atom_type[atom])
	atomRadius = atom_radius[atomType]
	z_center = sim_data.positions[atom,2]
	if z_center - atomRadius <= aqueousInterfaceBound and z_center + atomRadius >= -aqueousInterfaceBound:	
		x_pos = sim_data.positions[atom,0]
		y_pos = sim_data.positions[atom,1]
		# Calculate the surface array index bounds for the atom radius, the ceil function returns the smallest intege >= to x 
		x_minIndex = math.ceil((x_pos - atomRadius)*inverse_cell_dimension)
		x_maxIndex = math.ceil((x_pos + atomRadius)*inverse_cell_dimension)
		y_minIndex = math.ceil((y_pos - atomRadius)*inverse_cell_dimension)
		y_maxIndex = math.ceil((y_pos + atomRadius)*inverse_cell_dimension)
		for x_index in range(x_minIndex,x_maxIndex):
			for y_index in range(y_minIndex,y_maxIndex):
				#This calculates the intersection of the atom sphere with the line in the z direction at the surface index point. Uses always the larger intersection point
				zsquared = atomRadius**2 - (x_index*cell_dimension - x_pos)**2 - (y_index*cell_dimension - y_pos)**2
				if zsquared >= 0.0:
					z = z_center + np.sqrt(zsquared)
					# Test if the calculated height of atom at array position (x,y) is greater that current. If yes, then replace then height ant type.
					# Note: the arrays are 1-d for multiprocessing, therefore the index must be calculated.  Also, to accout for periodic boundaries
					# calculate the index with a %size -- i.e. the remainder.
					matrix_index=x_index%size + ((y_index%size)*size)
					if z > z_height_aqueous[matrix_index]:
						z_height_aqueous[matrix_index] = z
						z_type_aqueous[matrix_index] = atomType
	return True

def Test_Atom_Hydrophobic(atom):
	global atomLookup
	global sim_data
	# First test if atom is within +/- interfacebound (8 angstrom defined above) of interface 0	
	atomType = int(atom_type[atom])
	atomRadius = atom_radius[atomType]
	z_center = sim_data.positions[atom,2]
	if z_center - atomRadius <= hydrophobicInterfaceBound and z_center + atomRadius >= -hydrophobicInterfaceBound:	
		x_pos = sim_data.positions[atom,0]
		y_pos = sim_data.positions[atom,1]
		# Calculate the surface array index bounds for the atom radius	
		x_minIndex = math.ceil((x_pos - atomRadius)*inverse_cell_dimension)
		x_maxIndex = math.ceil((x_pos + atomRadius)*inverse_cell_dimension)
		y_minIndex = math.ceil((y_pos - atomRadius)*inverse_cell_dimension)
		y_maxIndex = math.ceil((y_pos + atomRadius)*inverse_cell_dimension)
		for x_index in range(x_minIndex,x_maxIndex):
			for y_index in range(y_minIndex,y_maxIndex):
				#This calculates the intersection of the atom sphere with the line in the z direction at the surface index point. Uses always the larger intersection point
				zsquared = atomRadius**2 - (x_index*cell_dimension - x_pos)**2 - (y_index*cell_dimension - y_pos)**2
				if zsquared >= 0.0:
					z = z_center - np.sqrt(zsquared)
					# Test if the calculated height of atom at array position (x,y) is greater that current. If yes, then replace then height ant type.
					# Note: the arrays are 1-d for multiprocessing, therefore the index must be calculated.  Also, to accout for periodic boundaries
					# calculate the index with a %size -- i.e. the remainder.
					matrix_index=x_index%size + ((y_index%size)*size)
					if z < z_height_hydrophobic[matrix_index]:
						z_height_hydrophobic[matrix_index] = z
						z_type_hydrophobic[matrix_index] = atomType
	return True

# End of the multiprocess core function

#######
# appends aqueous composition count to array
def append_aqueous_composition(frame):
	# Iterate through the atom types: Print frame result and save to interface composition 2-D array
	i = 0
	for atom in range ((len(aqueousTypes))):
		n= (z_type_aqueous == atom).sum()
		if not n == 0:			 									#check if atom is found in the file 
			print("{} count = {:d}".format(aqueousTypes[atom],n)) 
		aqueousComposition[frame,i]=n
		i+=1
	return True

#######
# appends hydrophobic composition count to array
def append_hydrophobic_composition(frame):
	# Iterate through the atom types: Print frame result and save to interface composition 2-D array
	i = 0
	#atomType = int(atom_type[atom])
	for atom in range ((len(hydrophobicTypes))):
		n= (z_type_hydrophobic == atom).sum()
		if not n == 0: 
			print("{} count = {:d}".format(hydrophobicTypes[atom],n))  #check if atom is found in the file
		hydrophobicComposition[frame,i]=n
		i+=1
	return True
	
#####	
# Completes 2-D fft, extracts diagonal and appends to array	
# calculates a power spectrum with a Hanning nwindow to reduce artifacts, then does approximate amplitude correction & normalization.
#Data saved for each frame in array.  File write  average, std calculation
def Analysis_FFT(frame):
	z_square = np.reshape(z_height_aqueous,(size,size))
	ps = np.abs(np.fft.rfft2(z_square * hanning2d))*hanning_correct/size**2
	for i in range(size//2):
		interface_fft[frame, i] = ps[i+1,i+1]
	return True
	
#####
# Do the fractal dimension analysis using pyramid analysis	This is a revision of the Clark method:  The eight panel method of W. Sun
def Analysis_FD(frame):
	# Some variables precalculated for efficiency	
	global x0, xm, x1, l2, l4, box_len
	#loop linearly through multiple of starting box length
	for box_len in range(box_increment, size +1, box_increment):
		Sum_FD_Area=0.0   # Collects calculated area in a row -- from MP function 
		count=len(range (0, size, box_len))  # need a counter to know the number of cells in a column -- for normalization
		
		# The base of the real pyramid is invarient for a given box length.  Center of pyramic is at x=0,y=0. Corners +-box_len/2
		corner = box_len//2*cell_dimension
		l2=corner*corner
		l4=l2*l2

		# Loop through columns spaced by the length of box  0 and 1 are coreners  m is midpoint  data wrapping requires %size
		for x0 in range (0, size, box_len):
			xm=(x0+box_len//2)%size
			x1=(x0+box_len)%size
			for y0 in range(0,size,box_len):
				Sum_FD_Area +=Sum_Cell_Aqueous(y0)
			
		#Save the sum and normalize by the area of the 'flat' base beneath it.
		interface_fd[frame,(box_len//box_increment)-1]= Sum_FD_Area/((box_len*cell_dimension*count)**2)
	return True

# Sum cell for aqueous interface 
def Sum_Cell_Aqueous (ystart):
	y0=ystart*size
	ym=((ystart+box_len//2)%size)*size
	y1=((ystart+box_len)%size)*size
	# calculate the coordinates of the 4 corners and four edge points going around the box of the current analysis box	
	# Note: center point is 0,0,0 and all z heights adjusted to this value.  indices precalculated for efficiency
	# This method uses the algebraic calculation of the 8 prism area, simplified and made efficient.  math.sqrt is 5x faster than np.sqrt!  
	# Sage was used to calculate the return expression		
	center= z_height_aqueous[xm + ym]
	z00= z_height_aqueous[x0 + y0]-center
	z01= z_height_aqueous[x0 + ym]-center
	z02= z_height_aqueous[x0 + y1]-center
	z12= z_height_aqueous[xm + y1]-center
	z22= z_height_aqueous[x1 + y1]-center
	z21= z_height_aqueous[x1 + ym]-center
	z20= z_height_aqueous[x1 + y0]-center
	z10= z_height_aqueous[xm + y0]-center
	return 0.5*(math.sqrt(l4 + l2*(2*(z01*z01 - z00*z01) + z00*z00)) + math.sqrt(l4 + l2*(2*(z01*z01 - z01*z02) + z02*z02)) + math.sqrt(l4 + l2*(2*(z10*z10 - z00*z10) + z00*z00)) + math.sqrt(l4 + l2*(2*(z10*z10 - z10*z20) + z20*z20)) + math.sqrt(l4 + l2*(2*(z12*z12 - z02*z12) + z02*z02)) + math.sqrt(l4 + l2*(2*(z12*z12 - z12*z22) + z22*z22)) + math.sqrt(l4 + l2*(2*(z21*z21 - z20*z21) + z20*z20)) + math.sqrt(l4 + l2*(2*(z21*z21 - z21*z22) + z22*z22)))

####### Output functions

# Aqueous Interface Composition, Average & Standard Deviation
def write_aqueous_composition():
	f=open(filename+"_aqueous_interface_comp.txt",'w')
	
	a_comp_ave = np.average (aqueousComposition, axis = 0)
	a_comp_std = np.std (aqueousComposition, axis = 0)
	
	f.write ('Aqueous Interface Composition\n')
	for atom in range ((len(aqueousTypes))):
		if not a_comp_std[atom] == 0: 
			f.write("\t{} ".format(aqueousTypes[atom])) 
		
	f.write ('\nAqueous Interface Average')
	for atom in range ((len(aqueousTypes))):
		if not a_comp_ave[atom] == 0:
			f.write ('\t{:.1f}'.format(a_comp_ave[atom]))
		
	f.write ('\nAqueous Interface Standard Deviation')
	for atom in range ((len(aqueousTypes))):
		if not a_comp_std[atom] == 0:
			f.write ('\t{:.1f}'.format(a_comp_std[atom]))
		
	f.write ('\nSampleCount = {:d}\n'.format(Analysis_Count))
	f.close ()
		
	return True

# Hydrophobic Interface Composition, Average & Standard Deviation 
def write_hydrophobic_composition():
	f=open(filename+"_hydrophobic_interface_comp.txt",'w')
	
	h_comp_ave = np.average (hydrophobicComposition, axis = 0)
	h_comp_std = np.std (hydrophobicComposition, axis = 0)
	
	f.write ('Hydrophobic Interface Composition\n')
	for atom in range ((len(hydrophobicTypes))):
		if not h_comp_std[atom] == 0:
			f.write ("\t{} ".format(hydrophobicTypes[atom]))
		
	f.write ('\nHydrophobic Interface Average')
	for atom in range(len(hydrophobicTypes)):
		if not h_comp_ave[atom] == 0.0: 
			f.write ('\t{:.1f}'.format(h_comp_ave[atom]))
		
	f.write ('\nHydrophobic Interface Standard Deviation')
	for atom in range (len(hydrophobicTypes)):
		if not h_comp_std[atom] == 0.0:
			f.write ('\t{:.1f}'.format(h_comp_std[atom]))
		
	f.write ('\nSampleCount = {:d}\n'.format(Analysis_Count))
	f.close ()
	return True

# Aqueous Interface FFT
def write_FFT():
	f=open(filename+"_aqueous_interface_FFT.txt",'w')
	
	FFT_ave = np.average(interface_fft, axis = 0)
	FFT_std = np.std(interface_fft, axis = 0)
	
	f.write('Aqueous Interface FFT\n')
	
	f.write('x(nm^-1)')
	for i in range(size//2-1):
		f.write('\t{:.3e}'.format(space_dim/(i+1)))
		
	f.write('\naverage')
	for i in range(size//2-1):
		f.write('\t{:.3e}'.format(FFT_ave[i]))
		
	f.write('\nstdev')
	for i in range(size//2-1):
		f.write('\t{:.3e}'.format(FFT_std[i]))
		
	f.write('\nSampleCount = {:d}\n'.format(Analysis_Count))
	f.close()
	return True

# Aqueous Interface FD 
def write_FD():
	f=open(filename+"_aqueous_interface_FD.txt",'w')
	
	FD_x_values = np.arange(box_increment,size+1,box_increment)* cell_dimension
	FD_ave = np.average(interface_fd, axis = 0)
	FD_std = np.std(interface_fd, axis = 0)
	
	f.write('Interface Fractal Dimension')
	
	f.write('\nx(nm)')
	for i in range(size//box_increment):
		f.write('\t{:.3e}'.format(FD_x_values[i]))
		
	f.write('\naverage')
	for i in range(size//box_increment):
		f.write('\t{:.3e}'.format(FD_ave[i]))
		
	f.write('\nstdev')
	for i in range(size//box_increment):
		f.write('\t{:.3e}'.format(FD_std[i]))
		
	f.write('\nSampleCount = {:d}\n'.format(Analysis_Count))
	f.close()
	return True

########
def display_time(seconds, granularity=2):
	result = []
	for name, count in intervals:
		value = seconds // count
		if value:
			seconds -= value * count
			if value == 1:
				name = name.rstrip('s')
			result.append("{} {}".format(value, name))
	return ', '.join(result[:granularity])
	
#########
if __name__ == '__main__':  
	mp.set_start_method ('fork')
	main ()
	
	
