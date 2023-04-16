"""

FrameDataClass.py

	Purpose is to read and return simulaion trajectory data frame by frame.
	This is an class with the following characteistics:

		1. In init:
			Function input is file base and simulation #
			Open and parse the name.siminfo and name.simtopo files
				extract: # frames and real dimensions
				
			Parse the name.simtopo for atoms
				
			Open the binary frame file for reading
				
		2. function self.next()
			Reads one frame of data & increments count.  Print the frame number
		3. Variable self.numframes = number of frames
		4. Variable self.topo = numpy array of topology onfo
		5. Variable self.positions = numpy array of atom positions as x, y, z, positions (in nm)
		6. Variable self.forces = numpy array of atom forces (in kJ/mol/nm)
		
		



import sys

import os.path
# import numpy for extensive array and math functions
import numpy as np
#from openmmtools import integrators
# Garbage Collector!
import gc
"""

# String regular expressions
import re
# import numpy for extensive array and math functions
import numpy as np


class FrameData:
	def __init__(self, basefile, simnumber):

# First open and parse the .siminfo file for the number of frames

		with open(basefile+'.siminfo','r') as f:
			for line in f:
				if not line:
					print('End of file {} reached and no frame count found'.format(basefile+'.siminfo'))
					sys.exit(1)
				test = re.search('(binary)',line)
				if not (test == None):
					self.framecount = int(line[test.end()+2:])
					break
					
# Now open the .simtopo file and parse.  First size, then atoms

		with open(basefile+'.simtopo','r') as f:
			line=f.readline()
			contents=line.replace(',',' ').replace(':','').replace('Vec3(x=','').replace(')','').replace('y=','').replace('z=','').split()
			self.simdimensions=SimulationDimensions(contents[3],contents[4],contents[5],contents[6])
			for line in f:
				if not line:
					print('End of file {} reached and no frame count found'.format(basefile+'.siminfo'))
					sys.exit(1)
				test = re.search('Total Atom Count ',line)
				if not test == None:
					self.atomcount = int(line[test.end():])
					break	
	
# Create a numpy array with the correct dimensions
			self.atominfo=np.empty([self.atomcount],dtype=[('type','S4'),('atomicno','int'),('molecule','S4'),('moleculeindex','int')])

# Create empty dictionary		
			self.atomdict={}	
	
#Return to beginning					
			f.seek(0)
		
# Find the first line of atom data
			for line in f:
				if not re.search('atom',line)==None:
					break
	
# Read in atom number of lines with atom data
	
			for index in range(self.atomcount):
				line = f.readline()
				contents=line.split()
				atno=int(contents[3])
				self.atominfo[index]=(contents[1],atno,contents[6],int(contents[5]))
				
				if not (atno in self.atomdict):
					self.atomdict[atno]=float(contents[4])
					
# Open the binary trajectory file for reading

		self._binaryfile = open("{}_{:02d}_traj.framebin".format(basefile,simnumber),'rb')

# Create the frame array

		self._framedata=np.empty([self.atomcount,6],float)
		self._framesread = 0

# Init complete		

# Class function to close the binary file
	def done(self):
		self._binaryfile.close()
		
# Class function to load the next data set.

	def next(self):
		if self._framesread < self.framecount:
			self._framedata=np.reshape(np.fromfile(self._binaryfile,count=6*self.atomcount),[self.atomcount,6])
			self._framesread +=1
			
			self.positions = self._framedata[:,0:3]
			self.forces = self._framedata[:,3:6]
		else:
			print("Attempt to read past end of frame data\nForced exit")
			sys.exit(1)
					
				
# Class for simulation dimension
class SimulationDimensions:
	def __init__(self, x,y,z, units):
		self.x= float(x)
		self.y= float(y)
		self.z=float(z)
		self.units= str(units)

