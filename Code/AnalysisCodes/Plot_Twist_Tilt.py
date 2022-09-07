"""
Plot_Twist_Tilt.py
Kevin Johnson Pacific University
Summer 2022

This Program will read in output from Twist_Tilt analysis and generate semicircular contournplots
"""
import sys

import os.path
# String regular expressions

import numpy as np

import math

from matplotlib import pyplot as plt

import matplotlib.tri as tri

import matplotlib.colors as colors

		
def main():
	print("Welcome To Twist and Tilt Plotter")
	
#This is the size the analysis bins

	depthBins = 12
	
	fig_rows=3
	fig_cols=4
	
	tiltBins = 40
	tiltValues=np.linspace(1,-1,tiltBins)
	
	twistBins = 20
	twistValues=np.linspace(0,1,twistBins)
	
	X_twist = np.empty(tiltBins*twistBins)
	Y_tilt = np.empty(tiltBins*twistBins)
	Z = np.empty([depthBins,tiltBins*twistBins])
	labels=[]
	
# Fill linear X abnd Y arrays with the plot values.  Note twist values generate a semicircle.	
	for y_index in range(tiltBins):
		for x_index in range(twistBins):
			X_twist[y_index*twistBins+ x_index]=twistValues[x_index]*math.sqrt(1-tiltValues[y_index]**2) 
			Y_tilt[y_index*twistBins+ x_index]=tiltValues[y_index]

# This triangulates the x-y values for plotting
	traing=tri.Triangulation(X_twist,Y_tilt)

# Check for presence of two parameters in input line
	global filename
	filename=''
	for i in range(len(sys.argv)-1):
		filename=filename+sys.argv[i+1]+' '
	filename=filename[0:len(filename)-1]
	if filename[len(filename)-4:]==".pdb" or filename[len(filename)-4:]==".chk": filename=filename[0:len(filename)-4]
	if filename[len(filename)-14:]!="_TwistTilt.txt" : filename=filename+"_TwistTilt.txt"
	
#
	if not os.path.isfile(filename):
		print("Simulation info file "+filename+" does not exist.  Exiting!")
		sys.exit(1)

# Open file and read first informational line
	f=open(filename,'r')
	line=f.readline()
	plotTitle=line[1+line.rfind("/"):len(line)]
	dataMax=float(f.readline().rstrip('\n'))

# Read label and each data array
	for count in range(depthBins):
		labels.append(f.readline())	

#read data line by line
		for y_index in range(tiltBins):
			rowvalues=f.readline().split()
			for x_index in range(twistBins):
# Assign Z values  note the normalization constant (from Holte paper) this should acount for the reduced probability space at tilt or twist approaching ±1

				Z[count,y_index*twistBins+ x_index]=float(rowvalues[x_index]) *math.sqrt(1-X_twist[x_index]**2-Y_tilt[y_index]**2)

	f.close()
					
# initiate a figure with a grid of plots
					
	fig, axs = plt.subplots(nrows=fig_rows, ncols=fig_cols,sharex=True,sharey=True,figsize=[4.8,6.4])
	fig.suptitle(plotTitle, fontsize=12)
	cmap = plt.get_cmap('gist_heat').reversed()
#	cmap = plt.get_cmap('turbo')
	norm =colors.SymLogNorm(linthresh=0.01, vmin=np.amin(Z), vmax=np.amax(Z))
					
					
	for plotrow in range(fig_rows):
		for plotcol in range(fig_cols):
# Generate each of the subplots
			axs[plotrow,plotcol].set_title(labels[(fig_rows -1 -plotrow)*fig_cols+plotcol], fontsize=7)
			axs[plotrow,plotcol].set_aspect('equal')
			if plotrow == fig_rows-1:
				axs[plotrow,plotcol].set_xlabel('Twist')
			if plotcol == 0 :
				axs[plotrow,plotcol].set_ylabel('Tilt')
			
					# Next construct the plot
			axs[plotrow,plotcol].tripcolor(traing,Z[(fig_rows -1 -plotrow)*fig_cols+plotcol,:], shading='gouraud', cmap=cmap, norm=norm)
			
	
	fig.tight_layout()
	plt.subplots_adjust(top=0.875,left=0.1,wspace=0.1)

# Must chose one or the other: display or save
	plt.show()
#	plt.savefig(filename[:-4]+'.png')
	

	return
## end Main



if __name__ == '__main__':  
		main()
	