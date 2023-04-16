#!/usr/bin/env python3
"""
TrajAverager.py

Purpose:	Reads interface analysis reqults for a single trajectory and averages. Write average file.
Types:		Mass Density, Charge Density, Force Profile

Interface analysis: For Mass density the water density is fit to a hyperbolic tangent function for each trajectory. fitting parameters and covarience are waved to a file for each tranjectory and for the average

input format: $python Interface_Traj_Average.py basefilename

Code will count the number of each type of file with the basefilename followed by file type.

"""
import sys
import numpy as np
import scipy 
import os
from scipy import optimize 

def main():
	
# Define the file type name extensions:
	global massdensity
	massdensity="_mass_density.txt"
	chargedensity="_charge_density.txt"
	forceprofile="_force_profile.txt"
	surfFFT="_interface_FFT.txt"
	surfFD="_interface_FD.txt"
	surfComp="_interface_comp.txt"
	
	
# This value will correct the depth to the interface tangent fit
	global zero_correct
	zero_correct = 0.0
	
	
# Check the input line
	if len(sys.argv)!=2:
		print("One input argument are expected: The input base file name for the simulation.")
		sys.exit(1)
		
	basefile = sys.argv[1]
	if basefile[len(basefile)-4:]==".pdb": basefile=basefile[0:len(basefile)-4]
	
	
# Count the  number of each type of 'density' file. Then average the densoity files and write out result
	masscount = 0
	while os.path.isfile("{}_{:02d}{}".format(basefile,masscount+1,massdensity)):
		masscount+=1
	
	if masscount>0:
		print("Averaging {} Mass Density Files".format(masscount))
		densityaverage(basefile,massdensity,masscount)

###	
	chargecount = 0
	while os.path.isfile("{}_{:02d}{}".format(basefile,chargecount+1,chargedensity)):
		chargecount+=1
		
	if chargecount>0:
		print("Averaging {:d} Charge Density Files".format(chargecount))
		densityaverage(basefile,chargedensity,chargecount)
		
###	
	forcecount = 0
	while os.path.isfile("{}_{:02d}{}".format(basefile,forcecount+1,forceprofile)):
		forcecount+=1
			
	
	if forcecount>0:
		print("Averaging {} Force profile Files".format(forcecount))
		densityaverage(basefile,forceprofile,forcecount)
		
###	
	FFTcount = 0
	while os.path.isfile("{}_{:02d}{}".format(basefile,FFTcount+1,surfFFT)):
		FFTcount+=1
		
		
	if FFTcount>0:
		print("Averaging {} FFT Interface Files".format(FFTcount))
		interfaceaverage(basefile,surfFFT,FFTcount)
		
###	
	FDcount = 0
	while os.path.isfile("{}_{:02d}{}".format(basefile,FDcount+1,surfFD)):
		FDcount+=1
		
		
	if FDcount>0:
		print("Averaging {} FD Interface Files".format(FDcount))
		interfaceaverage(basefile,surfFD,FDcount)

###	
	Compcount = 0
	while os.path.isfile("{}_{:02d}{}".format(basefile,Compcount+1,surfComp)):
		Compcount+=1
		
		
	if Compcount>0:
		print("Averaging {} Interface Composition Files".format(Compcount))
		compaverage(basefile,surfComp,Compcount)
		
###		
	
	return

#########

## Function to read and average data files.  Checks ti see if tis is a mass file

def densityaverage(basename, typename, filecount):

	global zero_correct

# First check if this is a force density
	InterfaceFit = typename == massdensity
	
	count=1
	
# Open first file then read and parse the header. Read the data into a 2-d array.
	
	
	with open("{}_{:02d}{}".format(basename,count,typename),"r") as density_file:
		header_line=density_file.readline()
		data_Array = np.loadtxt (density_file)
		
	header_items = header_line[2:].split()
	columns=len(header_items)
	
# Need the size to dimension the accumulated array.
	
	rows=data_Array.shape[0]
	
# Declare the array to save the file data
	
	full_array=np.zeros((filecount,rows,columns),dtype=float)
	
# Copy over the density data from the first file
	
	full_array[count-1]=data_Array
	
# Open, read, and copy remaining data files
	
	while os.path.isfile("{}_{:02d}{}".format(basename,count+1,typename)):
		count+=1
		with open("{}_{:02d}{}".format(basename,count,typename),"r") as density_file:
			header_line=density_file.readline()
			data_Array = np.loadtxt (density_file)
		full_array[count-1]=data_Array

# If this is a mass density average and water is in the mass listm then fit the mass density of water to a hyperbolic tangent.  The function returns a zero_correction
	
	
	if InterfaceFit and ('H2O' in header_items):
		zero_correct = interface_fit(basename,full_array,filecount,header_items)
		
	full_array[:,:,header_items.index('nm')]-=zero_correct
		
# Write out the average data 
		
	np.savetxt("{}_Average{}".format(basename,typename),np.average(full_array,axis=0),delimiter='\t', newline='\n', header=header_line[2:])
		
		
# writen out standard deviation data if more than one file analyzed		
	
	if filecount>1:
		np.savetxt("{}_StDev{}".format(basename,typename),np.std(full_array,axis=0),delimiter='\t', newline='\n', header=header_line[2:])
	
# If this is a mass density average then Calculate the water interface parameters
	

			
	return

### end density average function


## Interface fit complete the hyperbolic tangent fit to water and returns the depth zero correction as the fit inflection point
def interface_fit(basename,full_array,filecount,header_items):
	interface_Header = "\tmax\tinflection\twidth\n"
	maxHeight = 1.
	inflection = 0.
	width = 0.4
	guess=[maxHeight,inflection,width]

# Only fit data between depts in nm as analysis_start and Analysis_end.  Determine array indicies for these z values
	analysis_start=-1.0
	analysis_start_index=np.where(full_array[0][:,header_items.index ('nm'):] == analysis_start)[0][0]
	analysis_end = 1.0
	analysis_end_index=np.where(full_array[0][:,header_items.index ('nm'):] == analysis_end)[0][0]+1

# Extract the depth array

	Depth_Array=full_array[0,analysis_start_index:analysis_end_index,header_items.index ('nm')]


# Fit the data for each trajectory to the yyperbolic tangent function, and store data in arrays
	opt_values=np.empty([filecount,3],dtype=float)
	std_values=np.empty([filecount,3],dtype=float)

	for i in range(filecount):
		Interface_Array = full_array[i,analysis_start_index:analysis_end_index,header_items.index('H2O')]
		popt, cov = scipy.optimize.curve_fit (hyperbolicTanFunc, Depth_Array, Interface_Array,p0=guess)
		
		opt_values[i]=popt
		guess=popt
		std_values[i]= np.sqrt(np.diag (cov))
		
# Write to InterfaceFitting file		
		
	with open("{}_Water_Interface_Fit.txt".format(basename),"w") as fitfile:
		for i in range(filecount):
			fitfile.write("Water Interface Mass Density fit to hyperbolic tangent function for n={:d} data points\n".format(Depth_Array.shape[0]))
			fitfile.write("Trajectory {0:2d}\n".format(i+1))
			fitfile.write(interface_Header)
			fitfile.write("Fit:")
			for j in range(3): fitfile.write("\t{:.4f}".format(opt_values[i,j]))
			fitfile.write("\nStDev:")
			for j in range(3): fitfile.write("\t{:.4f}".format(std_values[i,j]))
			fitfile.write("\n\n")
			
			
		if filecount>1:
#First the aferage of the fits
			fitfile.write("\nStatistics for fit to {:d} trajectories:\n".format(filecount))
			fitfile.write("Average_Fit:")
			tmp=np.average(opt_values,axis=0)
			for j in range(3): fitfile.write("\t{:.4f}".format(tmp[j]))
			fitfile.write("\nStDev_Fit:")
			tmp=np.average(std_values,axis=0)
			for j in range(3): fitfile.write("\t{:.4f}".format(tmp[j]))
			fitfile.write("\n\n")

#Now the fit to the average
			Interface_Array = np.average(full_array,axis=0)[analysis_start_index:analysis_end_index,header_items.index('H2O')]
			popt, cov = scipy.optimize.curve_fit (hyperbolicTanFunc, Depth_Array, Interface_Array,p0=guess)
			fitfile.write("Fit parameters to Average Water Interface Density\n")
			fitfile.write("Fit:")
			for j in range(3): fitfile.write("\t{:.4f}".format(popt[j]))
			fitfile.write("\nStDev:")
			for j in range(3): fitfile.write("\t{:.4f}".format(np.sqrt(np.diag(cov))[j]))
			fitfile.write("\n\n")
			
# interface midpoint is either from the fit to the averate, or for a single trajectory

	return popt[1]

## End Interface_Fit function
	
## This is the interface fitting hyperbolic tangent function
def hyperbolicTanFunc (m, a, b, c):
	return (a/2)*(1-np.tanh ((m-b)/(c)))


# Function to average FFT and FD interface analysis results
def interfaceaverage(basefile,ftype,count):
	
# Read in the file data. Skip first line	
	for i in range(count):
		with open("{}_{:02d}{}".format(basefile,i+1,ftype),"r") as f:
			f.readline()
			xline=f.readline()
			yline=f.readline()
		
		
		if i==0:
			xlist=xline.split()
			labels=xlist[0]
			xlist=xlist[1:]
			numpoints=len(xlist)
			surfarray=np.empty([count,numpoints],dtype=float)
		
		surfarray[i,:]=np.asarray(yline.split()[1:],dtype=float)
		
# Now Construct the output array -- a transposition of the data such that data are in columns instead of rows
		
	if count>1:
		labels=labels+"\taverage\tStDev\n"
		resultarray=np.empty([numpoints,3],dtype=float)
		resultarray[:,2]=np.std(surfarray,axis=0)
	else:
		labels=labels+"\taverage\n"
		resultarray=np.empty([numpoints,2],dtype=float)
		
	resultarray[:,0]=np.asarray(xlist,dtype=float)
	resultarray[:,1]=np.average(surfarray,axis=0)

# Write the result file
	np.savetxt("{}_Average{}".format(basefile,ftype),resultarray,delimiter='\t', newline='\n', header=labels)
	
	
	return
# End Interface average

# Interface Composition Average
def compaverage(basefile,ftype,count):
	
	atomlist=[]
	for i in range(count):
		with open("{}_{:02d}{}".format(basefile,i+1,ftype),"r") as f:
			f.readline()
			headerline=f.readline()
			dataline=f.readline()
			
		if len(atomlist)==0:  # This is the first file read
			atomlist=headerline.split()
			columns=len(atomlist)
			comparray=np.empty([count,columns],dtype=float)
		
		comparray[i,:]=np.asarray(dataline.split()[1:],dtype=float)
	
# Now writeout data filewith header, average, stdef, and %
		
	compaverage=np.average(comparray,axis=0)
	sumaverage=np.sum(compaverage)
	compstdev=np.std(comparray,axis=0)
	comppercent=compaverage/sumaverage
	percentstd=compstdev/sumaverage
	
	with open("{}_Average{}".format(basefile,ftype),"w") as f:
		for i in range(columns): f.write("\t{}".format(atomlist[i]))
		f.write("\nAverage:")
		for v in compaverage: f.write("\t{:.1f}".format(v))
		f.write("\nStDev:")
		for v in compstdev: f.write("\t{:.1f}".format(v))
		f.write("\n\nratio:")
		for v in comppercent: f.write("\t{:.3%}".format(v))
		f.write("\nStDev:")
		for v in percentstd: f.write("\t{:.3%}".format(v))
		f.write("\n")

	
	return
# End Interface Compositon Averager

#########

if __name__ == '__main__':  
	main()						
	
	
