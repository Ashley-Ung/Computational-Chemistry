#!/usr/bin/python
"""

imulationTester.py
	
K.E. Johnson, Pacific University
Summer 2022

	Purpose is to run simulation on a given surfactant system.
	Thermal processing is:
		1. Minimize
		2. heat to  temp
		3. sumulate
		
	Reporters:
		DCD for animation and VMD analysis
		Progress / energy for meta analysis
		Frame data -- custom reporter that saves z- adjusted position and atom force data in binary format
		
		
"""


import sys
import traceback
import os.path
import datetime
import time
# import numpy for extensive array and math functions
import numpy as np
# openMM Imports
from openmm.app import *
from openmm import *
from openmm.unit import *
#from openmmtools import integrators
# Garbage Collector!
import gc

def main():
#main will read the input file base then run simulation as planned


# These two globalvariables are used to count simulation errors that should be trapped by 'try:'
# Used for each time simulation steps are run

# Check for presence of the pdb file
	if len(sys.argv) != 2:
		print("One input arguments is expected: The input pdb file for the test simulation.")
		sys.exit(1)
		
	basefile = sys.argv[1]
	
	if basefile[len(basefile)-4:]==".pdb": basefile=basefile[0:len(basefile)-4]
		
	if not os.path.isfile(basefile+".pdb"):
		print("File {} does not exist.  Exiting!".format(basefile+".pdb"))
		sys.exit(1)
		
	
##################
#  Simulation Constants

# Simulation Step Size
	sim_step= 2.0 * femtoseconds

# temperature setpoints
	T_Simulation=300
	
	
# Trajectory Length & Repeat

	Trajectory_time = 0.02 * nanoseconds  
	Trajectory_steps = round (Trajectory_time / sim_step)
		
# Reporting 	
	AnimationPoints = 200  
	ReportingPoints = 50  

	TrajReporterInterval=round(Trajectory_steps/ReportingPoints)
	TrajAnimationInterval=round(Trajectory_steps/AnimationPoints)
	
	


# Initialize the simulation 

	mySimulation=InitSimulation(basefile,T_Simulation,sim_step)

# minimize forces

	print('     Minimizing System Energy' )
			
	ExecuteSimulation(mySimulation, 0)
					
	print ('     Minimization Complete')
					

# Set Temp 
	print("Setting Temp at {}K".format(T_Simulation))
	SetSimTemp(mySimulation,T_Simulation,50)
		
	
	print("Starting data collection")
	CollectTrajectory(basefile+"_test", mySimulation, T_Simulation, Trajectory_steps, TrajReporterInterval, TrajAnimationInterval)
		
			
	print("Simulation complete.  Congratulations!")


	return(True)
	
#  END main



#########################################################	
def InitSimulation(FileBase,T_Start,step_size):
#Initialize and minimize simulation


	print("Initializing simulation.")	

	Sim_Temp=T_Start * kelvin

	
# Read in pdb file
	pdb_in=FileBase+".pdb"

	print ("    Loading initial pdb file")
	pdb = PDBFile(pdb_in)
	pdb.topology.setUnitCellDimensions((4.0, 4.0, 12.)*nanometers)
	
# Set up for system, NoCutof for Amoeba, Constrain only x-H bond lengths
	print ('    Loading Force Field')
	forcefield = ForceField('../AMOEBAff/solvents.xml','../AMOEBAff/surfactantsTest.xml','../AMOEBAff/ions.xml' )
	
	print ('    Setting up system')

# can consider 'direct' polarization calculation method
	system = forcefield.createSystem(	pdb.topology, 
										nonbondedMethod = PME,
										nonbondedCutoff=1.0*nanometer,
										vdwCutoff=1.2*nanometer,
										polarization='extrapolated', 
										constraints = HBonds)
										

# Set up integrator
	print ('    Setting up integrator')
	integrator_Friction = 1.0/picoseconds	# friction is 1.0 per ps, Future possibliltuy more frequent collisions
	
#Langevin Integrator. May consider experimental integrator in next iteration
#	integrator = LangevinIntegrator(Sim_Temp, integrator_Friction, step_size)
	integrator = LangevinMiddleIntegrator(Sim_Temp, integrator_Friction, step_size)

# setup to run CUDA at mixed precision. Single precision for forces, integrate @ double precision.
	print ('    Setting Platform to CUDA')
	platform = Platform.getPlatformByName('CUDA')
	properties = {'DeviceIndex': '0', 'CudaPrecision': 'mixed'}

# simulation is the object that accomplishes the MD work
	print ('Creating Simulation')
	mySimulation = Simulation(pdb.topology, 
							system, 
							integrator, 
							platform, 
							properties)

	mySimulation.context.setPositions(pdb.positions)
	
	return mySimulation	
	
# End InitSimulation






#########################################################	
def SetSimTemp(mySimulation,Temp,Dwell):
# Bring simulation to desired Temperature

	print("Bringing Simulation to set Temperature: {} K".format(Temp))
	
	Sim_Temp = Temp*kelvin
	
# Calculate Temperature Conversion. Used to calculate Temperature
	BOLTZ = MOLAR_GAS_CONSTANT_R
# calculate degrees of freedom in the system.  This used to calculate simulation temperature
	particles = mySimulation.system.getNumParticles()
	constraints = mySimulation.system.getNumConstraints()
	forces = [mySimulation.system.getForce(i) for i in range(mySimulation.system.getNumForces())]
	usesCMMR = any(isinstance(f, CMMotionRemover) for f in forces)
	DOF = 3*particles-constraints 
	if usesCMMR: DOF -= 3

# Set minimum error in sim temperature error reltive to set temperature
	Temp_tolerance = 0.008

# reset time and step count to zero
	mySimulation.context.setTime(0)
	mySimulation.currentStep = 0   #  Does this work?
# Set both integrator  temperature
	mySimulation.integrator.setTemperature(Sim_Temp)

# This function set's velocities independent of integrator.  Likely not appropriate for reaching thermodynamicequilibrium.  Coment out for production runs.
#	mySimulation.context.setVelocitiesToTemperature(Sim_Temp,int(time.time()))



# Begin Cycle of testing simulation temperature to set point temperature
	Equlib_Cycle = 0
	Current_Temp=mySimulation.context.getState(getEnergy=True).getKineticEnergy()/(0.5*DOF*BOLTZ)
	print ('  Current simulation temperature = {:.1f}K'.format(Current_Temp/kelvin))
	Temp_Err = Current_Temp/Sim_Temp - 1

	while (abs(Temp_Err) > Temp_tolerance):
		Equlib_Cycle +=1
		print ("     Starting Thermal Equlib Step: {}".format(Equlib_Cycle))
		
		ExecuteSimulation(mySimulation, Dwell)
		
		Current_Temp=mySimulation.context.getState(getEnergy=True).getKineticEnergy()/(0.5*DOF*BOLTZ)
		Temp_Err = Current_Temp/Sim_Temp - 1
		print("         Simulation T=","{:.1f}".format(Current_Temp/kelvin),"K.  Temp Error = ", "{:+.2f}".format(Temp_Err*100),"%")
	print ('    Equilibrated Temp = {:.1f}K'.format(Current_Temp/kelvin))
	print ('	Thermal Equilibrated for ', mySimulation.currentStep, ' steps')


	return
# End SetSimTemp



#########################################################	
def CollectTrajectory(basefile,mySimulation,T_Collect,Steps,ReportInt,AnimateInt):
# Simulate at T and save .dcd amd .log file



	
# Reset simulation
	mySimulation.context.setTime(0)
	mySimulation.currentStep = 0  #This does not work!
	

# Initialize Reporters
	mySimulation.reporters.append(DCDReporter(basefile+'test_traj.dcd' , AnimateInt))
	mySimulation.reporters.append(StateDataReporter(basefile + '_test.log', ReportInt, step=True, 
		time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, 
		temperature=True,  density=True, volume=True, progress=True, 
		remainingTime=True, speed=True, totalSteps=Steps, separator='\t'))
	

# Start the iterated simulation!
	print('Running Trajectory collection: ')
	print(  'Trajectory length: ', Steps)
	print('Startingtime: ',time.ctime())

	ExecuteSimulation(mySimulation, Steps)


#stop reporting progress and trajectory
	mySimulation.reporters.pop()
	mySimulation.reporters.pop()

					
# End CollectTempTrajectory

########################
# Code for general simulation.  if step= 0 then minimize

def ExecuteSimulation(sim,steps):

	if steps==0:
		sim.minimizeEnergy()
	else:
		sim.step(steps)


	return(True)


#####################
			
if __name__ == '__main__':  
	main()						
									
				
	