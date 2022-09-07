#!/usr/bin/python
"""

InterfaceSimulation.py
	
K.E. Johnson, Pacific University
Summer 2019

	Purpose is to run simulation on a given surfactant system.
	Thermal processing is:
		1. Minimize
		2. heat to anneal temp
		3. initial anneal
		save chkpoint file
		
		Repeat count times
		   Read chkpoint file
		4. Anneal
		5. cool to simulate temp
		6. Equilibrate
		7. Collect
		
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
	global maxerrors
	maxerrors = 5

# Check for presence of the pdb file
	if len(sys.argv) != 4:
		print("Three input arguments are expected: The input pdb file for the simulations, the integer repeat count, and the GPU identifier ('0' or '1').")
		sys.exit(1)
		
	basefile = sys.argv[1]
	
	if not os.path.isfile(basefile+".pdb"):
		print("File {} does not exist.  Exiting!".format(basefile+".pdb"))
		sys.exit(1)
		
	Trajectory_Repeat = int(sys.argv[2])
	if Trajectory_Repeat < 1 or Trajectory_Repeat>10:
		print("Simulation repeat must be an integer from 1 to 10")
		sys.exit(1)
		
	global gpu_id
	gpu_id = sys.argv[3]
	if (gpu_id !='0') and (gpu_id !='1'):
		print("GPU id must be '0' or '1'.")
		sys.exit(1)
	
##################
#  Simulation Constants

# Simulation Step Size
	sim_step= 2.0 * femtoseconds

# temperature setpoints
	T_Simulation=300
	T_Anneal=325
	Ramp_Delta_T=1  # 1
	
# Dwell times
	Dwell_T_change = 100. *picoseconds # 200. 
	Dwell_init_anneal = 2. * nanoseconds  # 2.0
	Dwell_anneal = 1. * nanoseconds  # 2.0
	Dwell_equilibrate = 2. * nanoseconds  # 5.0
	
# Convert to steps	
	Dwell_T_steps = round(Dwell_T_change / sim_step)
	Dwell_init_steps = round(Dwell_init_anneal / sim_step)
	Dwell_anneal_steps = round(Dwell_anneal / sim_step)
	Dwell_equilibrate_steps = round(Dwell_equilibrate / sim_step)
	
# Trajectory Length & Repeat

	Trajectory_time = 4.0 * nanoseconds  # 5.0 
	Trajectory_steps = round (Trajectory_time / sim_step)
		
# Reporting 	
	AnimationPoints = 1000  #1000
	ReportingPoints = 200  #200
	FrameSavePoints = 1000  # 1000

	
	TrajReporterInterval=round(Trajectory_steps/ReportingPoints)
	TrajAnimationInterval=round(Trajectory_steps/AnimationPoints)
	TrajFrameSaveInterval=round(Trajectory_steps/FrameSavePoints)
	
	AnnealReporting = 1000
	
# These constants used to determine z- location of the interface.

	global Interface_Res
	Interface_Res='HOH'
	global 	Interface_Element
	Interface_Element = 8

# This value is the count of interface atoms.   estimate of the top monolayer at 124.  Shamay reference	
# Empirial -- from an analysis pf density it appears the correct value is 80
	global Interface_use
	Interface_use = 80
	

# Initialize the simulation 

	mySimulation=InitSimulation(basefile,T_Anneal,sim_step)


#  Here check if a checkpoint file has previously been saved.  If not prepare and save it. If yes, skip to the simulation loop

	if not os.path.isfile(basefile+".chk"):

					
# Need to write info file and a topology data file for future analysis

		writeInfo(basefile, mySimulation, T_Anneal,T_Simulation, Dwell_init_anneal, Dwell_anneal, Dwell_equilibrate, Trajectory_time, AnimationPoints, FrameSavePoints)
	
		writeTopo(basefile,mySimulation)

# minimize forces

		print('     Minimizing System Energy' )
			
		ExecuteSimulation(mySimulation, 0, 0)
					
		print ('     Minimization Complete')
					

# Set Temp & Equilibrate  Save a log for analysis
		print("Initial Anneal at {}K".format(T_Anneal))
		SetSimTemp(mySimulation,T_Anneal,Dwell_init_steps, logit=True, logfile=basefile + '_00_anneal')
		
# Save .chk file
		mySimulation.saveCheckpoint(basefile+'.chk')

	
# This secion takes the .chk file through the simulation a number of times depending on the commandline paameter
# Checks for existing binary trajector files to start at the correct file index


# Determine the next binary file number to write
	Out_start = 1
	while os.path.isfile(basefile+"_{:02d}_traj.framebin".format(Out_start)):
		Out_start +=1

	for simcount in range(Out_start,Out_start+Trajectory_Repeat):
		print("Starting data collection for number {} of {} trajectories".format(simcount-Out_start+1,Trajectory_Repeat))
		print("Loading checkpoint file")
		mySimulation.loadCheckpoint(basefile+'.chk')

# Important here to randomoize the system, otherwise all simulations will be identical.  Proceed to anneal.
		print("Randomizing simulation velocities at anneal temperature")
		mySimulation.context.setVelocitiesToTemperature(T_Anneal,int(time.time()))
		
		print("Annealing at {}K for {} steps".format(T_Anneal, Dwell_anneal_steps))	
		
			
		SetSimTemp(mySimulation,T_Anneal,Dwell_anneal_steps, logit=True, logfile=basefile+"_{:02d}_anneal".format(simcount))
		
		print("Cooling to {}K".format(T_Simulation))
		RampTemp(mySimulation,T_Anneal,T_Simulation,Ramp_Delta_T,Dwell_T_steps)
		
		print("Equilibrate Dwell for {} steps".format(Dwell_equilibrate_steps))
		SetSimTemp(mySimulation,T_Simulation,Dwell_equilibrate_steps, logit=True, logfile=basefile+"_{:02d}_equilibrate".format(simcount))

		print("Beginning Data Collection Section")
		CollectTrajectory(basefile+"_{:02d}".format(simcount), mySimulation, T_Simulation, Trajectory_steps, TrajReporterInterval, TrajAnimationInterval, TrajFrameSaveInterval)
		
			
	print("{} Simulations complete.  Congratulations!".format(Trajectory_Repeat))


	return(True)
	
#  END main



#########################################################	
def InitSimulation(FileBase,T_Start,step_size):
#Initialize and minimize simulation

	errorcount=0
	global maxerrors

	print("Initializing simulation.")	

	Sim_Temp=T_Start * kelvin

	
# Read in pdb file
	pdb_in=FileBase+".pdb"

	print ("    Loading initial pdb file")
	pdb = PDBFile(pdb_in)
	pdb.topology.setUnitCellDimensions((4.0, 4.0, 12.)*nanometers)
	
# Set up for system, NoCutof for Amoeba, Constrain only x-H bond lengths
	print ('    Loading Force Field')
	forcefield = ForceField('../AMOEBAff/solvents.xml', '../AMOEBAff/ions.xml', '../AMOEBAff/surfactants.xml')
	
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
	integrator_Friction = 2.0/picoseconds	# friction is 1.0 per ps, Future possibliltuy more frequent collisions
	
#Langevin Integrator. May consider experimental integrator in next iteration
#	integrator = LangevinIntegrator(Sim_Temp, integrator_Friction, step_size)
	integrator = LangevinMiddleIntegrator(Sim_Temp, integrator_Friction, step_size)

# setup to run CUDA at mixed precision. Single precision for forces, integrate @ double precision.
	print ('    Setting Platform to CUDA')
	platform = Platform.getPlatformByName('CUDA')
	properties = {'DeviceIndex': gpu_id, 'CudaPrecision': 'mixed'}

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
def writeInfo(basefile,Sim,Tannl, Ttraj, timeAnneal1,timeAnneal2,timeEqulib, timesim, animCount, frameCount):

	print("Writing simulation conditions information file.")
	with open(basefile+".siminfo","w") as f:
		
		
# Write file name and data started
		f.write("Simulation started: {}\n\n".format(datetime.datetime.now()))

# Write simulation parameters

# Write Step duration
		f.write("Simulation Step Size: {}\n\n".format(Sim.integrator.getStepSize()))

# write temperatures

		f.write("Temp Anneal: {}K\n".format(Tannl))
		f.write("Temp Trajectory: {}K\n\n".format(Ttraj))


# write anneam time ramp time, equilibrate time, simulation time

		f.write("Time initial anneal: {}\n".format(timeAnneal1))
		f.write("Time sampling anneal: {}\n".format(timeAnneal2))
		f.write("Time equilibrate: {}\n".format(timeEqulib))
		f.write("Trajectory duration: {}\n".format(timesim))
		f.write("Trajectory number of frames saved (binary): {}\n\n".format(frameCount))
		
# Information about reporters
		f.write("Animation time step: {}\n".format(timesim/animCount))
		f.write("Data Frame time step: {}\n\n".format(timesim/frameCount))
		
	return(True)

# End writeInfo

#########################################################
def writeTopo(basefile,Sim):
	
	print("Writing simulation topology information file.")
	with open(basefile+".simtopo","w") as f:

# write simulation dimensions
		f.write("Simulation Box dimensions: {}\n\n".format(Sim.topology.getUnitCellDimensions()))
		
# write Topology header
		f.write("atom index\tatom name\tElement\tatomic number\tmolar mass\tresidue index\tresidue name\n")

# Write Topology data, and in the same loop count the residues, and atoms
		
		atomlist={}
		counter=0

		for a in Sim.topology.atoms():
			counter +=1
			f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(a.index, a.name, a.element.symbol, a.element.atomic_number, a.element.mass.value_in_unit(dalton), a.residue.index, a.residue.name))
			if a.element.symbol in atomlist:
				atomlist[a.element.symbol] +=1
			else:
				atomlist[a.element.symbol] =1	
					
		f.write("\n")
# Write list of elements and count
		f.write("element\tcount\n")
		for key, value in sorted(atomlist.items()):
			f.write("{}\t{}\n".format(key,value))
		f.write("\nTotal Atom Count {}\n\n".format(counter))

# Count and list the number of each residue (molecule)		

		reslist={}
		counter = 0
		
		for r in Sim.topology.residues():
			counter +=1
			if r.name in reslist:
				reslist[r.name] +=1
			else:
				reslist[r.name] =1

		f.write("molecule\tcount\n")
		for key, value in sorted(reslist.items()):
			f.write("{}\t{}\n".format(key,value))
		f.write("\nTotal Molecule Count {}".format(counter))
		
	return(True)

# End writeTopo


#########################################################
def RampTemp(sim,Tempstart,Tempend,Tempstep,Dwell):

	print("Temperature ramp from {}K to {}K".format(Tempstart,Tempend))
	direction= 2*int(Tempstart<Tempend)-1 		
	for T in range(Tempstart,Tempend+direction,direction*Tempstep):
		SetSimTemp(sim,T,Dwell)
		
	print("Temperature ramp commplete")

	return(True)
# End TampTemp



#########################################################	
def SetSimTemp(mySimulation,Temp,Dwell,logit=False,logfile="log"):
# Bring simulation to desired Temperature

	errorcount=0
	global maxerrors

	print("Bringing Simulation to set Temperature: {} K".format(Temp))
	
	Sim_Temp = Temp*kelvin
	Equlib_Steps = 1000
	
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
		
		ExecuteSimulation(mySimulation, Equlib_Steps, Temp)
		
		Current_Temp=mySimulation.context.getState(getEnergy=True).getKineticEnergy()/(0.5*DOF*BOLTZ)
		Temp_Err = Current_Temp/Sim_Temp - 1
		print("         Simulation T=","{:.1f}".format(Current_Temp/kelvin),"K.  Temp Error = ", "{:+.2f}".format(Temp_Err*100),"%")
	print ('    Equilibrated Temp = {:.1f}K'.format(Current_Temp/kelvin))
	print ('	Thermal Equilibrated for ', mySimulation.currentStep, ' steps')

# Continue with simulation 'wait' more setps

	if Dwell>0:
		print ("    Simulate for {} additional steps.".format(Dwell))

# reset time and step count to zero
		mySimulation.context.setTime(0)
		mySimulation.currentStep = 0   

		if logit:
			mySimulation.reporters.append(StateDataReporter(logfile + '.log', round(0.02*Dwell), step=True, 
				time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, 
				temperature=True,  density=True, volume=True, progress=True, 
				remainingTime=True, speed=True, totalSteps=Dwell, separator='\t'))
		
		ExecuteSimulation(mySimulation, Dwell, Temp)

		if logit:
			mySimulation.reporters.pop()


	return
# End SetSimTemp



#########################################################	
def CollectTrajectory(basefile,mySimulation,T_Collect,Steps,ReportInt,AnimateInt,FrameInt):
# Simulate at T and save .dcd amd .log file


#inteface information
	global Interface_Res
	global Interface_Element
	global Interface_check
	global Interface_use

	
# Reset simulation
	mySimulation.context.setTime(0)
	mySimulation.currentStep = 0  #This does not work!
	

# Initialize Reporters
	mySimulation.reporters.append(DCDReporter(basefile+'_traj.dcd' , AnimateInt))
	mySimulation.reporters.append(StateDataReporter(basefile + '_sim.log', ReportInt, step=True, 
		time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, 
		temperature=True,  density=True, volume=True, progress=True, 
		remainingTime=True, speed=True, totalSteps=Steps, separator='\t'))
	mySimulation.reporters.append(InterfaceReporter(basefile+"_traj" , FrameInt, mySimulation, Interface_Res, Interface_Element, Interface_use))
	

# Start the iterated simulation!
	print('Running Trajectory collection: ',basefile)
	print(  'Trajectory length: ', Steps)
	print('Startingtime: ',time.ctime())

	ExecuteSimulation(mySimulation, Steps, T_Collect)


#stop reporting progress and trajectory
	mySimulation.reporters.pop()
	mySimulation.reporters.pop()
	mySimulation.reporters.pop()

					
# End CollectTempTrajectory

########################
# Code for general simulation.  if step= 0 then minimize

def ExecuteSimulation(sim,steps,temp):

	errorcount = 0
	global maxerrors

	done=False
	while errorcount <= maxerrors and not done:
		try:
			if steps==0:
				sim.minimizeEnergy()
			else:
				sim.step(steps)
			done=True
		except:
			print('Simulation error in Simuation call at ', time.ctime())
			traceback.print_exc()
			errorcount += 1
			if errorcount > maxerrors:
				print('Maximum number of errors reached.  Terminating!')
				break
#					sys.exit(1)
			else:
				print('Resetting Simulation')
#  first get the positions (translated into the box).
# Then reset, set positons, set velocities to new randomized temperature/ Boltzman distribution
				statePositions = sim.context.getState(getPositions=True, enforcePeriodicBox=True)
				sim.context.reinitialize()
				sim.context.setState(statePositions)
				if steps>0:
					sim.context.setVelocitiesToTemperature(temp*kelvin ,int(time.time()))
					steps=steps-sim.currentStep
				print('  RESTARTING Simulation')

	return(True)


"""
Code for a custom OpenMM reporter
Purpose is to write a binary datafile with:
	xpos, ypos, (zpos - interface location), xforce, yforce, zforce

To use numpy array functions for efficiency

Interface location requires that we initially identify the atom indicies for the top 1/4th of the water O atoms.
	At each report we find the z positions of the idetified O atoms.  Sort to determine the top monolayer, then calculate the median value of that list.  That becomes the defined z = 0 / interface of the simulation report.



"""

class InterfaceReporter(object):
	def __init__(self, file, reportInterval, simsystem, resname, element_number, test_number):
		self._out = open(file+'.framebin', 'wb')
		self._intfc = open(file+'.intfloc','w')
		self._reportInterval = reportInterval
# here create a matrix for the frame data.  Filled in report
		natoms= simsystem.topology.getNumAtoms()
		self._frame=np.empty([natoms,6],float)

# This arrray used to find the atoms to test for interface location
# First index array and add z positions		
		z_list=np.zeros([natoms,3],float)
		z_list[:,0]=np.arange(0,natoms)
		

		pos = simsystem.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(nanometers)
		z_list[:,1]=pos[:,2]

#next determine if the atom is Oxygen in water
		for a in simsystem.topology.atoms():
			if a.residue.name==resname and a.element.atomic_number==element_number:
				z_list[a.index,2]=1.
		

# extract O-water atoms and sort and reverse order -- assume interface of interest has greater z
		z_list=z_list[z_list[:,2]>0][:,0:2]
		

#Extract the atom numbers in an array	
		self._testarray=np.sort(z_list[:,0].astype(int))
			
# Save the number of atoms to use in interface zero dertermination
		self._interfacecount = test_number
	
		

	def __del__(self):
		self._out.close()
		self._intfc.close()

	def describeNextReport(self, simulation):
		steps = self._reportInterval - simulation.currentStep%self._reportInterval
# boolian valuess indicate: position, velocity, force, enforce periodic boundaries
		return (steps, True, False, True, True)

	def report(self, simulation, state):
# get raw frame data		
		self._frame[:,3:6]=state.getForces(asNumpy=True).value_in_unit(kilojoules/mole/nanometer)
			
		self._frame[:,0:3]=state.getPositions(asNumpy=True).value_in_unit(nanometers)

# now extract the z positions of the top Oxygen-water atoms, calculate median, adjust the z position in frame accordingly

		adjust=np.mean((np.sort(self._frame[self._testarray,2]))[-self._interfacecount:])
		self._intfc.write('{}\n'.format(adjust))
		self._frame[:,2] -= adjust
		
# Write the binary data
		self._out.write(self._frame)


# END custom reporter

#####################
			
if __name__ == '__main__':  
	main()						
									
				
	