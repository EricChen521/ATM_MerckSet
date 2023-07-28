from __future__ import print_function

import openmm as mm
from openmm.app import *
from openmm import *
from simtk.unit import *
from sys import stdout
import os, re,time, shutil, math
from datetime import datetime

from atmmetaforce import *


# the multiple-time step integrator does not have a setTemperature() method
def setTemperature(self, temperature):
    self.setGlobalVariableByName('kT', MOLAR_GAS_CONSTANT_R*temperature)
MTSLangevinIntegrator.setTemperature = setTemperature

print("Started at: " + str(time.asctime()))
start=datetime.now()

jobname = '5a~5n'


displ = [40.48,1.58,4.62]
displacement      = [  displ[i] for i in range(3) ] * angstrom
lig1_restr_offset = [  0.       for i in range(3) ] * angstrom
lig2_restr_offset = [  displ[i] for i in range(3) ] * angstrom

lig1_atoms = [3297,3298,3299,3300,3301,3302,3303,3304,3305,3306,3307,3308,3309,3310,3311,3312,3313,3314,3315,3316,3317,3318,3319,3320,3321,3322,3323,3324,3325,3326,3327,3328,3329,3330,3331,3332,3333]
lig2_atoms = [3334,3335,3336,3337,3338,3339,3340,3341,3342,3343,3344,3345,3346,3347,3348,3349,3350,3351,3352,3353,3354,3355,3356,3357,3358,3359,3360,3361,3362,3363,3364,3365,3366,3367,3368,3369,3370,3371,3372,3373,3374,3375,3376,3377]
refatoms_lig1 = [15,14,0]
refatoms_lig2 = [5,6,15]
rcpt_cm_atoms = [1307,1324,1331,1598,1732,1753,1773,1856,1898,2878]
restrained_atoms = [ 9,16,30,49,68,87,99,118,137,143,155,167,189,204,224,241,252,268,283,298,313,330,347,358,372,388,412,427,444,468,480,487,494,511,521,528,535,554,574,588,612,633,647,666,685,707,726,743,765,781,792,806,828,850,869,893,908,932,953,967,984,1008,1032,1054,1069,1085,1096,1111,1126,1140,1157,1171,1188,1198,1212,1227,1251,1268,1287,1307,1324,1331,1350,1356,1376,1392,1406,1416,1435,1454,1471,1493,1500,1520,1532,1547,1571,1588,1598,1619,1638,1645,1652,1669,1689,1696,1706,1713,1732,1753,1773,1783,1798,1812,1823,1834,1856,1867,1881,1898,1919,1935,1956,1963,1982,1989,1996,2003,2017,2024,2043,2049,2065,2082,2104,2116,2140,2151,2162,2183,2202,2213,2230,2254,2271,2290,2309,2329,2340,2364,2380,2394,2413,2420,2442,2453,2473,2492,2509,2529,2540,2562,2579,2589,2606,2625,2639,2645,2652,2669,2686,2697,2713,2727,2734,2766,2772,2783,2799,2813,2820,2839,2849,2868,2878,2893,2914,2930,2949,2970,2994,3001,3016,3033,3043,3072,3078,3093,3114,3133,3152,3166,3187,3204,3223,3240,3272,3278 ]

# load system
prmtop = AmberPrmtopFile("complex.prmtop")
inpcrd = AmberInpcrdFile("complex.inpcrd")
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.9*nanometer,
                             constraints=HBonds)

# load the ATM Meta Force facility.
atm_utils = ATMMetaForceUtils(system)

# Vsite restraints
lig1_cm_atoms = lig1_atoms
lig2_cm_atoms = lig2_atoms

# force constant for Vsite CM-CM restraint
kf = 0.0 * kilocalorie_per_mole/angstrom**2
# radius of Vsite sphere
r0 = 5.0 * angstrom #radius of Vsite sphere
atm_utils.addRestraintForce(lig_cm_particles = lig1_cm_atoms,
                            rcpt_cm_particles = rcpt_cm_atoms,
                            kfcm = kf,
                            tolcm = r0,
                            offset = lig1_restr_offset)

atm_utils.addRestraintForce(lig_cm_particles = lig2_cm_atoms,
                            rcpt_cm_particles = rcpt_cm_atoms,
                            kfcm = kf,
                            tolcm = r0,
                            offset = lig2_restr_offset)

# alignment restraint

lig1_ref_atoms  = [ refatoms_lig1[i]+lig1_atoms[0] for i in range(3)]
lig2_ref_atoms  = [ refatoms_lig2[i]+lig2_atoms[0] for i in range(3)]

atm_utils.addAlignmentForce(liga_ref_particles = lig1_ref_atoms,
                            ligb_ref_particles = lig2_ref_atoms,
                            kfdispl = 2.5 * kilocalorie_per_mole/angstrom**2,
                            ktheta =  50.0 * kilocalorie_per_mole,
                            kpsi =  50.0 * kilocalorie_per_mole,
                            offset = lig2_restr_offset)

# restrain positions of specified atoms
fc = 25.0 * kilocalorie_per_mole/angstrom**2
tol = 1.5 * angstrom
posrestr_atoms = []
last_not_solvent_atom = lig2_atoms[len(lig2_atoms)-1]
for at in prmtop.topology.atoms():
    if at.index <= last_not_solvent_atom:
        posrestr_atoms.append(at.index)
atm_utils.addPosRestraints(posrestr_atoms, inpcrd.positions, fc, tol)

#places non-bonded force in a separate group
nonbonded_force_group = 1
atm_utils.setNonbondedForceGroup(nonbonded_force_group)

# Set up Langevin integrator
initial_temperature = 50 * kelvin
final_temperature = 298.15 * kelvin
temperature = initial_temperature
frictionCoeff = 0.5 / picosecond
MDstepsize = 0.001 * picosecond
barostat = MonteCarloBarostat(1*bar, final_temperature)
saved_barostat_frequency = barostat.getFrequency()
barostat.setFrequency(900000000)#disabled
system.addForce(barostat)
integrator = MTSLangevinIntegrator(temperature, frictionCoeff, MDstepsize, [(0,1), (nonbonded_force_group,1)])
integrator.setConstraintTolerance(0.00001)

# platform_name = 'OpenCL'
platform_name = 'CUDA'
platform = Platform.getPlatformByName(platform_name)
properties = {}
properties["Precision"] = "mixed"

simulation = Simulation(prmtop.topology, system, integrator,platform, properties)
print ("Using platform %s" % simulation.context.getPlatform().getName())
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

print("Potential energy before minimization =",
simulation.context.getState(getEnergy = True).getPotentialEnergy())

print("Energy minimizing the system ...")
simulation.minimizeEnergy()

print("Potential energy after minimization =",
simulation.context.getState(getEnergy = True).getPotentialEnergy())

print("Thermalization ...")

totalSteps = 100000
steps_per_cycle = 5000
number_of_cycles = int(totalSteps/steps_per_cycle)
delta_temperature = (final_temperature - initial_temperature)/number_of_cycles
simulation.reporters.append(
    StateDataReporter(stdout, steps_per_cycle,
     step=True, potentialEnergy = True, temperature=True, volume=True))

# MD with temperature ramp
for i in range(number_of_cycles):
    simulation.step(steps_per_cycle)
    # prepare system for new temperature
    temperature = temperature + delta_temperature
    integrator.setTemperature(temperature)

print("NPT equilibration ...")
barostat.setFrequency(saved_barostat_frequency)
# MD at constant pressure
for i in range(number_of_cycles):
    simulation.step(steps_per_cycle)

print("NVT equilibration ...")
barostat.setFrequency(900000000)#disabled
for i in range(number_of_cycles):
    simulation.step(steps_per_cycle)

print( "SaveState: equil.xml.")
# saves a checkpoint file and a pdb file
simulation.saveState('equil.xml')
positions = simulation.context.getState(getPositions=True).getPositions()
boxsize = simulation.context.getState().getPeriodicBoxVectors()
simulation.topology.setPeriodicBoxVectors(boxsize)
with open('equil.pdb', 'w') as output:
  PDBFile.writeFile(simulation.topology, positions, output)

end=datetime.now()
elapsed=end - start
print("elapsed time for equilbration ="+str(
    elapsed.seconds+elapsed.microseconds*1e-6)+"s")
