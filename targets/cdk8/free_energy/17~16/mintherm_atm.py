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

jobname = '17~16'


displ = [66.17,3.13,0.41]
displacement      = [  displ[i] for i in range(3) ] * angstrom
lig1_restr_offset = [  0.       for i in range(3) ] * angstrom
lig2_restr_offset = [  displ[i] for i in range(3) ] * angstrom

lig1_atoms = [5972,5973,5974,5975,5976,5977,5978,5979,5980,5981,5982,5983,5984,5985,5986,5987,5988,5989,5990,5991,5992,5993,5994,5995,5996,5997,5998,5999,6000,6001,6002,6003,6004,6005,6006,6007,6008,6009,6010,6011,6012]
lig2_atoms = [6013,6014,6015,6016,6017,6018,6019,6020,6021,6022,6023,6024,6025,6026,6027,6028,6029,6030,6031,6032,6033,6034,6035,6036,6037,6038,6039,6040,6041,6042,6043,6044,6045,6046,6047,6048,6049,6050]
refatoms_lig1 = [4,9,12]
refatoms_lig2 = [4,9,12]
rcpt_cm_atoms = [489,505,859,1632,1652,1664,1685,1727,2625,2635,2668,2890,2900]
restrained_atoms = [ 9,21,33,55,72,84,105,117,137,159,175,197,216,227,238,253,277,292,316,332,347,359,378,398,413,434,449,456,467,489,505,512,536,543,557,578,585,602,618,639,661,671,693,717,739,751,758,780,792,804,826,838,859,869,888,910,927,946,961,968,982,989,1008,1019,1036,1047,1057,1068,1092,1107,1126,1136,1155,1174,1198,1213,1232,1254,1279,1285,1299,1315,1334,1345,1364,1381,1403,1419,1439,1458,1469,1486,1496,1508,1532,1554,1570,1594,1613,1632,1652,1664,1685,1695,1710,1727,1739,1758,1782,1799,1818,1837,1859,1879,1896,1920,1930,1941,1963,1973,1987,2009,2039,2045,2061,2078,2105,2111,2135,2142,2159,2175,2197,2208,2227,2246,2267,2284,2303,2322,2334,2341,2360,2377,2398,2417,2434,2444,2458,2482,2498,2517,2534,2558,2570,2589,2619,2625,2635,2649,2668,2687,2703,2720,2727,2742,2757,2763,2778,2802,2809,2833,2849,2871,2890,2900,2912,2929,2936,2956,2966,2990,3009,3029,3043,3062,3068,3087,3117,3123,3142,3152,3164,3183,3203,3209,3225,3241,3257,3271,3291,3315,3336,3360,3378,3384,3399,3418,3437,3456,3463,3473,3497,3514,3535,3549,3571,3581,3600,3612,3631,3655,3665,3684,3691,3702,3721,3741,3751,3766,3785,3804,3818,3829,3852,3858,3877,3897,3914,3925,3949,3966,3981,3993,4012,4034,4048,4059,4081,4087,4108,4125,4142,4154,4171,4190,4202,4226,4245,4265,4279,4295,4312,4319,4347,4353,4363,4375,4397,4409,4433,4448,4460,4479,4501,4523,4548,4554,4569,4586,4597,4611,4630,4647,4669,4681,4701,4725,4749,4763,4777,4798,4812,4826,4837,4848,4867,4886,4908,4929,4946,4961,4983,5000,5022,5038,5068,5074,5086,5097,5119,5129,5149,5166,5185,5204,5221,5243,5262,5281,5295,5312,5332,5338,5357,5379,5403,5422,5436,5447,5462,5479,5489,5506,5523,5543,5549,5570,5590,5609,5624,5644,5650,5677,5683,5697,5708,5720,5736,5756,5766,5773,5784,5801,5828,5834,5863,5869,5891,5915,5930,5950 ]

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
