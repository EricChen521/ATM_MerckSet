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

jobname = '23~14'


displ = [66.17,3.13,0.41]
displacement      = [  displ[i] for i in range(3) ] * angstrom
lig1_restr_offset = [  0.       for i in range(3) ] * angstrom
lig2_restr_offset = [  displ[i] for i in range(3) ] * angstrom

lig1_atoms = [5972,5973,5974,5975,5976,5977,5978,5979,5980,5981,5982,5983,5984,5985,5986,5987,5988,5989,5990,5991,5992,5993,5994,5995,5996,5997,5998,5999,6000,6001,6002,6003,6004,6005,6006,6007,6008,6009,6010,6011,6012,6013,6014,6015,6016,6017,6018,6019,6020]
lig2_atoms = [6021,6022,6023,6024,6025,6026,6027,6028,6029,6030,6031,6032,6033,6034,6035,6036,6037,6038,6039,6040,6041,6042,6043,6044,6045,6046,6047,6048,6049,6050,6051,6052,6053,6054]
refatoms_lig1 = [4,9,10]
refatoms_lig2 = [4,9,12]
rcpt_cm_atoms = [489,505,859,1632,1652,1664,1685,1727,2625,2635,2668,2890,2900]
restrained_atoms = [ 9,21,33,55,72,84,105,117,137,159,175,197,216,227,238,253,277,292,316,332,347,359,378,398,413,434,449,456,467,489,505,512,536,543,557,578,585,602,618,639,661,671,693,717,739,751,758,780,792,804,826,838,859,869,888,910,927,946,961,968,982,989,1008,1019,1036,1047,1057,1068,1092,1107,1126,1136,1155,1174,1198,1213,1232,1254,1279,1285,1299,1315,1334,1345,1364,1381,1403,1419,1439,1458,1469,1486,1496,1508,1532,1554,1570,1594,1613,1632,1652,1664,1685,1695,1710,1727,1739,1758,1782,1799,1818,1837,1859,1879,1896,1920,1930,1941,1963,1973,1987,2009,2039,2045,2061,2078,2105,2111,2135,2142,2159,2175,2197,2208,2227,2246,2267,2284,2303,2322,2334,2341,2360,2377,2398,2417,2434,2444,2458,2482,2498,2517,2534,2558,2570,2589,2619,2625,2635,2649,2668,2687,2703,2720,2727,2742,2757,2763,2778,2802,2809,2833,2849,2871,2890,2900,2912,2929,2936,2956,2966,2990,3009,3029,3043,3062,3068,3087,3117,3123,3142,3152,3164,3183,3203,3209,3225,3241,3257,3271,3291,3315,3336,3360,3378,3384,3399,3418,3437,3456,3463,3473,3497,3514,3535,3549,3571,3581,3600,3612,3631,3655,3665,3684,3691,3702,3721,3741,3751,3766,3785,3804,3818,3829,3852,3858,3877,3897,3914,3925,3949,3966,3981,3993,4012,4034,4048,4059,4081,4087,4108,4125,4142,4154,4171,4190,4202,4226,4245,4265,4279,4295,4312,4319,4347,4353,4363,4375,4397,4409,4433,4448,4460,4479,4501,4523,4548,4554,4569,4586,4597,4611,4630,4647,4669,4681,4701,4725,4749,4763,4777,4798,4812,4826,4837,4848,4867,4886,4908,4929,4946,4961,4983,5000,5022,5038,5068,5074,5086,5097,5119,5129,5149,5166,5185,5204,5221,5243,5262,5281,5295,5312,5332,5338,5357,5379,5403,5422,5436,5447,5462,5479,5489,5506,5523,5543,5549,5570,5590,5609,5624,5644,5650,5677,5683,5697,5708,5720,5736,5756,5766,5773,5784,5801,5828,5834,5863,5869,5891,5915,5930,5950 ]

# define the thermodynamic/alchemical state
# the system is prepared at the alchemical intermediate state at lambda=0
temperature = 298.15 * kelvin
lmbd = 0.50
lambda1 = lmbd
lambda2 = lmbd
alpha = 0.0 / kilocalorie_per_mole
u0 = 0.0 * kilocalorie_per_mole
w0coeff = 0.0 * kilocalorie_per_mole
umsc =  1000.0 * kilocalorie_per_mole
ubcore = 500.0 * kilocalorie_per_mole
acore = 0.062500
direction = 1

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

# restrain the positions of specified atoms
fc = 25.0 * kilocalorie_per_mole/angstrom**2
tol = 1.5 * angstrom
atm_utils.addPosRestraints(restrained_atoms, inpcrd.positions, fc, tol)

# create ATM Force (direction is 1 by default)
atmforcegroup = 2
nonbonded_force_group = 1
atm_utils.setNonbondedForceGroup(nonbonded_force_group)
atmvariableforcegroups = [nonbonded_force_group]
atmforce = ATMMetaForce(
    lambda1, lambda2,  alpha * kilojoules_per_mole, u0/kilojoules_per_mole,
     w0coeff/kilojoules_per_mole, umsc/kilojoules_per_mole,
      ubcore/kilojoules_per_mole, acore, direction, atmvariableforcegroups )
# adds all atoms to the force with zero displacement
for at in prmtop.topology.atoms():
    atmforce.addParticle(at.index, 0., 0., 0.)
# the ligand atoms get displaced, ligand 1 from binding site to the solvent bulk
# and ligand 2 from the bulk solvent to the binding site
for i in lig1_atoms:
    atmforce.setParticleParameters(i, i,
    displ[0] * angstrom, displ[1] * angstrom, displ[2] * angstrom)
for i in lig2_atoms:
    atmforce.setParticleParameters(i, i,
    -displ[0] * angstrom, -displ[1] * angstrom, -displ[2] * angstrom)

atmforce.setForceGroup(atmforcegroup)
system.addForce(atmforce)

# Set up Langevin integrator
frictionCoeff = 0.5 / picosecond
MDstepsize = 0.001 * picosecond
barostat = MonteCarloBarostat(1*bar, temperature)
barostat.setFrequency(900000000)#disabled
system.addForce(barostat)
integrator = MTSLangevinIntegrator(temperature, frictionCoeff, MDstepsize, [(0,1), (atmforcegroup,1)])
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

state = simulation.context.getState(getEnergy = True, groups = {0,atmforcegroup})
pote = state.getPotentialEnergy()

print( "LoadState mdlambda.xml")
simulation.loadState('mdlambda.xml')

#override ATM parameters
simulation.context.setParameter(atmforce.Lambda1(), lambda1)
simulation.context.setParameter(atmforce.Lambda2(), lambda2)
simulation.context.setParameter(atmforce.Alpha(), alpha *kilojoules_per_mole)
simulation.context.setParameter(atmforce.U0(), u0 /kilojoules_per_mole)
simulation.context.setParameter(atmforce.W0(), w0coeff /kilojoules_per_mole)
simulation.context.setParameter(atmforce.Umax(), umsc /kilojoules_per_mole)
simulation.context.setParameter(atmforce.Ubcore(), ubcore /kilojoules_per_mole)
simulation.context.setParameter(atmforce.Acore(), acore)
simulation.context.setParameter(atmforce.Direction(), direction)

print("Potential energy =",
 simulation.context.getState(getEnergy = True, groups = {0,atmforcegroup}).getPotentialEnergy())

totalSteps = 150000
steps_per_cycle = 1000
simulation.reporters.append(
    StateDataReporter(stdout, steps_per_cycle,
     step=True, potentialEnergy = True, temperature=True, volume=True))

print("Equilibration at lambda=1/2...")
simulation.step(totalSteps)

print( "SaveState: complex_0.xml.")

# saves a checkpoint file and a pdb file at lambda=1/2.
# for historical reasons it has a "0" in the file name
#FIXME: the name is hard coded.
simulation.saveState('complex_0.xml')
positions = simulation.context.getState(getPositions=True).getPositions()
boxsize = simulation.context.getState().getPeriodicBoxVectors()
simulation.topology.setPeriodicBoxVectors(boxsize)
with open('complex_0.pdb', 'w') as output:
  PDBFile.writeFile(simulation.topology, positions, output)

end=datetime.now()
elapsed=end - start
print("elapsed time for equilbration ="+str(
    elapsed.seconds+elapsed.microseconds*1e-6)+"s")
