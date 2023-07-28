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

jobname = 'CHEMBL3402744-300~CHEMBL3402745-200'


displ = [1.91,43.96,7.25]
displacement      = [  displ[i] for i in range(3) ] * angstrom
lig1_restr_offset = [  0.       for i in range(3) ] * angstrom
lig2_restr_offset = [  displ[i] for i in range(3) ] * angstrom

lig1_atoms = [4577,4578,4579,4580,4581,4582,4583,4584,4585,4586,4587,4588,4589,4590,4591,4592,4593,4594,4595,4596,4597,4598,4599,4600,4601,4602,4603,4604,4605,4606,4607,4608,4609,4610,4611,4612,4613,4614,4615,4616,4617,4618]
lig2_atoms = [4619,4620,4621,4622,4623,4624,4625,4626,4627,4628,4629,4630,4631,4632,4633,4634,4635,4636,4637,4638,4639,4640,4641,4642,4643,4644,4645,4646,4647,4648,4649,4650,4651,4652,4653,4654,4655,4656,4657,4658,4659,4660,4661,4662,4663]
refatoms_lig1 = [15,18,10]
refatoms_lig2 = [16,19,11]
rcpt_cm_atoms = [347,366,708,1467,1494,1500,1521,1538,1560,1577,1584,2303,2327,2525,2535,2594,2657,2678]
restrained_atoms = [ 5,24,40,57,67,83,100,117,133,149,168,183,189,200,211,230,249,265,282,302,316,331,347,366,373,397,404,421,441,448,459,475,496,513,520,534,553,572,584,598,610,617,639,661,680,697,708,718,734,756,767,786,800,824,843,857,869,888,895,910,926,937,954,974,993,1007,1022,1029,1048,1067,1084,1106,1118,1138,1149,1174,1180,1194,1210,1229,1240,1259,1278,1285,1304,1315,1334,1358,1369,1384,1391,1410,1416,1435,1451,1467,1494,1500,1521,1538,1560,1577,1584,1596,1615,1639,1653,1673,1692,1716,1730,1745,1759,1776,1798,1804,1818,1834,1856,1868,1887,1906,1913,1933,1940,1959,1976,1992,2002,2024,2031,2048,2070,2091,2110,2120,2131,2153,2175,2195,2211,2228,2252,2264,2283,2293,2303,2327,2341,2352,2369,2388,2400,2415,2437,2457,2471,2487,2509,2525,2535,2548,2568,2575,2594,2604,2628,2640,2657,2678,2690,2712,2727,2748,2769,2780,2796,2813,2827,2849,2863,2870,2880,2902,2929,2935,2951,2973,2997,3014,3024,3043,3058,3069,3088,3105,3119,3136,3158,3178,3192,3206,3228,3239,3251,3267,3291,3302,3322,3329,3345,3364,3383,3407,3422,3441,3458,3472,3496,3503,3521,3535,3541,3570,3576,3588,3604,3618,3632,3652,3664,3683,3697,3713,3734,3753,3772,3789,3796,3820,3844,3863,3882,3907,3913,3928,3949,3968,3974,3994,4000,4019,4040,4055,4071,4088,4107,4129,4140,4164,4189,4195,4217,4227,4242,4259,4291,4297,4308,4328,4339,4354,4373,4389,4400,4424,4443,4454,4464,4483,4503,4514,4528,4548,4567 ]

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
