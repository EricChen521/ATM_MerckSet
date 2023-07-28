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

jobname = 'SHP099-1-Example-7~Example-2'


displ = [-7.06,49.83,4.57]
displacement      = [  displ[i] for i in range(3) ] * angstrom
lig1_restr_offset = [  0.       for i in range(3) ] * angstrom
lig2_restr_offset = [  displ[i] for i in range(3) ] * angstrom

lig1_atoms = [8363,8364,8365,8366,8367,8368,8369,8370,8371,8372,8373,8374,8375,8376,8377,8378,8379,8380,8381,8382,8383,8384,8385,8386,8387,8388,8389,8390,8391,8392,8393,8394,8395,8396,8397,8398,8399,8400,8401,8402,8403,8404,8405]
lig2_atoms = [8406,8407,8408,8409,8410,8411,8412,8413,8414,8415,8416,8417,8418,8419,8420,8421,8422,8423,8424,8425,8426,8427,8428,8429,8430,8431,8432,8433,8434,8435,8436,8437,8438,8439,8440,8441,8442,8443,8444,8445,8446,8447,8448,8449,8450,8451]
refatoms_lig1 = [0,13,11]
refatoms_lig2 = [0,13,11]
rcpt_cm_atoms = [1646,1660,1671,1686,1710,1734,1754,3383,3877,3892,3907,3942,3956,7779,7785]
restrained_atoms = [ 9,33,57,81,101,126,132,146,165,179,186,202,217,227,242,256,275,294,313,327,351,358,374,386,393,404,424,443,453,485,491,502,524,535,557,563,570,582,602,616,635,646,662,686,710,724,731,741,757,771,788,807,829,848,865,879,893,900,912,933,954,966,985,1006,1013,1020,1035,1057,1077,1087,1101,1120,1130,1145,1164,1180,1197,1218,1239,1256,1271,1288,1305,1312,1329,1348,1370,1385,1407,1421,1428,1440,1456,1475,1490,1509,1531,1560,1566,1585,1599,1610,1620,1640,1646,1660,1671,1686,1710,1734,1754,1771,1778,1795,1814,1825,1832,1854,1869,1879,1894,1916,1935,1954,1968,1983,2005,2012,2034,2051,2058,2069,2089,2108,2124,2148,2163,2174,2191,2202,2227,2233,2240,2252,2272,2288,2307,2318,2334,2358,2372,2379,2391,2403,2425,2432,2447,2458,2472,2484,2491,2513,2524,2546,2562,2576,2593,2609,2626,2645,2669,2680,2697,2712,2731,2753,2774,2786,2802,2809,2816,2823,2838,2862,2882,2894,2905,2924,2938,2950,2969,2985,3000,3017,3038,3060,3082,3104,3110,3127,3143,3158,3172,3191,3198,3212,3228,3247,3264,3283,3305,3330,3336,3355,3369,3383,3397,3421,3440,3454,3464,3474,3489,3508,3523,3534,3558,3574,3598,3613,3632,3643,3665,3684,3694,3709,3723,3737,3749,3771,3787,3809,3826,3833,3853,3877,3892,3907,3927,3942,3956,3975,3992,4009,4026,4041,4052,4074,4093,4112,4133,4144,4168,4190,4205,4212,4229,4253,4270,4285,4299,4321,4335,4357,4371,4395,4416,4438,4452,4471,4498,4504,4524,4536,4553,4567,4591,4607,4623,4642,4659,4671,4678,4698,4704,4718,4741,4747,4763,4774,4786,4807,4826,4840,4850,4864,4883,4902,4927,4933,4948,4968,4983,4997,5019,5030,5044,5058,5069,5099,5105,5127,5149,5160,5181,5200,5210,5224,5241,5248,5259,5278,5295,5309,5323,5339,5353,5365,5385,5409,5433,5450,5466,5486,5503,5518,5532,5543,5567,5583,5602,5618,5635,5649,5663,5685,5700,5716,5731,5755,5762,5784,5795,5817,5828,5844,5866,5887,5919,5925,5937,5952,5973,5983,6002,6024,6039,6060,6067,6083,6100,6124,6140,6164,6178,6194,6216,6231,6242,6252,6262,6279,6291,6312,6326,6345,6369,6384,6403,6425,6444,6455,6477,6493,6500,6517,6524,6538,6552,6567,6591,6605,6621,6645,6662,6683,6700,6720,6744,6758,6790,6796,6808,6825,6832,6856,6862,6873,6893,6899,6906,6913,6929,6948,6960,6980,6999,7014,7029,7045,7062,7079,7101,7118,7133,7144,7163,7180,7192,7202,7217,7223,7239,7255,7271,7288,7299,7310,7320,7327,7346,7353,7377,7391,7398,7412,7432,7451,7467,7486,7498,7517,7536,7555,7567,7586,7605,7629,7644,7666,7673,7689,7701,7712,7724,7743,7755,7779,7785,7807,7821,7840,7857,7874,7890,7914,7925,7942,7966,7977,7984,8001,8017,8034,8048,8063,8073,8090,8111,8135,8155,8174,8195,8212,8222,8238,8255,8272,8293,8312,8327,8341 ]

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
