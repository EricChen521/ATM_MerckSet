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

jobname = 'CHEMBL1077227~CHEMBL1086410'


displ = [-12.08,0.62,51.25]
displacement      = [  displ[i] for i in range(3) ] * angstrom
lig1_restr_offset = [  0.       for i in range(3) ] * angstrom
lig2_restr_offset = [  displ[i] for i in range(3) ] * angstrom

lig1_atoms = [5474,5475,5476,5477,5478,5479,5480,5481,5482,5483,5484,5485,5486,5487,5488,5489,5490,5491,5492,5493,5494,5495,5496,5497,5498,5499,5500,5501,5502,5503,5504,5505,5506,5507,5508,5509,5510,5511,5512,5513,5514,5515,5516,5517,5518,5519,5520,5521,5522,5523,5524,5525,5526,5527,5528,5529,5530]
lig2_atoms = [5531,5532,5533,5534,5535,5536,5537,5538,5539,5540,5541,5542,5543,5544,5545,5546,5547,5548,5549,5550,5551,5552,5553,5554,5555,5556,5557,5558,5559,5560,5561,5562,5563,5564,5565,5566,5567,5568,5569,5570,5571,5572,5573,5574,5575,5576,5577,5578,5579,5580,5581,5582,5583,5584,5585,5586,5587,5588,5589,5590]
refatoms_lig1 = [7,17,3]
refatoms_lig2 = [7,18,3]
rcpt_cm_atoms = [1564,1579,1586,1601,1625,1644,1729,1783,1809,1828,1864,1891,3146,3165,3202,3209]
restrained_atoms = [ 9,16,38,52,71,88,104,120,136,160,171,203,209,229,243,262,272,287,311,333,343,354,364,381,392,411,427,442,453,473,479,495,519,541,556,572,583,599,623,637,644,651,670,680,692,714,725,736,760,782,796,817,831,851,863,880,896,916,923,933,944,958,980,997,1016,1028,1044,1065,1089,1100,1116,1132,1151,1157,1176,1195,1207,1222,1238,1257,1274,1281,1302,1316,1327,1341,1360,1380,1390,1411,1418,1435,1449,1456,1470,1477,1499,1513,1533,1547,1564,1579,1586,1601,1625,1644,1650,1664,1679,1694,1715,1729,1753,1768,1783,1803,1809,1828,1838,1845,1864,1891,1897,1921,1935,1954,1971,1988,2007,2027,2042,2064,2083,2097,2109,2123,2130,2144,2159,2179,2190,2206,2228,2244,2255,2274,2293,2308,2327,2348,2362,2377,2392,2411,2431,2443,2462,2481,2503,2509,2520,2531,2543,2559,2570,2585,2609,2628,2645,2662,2682,2694,2714,2720,2744,2758,2780,2804,2811,2827,2846,2865,2887,2894,2913,2928,2943,2962,2976,2992,3010,3024,3046,3058,3073,3089,3110,3127,3146,3165,3180,3202,3209,3219,3229,3251,3275,3289,3303,3313,3323,3337,3356,3373,3387,3397,3418,3429,3440,3464,3475,3492,3503,3519,3539,3550,3566,3580,3599,3616,3633,3655,3670,3684,3698,3717,3729,3736,3751,3766,3785,3801,3823,3842,3849,3871,3890,3904,3923,3939,3951,3970,3980,3987,3998,4013,4027,4046,4053,4077,4088,4095,4105,4121,4133,4155,4179,4189,4213,4228,4238,4245,4259,4278,4292,4309,4320,4339,4358,4372,4391,4398,4422,4438,4457,4471,4481,4500,4516,4531,4555,4577,4583,4600,4624,4630,4651,4675,4690,4701,4723,4742,4756,4780,4799,4818,4835,4847,4858,4877,4884,4891,4915,4929,4953,4967,4978,4997,5016,5026,5040,5059,5078,5084,5094,5105,5124,5138,5157,5172,5187,5201,5220,5231,5245,5264,5279,5300,5310,5327,5351,5361,5383,5397,5416,5435,5449 ]

# define the thermodynamic/alchemical state
# the system is prepared at the alchemical intermediate state at lambda=0
temperature = 298.15 * kelvin
lmbd = 0.00
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
posrestr_atoms = []
last_not_solvent_atom = lig2_atoms[len(lig2_atoms)-1]
for at in prmtop.topology.atoms():
    if at.index <= last_not_solvent_atom:
        posrestr_atoms.append(at.index)
atm_utils.addPosRestraints(posrestr_atoms, inpcrd.positions, fc, tol)

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

print( "LoadState equil.xml")
simulation.loadState('equil.xml')

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

totalSteps = 250000
steps_per_cycle = 1000
number_of_cycles = int(totalSteps/steps_per_cycle)
deltalambda = (0.5 - 0.0)/float(number_of_cycles)

simulation.reporters.append(
    StateDataReporter(stdout, steps_per_cycle,
     step=True, potentialEnergy = True, temperature=True, volume=True))
# binding energy values and other parameters are recorded in this file
with open("mdlambda.out", 'w') as fh:
    for i in range(number_of_cycles):
        simulation.step(steps_per_cycle)
        state = simulation.context.getState(getEnergy = True, groups = {0,atmforcegroup})
        pot_energy = (state.getPotentialEnergy()).value_in_unit(kilocalorie_per_mole)
        pert_energy = (atmforce.getPerturbationEnergy(simulation.context)).value_in_unit(kilocalorie_per_mole)
        l1 = simulation.context.getParameter(atmforce.Lambda1())
        l2 = simulation.context.getParameter(atmforce.Lambda2())
        a = simulation.context.getParameter(atmforce.Alpha()) / kilojoules_per_mole
        umid = simulation.context.getParameter(atmforce.U0()) * kilojoules_per_mole
        w0 = simulation.context.getParameter(atmforce.W0()) * kilojoules_per_mole
        print("%f %f %f %f %f %f %f %f %f" % (temperature/kelvin,lmbd, l1, l2, a*kilocalorie_per_mole, umid/kilocalorie_per_mole, w0/kilocalorie_per_mole, pot_energy, pert_energy), file=fh )
        fh.flush()
        lmbd += deltalambda
        lambda1 += deltalambda
        lambda2 += deltalambda
        simulation.context.setParameter(atmforce.Lambda1(), lambda1)
        simulation.context.setParameter(atmforce.Lambda2(), lambda2)


print( "SaveState: mdlambda.xml.")
simulation.saveState('mdlambda.xml')

positions = simulation.context.getState(getPositions=True).getPositions()
boxsize = simulation.context.getState().getPeriodicBoxVectors()
simulation.topology.setPeriodicBoxVectors(boxsize)
with open('mdlambda.pdb', 'w') as output:
  PDBFile.writeFile(simulation.topology, positions, output)

end=datetime.now()
elapsed=end - start
print("elapsed time for equilbration ="+str(
    elapsed.seconds+elapsed.microseconds*1e-6)+"s")
