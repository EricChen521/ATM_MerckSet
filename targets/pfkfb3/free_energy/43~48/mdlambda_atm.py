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

jobname = '43~48'


displ = [51.6,12.84,-4.99]
displacement      = [  displ[i] for i in range(3) ] * angstrom
lig1_restr_offset = [  0.       for i in range(3) ] * angstrom
lig2_restr_offset = [  displ[i] for i in range(3) ] * angstrom

lig1_atoms = [6749,6750,6751,6752,6753,6754,6755,6756,6757,6758,6759,6760,6761,6762,6763,6764,6765,6766,6767,6768,6769,6770,6771,6772,6773,6774,6775,6776,6777,6778,6779,6780,6781,6782,6783,6784,6785,6786,6787,6788,6789,6790,6791]
lig2_atoms = [6792,6793,6794,6795,6796,6797,6798,6799,6800,6801,6802,6803,6804,6805,6806,6807,6808,6809,6810,6811,6812,6813,6814,6815,6816,6817,6818,6819,6820,6821,6822,6823,6824,6825,6826,6827,6828,6829,6830,6831,6832,6833,6834,6835,6836,6837,6838]
refatoms_lig1 = [22,4,10]
refatoms_lig2 = [22,4,10]
rcpt_cm_atoms = [196,220,1982,1993,2009,2125,2175,2190,2970,2986,3344,3363,3394,3413,3430,6369]
restrained_atoms = [ 9,23,42,48,62,78,97,113,130,146,153,180,186,196,220,227,249,263,284,303,314,336,358,377,391,415,436,455,469,493,512,519,543,549,563,585,601,621,635,651,658,673,694,718,742,757,767,783,805,822,843,854,865,886,900,920,940,972,978,990,1004,1019,1034,1044,1061,1083,1099,1123,1145,1162,1173,1183,1202,1212,1222,1241,1265,1277,1293,1315,1326,1347,1366,1376,1398,1413,1420,1427,1444,1463,1473,1489,1509,1521,1531,1545,1559,1573,1587,1611,1626,1650,1674,1691,1708,1727,1746,1763,1783,1793,1815,1830,1844,1856,1876,1898,1908,1928,1948,1967,1982,1993,2009,2020,2032,2052,2058,2072,2088,2104,2114,2125,2139,2158,2175,2190,2206,2228,2247,2258,2277,2283,2295,2316,2338,2350,2361,2375,2386,2396,2411,2421,2438,2450,2462,2482,2499,2521,2545,2564,2575,2586,2607,2622,2632,2643,2664,2689,2695,2714,2734,2740,2752,2774,2785,2797,2821,2833,2852,2863,2882,2901,2923,2939,2958,2970,2986,2993,3017,3041,3061,3080,3096,3110,3134,3150,3167,3179,3196,3215,3232,3243,3267,3286,3302,3323,3344,3363,3380,3394,3413,3430,3446,3471,3477,3501,3515,3534,3555,3574,3585,3609,3626,3633,3648,3662,3677,3694,3708,3727,3744,3751,3775,3794,3801,3808,3820,3831,3838,3857,3868,3879,3903,3910,3932,3954,3974,3984,3995,4005,4024,4035,4057,4077,4093,4108,4123,4140,4154,4173,4195,4207,4226,4250,4266,4290,4304,4315,4332,4351,4373,4384,4398,4417,4434,4448,4458,4473,4483,4502,4526,4553,4559,4580,4595,4612,4636,4658,4668,4687,4701,4716,4735,4747,4757,4764,4780,4791,4806,4821,4840,4854,4875,4890,4905,4924,4948,4960,4974,5003,5009,5024,5039,5060,5070,5089,5113,5128,5145,5157,5179,5200,5221,5242,5266,5295,5301,5315,5322,5337,5348,5369,5386,5398,5417,5433,5450,5474,5493,5516,5522,5538,5557,5574,5589,5608,5623,5647,5664,5679,5693,5709,5728,5744,5763,5774,5791,5808,5818,5834,5853,5877,5888,5907,5926,5936,5957,5977,5996,6008,6030,6041,6051,6066,6081,6106,6112,6133,6152,6174,6193,6199,6218,6235,6249,6265,6284,6306,6325,6347,6353,6369,6379,6400,6407,6418,6442,6458,6473,6484,6503,6524,6543,6557,6573,6588,6599,6615,6626,6640,6657,6681,6696,6720,6731 ]

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
