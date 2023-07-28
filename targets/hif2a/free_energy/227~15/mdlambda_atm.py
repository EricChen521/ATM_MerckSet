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

jobname = '227~15'


displ = [43.06,0.54,-0.88]
displacement      = [  displ[i] for i in range(3) ] * angstrom
lig1_restr_offset = [  0.       for i in range(3) ] * angstrom
lig2_restr_offset = [  displ[i] for i in range(3) ] * angstrom

lig1_atoms = [1758,1759,1760,1761,1762,1763,1764,1765,1766,1767,1768,1769,1770,1771,1772,1773,1774,1775,1776,1777,1778,1779,1780,1781,1782,1783,1784,1785,1786,1787,1788,1789,1790,1791,1792,1793,1794,1795]
lig2_atoms = [1796,1797,1798,1799,1800,1801,1802,1803,1804,1805,1806,1807,1808,1809,1810,1811,1812,1813,1814,1815,1816,1817,1818,1819,1820,1821,1822,1823,1824,1825,1826,1827,1828,1829,1830,1831,1832,1833]
refatoms_lig1 = [6,9,7]
refatoms_lig2 = [6,9,7]
rcpt_cm_atoms = [209,617,663,683,802,855,866,1036,1116,1310,1341]
restrained_atoms = [ 9,28,40,51,73,87,107,126,137,152,169,180,197,209,226,248,268,282,303,314,326,338,362,381,395,410,429,448,455,476,501,507,522,537,556,575,582,606,617,627,648,663,683,704,721,731,750,762,773,788,802,819,833,855,866,883,900,914,933,944,958,980,987,1004,1020,1036,1047,1054,1071,1092,1116,1133,1152,1162,1184,1201,1208,1215,1236,1252,1276,1295,1310,1324,1341,1348,1362,1378,1397,1418,1440,1446,1470,1484,1503,1528,1534,1551,1562,1581,1598,1609,1625,1639,1660,1676,1695,1706,1721,1740 ]

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
