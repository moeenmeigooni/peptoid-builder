from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout, exit, stderr, argv
import numpy as np
_, rep = argv

# file names
system_basename = 'peptoid'
crd_file = f'{system_basename}_rep{rep}.crd'
prm_file = f'{system_basename}_rep{rep}.parm7'

# load structure and topology
inpcrd = AmberInpcrdFile(crd_file)
prmtop = AmberPrmtopFile(prm_file)

# create system, integrator
system = prmtop.createSystem(nonbondedMethod=NoCutoff)
integrator = LangevinMiddleIntegrator(298*kelvin, 1/picosecond, 1*femtoseconds)

# set up platform for GPU
platform = Platform.getPlatformByName('CUDA')
properties = {'DeviceIndex': '0', 'Precision': 'mixed'}
simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
# set starting coordinates
simulation.context.setPositions(inpcrd.positions)

# do quick minimization
simulation.minimizeEnergy()

# save minimized structure
state = simulation.context.getState(getPositions=True, getEnergy=True)
with open(f'{system_basename}_rep{rep}_min.pdb', 'w') as outfile:
    PDBFile.writeFile(simulation.topology, state.getPositions(), outfile)
