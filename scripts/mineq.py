from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout, exit, stderr, argv
import numpy as np
import mdtraj as md
_, rep, stage = argv

# stage dict
stagemap = {'6': 0.6,
            '9': 0.9,
            '12': 1.2,
            '15': 1.5}

# file names
system_basename = 'peptoid'
pdb_file = f'{system_basename}_water_rep{rep}.pdb'
prm_file = f'{system_basename}_water_rep{rep}.parm7'
state_freq = 1000

# load structure and topology
pdb = PDBFile(pdb_file)
prmtop = AmberPrmtopFile(prm_file)

# get relevant atom indices for custom forces
t = md.load(pdb_file) 
unitcell = [Vec3(*vec) for vec in t.unitcell_vectors[0]]
#sequence = [residue.name for residue in t.top.residues if residue.chain.index == 0]
sequence = [residue.name for residue in t.top.residues if not residue.is_water]
n_residues = len(sequence)
n_atoms = t.n_atoms
protein_ind      = [int(ind) for ind in t.top.select('not water')] # indices must be int, not np.int64
s_ind1, s_ind2   = [int(ind) for ind in t.top.select('name S1 and resid 0 %i'%(n_residues - 2))] # last res is NH2 cap
cg_ind1, cg_ind2 = [int(ind) for ind in t.top.select('name C2 and resid 0 %i'%(n_residues - 2))] # first and second-to-last are N-methionine
ce_ind1, ce_ind2 = [int(ind) for ind in t.top.select('name C3 and resid 0 %i'%(n_residues - 2))] # so get anchor indices from resid 0 and n-2

# create system, integrator
system = prmtop.createSystem(nonbondedMethod=PME,
                             nonbondedCutoff=1.2*nanometer, switchDistance=1.0*nanometer,
                             constraints=HBonds, rigidWater=True)
integrator = LangevinMiddleIntegrator(298*kelvin, 1/picosecond, 2*femtoseconds)

# custom forces
# add S-S force via CustomCompoundBondForce
sspull_constant    = 418.4 # units: kJ/(mol*nm^2)
sspull_target_dist = stagemap[stage]   # units: nm
sspull = CustomCompoundBondForce(2, '0.5*k_ssp*((z2 - z1) - targetdist)^2')
sspull.addGlobalParameter('k_ssp', sspull_constant)
sspull.addGlobalParameter('targetdist', sspull_target_dist)
sspull.addBond([s_ind1, s_ind2], [])
system.addForce(sspull)

# add force orienting sulfur lone pair to z axis 
ssorient_constant = 4.184 * 10. # 10 kcal/mol when unaligned, -10 kcal/mol when aligned
ssorient = CustomCompoundBondForce(3, 'dir*k_sso*(z1 - (z2+z3)/2)/pointdistance(x1, y1, z1, (x2+x3)/2, (y2+y3)/2, (z2+z3)/2)')
ssorient.addPerBondParameter('dir')
ssorient.addGlobalParameter('k_sso', ssorient_constant)
ssorient.addBond([s_ind1, cg_ind1, ce_ind1], [1.])
ssorient.addBond([s_ind2, cg_ind2, ce_ind2], [-1.])
system.addForce(ssorient)

# add electric field
applied_potential   = 0.1       # Volts
separation_distance = sspull_target_dist  # nm
au_s_bond_length    = 0.25      # nm
avogadro_number     = 6.022e23
electron_volt       = 1.602e-19 # C/e-
J_to_kJ             = 0.001
efield_strength     = applied_potential * electron_volt * avogadro_number * J_to_kJ / (separation_distance + 2. * au_s_bond_length)
efield = CustomExternalForce('q*Ez*z')
efield.addPerParticleParameter('q')
efield.addGlobalParameter('Ez', efield_strength) # kJ/(mol*nm*e-)
nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
for i in range(n_atoms):
    charge, sigma, epsilon = nonbonded.getParticleParameters(i)
    efield.addParticle(i, [charge])
system.addForce(efield)

# add PBC and barostat
system.setDefaultPeriodicBoxVectors(*unitcell)
system.addForce(MonteCarloBarostat(1*atmosphere, 298*kelvin))

# set up platform for GPU
platform = Platform.getPlatformByName('CUDA')
properties = {'DeviceIndex': '0', 'Precision': 'mixed'}
simulation = Simulation(prmtop.topology, system, integrator, platform, properties)

# set starting coordinates
simulation.context.setPositions(pdb.positions)

# create log reporter
simulation.reporters.append(StateDataReporter(stdout, state_freq, step=True, potentialEnergy=True,
                                              temperature=True, speed=True, elapsedTime=True))

# do quick minimization
simulation.minimizeEnergy()

# equilibrate for 1 ns
simulation.step(500000)

# save minimized structure
state = simulation.context.getState(getPositions=True, getEnergy=True)
with open(f'{system_basename}_water_rep{rep}_hold{stage}_mineq.pdb', 'w') as outfile:
    PDBFile.writeFile(simulation.topology, state.getPositions(), outfile)

