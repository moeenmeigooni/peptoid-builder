from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout, exit, stderr, argv
import numpy as np
import argparse
from pathlib import Path

def parse_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-in', '--input_path', default='long_bonds.csv', help='Path to csv file containing unresolved ring piercings')
    args = parser.parse_args()
    return args

def main(args):
    # parse input csv, return if does not exist
    input_file = Path(args.input_path)
    if not input_file.is_file():
        return
    
    pdbs, tops, terminal_ids = [], [], []
    with open(input_file.resolve(), 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line is not '':
            spl = line.split(',')
            pdbi = spl[0]
            topi = spl[1]
            idsi = spl[2][:-1] # remove endl
            pdbs.append(pdbi)
            tops.append(topi)
            terminal_ids.append([int(ids) for ids in idsi.split(' ')])
    assert len(pdbs) == len(tops) == len(terminal_ids), 'Big oops'
    n_sims = len(pdbs)
    
    for sim in range(n_sims):
        # load structure and topology
        pdb = PDBFile(pdbs[sim])
        prmtop = AmberPrmtopFile(tops[sim])

        # create system, integrator
        system = prmtop.createSystem(nonbondedMethod=NoCutoff)
        integrator = LangevinMiddleIntegrator(298*kelvin, 1/picosecond, 1*femtoseconds)

        # set vdw eps to zero for some ring atoms and save old params for later
        nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
        unmodified_epsilons = []
        unmodified_sigmas = []
        unmodified_charges = []
        for ind in terminal_ids[sim]:
            charge, sigma, epsilon = nonbonded.getParticleParameters(ind)
            unmodified_epsilons.append(epsilon)
            unmodified_sigmas.append(sigma)
            unmodified_charges.append(charge)
            nonbonded.setParticleParameters(ind, 0, sigma, 0)

        # set up platform for GPU
        platform = Platform.getPlatformByName('CUDA')
        properties = {'DeviceIndex': '0', 'Precision': 'mixed'}
        simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
        # set starting coordinates
        simulation.context.setPositions(pdb.positions)

        # do minimization while zeroing out NB forces on terminal ring atoms
        simulation.minimizeEnergy()
        
        # restore NB params of terminal ring atoms
        for i,ind in enumerate(terminal_ids[sim]):
            nonbonded.setParticleParameters(ind, unmodified_charges[i], unmodified_sigmas[i], unmodified_epsilons[i])
        nonbonded.updateParametersInContext(simulation.context)

        # re-minimize after restoring NB forces
        simulation.minimizeEnergy()

        # save minimized structure
        state = simulation.context.getState(getPositions=True, getEnergy=True)
        outname = pdbs[sim]# resave as same name #.split('.')[0] + '_unpierced.pdb'
        with open(outname, 'w') as outfile:
            PDBFile.writeFile(simulation.topology, state.getPositions(), outfile)

if __name__ == '__main__':
    args = parse_arguments()
    main(args)
