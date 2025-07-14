import numpy as np
import argparse
from pathlib import Path
import mdtraj as md

class LongBondChecker:
    """
    Measures bond distances to ensure all bonds are below a maximum length.
    Long bonds are indicators for ring piercings. 
    Must use minimized structures for meaningful results.

    Parameters
    ----------
    pdb : str
          Path to minimized structure in pdb format.
    top : str
          Path to topology file corresponding to `pdbfile`.
    output : str
             Path for output csv file
    distance_cutoff : float
                      Distance above which a bond will be classified as long
    verbose : bool
              Unused
    """
    def __init__(self, pdb, top, output, distance_cutoff, verbose):
        pdb_file = Path(pdb)
        assert pdb_file.exists(), 'Invalid path to pdb file'
        self.pdb_path = pdb_file.resolve()
        
        top_file = Path(top)
        assert top_file.exists(), 'Invalid path to top file'
        self.top_path = top_file.resolve()

        self.traj = md.load(pdb, top)
        assert self.traj.n_frames == 1, 'Input pdb must contain only one frame'
        
        assert type(output) is str, 'output path must be string'
        self.output = output

        assert type(distance_cutoff) is float, 'Distance cutoff must be float'
        self.distance_cutoff = distance_cutoff # maximum allowable bond length (nm)
    
    def check_for_long_bonds(self):
        # get list of all bonded atom pairs
        pair_list = []
        for bond in self.traj.top.bonds:
            pair_list.append([bond.atom1.index, bond.atom2.index])
        
        # calculate distance
        self.bond_distances = md.compute_distances(self.traj, pair_list)[0] # assume 1 frame traj

        # iterate over bonds, check if any are above cutoff:
        self.n_long_bonds = 0
        self.long_bond_pairs = []
        for i,bond_distance in enumerate(self.bond_distances):
            if bond_distance > self.distance_cutoff:
                self.n_long_bonds += 1
                self.long_bond_pairs.append(pair_list[i])

    def identify_terminal_ring_atoms(self):
        """
        Given list of atom indices participating in long bonds, identify which rings
        are pierced and atom indices of their terminal atoms.
        This function should only be called if `pair_list` is not an empty list.
        """
        # define terminal and ring atoms for each residue of interest
        terminal_atoms_by_res = {'R20': ['C4', 'C5', 'C6', 'H4', 'H5', 'H6'],
                                 'R02': ['C5', 'C6', 'C7', 'H6', 'H7', 'H8']}
        ring_atoms_by_res = {'R20': ['C2', 'C3', 'C4', 'C5', 'C6', 'C7'],
                             'R02': ['C3', 'C4', 'C5', 'C6', 'C7', 'C8']}
        
        # get residues indices and names
        resids = [residue.index for residue in self.traj.top.residues]
        resnames = [residue.name for residue in self.traj.top.residues]

        # get residues with rings
        ring_resids = [resids[i] for i,resname in enumerate(resnames) if resname in ring_atoms_by_res.keys()]
        
        # get atom indices of each ring
        ring_ids = []
        for i,resid in enumerate(ring_resids):
            ring_ids.append(self.traj.top.select(f'resid {resid} and name {" ".join(ring_atoms_by_res[resnames[resid]])}'))
        
        # calculate COM of each ring
        ring_coms = [np.mean(self.traj.xyz[0,ind],axis=0) for ind in ring_ids]

        # calculate midpoints of long bonds
        long_bond_midpoints = [np.mean(self.traj.xyz[0,pair], axis=0) for pair in self.long_bond_pairs]

        # calculate which rings are being pierced by which rings
        pierced_ring_resids = [] # relative to ring_resids
        for midpoint in long_bond_midpoints:
            dists = np.sqrt(np.sum(np.square(ring_coms - midpoint),axis=1)) 
            if dists.min() < 0.1: # 1 Angstrom cutoff
                pierced_ring_resids.append(dists.argmin())
        
        # get terminal atoms of pierced rings
        self.terminal_atoms = []
        for pierced_ring_resid in pierced_ring_resids:
            real_resid = ring_resids[pierced_ring_resid]
            self.terminal_atoms.append(self.traj.top.select(f'resid {real_resid} and name {" ".join(terminal_atoms_by_res[resnames[real_resid]])}'))


def parse_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-pdb', '--pdb_file', help='Path to pdb file containing one frame')
    parser.add_argument('-top', '--top_file', help='Path to topology file (e.g. psf) corresponding to input pdb file')
    parser.add_argument('-out', '--output', default='long_bonds.csv', help='Path to output file to save atom indices of terminal ring atoms')
    parser.add_argument('-cutoff', '--distance_cutoff', default=0.25, help='Cutoff distance for defining long bonds')
    parser.add_argument('-v', '--verbose', default=True, help='verbose=True prints some output to stdout')
    args = parser.parse_args()
    return args

def main(args):
    checker = LongBondChecker(args.pdb_file, args.top_file, args.output, args.distance_cutoff, args.verbose)
    checker.check_for_long_bonds()
    print(f'{checker.n_long_bonds} long bond(s) found')

    # check for long bonds, proceed only if they exist
    if checker.n_long_bonds > 0:
        for pair in checker.long_bond_pairs:
            print(f'\t long bond between atoms: {*pair,}')
        
        # write output csv only if long bonds are matched to pierced rings
        checker.identify_terminal_ring_atoms()
        if len(checker.terminal_atoms) > 0:
            for terminal in checker.terminal_atoms:
                print(f'\t terminal atom indices: {*terminal,}')
        
            # remove duplicates from terminal atom indices, if any
            all_terminal_ids = list(set([idx for idy in checker.terminal_atoms for idx in idy]))
            with open(args.output, 'a') as f:
                f.write(f'{checker.pdb_path},{checker.top_path},{" ".join([str(ind) for ind in all_terminal_ids])}\n')

if __name__ == '__main__':
    args = parse_arguments()
    main(args)
