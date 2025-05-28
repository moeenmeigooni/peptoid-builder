import argparse
import numpy as np
import os

def random_value(a, b, n=1):
    # returns n random values in the interval [a,b)
    rng = np.random.default_rng()
    return (b - a) * rng.random(n) + a

def random_cistrans_angle(n=1):
    # returns random choice of 0. or 180. with equal probability
    rng = np.random.default_rng()
    return rng.integers(0,2,n) * 180

def write_file(output_name, lines):
    with open(output_name, 'w') as f:
        for line in lines:
            f.write(line + '\n')
class PeptoidLeaper:
    """
    Writes input files to generate and solvate peptoid structures in tleap
    """
    def __init__(self,
                 sequence,
                 replicate_ids,
                 output_directory=None,
                 system_basename='peptoid',
                 stepsff_pathname='/Scr/meigoon2-new/peptoid/STEPs',
                 solvate_input_suffix='min_align'):

        # check output dir
        if output_directory is not None:
            assert type(output_directory) is str, 'output_directory must be str type'
            self.output_directory = output_directory
        else:
            self.output_directory = os.getcwd()

        # check sequence
        if sequence is not None:
            assert type(sequence) is list, 'sequence must be list type'
            for res in sequence:
                assert type(res) is str, 'sequence elements must be int type'
            self.sequence = sequence

        # check replicate ids
        if replicate_ids is not None:
            assert type(replicate_ids) is list, 'replicate_ids must be list type'
            for rep in replicate_ids:
                assert type(rep) is int, 'replicate_ids elements must be int type'
            self.replicate_ids = replicate_ids
        else:
            raise ValueError('`replicate_ids` input not valid')

        # other class attributes
        self.system_basename = system_basename   # base name prefix for output structure/topology files
        self.stepsff_pathname = stepsff_pathname # path to STEPs forcefield directory
        self.solvate_input_suffix = solvate_input_suffix # input file suffix after minimizing, aligning

    def write_generate_inputs(self):
        omega_angles = random_cistrans_angle(n=len(self.replicate_ids)) # set omega to 0 or 180
        phi_angles = random_value(-180,180,n=len(self.replicate_ids))   # set phi to random [-180, 180)
        psi_angles = random_value(-180,180,n=len(self.replicate_ids))   # set psi to random [-180, 180)
        for replicate_id in self.replicate_ids:
            # initialize, load STEPs forcefield
            lines = []
            lines.append(f'logfile build_rep{replicate_id}.log')
            lines.append(f'source {self.stepsff_pathname}/load_peptoids.leap')
            # build sequence
            lines.append(f'{self.system_basename} = sequence {{ {" ".join(self.sequence)} }}')
            # determine number of dihedrals to rotate based on sequence length
            # assume uncapped N-terminus and NH2 cap on C-terminus
            n_res = len(self.sequence) - 1
            # randomize dihedrals
            omega_angles = random_cistrans_angle(n=n_res-1) # set omega to 0 or 180
            phi_angles = random_value(-180,180,n=n_res-1)   # set phi to random [-180, 180)
            psi_angles = random_value(-180,180,n=n_res-1)   # set psi to random [-180, 180)
            for i in range(n_res - 1):
                rand_omega = omega_angles[i]
                rand_phi = phi_angles[i]
                rand_psi = psi_angles[i]
                lines.append(f'impose {self.system_basename} {{ {i+1} }} {{{{ "O" "C" "N" "C1" {rand_omega:.2f} }}}}')
                lines.append(f'impose {self.system_basename} {{ {i+1} }} {{{{ "C" "N" "CA" "C" {rand_phi:.2f} }}}}')
                lines.append(f'impose {self.system_basename} {{ {i+1} }} {{{{ "N" "CA" "C" "N" {rand_psi:.2f} }}}}')
            # save parm7, inpcrd, pdb
            lines.append(f'saveamberparm {self.system_basename} {self.system_basename}_rep{replicate_id}.parm7 {self.system_basename}_rep{replicate_id}.crd')
            lines.append(f'savepdb {self.system_basename} {self.system_basename}_rep{replicate_id}.pdb')
            # quit
            lines.append('q')
            # save leap file
            write_file(f'{self.output_directory}/build_rep{replicate_id}.leap', lines)

    def write_solvate_inputs(self):
        for replicate_id in self.replicate_ids:
            # load precomputed padding
            with open(f'pad_rep{replicate_id}.dat', 'r') as f:
                dat = f.readlines()
            padx,pady,padz = [f'{float(x):.4f}' for x in dat[0][:-1].split()]

            # initialize, load forcefields
            lines = []
            lines.append(f'logfile solvate_rep{replicate_id}.log')
            lines.append(f'source {self.stepsff_pathname}/load_peptoids.leap')
            lines.append('source leaprc.water.tip3p')
            # remap names of hydrogens in NH2 C-term cap because VMD renamed it :(
            lines.append('addPdbAtomMap {{HN1 H1} {HN2 H2}}')
            # load aligned pdb
            lines.append(f'{self.system_basename} = loadpdb {self.system_basename}_rep{replicate_id}_{self.solvate_input_suffix}.pdb')
            # solvate and save
            lines.append(f'solvateBox {self.system_basename} TIP3PBOX {{ {padx}, {pady}, {padz} }}')
            lines.append(f'saveamberparm {self.system_basename} {self.system_basename}_water_rep{replicate_id}.parm7 {self.system_basename}_water_rep{replicate_id}.crd')
            lines.append(f'savepdb {self.system_basename} {self.system_basename}_water_rep{replicate_id}.pdb')
            lines.append('quit')
            write_file(f'{self.output_directory}/solvate_rep{replicate_id}.leap', lines)

def parse_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-seq', '--sequence', help='Comma-separated residue names (e.g. "R59,R20,R20,R59,NH2"')
    parser.add_argument('-reps', '--replicate_ids', help='Comma-separated replicate indices (e.g. "0,1,2,3", "8,9")')
    parser.add_argument('-mode', '--run_mode', help='Choice of "generate" or "solvate"', choices=['generate','solvate'])
    args = parser.parse_args()
    return args

def main(args):
    sequence = args.sequence.split(',')
    replicate_ids = [int(ind) for ind in args.replicate_ids.split(',')]
    builder = PeptoidLeaper(sequence, replicate_ids)
    if args.run_mode == 'generate':
        builder.write_generate_inputs()
    elif args.run_mode == 'solvate':
        builder.write_solvate_inputs()

if __name__ == '__main__':
    args = parse_arguments()
    main(args)
