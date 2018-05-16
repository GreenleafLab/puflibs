#!/usr/bin/env python
""" Fit on or off rates.

Sarah Denny """

##### IMPORT #####
import numpy as np
import pandas as pd
import sys
import os
import argparse
import itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from puflibs import seqmodel
from Bio import SeqIO
################ Parse input parameters ################

#set up command line argument parser
parser = argparse.ArgumentParser(description='bootstrap on off rate fits')
parser.add_argument('-s', '--sequence', metavar='ACGU', 
                   help='sequence to test')
parser.add_argument('-fi', '--infasta',
                    help='input fasta files containing sequences to test')
parser.add_argument('-o', '--output', default='test.dat',
                    help='output filename. Only required if input is fasta')

parser.add_argument('-m', '--model' , default='annotations/RNAmap/qMotif_20180302_',
                   help='basename of files giving the model parameters')
parser.add_argument('-t', '--temperature', type=float, default=37,
                   help='temperature in Celsius')
parser.add_argument('-p', '--plot', action='store_true', 
                   help='if flagged, make plot of occupancy')

def check_sequence(sequence):
    """Make sure sequence is all uppercase and has only A,C,G,U"""
    # check sequence
    acceptable_bases = ['A', 'C', 'G', 'U']
    if any([s=='T' for s in sequence]):
        sys.stderr.write("replacing T's in sequence with U's")
        sequence = sequence.replace('T', 'U')
    if not all([s in acceptable_bases for s in sequence]):
        if not all([s in acceptable_bases for s in sequence.upper()]):
            raise ValueError('Sequence improperly formatted')
        else:
            sys.stderr.write("replacing lowercase letters with uppercase")
            sequence = sequence.upper()
    return sequence

def find_occupancy(sequence, base_params, coupling_params, flip_params, dflip_params, temperature):
    """Find dG and then occupancy across sequence"""
    ddGs = seqmodel.find_energy_for_1nt_sequence_registers(sequence, base_params, coupling_params, flip_params, dflip_params, temperature)
    # return occupancy
    occupancy = np.exp(ddGs/seqmodel.get_ddG_conversion(temperature))
    return occupancy

def make_plot(occupancy, kind=None, sequence=None):
    """Plot the occupancy across the sequence"""
    plt.figure(figsize=(10,3));
    xvalues = occupancy.index.tolist()
    length_cutoff = 100
    if kind is None:
        # if the len of the region of interest is very long, don't plot full binding region, just value at center of register
        if len(occupancy) > length_cutoff:
            kind='plot'
        else:
            kind='bar'
        
    if kind=='plot':
        plt.plot(xvalues, occupancy, color='0.5', linewidth=1)
    elif kind=='bar':
        plt.bar(xvalues, occupancy, color='0.5', linewidth=1, width=0.8, alpha=0.4)
        
    plt.xlim(xvalues[0]-0.5, xvalues[-1]+0.5)
    plt.ylim(0, 1.5)
    plt.xlabel('location of center of 11nt register (nt)')
    plt.ylabel('occupancy')
    plt.tight_layout()
    
    # annotate the sequence at the top of the plot
    if sequence is not None and len(sequence)==len(occupancy) and len(sequence) < length_cutoff:
        for x, s in zip(xvalues, sequence):
            plt.annotate(s, (x, 1.5), horizontalalignment='center', verticalalignment='top', )

    

if __name__ == '__main__':
    args = parser.parse_args()
    
    # load model parameters
    temperature = args.temperature
    flip_params, base_params, coupling_params, dflip_params = seqmodel.load_params(args.model)
    
    # load sequences
    if args.sequence is not None:
        # proceed with just a single sequence
        sequence = check_sequence(args.sequence)
        occupancy = find_occupancy(sequence, base_params, coupling_params, flip_params, dflip_params, temperature)
        output = '\n'.join(['%.3e'%occupancy.loc[i] if i in occupancy.keys() else 'nan' for i in range(len(sequence))])
        sys.stdout.write(output + '\n')
        if args.plot:
            make_plot(occupancy.loc[range(len(sequence))], sequence=sequence)
            plt.savefig(os.path.splitext(args.output)[0] + '.pdf')
            
    elif args.infasta is None:
        raise InputError('Need to provide sequence or input fasta file')
    else:
        # proceed with fasta file
        fasta_sequences = SeqIO.parse(open(args.infasta),'fasta')
        with open(args.output, 'w') as f:
            for fasta in fasta_sequences:
                sequence = check_sequence(fasta.seq.tostring())
                occupancy = find_occupancy(sequence, base_params, coupling_params, flip_params, dflip_params, temperature)
                f.write('>%s\n'%fasta.id)
                for i in range(len(sequence)):
                    if i in occupancy.keys():
                        f.write('%.3e\n'%occupancy.loc[i])
                    else:
                        f.write('nan\n')
                if args.plot:
                    make_plot(occupancy, sequence=sequence)
                    plt.savefig(os.path.splitext(args.output)[0] + '.pdf')
    
    
