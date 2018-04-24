import luigi
import sciluigi
import os
import subprocess
import logging
import itertools
import pandas as pd
import numpy as np
import ipdb
import matplotlib.pyplot as plt
import seaborn as sns
from puflibs import variables, processing


class FilenameToTaskOutput(sciluigi.Task):
    """Take an input and make that the output"""
    # parameters
    filename = luigi.Parameter()
    
    # outputs
    def out_file(self):
        return sciluigi.TargetInfo(self, self.filename)
    
    # run
    def run(self):
        pass
    
class CombineBeds(sciluigi.Task):
    """Filter the bed output from HOMER for strandedness, protein-coding, etc."""
    # parameters
    outdir = luigi.Parameter()
    
    # inputs
    in_beds = None

    
    # outputs
    def out_bed(self):

        basenames = processing.combine_filenames([target().path for target in self.in_beds])
        basename = '.'.join(basenames)
        
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, basename + '.bed'))    
    # run
    def run(self):
        # make out directory if it doesn't exist
        dirname = os.path.dirname(self.out_bed().path)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
            
        combined_data = pd.concat([processing.load_bed(target().path) for target in self.in_beds])
        processing.save_bed(combined_data.sort_values(['chrm', 'start']), self.out_bed().path)


      
class DivideBed2(sciluigi.Task):
    """Divide Bed into two beds."""
    # parameters
    outdir = luigi.Parameter()
    num = 2.
    
    # inputs
    in_bed = None

    ## outputs
    def out_bed0(self):
        return sciluigi.TargetInfo(self, self.get_outval(0))
    
    def out_bed1(self):
        return sciluigi.TargetInfo(self, self.get_outval(1))

    def get_prefix(self):
        outdir = self.outdir
        basename = os.path.splitext(os.path.basename(self.in_bed().path))[0]
        return os.path.join(outdir, basename)

    def get_outval(self, i):
        prefix = self.get_prefix()
        return '%s%02d.bed'%(prefix, i)  
    
    # run
    def run(self):
        # make out directory if it doesn't exist
        dirname = os.path.dirname(self.get_prefix())
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        
        # divide into n lines
        num_total_lines = int(subprocess.check_output("wc -l %s | awk '{print $1}'"%self.in_bed().path, shell=True).strip()) - 1
        num_lines_per_file = int(np.ceil(num_total_lines/self.num))

        task_call = ('awk \'{if (substr($1, 1, 1)!="#") print}\' %s | split -d -l %d --additional-suffix .bed - %s'
                     %(self.in_bed().path, num_lines_per_file, self.get_prefix()))
        subprocess.call(task_call, shell=True) 
        
        
class DivideBed10(sciluigi.Task):
    """Filter the bed output from HOMER for strandedness, protein-coding, etc."""
    # parameters
    outdir = luigi.Parameter()
    num = 10.
    
    # inputs
    in_bed = None

    ## outputs
    def out_bed0(self):
        return sciluigi.TargetInfo(self, self.get_outval(0))
    def out_bed1(self):
        return sciluigi.TargetInfo(self, self.get_outval(1))
    def out_bed2(self):
        return sciluigi.TargetInfo(self, self.get_outval(2))
    def out_bed3(self):
        return sciluigi.TargetInfo(self, self.get_outval(3))
    def out_bed4(self):
        return sciluigi.TargetInfo(self, self.get_outval(4))
    def out_bed5(self):
        return sciluigi.TargetInfo(self, self.get_outval(5))
    def out_bed6(self):
        return sciluigi.TargetInfo(self, self.get_outval(6))
    def out_bed7(self):
        return sciluigi.TargetInfo(self, self.get_outval(7))
    def out_bed8(self):
        return sciluigi.TargetInfo(self, self.get_outval(8))
    def out_bed9(self):
        return sciluigi.TargetInfo(self, self.get_outval(9))

    def get_prefix(self):
        outdir = self.outdir
        basename = os.path.splitext(os.path.basename(self.in_bed().path))[0]
        return os.path.join(outdir, basename)

    def get_outval(self, i):
        prefix = self.get_prefix()
        return '%s.%02d.bed'%(prefix, i)  
    
    # run
    def run(self):
        # make out directory if it doesn't exist
        dirname = os.path.dirname(self.get_prefix())
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        
        # divide into n lines
        num_total_lines = int(subprocess.check_output("wc -l %s | awk '{print $1}'"%self.in_bed().path, shell=True).strip()) - 1
        num_lines_per_file = int(np.ceil(num_total_lines/self.num))

        task_call = ('awk \'{if (substr($1, 1, 1)!="#") print}\' %s | split -d -l %d --additional-suffix .bed - %s.'
                     %(self.in_bed().path, num_lines_per_file, self.get_prefix()))
        subprocess.call(task_call, shell=True)
        
        

class FindSequence(sciluigi.Task):
    """Find the sequence of a bed file."""
    # parameters
    genome_fasta = luigi.Parameter()
    window_size = luigi.IntParameter()
    outdir = luigi.Parameter()
    #input
    in_bed = None
    
    def out_fasta(self):
        outfile = self.get_basename() + '.%d.fasta'%self.window_size
        return sciluigi.TargetInfo(self, outfile)
    
    def get_basename(self, ):
        basename = os.path.splitext(os.path.basename(self.in_bed().path))[0] 
        return os.path.join(self.outdir, basename)
    
    
    def run(self):
        # make the fasta file
        outfile = self.out_fasta().path
        if not os.path.exists(os.path.dirname(outfile)):
            os.makedirs(os.path.dirname(outfile))
        num_up_bases = int(self.window_size/2.)
        fasta_call = ('awk -v ws=%d \'{OFS="\\t"}{print $1, $2-ws, $3+ws, $4, $5, $6}\' %s | '
                      'bedtools getfasta -s -name -fi %s -bed stdin -fo %s')%(num_up_bases, self.in_bed().path, self.genome_fasta, outfile)
        print fasta_call
        subprocess.call(fasta_call, shell=True)


class FindSSEnergy(sciluigi.Task):
    """Calculate the free energy change with and without constraint"""
    seq_length = luigi.IntParameter()
    window_size = luigi.IntParameter()
    temperature = luigi.IntParameter()
    constraint = luigi.BoolParameter()
    outdir = luigi.Parameter()
    in_fasta = None
    
    def out_rnafold(self, ):
        basename = os.path.splitext(os.path.basename(self.in_fasta().path))[0]
        if self.constraint:
            outfile = os.path.join(self.outdir, basename + '.constraint.dG.dat' )
        else:
            outfile = os.path.join(self.outdir, basename + '.dG.dat' )
        return sciluigi.TargetInfo(self, outfile)
    
    def run(self, ):
        # initialize out files
        outdir = self.outdir
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        basename = os.path.splitext(os.path.basename(self.in_fasta().path))[0]
        
        if self.constraint:
            fasta_outfile = os.path.join(outdir, basename + '.constraint.fasta')
            rnafold_outfile = os.path.join(outdir, basename + '.constraint.rnafold')
        else:
            fasta_outfile = self.in_fasta().path
            rnafold_outfile = os.path.join(outdir, basename + '.rnafold')
        
        # process params
        num_up_bases = int(self.window_size/2.)
        
        if self.constraint:
            # from fasta file, generate a constraint file.
            with open(self.in_fasta().path) as f:
                lines = f.readlines()
            # save to a file
            with open(fasta_outfile, 'w') as f:
                for name, seq in zip(lines[:-1:2], lines[1::2]):
                    key = name.strip()[1:]
                    seq = seq.strip()
                    # constraint
                    constraint_string = ''.join(['.']*num_up_bases +
                                         ['x']*self.seq_length +
                                         ['.']*(len(seq) - num_up_bases - self.seq_length))
                    f.write('>%s\n'%key)
                    f.write('%s\n'%seq)
                    f.write('%s\n'%constraint_string)
            
        # run fasta file through RNAfold
        if self.constraint:
            rnafold_call = 'cat %s | RNAfold --noPS -p0 -C -T%d > %s'%(fasta_outfile, self.temperature, rnafold_outfile)
        else:
            rnafold_call = 'cat %s | RNAfold --noPS -p0 -T%d > %s'%(fasta_outfile, self.temperature, rnafold_outfile)

        subprocess.call(rnafold_call, shell=True)

        consolidate_call = ('paste <(awk \'{if ((NR-1)%%5==0) print substr($1, 2, length($1))}\' %s) '
         '<(awk \'{n=index($0, "="); if ((NR-4)%%5==0) print substr($0, n+1, 7)}\' %s) '
         '> %s')%(rnafold_outfile, rnafold_outfile, self.out_rnafold().path)
        subprocess.call(consolidate_call, shell=True, executable='/bin/bash')