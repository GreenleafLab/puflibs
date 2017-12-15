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
        