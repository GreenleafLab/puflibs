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
from puflibs import variables, processing, scltasks, seqmodel

log = logging.getLogger('sciluigi-interface')

class MyWorkflow(sciluigi.WorkflowTask):
    # only required parameter
    outdir = luigi.Parameter()
    cores =  luigi.IntParameter(default=1)
    
    # genome data
    genome = luigi.Parameter(default='hg38')
    genome_fasta = luigi.Parameter(default='annotations/fasta/hg38.fa')

    
    # input data
    input_bed = luigi.Parameter()
    
    # parameters
    temperature = luigi.IntParameter(default=0)
    ss_window_size = luigi.ListParameter(default=[2, 4, 6, 8, 10])
    len_consensus_seq = luigi.IntParameter(default=8)

    def workflow(self):
        ####### DEfINE inPUT ########
        bedtask = self.new_task('getbed', scltasks.FilenameToTaskOutput, filename=self.input_bed)
        
 
    
        ####### FIND SECONDARY STRUCTURE AT MOTIF SITES ########
    
        # find the secondary structure energy of the non-random areas
        findssenergy = {}
        
        for window_size in self.ss_window_size:
            outdir_ss = os.path.join(self.outdir, 'sec_structure', 'temp_%d'%(self.temperature), 'window_size_%d'%window_size)
            findsequence = self.new_task('findsequence_ss_%d'%window_size, scltasks.FindSequence,
                                         genome_fasta=self.genome_fasta,
                                         window_size=window_size,
                                         outdir=outdir_ss)
                                         
            findsequence.in_bed = bedtask.out_file
    
            
            

            findssenergy[window_size] = {}
            for constraint in [False, True]:
                findssenergy[window_size][constraint] = (
                    self.new_task('findssenergy_%d_%d'%(window_size, constraint),
                                  scltasks.FindSSEnergy, seq_length=self.len_consensus_seq,
                                  window_size=window_size, temperature=self.temperature,
                                  outdir=outdir_ss,
                                  constraint=constraint))
                findssenergy[window_size][constraint].in_fasta = findsequence.out_fasta

        ####### COMBINE INFO AT MOTIF SITES ########
        outdir_ss = os.path.join(self.outdir, 'sec_structure', 'temp_%d'%(self.temperature))

        # combine data in meaningful way
        combinedata = self.new_task('combinedata', CombineDataAll,
                                    outdir=outdir_ss,
                                    temperature=self.temperature)

        combinedata.in_bed = bedtask.out_file
        combinedata.in_secstructure = {window_size:
            {constraint:findssenergy[window_size][constraint].out_rnafold
             for constraint in [False, True]}
            for window_size in self.ss_window_size}

            
        
        #combinedata.in_secstructure = {key:target.out_rnafold for key, target in findssenergy.items()}
        


        return combinedata



class CombineDataAll(sciluigi.Task):
    """Combine al the relevant outputs into a table."""

    in_bed = None
    in_secstructure = None
    outdir = luigi.Parameter()
    temperature = luigi.FloatParameter()
    
    def out_table(self, ):
        filenames = [os.path.basename(subdict.values()[0]().path) for windowsize, subdict in self.in_secstructure.items()]
        #filenames = [os.path.basename(target().path) for target in self.in_secstructure.values()]
        filename = processing.combine_filenames_split(filenames, avoid_elements=['dat', 'constraint'])
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, filename + '.combined_data.gz'))


    
    def run(self, ):
        # make directory f it doesn't exist
        outdir = self.outdir
        if not os.path.exists(outdir):
            os.makedirs(outdir)

    
        # load bed data
        beddata = processing.load_bed(self.in_bed().path).set_index('name')
        
        
        # load ss energy
        ss_dG_data_all = {}
        for windowsize, subdict in self.in_secstructure.items():
            ss_dG_data = {}
            for constraint, target in subdict.items():
                ss_dG_data[constraint] = pd.read_table(target().path, header=None, index_col=0, squeeze=True, names=['name', 'dG'])
            ss_dG_data = pd.concat(ss_dG_data, names=['constraint']).unstack(level=0)
            ss_dG_data_all['ss_ddG_%d'%windowsize] = (ss_dG_data.loc[:, True] - ss_dG_data.loc[:, False]).rename('ss_ddG')
        ss_dG_data_all = pd.concat(ss_dG_data_all, axis=1)
        
        # combine
        out_data = pd.concat([beddata, ss_dG_data_all], axis=1).reset_index().loc[:, variables.bed_fields + ss_dG_data_all.columns.tolist()]

        out_data.to_csv(self.out_table().path, sep='\t', compression='gzip')


if __name__ == '__main__':
    luigi.run(local_scheduler=True, main_task_cls=MyWorkflow)