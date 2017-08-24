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
from puflibs import variables

log = logging.getLogger('sciluigi-interface')

class MyWorkflow(sciluigi.WorkflowTask):
    # only required parameter
    outdir = luigi.Parameter()
    cores =  luigi.IntParameter(default=1)
    
    # genome data
    genome = luigi.Parameter(default='hg38')
    genome_bowtie2 = luigi.Parameter(default='/shr/genomes/bowtie2/hg38/hg38')
    genome_fasta = luigi.Parameter(default='/shr/genomes/fasta/hg38/hg38.fa')
    genome_size = luigi.Parameter(default='/shr/gSizes/hg38.genomsize')
    
    # trasncriptome input data
    # after concatenating /lab/curtis/20160325_Pum2_HumanTranscriptome/Transcribable_S2_L001_R1_001.fastq.gz and /lab/curtis/20160325_Pum2_HumanTranscriptome/Undetermined_S0_L001_R1_001.fastq.gz
    read1_fastq = luigi.Parameter(default='transcriptome/data/transcribable_undet_R1_001.fastq.gz') 
    read2_fastq = luigi.Parameter(default='transcriptome/data/transcribable_undet_R2_001.fastq.gz')
    
    # transcriptome processing inputs
    ad1_seq = luigi.Parameter(default='CTGTCTCTTATACACATCTCCGAGCCCACGAGTCATCTCGTATGCCGTCTTCTGCTTGAA')
    ad2_seq = luigi.Parameter(default='CTGTCTCTTATACACATCTGACGCTGCCGACGTGGATCCAGGAACGTCTTCCATACAACC')
    consensus_seq = luigi.Parameter(default='TGTATATA')
    motif_name = luigi.Parameter(default='hPUM2')
    temperature = luigi.IntParameter(default=25)
    
    # fitting inputs
    # after concatenating /lab/curtis/20160325_Pum2_HumanTranscriptome/binding_curves_Pum2Pos/ANMLV_ALL_Bottom_filtered_reduced.CPfitted.pkl and
    # /lab/curtis/20160325_Pum2_HumanTranscriptome/binding_curves_allTranscribable/ANMLV_ALL_Bottom_filtered_reduced.CPfitted.pkl
    fitted_vals_filename = luigi.Parameter('transcriptome/data/transcribable_pum2pos.CPfitted.pkl')
    
    def workflow(self):
        
        # combind undetermined and transcribable 
        
        # cut adapters from fastq data
        cutadapt = self.new_task('cutadapt', CutAdapt, outdir=os.path.join(self.outdir, 'trimmed_fastqs'),
                                 read1_fastq=self.read1_fastq,
                                 read2_fastq=self.read2_fastq,
                                 ad1_seq=self.ad1_seq,
                                 ad2_seq=self.ad2_seq)
    
        # align fastq to genome to determine insert sequence
        align_fastqs = self.new_task('align_fastqs', AlignFastqs, outdir=os.path.join(self.outdir, 'aligned_genome'),
                                 genome=self.genome_bowtie2)
        align_fastqs.in_read1_fastq = cutadapt.out_read1_fastq
        align_fastqs.in_read2_fastq = cutadapt.out_read2_fastq

        # make reads from bed file into bam file
        bamtobed = self.new_task('bamtobed', BamtoBed, outdir=os.path.join(self.outdir, 'aligned_genome'))
        bamtobed.in_bam = align_fastqs.out_bam

        # filter bed file for bottom side of chip only
        #filterbedbottom = self.new_task('filterbedbottom_%d'%i, FilterBedBottom, outdir=os.path.join(self.outdir, 'trimmed_fastqs'))
        #filterbedbottom.in_bed = bamtobed.out_bed
        
        # find fasta file in regions
        fastafrombed = self.new_task('fastafrombed', FastaFromBed, outdir=os.path.join(self.outdir, 'aligned_genome'), genome=self.genome_fasta)
        fastafrombed.in_bed = bamtobed.out_bed
        
        # also make file of sequences from pear
        alignreads = self.new_task('seqfromreads', SeqFromReads, outdir=os.path.join(self.outdir, 'aligned_reads'))
        alignreads.in_read1_fastq = cutadapt.out_read1_fastq
        alignreads.in_read2_fastq = cutadapt.out_read2_fastq
        
        # filter for bottom only
        filterreads = self.new_task('filterreads', FilterBottom, outdir=os.path.join(self.outdir, 'aligned_reads'), genome=self.genome_fasta)
        filterreads.in_file = alignreads.out_annot
        filterref = self.new_task('filterref', FilterBottom, outdir=os.path.join(self.outdir, 'aligned_genome'), genome=self.genome_fasta)
        filterref.in_file = fastafrombed.out_annot        
        
        # combine
        combineseqdata = self.new_task('combineseqdata',CombineSeqData, outdir=os.path.join(self.outdir, 'final_seqs'))
        combineseqdata.in_annot_ref = filterref.out_file
        combineseqdata.in_annot_reads = filterreads.out_file
       
        # from fasta file, find consensus seqs
        #annotvariants = self.new_task('annotvariants', AnnotClusters, outdir=os.path.join(self.outdir, 'trimmed_fastqs'), consensus_seq=self.consensus_seq)
        #annotvariants.in_annot = fastafrombed.out_annot
        
        # combine 
        
        return combineseqdata
    
class CatFastq(sciluigi.Task):
    """join two sets of fastq files"""
    # parameters
    outdir = luigi.Parameter()
    read1_fastqs = luigi.ListParameter()
    read2_fastqs = luigi.ListParameter()

    def out_read1_fastq(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, 'alldata_R1.fastq.gz'))

    def out_read2_fastq(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, 'alldata_R2.fastq.gz'))
    
    def run(self):
        # make out directory if it doesn't exist
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
            
        # cat fastq
        task_call1 = 'cat %s > %s'%(' '.join([target().path for target in self.read1_fastqs]), self.out_read1_fastq().path)
        task_call2 = 'cat %s > %s'%(' '.join([target().path for target in self.read2_fastqs]), self.out_read2_fastq().path)
        
        print task_call1
        subprocess.check_call(task_call1, shell=True)
        print task_call2
        subprocess.check_call(task_call2, shell=True)       

class CutAdapt(sciluigi.Task):
    """cut adapters from paired-end fastq """
    # parameters
    outdir = luigi.Parameter()
    read1_fastq = luigi.Parameter()
    read2_fastq = luigi.Parameter()
    ad1_seq = luigi.Parameter()
    ad2_seq = luigi.Parameter()
    
    # input
    # no inputs
    
    # outputs
    def out_read1_fastq(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, self.get_basename(self.read1_fastq)))

    def out_read2_fastq(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, self.get_basename(self.read2_fastq)))    
    
    def get_basename(self, filename):
        filebasename = os.path.basename(filename)
        first_part, second_part = (filebasename[:filebasename.rfind('fastq')], filebasename[filebasename.rfind('fastq'):])
        return first_part + 'trim.' + second_part
    
    # run
    def run(self):
        # make out directory if it doesn't exist
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        
        # run cutadapt
        task_call = ('cutadapt -m 10 -q 20 -a {ad1} -A {ad2} -o {out_fastq1} -p {out_fastq2} {fastq1} {fastq2} 2>{log}'.format(
            ad1=self.ad1_seq, ad2=self.ad2_seq, out_fastq1=self.out_read1_fastq().path, out_fastq2=self.out_read2_fastq().path,
            fastq1=self.read1_fastq, fastq2=self.read2_fastq, log=os.path.join(self.outdir, 'cutadapt_out.log')))
        subprocess.check_call(task_call, shell=True)

class AlignFastqs(sciluigi.Task):
    """ align fastqs to genome."""
    # parameters
    outdir = luigi.Parameter()    
    genome = luigi.Parameter()
    
    # input
    in_read1_fastq = None
    in_read2_fastq = None
    
    #outputs
    def out_bam(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, self.get_basename(self.in_read1_fastq().path)))

    def get_basename(self, filename):
        filebasename = os.path.basename(filename)
        return filebasename[:filebasename.find('_R1')] + '.bam'
    
    # run
    def run(self):
        # make out directory if it doesn't exist
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        
        # run bowtie 2
        outsam = os.path.splitext(self.out_bam().path)[0]+'.sam'
        align_call = 'bowtie2 -X2000 -p 8 {genome} -1 {fastq1} -2 {fastq2} > {outsam} 2>{log}'.format(
            genome=self.genome, fastq1=self.in_read1_fastq().path, fastq2=self.in_read2_fastq().path, outsam=outsam, log=os.path.join(self.outdir, 'bowtie2_out.log'))
        print align_call
        subprocess.check_call(align_call, shell=True)
        
        make_bam_call = 'samtools view -bS {outsam} -o {outbam}'.format(outsam=outsam, outbam=self.out_bam().path)
        print make_bam_call
        subprocess.check_call(make_bam_call, shell=True)
        
        subprocess.call('rm {outsam}'.format(outsam=outsam), shell=True)
        
class BamtoBed(sciluigi.Task):
    """make a bed file of inserts from paired-end bam"""
    
    # parameters
    outdir = luigi.Parameter()    
    
    # input
    in_bam = None
    
    #output
    def out_bed(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, self.get_basename(self.in_bam().path) + '.bed'))
    
    def get_basename(self, filename):
        return os.path.splitext(os.path.basename(filename))[0]
    #run
    def run(self):
        
        task_call = ('samtools view -b -f 3 -q 30 %s | bedtools bamtobed -bedpe -mate1 -i - | '
        'awk \'{OFS="\t"}{if ($1==$4 ) if ($9=="+") {print $1, $2, $6, $7, $8, $9} else if ($9=="-") print $1, $5, $3, $7, $8, $9}\''
        ' > %s')%(self.in_bam().path, self.out_bed().path)
        print task_call
        subprocess.call(task_call, shell=True)
        
class FilterBottom(sciluigi.Task):
    """filter regions for those on bottom surface of chip """
    
    # parameters
    outdir = luigi.Parameter()    
    col = luigi.IntParameter(default=1)
    # input
    in_file = None
    
    #output
    def out_file(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, get_basename(self.in_file().path) + '_bottom' + os.path.splitext(self.in_file().path)[1]))
    
    #run
    def run(self):
 
        task_call = ('cat %s | awk \'{n=split($%d, a, ":"); if (substr(a[5], 1, 1)==2) print}\' > %s')%(self.in_file().path, self.col, self.out_file().path)

        print task_call
        subprocess.call(task_call, shell=True)        

class FastaFromBed(sciluigi.Task):
    """filter regions for those on bottom surface of chip """
    
    # parameters
    outdir = luigi.Parameter()
    genome = luigi.Parameter()
    
    # input
    in_bed = None
    
    #output
    def out_annot(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, self.get_basename(self.in_bed().path) + '_seqs.CPannot'))
    
    def get_basename(self, filename):
        return os.path.splitext(os.path.basename(filename))[0]
    #run
    def run(self):
        
        task_call = ('bedtools getfasta -fi {reffasta} -bed {bedfile} -fo {outfasta} -name -tab -s').format(reffasta=self.genome, bedfile=self.in_bed().path,
                                                                                                            outfasta=self.out_annot().path)

        print task_call
        subprocess.call(task_call, shell=True)

class SeqFromReads(sciluigi.Task):
    """assemble sequences using PEAR"""
    
    # parameters
    outdir = luigi.Parameter()
    
    # input
    in_read1_fastq = None
    in_read2_fastq = None
    
    #output
    def out_annot(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, self.get_basename(self.in_read1_fastq().path) + '_reads.CPannot'))

    def get_basename(self, filename):
        filebasename = os.path.basename(filename)
        return filebasename[:filebasename.find('_R1')] 
    
    #run
    def run(self):
        # make out directory if it doesn't exist
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        
        # run pear
        assembled_basename = os.path.join(self.outdir, self.get_basename(self.in_read1_fastq().path))
        task_call = ('pear -j 30 -n 15 -q 20  -f %s -r %s -o %s')%(self.in_read1_fastq().path, self.in_read2_fastq().path, assembled_basename)

        print task_call
        subprocess.call(task_call, shell=True)
        
        # change fastq to data filee
        assembled_filename = assembled_basename + '.assembled.fastq'
        
        reshape_call = ('paste <(cat %s | awk \'{if ((NR-1)%%4==0) print substr($1, 2)}\') <(cat %s | awk \'{if ((NR-2)%%4==0) print $1}\') > %s')%(assembled_filename, assembled_filename, self.out_annot().path)
        print reshape_call
        subprocess.call(reshape_call, shell=True, executable='/bin/bash')
        
        # also gzip all of the fastq files that were generated
        gzip_call = 'gzip %s/*fastq'%assembled_basename
        print gzip_call
        subprocess.call(gzip_call, shell=True)
        
class CombineSeqData(sciluigi.Task):
    """From two sequences, decide which one to take."""
    
    # parameters
    outdir = luigi.Parameter()
    threshold = luigi.IntParameter(default=100) # above this value, use the genome ref seq 
    # input
    in_annot_ref = None
    in_annot_reads = None
    
    # output
    def out_annot(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, get_basename(self.in_annot_reads().path).split('_reads')[0] + '_comb.CPannot'))
    # run
    def run(self):
        # make out directory if it doesn't exist
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
            
        # open the two refs in chunks and compare
        data = pd.concat([pd.read_table(target().path, header=None, names=['clusterID', 'sequence'], index_col=0, squeeze=True).str.upper()
                          for target in [self.in_annot_ref, self.in_annot_reads]], axis=1)
        try:
            seqs_one = pd.Series({idx:s1 if pd.isnull(s2) else s2 for idx, s1, s2 in data.loc[data.isnull().sum(axis=1)==1].itertuples()})
            genomic = pd.Series({idx:True if pd.isnull(s2) else False for idx, s1, s2 in data.loc[data.isnull().sum(axis=1)==1].itertuples()})
            # if one is not defined, use the otehr
            threshold = self.threshold
            seqs_both = {}
            for idx, seq_ref, seq_reads in data.dropna().itertuples():
    
                    
                # if the two sequences are the same, use either
                if seq_ref==seq_reads:
                    seqs_both[idx] = seq_ref 
                
                # if they are different, and the read is more than 100 bps, use 
                else:
                    if len(seq_ref) > threshold:
                        seqs_both[idx] = seq_ref
                    else:
                        seqs_both[idx] = seq_reads
            seqs_both = pd.Series(seqs_both)
            genomic = pd.concat([genomic, pd.Series(True, index=seqs_both.index)])
            all_seqs = pd.concat([seqs_one, seqs_both])
            all_seqs.to_csv(os.path.splitext(self.out_annot().path)[0] + '_seqs.CPannot', sep='\t', header=False)
            
            # also number variants
            seq_data = pd.concat([all_seqs.rename('sequence'), genomic.rename('genomic')], axis=1)
            variant_nums = all_seqs.value_counts().rename('num_clusters')
            variant_ids = pd.Series(np.arange(len(variant_nums)), index=variant_nums.index).rename('variant_number')
            variant_genomic = seq_data.groupby('sequence')['genomic'].mean()
            variant_data = pd.concat([variant_ids, variant_nums, variant_genomic], axis=1).sort_values('num_clusters', ascending=False).reset_index().rename(columns={'index':'sequence'})
            
            variant_data.loc[:, ['variant_number', 'sequence', 'num_clusters', 'genomic']].to_csv(os.path.splitext(self.out_annot().path)[0] + '.libChar', sep='\t', index=False)
            all_variant_numbers = pd.Series(variant_data.set_index('sequence').loc[all_seqs].variant_number.values, index=all_seqs.index).rename('variant_number')
            all_variant_numbers.to_csv(self.out_annot().path, sep='\t', header=False)
            all_variant_numbers.to_pickle(self.out_annot().path + '.pkl')

        except:
            ipdb.set_trace()

class AnnotClusters(sciluigi.Task):
    """find consensus seqs and mutants in fasta sequences"""
    
    # parameters
    outdir = luigi.Parameter()
    consensus_seq = luigi.Parameter()
    
    # input
    in_annot = None
    
    #output
    def out_annot(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, get_basename(self.in_annot().path) + '_%s_muts.CPannot'%self.consensus_seq))
    
    #run
    def run(self):
        
        task_call = ('python /home/sarah/puflibs/process_transcriptome_seqs.py -c {consensus} --seq {annot_seqs} --out {out}').format(
            consensus=self.consensus_seq, annot_seqs=self.in_annot().path, out=self.out_annot().path)

        print task_call
        subprocess.call(task_call, shell=True)

def get_basename(filename):
    return os.path.splitext(os.path.basename(filename))[0]

if __name__ == '__main__':
    luigi.run(local_scheduler=True, main_task_cls=MyWorkflow)