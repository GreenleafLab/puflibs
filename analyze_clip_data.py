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
    genome_fasta = luigi.Parameter(default='/shr/genomes/fasta/hg38/hg38.fa')
    genome_size = luigi.Parameter(default='/shr/gSizes/hg38.genomsize')
    
    # CLIP input data
    input_bam = luigi.Parameter(default='CLIP/hPUM2/bams/input.ENCFF786ZZB.bam')
    rep1_bam = luigi.Parameter(default='CLIP/hPUM2/bams/rep1.ENCFF231WHF.bam')
    rep2_bam = luigi.Parameter(default='CLIP/hPUM2/bams/rep2.ENCFF732EQX.bam')
    rep1_peaks = luigi.Parameter(default='CLIP/hPUM2/peaks/ENCFF372VPV.rep1.st.bed')
    rep2_peaks = luigi.Parameter(default='CLIP/hPUM2/peaks/ENCFF141SVY.rep2.st.bed')
    
    # CLIP processing inputs
    len_consensus_seq = luigi.IntParameter(default=11)
    check_for_seq = luigi.Parameter(default='TGTA')
    motif_file = luigi.Parameter(default='annotations/RNAmap/hPUM2_all.motif')
    num_muts = luigi.IntParameter(default=1)
    window_size = luigi.IntParameter(default=500)
    temperature = luigi.IntParameter(default=0)
    num_random = luigi.IntParameter(default=5000000)
    randomseed = luigi.FloatParameter(default=np.nan)
    
    # sec structure processing inputs
    ss_window_size = luigi.IntParameter(default=100)
    
    # RNAMap input data
    model_param_basename = luigi.Parameter(default='annotations/RNAmap/qMotif_table_05_012318_1_')
    
    # transcript data
    tpm_cutoff = luigi.FloatParameter(default = 0.01)
    #tpm_file = luigi.Parameter(default='RNAseq/transcript_quant/rna_seq_combined.tpm.above_0.01_both.dat')
    #rnaseq_file1 = luigi.Parameter(default='RNAseq/transcript_quant/ENCFF272HJP.rep1.tsv')
    #rnaseq_file2 = luigi.Parameter(default='RNAseq/transcript_quant/ENCFF471SEN.rep2.tsv')
    #regions = luigi.Parameter(default='RNAseq/transcript_quant/exons.st.merge_transcript.above_0.01_both.bed') # the regions in which to look for motifs
    transcript_bed = luigi.Parameter(default='annotations/refseq/hg38_refGene.transcripts.st.bed')
    biomart_file = luigi.Parameter(default='annotations/ensemble_gene_converter_biomart.txt')
    
    def workflow(self):
        ####### CLIP ########
        # download CLIP data
        # TODO
        
        # process CLIP data bams
        # get the bam file of the clip data
        processclipbams = {}
        findtotalreads = {}
        outdir_bams = os.path.join(self.outdir, 'bams')
        for key, bamfile in zip(['rep1', 'rep2', 'input'], [self.rep1_bam, self.rep2_bam, self.input_bam]):
             processclipbams[key] = self.new_task('processclipbam_%s'%key, ProcessRawClipBam, bamfile=bamfile, outdir=outdir_bams)
             findtotalreads[key] = self.new_task('findtotalreads_%s'%key, FindTotalReads)
             findtotalreads[key].in_bam = processclipbams[key].out_bam

        # make bed graph of each strand of clip data
        getbedgraphs = {}
        outdir_clips = os.path.join(self.outdir, 'clip', 'bedgraphs')
        for key, processclipbam in processclipbams.items():
             getbedgraphs[key] = self.new_task('getbedgraphs_%s'%key, GetBedGraphFromBam, outdir=outdir_clips, genome_size=self.genome_size)
             getbedgraphs[key].in_bam = processclipbam.out_bam
        
        # process CLIP data peaks
        processclipbeds = self.new_task('processclipbeds%s'%key, ProcessClipPeaks, outdir=os.path.join(self.outdir, 'beds'), bedfile1=self.rep1_peaks, bedfile2=self.rep2_peaks)
        
        ####### FIND REGIONS ########
        # load RNA seq data
        downloadrna = self.new_task('downloadrna', DownloadRNAseq, outdir=os.path.join(self.outdir, 'expression'))

        # download and process gencode data
        downloadgencode = self.new_task('downloadgencode', DownloadGencode, outdir=os.path.join(self.outdir, 'annotations'))
        findtranscribedregions = self.new_task('findtranscribedregions', FindTranscribedRegions, outdir=os.path.join(self.outdir, 'annotations'))
        findtranscribedregions.in_annotations = downloadgencode.out_annotations

        # make random bed
        makerandombed = self.new_task('makerandombed', MakeRandomBed, num_random=self.num_random, outdir=os.path.join(self.outdir, 'beds'), window_size=self.len_consensus_seq, seed=self.randomseed)
        makerandombed.in_regions = findtranscribedregions.out_bed
        
        ####### DIVIDE REGIONS FOR MULTIPROCESSING ########
        # divide regions
        divideregions = self.new_task('divideregions', scltasks.DivideBed10, outdir=os.path.join(self.outdir, 'beds/split'))
        divideregions.in_bed = findtranscribedregions.out_bed

        # divide random bed
        dividebedrandom = self.new_task('dividebedrandom', scltasks.DivideBed10, outdir=os.path.join(self.outdir, 'beds/split'))
        dividebedrandom.in_bed = makerandombed.out_bed
        
        combinedataall = {}
        filterbedall = {}
        ####### FIND MOTIFS IN REGIONS ########
        # iterate over all the region beds
        num_iter = 10
        iter_values = range(num_iter+1)
        iter_values = [9,10]
        for i in iter_values:
            # normal workflow, i = 1:10:
            if i < num_iter:
                # make motif bed
                makemotifbed = self.new_task('makemotifbed_%d'%i, MakeBedFileHomer, genome=self.genome, seq_length=self.len_consensus_seq, outdir=os.path.join(self.outdir, 'beds/split/%d'%i), motif=self.motif_file)
                makemotifbed.in_regions = getattr(divideregions, 'out_bed%d'%i)
                    
                # combine motif and random bed, annotate, and filter
                combinebed = self.new_task('combinebeds_%d'%i, scltasks.CombineBeds, outdir=os.path.join(self.outdir, 'beds/split/%d'%i))
                combinebed.in_beds = [makemotifbed.out_bed, getattr(dividebedrandom, 'out_bed%d'%i)]
                filter_genetype = True
            
            elif i == num_iter:
                ###### APPEND THE CLIP PEAKS ######
                makemotifbed = self.new_task('makemotifbed_%d'%i, MakeBedFileHomer, genome=self.genome, seq_length=self.len_consensus_seq, outdir=os.path.join(self.outdir, 'beds/split/%d'%i), motif=self.motif_file)
                makemotifbed.in_regions = processclipbeds.out_bed ### CLIP PEAK REGIONS
                combinebed = makemotifbed
                filter_genetype = False
                
            # annotate and filter the clip peak motif bed
            annmotifbed = self.new_task('annmotifbed_%d'%i, AnnBedFile, genome=self.genome)
            annmotifbed.in_bed = combinebed.out_bed
            filterbed = self.new_task('filtermotifbed_%d'%i, ApplyFilterBedFile, transcript_bed=self.transcript_bed, filter_genetype=filter_genetype )
            filterbed.in_filt_dat = annmotifbed.out_filt_dat               
            filterbedall[i] = filterbed

            # split bed file into strands
            splitbedfile = self.new_task('splitbedfile_%d'%i, DividBedByStrand)
            splitbedfile.in_bed_file = filterbed.out_filt_bed
                
            # go through bedgraph files and run all clip commands
            outdir_clips = os.path.join(self.outdir, 'clip', 'split_%d'%i, 'strands')
            combinestrandsall = {}
            for key, getbedgraph in getbedgraphs.items():
                # find signal in plus strand
                clipsignalplus = self.new_task('getclipsignalplus_%s_%d'%(key, i), GetClipSignal, window_size=self.window_size, genome_size=self.genome_size, outdir=outdir_clips)
                clipsignalplus.in_bed_file = splitbedfile.out_bed_plus
                clipsignalplus.in_bg_file = getbedgraph.out_bg_plus
    
                # find signal in minus strand
                clipsignalminus = self.new_task('getclipsignalminus_%s_%d'%(key, i), GetClipSignal, window_size=self.window_size, genome_size=self.genome_size, outdir=outdir_clips)
                clipsignalminus.in_bed_file = splitbedfile.out_bed_minus
                clipsignalminus.in_bg_file = getbedgraph.out_bg_minus
                
                # combine the two
                combinestrands = self.new_task('combinestrands_%s_%d'%(key, i), CombineStrandData, outdir=os.path.join(self.outdir, 'clip', 'split_%d'%i))
                combinestrands.in_datafiles = [clipsignalplus.out_signal, clipsignalminus.out_signal]
                combinestrandsall[key] = combinestrands

            ####### FIND EXPRESSION OF TRANSCRIPTS AT MOTIF SITES ########
    
            # find the transcript count per motif site based on the annotated refseq gene and the rnaseq data
            outdir_tpm = os.path.join(self.outdir, 'expression', 'split_%d'%i)
            findmotiftpm = self.new_task('findmotiftpm_%d'%i, ProcessRNASeq, biomart_file=self.biomart_file, outdir=outdir_tpm)
            findmotiftpm.in_bed = filterbed.out_filt_bed
            findmotiftpm.in_rna1 = downloadrna.out_rna1
            findmotiftpm.in_rna2 = downloadrna.out_rna2

            ####### FIND SEQUENCE AT MOTIF SITES ########
    
            # find sequence of intervals
            outdir_seq = os.path.join(self.outdir, 'sequences', 'split_%d'%i)
            findsequence = self.new_task('findsequence_%d'%i, FindSequence, genome_fasta=self.genome_fasta,
                                         window_size=self.window_size, outdir=outdir_seq)
            findsequence.in_bed = filterbed.out_filt_bed
            
            find_seqdata = self.new_task('findseqdata_%d'%i, FindMotifSequenceData, seq_length=self.len_consensus_seq, check_for_seq=self.check_for_seq,
                                         window_size=self.window_size, outdir=outdir_seq)
            find_seqdata.in_fasta = findsequence.out_fasta
            
            ####### PREDICT EFFECTS AT MOTIF SITES #######
            outdir_model = os.path.join(self.outdir, 'effects', 'temp_%d'%(self.temperature), 'split_%d'%i)
            find_effect = self.new_task('findeffect_%d'%i, FindPredictedSeqEffect, outdir=outdir_model, model_param_basename=self.model_param_basename,
                                        temperature=self.temperature)
            find_effect.in_seqdata = find_seqdata.out_seqdata   

            ####### FIND SECONDARY STRUCTURE AT MOTIF SITES ########
        
            # find the secondary structure energy of the non-random areas
            outdir_ss = os.path.join(self.outdir, 'sec_structure', 'temp_%d'%(self.temperature), 'split_%d'%i)
            findsequence = self.new_task('findsequence_ss_%d'%i, FindSequence, genome_fasta=self.genome_fasta,
                                         window_size=self.ss_window_size, outdir=outdir_ss)
            findsequence.in_bed = filterbed.out_filt_bed
            
            findssenergy = {}
            for constraint in [False, True]:
                findssenergy[constraint] = self.new_task('findssenergy_%d_%d'%(constraint, i), FindSSEnergy, seq_length=self.len_consensus_seq,
                                                         window_size=self.ss_window_size, temperature=self.temperature, outdir=outdir_ss,
                                             constraint=constraint)
                findssenergy[constraint].in_fasta = findsequence.out_fasta

            ####### COMBINE INFO AT MOTIF SITES ########
            
            # combine data in meaningful way
            combinedata = self.new_task('combinedata_%d'%i, CombineDataAll, window_size=self.window_size, outdir=os.path.join(self.outdir, 'output', 'split_%d'%i))
            combinedata.in_counts = {key:target.out_signal for key, target in combinestrandsall.items()}
            combinedata.in_seq = find_seqdata.out_seqdata
            combinedata.in_tpm = findmotiftpm.out_motif_tpm
            combinedata.in_bed = filterbed.out_filt_bed
            combinedata.in_effect = find_effect.out_seqdata
            combinedata.in_secstructure = {key:target.out_rnafold for key, target in findssenergy.items()}
            combinedataall[i] = combinedata
        
        combinesplits = self.new_task('combinesplits', CombineSplitData, outdir = os.path.join(self.outdir, 'output'))
        combinesplits.in_data = {key:target.out_table for key, target in combinedataall.items()}
        return combinesplits




class DownloadRNAseq(sciluigi.Task):
    """Use wget to download RNA seq data."""
    # parameters
    outdir = luigi.Parameter()

    # input

    # no inputs
    
    # outputs
    def out_rna1(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, 'rna_seq_rep1.dat'))
    
    def out_rna2(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, 'rna_seq_rep2.dat'))
    
    # run
    def run(self):
        # make out directory if it doesn't exist
        dirname = os.path.dirname(self.out_rna1().path)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
            
        # load rep1 and rep2 of the RNA seq data
        subprocess.call('wget https://www.encodeproject.org/files/ENCFF272HJP/@@download/ENCFF272HJP.tsv', shell=True)
        subprocess.call('mv ENCFF272HJP.tsv %s'%self.out_rna1().path, shell=True)
        
        # load rep1 and rep2 of the RNA seq data
        subprocess.call('wget https://www.encodeproject.org/files/ENCFF471SEN/@@download/ENCFF471SEN.tsv', shell=True)
        subprocess.call('mv ENCFF471SEN.tsv %s'%self.out_rna2().path, shell=True)


class DownloadGencode(sciluigi.Task):
    """Use wget to download RNA seq data."""
    # parameters
    outdir = luigi.Parameter()

    # input

    # no inputs
    
    # outputs
    def out_annotations(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, 'gencode.v24.primary_assembly.annotation.gtf.gz'))

    # run
    def run(self):
        # make out directory if it doesn't exist
        dirname = self.outdir
        if not os.path.exists(dirname):
            os.makedirs(dirname)
            
        # download the gencode annotation file
        task1 = 'wget https://www.encodeproject.org/files/gencode.v24.primary_assembly.annotation/@@download/gencode.v24.primary_assembly.annotation.gtf.gz'
        log.info(task1)
        subprocess.call(task1, shell=True)
        subprocess.call('mv gencode.v24.primary_assembly.annotation.gtf.gz %s'%self.out_annotations().path, shell=True)


class FindTranscribedRegions(sciluigi.Task):
    """Download the exons, st, merge. """
    # parameters
    outdir = luigi.Parameter()

    # input
    in_annotations = None

    def out_bed_unprocessed(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, 'all_unprocessed.bed'))
    
    def out_bed_sorted(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, 'all_unprocessed_st.bed'))

    def out_bed(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, 'all_unprocessed_st_merged.bed'))
    
    def run(self):
        # make out directory if it doesn't exist
        dirname = os.path.dirname(self.out_bed().path)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
            
        # process to save exons
        task2 = ('gunzip -c %s | awk \'{if (substr($1, 1, 3)=="chr") print}\' | '
                 'awk \'{OFS="\\t"}{n=split($0, a, "\\t"); print $1, $4-1, $5, $3"_"(NR-1), $6, $7, a[n] }\' > %s')%(self.in_annotations().path, self.out_bed_unprocessed().path)
        log.info(task2)
        subprocess.call(task2, shell=True)
        
        # merge bed file
        task3 = 'bedtools sort -i %s > %s'%(self.out_bed_unprocessed().path, self.out_bed_sorted().path)
        log.info(task3)
        subprocess.call(task3, shell=True)
        
        task4 = 'bedtools merge -i %s -c 4 -o first |  awk \'{print}\'> %s'%(self.out_bed_sorted().path, self.out_bed().path)
        log.info(task4)
        subprocess.call(task4, shell=True)   

class MakeHomerMotif(sciluigi.Task):
    """Make a motif using HOMER and consensus seq."""
    # parameters
    seq = luigi.Parameter()
    num_muts = luigi.IntParameter()
    outdir = luigi.Parameter()
    motif_name = luigi.Parameter()
    
    # input
    # no inputs
    
    # outputs
    def out_homer_motif(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, 'motifs', '%s.%d.motif'%(self.motif_name, self.num_muts)))
    
    # run
    def run(self):
        # make out directory if it doesn't exist
        dirname = os.path.dirname(self.out_homer_motif().path)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        task_call = 'seq2profile.pl %s %d %s > %s'%(self.seq, self.num_muts,  self.motif_name, self.out_homer_motif().path)
        subprocess.call(task_call, shell=True)
        
class MakeBedFileHomer(sciluigi.Task):
    """Make a bed file of motif sites using HOMER."""
    # parameters
    genome = luigi.Parameter()
    outdir = luigi.Parameter()
    motif = luigi.Parameter()
    seq_length = luigi.IntParameter()
    
    # input
    #in_homer_motif = None
    in_regions = None
    
    # outputs
    def out_bed(self):
        basename1 = os.path.splitext(os.path.basename(self.in_regions().path))[0]
        basename2 = os.path.splitext(os.path.basename(self.motif))[0]
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, '%s.%s.bed'%(basename1, basename2)))
    
    # run
    def run(self):
        # make out directory if it doesn't exist
        dirname = os.path.dirname(self.out_bed().path)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
            
        # find the motif
        output1 = self.out_bed().path + '.tmp'
        task_call = 'annotatePeaks.pl %s %s -m %s -mbed %s -noann -nogene > /dev/null'%(self.in_regions().path, self.genome, self.motif, output1)
        log.info(task_call)
        subprocess.call(task_call, shell=True)
        
        # take out trackline and sort 
        output2 = self.out_bed().path + '.st.tmp'
        sort_call = ('tail -n+2 %s | bedtools sort -i stdin | '
                     'bedtools merge -s -d -%d -c 4,5 -o last,min -s -i stdin | '
                     'awk \'{OFS="\\t"}{print $1, $2, $3, $5 "_"NR-1, $6, $4}\' > %s')%(output1, self.seq_length, output2)
        log.info(sort_call)
        subprocess.call(sort_call, shell=True)
        
        # make sure its the same strand as the parent bed file, if strand info is given.
        output3 = self.out_bed().path 
        strand_call =  ('bedtools closest -d -a %s -b %s | '
                        'awk -F "\\t" \'{OFS="\\t"}{'
                        'distance=$NF; '
                        'motifstrand=$6; '
                        'peakstrand=$(NF-1); '
                        'if (peakstrand=="+" || peakstrand == "-") {if (peakstrand==motifstrand) print} '
                        'else print}\' | cut -f 1-6 > %s')%(output2, self.in_regions().path, output3)
        log.info(strand_call)
        subprocess.call(strand_call, shell=True)
        
        # load bed and make sure only unique sites are reported

        

class AnnBedFile(sciluigi.Task):
    """Filter the bed output from HOMER for strandedness, protein-coding, etc."""
    # parameters
    genome = luigi.Parameter()
    
    # inputs
    in_bed = None
    
    # outputs
    def out_filt_dat(self):
        outdir = os.path.dirname(self.in_bed().path)
        basename = os.path.splitext(os.path.basename(self.in_bed().path))[0]
        return sciluigi.TargetInfo(self, os.path.join(outdir, '%s.ann.dat'%(basename)))    
    # run
    def run(self):
        bedfile = self.in_bed().path
        annfile = self.out_filt_dat().path
        
        # filter
        annotate_call = 'annotatePeaks.pl %s %s > %s'%(bedfile, self.genome, annfile)

        # do calls
        subprocess.call(annotate_call, shell=True)


class ApplyFilterBedFile(sciluigi.Task):
    """Filter the bed output from HOMER for strandedness, protein-coding, etc."""
    # parameters
    transcript_bed = luigi.Parameter()
    filter_genetype = luigi.BoolParameter(default=True)
    # inputs
    in_filt_dat = None
    
    # outputs
    def out_filt_bed(self):
        outdir = os.path.dirname(self.in_filt_dat().path)
        basename = os.path.splitext(os.path.basename(self.in_filt_dat().path))[0]
        return sciluigi.TargetInfo(self, os.path.join(outdir, '%s.filt.bed'%(basename)))    
    # run
    def run(self):
        annfile = self.in_filt_dat().path
        interfile_basename = os.path.splitext(self.in_filt_dat().path)[0]
        filtfile = interfile_basename + '.filt.bed.tmp'
        
        # filter
        process_annotations = ('awk \'BEGIN {FS="\\t"}{OFS="\\t"}{'
                       'if ($NF=="") {gene_type="NA"} else {gene_type=$NF}; '
                       'n=index($8, " ("); '
                       'if (n>0) {ann=substr($8, 1, n-1); notann=substr($8, n+2, length($8)-n-2);} '
                       'else {ann=$8; notann=""}; '
                       'm=index(notann, ","); '
                       'if (m>0) {gene=substr(notann, 1, m-1); exon=substr(notann, m+2, length(notann));} '
                       'else {gene=notann; exon=""}; '
                       'print $2, $3-1, $4, $1, $6, $5, ann, gene, exon, gene_type}\'')
        
        # only include protein coding genes or the NORAD subset
        #apply_filter = ('awk \'BEGIN {FS="\\t"}{OFS="\\t"}{'
        #                 'if (($NF=="protein-coding" &&($7=="exon" || index($7, "UTR")==4)) || ($NF=="ncRNA" && $8=="%s")) print}\'')%self.nc_gene
        if self.filter_genetype:
            apply_filter = ('awk \'BEGIN {FS="\\t"}{OFS="\\t"}{'
                             'if (($NF=="protein-coding" &&($7=="exon" || index($7, "UTR")==4)) || '
                                 '($NF=="ncRNA" && $7=="non-coding")) print}\'')
            filter_call_sub = process_annotations + ' | ' + apply_filter
        else:
            filter_call_sub = process_annotations
                         
                         
        # process and apply filter
        filter_call = ('tail -n+2 %s | ' + filter_call_sub +
                       ' | grep -v chrUn | grep -v random | bedtools sort -i stdin > %s')%(annfile, filtfile)
        
        log.info(filter_call)
        subprocess.call(filter_call, shell=True)
        
        
        # find closest transcript and only keep that which aligns
        strand_call = (('bedtools closest -d -a %s -b %s | '
                        'awk -F "\\t" \'{OFS="\\t"}{'
                        'distance=$NF; '
                        'genemotif=$8; '
                        'geneclosest=$14; '
                        'strandmotif=$6; '
                        'strandgene=$16; '
                        'if ((genemotif==geneclosest && strandmotif==strandgene) || genemotif=="") print}\' | '
                       'cut -f 1-10 | awk \'{print}\'> %s')%
            (filtfile, self.transcript_bed, self.out_filt_bed().path))
        
        # do calls

        log.info(strand_call)
        subprocess.call(strand_call, shell=True)
        #subprocess.call('rm %s'%filtfile, shell=True)
        
class MakeRandomBed(sciluigi.Task):
    """Make bed file of intervals within regions of a certain length."""
    # parameters
    num_random = luigi.IntParameter()
    outdir = luigi.Parameter()
    window_size = luigi.IntParameter()
    in_regions = None
    seed = luigi.FloatParameter()
    
    # output
    def out_bed(self):
        basename1 = os.path.splitext(os.path.basename(self.in_regions().path))[0]
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, '%s.random_%.0e.bed'%(basename1, self.num_random)))
    
    # run
    def run(self):

        # make out directory if it doesn't exist
        dirname = os.path.dirname(self.out_bed().path)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
            
        # load regions
        regions = processing.load_bed(self.in_regions().path)        
        regions.loc[:, 'int_size'] = regions.stop - regions.start
        regions.loc[:, 'cum_int_size'] = np.cumsum(regions.int_size)
        
        # find n start sites
        n = self.num_random
        if not pd.isnull(self.seed):
            np.random.seed(self.seed)
        vals = np.sort(np.random.choice(np.arange(regions.int_size.sum()), size=n))
        locations = np.searchsorted(regions.cum_int_size, vals)
        strands = np.random.choice(['+', '-'], size=n)

        # go through each location and get the interval
        new_regions = []
        for i, (loc, val) in enumerate(zip(locations, vals)):
            region = regions.iloc[loc]
            diff = region.cum_int_size - val
            new_start = region.start + diff 
            new_stop = new_start + self.window_size
            new_region = pd.Series({'chrm':region.chrm,
                                    'start':new_start,
                                    'stop':new_stop,
                                    'name':'%s_%d'%(region.loc['name'], i),
                                    'strand':strands[i],
                                    'score':'.'})
            new_regions.append(new_region)
        new_regions = pd.concat(new_regions, axis=1).transpose().loc[:, variables.bed_fields]
        processing.save_bed(new_regions, self.out_bed().path)

         

class DividBedByStrand(sciluigi.Task):
    """Filter the bed output from HOMER for strandedness, protein-coding, etc."""
    # inputs
    in_bed_file = None

    
    # outputs
    def out_bed_plus(self):
        outdir = os.path.dirname(self.in_bed_file().path)
        basename = os.path.splitext(os.path.basename(self.in_bed_file().path))[0]    
        return sciluigi.TargetInfo(self, os.path.join(outdir, '%s.plus.bed'%(basename)))   
    
    def out_bed_minus(self):
        outdir = os.path.dirname(self.in_bed_file().path)
        basename = os.path.splitext(os.path.basename(self.in_bed_file().path))[0]    
        return sciluigi.TargetInfo(self, os.path.join(outdir, '%s.minus.bed'%(basename)))   

    # run
    def run(self):
        # divide into plus and minus strands
        bedfile_plus, bedfile_minus = [target().path for target in [self.out_bed_plus, self.out_bed_minus]]
        for strand, bedfile in zip(["+", "-"], [bedfile_plus, bedfile_minus]):
            task_call = ('awk \'{if ($6=="%s") print}\' %s > %s'%
                         (strand, self.in_bed_file().path, bedfile))
            subprocess.call(task_call, shell=True)  

class MakeFootprint(sciluigi.Task):
    """Given the alignments from ENCODE, process to get plus and minus strand bedgraph."""
    # inputs
    bamfile = luigi.Parameter()
    outdir = luigi.Parameter()
    cores = luigi.IntParameter()
    in_bed = None
    
    def out_vplot(self, ):
        outname = os.path.join(self.outdir, os.path.splitext(os.path.basename(self.bamfile))[0])
        return sciluigi.TargetInfo(self, outname + '.VMat')
    
    def run(self, ):
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        if not os.path.exists(self.bamfile + '.bai'):
            subprocess.call('samtools index %s'%self.bamfile, shell=True)
        outname = os.path.join(self.outdir, os.path.splitext(os.path.basename(self.bamfile))[0])

        call = 'pyatac vplot --bed %s --bam %s --not_atac --out %s --cores %d'%(self.in_bed().path, self.bamfile, outname, self.cores)
        subprocess.call(call, shell=True)

class ProcessClipPeaks(sciluigi.Task):
    """Given the two peak calls from ENCODE, processs to get a single combined file."""
    # inputs
    outdir = luigi.Parameter()
    bedfile1 = luigi.Parameter()
    bedfile2 = luigi.Parameter()
    
    def out_bed(self, ):
        filenames = [self.bedfile1, self.bedfile2]
        outname = processing.combine_filenames_split(filenames, avoid_elements=['st', 'bed']) + '.bed'
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, outname))
    
    def run(self, ):
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        
        # combine beds, sort and merge
        call = ('cat %s %s | '
                'bedtools sort -i stdin | '
                'bedtools merge -s -c 4,8 -o distinct,max -i stdin | '
                'awk \'{OFS="\\t"}{print $1, $2, $3, $5, $6, $4}\' > %s')%(self.bedfile1, self.bedfile2, self.out_bed().path)
        log.info(call)
        subprocess.call(call, shell=True)
        

class ProcessRawClipBam(sciluigi.Task):
    """Given the alignments from ENCODE, process to get plus and minus strand bedgraph."""
    # inputs
    bamfile = luigi.Parameter()
    outdir = luigi.Parameter()
    
    def out_bam(self, ):
        outname = os.path.join(self.outdir, os.path.splitext(os.path.basename(self.bamfile))[0] + '.R2.bam')
        return sciluigi.TargetInfo(self, outname)
    
    def run(self, ):
        outfile = self.out_bam().path
        if not os.path.exists(os.path.dirname(outfile)):
            os.makedirs(os.path.dirname(outfile))
        call = 'samtools view -bh -f 128 %s > %s'%(self.bamfile, self.out_bam().path)
        subprocess.call(call, shell=True)

class GetBedGraphFromBam(sciluigi.Task):
    """Given the alignments from ENCODE, process to get plus and minus strand bedgraph."""
    # inputs
    outdir = luigi.Parameter()
    genome_size = luigi.Parameter()

    in_bam = None
    
    def out_bg_plus(self, ):
        outname = self.get_basename() + '.plus.bedGraph.gz' 
        return sciluigi.TargetInfo(self, outname)

    def out_bg_minus(self, ):
        outname = self.get_basename() + '.minus.bedGraph.gz' 
        return sciluigi.TargetInfo(self, outname)
    
    def get_basename(self, ):
        return os.path.join(self.outdir, os.path.basename(os.path.splitext(self.in_bam().path)[0]))
        
    
    def run(self, ):
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
            
        for strand, target in zip(['+', '-'], [self.out_bg_plus, self.out_bg_minus]):
            extension = '.bedGraph.gz' 
            outbase= target().path[:-len(extension)]
            outcov_file = outbase + '.bg.tmp'
            cov_call = 'bedtools genomecov -ibam %s -strand %s -bg -5 > %s'%(self.in_bam().path, strand, outcov_file)
            log.info(cov_call)
            subprocess.call(cov_call, shell=True)
            
            log.info('making sure all chromosomes get tabix')
            # make a dummy bedgraph file of only chromosomes
            outdummy_file = 'temp.txt'
            dummy_bed_call = 'awk \'{OFS="\\t"}{print $1, "0", "1", "0"}\' %s > %s'%(self.genome_size, outdummy_file )
            subprocess.call(dummy_bed_call, shell=True)
            
            # concatenate
            cat_call = 'cat %s %s | bedtools sort -i stdin | bgzip > %s'%(outcov_file, outdummy_file, outbase+extension)
            log.info(cat_call)
            subprocess.call(cat_call, shell=True)
            
            # tabix index
            subprocess.call('tabix -p bed %s'%outbase+extension, shell=True)
            subprocess.call('rm %s'%outcov_file, shell=True)
               
    
class FindTotalReads(sciluigi.Task):
    """Find the total number of reads per clip bam."""
    in_bam = None
    def out_bamcount(self):
        outfile = os.path.splitext(self.in_bam().path)[0] + '.bamcount'
        return sciluigi.TargetInfo(self, outfile)
    def run(self, ):
        num_million_reads = int(subprocess.check_output('samtools view  %s | wc -l'%self.in_bam().path, shell=True).strip())/1E6
        np.savetxt(self.out_bamcount().path, [num_million_reads])
    

 
class GetClipSignal(sciluigi.Task):
    """Get the CLIP signal around motif sites."""
    # parameters
    window_size = luigi.IntParameter()
    genome_size = luigi.Parameter()
    outdir = luigi.Parameter()
    # input
    in_bed_file = None # where to find the signal
    in_bg_file = None
    # outputs
    def out_signal(self):
        return sciluigi.TargetInfo(self, '%s.tracks.txt.gz'%self.get_filename_no_ext())
    
    def get_filename_no_ext(self):
        basename = os.path.splitext(os.path.basename(self.in_bed_file().path))[0]
        basename2 = os.path.splitext(os.path.basename(self.in_bg_file().path))[0]
        return os.path.join(self.outdir, '%s.%s.%d'%(basename, basename2, self.window_size))
                                   
    def run(self):
        # make directory if it doesn't exist
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        
        # run pyatac signal per strand
        extension = '.tracks.pkl'
        outfile = self.get_filename_no_ext()
        up_amount = int(self.window_size/2.)
        call = 'pyatac signal --bed %s --bg %s --sizes %s --out %s --all --up %d --down %d --strand 6 --no_agg'%(self.in_bed_file().path, self.in_bg_file().path, self.genome_size, outfile, up_amount, up_amount)
        log.info(call)
        subprocess.call(call, shell=True)
        
        # load and add index from bed file
        bed_names = pd.read_table(self.in_bed_file().path, usecols=(3,), squeeze=True, header=None).rename('motif')
        clip_data = pd.read_csv(self.out_signal().path, compression='gzip', header=None).fillna(0).astype(float)
        clip_data.index = bed_names
        clip_data.to_csv(self.out_signal().path, compression='gzip')
        
class CombineStrandData(sciluigi.Task):
    """combine the clip signal from multiple instances of getclipsignal."""
    # parameters
    outdir = luigi.Parameter()
    
    # input
    in_datafiles = None
    
    # output
    def out_signal(self):
        basename = self.get_basename()
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, basename))

    def get_basename(self):
        name1, name2 = [os.path.basename(target().path) for target in self.in_datafiles]
        return '.'.join([s for s in name1.split('.') if s in name2.split('.')])
    
    # run
    def run(self):
        data = pd.concat([np.abs(pd.read_csv(target().path, compression='gzip', index_col=0)) for target in self.in_datafiles])
        data.to_csv(self.out_signal().path, compression='gzip')
        
class ProcessRNASeq(sciluigi.Task):
    """Map the RNAseq data to the motif data."""
    # parameters
    biomart_file = luigi.Parameter()
    outdir = luigi.Parameter()
    
    # input
    in_bed = None
    in_rna1 = None
    in_rna2 = None
    
    # output
    def out_motif_tpm(self):
        
        outfile = os.path.splitext(os.path.basename(self.in_bed().path))[0] + '.tpm.peak_id'
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, outfile))
    
    def run(self):
        # check directory exists
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        # load tpm files       
        rep1 = pd.read_table(self.in_rna1().path)
        rep2 = pd.read_table(self.in_rna2().path)
        tpm_data = pd.concat([rep1.set_index('transcript_id').TPM.rename('rep1'), rep2.set_index('transcript_id').TPM.rename('rep2')], axis=1).reset_index().rename(columns={'transcript_id':'transcript_idx'})
        tpm_data.index = [s.split('.')[0] for s in tpm_data.transcript_idx]
        tpm_combined = np.exp(np.log(tpm_data.loc[:, ['rep1', 'rep2']]).mean(axis=1))
        
        # load the biomart data mapping transcript id to refseq id
        biomart_data = pd.read_table(self.biomart_file, names=['gene_id', 'transcript_id', 'gene_name', 'refseq_id', 'refseq_nc'], header=0)
        # process to add nc to id column if id column is nan
        biomart_data.loc[:, 'refseq_comb'] = [refseq_id if not str(refseq_id)=='nan' else refseq_nc for idx, refseq_id, refseq_nc in biomart_data.loc[:, ['refseq_id', 'refseq_nc']].itertuples()]
        
        # annotate tpm data with refseq id
        biomart_data.loc[:, 'tpm'] = tpm_combined.loc[biomart_data.transcript_id].values
        # take whichever refseq id has the most tpm (or the most)
        tpm_refseq = biomart_data.groupby('refseq_id')['tpm'].max()
        
        # load bed data
        bed_data = processing.load_bed(self.in_bed().path, additional_cols=variables.motif_fields_additional)
        index = ~bed_data.refseq_id.isnull()
        bed_data.loc[index, 'tpm'] = tpm_refseq.loc[bed_data.loc[index].refseq_id].values
        log.info('%d out of %d motif sites had no TPM data'%(bed_data.tpm.isnull().sum(), len(bed_data)))
        bed_data.loc[:, ['name', 'tpm']].to_csv(self.out_motif_tpm().path, sep='\t', index=False)
        
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
        subprocess.call(fasta_call, shell=True)

class FindPredictedSeqEffect(sciluigi.Task):
    """From the file with sequence info, and predict relative affinity of that site."""
    model_param_basename = luigi.Parameter()
    outdir = luigi.Parameter()
    temperature = luigi.IntParameter()
    in_seqdata = None
    
    def out_seqdata(self):
        outfile = os.path.join(self.outdir, self.get_basename() + '.affinity.gz')
        return sciluigi.TargetInfo(self, outfile)
    
    def get_basename(self, ):
        basename = os.path.basename(os.path.splitext(self.in_seqdata().path)[0] )
        return basename

    def run(self):
        # initialize out files
        outdir = self.outdir
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        
        # load seqdata
        seqdata = pd.read_table(self.in_seqdata().path, compression='gzip', index_col=0)
        seqs = seqdata.seq.str.replace('T', 'U')
        
        # run the model
        base_params     = pd.read_csv(self.model_param_basename + 'term1.csv', index_col=0) #.stack().values
        flip_params     = pd.read_csv(self.model_param_basename + 'term2_single.csv', index_col=0) #.stack().values
        dflip_params    = pd.read_csv(self.model_param_basename + 'term2_double.csv', index_col=0, squeeze=True) #.stack().values, 10])
        coupling_params = pd.read_csv(self.model_param_basename + 'term3.csv', index_col=0, squeeze=True) #.stack().values
        
        output = [seqmodel.additive_PUF_flip_model(seq, flip_params, base_params, coupling_params, dflip_params, self.temperature, return_ensemble=True)
                  for seq in seqs]
        ddGs = pd.Series({idx:val[0] for val, idx in zip(output, seqs.index.tolist())})
        ddG_ensembles = pd.concat({idx:val[1] for val, idx in zip(output, seqs.index.tolist())}).unstack()
        
        # save
        pd.concat([ddGs.rename('ddG'), ddG_ensembles], axis=1).to_csv(self.out_seqdata().path, sep='\t', compression='gzip')  
        

        

class FindMotifSequenceData(sciluigi.Task):
    """From the fasta file, find data about the motif site."""
    seq_length = luigi.IntParameter()
    check_for_seq = luigi.Parameter()
    window_size = luigi.IntParameter()
    in_fasta = None
    
    def out_seqdata(self):
        outfile = self.get_basename() + '.seqdata.gz'
        return sciluigi.TargetInfo(self, outfile)
    
    def get_basename(self, ):
        basename = os.path.splitext(self.in_fasta().path)[0] 
        return basename

    def run(self):
        # make the fasta file
        num_up_bases = int(self.window_size/2.)
        
        # load fasta file. save the middle sequence, and save whetehr 'check for seq' is anywhere in the window'
        with open(self.in_fasta().path) as f:
            lines = f.readlines()
        # save to a dict
        fasta_dict = {}
        for name, seq in zip(lines[:-1:2], lines[1::2]):
            key = name.strip()[1:]
            fasta_dict[key] = seq.strip().upper()
        
        # now locate the sequence of the motif site itself
        motif_seqs = {}
        has_seq = {}
        where_is_seq_upstream = {}
        where_is_seq_downstream = {}
        for name, seq in fasta_dict.items():
            motif_seqs[name] = seq[num_up_bases:][:self.seq_length]
            has_seq[name] = seq.find(self.check_for_seq)>=0
            
            # downstream distance
            where_is_seq_upstream[name] = seq[:num_up_bases][::-1].find(self.check_for_seq[::-1])
            where_is_seq_downstream[name] = seq[num_up_bases+self.seq_length:].find(self.check_for_seq)              
                    
            
        motif_seqs = pd.Series(motif_seqs).rename('seq')
        where_is_seq_upstream = pd.Series(where_is_seq_upstream).replace(-1, num_up_bases+1).rename('upstream_bases_to_%s'%self.check_for_seq)
        where_is_seq_downstream = pd.Series(where_is_seq_downstream).replace(-1, num_up_bases+1).rename('downstream_bases_to_%s'%self.check_for_seq)
        pd.concat([motif_seqs, where_is_seq_upstream, where_is_seq_downstream], axis=1).to_csv(self.out_seqdata().path, sep='\t', compression='gzip')  

    

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

class CombineBamCounts(sciluigi.Task):
    """Combine al the relevant outputs into a table."""

    in_bamcounts = None

    outdir = luigi.Parameter()

    def out_bamcounts(self, ):
        return sciluigi.TargetInfo(self, self.get_basename() + '.total_bamcounts.pkl')
    
    def get_basename(self):
        name_list = list(itertools.chain(*[os.path.basename(target().path).split('.') for key, target in self.in_counts.items()]))
        _, idx = np.unique(name_list, return_index=True)
        name_unique = [name_list[i] for i in np.sort(idx)]
        basename = '.'.join([s for s in name_unique if s not in ['ann', 'filt', 'bedGraph', 'tracks', 'txt', 'gz', 'counts', 'pkl']])

        return os.path.join(self.outdir, basename)
    
    def run(self, ):
        # make directory f it doesn't exist
        outdir = self.outdir
        if not os.path.exists(outdir):
            os.makedirs(outdir)
            
        # load total counts
        bamcounts = {}
        for key, target in self.in_bamcounts.items():
            bamcounts[key] = np.loadtxt(target().path)
        bamcounts = pd.Series(bamcounts)
        bamcounts.to_pickle(self.out_bamcounts().path)       
       

class CombineDataAll(sciluigi.Task):
    """Combine al the relevant outputs into a table."""
    in_counts = None
    in_seq = None
    in_tpm = None
    in_bed = None
    in_secstructure = None
    in_effect = None
    outdir = luigi.Parameter()
    window_size = luigi.IntParameter()
    
    def out_table(self, ):
        filenames = [target().path for target in self.in_counts.values()]
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, processing.combine_filenames_split(filenames) + '.combined_data.gz'))


    
    def run(self, ):
        # make directory f it doesn't exist
        outdir = self.outdir
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        # load counts
        interval_radius = 40 #1/2 width to sum over
        offset = 15
        start_loc = int(self.window_size/2-interval_radius-offset)
        end_loc = start_loc + interval_radius*2
        counts = {}

        for key, target in self.in_counts.items():

            counts['%s'%(key)] = pd.read_csv(target().path, compression='gzip', index_col=0).iloc[:, start_loc:end_loc].sum(axis=1)
        counts = pd.concat(counts).unstack(level=0)
        
        # load seqdata
        seqdata = pd.read_table(self.in_seq().path, compression='gzip', index_col=0)

        # load seqdata_effects
        seqeffect =  pd.read_table(self.in_effect().path, compression='gzip', index_col=0)
        seqeffect.loc[:, 'flag'] = pd.Series({idx:seqmodel.flag_ensemble(row.drop('ddG')-row.ddG) for idx, row in seqeffect.iterrows()})

        kT = seqmodel.get_ddG_conversion(temperature)
        noflip_cols = [idx for idx in seqeffect if idx.find('noflip')==0]
        seqeffect.loc[:, 'ddG_noflip'] = kT*np.log(np.exp(seqeffect.loc[:, noflip_cols]/kT).sum(axis=1))

        # load bed data
        beddata = processing.load_bed(self.in_bed().path, additional_cols=variables.motif_fields_additional).set_index('name')
        
        # load tpm
        expression = pd.read_table(self.in_tpm().path, index_col=0, squeeze=True)
    
        # load ss energy
        ss_dG_data = {}
        for constraint, target in self.in_secstructure.items():
            ss_dG_data[constraint] = pd.read_table(target().path, header=None, index_col=0, squeeze=True, names=['motif', 'dG'])
        ss_dG_data = pd.concat(ss_dG_data, names=['constraint']).unstack(level=0)
        ss_dG_diff = (ss_dG_data.loc[:, True] - ss_dG_data.loc[:, False]).rename('ss_ddG')
        
        # combine
        out_data = pd.concat([beddata, counts, expression, seqdata, seqeffect.ddG, seqeffect.flag, ss_dG_diff], axis=1)
        out_data.loc[:, 'clip_signal_per_tpm'] = (out_data.rep1 + out_data.rep2)/out_data.tpm
        out_data.loc[:, 'clip_input_per_tpm'] = (out_data.input)/out_data.tpm

        out_data.to_csv(self.out_table().path, sep='\t', compression='gzip')


class CombineSplitData(sciluigi.Task):
    """Combine al the relevant outputs into a table."""
    in_data = None

    outdir = luigi.Parameter()

    
    def out_data(self, ):
        filenames = [target().path for target in self.in_data.values()]
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, processing.combine_filenames_split(filenames) + '.combined_data.gz'))
    
    def run(self, ):
        # make directory f it doesn't exist
        outdir = self.outdir
        if not os.path.exists(outdir):
            os.makedirs(outdir)
            
        # load files
        filename_dict = {key:target().path for key, target in self.in_data.items()}
        table = pd.concat({key:pd.read_table(filename, compression='gzip', index_col=0) for key, filename in filename_dict.items()}, names=['split_id', 'name'])
        table.reset_index().to_csv(self.out_data().path, compression='gzip', sep='\t')

    
if __name__ == '__main__':
    luigi.run(local_scheduler=True, main_task_cls=MyWorkflow)