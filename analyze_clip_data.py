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
from puflibs import variables, processing, scltasks

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
    
    # CLIP processing inputs
    consensus_seq = luigi.Parameter(default='TGTATATA')
    motif_name = luigi.Parameter(default='hPUM2')
    num_muts = luigi.IntParameter(default=1)
    window_size = luigi.IntParameter(default=200)
    temperature = luigi.IntParameter(default=0)
    num_random = luigi.IntParameter(default=1000000)
    
    # transcript data
    tpm_cutoff = luigi.FloatParameter(default = 0.01)
    tpm_file = luigi.Parameter(default='RNAseq/transcript_quant/rna_seq_combined.tpm.above_0.01_both.dat')
    rnaseq_file1 = luigi.Parameter(default='RNAseq/transcript_quant/ENCFF272HJP.rep1.tsv')
    rnaseq_file2 = luigi.Parameter(default='RNAseq/transcript_quant/ENCFF471SEN.rep2.tsv')
    regions = luigi.Parameter(default='RNAseq/transcript_quant/exons.st.merge_transcript.above_0.01_both.bed') # the regions in which to look for motifs
    transcript_bed = luigi.Parameter(default='annotations/refseq/hg38_refGene.transcripts.st.bed')
    biomart_file = luigi.Parameter(default='annotations/ensemble_gene_converter_biomart.txt')
    
    def workflow(self):
        # load RNA seq data
        downloadrna = self.new_task('downloadrna', DownloadRNAseq, outdir=self.outdir)
        getexpressedtranscripts = self.new_task('getexpressedtranscripts', GetExpressedGenes, outdir=self.outdir, tpm_cutoff=self.tpm_cutoff)
        getexpressedtranscripts.in_rna1 = downloadrna.out_rna1
        getexpressedtranscripts.in_rna2 = downloadrna.out_rna2
        
        # download and process gencode data
        downloadgencode = self.new_task('downloadgencode', DownloadGencode, outdir=self.outdir)
        findtranscribedregions = self.new_task('findtranscribedregions', FindTranscribedRegions, outdir=self.outdir)
        findtranscribedregions.in_annotations = downloadgencode.out_annotations
        
        # download and process CLIP data
        downloadCLIP = self.new_task('downloadgencode', DownloadGencode, outdir=self.outdir)
        
        
        # make homer motif
        makehomermotif = self.new_task('makehomermotif', MakeHomerMotif, num_muts=self.num_muts, outdir=self.outdir, seq=self.consensus_seq, motif_name=self.motif_name)
        
        # make motif bed
        makemotifbed = self.new_task('makemotifbed', MakeBedFileHomer, genome=self.genome, outdir=os.path.join(self.outdir, 'beds'))
        makemotifbed.in_homer_motif = makehomermotif.out_homer_motif
        makemotifbed.in_regions = findtranscribedregions.out_bed
        
        # make random bed
        makerandombed = self.new_task('makerandombed', MakeRandomBed, num_random=self.num_random, outdir=os.path.join(self.outdir, 'beds'))
        makerandombed.in_regions = findtranscribedregions.out_bed
        
        # combine motif and random bed, annotate, and filter
        combinebed = self.new_task('combinebeds', scltasks.CombineBeds, outdir=os.path.join(self.outdir, 'beds'))
        combinebed.in_beds = [makemotifbed.out_bed, makerandombed.out_bed]
        
        # annotate and filter the combined motif bed
        annmotifbed = self.new_task('annmotifbed', AnnBedFile, genome=self.genome)
        annmotifbed.in_bed = combinebed.out_bed
        filterbed = self.new_task('filtermotifbed', ApplyFilterBedFile, transcript_bed=self.transcript_bed )
        filterbed.in_filt_dat = annmotifbed.out_filt_dat
        
        """
        # make footprint
        makefootprintall = {}
        outdir_footprint = os.path.join(self.outdir, 'footprints')
        for key, bamfile in zip(['rep1', 'rep2', 'input'], [self.rep1_bam, self.rep2_bam, self.input_bam]):            
            makefootprint = self.new_task('makefootprint_%s'%key, MakeFootprint, bamfile=bamfile,
                                          outdir=outdir_footprint, cores=self.cores)
            makefootprint.in_bed = filtermotifbed.out_filt_bed
            makefootprintall[key] = makefootprint
        """

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


        # split bed file into strands
        splitbedfile = self.new_task('splitbedfile', DividBedByStrand)
        splitbedfile.in_bed_file = filterbed.out_filt_bed
            
        # go through bedgraph files and run all clip commands
        outdir_clips = os.path.join(self.outdir, 'clip', 'strands')
        combinestrandsall = {}
        for key, getbedgraph in getbedgraphs.items():
            # find signal in plus strand
            clipsignalplus = self.new_task('getclipsignalplus_%s'%(key), GetClipSignal, window_size=self.window_size, genome_size=self.genome_size, outdir=outdir_clips)
            clipsignalplus.in_bed_file = splitbedfile.out_bed_plus
            clipsignalplus.in_bg_file = getbedgraph.out_bg_plus

            # find signal in minus strand
            clipsignalminus = self.new_task('getclipsignalminus_%s'%(key), GetClipSignal, window_size=self.window_size, genome_size=self.genome_size, outdir=outdir_clips)
            clipsignalminus.in_bed_file = splitbedfile.out_bed_minus
            clipsignalminus.in_bg_file = getbedgraph.out_bg_minus
            
            # combine the two
            combinestrands = self.new_task('combinestrands_%s'%(key), CombineStrandData, outdir=os.path.join(self.outdir, 'clip'))
            combinestrands.in_datafiles = [clipsignalplus.out_signal, clipsignalminus.out_signal]
            combinestrandsall[key] = combinestrands

        """
        # find the transcript count per motif site based on the annotated refseq gene and the rnaseq data
        outdir_tpm = os.path.join(self.outdir, 'expression')
        findmotiftpm = self.new_task('findmotiftpm_%s'%key_rand, ProcessRNASeq, rnaseq_file1=self.rnaseq_file1,
                                     rnaseq_file2=self.rnaseq_file2, biomart_file=self.biomart_file, outdir=outdir_tpm)
        findmotiftpm.in_bed = filterbed.out_filt_bed
        motiftpmboth[key_rand] = findmotiftpm
        
        # find sequence of intervals
        outdir_seq = os.path.join(self.outdir, 'sequences')
        findsequence = self.new_task('findsequence_%s'%key_rand, FindSequence, genome_fasta=self.genome_fasta,
                                     consensus_seq=self.consensus_seq, check_for_seq=self.consensus_seq[:4],
                                     window_size=self.window_size, outdir=outdir_seq)
        findsequence.in_bed = filterbed.out_filt_bed
        find_seqdata = self.new_task('findseqdata_%s'%key_rand, FindMotifSequenceData, consensus_seq=self.consensus_seq, check_for_seq=self.consensus_seq[:4],
                                     window_size=self.window_size, outdir=outdir_seq)
        find_seqdata.in_fasta = findsequence.out_fasta

        
        # find the secondary structure energy of the non-random areas
        secstructuremotif = {}
        outdir_ss = os.path.join(self.outdir, 'sec_structure_%d'%self.temperature)
        for window_size in [50, 75, 100]:
            findsequence = self.new_task('findsequence_%d'%(window_size), FindSequence, genome_fasta=self.genome_fasta,
                                         consensus_seq=self.consensus_seq, check_for_seq=self.consensus_seq[:4],
                                         window_size=window_size, outdir=outdir_ss)
            findsequence.in_bed = filtermotifbed.out_filt_bed
            
            for constraint in [False, True]:
                findssenergy = self.new_task('findssenergy_%d_%d'%(constraint, window_size), FindSSEnergy, consensus_seq=self.consensus_seq, window_size=window_size, temperature=self.temperature, outdir=outdir_ss,
                                             constraint=constraint)
                findssenergy.in_fasta = findsequence.out_fasta
                secstructuremotif[(window_size, constraint)] = findssenergy

        
        # combine data in meaningful way
        combinedata = self.new_task('combinedata', CombineDataAll, outdir=os.path.join(self.outdir, 'output'))
        combinedata.in_counts      = {key:target.out_signal for key, target in motifcountsboth['motif'].items()}
        combinedata.in_counts_rand = {key:target.out_signal for key, target in motifcountsboth['random'].items()}
        combinedata.in_seq      = sequenceboth['motif'].out_seqdata
        combinedata.in_seq_rand = sequenceboth['random'].out_seqdata
        combinedata.in_tpm      = motiftpmboth['motif'].out_motif_tpm
        combinedata.in_tpm_rand = motiftpmboth['random'].out_motif_tpm
        combinedata.in_bamcounts =  {key:target.out_bamcount for key, target in findtotalreads.items()}
        combinedata.in_bed = filtermotifbed.out_filt_bed
        combinedata.in_bed_rand = filterrandombed.out_filt_bed
        combinedata.in_secstructure = {key:target.out_rnafold for key, target in secstructuremotif.items()}
        ""
        # make plots
        makeplots = self.new_task('makeplots', MakePlots, outdir=os.path.join(self.outdir, 'plots'))
        makeplots.in_counts = {key:target.out_signal for key, target in motifcountsboth['motif'].items()}
        makeplots.in_counts_rand = {key:target.out_signal for key, target in motifcountsboth['random'].items()}
        makeplots.in_signal = {key:target.out_signal for key, target in motifsignalboth['motif'].items()}
        makeplots.in_signal_rand = {key:target.out_signal for key, target in motifsignalboth['random'].items()}
        makeplots.in_sequence = filtermotifbed.out_filt_bed
        makeplots.in_bed_rand = filterrandombed.out_filt_bed
        ""
        

        #makeplots.in_signal = [target.out_signal for target in motifsignalboth['motif']]
        #makeplots.in_signal_rand = [target.out_signal for target in motifsignalboth['random'].values()]
        
        #return combinedata
        """
        return combinestrandsall


class DownloadRNAseq(sciluigi.Task):
    """Use wget to download RNA seq data."""
    # parameters
    outdir = luigi.Parameter()

    # input

    # no inputs
    
    # outputs
    def out_rna1(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, 'expression_data', 'rna_seq_rep1.dat'))
    
    def out_rna2(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, 'expression_data', 'rna_seq_rep2.dat'))
    
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
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, 'annotations', 'gencode.v24.primary_assembly.annotation.gtf.gz'))

    # run
    def run(self):
        # make out directory if it doesn't exist
        dirname = os.path.dirname(self.out_rna1().path)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
            
        # download the gencode annotation file
        task1 = 'wget https://www.encodeproject.org/files/gencode.v24.primary_assembly.annotation/@@download/gencode.v24.primary_assembly.annotation.gtf.gz'
        log.info(task1)
        subprocess.call(task1, shell=True)
        subprocess.call('mv gencode.v24.primary_assembly.annotation.gtf.gz %s'%self.out_annotations().path, shell=True)
        
  
class GetExpressedGenes(sciluigi.Task):
    """Process RNA seq expression data and save transcripts with some expression"""
    # parameters
    outdir = luigi.Parameter()
    tpm_cutoff = luigi.FloatParameter()
    
    # input
    in_rna1 = None
    in_rna2 = None
    
    # outputs
    def out_tpm(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, 'expression_data', 'transcripts_tpm_%.2f.dat'%self.tpm_cutoff))
    
    # run
    def run(self):
        # make out directory if it doesn't exist
        dirname = os.path.dirname(self.out_tpm().path)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
            
        # load rep1 and rep2 of the RNA seq data
        rep1 = pd.read_table(self.in_rna1().path)
        rep2 = pd.read_table(self.in_rna2().path)
        tpm = pd.concat([rep1.set_index('transcript_id').TPM.rename('rep1'),
                         rep2.set_index('transcript_id').TPM.rename('rep2')], axis=1);
        eps = self.tpm_cutoff/10.
        tpm.loc[(np.log10(tpm + eps) > np.log10(self.tpm_cutoff)).all(axis=1)].to_csv(self.out_tpm().path, sep='\t')

class FindExonicRegions(sciluigi.Task):
    """Download the exons, st, merge. """
    # parameters
    outdir = luigi.Parameter()

    # input
    

    def out_bed_unprocessed(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, 'annotations', 'exons_unprocessed.bed'))

    def out_bed_processed(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, 'annotations', 'exons_processed.bed'))
    
    def out_bed_sorted(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, 'annotations', 'exons_processed_st.bed'))

    def out_bed(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, 'annotations', 'exons_processed_st_merged.bed'))
    
    def run(self):
        # make out directory if it doesn't exist
        dirname = os.path.dirname(self.out_annotations().path)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
            

        # process to save exons
        task2 = ('gunzip -c %s | awk \'{if ($3=="exon") print}\' | '
                 'awk \'{OFS="\\t"}{n=split($0, a, "\\t"); print $1, $4-1, $5, $3"_"(NR-1), $6, $7, a[n] }\' > %s')%(self.out_annotations().path, self.out_bed_unprocessed().path)
        log.info(task2)
        subprocess.call(task2, shell=True)
        
        ##### process to make annotations more reasonable. #####
        with open(self.out_bed_unprocessed().path) as f:
            lines = f.readlines()
        
        # go through lines and save fields
        log.info('processing annotation lines')
        line_dict = {}
        for i, line in enumerate(lines):
            line_series = pd.Series(line.split('\t')[:-1], index=variables.bed_fields)
            line_annots = pd.Series({val.split()[0]:val.split()[1].replace('"', '') for val in line.split('\t')[-1].strip('\n;').split('; ')})
            line_dict[i] = pd.concat([line_series, line_annots])
        
        # concat fields
        log.info('concatenating annotations')
        line_table = pd.concat(line_dict).unstack()
        
        # reorder columns
        cols = line_table.columns.tolist()
        cols_ordered = variables.bed_fields + [col for col in cols if col not in variables.bed_fields]
        line_table = (line_table.loc[:, cols_ordered].sort_values(['chrm', 'start', 'stop'])).copy()
        
        # save output with fields prepended by comment for use in bed
        line_table.columns = ['#'+col if i==0 else col for i, col in enumerate(cols_ordered)]
        line_table.to_csv(self.out_bed_processed().path, sep='\t', index=False)
        
        # merge bed file
        task3 = 'bedtools sort -i %s > %s'%(self.out_bed_processed().path, self.out_bed_sorted().path)
        log.info(task3)
        subprocess.call(task3, shell=True)
        
        task4 = 'bedtools merge -i %s -c 20 -o distinct |  awk \'BEGIN {print "%s"} {print}\'> %s'%(self.out_bed_sorted().path, '#' + '\\t'.join(variables.bed_fields), self.out_bed().path)
        log.info(task4)
        subprocess.call(task4, shell=True)        

class FindTranscribedRegions(sciluigi.Task):
    """Download the exons, st, merge. """
    # parameters
    outdir = luigi.Parameter()

    # input
    in_annotations = None

    def out_bed_unprocessed(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, 'annotations', 'all_unprocessed.bed'))
    
    def out_bed_sorted(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, 'annotations', 'all_unprocessed_st.bed'))

    def out_bed(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, 'annotations', 'all_unprocessed_st_merged.bed'))
    
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
        
        task4 = 'bedtools merge -i %s -c 4 -o first |  awk \'BEGIN {print "%s"} {print}\'> %s'%(self.out_bed_sorted().path, '#' + '\\t'.join(variables.bed_fields), self.out_bed().path)
        log.info(task4)
        subprocess.call(task4, shell=True)   


class FindExpressedRegions(sciluigi.Task):
    """Download the exons, st, merge. """
    # parameters
    outdir = luigi.Parameter()

    # input
    in_bed = None
    in_tpm = None
    
    def out_bed(self):
        basename2 = os.path.splitext(os.path.basename(self.in_tpm().path))[0]
        basename = os.path.splitext(os.path.basename(self.in_bed().path))[0]
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, basename + '_' + basename2 + '.bed'))


    
    def run(self):
        # make out directory if it doesn't exist
        dirname = os.path.dirname(self.out_bed().path)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        
        # load merged exonic regions
        regions = processing.load_bed(self.in_bed().path)
        regions.rename(columns={'name':'transcript_ids'}, inplace=True);
        regions.loc[:, 'name'] = ['exon_%d'%i for i in range(len(regions))];
        names = regions.transcript_ids.str.split(',', expand=True).replace([None], np.nan);
        
        # load transcripts with tpm > cutoff
        a = pd.read_table(self.in_tpm().path, index_col=0);
        
        # find those regions with transcript ids in a
        names_match = pd.concat({col:pd.Series(np.in1d(names.loc[:, col].dropna(), a.index.tolist()), index=names.loc[:, col].dropna().index) for col in names}, axis=1);
        subset = names_match.fillna(False).any(axis=1);
        cols_ordered = variables.bed_fields + ['transcript_ids'];
        sub_bed = regions.loc[subset, cols_ordered]
        processing.save_bed(sub_bed.fillna('.'), self.out_bed().path)

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
    
    # input
    in_homer_motif = None
    in_regions = None
    
    # outputs
    def out_bed(self):
        basename1 = os.path.splitext(os.path.basename(self.in_regions().path))[0]
        basename2 = os.path.splitext(os.path.basename(self.in_homer_motif().path))[0]
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, '%s.%s.bed'%(basename1, basename2)))
    
    # run
    def run(self):
        # make out directory if it doesn't exist
        dirname = os.path.dirname(self.out_bed().path)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
            
        # find the motif
        output1 = self.out_bed().path + '.tmp'
        task_call = 'annotatePeaks.pl %s %s -m %s -mbed %s -noann -nogene > /dev/null'%(self.in_regions().path, self.genome, self.in_homer_motif().path, output1)
        subprocess.call(task_call, shell=True)
        
        # take out trackline and sort
        output2 = self.out_bed().path 
        sort_call = 'tail -n+2 %s | bedtools sort -i stdin | awk \'BEGIN {print "%s"} {OFS="\t"}{center=int(($2+$3)*0.5); print $1, center, center+1, $4 "_"NR-1, $5, $6}\' > %s'%(output1, '#' + '\\t'.join(variables.bed_fields), output2)
        subprocess.call(sort_call, shell=True)
        

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
        filter_call = ('tail -n+2 %s | '
                       'awk \'BEGIN {FS="\\t"}{OFS="\\t"}{n=index($8, " ("); ann=substr($8, 1, n-1); notann=substr($8, n+2, length($8)-n-2); m=index(notann, ","); gene=substr(notann, 1, m-1); exon=substr(notann, m+2, length(notann)); '
                        'if ($NF=="protein-coding" &&(ann=="exon" || index(ann, "UTR")==4))  print $2, $3-1, $4, $1, $6, $5, ann, gene, exon}\' | '
                        'grep -v chrUn | grep -v random | bedtools sort -i stdin > %s')%(annfile, filtfile)
        
        # find closest transcript and only keep that which aligns
        strand_call = ('bedtools closest -d -a %s -b %s | awk -F "\t" \'{OFS="\t"}{if (($8==$13) && ($6==$15)) print}\' | '
                       'cut -f 1-9 > %s')%(filtfile, self.transcript_bed, self.out_filt_bed().path)
        
        # do calls
        subprocess.call(filter_call, shell=True)
        subprocess.call(strand_call, shell=True)
        subprocess.call('rm %s'%filtfile, shell=True)
        
class MakeRandomBed(sciluigi.Task):
    """Make bed file of intervals within regions of a certain length."""
    # parameters
    num_random = luigi.IntParameter()
    outdir = luigi.Parameter()
    
    in_regions = None
    
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
        np.random.seed(0)
        vals = np.sort(np.random.choice(np.arange(regions.int_size.sum()), size=n))
        locations = np.searchsorted(regions.cum_int_size, vals)
        strands = np.random.choice(['+', '-'], size=n)

        # go through each location and get the interval
        new_regions = []
        for i, (loc, val) in enumerate(zip(locations, vals)):
            region = regions.iloc[loc]
            diff = region.cum_int_size - val
            new_start = region.start + diff 
            new_stop = new_start + 1
            new_region = pd.Series({'chrm':region.chrm,
                                    'start':new_start,
                                    'stop':new_stop,
                                    'name':'%s_%d'%(region.loc['name'], i),
                                    'strand':strands[i],
                                    'score':region.score})
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
            task_call = 'awk \'{if ($6=="%s") print}\' %s > %s'%(strand, self.in_bed_file().path, bedfile)
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
            dummy_bed_call = 'awk \'{OFS="\t"}{print $1, "0", "1", "0"}\' %s > %s'%(self.genome_size, outdummy_file )
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

        basename = os.path.splitext(os.path.basename(self.in_bed_file().path))[0]
        basename2 = os.path.splitext(os.path.basename(self.in_bg_file().path))[0]
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, '%s.%s.%d.tracks.pkl'%(basename, basename2, self.window_size)))
                                   
    def run(self):
        # make directory if it doesn't exist
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        
        # run pyatac signal per strand
        extension = '.tracks.pkl'
        outfile = self.out_signal().path[:-len(extension)]
        up_amount = int(self.window_size/2.)
        call = 'pyatac signal --bed %s --bg %s --sizes %s --out %s --all --up %d --down %d --strand 6 --no_agg'%(self.in_bed_file().path, self.in_bg_file().path, self.genome_size, outfile, up_amount, up_amount)
        log.info(call)
        subprocess.call(call, shell=True)
        
        # load and add index from bed file
        bed_names = pd.read_table(self.in_bed_file().path, usecols=(3,), squeeze=True, header=None).rename('motif')
        clip_data = pd.read_csv(outfile + '.tracks.txt.gz', compression='gzip', header=None).fillna(0).astype(float)
        clip_data.index = bed_names
        clip_data.to_pickle(self.out_signal().path)
        
class CombineStrandData(sciluigi.Task):
    """combine the clip signal from multiple instances of getclipsignal."""
    # parameters
    outdir = luigi.Parameter()
    
    # input
    in_datafiles = None
    
    # output
    def out_signal(self):
        basename = self.get_basename()
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, basename + '.counts.pkl'))

    def get_basename(self):
        name1, name2 = [os.path.basename(target().path) for target in self.in_datafiles]
        return '.'.join([s for s in name1.split('.')[:-1] if s in name2.split('.')])
    
    # run
    def run(self):
        data = pd.concat([np.abs(pd.read_pickle(target().path)) for target in self.in_datafiles])
        data.to_pickle(self.out_signal().path)
        
class ProcessRNASeq(sciluigi.Task):
    """Map the RNAseq data to the motif data."""
    # parameters
    rnaseq_file1 = luigi.Parameter()
    rnaseq_file2 = luigi.Parameter()
    biomart_file = luigi.Parameter()
    outdir = luigi.Parameter()
    
    # input
    in_bed = None
    
    # output
    def out_motif_tpm(self):
        
        outfile = os.path.splitext(os.path.basename(self.in_bed().path))[0] + '.tpm.peak_id'
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, outfile))
    
    def run(self):
        # check directory exists
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        # load tpm files       
        rep1 = pd.read_table(self.rnaseq_file1)
        rep2 = pd.read_table(self.rnaseq_file2)
        tpm_data = pd.concat([rep1.set_index('transcript_id').TPM.rename('rep1'), rep2.set_index('transcript_id').TPM.rename('rep2')], axis=1);
        tpm_data = tpm_data.loc[(tpm_data >= 0.01).all(axis=1)].reset_index().rename(columns={'transcript_id':'transcript_idx'})
        tpm_data.index = [s.split('.')[0] for s in tpm_data.transcript_idx]
        tpm_combined = np.exp(np.log(tpm_data.loc[:, ['rep1', 'rep2']]).mean(axis=1))
        
        # load the biomart data mapping transcript id to refseq id
        biomart_data = pd.read_table(self.biomart_file, names=['gene_id', 'transcript_id', 'gene_name', 'refseq_id', 'refseq_nc'], header=0)
        # process to add nc to id column if id column is nan
        biomart_data.loc[:, 'refseq_id'] = [refseq_id if not str(refseq_id)=='nan' else refseq_nc for idx, refseq_id, refseq_nc in biomart_data.loc[:, ['refseq_id', 'refseq_nc']].itertuples()]
        
        # annotate tpm data with refseq id
        biomart_data.loc[:, 'tpm'] = tpm_combined.loc[biomart_data.transcript_id].values
        # take whichever refseq id has the most tpm (or the most)
        tpm_refseq = biomart_data.groupby('refseq_id')['tpm'].max()
        
        # load bed data
        bed_data = pd.read_table(self.in_bed().path, header=None, names=variables.bed_fields + variables.motif_fields_additional)
        bed_data.loc[:, 'tpm'] = tpm_refseq.loc[bed_data.refseq_id].values
        log.info('%d out of %d motif sites had no TPM data'%(bed_data.tpm.isnull().sum(), len(bed_data)))
        bed_data.loc[:, ['name', 'tpm']].to_csv(self.out_motif_tpm().path, sep='\t', index=False)
        


class NormalizeClip(sciluigi.Task):
    """Sum the read counts around motif sites and normalize by tpm"""
    # parameters

    # input
    in_signal = None
    in_tpm = None
    in_bam = None
    
    def out_signal(self):
        outfile = os.path.splitext(self.in_signal().path)[0] +'.signal.pkl'
        return sciluigi.TargetInfo(self, outfile)
    
    def run(self):
        # load signal
        counts = pd.read_pickle(self.in_signal().path).sum(axis=1)
        tpm = pd.read_table(self.in_tpm().path, index_col=0, squeeze=True)
        
        # num unique reads
        log.info('Finding the number of total clip reads')
        num_million_reads = int(subprocess.check_output('samtools view  %s | wc -l'%self.in_bam().path, shell=True).strip())/1E6
        signal = np.divide(counts, tpm.loc[counts.index.tolist()]).astype(float)
        min_signal = signal.replace(0, np.nan).dropna().min()
        log_signal = np.log2(signal+min_signal/2.)
        log_signal.to_pickle(self.out_signal().path)

class NormalizeClipSMinput(sciluigi.Task):
    """Sum the read counts around motif sites and normalize by tpm"""
    # parameters

    # input
    in_signal = None
    in_signal_input = None
    in_bam = None
    in_bam_input = None
    
    def out_signal(self):
        outfile = os.path.splitext(self.in_signal().path)[0] +'.input_norm.signal.pkl'
        return sciluigi.TargetInfo(self, outfile)
    
    def run(self):
        # load signal
        counts = pd.read_pickle(self.in_signal().path).sum(axis=1)
        counts_input = pd.read_pickle(self.in_signal_input().path).sum(axis=1)
        
        # num unique reads
        log.info('Finding the number of total clip reads')
        num_million_reads = int(subprocess.check_output('samtools view  %s | wc -l'%self.in_bam().path, shell=True).strip())/1E6
        
        # num unique reads
        log.info('Finding the number of total clip reads in input')
        num_million_reads_input = int(subprocess.check_output('samtools view  %s | wc -l'%self.in_bam_input().path, shell=True).strip())/1E6
                
        
        signal = np.divide(counts/num_million_reads, counts_input/num_million_reads_input).astype(float).replace(np.inf, np.nan)
        min_signal = signal.replace(0, np.nan).dropna().min()
        log_signal = np.log2(signal+min_signal/2.)
        log_signal.to_pickle(self.out_signal().path)
        
class FindSequence(sciluigi.Task):
    """Find the sequence of a bed file."""
    # parameters
    genome_fasta = luigi.Parameter()
    consensus_seq = luigi.Parameter()
    check_for_seq = luigi.Parameter()
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

class FindMotifSequenceData(sciluigi.Task):
    """From the fasta file, find data about the motif site."""
    consensus_seq = luigi.Parameter()
    check_for_seq = luigi.Parameter()
    window_size = luigi.IntParameter()
    in_fasta = None
    
    def out_seqdata(self):
        outfile = self.get_basename() + '.%d.seqdata.pkl'%self.window_size
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
        for name, seq in fasta_dict.items():
            motif_seqs[name] = seq[num_up_bases:][:len(self.consensus_seq)]
            has_seq[name] = seq.find(self.check_for_seq)>=0
        motif_seqs = pd.Series(motif_seqs).rename('seq')
        has_seq = pd.Series(has_seq).rename('seq_bool')
        pd.concat([motif_seqs, has_seq], axis=1).to_pickle(self.out_seqdata().path)  

    

class FindSSEnergy(sciluigi.Task):
    """Calculate the free energy change with and without constraint"""
    consensus_seq = luigi.Parameter()
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
                                         ['x']*len(self.consensus_seq) +
                                         ['.']*(len(seq) - num_up_bases -len(self.consensus_seq)))
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

        
       

class CombineDataAll(sciluigi.Task):
    """Combine al the relevant outputs into a table."""
    in_counts = None
    in_counts_rand = None
    in_seq = None
    in_seq_rand = None
    in_tpm = None
    in_tpm_rand = None
    in_bamcounts = None
    in_bed = None
    in_bed_rand = None
    in_secstructure = None
    outdir = luigi.Parameter()
    
    def out_table(self, ):
        return sciluigi.TargetInfo(self, self.get_basename() + '.combined_data.pkl')

    def out_bamcounts(self, ):
        return sciluigi.TargetInfo(self, self.get_basename() + '.total_bamcounts.pkl')
    
    def get_basename(self):
        name_list = list(itertools.chain(*[os.path.basename(target().path).split('.') for key, target in self.in_counts.items()]))
        _, idx = np.unique(name_list, return_index=True)
        name_unique = [name_list[i] for i in np.sort(idx)]
        basename = '.'.join([s for s in name_unique if s not in ['ann', 'filt', 'bedGraph', 'tracks', 'counts', 'pkl']])

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
            
        # load counts
        counts = {}
        for key_rand, target_dict in zip(['motif', 'random'], [self.in_counts, self.in_counts_rand]):
            for key, target in target_dict.items():
                counts[(key_rand, key)] = pd.read_pickle(target().path).sum(axis=1)
        counts = pd.concat(counts).unstack(level=1)
        
        # load seqdata
        seqdata = {}
        for key_rand, target in zip(['motif', 'random'], [self.in_seq, self.in_seq_rand]):
            seqdata[key_rand] = pd.read_pickle(target().path)
        seqdata = pd.concat(seqdata)

        # load bed data
        beddata = {}
        for key_rand, target in zip(['motif', 'random'], [self.in_bed, self.in_bed_rand]):
            beddata[key_rand] = pd.read_table(target().path, header=None, names=variables.bed_fields + variables.motif_fields_additional, index_col='name')
        beddata = pd.concat(beddata)
        
        # load tpm
        expression = {}
        for key_rand, target in zip(['motif', 'random'], [self.in_tpm, self.in_tpm_rand]):
            expression[key_rand] = pd.read_table(target().path, index_col=0, squeeze=True)
        expression = pd.concat(expression)
        
        # load ss energy
        ss_dG_data = {}
        for (window_size, constraint), target in self.in_secstructure.items():
            ss_dG_data[(window_size, constraint)] = pd.read_table(target().path, header=None, index_col=0, squeeze=True, names=['motif', 'dG'])
        ss_dG_data = pd.concat(ss_dG_data, names=['window_size', 'constraint']).unstack(level=1)
        ss_dG_diff = (ss_dG_data.loc[:, True] - ss_dG_data.loc[:, False]).unstack(level=0)
        ss_dG_diff.rename(columns={s:'ss_%d'%s for s in ss_dG_diff}, inplace=True)
        
        ssdata = pd.concat({'motif':ss_dG_diff, 'random':pd.DataFrame(index=beddata.loc['random'].index, columns=ss_dG_diff.columns)})
        
        out_data = pd.concat([beddata, counts, expression, seqdata, ssdata], axis=1)
        out_data.to_pickle(self.out_table().path)
            
        
   
    

class MakePlots(sciluigi.Task):
    """Make the plots you want to make!"""
    # parameters
    outdir = luigi.Parameter()
    # inputs
    in_counts = None
    in_counts_rand = None
    
    in_signal = None
    in_signal_rand = None
    in_bed = None
    in_bed_rand = None


    
    def out_plot_names(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, 'normalized_signal_consensus_versus_random_no_neg.pdf'))
    
    def run(self):
        outdir = self.outdir
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        # load the counts   
        data_agg = pd.concat([pd.concat({('motif', key):pd.read_pickle(target().path).mean() for key, target in self.in_counts.items()},names=['site_type', 'sample']),
                              pd.concat({('random', key):pd.read_pickle(target().path).mean() for key, target in self.in_counts_rand.items()},names=['site_type', 'sample'])])
            
        #data_agg = pd.concat({key:pd.read_pickle(filename).mean() for key, filename in count_filenames.items()}, names=['site_type', 'sample'])
        data_agg.unstack().transpose().plot(linewidth=1, figsize=(3,3))
        plt.xlabel('distance from motif center (bp)')
        plt.ylabel('mean read counts per site')
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, 'aggregate_counts_per_site.pdf'))
        
        # laod the signal
        data_agg = pd.concat([pd.concat({('motif', key):pd.read_pickle(target().path) for key, target in self.in_signal.items()},names=['site_type', 'sample']),
                              pd.concat({('random', key):pd.read_pickle(target().path) for key, target in self.in_signal_rand.items()},names=['site_type', 'sample'])])
        data_agg = data_agg.unstack('sample')
        
        # laod bed file
        bed_data = pd.concat({'motif':pd.read_table(self.in_bed().path, header=None, names=variables.bed_fields + variables.motif_fields_additional, index_col='name'),
                              'random':pd.read_table(self.in_bed_rand().path, header=None, names=variables.bed_fields + variables.motif_fields_additional, index_col='name')}, names=['site_type', 'motif'])
        data_agg.loc[:, 'score'] = bed_data.score
        
        data_melt = pd.melt(data_agg, id_vars=['score']).replace(np.inf, np.nan).dropna(subset=['value'])
        g = sns.FacetGrid(data=data_melt, col='sample', hue='score'); g.map(sns.distplot, 'value', kde=False, norm_hist=True)
        g.add_legend()
        plt.savefig(os.path.join(outdir, 'normalized_signal_consensus_versus_random.pdf'))

        g = sns.FacetGrid(data=data_melt, col='sample', hue='score'); g.map(sns.distplot, 'value', norm_hist=True, bins=np.linspace(-3, 10, 15), kde_kws={'clip':[-3, 10]})
        g.add_legend()
        plt.savefig(os.path.join(outdir, 'normalized_signal_consensus_versus_random_no_neg.pdf'))
        


    
if __name__ == '__main__':
    luigi.run(local_scheduler=True, main_task_cls=MyWorkflow)