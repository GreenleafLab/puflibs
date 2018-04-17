# load per-transcript stuff.
cd '/lab/sarah/pufProject/K562_encode/RNAseq/transcript_quant'
rep1 = pd.read_table('ENCFF664LYH.tsv'); rep2 =  pd.read_table('ENCFF855OAF.tsv');
tpm = pd.concat([rep1.set_index('transcript_id').TPM.rename('rep1'), rep2.set_index('transcript_id').TPM.rename('rep2')], axis=1); tpm.loc[(np.log10(tpm + 1E-3) > -2).all(axis=1)].to_csv('rna_seq_combined.tpm.above_0.01_both.dat', sep='\t')

# load
cd '/lab/sarah/pufProject/K562_encode'
a = pd.read_table('RNAseq/transcript_quant/rna_seq_combined.tpm.above_0.01_both.dat', index_col=0);
bedFields = ['chrm', 'start', 'stop', 'name', 'strand', 'score'];
bed_file = pd.read_table('annotations/gencode/exon.st.merge_transcript.bed', header=None, names=bedFields);
bed_file.rename(columns={'name':'transcript_ids'}, inplace=True);
bed_file.loc[:, 'name'] = ['exon_%d'%i for i in range(len(bed_file))];
names = bed_file.transcript_ids.str.split(',', expand=True).replace([None], np.nan);
names_match = pd.concat({col:pd.Series(np.in1d(names.loc[:, col].dropna(), a.index.tolist()), index=names.loc[:, col].dropna().index) for col in names}, axis=1);
subset = names_match.fillna(False).any(axis=1);
cols_ordered = bedFields + ['transcript_ids'];
sub_bed = bed_file.loc[subset, cols_ordered];
sub_bed.columns = ['#' + col if i==0 else col for i, col in enumerate(sub_bed.columns.tolist())];
sub_bed.fillna('.').to_csv('RNAseq/transcript_quant//exons.st.merge.above_0.01_both.bed', index=False, sep='\t')

# find the fasta and make constraint string
ws = 20
motif_length = 8
constraint = ''.join(['.']*(ws-1) + ['x']*motif_length + ['.']*ws)
fastafile = 'analysis/hPUM2_0muts/sequence/hPUM2.1_muts.exons.st.merge_transcript.above_0.01_both.st.ann.filt.wc_%d.fasta'%ws
with open(fastafile) as f:
    lines = f.readlines()
fasta_dict = {key[1:].strip():val.strip() for key, val in zip(lines[::2], lines[1::2])}

fastafile_out = 'analysis/hPUM2_0muts/sequence/hPUM2.1_muts.exons.st.merge_transcript.above_0.01_both.st.ann.filt.wc_%d.constraint.fasta'%ws
with open(fastafile_out, 'w') as f:
    for key, val in fasta_dict.items():
        f.write('>' + key + '\n')
        f.write(val + '\n')
        f.write(constraint + '\n')
        
        
    