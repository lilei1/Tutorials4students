{
    # Maximum intron length. If unsure, use the Iso-Seq alignments for an estimate.
    intron_len=>30000,

    # Whether the reads are paired-end.
    paired=>1,
    # Whether the paired reads are in the same file or in two separate files.
    # Files for paired-end mates in two files must have identical names with letters
    # "1" and "2" indicating the correspondence. E.g., "reads_1.fastq, reads_2.fastq".
    interleave=>1,

    # Whether the reads are stranded.
    stranded=>1,

    genome=>'/pscratch/sd/m/mcampos/Annotation/Genomes/Bhybridum/Bhyb7.fasta',
    read_root=>'/pscratch/sd/m/mcampos/Annotation/Transcriptomes/Bhybridum/Bhyb7',

    # Pseudohash specifying which reads belong to which tissue (or condition). The read locations
    # are relative to the read_root.
    reads=>{
       drought =>[qw(Bhyb7_1_T*.fastq.gz Bhyb7_2_T*.fastq.gz Bhyb7_3_T*.fastq.gz Bhyb7_4_T*.fastq.gz Bhyb7_5_T*.fastq.gz Bhyb7_6_T*.fastq.gz)],
       nutrition =>[qw(Bhyb7_*_C.fastq.gz Bhyb7_*_P.fastq.gz Bhyb7_*_N.fastq.gz)],
        #tissue_n=>[qw(read_c_*.fastq read_d_*.fastq read_e_*.fastq)],
    },

    # Advanced options below. The defaults are fine for the majority of inputs.

    # Increase this value for better specificity if the RNA coverage is high.
    minsplicedintrondepth=>1,

    # How many hits to the genome are allowed. Incerase for polyploids.
    n_hit=>8,

    # If stranded and paired, this means first mate (single end read) is anti-sense (reverse completement of mRNA).
    firstReadStrandIsmRNA=>0,

    # For cases where some trimming was done in fastq files. This is used for estimating read lengths
    variable_length_in_file=>0
}
