This is the tutorial for running JAWS_PERTRAN 

Steps:

### 1> Create a dir called JAWS_PERTRAN and in this file, creat a file with accession name

```
mkdir JAWS_PERTRAN
cd JAWS_PERTRAN
mkdir Bhyb7
```
### 2> Revise the pertran.spec according to your accessions informations, basically the this section

```
 genome=>'/pscratch/sd/m/mcampos/Annotation/Genomes/Bhybridum/Bhyb7.fasta',
    read_root=>'/pscratch/sd/m/mcampos/Annotation/Transcriptomes/Bhybridum/Bhyb7',

    # Pseudohash specifying which reads belong to which tissue (or condition). The read locations
    # are relative to the read_root.
    reads=>{
       drought =>[qw(Bhyb7_1_T*.fastq.gz Bhyb7_2_T*.fastq.gz Bhyb7_3_T*.fastq.gz Bhyb7_4_T*.fastq.gz Bhyb7_5_T*.fastq.gz Bhyb7_6_T*.fastq.gz)],
       nutrition =>[qw(Bhyb7_*_C.fastq.gz Bhyb7_*_P.fastq.gz Bhyb7_*_N.fastq.gz)],
        #tissue_n=>[qw(read_c_*.fastq read_d_*.fastq read_e_*.fastq)],
    },
```

### 3> Run perl script to generate the `pertran.input.json` file 

```
/global/cfs/cdirs/plantbox/compgen/compute_farm/saps/tools/makePERTRANinputSpec.pl -spec pertran.spec >pertran.input.json
n input fastq groups/tissues: 2
Done
```

### 4> Then submit the jaws job 

```
module load jaws
jaws submit /global/cfs/cdirs/plantbox/compgen/compute_farm/WDL/ptr_Assemble.wdl  pertran.input.json dori
```
