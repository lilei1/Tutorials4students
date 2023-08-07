# Aims: this Tutorial is to teach people with beginner level to run Shu's gene annotation pipeline in parralle way and also in Perlmutter!

The basical steps are the same as the IGC with the single genome in [gene_annotation_IGS_singleGenome.md](https://github.com/lilei1/Tutorials4students/blob/main/Gene_anno/gene_annotation_IGS_singleGenome.md).

## Step1: Run Pertran one by one

Since run Pertran in paralle way will cause problems and it will be hard to debug, we still run Pertran one by one. 

In the below example, I am going to run three strains: Bd2-2-2, Bstat7,and Bstat8. Bstat7 and Bstat8 both have Illumina RNAseq data from the nutriant deficiency experiment; however,  Bd2-2-2 did not have.

### 1. Download the RNAseq data from jamo.

RNAseq data information:

```
#Illumina RNAseq:
$spid= 1386228

Bstat7
=====================
 PROJECT INFORMATION
=====================
Scientific Program: Plant
Organism Name: Brachypodium hybridum Bsta7
Principal Investigator: Pilar Catalan
Proposal ID: 503504
Sequencing Project ID: 1386228
Final Deliverable Project ID: 1386201
Proposal Name: Genomic Characterization of the Brachypodium Polyploid Model to Unravel Bases of Success of Polyploidy in Flowering Plants
Sequencing Project Name: Brachypodium hybridum Bsta7 Nutrient Deficit Gene Expression Profiling

=====================
 SAMPLE SUMMARY     
=====================
1=libraryName   2=SampleId      3=rawReads      4=filteredReads 5=sampleName    6=conditionNumber       7=groupName     8=sequencerType 9=runType       10=fileUsed
NNPBS   305601  87154218        77506394        Bsta7_1_C       1       Bsta7_C NovaSeq S4      2x151   52748.4.439898.GTTGTAGC-GTTGTAGC.filter-RNA.fastq.gz
NNPBT   305602  154960430       133790078       Bsta7_2_C       1       Bsta7_C NovaSeq S4      2x151   52748.4.439898.ACCACGAT-ACCACGAT.filter-RNA.fastq.gz
NNPBU   305603  99526268        86240706        Bsta7_3_C       1       Bsta7_C NovaSeq S4      2x151   52748.4.439898.CTGAAGCT-CTGAAGCT.filter-RNA.fastq.gz
NNPBW   305604  85355444        76372150        Bsta7_4_C       1       Bsta7_C NovaSeq S4      2x151   52748.4.439898.TCCTTAGC-TCCTTAGC.filter-RNA.fastq.gz
NNPBX   305605  104778498       95523446        Bsta7_1_N       2       Bsta7_N NovaSeq S4      2x151   52748.4.439898.ACGACTTG-ACGACTTG.filter-RNA.fastq.gz
NNPBY   305606  84771870        71351328        Bsta7_2_N       2       Bsta7_N NovaSeq S4      2x151   52748.4.439898.TGCTTCCA-TGCTTCCA.filter-RNA.fastq.gz
NNPBZ   305607  79987886        73964932        Bsta7_3_N       2       Bsta7_N NovaSeq S4      2x151   52748.4.439898.GCATACAG-GCATACAG.filter-RNA.fastq.gz
NNPCA   305608  89118852        76302472        Bsta7_4_N       2       Bsta7_N NovaSeq S4      2x151   52748.4.439898.CTTCGTTC-CTTCGTTC.filter-RNA.fastq.gz
NNPCB   305609  102125846       91736588        Bsta7_1_P       3       Bsta7_P NovaSeq S4      2x151   52748.4.439898.TGCGAACT-TGCGAACT.filter-RNA.fastq.gz
NNPCC   305610  89841346        79452110        Bsta7_2_P       3       Bsta7_P NovaSeq S4      2x151   52748.4.439898.CATCGTGA-CATCGTGA.filter-RNA.fastq.gz
NNPCG   305611  107648402       99836764        Bsta7_3_P       3       Bsta7_P NovaSeq S4      2x151   52748.4.439898.TCTCCGAT-TCTCCGAT.filter-RNA.fastq.gz
NNPCH   305612  78176248        46969284        Bsta7_4_P       3       Bsta7_P NovaSeq S4      2x151   52748.4.439898.AGAGTAGC-AGAGTAGC.filter-RNA.fastq.gz

=====================
MAPPING STATISTICS
=====================
1=libraryName   2=totalFragments        3=mappedFragments       4=assignedFragments     5=unassignedAmbiguous   6=unassignedNoFeatures  7=unassignedSecondaryHits       8=ratioStrandedness
NNPBS   38753197        38049887        36246676        367022  1436189 0       0.9790
NNPBT   66895039        65591490        62721423        628213  2241854 0       0.9792
NNPBU   43120353        42290560        40572430        361213  1356917 0       0.9777
NNPBW   38186075        37407417        35766193        272060  1369164 0       0.9788
NNPBX   47761723        46882846        44807716        444272  1630858 0       0.9791
NNPBY   35675664        35022415        33606668        324248  1091499 0       0.9775
NNPBZ   36982466        36183352        34484589        374837  1323926 0       0.9783
NNPCA   38151236        37441600        35835820        347929  1257851 0       0.9769
NNPCB   45868294        45037321        42945143        457366  1634812 0       0.9779
NNPCC   39726055        39024934        37355043        352947  1316944 0       0.9775
NNPCG   49918382        49021428        46884196        492373  1644859 0       0.9768
NNPCH   23484642        22923418        22032021        103962  787435  0       0.9773

NOTE: mappedFragments - fragments aligned to the reference genome.
assignedFragments - aligned fragments assigned to the gene features.                 
unassignedAmbiguous - aligned fragments overlapping two or more features.
unassignedNoFeatures - aligned fragments not overlapping with any features included in the annotation.                 
unassignedSecondaryHits - fragments marked as second alignment in the FLAG field in SAM/BAM input.

=====================
 RELEASE FILES       
=====================
Portal:
    http://genome.jgi.doe.gov/portal/lookup?keyName=jgiProjectId&keyValue=1386201&app=Info
   Note: If you have problems logging in, please contact your Project Manager.

#Bsta 8:
Scientific Program: Plant
Organism Name: Brachypodium hybridum Bsta8
Principal Investigator: Pilar Catalan
Proposal ID: 503504
Sequencing Project ID: 1386229
Final Deliverable Project ID: 1386203
Proposal Name: Genomic Characterization of the Brachypodium Polyploid Model to Unravel Bases of Success of Polyploidy in Flowering Plants
Sequencing Project Name: Brachypodium hybridum Bsta8 Nutrient Deficit Gene Expression Profiling

=====================
 SAMPLE SUMMARY     
=====================
1=libraryName   2=SampleId      3=rawReads      4=filteredReads 5=sampleName    6=conditionNumber       7=groupName     8=sequencerType 9=runType       10=fileUsed
NNPCN   305613  111699942       101173432       Bsta8_1_C       1       Bsta8_C NovaSeq S4      2x151   52748.4.439898.CAACACCT-CAACACCT.filter-RNA.fastq.gz
NNPCO   305614  96890558        88578998        Bsta8_2_C       1       Bsta8_C NovaSeq S4      2x151   52748.4.439898.TTGTGTGC-TTGTGTGC.filter-RNA.fastq.gz
NNPCP   305615  93769606        86078604        Bsta8_1_N       2       Bsta8_N NovaSeq S4      2x151   52748.4.439898.ACTGCTAG-ACTGCTAG.filter-RNA.fastq.gz
NNPCS   305616  85536710        76684100        Bsta8_2_N       2       Bsta8_N NovaSeq S4      2x151   52748.4.439898.GTCATCGA-GTCATCGA.filter-RNA.fastq.gz
NNPCT   305617  110427398       98242690        Bsta8_3_N       2       Bsta8_N NovaSeq S4      2x151   52748.4.439898.GAAGAGGT-GAAGAGGT.filter-RNA.fastq.gz
NNPCU   305618  119240876       110268632       Bsta8_4_N       2       Bsta8_N NovaSeq S4      2x151   52748.4.439898.AGTCGACA-AGTCGACA.filter-RNA.fastq.gz
NNPCW   305619  89589370        83322736        Bsta8_1_P       3       Bsta8_P NovaSeq S4      2x151   52748.4.439898.AAGCACTG-AAGCACTG.filter-RNA.fastq.gz
NNPCX   305620  105696094       95964012        Bsta8_2_P       3       Bsta8_P NovaSeq S4      2x151   52748.4.439898.GTCGGTAA-GTCGGTAA.filter-RNA.fastq.gz
NNPCY   305621  97923044        57429112        Bsta8_3_P       0       Bsta8_P NovaSeq S4      2x151   52748.4.439898.CTCCTAGA-CTCCTAGA.filter-RNA.fastq.gz
NNPCZ   305622  127549544       111279650       Bsta8_4_P       3       Bsta8_P NovaSeq S4      2x151   52751.4.439598.TCTGAGAG-TCTGAGAG.filter-RNA.fastq.gz

=====================
MAPPING STATISTICS
=====================
1=libraryName   2=totalFragments        3=mappedFragments       4=assignedFragments     5=unassignedAmbiguous   6=unassignedNoFeatures  7=unassignedSecondaryHits       8=ratioStrandedness
NNPCY   28714556        27909920        25907297        208091  1794532 0       0.9758
NNPCN   50586716        48977790        46764375        419642  1793773 0       0.9784
NNPCO   44289499        43181170        41245190        412349  1523631 0       0.9785
NNPCP   43039302        42026213        40123755        430016  1472442 0       0.9789
NNPCS   38342050        37312802        35783239        403070  1126493 0       0.9792
NNPCT   49121345        47921327        45824650        532523  1564154 0       0.9788
NNPCU   55134316        53954938        51461507        707887  1785544 0       0.9782
NNPCW   41661368        40868837        38831571        516651  1520615 0       0.9746
NNPCX   47982006        47078126        44877269        542626  1658231 0       0.9733
NNPCZ   55639825        53477095        45743594        320703  7412798 0       0.9621

```

Install python 3.8  because spidFromDataWarehouse.py need this specific python version

```
module load python
mamba create --yes --name python3_8 python=3.8 -c conda-forge
mamba activate python3_8

```

Download the data:
```
#download Bstat7
cd /global/u2/l/llei2019/plantbox/annotation/hybridum/transcriptomes/Illumina 
spid=1386228
mkdir Bstat7

pip install psycopg2-binary
python3 /global/homes/t/tbruna/compgen/data_warehouse/spidFromDataWarehouse.py $spid > spid.$spid.json

 grep -A2 ' "lib_protocol": "Low Input (RNA)"' spid.$spid.json|grep "lib_name"|cut -f2 -d ":"|tr -d "\","|sort|uniq |tr "\n" " " >libs.$spid.txt

module load jamo
while read lib; do   export lib;   jamo fetch filtered library $lib | cut -f1 -d " " | xargs -I {} bash -c 'ln -s {} ${lib}.$(basename {})'; done < <(tr " " "\n" < libs.$spid.txt)


#download Bstat8
spid=1386229
python3 -m pip install pymysql
python3 /global/homes/t/tbruna/compgen/data_warehouse/spidFromDataWarehouse.py $spid > spid.$spid.json

 grep -A2 ' "lib_protocol": "Low Input (RNA)"' spid.$spid.json|grep "lib_name"|cut -f2 -d ":"|tr -d "\","|sort|uniq |tr "\n" " " >libs.$spid.txt

module load jamo
while read lib; do   export lib;   jamo fetch filtered library $lib | cut -f1 -d " " | xargs -I {} bash -c 'ln -s {} ${lib}.$(basename {})'; done < <(tr " " "\n" < libs.$spid.txt)

```

### 2 run pertran one by one for Bstat7 and Bstat8:

The reference genome:
`/global/cfs/cdirs/plantbox/llei/annotation/hybridum/genome/assembliesForLi_May18.23/Brachypodium_hybridum_var_Bsta7.fasta`

`/global/cfs/cdirs/plantbox/llei/annotation/hybridum/genome/assembliesForLi_May18.23/Brachypodium_hybridum_var_Bsta8.fasta`

The Ilumina RNAseq:

`/global/cfs/cdirs/plantbox/llei/annotation/hybridum/transcriptomes/Illumina/Bstat7`

`/global/cfs/cdirs/plantbox/llei/annotation/hybridum/transcriptomes/Illumina/Bstat8`

Bstat7:

#copy the upertran configure file into the current dir:
```
cd ~/plantbox/annotation/hybridum/annotation/migc/Bstat7_V3
mkdir PERTRAN

cp /global/cfs/cdirs/plantbox/llei/annotation/hybridum/annotation/3-7-2/v1.0t/PERTRAN/upertran.addon.conf . 

# then revise the configure file, the revised configure file can be seen [here](https://github.com/lilei1/Tutorials4students/blob/main/Gene_anno/configures/Pertran/Bstat7/upertran.addon_perlmutter.txt):

!!!! keep in mind the "-name" should be unique for each strain!!!

nohup ~/saps/scripts/go.sh -d jpa9 -addon ./upertran.addon_perlmutter.conf -name upertran_Bstat7 &> log.config_run &
```

#I've got the error and ask me to install perl module of DBD::Pg module, I asked Tomas and he told me that I need to revise my bashrc.ext file:

```
if [[ $LMOD_SYSHOST == "perlmutter" ]]; then
   export PERL5LIB=/global/dna/projectdirs/plant/tools/compgen/lib/perl5.26.perlmutter/lib/perl5:$PERL5LIB
elif [[ $NERSC_HOST == "denovo" ]]; then
   export PERL5LIB=/global/dna/projectdirs/plant/tools/compgen/lib/perl5.24:$PERL5LIB
   shopt -s direxpand
   module load emacs
elif [[ $NERSC_HOST == "cori" ]]; then
   export PERL5LIB=/global/dna/projectdirs/plant/tools/compgen/lib/perl5.26/lib/perl5:$PERL5LIB
   shopt -s direxpand
elif [[ $NERSC_HOST == "datatran" ]]; then
   export PERL5LIB=/global/dna/projectdirs/plant/tools/compgen/lib/perl5.26/lib/perl5:$PERL5LIB
else
   export PERL5LIB=/global/dna/projectdirs/plant/tools/compgen/lib/perl5.16:$PERL5LIB
fi
```

```
souce ~/.bashrc.ext

```
####  Rerun the command line 

```
nohup ~/saps/scripts/go.sh -d jpa9 -addon ./upertran.addon_perlmutter.conf -name upertran_Bstat7 &> log.config_run &

```

It works!!!

#### Check job status:

```
jobs
```

```
#Then find the job from the log file:
BATCH: 173 (name:upertran_Bstat7) in jpa9

```

#Then run this:

```
nohup ~/saps/scripts/go.sh -d jpa9 -b 173 -wait &> run.log &
[2] 346818
```


#check status!!!

```
~/saps/scripts/go.sh -d jpa9 -b 173

```

### Bstat8

# copy the configure file from Bstat7 to this dir, then revise the configure file, the revised configure file can be seen [here](https://github.com/lilei1/Tutorials4students/blob/main/Gene_anno/configures/Pertran/Bstat8/upertran.addon_perlmutter.txt):

```
nohup ~/saps/scripts/go.sh -d jpa9 -addon ./upertran.addon_perlmutter.conf -name upertran_Bstat8 &> log.config_run &

```

######  BATCH: 165 (name:upertran_Bstat8) in jpa9


```
nohup ~/saps/scripts/go.sh -d jpa9 -b 165 -wait &> run.log & 
```

##### check the job status

```
~/saps/scripts/go.sh -d jpa9 -b 165
```

#Create the concatenate est with the RNAseq and high confident transcripts

```
cat /global/cfs/cdirs/plantbox/annotation/Bhybridum/v2.0/annotv2.1/hiConf.transcripts.fa /global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/migc/Bstat7/PERTRAN/results/mergefasta_0_5169371_jpa9.fa  >/global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/migc/Bstat7/BhybridumBstat7.illuminaPlusHCest.pertran.TAs.fa

cat /global/cfs/cdirs/plantbox/annotation/Bhybridum/v2.0/annotv2.1/hiConf.transcripts.fa /global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/migc/Bstat8/PERTRAN/results/mergefasta_0_4774351_jpa9.fa  >/global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/migc/Bstat8/BhybridumBstat8.illuminaPlusHCest.pertran.TAs.fa

```


##### Check if all of the scafoold belong to the corresponding genome!!!

```
grep ">" /global/u2/l/llei2019/plantbox/annotation/hybridum/genome/assembliesForLi_May18.23/Brachypodium_hybridum_var_Bsta8.fasta >geno

samtools idxstats leaf_tissue_control_jpa9_628335.bam >bam1
 
diff -y bam1 geno

```


### 3 run MIGC 1st round

# you need to have three configure files:

- [common.igc.spec](https://github.com/lilei1/Tutorials4students/blob/main/Gene_anno/configures/migc/common.igc.spec)

- [genefilter.spec](https://github.com/lilei1/Tutorials4students/blob/main/Gene_anno/configures/migc/genefilter.spec)

- [test3strains.migc.spec](https://github.com/lilei1/Tutorials4students/blob/main/Gene_anno/configures/migc/test3strains.migc.spec)

#edit them 

#run validation

```
/global/homes/l/llei2019/saps/scripts/migc -d jpa9 -single -spec test3strains.migc.spec -validate
```
##### run 1st round

```
/global/homes/l/llei2019/saps/scripts/migc -spec test3strains.migc.spec -single &>log_1stRun&

```

##### #run PAV mode for migc!!!

cp test3strains.migc.spec [m.test3strains.migc.spec](https://github.com/lilei1/Tutorials4students/blob/main/Gene_anno/configures/migc/m.test3strains.migc.spec)

#then uncommented the lines below:

```
##value can be migc-workdir-full-path/*/hiConf.proteome.nr.fa
    makenonselfpanproteomefrom=>{default=>'/global/cfs/cdirs/plantbox/llei/annotation/hybridum/annotation/migc/*/hiConf.proteome.nr.fa'},

    high_confidence_transcript_file=>{default=>'./hiConf.ovlps'}, #relative and in each IGC own work_dir

```

```
/global/homes/l/llei2019/saps/scripts/migc -spec m.test3strains.migc.spec -single &>log_2ndRun&
[3] 32425
```

##### Actually I need to run the liftoff to double check if the LQ genes are good!!


##### run the  inactivation step (filtering the LQ genes):

cp m.test3strains.migc.spec [in.test3strains.migc.spec](https://github.com/lilei1/Tutorials4students/blob/main/Gene_anno/configures/migc/in.test3strains.migc.spec)

#commented the lines below:

```
#makenonselfpanproteomefrom=>{default=>'/global/cfs/cdirs/plantbox/llei/annotation/hybridum/annotation/migc/*/hiConf.proteome.nr.fa'},
    #high_confidence_transcript_file=>{default=>'./hiConf.ovlps'}, #relative and in each IGC own work_dir
```

#uncommented the lines:

```
inactive_transcript=>{default=>'./lowQ.defaultParams.plus.rRNA.simple.repetitive.tr.txt'}
```

```
/global/homes/l/llei2019/saps/scripts/migc -spec in.test3strains.migc.spec -single &>log_5ndRun&

```

!!!!Keep in mind:

#Tomas suggested I need to use the unique name for each run
