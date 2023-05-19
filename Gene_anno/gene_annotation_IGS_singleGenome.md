# Aims: this Tutorial is to teach people with beginner level to run Shu's gene annotation pipeline step by step


## Introduction

1. General steps

	- 1. run repeat modeler (Shu is done for Brachypodium)
	- 2. run repeat masker (annotation pipeline will automatically do)
	- 3. align RNA back to the assembled genome - do isoseq and Illumina based mRNA-seq seperately
	- 4. run IGC
    - 5. Pnagenome project, need to run pan-genome mode after annotate the individual pipeline

2. Some terminologies associated with Shu's recoded vedioe or the slides


## Preparation 

- Installed SAPS in cori

```
 ln -s /global/dna/projectdirs/plant/tools/compgen/compute_farm/saps ~/saps
```

- Test if SAP install well by print the hep message

```
~/saps/scripts/go.sh
```

- Modify your .bashrc

```
vi ~/.bashrc.ext

if [[ $NERSC_HOST == "denovo" ]]; then
   module load perl/5.24.0
   export PERL5LIB=/global/dna/projectdirs/plant/tools/compgen/lib/perl5.24:$PERL5LIB
elif [[ $NERSC_HOST == "cori" ]]; then
   export PERL5LIB=/global/dna/projectdirs/plant/tools/compgen/lib/perl5.26/lib/perl5:$PERL5LIB
else
   export PERL5LIB=/global/dna/projectdirs/plant/tools/compgen/lib/perl5.16:$PERL5LIB
fi

source ~/.bashrc.ext

```

- Get a copy of addon template from ~/saps/scripts/config/misc.conf/ and edit it to suit your case


## Start to do analysis 

Take Strain 3-7-2 as an example, we did not have isoseq, and also did not have other illumina RNAseq data, but we can use other strain's illumina data, for example, $spid= 1386221 to do the test.

### 1. Download the Illumina data from jamo with Tomas's scripts

```
mkdir Illumina

conda install python=3.8

pip install pymysql

pip install psycopg2-binary

/global/homes/t/tbruna/compgen/data_warehouse/spidFromDataWarehouse.py $spid > spid.$spid.json

grep -A2 ' "lib_protocol": "Low Input (RNA)"' spid.$spid.json|grep "lib_name"|cut -f2 -d ":"|tr -d "\","|sort|uniq |tr "\n" " " >libs.$spid.txt

cd Illumina

module load jamo

# The loop may seem too complicated -- it handles a case when multiple files are deposited

# under the same library (which does not occur here).

while read lib; do
  export lib
  jamo fetch filtered library $lib | cut -f1 -d " " | xargs -I {} bash -c 'ln -s {} ${lib}.$(basename {})'
done < <(tr " " "\n" < libs.$spid.txt)

```

### 2 Map the reads back to the assembled genome. The assembled genome can be avaible when you required from Adam


#### The reference genome path:
`/global/u2/l/llei2019/plantbox/annotation/hybridum/genome/Bhybridums_testAnnotation/Brachypodium_hybridum_var_3-7-2.fasta`

#### Edit the configure file: `upertran.addon.conf `

#### Run configure 

```
nohup ~/saps/scripts/go.sh -d jpa8 -addon ./upertran.addon.conf -name upertran &> log.config_run &
```

#### Once the configure file is done, then find the number of the batch, and then run actual job fron the log.config file

BATCH: 1379 (name:upertran) in jpa8

#### Run actual job

```
nohup ~/saps/scripts/go.sh -d jpa8 -b 1379 -wait &> run.log &
```

#### Check status:

```
~/saps/scripts/go.sh -d jpa8 -b 1379
```

#### The results file can be seen like below:

In `/global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/3-7-2/v1.0t/PERTRAN/results`

```
$ ls
total 4636274
-rw-rw-r-- 1 llei2019 llei2019       5728 Apr 26 09:25 chromgff_1290663_86197211_jpa8.gff
-rw-rw-r-- 1 llei2019 llei2019  178108592 Apr 27 00:30 pileupreadcontig_8619788_jpa8.gff.gz
-rw-rw-rw- 1 llei2019 llei2019       2806 Apr 27 00:30 mergecountread_0_86197851_jpa8.txt
-rw-rw-rw- 1 llei2019 llei2019        710 Apr 27 00:34 expresseddepth_8619789_jpa8.raw
-rw-rw-rw- 1 llei2019 llei2019  211360201 Apr 27 01:05 merge_enmall_0_86197681_jpa8.gff.gz
-rw-rw-r-- 1 llei2019 llei2019        163 Apr 27 01:09 intronlen_1295844_87116271_jpa8.raw
-rw-rw-r-- 1 llei2019 llei2019 3469845636 Apr 27 01:10 bedgraph_1295834_86197900_jpa8.bg
-rw-rw-r-- 1 llei2019 llei2019       7193 Apr 27 01:15 mexonhist_1295841_86197700_jpa8.hist
-rw-rw-rw- 1 llei2019 llei2019  299291343 Apr 27 01:17 mergefasta_0_86197731_jpa8.fa
-rw-rw-r-- 1 llei2019 llei2019       7194 Apr 27 01:17 exonhist_1295842_86197690_jpa8.hist
-rw-rw-r-- 1 llei2019 llei2019       7204 Apr 27 01:17 mexongaphist_1295843_86197720_jpa8.hist
-rw-rw-r-- 1 llei2019 llei2019  588828587 Apr 27 01:18 bigwig_1295847_87116280_jpa8.bw
-rw-rw-rw- 1 llei2019 llei2019        403 Apr 27 01:26 strandednesscheck_0_86197861_jpa8.txt
(base) 

```

#### Here I need to double check the files from the results


##### Last time, I met with Shu and discuss about the stratagies, and we think we need to treat the high confident gene model as EST and combined the EST with the output of the PERTRAN


#### Mannually combine the file with the high confident EST. High confident EST, Tomas is going to do that for me.

```
cat /global/cfs/cdirs/plantbox/annotation/Bhybridum/v2.0/annotv2.1/hiConf.transcripts.fa /global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/3-7-2/v1.0t/PERTRAN/results/mergefasta_0_86197731_jpa8.fa >/global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/3-7-2/v1.0t/annotv1.0t/Bhybridum.illuminaPlusHCest.pertran.TAs.fa

```

### 3. Run IGC
#copy Tomas's configure file into my annotation dir

```
llei2019@cori14 09:51:42 /global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/3-7-2/v1.0t/annotv1.0t 
$ cp  /global/cfs/cdirs/plantbox/annotation/Bhybridum/v2.0/annotv2.1/igc.spec .

```

#### commented this line for the first time:

```
 inactive_transcript=>'/global/cfs/cdirs/plantbox/annotation/Bhybridum/v2.0/annotv2.1/lowQ.defaultParams.rhitcnt.20.rRNA.simple.tr.txt'
```

```
 cp /global/cfs/cdirs/plantbox/annotation/Bhybridum/v2.0/annotv2.1/Jascendens.Aamericanus.Aofficinalis.Acomosus.Macuminata.Eguineensis.nr95.fa .
 1063  cp /global/cfs/cdirs/plantbox/annotation/Bhybridum/v2.0/annotv2.1/Lrigidum.Atauschii.nr95.fa .

```

```
~/saps/scripts/igc.sh -spec igc.spec -validate

Proper way to kill this one is on cori12
kill -HUP -23263
You may have to manually kill all child processes of this process: 23263
auto cmd: igc.sh -spec igc.spec
All args:
 -ppid 23263 -specfile igc.spec -validation_only 1

/global/dna/projectdirs/plant/tools/compgen/gmap-2019-09-12/bin/gmap
could query IGC run using: go.sh -d jpa8 -gn 3-7-2.v1.1
Validation done
(base) 

```
#### Here I need to deactivate the conda base enviroment and make sure the python is python2

```
 ~/saps/scripts/igc.sh -spec igc.spec
```

check job status

```
jobs
```

Kill the job

```
kill -9 jobID
```

### 4. mannually check and filter

#### 1) Creat the groundtrueth dataset with [liftoff](https://github.com/agshumate/Liftoff).

I will use the wdl scripts Tomas shared and run this in the JAW

#set up the jaws files:
```
cp /global/cfs/projectdirs/jaws/jaws-prod/jaws.conf ~/jaws.conf
chmod 600 ~/jaws.conf
```

#add the token:
```
vi ~/jaws.conf
Token = zKWpQQtbIMrAxHZ5sJpU-Od74UYa6G9vLbVVbkzUnGYB6dDynfkul9hE1NPwjWckBgHMMSvAXoVJwxqJqi-Gyini6ECKfI7vhlllGGf0lEWrYUgrZWbJgUHTmj-w_KbpuLWoHl13vgolme6eVzJI_iHhxSQYHAIvieIg1YLTmw
```

#Put those content in the `/global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/3-7-2/v1.0t/inputs.jason`
```
{
  "liftoff.processingScript": "/global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/3-7-2/v1.0t/computePeptideLengthFractions.py",
  "liftoff.referenceFasta": "/global/cfs/cdirs/plantbox/annotation/Bhybridum/v2.0/Brachypodium_hybridum_var_ABR113.mainGenome.fasta",
  "liftoff.referenceGff": "/global/cfs/cdirs/plantbox/annotation/Bhybridum/v2.0/annotv2.1/DATA/Bhybridumvar.ABR113v2.1.gene_exons.gff3",
  "liftoff.targetFasta": "/global/cfs/cdirs/plantbox/llei/annotation/hybridum/annotation/3-7-2/v1.0t/Brachypodium_hybridum_var_3-7-2.fasta",
  "liftoff.n_cpu": "32",
  "liftoff.outputName": "/global/cfs/cdirs/plantbox/annotation/Bhybridum/v2.0/evaluate_3_7_2.gff3"
}
```

```
cd /global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/3-7-2/v1.0t/
jaws inputs liftoff.wdl
```

#run job
```
jaws submit liftoff.wdl inputs.json cori

WARNING: your docker image python is not using a sha tag. call-caching will not be enabled for task: processResult.
WARNING: The runtime parameter "time" was not found in the task "lift". You should consider adding it if you ever want to share your WDL.
WARNING: The runtime parameter "time" was not found in the task "processResult". You should consider adding it if you ever want to share your WDL.
Copying Bhybridumvar.ABR113v2.1.gene_exons.gff3
100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████| 113M/113M [00:00<00:00, 527MB/s]
Copied 112912815 bytes in 0.2 seconds.
Copying Brachypodium_hybridum_var_3-7-2.fasta
100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████| 531M/531M [00:00<00:00, 680MB/s]
Copied 530948739 bytes in 0.8 seconds.
Copying computePeptideLengthFractions.py
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 2.79k/2.79k [00:00<00:00, 1.03MB/s]
Copied 2786 bytes in 0.0 seconds.
Copying Brachypodium_hybridum_var_ABR113.mainGenome.fasta
100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████| 537M/537M [00:00<00:00, 621MB/s]
Copied 537128886 bytes in 0.9 seconds.
{
    "run_id": 64099
}
(cori-prod) 

```

#check status
```
jaws status 64099

jaws log 64099
#STATUS_FROM     STATUS_TO        TIMESTAMP            COMMENT  
created          upload queued    2023-05-10 10:48:10           
upload queued    upload complete  2023-05-10 10:48:10           
upload complete  ready            2023-05-10 10:48:39           
ready            submitted        2023-05-10 10:48:45           
submitted        queued           2023-05-10 10:48:55           
queued           running          2023-05-10 10:49:26           
(cori-prod) 
```

#get the output
```
jaws get 64099 /global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/3-7-2/v1.0t

```

#I can get the path of the output:

`/global/cfs/cdirs/plantbox/llei/annotation/hybridum/annotation/3-7-2/v1.0t/annotv1.0t/DATA/Bhybridumvar.3-7-2v1.1.gene_exons.gff3`

`/global/cfs/cdirs/plantbox/annotation/Bhybridum/v2.0/evaluate_3_7_2.gff3`


#### 2) Compare my results with the groundtruth with [gffcompare](https://github.com/gpertea/gffcompare)

```
git clone https://github.com/gpertea/gffcompare

cd /some/build/dir
  git clone https://github.com/gpertea/gffcompare
  cd gffcompare
  make release

git clone https://github.com/gpertea/gclib.git

  make release

./gffcompare -R -r /global/cfs/cdirs/plantbox/llei/annotation/hybridum/annotation/3-7-2/v1.0t/annotv1.0t/DATA/Bhybridumvar.3-7-2v1.1.gene_exons.gff3 -o IGC_vs_liftoff_cmp  /global/cfs/cdirs/plantbox/annotation/Bhybridum/v2.0/evaluate_3_7_2.gff3
  75445 reference transcripts loaded.
  26 duplicate reference transcripts discarded.
  75079 query transfrags loaded.
```


```
#Summary for dataset: /global/cfs/cdirs/plantbox/annotation/Bhybridum/v2.0/evaluate_3_7_2.gff3 
#     Query mRNAs :   75079 in   53222 loci  (60059 multi-exon transcripts)
#            (11768 multi-transcript loci, ~1.4 transcripts per locus)
# Reference mRNAs :   63643 in   50409 loci  (53440 multi-exon)
# Super-loci w/ reference transcripts:    49144
#-----------------| Sensitivity | Precision  |
        Base level:    89.9     |    83.5    |
        Exon level:    78.0     |    73.1    |
      Intron level:    93.2     |    90.8    |
Intron chain level:    67.0     |    59.6    |
  Transcript level:    65.6     |    55.6    |
       Locus level:    75.9     |    71.7    |

     Matching intron chains:   35810
       Matching transcripts:   41745
              Matching loci:   38267

          Missed exons:    7789/286167  (  2.7%)
           Novel exons:   15198/316923  (  4.8%)
        Missed introns:    8362/227866  (  3.7%)
         Novel introns:   11893/234015  (  5.1%)
           Missed loci:       0/50409   (  0.0%)
            Novel loci:    3178/53222   (  6.0%)

 Total union super-loci across all input datasets: 52322 
75079 out of 75079 consensus transcripts written in IGC_vs_liftoff_cmp.annotated.gtf (0 discarded as redundant)

```


#### 3)Then load the "evaluate_3_7_2.gff3" to the Jbrowser
```
Firstly download the file to the local computer

Then open the Jbrowser

The Jbrowse is here  https://phytozome-next.jgi.doe.gov/jbprivate/index.html?data=genomes/Bhybridumvar372v1_1
12:13
user/pass is brachys_v2/brachypodiums_genepred


Then click the file and open track file of the liftoff uotput and check:

/global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/3-7-2/v1.0t/annotv1.0t/lowQ.defaultParams.plus.rRNA.simple.repetitive.tr.txt

```

```
#Tomas suggested that 

Looks reasonable. You have 0% missed loci and 6% novel loci. The novel loci should ideally mostly correspond to the low confidence genes you will filtering out
12:47
you can also try to do this comparison for coding exons only (CDS) because right now, it includes UTRs which can be very variable.
12:48
to do that, just remove the “exon” features from gff, keep only “CDS” features.
12:48
But all that’s optional. Up to you how much you want to inspect this set.
```

#Then Tomas suggested use the cds to compare 

```
./gffcompare -R -r /global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/3-7-2/v1.0t/annotv1.0t/DATA/CDS_hybridumvar.3-7-2v1.1.gene_exons.gff3 -o IGC_vs_liftoff_cmp_cds  /global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/3-7-2/v1.0t/CDS_evaluate_3_7_2.gff3


#= Summary for dataset: /global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/3-7-2/v1.0t/CDS_evaluate_3_7_2.gff3 
#     Query mRNAs :   75079 in   53222 loci  (60063 multi-exon transcripts)
#            (11768 multi-transcript loci, ~1.4 transcripts per locus)
# Reference mRNAs :   63643 in   50409 loci  (53440 multi-exon)
# Super-loci w/ reference transcripts:    49144
#-----------------| Sensitivity | Precision  |
        Base level:    89.9     |    83.5    |
        Exon level:    78.0     |    73.1    |
      Intron level:    93.2     |    90.8    |
Intron chain level:    67.0     |    59.6    |
  Transcript level:    65.6     |    55.6    |
       Locus level:    75.9     |    71.7    |

     Matching intron chains:   35808
       Matching transcripts:   41743
              Matching loci:   38265

          Missed exons:    7791/286167  (  2.7%)
           Novel exons:   15207/316947  (  4.8%)
        Missed introns:    8361/227866  (  3.7%)
         Novel introns:   11900/234024  (  5.1%)
           Missed loci:       0/50409   (  0.0%)
            Novel loci:    3178/53222   (  6.0%)

 Total union super-loci across all input datasets: 52322 
75079 out of 75079 consensus transcripts written in IGC_vs_liftoff_cmp_cds.annotated.gtf (0 discarded as redundant)

```

