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