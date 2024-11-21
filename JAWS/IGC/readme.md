This is for running IGC with the samples withou RNAseq data.
Files in configure dir:

- genefilter.spec

- input.json

# Steps:
### 1> copy the file `genefilter.spec` in YOURPATH;

### 2> creat a new dir and name it with your accessions ID;

```
mkdir ABR137

```
### 3> copy the `input.json` to your target dir;

```
cp ./input.json ABR137
```

### 4> edit the `input.json`, bascially the firsly 11 lines and last line.

### 5> run jaws

```
cd ABR137
jaws submit /global/cfs/cdirs/plantbox/tbruna/compgen/compute_farm/WDL/IGC.wdl input.json dori 
```
### 5> record the job id and check the status later:

```
 jaws log 91250
 jaws status 91250
```