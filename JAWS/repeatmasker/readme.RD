#Run code for repeat RepeatMasker
```
module use /global/cfs/projectdirs/jaws/modulefiles/
module load jaws
mkdir JAWS_masking
cd JAWS_masking
mkdir strain_ID
cd strain_ID
touch input.json
##copy the template of the file content and adust based on your path of the file
##the template file is [here](https://github.com/lilei1/Tutorials4students/blob/main/JAWS/repeatmasker/input.jason)
jaws submit /global/cfs/cdirs/plantbox/tbruna/compgen/compute_farm/WDL/repeatMasker.wdl  input.json dori

##record the job_ID
##check your status with
jaws status job_ID
jaws log job_ID
```