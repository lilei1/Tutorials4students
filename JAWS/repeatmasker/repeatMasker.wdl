version 1.0

import "task/splitFasta.wdl" as splitFa
import "task/chunkChromosome.wdl" as chunk
import "task/utils.wdl" as util
import "taskParams.wdl" as taskParams

#-e(ngine) [crossmatch|rmblast]
workflow RepeatMasker {
  input {
   String compgen
   File   genomeFasta
   File   repeat_lib
   Int?   chunk_size = 1000000
   Int?	 num_chunk = 4
   Int?   overlap = 1000
   String? rpm_args = "-xsmall -gccalc -cutoff 400"
   # doejgi/repeatmasker:4.1.2.fixed.trf - works around the issue of a trf failure due to a long filepath
   String? repeatMaskerTag = "@sha256:d24945453f80523efc5acd723fb73e5a57a0a9922369fa8e8fdb405224833a12"
   String? compgenTag = "@sha256:26e87557d0111a81709403465fb41cfa4cd72add1916d768fe071b3033bc7ece"
   File?   optTaskParamFile
   Map[String,TaskParams]? tParams
   String?		   checkDataOut
  }

  if (!defined(tParams)) {
    call taskParams.LoadTaskParams as jP {
      input: optional=optTaskParamFile
    }
  }
  Map[String,TaskParams] jParams = select_first([tParams,jP.params])
  #1 cpu, 2h job
  TaskParams tinyJobParam = jParams["xs_short"]
  #2 cpu, 6h job
  TaskParams smallJobParam = jParams["s_medium"]
  #2 cpu, 12h job
  TaskParams longJobParam = jParams["s_long"]
  
  Int tmpMaxSeqLen = num_chunk * chunk_size

  if (!defined(checkDataOut)) {
    call check_repeat_lib {
      input: compgen=compgen, repeatLibFasta=repeat_lib, compgenTag=compgenTag
    }
  }
  call util.genomeGFF as gGFF {
    input: compgen=compgen, genomeFasta=genomeFasta, jParam=jParams["s_short"], compgenTag=compgenTag
  }
  
  call chunk.chunkFasta as chunkGenome {
    input: compgen=compgen, fastaFile=genomeFasta, chunk_size=chunk_size, overlap=overlap,
    jParam=smallJobParam, compgenTag=compgenTag
  }

  call splitFa.splitFasta as splitGenome {
    input: compgen=compgen, fastaFile=chunkGenome.chunkedFasta, job_seq_len=tmpMaxSeqLen,
    jParam=jParams["m_medium"], compgenTag=compgenTag
  }

  Array[Int] genomeChunkCounter = range(length(splitGenome.fastaFiles))
  scatter (j_n_genome_chunk in zip(genomeChunkCounter,splitGenome.fastaFiles)) {
    call repeatMask {
      input: compgen=compgen, chunkFasta=j_n_genome_chunk.right, job_id=j_n_genome_chunk.left, repeat_lib=repeat_lib, rpm_args=rpm_args,
      cpu=longJobParam.cpu, notUsed=select_first([check_repeat_lib.checkDataOut, checkDataOut]), jParam=longJobParam, repeatMaskerTag=repeatMaskerTag
    }
  }

  call util.merge as rptMaskGFF {
    input: files=repeatMask.chunkRptMaskGff, type="gff2"
  }

  call util.merge as rptMaskOut {
    input: files=repeatMask.chunkRptMaskOut, type="out"
  }
  
  call util.merge as rptMaskSeq {
    input: files=repeatMask.chunkMaskSeq, type="fa"
  }

  String genomeFileStem = basename(genomeFasta)
  call joinMaskedSeq {
    input: compgen=compgen, mergedMasked=rptMaskSeq.mergedFile, jParam=tinyJobParam, prefix=genomeFileStem, compgenTag=compgenTag
  }
  
  call transformRptGFF {
    input: compgen=compgen, mergedGFF=rptMaskGFF.mergedFile, chunkedGenomeFasta=chunkGenome.chunkedFasta, mergedOut=rptMaskOut.mergedFile,
    jParam=tinyJobParam, prefix=genomeFileStem, compgenTag=compgenTag
  }

  #16 cpu 6h, num of cpu could be based on genome size to be more efficient
  call util.exonHistogram as histo {
    input: compgen=compgen, genome=gGFF.genomeGFFfile, isGenomeGFF=true, featureGFF=transformRptGFF.GFF_and_i.left, parentType="similarity",
    jParam=jParams["xl_medium"], compgenTag=compgenTag
  }
  
  output {
    File repeatMaskedGenome = joinMaskedSeq.repeatMaskedGenome
    Pair[File,File] genomeRepeatMaskedGFF_and_i = transformRptGFF.GFF_and_i
    File repeatMaskOut = rptMaskOut.mergedFile
    File rpmHistogram = histo.histogram
  }
}

task check_repeat_lib {
  input {
    String      compgen
    File	repeatLibFasta
    String?	compgenTag = ":latest"
  }

  String outLog = "log.check"
  String myOut = "my.stdo.txt"
  
  command {
  set -e
  perl ~{compgen}/compute_farm/saps/tools/check_data4gene_call.pl -ml ~{repeatLibFasta} &> ~{outLog}
  echo "GOOD" > ~{myOut}
  sleep 10
  }

  output {
    String checkDataOut = read_string("~{myOut}")
  }

  #2 cpu is enough
  runtime {
    docker: "doejgi/compgen"  + compgenTag
    cpu:    2
    memory: "8G"
    time:   "00:30:00"
  } 
}
task repeatMask {
  input {
   String compgen
   File   chunkFasta
   File   repeat_lib
   Int	  job_id
   TaskParams jParam
   String? notUsed 
   String? rpm_args = "-xsmall -gccalc"
   String? repeatMaskerTag = ":4.1.2"
   Int?   cpu = jParam.cpu
  }
  
  String outFile = "~{job_id}.chunk.rptMask.gff"
  String outFile2 = "~{job_id}.chunk.rptMask.out"
  String outFile3 = "~{job_id}.chunk.fa"

  #cromwell has new input file for each job by hard link, softlink, etc
  #let make our input (can not be in output dir!) so we can get masked since it is in the same place as input
  command {
    set -e
    rm -rf ~{job_id}Dir
    mkdir ~{job_id}Dir
    ln -s ~{chunkFasta} ~{job_id}Dir/~{outFile3}
    >&2 echo "Starting rptmasker at $(date)"
    ~{compgen}/compute_farm/saps/tools/rpmaskerwrap.pl -seq ~{job_id}Dir/~{outFile3} -rpmp ~{compgen}/RepeatMasker -gff -dir ./ -pa ~{cpu} -keepmasked ~{rpm_args}  -lib ~{repeat_lib} -dotout ~{outFile2} > ~{outFile}
    >&2 echo "Done with rptmasker at $(date)"
  }

  output {
    File chunkRptMaskGff = "~{outFile}"
    File chunkRptMaskOut = "~{outFile2}"
    File chunkMaskSeq = "~{job_id}Dir/~{outFile3}.masked"
  }

  runtime {
    docker: "doejgi/repeatmasker"  + repeatMaskerTag
    cpu:    jParam.cpu
    memory: jParam.ram
    time:   jParam.time    
  }
}

task joinMaskedSeq {
  input {
   String compgen
   File   mergedMasked
   TaskParams jParam
   String     prefix
   String? compgenTag = ":latest"
  }
  
  String outFile = "~{prefix}.repeatMasked.fa"

  command {
    perl ~{compgen}/compute_farm/saps/tools/wrangleFasta.pl -action joinChunked ~{mergedMasked} > ~{outFile}
  }

  runtime {
    docker: "doejgi/compgen"  + compgenTag
    cpu:    jParam.cpu
    memory: jParam.ram
    time:   jParam.time    
  }

  output {
    File repeatMaskedGenome = "~{outFile}"
  }
}

task transformRptGFF {
  input {
   String compgen
   File   mergedGFF
   File   chunkedGenomeFasta
   File   mergedOut
   TaskParams jParam
   String     prefix
   String? compgenTag = ":latest"
  }
  
  String outFile = "~{prefix}.repeatMasked.gff"

  command {
  set -e
  (perl ~{compgen}/data_wrangling/GFFdbIO/utils/transformCoordsGFF.pl ~{mergedGFF} ~{chunkedGenomeFasta} -repeatout4class ~{mergedOut} 2> /dev/null | (grep -v '^#' || true) | sort -s -k1,1 > ~{outFile})
  sleep 10
  perl ~{compgen}/data_wrangling/GFFdbIO/utils/indexGFF.pl ~{outFile}
  }

  runtime {
    docker: "doejgi/compgen"  + compgenTag
    cpu:    jParam.cpu
    memory: jParam.ram
    time:   jParam.time    
  }

  output {
    Pair[File,File] GFF_and_i = ("~{outFile}","~{outFile}.gi")
  }
}

