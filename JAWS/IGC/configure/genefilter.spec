{
   ##mature chromosome level assembly and annotation will be submitted to genbank, below it fix/filter gap in gene model
   gapthreshold=>0.0000001,
   #for homol total>=1.8,not count Unknown repeat as repeat, more diverged from homology seed, could use lower threshold
   discountunknownrepeat=>1.8,
   ##filter out RNAseq-supported-only gene model which has a good gene model on the opposite strand
#   rnaseqprefix=>'enm', #pertran assembly id starting with enm
    #note all *type keys will be added/overwriten by IGC script
    #IGC will add these 2 if not in
    positivefilter=>{name=>'valid',synonym=>'valid',Valid=>1},
    altspliceselectfilter=>{name=>'valid',synonym=>'valid',Valid=>1},
    #IGC will properly add positivetype/esttype, but positivethreshold and valid_intron_only can be added
    positivethreshold=>0.01,
    valid_intron_only=>1,
    #use one with homology_total for genome where more truncated genes (cscore 1,coverage  0.4) are real like no very good homology seeds
    #addwiththreshold=>{coverage=>0.5,cscore=>0.5,homology_total=>1.4},
    #more homology available, increase a bit from: coverage: 0.5, cscore: 0.5
    addwiththreshold=>{coverage=>0.7,cscore=>0.6},
    andthreshold=>1,
    completeness=>1,
    #these for getting gene fragments (gene without start codon or stop codon or neither)
    #homology_total 1.6 for cases where coverage 0.9 and cscore 0.7 or other combinations but cscore+coverage >= 1.6
    addwiththreshold2=>{coverage=>0.7,cscore=>0.8,homology_total=>1.6},
    andthreshold2=>1,
    #for very fragmented genome assembly
    #addedgegenethreshold=>{distance=>1000,cscore=>0.6,coverage=>0.4,mincdslength=>90},
    #get rid of gene overlapping with repeats, by default gene filter will filter out repeat-ov gene
    negativethreshold=>0.2,
    #to ignore simple repeat when filtering gene by overlapping with repeats, turn flag below on when repeatmasking without -nolow option
    discountsimplerepeat=>1,
    #overwrite repeat filtering (subject to completeness check if specified), see above for homology_total
    addbackthreshold=>{cscore=>0.9,coverage=>0.7,homology_total=>1.6},
    addbacknegativethreshold=>0.4,
    #incr addbacknegativethreshold to 1/0.9/0.8 and let Pfam to find TE genes (and live with a few TE genes) since we have high homology requirement?
    #if these params still did not work when we have repeat over masking, manually get a set in gff3 format from filtered out PASA genes, say complete and homology total 1.8 and overlap with picked regardless strand, and use filtered_gff_to_add option in IGC. However we must redo_load2pac with -testredo and then delete mergecompleteness analysis before we can start IGC to make this work. So do this before inactivate transcript
    ## low quality intergenic gene filtering (only used in completeness analysis, not in genefilter analysis)
    #intergenichomology_total=>1.1,
    #intergenicpercentilecutoff=>0.1, #10 percentile of intergenic distance
    #intergenicsplicedcdslength=>500, #for not good homology support ones
    #intergenicunsplicedcdslength=>900, #for not good homology support ones
    removeinternalstopcodongene=>1,
    mincdslength=>200,
    #attributes and their value to overwrite mincdslength requirement: get short gene supported by homology
    overwritelength=>{cscore=>0.9,coverage=>0.8,homology_total=>1.7,Gcomplete=>1},
    #get shorter gene if EST supports it (number of ESTs and cds coverage by ESTs) but not overlap with repeat
    #could add expression_level in addition to depth: depth and expression_level is ORed, depth is number of ESTs.
    #expression_level is default to 10 for spliced, 20 for unspliced,
    #cdsthreshold is EST overlap threshold with CDS (aka ovthreshold)
    #should set expression_level using PERTRAN depth distribution data when having rnaseq_contig_gff in IGC spec
    #nohomology_expression_level under unsplicedest is CDS expression level and can be used to overwrite 0.3 cscore from attrthreshold4unsplicedest
    splicedest=>{depth=>3,cdsthreshold=>0.5,cdslength=>90,negativethreshold=>0.2},
    unsplicedest=>{depth=>8,cdsthreshold=>0.5,cdslength=>150,negativethreshold=>0.2,nohomology_expression_level=>15},
    attrthreshold4unsplicedest=>{cscore=>0.3},
    ##for promoting pasa assembly into a gene when no gene is called
    #promotedminimumcdsratio=>0.2
    #allow some CDS overlap on opposite strand
    maxcdsoverlap=>0.05,
    #allow a bit more CDS overlap on opposite strand if homology is high, translation start overextension problem
    allowcdsoverlap=>{homology_total=>'1.9',maxcdsoverlap=>'0.12'},
}
