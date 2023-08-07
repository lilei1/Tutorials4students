{
    farm=>'SLURM',
    ##be sure to have correct params related to farm
    ## cori
#    farm_host=>'cori',
#    whole_node=>'32',
#    disable_hyper_threading=>1,
#    farm_option=>'constraint=haswell qos=genepool account=plant',

    ## perlmutter
    farm_host=>'perlmutter',
    disable_hyper_threading=>1,
    farm_option=>'constraint=cpu qos=regular account=m342',

    igc_ppn=>8,
    monitor_ppn=>2,
    database=>'jpa9',
    work_dir=>'./',
    site=>'jgi',
    #single word, better to use upper case
    supergroupname=>'panbhybgcV3',
    totalcpu=>32,
    #use common spec to specify stuff for all strains/genomes, it is like IGC spec without specific part to individual strain/haplotype
    #two below must be in the migc working dir (work_dir above, not same as each IGC work dir)
    igc_spec=>'common.igc.spec',
    genemodel_filter_spec=>'genefilter.spec',

    ##all option below is to specify individual strain/haplotype, normally to use strain name as key in the hash, each option hash must have all
    ##the same keys (default means any thing other the specific one, blank default value also mean no such value
    ##for some options, use default instead of specifying to speed writing spec

    #to make hiconf proteome for PAV, after 1st round done, turn PAV pipeline on below: makenonselfpanproteomefrom, high_confidence_transcript_file
    #IGC will generate 2 files under each igc work-dir: hiConf.proteome.nr.fa, hiConf.ovlps
    hiconfproteomefasta=>{default=>1},
    #default can be used for all strains other than a specific one(s) (xyz below)
    assembly_version=>{default=>'v3.0'},
    annotation_version=>{
	default=>'v3.1'
	},
    ##to make up annotation_version with assembly version, but overwriten by annotation_version if any
    annotation_version_number=>1,
#    display_version=>{default=>'v3.1'},
    assembly_site=>{default=>'JGI'},
    ##value can be migc-workdir-full-path/*/hiConf.proteome.nr.fa
    #makenonselfpanproteomefrom=>{default=>'/global/cfs/cdirs/plantbox/llei/annotation/hybridum/annotation/migc/*/hiConf.proteome.nr.fa'},
    #high_confidence_transcript_file=>{default=>'./hiConf.ovlps'}, #relative and in each IGC own work_dir

    ##num of chr to make sure pasa has enough part to have max of 1 chr per part, and add comma (,) after the number
    n_chromosome=>10,

    ##intron spec
    intron_len=>{
	default=>10000, 
    },
    est_intron_len=>{
	default=>25000,
    },
    
    ##use proper key to distinquish from other IGC runs, not use just primary/alt as key!
    ##for this and anything that is unique to each strain/haplotype like name2genome (masked_genomic_fasta, repeat_gff, etc) must have all strain/haplotype
    species_name=>{
	Bhyb2_2_2_V3 => 'Bhrachypodium hybridum 2-2-2',
	Bstat7_V3 => 'Bhrachypodium hybridum Bstat7',
	Bstat8_V3 => 'Bhrachypodium hybridum Bstat8',
    },
    #key is each strain/haplotype and value is full path genome fasta or relative to this migc work dir
    name2genome=>{
	Bhyb2_2_2_V3 => '/global/u2/l/llei2019/plantbox/annotation/hybridum/genome/Bhybridums_testAnnotation/Brachypodium_hybridum_var_2-2-2.fasta',
	Bstat7_V3 => '/global/u2/l/llei2019/plantbox/annotation/hybridum/genome/assembliesForLi_May18.23/Brachypodium_hybridum_var_Bsta7.fasta',
	Bstat8_V3 => '/global/u2/l/llei2019/plantbox/annotation/hybridum/genome/assembliesForLi_May18.23/Brachypodium_hybridum_var_Bsta8.fasta',
    },

    taxid=>{default=>1071398},

    jbrowse=>{default=>1},
    ##uncomment below after ready to do inactive transcript
#    prepare_vetting_files=>{default=>1},
    inactive_transcript=>{default=>'./lowQ.defaultParams.plus.rRNA.simple.repetitive.tr.txt'}, #relative and they are in respective IGC work dir
    ##either way works or get from CCS pipeline
    #masked_lib=>{default=>'../Zmays.var.B73/liftover.RefGen_v4/Zea.mays.RepBase.repeats.remove.MuDR-9_ZM.fa'},   
	masked_lib=>{default=>'/global/cfs/cdirs/plantbox/annotation/Bhybridum/v2.0/repeatModeler/results/Bhybridum.denovo.repBase.repeats.fa'},
	
    #if genome is masked, specify it here and repeat_gff, but must have default=>'' for other genome without masked genome
    #masked_genomic_fasta=>{
	#default=>'',
    #},
    #repeat_gff=>{
	#default=>'',
    #},

    ##could make one from TA2est (pertran), CCS2est (CCS pipeline), and/or EST2est
    ##either specify est_fasta or merge from max of 3
    est_fasta=>{
	default=>'',
	Bhyb2_2_2_V3 => '/global/cfs/cdirs/plantbox/annotation/Bhybridum/v2.0/annotv2.1/hiConf.transcripts.fa',
	Bstat7_V3 => '/global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/migc/Bstat7_V3/BhybridumBstat7.illuminaPlusHCest.pertran.TAs.fa',
	Bstat8_V3 => '/global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/migc/Bstat8_V3/BhybridumBstat8.illuminaPlusHCest.pertran.TAs.fa',
    },

    #TA2est=>{
	#default=>'',
	#Brachy_strain_1 => '/global/cfs/cdirs/plantbox/annotation/Bhybridum/v2.0/annotv2.1/hiConf.transcripts.fa',
	#Brachy_strain_2 => '/global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/migc/Bstat7/BhybridumBstat7.illuminaPlusHCest.pertran.TAs.fa',
	#Brachy_strain_3 => 'global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/migc/Bstat8/BhybridumBstat8.illuminaPlusHCest.pertran.TAs.fa',
	
    #},
    #CCS2est=>{
	#default=>'',
    #},
    #EST2est=>{default=>''},
    
    ##below could get it from CCS pipeline or ta pipeline
    #fl_acc=>{default=>''},
    #count from rerun result using larger intron
    #transcript_cluster_count=>{
	#default=>'',
    #},
    ##space-separate list or wild carded files path, default here mean no bam_list
    ##this is to do self train augustus gene prediction
    bam_list=>{
	default=>'',
	Bhyb2_2_2_V3 => '',
	Bstat7_V3 => '/pscratch/sd/l/llei2019/Bstat7/pertran/gcomb/*.bam',
	Bstat8_V3 => '/pscratch/sd/l/llei2019/Bstat8/pertran/gcomb/*.bam',
    },
    rnaseq_contig_gff=>{
	default=>'',
	Bhyb2_2_2_V3 => '',
	Bstat7_V3 => '/global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/migc/Bstat7_V3/PERTRAN/results/pileupreadcontig_614394_jpa9.gff.gz ',
	Bstat8_V3 => '/global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/migc/Bstat8_V3/PERTRAN/results/pileupreadcontig_625340_jpa9.gff.gz ',
    },
    rnaseq_bedgraph=>{
        default=>'',
	Bhyb2_2_2_V3 =>'',
	Bstat7_V3 => '/global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/migc/Bstat7_V3/PERTRAN/results/bedgraph_23735_6143960_jpa9.bg',
	Bstat8_V3 => '/global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/migc/Bstat8_V3/PERTRAN/results/bedgraph_24159_6253420_jpa9.bg',
    },
    rnaseq_bigwig=>{
	default=>'',
	Bhyb2_2_2_V3 => '',
	Bstat7_V3 => '/global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/migc/Bstat7_V3/PERTRAN/results/bigwig_23747_6247160_jpa9.bw',
	Bstat8_V3 => '/global/u2/l/llei2019/plantbox/annotation/hybridum/annotation/migc/Bstat8_V3/PERTRAN/results/bigwig_24197_6333370_jpa9.bw',
    },

    ##uncomment to turn on when having fl (correct/cluster CCS) or expression (pertran)
    ##using single value for all, or try hash with default for most common ones
    ##fl alt filter params when have full length transcriptome like Iso-Seq CCS:
#    filter_fl_unsupported_alt_tr=>1,
    ##default should be 1, increase it when having reasonable amount of Iso-Seq CCS reads
#    num_fl_as_full_length=>2,

    ##full length expression level:
#    min_expression_full_length=>10,

    ##for CSP/JGI genome, we must have GA_meta (atid:NNN is good enough without spid)
    #TA_meta=>{default=>''},
    #GA_meta=>{XYZ=>'spid:1257909 atid:419181',QST=>'spid:1257909 atid:421387',default=>''},
    
}
