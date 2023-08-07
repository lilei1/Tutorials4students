{
 mail_to=>'lilei@lbl.gov', #use single quote or escape @ (\@)
 ## optionally specify email address where message will be promptly received and read, default to mail_to
 #fail_mail_to=>'',

# genecall_template=>'gene_call_AUG_GFF_file.master.tmpl',
 latesttemplate=>1,

 #R500: v1.2.1, unconventional
 okversion=>1,
 use_common_name4nameprefix=>1,
 #transcript_cluster_count=>'',
 est_identity=>95,
 est_coverage=>50,
 #sibling_est_fasta=>'./mg.related.EST.TC.TA.fa',
 #sibling_est_identity=>80,
 #sibling_est_coverage=>70,

 ##when sense_transcript is true, majority of ESTs/TAs should be sense transcripts
 ##since PASA will not merge single exon alignment overlapping with alignment on the opposite strand even if single exon is from three prime EST
 #R500 has both stranded and unstranded, no CCS, all have ESTs (3'), a good number of CCS clusters, selectively set sense_transcript
# sense_transcript=>0,
 sibling_sense_transcript=>0,

 ##get revcomplement TA ORF prediction for single exon TAs
 ##only useful if transcriptome assemblies from unstranded RNA-seq and no good homology seeds as this will bring noise at the same time
 #revcompsingleexon=>1,

 ##use PASA assembly GFF instead of blat PASA assembly fasta for loci location
 ##use this when have lots of ESTs/RNAseq data
 use_assembly_gff=>1,

 relative_proteome=>'',
 ## good idea to include swiss-prot (look for appropriate latest version)
 nonpac_proteomes=>[qw(/global/dna/projectdirs/plant/data/homology_seeds/Swiss-Prot/uniprot_sprot.eukaryote.2022_04.uni95.fa /global/dna/projectdirs/plant/data/homology_seeds/JGI/Paspalum_vaginatum/v3.1/Pvaginatum.v3.1.prot672.hiConf.pep.fa /global/dna/projectdirs/plant/data/homology_seeds/JGI/Panicum_hallii.FIL/v3.2/Phallii.FIL.v3.2.prot590.hiConf.pep.fa /global/dna/projectdirs/plant/data/homology_seeds/JGI/Zea_mays_NKH8431/v1.2/Zmays.NKH8431.v1.2.prot683.hiconf.pep.fa /global/dna/projectdirs/plant/data/homology_seeds/JGI/Sorghum_bicolor/v5.1/Sbicolor.BTx623.v5.1.prot730.hiConf.pep.fa /global/dna/projectdirs/plant/data/homology_seeds/JGI/Setaria_viridis/v2.1/Sviridis.v2.1.prot500.hiconf.tr.homol.pep.fa /global/dna/projectdirs/plant/data/homology_seeds/JGI/Oryza_sativa/JGI.v7.1b/Osativa.Nipponbare.JGI.v7.1b.GMI.prot535.hiConf.pep.fa /global/cfs/cdirs/plantbox/llei/annotation/hybridum/annotation/3-7-2/v1.0t/annotv1.0t/Lrigidum.Atauschii.nr95.fa /global/cfs/cdirs/plantbox/llei/annotation/hybridum/annotation/3-7-2/v1.0t/annotv1.0t/Jascendens.Aamericanus.Aofficinalis.Acomosus.Macuminata.Eguineensis.nr95.fa /global/dna/projectdirs/plant/data/homology_seeds/JGI/Populus_trichocarpa/v4.1/Ptrichocarpa.v4.1.prot533.hiConf.pep.fa /global/dna/projectdirs/plant/data/homology_seeds/JGI/Aquilegia_coerulea/v3.1/Acoerulea.v3.1.prot322.hiconf.pep.fa)],
#arabi,soybean,sorghum BTx642,rice(GMI),Sviridis,aquilegia,grape
proteome_ids=>[qw(680 690 167 508 457)],
 ##comparators are used to run blastp for cscore and coverage, default to nonpac_proteomes and proteome_ids
# comparators=>[qw(130 120 145)],
est_orf_homology_target=>'/global/dna/projectdirs/plant/data/homology_seeds/JGI/Oryza_sativa/JGI.v7.1b/Osativa.Nipponbare.JGI.v7.1b.GMI.prot535.hiConf.pep.fa',
# intron_len=>40000, #intron for homology (proteomes) alignment, use pertran intron distribution say 99.8% if any
# est_intron_len=>67000, #when not specified, est_intron_len will be estimated from EST fasta from PASA run, use 99.99% from pertran
 wiggle=>2000,
 ##Maize, Human, Arabidopsis, not necessary any more, leave it unchanged for now
 genomescan_matrix=>'Arabidopsis',
 ##be sure to change it to Monocots for monocot genome, Dicots is a bit better for magnoliids, probably for others
 fgenesh_param=>'/global/dna/projectdirs/plant/tools/compgen/FGENESH2011/params/Monocots',
 ##uncomment below for algae, if not, busco is on embryophyta_odb9, which may or may not be a bad thing but score is lower
# busco_profile=>'/global/dna/projectdirs/plant/tools/compgen/busco/eukaryota_odb9/',

 ##external prediction gff as additional gene predicitons
 #external_prediction_gff=>'',

 #example=>'Chr01',
 chr_id_regexp=>'^(Chr\d{1,2})', #binning gff for speed and ram usage, in example, we have max of 91 buckets if min 2 digits id, 10-99 plus 'rest', scaffold_10,scaffold_10* goes to bucket of scaffold_10

 ## extra tracks in gff
 #extra_gff=>[qw()],

 ##when not promote prior, use this to archive prior version gene models in gff2 so we can build a track in Zome release
 ##take alllifted_transcript* gff2 from liftover pipeline
 #priorversion_gff2archive=>'',

 ## promote not-so-good gene prediction overlapped with prior version either gene model on same assembly or aligned CDS/tr on different assembly
 ## it is better to use transcript in GFF2 from liftedover pipeline if on a different assembly and jbrowse will build a track for this
 #priorversion_cds_gff=>[qw()],
 #prior_proteome_id=>0,

 ## to get transcript expression level, use output of pileupreadcontig of pertran run, final merged one (name with pileupreadcontigCombined) if you have multiple runs
 #rnaseq_contig_gff=>'',
 ## in conjunction with prior version, promote if expression >= threshold (default 5, use 0 to overwrite), must have rnaseq_config_gff
 #minexpression=>5,#min expression level for gene model to be promoted based on prior version overlap, get from pertran depth distr 25%?
 ## decent expresion like 25 percentile depth or higher (min depth 15 to be safe?) to say RNAseq transcript assembly of this level or higher is full length like full length cDNA
 #get from ta pipeline min_expression_full_length=>10,

 ## difference b/n extra_* and explicit options is to store type info for explicit ones in tracking of archive
 #rnaseq_bedgraph=>[qw()],
 #rnaseq_bigwig=>[qw()],

 ##correct options to run cori, but must use pasa_work_dir_root under your $CSCRATCH to avoid PASA run timeout (SLOW /projectb on cori),
 ##same for staging_root. Can not switch to denovo since CSCRATCH is not accessible in denovo once starting to run on cori
 staging_root=>'${PSCRATCH}/pstaging',
 pasa_work_dir_root=>'${PSCRATCH}/PASA_RUNS',
 avoid_exclusive=>1,
 hiconfaup_cpu=>64,
 usergroup=>'compgen',
 smtp_server=>'smtp.lbl.gov',
 ##change or add email address, use comma as separator, to send archive (gene set finalize) email
 mail_to_when_finalize=>'lilei@lbl.gov',

 do_pasa_on_farm=>1,

  ##blastx/blastp matrix and parameters
 #blastx_options=>'-F "m S" -U  -b 15000 -v 15000 -K 20 -e 1e-5',
 #blastx_matrix=>'BLOSUM62',
 #blastp_options=>'-F F -C F -b 100 -v 100 -e 1e-3',
 #blastp_matrix=>'BLOSUM45',

 ## these 2 options are for exonerate: --model protein2genome --showtargetgff yes and they can not be changed'
 ## default additiona exonerate options below and intron_len is from est_intron_len option
 ## when overwrite additional default, make sure to use right value (literal value) for maxintron',
 #exonerate_options=>'--softmasktarget yes --maxintron $est_intron_len',

 ## num of peptides of loci in one exonerate job, testing indicate 1000 is about right and job should finish in <12 hours
 #exonerate_job_n_peptide=>1000,

 ## num of peptides allowed a for given locus, use smaller value 40 or even 20 for cases with lots of redundant seeds (lots of hits in one locus)
 #locus_max_peptide=>100,

 ## for non-main-genome translation, we must have map from chrom to BioPerl CodonTable Id, consult web or SOI::Visitor->translate for ID
 ## perldoc /global/dna/projectdirs/plant/tools/compgen/lib/perl/SOI/Visitor.pm
 ## make sure chromosome and its CodonTable Id are valid!
 #codontable_map=>{chrMt=>22,chrCp=>11},

 ## below is good weight when use exonerate
 pep_coverage_weight=>0.5, #blastp
 target_coverage_weight=>0.5, #exonerate
 exon_coverage_weight=>1, #exonerate
 gene_complete_weight=>0.5, #complete score 0.5, 0.33 having stopCodon, 0.17 having startCodon

 ##track look-and-feel hash: {trackDescription=>'bla bla',trackKey=>'Gene locus Expression of experimentName'}
 ##without look-and-feel-hash, IGC will auto fill these keys: trackLabel,trackKey,trackShortDescription,trackDescription
 ##search all fpkm options to see how these can be modified without specifying the hash for each fpkm file
 #gene_fpkm=>'', #value can be hash of key file and value is another hash: gene_fpkm=>{file1=>{key=>value,key2=>value2}},
 #transcript_fpkm=>'', #value can be hash
 ##see help on archive_fpkm when you want to archive fpkm files
 #archive_fpkm=>0,

##for dev only
 ##WARNING: must NOT use gmap-2016-06-09 which is good for Iso-Seq
 ##use new PASA (PASA_v2.0.2), PASApipeline2.0.2a (PASA_v2.0.2a) has JGI intron scoring
 ##PASApipeline-2.0.2 only supports GT/AG, GC/AG, AT/AC
 pasa_root=>'/global/dna/projectdirs/plant/tools/compgen/PASApipeline-2.0.2b',
 ##genepool
# gmap_root=>'/global/dna/projectdirs/plant/tools/compgen/gmap-2016-09-20',
 ## denovo, cori
# gmap_root=>'/global/dna/projectdirs/plant/tools/compgen/gmap-2016-09-20b',

# script_root=>'/global/dna/projectdirs/plant/tools/compgen/dev_pipe',
}
