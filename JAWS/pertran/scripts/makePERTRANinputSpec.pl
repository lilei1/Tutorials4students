eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use FindBin qw($RealBin);
use lib ($INC[0]=~/^\./)?"$RealBin/$INC[0]":"$RealBin/..";

use strict;
use File::Spec;
use FileHandle;
use Cwd 'abs_path';
use Getopt::Long;
use JSON::XS;
use SAPS::Utils;
use SAPS::Manager;

my $argh = {};
my $ok =
  GetOptions ($argh,
	      "genomefasta=s",
	      "read_root=s",
	      "reads=s%",
	      "is_single_end",
	      "has_use_shared_memory",
	      "stranded",
	      "pair_is_1_or_2=s",
	      "intron_length=s",
	      "interleave",
	      "n_read_a_job=i",
	      "variable_length_in_file",
	      "specfile=s",
	      "help",
	     );

sub help_msg {
  my $msg = shift;
  print "$0 to generate input json (to stdout) for WDL from, genome, intron-length read root and tissue specification (hash)\n";
  print "  -spec yourSpecFile which must have genomefasta, intron_length, root_root and reads options,\n";
  print "  other options: is_single_end (default false), interleave (default false), n_read_a_job (default 1M per gsnap job), pair_is_1_or_2 (default true/1),\n";
  print "    no effect if interleave, set to 0 if not 1 or 2)\n";
  print "  -has_use_shared_memory, for gsnap newer than May 2017, gsnap has use-shared-memory option and it is on by default and it is turned off in pertran\n";
  print "  -straned, default is false\n";
  print "    but one node could have jobs doing gsnap for different genome, we need to turn it off\n";
  print "  still others not listed but they are part of argument construction that can be changed in json before run WDL,\n";
}
if ($argh->{specfile}) {
  my $inh = get_hash_in_file($argh->{specfile});
  $argh = {%$inh, %$argh};
  $argh->{intron_length} = $argh->{intron_len} if ($argh->{intron_len});
  $argh->{genomefasta} = $argh->{genome} if ($argh->{genome});
  $argh->{is_single_end} = 0 if ($argh->{paired}); #in case using saps pertran addon which has 'paired' option
}
$argh->{pair_is_1_or_2} = 1 unless (defined($argh->{pair_is_1_or_2}));
if ($argh->{is_single_end}) {
  print STDERR "single end reads should be be interleave nor stranded, they are turned off if any\n" if ($argh->{stranded} || $argh->{interleave});
  delete $argh->{stranded};
  delete $argh->{interleave};
}
if ($argh->{help} || !($argh->{genomefasta} && $argh->{read_root} && $argh->{reads} && $argh->{intron_length})) {help_msg();exit}
add_defaults($argh);

my $inHash = get_input_hash($argh->{read_root},$argh->{reads},$argh);
my $nameArray = get_input_names($inHash);
#check tissue
my @notIn = grep{!exists $inHash->{$_}}keys %{$argh->{reads}};
die "some tissue has no or fail to get any fastq: ".join(", ",@notIn) if (@notIn);
my @extra = grep{!exists $argh->{reads}->{$_}}keys %$inHash;
die "get extra tissue when construct read map for PTR_Assemble WDL: ".join(", ", @extra) if (@extra);

my @fq;
map{push @fq, map{@$_}@{$inHash->{$_}}}keys %$inHash;
my $n_lib = scalar(map{@{$inHash->{$_}}}keys %$inHash);
#my $hash = {"PTR_Assemble.tissueNames"=>$nameArray, "PTR_Assemble.tissue2inputMap"=>$inHash};
#my $utf8_encoded_json_hash_text = encode_json($hash);
#print "$utf8_encoded_json_hash_text\n";
my $is_I_64 = &is_illumina_64($fq[0]);
my $r_len = &get_r_len($fq[0],$argh);
my $i_len = $argh->{intron_length};
#SAPS pertran default withut this:--maxsearch=100
my $align_args = sprintf("--gunzip -n %d -N 1 -w $i_len --novelend-splicedist=$i_len -K $argh->{min_hsp_len} -l $argh->{min_hsp_len} -B 5 %s %s %s", $argh->{n_hit}+1, mismatch_rate_arg($argh->{min_align_id}), $argh->{has_use_share_memory} ? "--use-shared-memory=0 --use-sarray=0" : "", $is_I_64 ? "--quality-protocol=illumina" : "");
#let wdl construct this, pairmax_arg($argh->{intron_length},$r_len,!$argh->{is_single_end})); as read length is per lib
my $parse_args = sprintf("-maxintron $i_len -minimumaligned $argh->{min_align_id} -minimumquality $argh->{minmapq} -minimumbasequality $argh->{minbaseq} -hsplength4realign $argh->{hsp_len_to_realign} -realignminintron $argh->{min_intron_len_to_realign} -gaplength $argh->{maxgaplength} -nhit4repetitive %d %s %s",$argh->{n_hit}+1, $argh->{stranded} ? ($argh->{firstReadStrandIsmRNA} ? "-flipstrand2" : "-flipstrand") : "", $argh->{is_single_end} ? "-single" : "");
#let wdl construct this for assemble, -readlength $r_len, as read length is per lib
#add strandedread ext_n_merge, assemble (readcontig, linkmerge) args
my $asmbl_args = sprintf("-threshold %s -minintrondepth %d -minsplicedintrondepth %d -minmutantsplicedintrondepth %d -minexonlength %d -maxintronlength %d -maxgaplength %d %s -minorminratio %s",$argh->{threshold},$argh->{mingappedintrondepth},$argh->{minsplicedintrondepth},$argh->{minmutantsplicedintrondepth},$argh->{minexonlength},$i_len,$argh->{maxgaplength},$argh->{stranded}?"-strandedread":"",$argh->{minorminratio});
my $minSpliced=$argh->{minsplicedintrondepth_in_all_assembly} ? $argh->{minsplicedintrondepth_in_all_assembly} : int(1.8 * $argh->{minsplicedintrondepth});
my $minGapped = $argh->{mingappedintrondepth_in_all_assembly} ? $argh->{mingappedintrondepth_in_all_assembly} : int(1.8 * $argh->{mingappedintrondepth});
my $asmbl_args4whole = sprintf("-threshold %s -minintrondepth %d -minsplicedintrondepth %d -minmutantsplicedintrondepth %d -minexonlength %d -maxintronlength %d -maxgaplength %d %s -minorminratio %s",$argh->{threshold},$minGapped,$minSpliced,$argh->{minmutantsplicedintrondepth},$argh->{minexonlength},$i_len,$argh->{maxgaplength},$argh->{stranded}?"-strandedread":"",$argh->{minorminratio});
my $ext_n_merge = sprintf("-minintronreadthroughexpressionratio $argh->{min_intron_read_through_expression_ratio} -intronreadthroughmindepth $argh->{min_intron_read_through_depth} $argh->{want_readthrough_arg} %s -minorminratio $argh->{minorminratio}",$argh->{stranded}?"-strandedread":"");
$align_args =~ s/\s+/ /g;
$parse_args =~ s/\s+/ /g;
$asmbl_args =~ s/\s+/ /g;
$asmbl_args4whole =~ s/\s+/ /g;
map{s/\s+$//}($align_args,$parse_args,$asmbl_args,$asmbl_args4whole);
print  "{\n";
print  "   \"PTR_Assemble.compgen\": \"/usr/local/compgen\",\n";
printf "   \"PTR_Assemble.isIllumina64\": %s,\n",$is_I_64 ? "true" : "false";
printf "   \"PTR_Assemble.variable_length_in_file\": true,\n" if ($argh->{variable_length_in_file});
printf "   \"PTR_Assemble.intronLength\": %d,\n", $i_len;
printf "   \"PTR_Assemble.filterBaseQuality\": %d,\n",$argh->{minbaseq};
printf "   \"PTR_Assemble.genomeFasta\": \"%s\",\n",full_path($argh->{genomefasta});
printf "   \"PTR_Assemble.isSingleEnd\": %s,\n",$argh->{is_single_end} ? "true" : "false";
printf "   \"PTR_Assemble.isStranded\": %s,\n",$argh->{stranded} ? "true" : "false";
printf "   \"PTR_Assemble.gsnap_args\": \"%s\",\n",$align_args;
printf "   \"PTR_Assemble.parse_args\": \"%s\",\n",$parse_args;
printf "   \"PTR_Assemble.assemble_args\": \"%s\",\n",$asmbl_args;
printf "   \"PTR_Assemble.assemble_args4whole\": \"%s\",\n",$asmbl_args4whole;
printf "   \"PTR_Assemble.extend_n_merge_args\": \"%s\",\n",$ext_n_merge;
printf "   \"PTR_Assemble.n_read_a_job\": \"%d\",\n",$argh->{n_read_a_job} || 1000000;
printf "   \"PTR_Assemble.n_library\": %d,\n",$n_lib;
print  "   \"PTR_Assemble.tissueNames\": ".encode_json($nameArray).",\n";
print  "   \"PTR_Assemble.tissue2inputMap\": ".encode_json($inHash)."\n";
print  "}\n";
print STDERR "Done\n";
exit;

sub get_hash_in_file {
  my $file = shift;

  open(F, $file) || die "can not open spec file: $file: $!";
  my $perlstr = join("", <F>);
  close(F);
  my $inh = eval $perlstr;
  if ($@) {
    die "Error in _get_hash_in_file: $@";
  }
  return $inh;
}

sub full_path {
  my $file = shift;
  die "file, $file, does not exist or readable" unless (-f $file && -r $file);
  return abs_path($file);
}
sub get_input_hash {
  my ($root,$fhash,$argh)=@_;
  my $paired = not $argh->{is_single_end};
  my $interleave = $argh->{interleave};
  my $pair_is_1_or_2 = $argh->{pair_is_1_or_2};
  $root=~s/\/$//;
  my @list;
  my @k=grep{$_}(keys %$fhash);
  my @v=grep{$_}(values %$fhash);
  die "\nreads value must be a hash (key value pairs)" unless (@k == @v);
  my $limit = 10;
  die "\nCan not run more than $limit samples/tissues. Please reorganize data with num of keys in reads hash <= $limit. You can also use 2 level key hash to organize data: first level key num <= $limit and for each key, you can have any number of keys underneath.\nIdeally keys for arrayref of read files are to partition read files in a reasonable biology fashion for differential expression compute, and second level keys are to group read files in a such a way that they have similar alternative transcription for a given first level key.\n" if (@k>$limit);
  if (ref($fhash->{$k[0]}) eq "ARRAY") {
    map{push @list,@{$fhash->{$_} || []}}keys %{$fhash}
  } else {
    die "does not support 2 level hash spec for reads and value must be array ref like tissue1=>[qw(xxx)],";
    #map{my $k=$_;map{push @list,@{$fhash->{$k}->{$_} || []}}keys %{$fhash->{$k}}}@k; #2 level hash spec
  }
  my $list=join(" ",@list);
  my @rlist;
  for (split/\s+/,$list) {
    if ($_=~/(\*+|\[\w+-\w+\]|\[\d+\]|\?)/) {
      my $wfile=/^\//?$_:"$root/$_";
      push @rlist,split/\n/,`ls $wfile`;
    } else {
      push @rlist,/^\//?$_:"$root/$_";
    }
  }
  map{die("file, $_, does not exist or not readable") unless (-f $_ && -r $_)}@rlist;
  my %uniqf;
  map{my @n=split/\//;$uniqf{$n[-1]}++}@rlist;
  my @dup=grep{$uniqf{$_}>1}keys %uniqf;
  die("Read file name must be unique and can NOT be in multiple tissue/condition groups which could be sub word in non-wild card in reads spec\nhere are duplicated: ".join(",",@dup)) if (@dup);
  return "" unless (length($list)>0);
  my %inHash;
  @rlist=sort{$a cmp $b}@rlist;
  if ($paired && not $interleave) {
    my ($pairedHash,@unpaired)=SAPS::Utils->pairing_pe_read_files(@rlist,$pair_is_1_or_2);
    %inHash=%{$pairedHash};
    die("file name for paired end mates in 2 fils must have only 1 letter different and they must be 1 or 2, here are files finding no pairing:\n".join("\n",@unpaired)) if (@unpaired || scalar(keys %inHash) != @rlist / 2)
  } else {
    map{my @n=split/\//,$_;push @{$inHash{$n[-1]}},$_}@rlist
  }
  die("input file name must not have space") if (grep{/\s+/}keys %inHash);
  die("reads spec keys must have not single quote ror double quote nor space") if (grep{/['" ]/}@k);
  #above is to validate
  #now do one tissue at a time
  my %in_hash;
  for my $k (@k) {
    my @list2;
    if (ref($fhash->{$k}) eq "ARRAY") {
      push @list2,@{$fhash->{$k}};
    }
    my $list = join(" ", @list2);
    my @rlist;
    for (split/\s+/,$list) {
      if ($_=~/(\*+|\[\w+-\w+\]|\[\d+\]|\?)/) {
	my $wfile=/^\//?$_:"$root/$_";
	push @rlist,split/\n/,`ls $wfile`;
      } else {
	push @rlist,/^\//?$_:"$root/$_";
      }
    }
    if ($paired && not $interleave) {
      my ($pairedHash,@unpaired)=SAPS::Utils->pairing_pe_read_files(@rlist,$pair_is_1_or_2);
      push @{$in_hash{$k}}, map{$_}values %{$pairedHash};
    } else {
      push @{$in_hash{$k}}, map{[$_]}@rlist;
    }
  }
  return \%in_hash;
}
sub get_input_names {
  my @all=@_;
  my $hashStr;
  my $hash;
  if (ref($all[0])) {
    $hash=$all[0]
  } else {
    $hashStr=join(", ",@all);
    $hash= eval $hashStr;
    die("$@\nhere is hash string\n$hashStr\n") if ($@)
  }
  die("KEYS: ".join(",",keys %{$hash})."$hashStr") unless (keys %{$hash});
  printf STDERR "n input fastq groups/tissues: %s\n",scalar(keys %{$hash});
  return [keys %{$hash}];
}

sub is_illumina_64 {
  my ($rf, $argh) = @_;
  if (ref($rf)) { die "first arg must be scalar (file)"}
  my $unc;
  if ($rf =~ /\.(gz|gzip)$/) { $unc = "gzip -dc"}
  if ($rf =~ /\.bz2$/) { $unc = "bzip2 -dc";}
  return if ($argh->{q_is_phred33});
  my $fh;
  if ($unc) {
    $fh = FileHandle->new("$unc $rf|") or die("Error making read file handle with $unc $rf|: $!")
  } else {
    $fh = FileHandle->new("$rf") or die("Error open $rf to read: $!")
  }
  die("empty fastq? as it is at eof") if ($fh->eof);
  while (1) {
    my $l1=<$fh>;my $l2=<$fh>;my $l3=<$fh>;my $q=<$fh>;
    die ("did not find quality line from $rf (compress command: $unc)\nhere is first 3 lines\n".join("\n","l1:$l1","l2:$l2","l3:$l3","l4:$q",<$fh>,$fh->eof)) unless ($q=~/\S+/);
    die("$rf is not fastq format, here is first line:$l1") unless ($l1=~/^\@/);
    if ($q=~/[!\"#$%&\'()*+,-.\/0123456789:]/) {return }
    elsif ($q=~/[KLMNOPQRSVWXYZ\[\]^_abcdefgh]/) {return "--quality-protocol=illumina"}
  }
  return;
}

sub add_defaults {
  my $argh = shift;

  my $defaults =
    {
     stranded=>0,
     #this should be default unless it is NOT dUPT protocol,take mate2 strand as mRNA's, in other words, mate is first strand of cDNA
     firstReadStrandIsmRNA=>0,

     #for transcriptome assembly
     #when paire-end: inner distance b/n 2 reads to infer intron when it is > 0 (larger than maxgaplenth implying intron), use max instead of avg
     #must not change this to non-0 as inferred intron is not accurate for modern long reads
     maxgaplength=>0,
     #for making contig
     threshold=>1.5,
     #spliced + gap implied: assembly option name: minintrondepth, meaningless when maxgaplength set to 0
     mingappedintrondepth=>5,

     #increase it for better specificity if lots of reads per sample/tissue, or reduce 2 below
     minsplicedintrondepth=>1,
     minmutantsplicedintrondepth=>5,

     minexonlength=>18, #still can get smaller exon like 6 bp exon if spliced aligned reads high, 6 (old microexon: 3bp, 10 spliced)
     #when assembly by tissue/sample, need to have higher threshold
     minsplicedintrondepth_in_all_assembly=>0, #use default: 1.8 * minsplicedintrondepth
     mingappedintrondepth_in_all_assembly=>0, #use default: 1.8 * mingappedintronthreshold
     #maxintronlength=>15000,#made from intron_len

     #extend and merge
     minorminratio=>0.1,
     min_intron_read_through_expression_ratio=>0.9,
     min_intron_read_through_depth=>10,
     want_readthrough_arg=>'-wantreadthrough', #blank for not to have intron read-thru

     variable_length_in_file=>0, #for cases where some trimming done per fastq file
     #compiled gsnap capability
     gsnap_handle_length=>250,
     max_read_length=>251, #some produce a bit longer than the handled (250), we trim anything longer than the handled, this is to check if we trim too much

     #for alignment
     n_read_a_job=>2000000,
     n_hit=>5,
     min_align_id=>0.95,
     minmapq=>1, #alignment quality, use n_hit instead as multi hit get very poor map score
     minbaseq=>17, #for map filtering when calculating good base length, 17=2% error, 20=1% error
     min_hsp_len=>18,
     #more_gsnap_options=>'',#should only for ones not in gsnap analysis args, e.g. --clip-overlap -m 0.05

     #parsing
     #non 0 to trigger realignment to possible reduce intron size
     #small hsp that could be coming from neighbouring tandem duplication and is short
     #looks for exact match at shorter distance
     #bigger hsp_len_to_realign and/or smaller min_intron_len_to_realign triggers more realignment
     #in expense of more compute cycle/time
     hsp_len_to_realign=>20,
     min_intron_len_to_realign=>800, #no effect unless hsp_len_to_realign
    };
  for my $k (keys %$defaults) {
    $argh->{$k} = $defaults->{$k} unless (defined($argh->{$k}));
  }

  return $argh;
}

sub merge_spec_arg {
  my ($rlength,$argh)=@_;
  return sprintf("-threshold %s -minintrondepth %d -minsplicedintrondepth %d -minmutantsplicedintrondepth %d -minexonlength %d -maxintronlength %d -maxgaplength %d -readlength %d",$argh->{threshold},$argh->{mingappedintrondepth},$argh->{minsplicedintrondepth},$argh->{minmutantsplicedintrondepth},$argh->{minexonlength},$argh->{intron_length},$argh->{maxgaplength},$rlength);
}
sub merge_spec_arg4whole {
  my ($rlength,$argh)=@_;
  my $minSpliced=$argh->{minsplicedintrondepth_in_all_assembly} ? $argh->{minsplicedintrondepth_in_all_assembly} : int(1.8 * $argh->{minsplicedintrondepth});
  my $minGapped = $argh->{mingappedintrondepth_in_all_assembly} ? $argh->{mingappedintrondepth_in_all_assembly} : int(1.8 * $argh->{mingappedintrondepth});
  return sprintf("-threshold %s -minintrondepth %d -minsplicedintrondepth %d -minmutantsplicedintrondepth %d -minexonlength %d -maxintronlength %d -maxgaplength %d -readlength %d",$argh->{threshold},$minGapped,$minSpliced,$argh->{minmutantsplicedintrondepth},$argh->{minexonlength},$argh->{intron_length},$argh->{maxgaplength},$rlength);
}
#tighter mismatch criteria?
sub mismatch_rate_arg {
  my $id=shift;
  my $mr=sprintf("%1.3f",1-$id-0.001);$mr=~s/0+$//;
  die("min align id must be a fraction, you have $id") unless ($id <=1 && $id=>0.8);
  return "--max-mismatches=$mr"
}
sub get_r_len {
  my ($rf,$argh)=@_;
  my $unc;
  if ($rf =~ /\.(gz|gzip)$/) { $unc = "gzip -dc"}
  if ($rf =~ /\.bz2$/) { $unc = "bzip2 -dc";}
  my $can_len=$argh->{gsnap_handle_length};
  my $max_len=$argh->{max_read_length};
  my $variableL=$argh->{variable_length_in_file};
  my $paired=not $argh->{is_single_end};
  my $gapL=$argh->{maxgaplength};
  if($unc){open R,"$unc $rf|" or die("Error making read file handle with $unc $rf|: $!")}else{open R,"$rf" or confess("Error open $rf to read: $!")};
  my $rLen;
  if ($variableL && $paired && $gapL){
    my $c=0;
    my @lens;
    while(<R>){my $r=<R>;<R>,<R>;chomp $r;push @lens,length($r);last if ($c++>1000)};
    @lens=sort{$b<=>$a}@lens;$rLen=$lens[0]
  } else {
    <R>;
    my $r=<R>;chomp $r;$rLen=length($r)
  }
  close R;
  die("compiled gsnap can only handle up to $can_len, need to recompile for longer length, this one length is ".$rLen) if ($rLen > $max_len);
  return $rLen;
}
#avg from sampling first 10K if variable
sub get_avg_r_len {
  my ($rf,$argh)=@_;
  my $unc;
  if ($rf =~ /\.(gz|gzip)$/) { $unc = "gzip -dc"}
  if ($rf =~ /\.bz2$/) { $unc = "bzip2 -dc";}
  my $can_len=$argh->{gsnap_handle_length};
  my $variableL=$argh->{variable_length_in_file};
  my $paired=not $argh->{is_single_end};

  if($unc){open R,"$unc $rf|" or confess("Error making read file handle with $unc $rf|: $!")}else{open R,"$rf" or confess("Error open $rf to read: $!")};
  my $rLen;
  if ($variableL){
    my $c=0;
    my @lens;
    my $t=0;
    while(<R>) {
      my $r=<R>;<R>,<R>;chomp $r;push @lens,length($r);
      if ($c++>10000) {
	map{$t += $_}@lens;$rLen=int($t/scalar(@lens));
	last;
      }
    }
  } else {
    <R>;my $r=<R>;chomp $r;$rLen=length($r)
  }
  close R;
  return $rLen
}

sub pairmax_arg {
  my ($intronl,$readl,$paired)=@_;
  return "" unless ($paired);
  my $l=$intronl+4*$readl;
  return "--pairmax-rna=$l"
}
#each align handle this in WDL
sub sam_out_opt_arg {
  my ($in,$root,$fhash,$PL)=@_;
  my (@id)=split/\//,$in;
  my $ID=$id[-1];$ID =~ s/\.(gz|bz1)$//;$ID =~ s/\.(fastq|fq)$//;
  my (%sm_h,%lb_h);
  map{
    my $sm=$_;
    if(ref($fhash->{$sm}) eq "ARRAY"){
      map{
	my $fn="$root/$_";
	map{my $fpath=$_;$fpath =~ s/\/+/\//g;$sm_h{$fpath}=$sm}split/\n/,`ls $fn`
      }@{$fhash->{$sm}}
    }else{
      map{
	my $k=$_;
	map{my $fn="$root/$_";
	    map{my $fpath=$_;$fpath =~ s/\/+/\//g;$sm_h{$fpath}=$sm}split/\n/,`ls $fn`
	  }@{$fhash->{$sm}->{$k}}
	}keys %{$fhash->{$sm}}}
  }keys %$fhash;$in =~ s/\/+/\//g;
  my $SM=$sm_h{$in};
  die ("did not get tissue name") unless ($SM);
  my $arg = "--read-group-id=$ID --read-group-name=$SM --read-group-library=$SM --read-group-platform=$PL";
  return $arg
}
