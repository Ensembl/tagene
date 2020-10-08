
package LoutreWrite::IntronFilter;

use strict;
use warnings;
use LoutreWrite::Default;
use LoutreWrite::Config;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::DB::Fasta;


my $host = `hostname`;
chomp $host;

sub new {
  my $package = shift;
  return bless({}, $package);
}


sub predict_outcome {
    my ($self, $intron, $transcript, $only_novel) = @_;
    my $rel_score_cutoff = -7;
    my $length_cutoff = 50;
    my $source;
    foreach my $att (@{$transcript->get_all_Attributes('remark')}){
        if ($att->value =~ /^Assembled from PacBio CLS reads - (\w+)$/){
            $source = $1;
        }
    }

    my $is_novel = is_novel($intron);
    if ($only_novel and $is_novel eq "no"){
        return 1;  #Skip filtering for existing introns
    }
    my $ss_sequence = get_ss_seq($intron);
    my $antisense_ovlp = get_ss_antisense_overlap($intron);
    my $repeat_ovlp = get_ss_repeat_overlap($intron);
    my $intropolis_score = get_intropolis_score($intron);
    my $rel_score = get_relative_intron_score($intron, $transcript);
    my $intron_length = get_intron_length($intron);
    my $exonerate_ali = get_exonerate_alignment_support($intron, $transcript);
    print "INTRON: ".$intron->seq_region_start."-".$intron->seq_region_end."  NOV=$is_novel SS=$ss_sequence AS=$antisense_ovlp REP=$repeat_ovlp IP=$intropolis_score REL=$rel_score LEN=$intron_length EX=$exonerate_ali\n";

    if ($is_novel eq "yes" and
        ($ss_sequence ne "GT..AG" or
        $antisense_ovlp eq "yes" or
        $repeat_ovlp =~ /SINE|Tandem_repeats/ or
        $intropolis_score == 0 or
        ($rel_score ne "NA" and $rel_score < $rel_score_cutoff) or
        $intron_length < $length_cutoff)){
        return 0;
    }
    else{
        return 1;
    }
}



sub exonerate_support {
    my ($self, $intron, $transcript) = @_;
    my $exonerate_ali = get_exonerate_alignment_support($intron, $transcript);
    if ($exonerate_ali eq "yes"){
        return 1;
    }
    else{
        return 0;
    }
}



# #Connect to the loutre database
# sub get_loutre_db_adaptor {
# #   my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
# #     -host   => 'mysql-ens-havana-prod-1',
# #     -port   => 4581,
# #     -user   => 'ensro',
# #     -pass   => undef,
# #     -dbname => 'loutre_human_tagene_test_6', #this should point to the database being modified - CHANGE THIS!
# #     -driver => 'mysql',
# #   );
#   my $db = $DBA{'otter'};
#   $db->dbc->disconnect_when_inactive(1);
#   return $db;
# }
# 
# # my $l_db = get_loutre_db_adaptor();
# # my $l_sa = $l_db->get_SliceAdaptor;
# 
# 
# 
# #Connect to the Otter pipeline database
# sub get_pipe_db_adaptor { 
#   my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
#     -host   => 'mysql-ens-havana-prod-1',
#     -port   => 4581,
#     -user   => 'ensro',
#     -pass   => undef,
#     -dbname => 'pipe_human',
#     -driver => 'mysql',
#   );
#   $db->dbc->disconnect_when_inactive(1);
#   return $db;
# }
# 
# my $p_db = get_pipe_db_adaptor();
# my $p_sa = $p_db->get_SliceAdaptor;
# 
# 
# #Connect to the Intropolis database
# sub get_intropolis_db_adaptor { 
#   my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
#     -host   => 'mysql-ens-havana-prod-1',
#     -port   => 4581,
#     -user   => 'ensro',
#     -pass   => undef,
#     -dbname => 'gencode_sf5_human_introns',
#     -driver => 'mysql',
#   );
#   $db->dbc->disconnect_when_inactive(1);
#   return $db;
# }
# 
# my $i_db = get_intropolis_db_adaptor();
# my $i_sa = $i_db->get_SliceAdaptor;



sub is_novel {
  my $intron = shift;
  my $chr = $intron->seq_region_name;
  my $sa = $DBA{'otter'}->get_SliceAdaptor();
  my $intron_slice = $sa->fetch_by_region("chromosome", $chr, $intron->seq_region_start, $intron->seq_region_end);
  foreach my $tr (@{$intron_slice->get_all_Transcripts}){
    if ($tr->source =~ /(ensembl|havana)/ and $tr->biotype ne "artifact" and
        scalar(grep {$_->value eq "not for VEGA"} @{$tr->get_all_Attributes('remark')}) == 0){
      foreach my $int (@{$tr->get_all_Introns}){
        if ($intron->seq_region_start == $int->seq_region_start and $intron->seq_region_end == $int->seq_region_end and $intron->seq_region_strand == $int->seq_region_strand){
          return "no";
        }
      }
    }
  }
  return "yes";
}


sub get_ss_seq {
  my $intron = shift;
  my $sa = $DBA{'otter'}->get_SliceAdaptor();
  my $intron_slice = $sa->fetch_by_region("chromosome", $intron->seq_region_name, $intron->seq_region_start, $intron->seq_region_end);
  my $donor = substr($intron_slice->seq, 0, 2);
  my $acceptor = substr($intron_slice->seq, -2);
  my $ssite = "$donor..$acceptor";
  if ($intron->seq_region_strand == -1){
    $ssite = reverse($ssite);
    $ssite =~ tr/ACGTacgt/TGCAtgca/;
  }
  return $ssite;
}


sub get_ss_antisense_overlap {
  my $intron = shift;
  my @as_exon_overlaps;
  my $padding = 10; #Number of nucleotides each side of the splice site that will define the slices 
  my $donor_as = 0;
  my $sa = $DBA{'otter'}->get_SliceAdaptor();
  my $donor_slice = $sa->fetch_by_region("chromosome", $intron->seq_region_name, $intron->seq_region_start-$padding, $intron->seq_region_start+$padding-1);
  if (scalar(grep {$_->seq_region_strand != $intron->seq_region_strand} @{$donor_slice->get_all_Exons})){
    $donor_as = 1;
    #Check that the same splice site does not exist in annotation
    TR:foreach my $tr (grep {$_->seq_region_strand == $intron->seq_region_strand} @{$donor_slice->get_all_Transcripts}){
      next unless $tr->source =~ /(ensembl|havana)/;
      next if scalar(grep {$_->value eq "not for VEGA"} @{$tr->get_all_Attributes('remark')}) == 1;
      next if $tr->biotype eq "artifact";
      foreach my $int (@{$tr->get_all_Introns}){
        if ($int->seq_region_start == $intron->seq_region_start){
          $donor_as = 0;
          last TR;
        }
      }
    }
  }
  my $acceptor_as = 0;
  my $acceptor_slice = $sa->fetch_by_region("chromosome", $intron->seq_region_name, $intron->seq_region_end-$padding+1, $intron->seq_region_end+$padding);
  if (scalar(grep {$_->seq_region_strand != $intron->seq_region_strand} @{$acceptor_slice->get_all_Exons})){
    $acceptor_as = 1;
    #Check that the same splice site does not exist in annotation
    TR:foreach my $tr (grep {$_->seq_region_strand == $intron->seq_region_strand} @{$acceptor_slice->get_all_Transcripts}){
      next unless $tr->source =~ /(ensembl|havana)/;
      next if scalar(grep {$_->value eq "not for VEGA"} @{$tr->get_all_Attributes('remark')}) == 1;
      next if $tr->biotype eq "artifact";
      foreach my $int (@{$tr->get_all_Introns}){
        if ($int->seq_region_end == $intron->seq_region_end){
          $acceptor_as = 0;
          last TR;
        }
      }
    }
  }
  if ($donor_as or $acceptor_as){
    return "yes";
  }
  else{
    return "no";
  }
}


sub get_ss_repeat_overlap {
  my $intron = shift;
  my $sa = $DBA{'loutre'}->get_SliceAdaptor();
  my $intron_slice = $sa->fetch_by_region("chromosome", "chr".$intron->seq_region_name."-38", $intron->seq_region_start-2, $intron->seq_region_end+2);
  #Project slice to clone level as repeat features are stored in clone coordinates
  my @rfs;
  my @projs = @{$intron_slice->project("contig")};
  foreach my $proj (@projs){
    my $clone_slice = $proj->to_Slice();
    push(@rfs, @{$clone_slice->get_all_RepeatFeatures('RepeatMasker')}, @{$clone_slice->get_all_RepeatFeatures('trf')});
  }
  if (scalar(@projs==0)){
    print "NO PROJ\t".$intron->seq_region_name.":".$intron->seq_region_start."-".$intron->seq_region_end."\n";
  }
  my $overlaps_repeat;
  my %intron_repeat_overlaps;
  foreach my $rf (@rfs){
    $rf = $rf->transform("chromosome", "Otter", $intron_slice);
    next unless $rf;
    #overlap with splice site: -2/+2 nt
    if (($rf->seq_region_start <= $intron->seq_region_start + 1 and $rf->seq_region_end >= $intron->seq_region_start - 2) or
        ($rf->seq_region_start <= $intron->seq_region_end + 2 and $rf->seq_region_end >= $intron->seq_region_end - 1)){
      $overlaps_repeat = 1;
      my $repeat_type = get_repeat_type($rf); #pipe_human does not store repeat types, unlike ensembl_core
      $repeat_type =~ s/ /_/g;
      $intron_repeat_overlaps{$repeat_type}++;
    }
  }
  if ($overlaps_repeat){
    return join(",", keys %intron_repeat_overlaps);
  }
  else{
    return "No_overlap";
  }
}


sub get_intropolis_score {
  my $intron = shift;
  my $sa = $DBA{'intron'}->get_SliceAdaptor();
  my $intron_slice = $sa->fetch_by_region("chromosome", $intron->seq_region_name, $intron->seq_region_start, $intron->seq_region_end);
  foreach my $sf (@{$intron_slice->get_all_SimpleFeatures}){
    if ($sf->seq_region_start==$intron->seq_region_start and $sf->seq_region_end==$intron->seq_region_end and $sf->seq_region_strand==$intron->seq_region_strand){
      return $sf->score;
    }
  }
  return "0";
}


sub get_relative_intron_score {
  my ($intron, $tr) = @_;
  my (@scores, $pos, @novel);
  my $n = 0;
  foreach my $int (@{$tr->get_all_Introns}){
    my $score = get_intropolis_score($intron);
    push(@scores, $score);
    my $is_novel = is_novel($intron);
    push(@novel, $is_novel);
    if ($intron->seq_region_start == $int->seq_region_start and $intron->seq_region_end == $int->seq_region_end){
      $pos = $n;
    }
    $n++;
  }
  my $other_score_sum = 0;
  my $other_score_n = 0;
  for (my $j=0; $j<scalar(@scores); $j++){
    next if $novel[$j] eq "yes";
    next if $j==$pos;
    $other_score_sum += $scores[$j];
    $other_score_n++;
  }

  my $log_fold_change;
  my $score = $scores[$pos];
  $score++ if $score == 0;
  if ($other_score_n > 0 and $other_score_sum > 0){
    my $avg_score = $other_score_sum/$other_score_n;
    $log_fold_change = log($score/$avg_score);
  }
  else{
    $log_fold_change = "NA";
  }
  return $log_fold_change;
}


sub get_repeat_type {
  my $rf = shift;
  my $type;
  my $class = $rf->repeat_consensus->repeat_class;
  if ($class =~ /^DNA/){
    $type = "Type II Transposons";
  }
  elsif ($class =~ /^LINE/){
    $type = "Type I Transposons/LINE";
  }
  elsif ($class =~ /^Low_complexity/){
    $type = "Low complexity regions";
  }
  elsif ($class =~ /^LTR/){
    $type = "LTRs";
  }
  elsif ($class =~ /^SINE/){
    $type = "Type I Transposons/SINE";
  }
  elsif ($class =~ /^Satellite/){
    $type = "Satellite repeats";
  }
  elsif ($class =~ /RNA/){
    $type = "RNA repeats";
  }
  elsif ($class eq "trf"){
    $type = "Tandem repeats";
  }
  else{
    $type = "Unknown";
  }
  return $type;
}


sub get_intron_length {
  my $intron = shift;
  return $intron->length;
}



sub get_exonerate_alignment_support {
  my ($intron, $tr) = @_;
  my $read_seq = "";
  foreach my $att (@{$tr->get_all_Attributes('hidden_remark')}){
    if ($att->value =~ /^(pacbio_capture_seq|pacbio_raceseq|SLR-seq|postfilter)_\w+/){ #TO DO: COVER OTHER CASES!
      $read_seq .= get_read_sequences($att->value);
      last;
    }
  }
  my @all_introns;
  
  #generate random file name for the query and target sequences
  my $filename = join("", map{('a'..'z','A'..'Z',0..9)[rand 62]} 0..8);
  print_query_seq($read_seq, $filename);
  print_target_seq($tr, $filename);

      #my $cmd = "exonerate -q query.fa -t target.fa -m est2genome --exhaustive --geneseed 250 -n 1 > z_ex_out";
      #my $cmd = "exonerate -q query.fa -t target.fa -m est2genome --refine region --geneseed 250 -n 1 > z_ex_out";
      #my $cmd = "exonerate -q query.fa -t target.fa -m est2genome --geneseed 250 -n 1 --forcegtag yes > z_ex_out";
  
  my $exonerate = $host eq "havana-06" ? "exonerate" : "/nfs/ensembl/bin/exonerate-2.2.0/exonerate";
  my $pssm_dir = "/homes/jmgonzalez/work/long_read_pipeline/annotation_exercise_results/third_round/exonerate";
  my $dir = "/tmp";
  if (`wc -l $dir/$filename.query.fa | cut -d' ' -f1` < 2){
    return "NA";
  }
  my $cmd = "$exonerate -q $dir/$filename.query.fa -t $dir/$filename.target.fa -m est2genome --splice5 $pssm_dir/s5.pssm --splice3 $pssm_dir/s3.pssm --geneseed 250 -n 1 > $dir/$filename.ali";
  my $err = system($cmd);
  unless ($err){
    my @vulgars = `grep ^vulgar $dir/$filename.ali`;
    chomp @vulgars;
    foreach my $vulgar (@vulgars){
      my $introns = parse_vulgar($vulgar);
      push(@all_introns, @$introns);
    }
  }
  system("rm $dir/$filename.*");
  my %all_introns = map {$_ => 1} @all_introns;
  #print OUT $slice->seq_region_name."\n".$tr->stable_id."\n".join("\n", keys %all_introns)."\n";
  my $e_strand = $intron->seq_region_strand == 1 ? "+" : "-";
  my $intron_str = $intron->seq_region_name.":".$intron->start.":".$intron->end.":".$e_strand."_".get_ss_seq($intron);
  if ($all_introns{$intron_str}){
    return "yes";
  }
  else {
    return "no";
  }
}

sub get_read_sequences {
  my $read_names = shift;
  my $fasta_dir = $host eq "havana-06" ? "/ebi/teams/ensembl/jmgonzalez/lrp" : "/nfs/production/panda/ensembl/havana/lrp";
  my %db = ('SLRseq'  => Bio::DB::Fasta->new("$fasta_dir/SLRseq_merged.fasta"),
            'CLS'     => Bio::DB::Fasta->new("$fasta_dir/Captureseq_merged.fasta"),
            'RACEseq' => Bio::DB::Fasta->new("$fasta_dir/RACEseq_merged.fasta"),
            'PB_test' => Bio::DB::Fasta->new("$fasta_dir/PB_test.fasta"),
           );
  my %fasta_seqs;
  my @warnings;
  #pacbio_raceseq_testis : m160112_015449_42182_c100967152550000001823202805121624_s1_p0_34050_ccs, m160115_185856_42182_c100967052550000001823202805121653_s1_p0_130345_ccs ; pacbio_raceseq_lung : m160322_092041_42182_c100987092550000001823224907191673_s1_p0_45294_ccs ; pacbio_raceseq_liver : m160321_201548_42182_c100987182550000001823224907191657_s1_p0_25946_ccs ;
  my @reads_by_sample = split(/\s?;\s?/, $read_names);
  foreach my $sample (@reads_by_sample){
    my ($sample_name, $sample_read_names) = split(" : ", $sample);
    #if no sample name provided, reassign read names to the right variable
    if ($sample_name and !($sample_read_names)){
      ($sample_name, $sample_read_names) = ("none", $sample_name);
    }
    foreach my $read_name (split(", ", $sample_read_names)){
      if ($read_name =~ /\%3D/){
        $read_name =~ s/\%3D/=/g;
      }
      my $seq;

      if ($sample_name =~ /SLR-seq/){
        $seq = $db{'SLRseq'}->seq($read_name);
      }
      elsif ($sample_name =~ /pacbio_capture_seq/){
        $seq = $db{'CLS'}->seq($read_name);
      }
      elsif ($sample_name =~ /pacbio_raceseq/){
        $seq = $db{'RACEseq'}->seq($read_name);
      }
      elsif ($sample_name =~ /postfilter/){
        $seq = $db{'PB_test'}->seq($read_name);
      }
      elsif ($sample_name eq "none"){
        foreach my $dataset (qw(SLRseq CLS RACEseq)){
          $seq = $db{$dataset}->seq($read_name);
          last if $seq;
        }
      }

      if ($seq){
        $read_name =~ s/=/_/g; #Replace "=" character as it would be transformed to "%3D" by OTF, which would render it unusable by Blixem
        $fasta_seqs{$read_name} = $seq;
      }
      else{
        push(@warnings, "WARNING: no sequence found for read name $read_name in library $sample_name");
      }
    }
  }
  
  my @keys = keys(%fasta_seqs);
  my $multi_seq;
  foreach my $key (@keys){
    $multi_seq .= ">".$key."\n".$fasta_seqs{$key}."\n";
  }
  if (scalar @warnings){
    $multi_seq = join("\n", @warnings)."\n".$multi_seq;
  }
  return $multi_seq;
}


sub print_query_seq {
  my ($seq, $filename) = @_;
  my $dir = "/tmp";
  
  open (SEQ, ">$dir/$filename.query.fa") or die "Can't open $filename.query.fa: $!";
  print SEQ $seq."\n";
  close (SEQ);
}

sub print_target_seq {
  my ($transcript, $filename) = @_;
  #my $sa = $transcript->adaptor->db->get_SliceAdaptor;
  my $dir = "/tmp";
  open (SEQ, ">$dir/$filename.target.fa") or die "Can't open $filename.target.fa: $!";
  my $ext = 0;
  #my $slice = $sa->fetch_by_region("chromosome", $transcript->seq_region_name, $transcript->start - $ext, $transcript->end + $ext);
  my $tslice = $transcript->feature_Slice; 
  my $slice = $tslice->expand($ext, $ext);
  my $seq = $slice->seq;
  if ($transcript->seq_region_strand == -1){
    $seq = reverse($seq);
    $seq =~ tr/ACGTacgt/TGCAtgca/;
  }
  print SEQ ">".$transcript->seq_region_name.":".($transcript->seq_region_start-$ext).":".($transcript->seq_region_end+$ext)."\n".$seq."\n";
  close (SEQ);
}

sub parse_vulgar {
  my $vulgar = shift;
  my @introns;
  #vulgar: m150420_082649_42137_c100793212550000001823172010081515_s1_p0_100150_ccs 120 667 + 86720575-86722026 1448 2 - 2653 M 80 80 G 0 1 M 65 65 G 1 0 M 276 276 3 0 2 I 0 895 5 0 2 M 125 125
  $vulgar =~ s/^vulgar: //; #print OUT "\n".$vulgar."\n";
  my ($q_id, $q_start, $q_end, $q_strand, $t_id, $t_start, $t_end, $t_strand, $score, @elem) = split(/\s+/, $vulgar);
  my ($chr, $gnm_start, $gnm_end) = split(/:/, $t_id); #Genomic start and end
  my $target_coord = ($t_strand eq "+" ) ? $gnm_start + $t_start : $gnm_end - $t_start + 1;
  my $intron_coord;
  my ($has_5, $has_3);
  while (my ($type, $q_size, $t_size) = splice(@elem, 0, 3)){
    #We expect a VULGAR string like "5 0 2 I 0 X 3 0 2" or "3 0 2 I 0 X 5 0 2", usually the former for the plus strand and the latter for the minus strand but, strangely, not always.
    if ($type eq "5" and !$has_3 or
        $type eq "3" and !$has_5){
      my $intron_start = $target_coord;
      $intron_coord = "$chr:$intron_start:";
      $has_5 = 1;
      $has_3 = 1;
    }
    elsif ($type eq "3" and $has_5 or
           $type eq "5" and $has_3){
      my $intron_end = $target_coord + $t_size - 1;
      $intron_coord .= "$intron_end:$t_strand";
      $intron_coord = add_ss_seq($intron_coord);
      #push(@introns, "$intron_start:$intron_end");
      push(@introns, $intron_coord);
      $intron_coord = "";
      $has_5 = 0;
      $has_3 = 0;
    }
    $target_coord += $t_size;
  }
  return \@introns;
}


sub add_ss_seq {
  my $coord = shift;
  my ($chr, $start, $end, $strand) = split(/:/, $coord);
  my $sa = $DBA{'otter'}->get_SliceAdaptor();
  my $intron_slice = $sa->fetch_by_region("chromosome", $chr, $start, $end);
  my $donor = substr($intron_slice->seq, 0, 2);
  my $acceptor = substr($intron_slice->seq, -2);
  my $ssite = "$donor..$acceptor";
  if ($strand eq "-"){
    $ssite = reverse($ssite);
    $ssite =~ tr/ACGTacgt/TGCAtgca/;
  }
  return $coord."_".$ssite;
}



1;



