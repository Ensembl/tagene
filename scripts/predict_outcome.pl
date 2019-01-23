#!/usr/bin/env perl -w

#This script predicts the annotation outcome
#for a set of comp_pipe models in a database
#using the following filtering steps at the intron level:
# splice site sequence
# splice site overlap with antisense exon
# splice site overlap with repeat elements
# Intropolis intron score
# relative intron score
# intron length
#At the transcript level, a model is rejected if at least 
#one of its introns is rejected.


use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
$|=1;

my $tlist;
my $outfile;
my $dbhost;
my $dbport;
my $dbuser;
my $dbpass;
my $dbname;
my $analysis_name;
my $biotype;

&GetOptions(
            'tids=s'     => \$tlist,
            'out=s'      => \$outfile,
            'host=s'     => \$dbhost,
            'port=s'     => \$dbport,
            'user=s'     => \$dbuser,
            'pass=s'     => \$dbpass,
            'dbname=s'   => \$dbname,
            'analysis=s' => \$analysis_name,
            'biotype=s'  => \$biotype,
           );


#Connect to the query database
my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host   => $dbhost,
    -port   => $dbport,
    -user   => $dbuser,
    -pass   => $dbpass,
    -dbname => $dbname,
    -driver => 'mysql',
);
my $dbc = $db->dbc;
my $sa = $db->get_SliceAdaptor;


#Connect to the loutre database
my $l_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host   => 'mysql-ens-havana-prod-1',
    -port   => 4581,
    -user   => 'ensro',
    -pass   => undef,
    -dbname => 'loutre_human',
    -driver => 'mysql',
);
my $l_dbc = $l_db->dbc;
my $l_sa = $l_db->get_SliceAdaptor;


#Connect to the database that stores the long read transcripts
my $lr_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host   => 'mysql-ens-havana-prod-1',
    -port   => 4581,
    -user   => 'ensro',
    -pass   => undef,
    -dbname => 'gencode_long_read_pipeline',
    -driver => 'mysql',
);
my $lr_sa = $lr_db->get_SliceAdaptor;


#Connect to the Ensembl core database
my $e_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host   => 'mysql-ensembl-mirror',
    -port   => 4240,
    -user   => 'ensro',
    -pass   => undef,
    -dbname => 'homo_sapiens_core_95_38',
    -driver => 'mysql',
);
my $e_sa = $e_db->get_SliceAdaptor;


#Connect to the Otter pipeline database
my $p_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host   => 'mysql-ens-havana-prod-1',
    -port   => 4581,
    -user   => 'ensro',
    -pass   => undef,
    -dbname => 'pipe_human',
    -driver => 'mysql',
);
my $p_sa = $p_db->get_SliceAdaptor;


#Connect to the Intropolis database
my $i_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host   => 'mysql-ens-havana-prod-1',
    -port   => 4581,
    -user   => 'ensro',
    -pass   => undef,
    -dbname => 'gencode_sf5_human_introns',
    -driver => 'mysql',
);
my $i_sa = $i_db->get_SliceAdaptor;


my $rel_score_cutoff = -7;
my $length_cutoff = 50;

my %tids;
if ($tlist){
  open (LIST, $tlist) or die "Can't open $tlist: $!";
  while (<LIST>){
    chomp;
    $tids{$_} = 1;
  }
  close (LIST);
}


open (OUT, ">$outfile") or die "Can't open $outfile: $!";
print OUT "coords\ttranscript_id\tgene_name\tis_novel\tsplice_site\tss_antisense\trepeat_overlap\tscore\trel_int_score\tlength\tbroad_biotype\toutcome\treal_outcome\n";

my %all_introns;
my %rejected_introns;
my %rejected_transcripts;
my %rej_reasons;
my $tr_c = 0;
my %real_outcome;

#Fetch comp_pipe models
foreach my $slice (@{$sa->fetch_all("chromosome")}){
  foreach my $gene (@{$slice->get_all_Genes($analysis_name, undef, 1)}){
    my $gene_name;
    if (scalar(@{$gene->get_all_Attributes('name')})){
      $gene_name = $gene->get_all_Attributes('name')->[0]->value;
    }
    foreach my $tr (@{$gene->get_all_Transcripts}){
      if ($analysis_name){
        next unless $tr->analysis->logic_name eq $analysis_name;
      }
      if ($biotype){
        next unless $tr->biotype eq $biotype;
      }
      if ($tlist){
        next unless $tids{$tr->stable_id};
      }
      $tr_c++;
      foreach my $intron (@{$tr->get_all_Introns}){
        #Introns may be present in several transcripts - use coordinates to avoid redundancy
        my $intron_coord = $intron->seq_region_name.":".$intron->seq_region_start."-".$intron->seq_region_end;
        $all_introns{$intron_coord}{'is_novel'} = is_novel($intron);
        $all_introns{$intron_coord}{'ss_sequence'} = get_ss_seq($intron);
        $all_introns{$intron_coord}{'antisense_ovlp'} = get_ss_antisense_overlap($intron);
        $all_introns{$intron_coord}{'repeat_ovlp'} = get_ss_repeat_overlap($intron);
        $all_introns{$intron_coord}{'intropolis_score'} = get_intropolis_score($intron);
        $all_introns{$intron_coord}{'length'} = get_intron_length($intron);
        $all_introns{$intron_coord}{'rel_score'} = get_relative_intron_score($intron, $tr);
        $all_introns{$intron_coord}{'inferred_outcome'} = infer_intron_outcome($intron, $tr);

        if ($all_introns{$intron_coord}{'is_novel'} eq "yes" and
           ($all_introns{$intron_coord}{'ss_sequence'} ne "GT..AG" or
            $all_introns{$intron_coord}{'antisense_ovlp'} eq "yes" or
            $all_introns{$intron_coord}{'repeat_ovlp'} =~ /SINE|Tandem repeats/ or
            $all_introns{$intron_coord}{'intropolis_score'} == 0 or
            ($all_introns{$intron_coord}{'rel_score'} ne "NA" and $all_introns{$intron_coord}{'rel_score'} < $rel_score_cutoff) or
            $all_introns{$intron_coord}{'length'} < $length_cutoff)){
          $rejected_introns{$intron_coord} = 1;
          $rejected_transcripts{$tr->stable_id} = 1;
        }
        print OUT $intron_coord."\t".
                  $tr->stable_id."\t".
                  $gene_name."\t".
                  $all_introns{$intron_coord}{'is_novel'}."\t".
                  $all_introns{$intron_coord}{'ss_sequence'}."\t".
                  $all_introns{$intron_coord}{'antisense_ovlp'}."\t".
                  $all_introns{$intron_coord}{'repeat_ovlp'}."\t".
                  $all_introns{$intron_coord}{'intropolis_score'}."\t".
                  $all_introns{$intron_coord}{'rel_score'}."\t".
                  $all_introns{$intron_coord}{'length'}."\t".
                  broad_biotype($gene->biotype)."\t".
                  ($rejected_introns{$intron_coord} ? "rejected" : "accepted")."\t".
                  $all_introns{$intron_coord}{'inferred_outcome'}."\n";
      }
      #Infer transcript outcome if database is loutre
      my $outcome = "unknown";
      if ($dbname =~ /loutre/){
        if ($tr->biotype eq "comp_pipe"){
          #print $tr->stable_id."\n";
          if (scalar grep {$_->value eq "comp_pipe_rejected"} @{$tr->get_all_Attributes('remark')} or
              scalar grep {$_->value eq "comp_pipe_rejected"} @{$tr->get_all_Attributes('hidden_remark')}){
            $outcome = "rejected";
          }
        }
        else{
          $outcome = "accepted"; #although perhaps modified - CHECK THIS
        }
      }
      $real_outcome{$outcome}++;
      #my $source = map {$_->value} grep {$_->value =~ /^Assembled/} @{$tr->get_all_Attributes('remark')};
      my $source;
      foreach my $att (@{$tr->get_all_Attributes('remark')}){
        if ($att->value =~ /^Assembled/){
          $source = $att->value;
        }
      }
      print "\nTR\t".$source."\t".$outcome."\t".$tr->biotype."\n";
    }
  }
}

my $novel_intron_c = 0;
foreach my $ic (keys %all_introns){
  if ($all_introns{$ic}{'is_novel'} eq "yes"){
    $novel_intron_c++;
    if ($all_introns{$ic}{'ss_sequence'} ne "GT..AG"){
      $rej_reasons{'non-GT-AG'}++;
    }
    if ($all_introns{$ic}{'antisense_ovlp'} eq "yes"){
      $rej_reasons{'antisense_ovlp'}++;
    }
    if ($all_introns{$ic}{'repeat_ovlp'} =~ /SINE|Tandem repeats/){
      $rej_reasons{'repeat_ovlp'}++;
    }
    if ($all_introns{$ic}{'intropolis_score'} == 0){
      $rej_reasons{'no_intropolis_support'}++;
    }
    if ($all_introns{$ic}{'rel_score'} ne "NA" and $all_introns{$ic}{'rel_score'} < $rel_score_cutoff){
      $rej_reasons{'low_relative_intron_score'}++;
    }
    if ($all_introns{$ic}{'length'} < $length_cutoff){
      $rej_reasons{'short_intron_length'}++;
    }
  }
}

print OUT "\n\n\n";
print OUT "$novel_intron_c introns are novel out of ".scalar(keys %all_introns)."\n\n";
print OUT scalar(keys %rejected_introns)." novel introns would be rejected out of ".$novel_intron_c."\n\n";
print OUT scalar(keys %rejected_transcripts)." transcripts would be rejected out of $tr_c\n\n";

print OUT "Rejection reasons:\n";
foreach my $reason (keys %rej_reasons){
  print OUT $reason."\t".$rej_reasons{$reason}."/".$novel_intron_c."\t".sprintf("%.2f", $rej_reasons{$reason} / $novel_intron_c * 100)."%\n";
}

print OUT "\nReal outcome:\n";
foreach my $out (keys %real_outcome){
  print OUT "$out: ".$real_outcome{$out}." transcripts\n";
}
close (OUT);


####


sub is_novel {
  my $intron = shift; 
  
  my $chr = ($intron->seq_region_name =~ /^chr.+-38$/) ? $intron->seq_region_name : "chr".$intron->seq_region_name."-38";
  #print "$chr\t".$intron->seq_region_start."\t".$intron->seq_region_end."\n";
  my $intron_slice = $l_sa->fetch_by_region("chromosome", $chr, $intron->seq_region_start, $intron->seq_region_end);
          #if ($intron->seq_region_start == 1035981 and $intron->seq_region_end == 1036796){
          #  print "\nA - THIS\n";
          #}
  foreach my $tr (@{$intron_slice->get_all_Transcripts}){
    if ($tr->source eq "havana" and $tr->biotype ne "artifact" and !($tids{$tr->stable_id}) and
        scalar(grep {$_->value eq "not for VEGA"} @{$tr->get_all_Attributes('remark')}) == 0){
        #if ($intron->seq_region_start == 1035981 and $intron->seq_region_end == 1036796){
        #    print "\nB - THIS: ".$tr->stable_id."\n";
        #}
      foreach my $int (@{$tr->get_all_Introns}){
        if ($intron->seq_region_start == $int->seq_region_start and $intron->seq_region_end == $int->seq_region_end and $intron->seq_region_strand == $int->seq_region_strand){
          #if ($intron->seq_region_start == 1035981 and $intron->seq_region_end == 1036796){
          #  print "\nC - THIS: ".$tr->stable_id."\n";
          #}
          return "no";
        }
      }
    }
  }
  return "yes";
}


sub get_ss_seq {
  my $intron = shift;
  my $chr = ($intron->seq_region_name =~ /^chr.+-38$/) ? $intron->seq_region_name : "chr".$intron->seq_region_name."-38";
  my $intron_slice = $l_sa->fetch_by_region("chromosome", $chr, $intron->seq_region_start, $intron->seq_region_end);
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
  my $chr = ($intron->seq_region_name =~ /^chr.+-38$/) ? $intron->seq_region_name : "chr".$intron->seq_region_name."-38";
  my $donor_slice = $l_sa->fetch_by_region("chromosome", $chr, $intron->seq_region_start-$padding, $intron->seq_region_start+$padding-1);
  if (scalar(grep {$_->seq_region_strand != $intron->seq_region_strand} @{$donor_slice->get_all_Exons})){
    $donor_as = 1;
    #Check that the same splice site does not exist in annotation
    TR:foreach my $tr (grep {$_->seq_region_strand == $intron->seq_region_strand} @{$donor_slice->get_all_Transcripts}){
      next unless $tr->source eq "havana";
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
  my $acceptor_slice = $l_sa->fetch_by_region("chromosome", $chr, $intron->seq_region_end-$padding+1, $intron->seq_region_end+$padding);
  if (scalar(grep {$_->seq_region_strand != $intron->seq_region_strand} @{$acceptor_slice->get_all_Exons})){
    $acceptor_as = 1;
    #Check that the same splice site does not exist in annotation
    TR:foreach my $tr (grep {$_->seq_region_strand == $intron->seq_region_strand} @{$acceptor_slice->get_all_Transcripts}){
      next unless $tr->source eq "havana";
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
  my $chr = ($intron->seq_region_name =~ /^chr.+-38$/) ? $intron->seq_region_name : "chr".$intron->seq_region_name."-38";
  my $intron_slice = $p_sa->fetch_by_region("chromosome", $chr, $intron->seq_region_start-2, $intron->seq_region_end+2);
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
      $intron_repeat_overlaps{$repeat_type}++;
    }
  }
  if ($overlaps_repeat){
    return join(",", keys %intron_repeat_overlaps);
  }
  else{
    return "No overlap";
  }
}


sub get_intropolis_score {
  my $intron = shift;
  my $chr = $intron->seq_region_name;
  $chr =~ s/^chr(.+)-38$/$1/;
  my $intron_slice = $i_sa->fetch_by_region("chromosome", $chr, $intron->seq_region_start, $intron->seq_region_end);
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
    my $int_coord = $int->seq_region_name.":".$int->seq_region_start."-".$int->seq_region_end;
    my $score;
    if ($all_introns{$int_coord}{'intropolis_score'}){
      $score = $all_introns{$int_coord}{'intropolis_score'};
    }
    else{
      $score = get_intropolis_score($intron);
      $all_introns{$int_coord}{'intropolis_score'} = $score;
    }
    push(@scores, $score);
    my $is_novel;
    if ($all_introns{$int_coord}{'is_novel'}){
      $is_novel = $all_introns{$int_coord}{'is_novel'};
    }
    else{
      $is_novel = is_novel($intron);
      $all_introns{$int_coord}{'is_novel'} = $is_novel;
    }
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


sub broad_biotype {
  my $biotype = shift;
  if ($biotype eq "comp_pipe"){
    return "novel";
  }
  elsif ($biotype eq "protein_coding"){
    return "coding";
  }
  else{
    return "non-coding";
  }
}


sub gene_name {
  my $tr = shift;
  my $gene = $tr->get_Gene;
  my @gene_attribs = @{$gene->get_all_Attributes('name')};
  if (scalar(@gene_attribs)){
    return $gene_attribs[0]->value;
  }
  return "NA";
}


sub infer_intron_outcome {
  my ($intron, $tr) = @_;
  my $outcome = "unknown";
  
  #Find other transcripts that may contain the intron
  my @otrs;
  my $chr = ($intron->seq_region_name =~ /^chr.+-38$/) ? $intron->seq_region_name : "chr".$intron->seq_region_name."-38";
  my $intron_slice = $l_sa->fetch_by_region("chromosome", $chr, $intron->seq_region_start, $intron->seq_region_end);
  foreach my $otr (grep {$_->stable_id ne $tr->stable_id} @{$intron_slice->get_all_Transcripts}){
    if ($otr->source eq "havana" and $otr->biotype ne "artifact" and
        scalar(grep {$_->value eq "not for VEGA"} @{$otr->get_all_Attributes('remark')}) == 0){
      foreach my $int (@{$otr->get_all_Introns}){
        if ($intron->seq_region_start == $int->seq_region_start and $intron->seq_region_end == $int->seq_region_end and $intron->seq_region_strand == $int->seq_region_strand){
          push(@otrs, $otr);
          last;
        }
      }
    }
  }
  #if comp_pipe transcript was rejected and intron not in other transcript...
  if (scalar grep {$_->value eq "comp_pipe_rejected"} @{$tr->get_all_Attributes('remark')} or 
      scalar grep {$_->value eq "comp_pipe_rejected"} @{$tr->get_all_Attributes('hidden_remark')}){
    if (scalar @otrs == 0){
      $outcome = "rejected";
    }
    else{
      $outcome = "accepted";
      print "INTRON ".$intron->seq_region_name.":".$intron->seq_region_start."-".$intron->seq_region_end." IN TR  ".$tr->stable_id."\n";
    }
  }
  #if comp_pipe transcript was accepted, confirm that it still contains the intron
  if ($tr->biotype ne "comp_pipe"){
    foreach my $int (@{$tr->get_all_Introns}){
      if ($intron->seq_region_start == $int->seq_region_start and $intron->seq_region_end == $int->seq_region_end and $intron->seq_region_strand == $int->seq_region_strand){
        $outcome = "accepted";
        last;
      }
    }
  }

  return $outcome;
}