#!/usr/bin/perl -w

#This script performs quality checks of the TAGENE annotation
#Usage: perl post_import_QC.pl -dataset human_test

# Transcripts that could be merged
# NMD transcripts with an ORF < 35 aa
# Fragmented loci


use strict;
use warnings;
use Getopt::Long;
use LoutreWrite::Default qw ( exon_novelty can_be_merged );
use Bio::Otter::Lace::Defaults;
use Bio::Otter::Server::Config;
$| = 1;


my $dataset_name;
my $only_chr;
my $only_gene;
my $outfile;

&GetOptions(
            'dataset=s'  =>  \$dataset_name,
            'chr=s'      =>  \$only_chr,
            'geneid=s'   =>  \$only_gene,
            'outfile=s'  =>  \$outfile,
            );

#Find chromosome name if it has been passed as a job array index
$only_chr ||= $ENV{LSB_JOBINDEX};
if ($only_chr){
  $only_chr = "X" if $only_chr eq "23";
  $only_chr = "Y" if $only_chr eq "24";
}

#Connect to loutre database
#DataSet interacts directly with an otter database
my $dataset = Bio::Otter::Server::Config->SpeciesDat->dataset($dataset_name);
my $otter_dba = $dataset->otter_dba or die "can't get db adaptor\n";
$otter_dba->dbc->reconnect_when_lost(1);
my $sa = $otter_dba->get_SliceAdaptor();
my $ga = $otter_dba->get_GeneAdaptor();
my $ta = $otter_dba->get_TranscriptAdaptor();

open (OUT, ">$outfile.$only_chr") or die "Can't open $outfile.$only_chr: $!";

foreach my $slice (@{$sa->fetch_all("chromosome", "Otter")}){
  next unless $slice->seq_region_name =~ /-38$/;
  if ($only_chr){
    next unless $slice->seq_region_name eq "chr$only_chr-38";
  }
  print $slice->seq_region_name."\n";
  foreach my $gene (@{$slice->get_all_Genes}){
    if ($only_gene){
      next unless $gene->stable_id eq $only_gene;
    }
    if (is_tagene($gene)){
      check_redundant_transcripts($gene);
      check_short_NMD_ORF($gene);
      check_fragmented_locus($gene);
    }
  }
}

close (OUT);


##########


sub is_tagene {
  my $gene = shift;
  foreach my $tr (@{$gene->get_all_Transcripts}){
    if (scalar grep {$_->value eq "TAGENE_transcript"} @{$tr->get_all_Attributes('remark')}){
      return 1;
    }
  }
  return 0;
}



sub check_redundant_transcripts {
  my $gene = shift;
  my @trs = @{$gene->get_all_Transcripts};
  for (my $i=0; $i<scalar(@trs); $i++){
    my $tr = $trs[$i];
    if (scalar grep {$_->value eq "TAGENE_transcript"} @{$tr->get_all_Attributes('remark')}){
      #Check for redundant transcripts (where one could be merged into the other)
      for (my $j=$i+1; $j<scalar(@trs); $j++){
        my $db_tr = $trs[$j];
        next if $db_tr->biotype eq "artifact";
        next if scalar(grep {$_->value eq "not for VEGA"} @{$db_tr->get_all_Attributes('remark')}) and $db_tr->biotype ne "comp_pipe";
        next if scalar(grep {$_->value eq "comp_pipe_rejected"} @{$db_tr->get_all_Attributes('remark')});
        next if scalar @{$db_tr->get_all_Introns} == 0;

        if ((exon_novelty($tr, $db_tr) or exon_novelty($db_tr, $tr)) and can_be_merged($tr, $db_tr)){
          print OUT join("\t", "Redundant_transcript", 
                                 $tr->stable_id, scalar(@{$tr->get_all_Exons}), $tr->biotype,
                                 $db_tr->stable_id, scalar(@{$db_tr->get_all_Exons}), $db_tr->biotype,
                                 $gene->stable_id, $gene->biotype)."\n";
        }
      }
    }
  }
}



sub check_short_NMD_ORF {
  my $gene = shift;
  foreach my $tr (@{$gene->get_all_Transcripts}){
    if (scalar grep {$_->value eq "TAGENE_transcript"} @{$tr->get_all_Attributes('remark')}){
      #Check for NMDs with an ORF < 35 aa
      if ($tr->biotype eq "nonsense_mediated_decay" and $tr->translation->length < 35){
        print OUT join("\t", "Short_NMD_ORF", $tr->stable_id, $tr->translation->length)."\n";
      }
    }
  }
}



sub check_fragmented_locus {
  my $gene = shift;
  my %tmp_gene_coords;
  foreach my $tr (sort {$a->seq_region_start <=> $b->seq_region_start} @{$gene->get_all_Transcripts}){
    next if $tr->biotype eq "artifact";
    next if scalar(grep {$_->value eq "comp_pipe_rejected"} @{$tr->get_all_Attributes('remark')});
    #Check for gaps in the gene (fragmented locus)
    if ($tmp_gene_coords{'end'} and $tr->seq_region_start > $tmp_gene_coords{'end'}){
      print OUT join("\t", "Fragmented_locus", $gene->stable_id)."\n";
    }
    unless ($tmp_gene_coords{'start'} and $tmp_gene_coords{'start'} < $tr->seq_region_start){
      $tmp_gene_coords{'start'} = $tr->seq_region_start;
    }
    unless ($tmp_gene_coords{'end'} and $tmp_gene_coords{'end'} > $tr->seq_region_end){
      $tmp_gene_coords{'end'} = $tr->seq_region_end;
    }
  }
  if ($tmp_gene_coords{'start'} ne $gene->seq_region_start or $tmp_gene_coords{'end'} ne $gene->seq_region_end){
    warn "Gene ".$gene->stable_id." coordinates don't match global transcript coordinates!\n".
         $tmp_gene_coords{'start'}." <> ".$gene->seq_region_start." || ".$tmp_gene_coords{'end'}." <> ".$gene->seq_region_end."\n";
  }
}




