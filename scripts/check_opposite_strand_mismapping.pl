#!/usr/bin/env perl

#For each CLS transcript model (input GTF), find if it overlaps a coding gene on the opposite strand
#and is contained in a window around the gene coordinates with an additional buffer (eg. 400 nt)
#and the start of the transcript overlaps the start exon of one of the transcripts in the basic set
#or the end of the transcript overlaps the end exon of one of the transcripts in the basic set.

#If those conditions are met, the filter is not passed unless the CLS transcript appears to be genuine,
#i.e. it overlaps an existing gene on the same strand with > 50% exonic overlap
#or >10% exonic overlap and one splice site matches an annotated splice site.

use strict;
use warnings;
use Getopt::Long;
use Bio::Otter::Server::Config;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use List::Util qw(min max);

$|=1;

my $dataset_name;
my $host;
my $port;
my $user;
my $dbname;
my $input_gtf;
my $only_chr;
my $outfile;

&GetOptions(
            'dataset=s'  => \$dataset_name,
            'host=s'     => \$host,
            'port=i'     => \$port,
            'user=s'     => \$user,
            'dbname=s'   => \$dbname,
            'gtf=s'      => \$input_gtf,
            'chr=s'      => \$only_chr,
            'out=s'      => \$outfile,
           );

if ($only_chr){
  $only_chr =~ s/^chr//;
}

my ($sa, $ta);
if ($dataset_name){
  my $dataset = Bio::Otter::Server::Config->SpeciesDat->dataset($dataset_name);
  my $otter_dba = $dataset->otter_dba or die "can't get db adaptor\n";
  $otter_dba->dbc->reconnect_when_lost(1);
  $sa = $otter_dba->get_SliceAdaptor();
  $ta = $otter_dba->get_TranscriptAdaptor();
}
elsif ($host and $port and $user and $dbname){
  my $core_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
      -host => $host,
      -port => $port,    
      -user => $user,
      -dbname => $dbname,
  );
  $sa = $core_db->get_SliceAdaptor();
  $ta = $core_db->get_TranscriptAdaptor();
}
else{
  die "Can't connect to an annotation database: $!";
}

print "Reading input GTF...\n";
my @tids;
my %tr_coords;
my %exon_coords;
if ($input_gtf =~ /\.gz$/){
  open (IN, "zcat $input_gtf | ") or die "Can't open input file $input_gtf: $!";
}
else{
  open (IN, "cat $input_gtf | ") or die "Can't open input file $input_gtf: $!";
}
while (<IN>){
  next if /^#/;
  chomp;
  my @cols = split(/\t/);
  if ($cols[2] eq "exon"){
    my $chr = $cols[0];
    $chr =~ s/^chr//;
    if ($only_chr){
      next unless $chr eq $only_chr;
    }
    my @atts = split(/; /, $cols[8]);
    my ($start, $end, $strand) = @cols[3,4,6];
    my $transcript_id;
    foreach my $att (@atts){
      my ($field, $value) = split(/\s/, $att);
      $value =~ s/\"//g;
      $value =~ s/;//;
      if ($field eq "transcript_id"){
        $transcript_id = $value;
        last;
      }
    }
    push(@{$exon_coords{$transcript_id}}, $start."-".$end);
    if ($tr_coords{$transcript_id}){
      if ($start < $tr_coords{$transcript_id}{'start'}){
        $tr_coords{$transcript_id}{'start'} = $start;
      }
      if ($end > $tr_coords{$transcript_id}{'end'}){
        $tr_coords{$transcript_id}{'end'} = $end;
      }     
    }
    else{
      $tr_coords{$transcript_id}{'chr'} = $chr;
      $tr_coords{$transcript_id}{'start'} = $start;
      $tr_coords{$transcript_id}{'end'} = $end;
      $tr_coords{$transcript_id}{'strand'} = $strand;
      push(@tids, $transcript_id);
    }
  }
}
close (IN);

#Calculate intron coordinates for each transcript
my %intron_coords;
foreach my $transcript_id (@tids){
  my @exon_starts;
  my @exon_ends;
  foreach my $exon_coord (@{$exon_coords{$transcript_id}}){
    my ($start, $end) = split(/-/, $exon_coord);
    push(@exon_starts, $start);
    push(@exon_ends, $end);
  }
  @exon_starts = sort {$a<=>$b} @exon_starts;
  @exon_ends = sort {$a<=>$b} @exon_ends;

  @{$intron_coords{$transcript_id}{'starts'}} = map {$_ + 1} @exon_ends[0..$#exon_ends-1];
  @{$intron_coords{$transcript_id}{'ends'}} = map {$_ - 1} @exon_starts[1..$#exon_starts];
}

print "Finding overlap with coding genes...\n";
open (OUT, ">$outfile") or die "Can't open $outfile\n";

my %terminal_exons;
foreach my $transcript_id (@tids){
  my $transcript_length = 0;
  foreach my $exon_coord (@{$exon_coords{$transcript_id}}){
    my ($start, $end) = split(/-/, $exon_coord);
    $transcript_length += $end - $start + 1;
  }
  my $slice = $sa->fetch_by_region("toplevel", $tr_coords{$transcript_id}{'chr'}, $tr_coords{$transcript_id}{'start'}, $tr_coords{$transcript_id}{'end'});
  my $overlapping_genes = $slice->get_all_Genes();
  next unless scalar(grep {$_->biotype eq "protein_coding"} @{$overlapping_genes});
  
  #Look for exonic overlap with a coding gene on the same strand
  my $sense_overlap = 0;
  G:foreach my $ov_gene (@{$overlapping_genes}){
    my $gene_strand = $ov_gene->seq_region_strand == 1 ? "+" : "-";
    if ($gene_strand eq $tr_coords{$transcript_id}{'strand'}){
      #my $overlap = min($gene->seq_region_end, $tr_coords{$transcript_id}{'end'}) - max($gene->seq_region_start, $tr_coords{$transcript_id}{'start'}) + 1;
      #my $tr_span = $tr_coords{$transcript_id}{'end'} - $tr_coords{$transcript_id}{'start'} + 1;
      #next G unless $overlap/$tr_span > 0.5; #Require that CLS transcript and coding gene overlap at least for half of transcript span
                                             #This is to avoid false negatives where the end of one barely overlaps the start of the other


      #my $max_exonic_overlap = 0;
      foreach my $ov_transcript (@{$ov_gene->get_all_Transcripts}){
        my $exonic_overlap = 0;
        foreach my $exon_coord (@{$exon_coords{$transcript_id}}){
          my ($start, $end) = split(/-/, $exon_coord);
          foreach my $ov_exon (@{$ov_transcript->get_all_Exons}){
            if ($ov_exon->seq_region_start <= $end and $ov_exon->seq_region_end >= $start){
              $exonic_overlap += min($end, $ov_exon->seq_region_end) - max($start, $ov_exon->seq_region_start);
            }
          }
        }
        #if ($exonic_overlap > $max_exonic_overlap){
        #  $max_exonic_overlap = $exonic_overlap;
        #}

        my $seen_splice_sites = 0;
        my $seen_introns = 0;
        foreach my $ov_intron (@{$ov_transcript->get_all_Introns}){
          for (my $i=0; $i<scalar(@{$intron_coords{$transcript_id}{'starts'}}); $i++){
            my $seen_start = 0;
            my $seen_end = 0;
            if ($intron_coords{$transcript_id}{'starts'}[$i] == $ov_intron->seq_region_start){
              $seen_splice_sites++;
              $seen_start = 1;
            }
            if ($intron_coords{$transcript_id}{'ends'}[$i] == $ov_intron->seq_region_end){
              $seen_splice_sites++;
              $seen_end = 1;
            }
            if ($seen_start and $seen_end){
              $seen_introns++;
            }
          }
        }
        if ($exonic_overlap / $transcript_length > 0.5  #If same sense transcript exonic overlap covers at least 50% of CLS transcript,
            or ($exonic_overlap / $transcript_length > 0.3 and $seen_splice_sites >= 1)
            or ($exonic_overlap / $transcript_length > 0.1 and $seen_introns >= 1 )){ #or this has an annotated intron or two annotated splice sites, assume that this alignment is genuine
          $sense_overlap = 1;
          last G;
        }
      }

      #foreach my $exon (@{$gene->get_all_Exons}){
        #foreach my $exon_coord (@{$exon_coords{$transcript_id}}){
          #my ($start, $end) = split(/-/, $exon_coord);
          #if ($exon->seq_region_start < $end and $exon->seq_region_end > $start){
            #$sense_overlap = 1;
            #last G;
          #}
        #}
      #}
    }
  }
  if ($sense_overlap){
    next; #Skip transcript if exonic overlap with a coding gene on the same strand
  }

  #If not, look for exonic overlap with terminal exons of a coding gene on the opposite strand
  my @overlapping_coding_genes = grep {$_->biotype eq "protein_coding"} @{$overlapping_genes};
  foreach my $gene (@overlapping_coding_genes){
    my $found = 0;
    my $gene_strand = $gene->seq_region_strand == 1 ? "+" : "-";
    my $buffer = 400;
    my $window_start = $gene->seq_region_start - $buffer;
    my $window_end = $gene->seq_region_end + $buffer;
    if ($gene_strand ne $tr_coords{$transcript_id}{'strand'} and $tr_coords{$transcript_id}{'start'} >= $window_start and $tr_coords{$transcript_id}{'end'} <= $window_end){
      unless ($terminal_exons{$gene->stable_id}){
        $terminal_exons{$gene->stable_id} = get_terminal_exons($gene);
      }
      foreach my $exon_pair (@{$terminal_exons{$gene->stable_id}}){
        my ($overlaps_start_exon, $overlaps_end_exon);
        if ($exon_pair->[0]->seq_region_start <= $tr_coords{$transcript_id}{'start'} and $exon_pair->[0]->seq_region_end >= $tr_coords{$transcript_id}{'start'}){
          $overlaps_start_exon = 1;
        }
        if ($exon_pair->[1]->seq_region_start <= $tr_coords{$transcript_id}{'end'} and $exon_pair->[1]->seq_region_end >= $tr_coords{$transcript_id}{'end'}){
          $overlaps_end_exon = 1;
        }      
        if ($overlaps_start_exon or $overlaps_end_exon){
          $found = 1;
        }
      }
    }
    if ($found){
      print OUT join("\t", $transcript_id,
                          $tr_coords{$transcript_id}{'chr'}.":".$tr_coords{$transcript_id}{'start'}."-".$tr_coords{$transcript_id}{'end'}.":".$tr_coords{$transcript_id}{'strand'},
                          $gene->stable_id,
                          $gene->seq_region_name.":".$gene->seq_region_start."-".$gene->seq_region_end.":".($gene->seq_region_strand == 1 ? "+" : "-"),
                    )."\n";
      last;
    }
  }
}
close (OUT);


sub get_terminal_exons {
  my $g = shift;
  my @a;
  foreach my $t (@{$g->get_all_Transcripts}){
    if ($t->gencode_basic){
      push(@a, [sort {$a->seq_region_start <=> $b->seq_region_start} ($t->start_Exon, $t->end_Exon)]);
    }
  }
  return \@a;
}


sub get_terminal_exons_OLD {
  #find longest terminal exons
  my $g = shift;
  my @start_exons = sort {$b->length <=> $a->length} grep {$_->seq_region_start == $g->seq_region_start} @{$g->get_all_Exons};
  my @end_exons = sort {$b->length <=> $a->length} grep {$_->seq_region_end == $g->seq_region_end} @{$g->get_all_Exons};  
  return [$start_exons[0], $end_exons[0]];
}

  
