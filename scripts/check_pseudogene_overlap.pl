#!/usr/bin/env perl

#For each CLS transcript model (input GTF), find if it overlaps a pseudogene on either strand and
#is encompassed by the pseudogene boundaries extended with the maximum 5' and 3'UTRs found in the gene annotation
#(to account for pseudogenized UTRs, which are not included in the pseudogene annotation).


use strict;
use warnings;
use Getopt::Long;
use Bio::Otter::Server::Config;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use List::Util qw( max sum );

$|=1;

my $dataset_name;
my $host;
my $port;
my $user;
my $dbname;
my $input_gtf;
my $only_chr;
my $utr5_length;
my $utr3_length;
my $outfile;

&GetOptions(
            'dataset=s'  => \$dataset_name,
            'host=s'     => \$host,
            'port=i'     => \$port,
            'user=s'     => \$user,
            'dbname=s'   => \$dbname,
            'gtf=s'      => \$input_gtf,
            'chr=s'      => \$only_chr,   
            'utr5=s'     => \$utr5_length,
            'utr3=s'     => \$utr3_length,
            'out=s'      => \$outfile,
           );


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
open (IN, "zcat $input_gtf | ") or die "Can't open input file $input_gtf: $!";
while (<IN>){
  next if /^#/;
  chomp;
  my @cols = split(/\t/);
  if ($cols[2] eq "exon"){
    if ($only_chr){
      unless ($cols[0] eq $only_chr or "chr".$cols[0] eq $only_chr or $cols[0] eq "chr".$only_chr){
        next;
      }
    }
    my @atts = split(/; /, $cols[8]);
    my ($chr, $start, $end, $strand) = @cols[0,3,4,6];
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


#Find max 5' and 3' UTR lengths
#$utr5_length = 0;
#$utr3_length = 0;

unless ($utr5_length and $utr3_length){
  print "Finding maximum UTR lengths in current annotation...\n";
  my @five_prime_utrs;
  my @three_prime_utrs;
  foreach my $transcript (@{$ta->fetch_all_by_biotype("protein_coding")}){
    if ($transcript->five_prime_utr){
      push(@five_prime_utrs, $transcript->five_prime_utr->length);
    }
    if ($transcript->three_prime_utr){
      push(@three_prime_utrs, $transcript->three_prime_utr->length);
    }
  }
  my $max_5utr = max(@five_prime_utrs);
  my $max_3utr = max(@three_prime_utrs);
  my $median_5utr = median(@five_prime_utrs);
  my $median_3utr = median(@three_prime_utrs);
  print "MAX 5'UTR=$max_5utr, MAX 3'UTR=$max_3utr\n";
  print "MEDIAN 5'UTR=$median_5utr, MEDIAN 3'UTR=$median_3utr\n";
  $utr5_length = $median_5utr;
  $utr3_length = $median_3utr;
}


print "Finding overlap with pseudogenes...\n";
open (OUT, ">$outfile") or die "Can't open $outfile\n";

foreach my $transcript_id (@tids){
  my $slice = $sa->fetch_by_region("toplevel", $tr_coords{$transcript_id}{'chr'}, $tr_coords{$transcript_id}{'start'}, $tr_coords{$transcript_id}{'end'});
  foreach my $gene (@{$slice->get_all_Genes}){
    my $overlap = 0;
    if ($gene->biotype =~ /^(processed_pseudogene|transcribed_processed_pseudogene|translated_processed_pseudogene)$/){
      if ($gene->seq_region_strand == 1){
        if ($tr_coords{$transcript_id}{'start'} >= $gene->seq_region_start-$utr5_length and $tr_coords{$transcript_id}{'end'} <= $gene->seq_region_end+$utr3_length){
          $overlap = 1;
        }
      }
      elsif ($gene->seq_region_strand == -1){
        if ($tr_coords{$transcript_id}{'start'} >= $gene->seq_region_start-$utr3_length and $tr_coords{$transcript_id}{'end'} <= $gene->seq_region_end+$utr5_length){
          $overlap = 1;
        }
      }
    }
    if ($overlap){
      print OUT join("\t", $transcript_id,
                          $tr_coords{$transcript_id}{'chr'}.":".$tr_coords{$transcript_id}{'start'}."-".$tr_coords{$transcript_id}{'end'}.":".$tr_coords{$transcript_id}{'strand'},
                          $gene->stable_id,
                          "chr".$gene->seq_region_name.":".$gene->seq_region_start."-".$gene->seq_region_end.":".($gene->seq_region_strand ? "+" : "-"),
                          $gene->biotype,
                    )."\n";
      last;
    }
  }
}
close (OUT);



sub median {
  my @sorted_array = sort {$a<=>$b} @_;
  if (scalar @_ % 2){
    return $sorted_array[int(scalar(@_)/2)];
  }
  else{
    return ($sorted_array[int(scalar(@_)/2)-1] + $sorted_array[int(scalar(@_)/2)])/2;
  }
}
