#!/usr/bin/env perl

#Check if the start and end exon sizes are smaller than a given threshold (default is 15nt).
#Both exons must have the minimum size for the transcript to pass this filter.

use strict;
use warnings;
use Getopt::Long;
use List::Util qw(uniq);

$|=1;

my $gtf_file;
my $id_list;
my $minimum_size = 15;
my $only_chr;
my $outfile;

&GetOptions(
            'gtf=s'      => \$gtf_file,
            'idlist=s'   => \$id_list,
            'chr=s'      => \$only_chr,
            'minsize=i'  => \$minimum_size,
            'out=s'      => \$outfile,
           );

if ($only_chr){
  $only_chr =~ s/^chr//;
}


#Read list of selected RNA-seq transcript ids, if any
my @selected_ids;
my %selected_id_list;
if ($id_list){
  open (IN, $id_list) or die "Can't open $id_list: $!";
  while (<IN>){
    chomp;
    push(@selected_ids, $_);
    $selected_id_list{$_} = 1;
  }
  close (IN);
}


#Read input GTF file
my %exon_coords_by_tid;
if ($gtf_file =~ /\.gz$/){
  open (GTF, "zcat $gtf_file |") or die "Can't open $gtf_file : $!";
}
else{
  open (GTF, "cat $gtf_file |") or die "Can't open $gtf_file : $!";
}
while (<GTF>){
  next if /^#/;
  chomp;
  #Read exon lines (do not assume that there are transcript lines)
  my @cols = split("\t");
  if ($cols[2] eq "exon"){
    my $chr = $cols[0];
    $chr =~ s/^chr//;
    if ($only_chr){
      next unless $chr eq $only_chr;
    }
    if ($cols[8] =~ /transcript_id \"(\S+)\";/){
      my $tid = $1;
      #Filter by ID list if any
      if ($id_list){
        next unless $selected_id_list{$tid};
      }
      else{
        push(@selected_ids, $tid);
      }
      #Store exon coordinates
      my ($start, $end, $strand) = @cols[3,4,6];
      $exon_coords_by_tid{$tid}{'chr'} = $chr;
      $exon_coords_by_tid{$tid}{'strand'} = $strand;
      push(@{$exon_coords_by_tid{$tid}{'starts'}}, $start);
      push(@{$exon_coords_by_tid{$tid}{'ends'}}, $end);  
    }
    else{
      die "No transcript_id field!";
    }
  }
}
close (GTF);


open (OUT, ">$outfile") or die "Can't open $outfile: $!";
print OUT "#".join("\t", "transcript_ID",
                        "passes_filter",
                        "start_exon_size",
                        "end_exon_size",
                  )."\n";


#Loop through transcripts
foreach my $tid (uniq @selected_ids){
  my $passed = 0;
  my @number_of_exons = @{$exon_coords_by_tid{$tid}{'starts'}};
  my @starts = sort {$a<=>$b} @{$exon_coords_by_tid{$tid}{'starts'}};
  my @ends = sort {$a<=>$b} @{$exon_coords_by_tid{$tid}{'ends'}};
  my ($start_exon_size, $end_exon_size);
  if ($exon_coords_by_tid{$tid}{'strand'} eq "+"){
    $start_exon_size = $ends[0] - $starts[0] + 1;
    $end_exon_size = $ends[$#number_of_exons] - $starts[$#number_of_exons] + 1;
  }
  else{
    $start_exon_size = $ends[$#number_of_exons] - $starts[$#number_of_exons] + 1;    
    $end_exon_size = $ends[0] - $starts[0] + 1;
  }  
  if ($start_exon_size >= $minimum_size and $end_exon_size >= $minimum_size){
    $passed = 1;
  }
  print OUT join("\t", $tid, ($passed ? "yes" : "no"), $start_exon_size, $end_exon_size)."\n";
}
close (OUT);

