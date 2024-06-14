#!/usr/bin/env perl

#For each CLS transcript model (input GTF), find if it has an intron that falls entirely within a single tandem repeat (trf) region


use strict;
use warnings;
use Getopt::Long;
use Bio::DB::HTS::Tabix;

$|=1;


my $repeat_file;
my $input_gtf;
my $only_chr;
my $outfile;

&GetOptions(
            'repeats=s'  => \$repeat_file,
            'gtf=s'      => \$input_gtf,
            'chr=s'      => \$only_chr,
            'out=s'      => \$outfile,
           );



#Open Tabix for the repeat file
my $repeat_tabix = Bio::DB::HTS::Tabix->new( filename => $repeat_file );


#Read GTF file with CLS models
print "Reading input GTF...\n";
my %exon_coords_by_tid;
my @tids;
open (GTF, "zcat $input_gtf |") or die "Can't open $input_gtf: $!";
while (<GTF>){
  next if /^#/;
  chomp;
  #Read exon lines (do not assume that there are transcript lines)
  my @cols = split("\t");
  if ($cols[2] eq "exon"){
    if ($only_chr){
      unless ($cols[0] eq $only_chr or "chr".$cols[0] eq $only_chr or $cols[0] eq "chr".$only_chr){
        next;
      }
    }
    if ($cols[8] =~ /transcript_id \"(\w+)\";/){
      my $tid = $1;
      #Add to ID list
      unless ($exon_coords_by_tid{$tid}){
        push(@tids, $tid);
      }
      #Store exon coordinates
      my ($chr, $start, $end, $strand) = @cols[0,3,4,6];
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

#Calculate intron coordinates for each transcript
print "Calculating intron coordinates...\n";
my %intron_coords_by_tid;
foreach my $tid (keys %exon_coords_by_tid){
  $intron_coords_by_tid{$tid}{'chr'} = $exon_coords_by_tid{$tid}{'chr'};
  $intron_coords_by_tid{$tid}{'strand'} = $exon_coords_by_tid{$tid}{'strand'};
  if (scalar @{$exon_coords_by_tid{$tid}{'starts'}} > 1){
    my @starts = sort {$a<=>$b} @{$exon_coords_by_tid{$tid}{'starts'}};
    my @ends = sort {$a<=>$b} @{$exon_coords_by_tid{$tid}{'ends'}};
    @{$intron_coords_by_tid{$tid}{'starts'}} = map {$_ + 1} @ends[0..$#ends-1];
    @{$intron_coords_by_tid{$tid}{'ends'}} = map {$_ - 1} @starts[1..$#starts];
  }
}


print "Finding overlap with repeats...\n";
open (OUT, ">$outfile") or die "Can't open $outfile\n";

my %seen_introns;
foreach my $tid (@tids){
  next unless $intron_coords_by_tid{$tid}{'starts'};
  my $trf_overlap = 0;
  my $failed_intron_coord;
  my $repeat_feat_coord = "na";
  INT:for (my $i=0; $i<scalar(@{$intron_coords_by_tid{$tid}{'starts'}}); $i++){
    my $coord = $intron_coords_by_tid{$tid}{'chr'}.":".$intron_coords_by_tid{$tid}{'starts'}[$i]."-".$intron_coords_by_tid{$tid}{'ends'}[$i];
    if (my $overlap = $seen_introns{$coord}){
      if ($overlap eq "yes"){
        $trf_overlap = 1;
        $failed_intron_coord = $coord;
      }
    }
    else{
      my $iter = $repeat_tabix->query($coord);
      my $overlap = "no";
      if ($iter){
        while (my $line = $iter->next){
          my @cols = split(/\t/, $line);
          if ($cols[1]+1 <= $intron_coords_by_tid{$tid}{'starts'}[$i] and $cols[2] >= $intron_coords_by_tid{$tid}{'ends'}[$i]){
            $overlap = "yes";
            $trf_overlap = 1;
            $failed_intron_coord = $coord;
            $repeat_feat_coord = $cols[3];
            last;
          }
        }
      }
      $seen_introns{$coord} = $overlap;
    }

    if ($trf_overlap){
      print OUT join("\t", $tid,
                          $failed_intron_coord,
                          $repeat_feat_coord,
                    )."\n";
      last INT;
    }
  }
}
close (OUT);



