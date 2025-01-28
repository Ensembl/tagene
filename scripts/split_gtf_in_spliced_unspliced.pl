#!/usr/bin/env perl

#Split an input gtf into a gtf with spliced transcripts and a gtf with single-exon transcripts

use strict;
use warnings;
use Getopt::Long;

$|=1;


my $input_gtf;
my $spliced_gtf;
my $unspliced_gtf;

&GetOptions(
            'infile=s'      => \$input_gtf,
            'spliced=s'     => \$spliced_gtf,
            'unspliced=s'   => \$unspliced_gtf,
           );


#Read input GTF file and store transcript ids of spliced and unspliced transcripts separately
print "Reading input GTF...\n";
my %exon_count_by_tid;
if ($input_gtf =~ /\.gz$/){
  open (GTF, "zcat $input_gtf |") or die "Can't open $input_gtf: $!";
}
else{
  open (GTF, $input_gtf) or die "Can't open $input_gtf: $!";
}
while (<GTF>){
  next if /^#/;
  chomp;
  #Read exon lines (do not assume that there are transcript lines)
  my @cols = split("\t");
  if ($cols[2] eq "exon"){
    my @atts = split(/;/, $cols[8]);
    foreach my $att (@atts){
      if ($att =~ /transcript_id "(.+)"/){
        my $tid = $1;
        $exon_count_by_tid{$tid}++;
      }
    }
  }
}
close (GTF);

#Count transcripts
my $spliced_tid_count;
my $unspliced_tid_count;
foreach my $tid (keys %exon_count_by_tid){
  if ($exon_count_by_tid{$tid} == 1){
    $unspliced_tid_count++;
  }
  else{
    $spliced_tid_count++;
  }
}
print "SPLICED=$spliced_tid_count; UNSPLICED=$unspliced_tid_count\n";


#Read input GTF file and print lines into spliced or unspliced GTF files
print "Writing output GTF files...\n";
if ($spliced_gtf =~ /.gz$/){
  open (OUTSP, "| gzip -c > $spliced_gtf") or die "Can't open $spliced_gtf\n";
}
else{
  open (OUTSP, ">$spliced_gtf") or die "Can't open $spliced_gtf\n";
}
if ($unspliced_gtf =~ /.gz$/){
  open (OUTUN, "| gzip -c > $unspliced_gtf") or die "Can't open $unspliced_gtf\n";
}
else{
  open (OUTUN, ">$unspliced_gtf") or die "Can't open $unspliced_gtf\n";
} 
  
if ($input_gtf =~ /\.gz$/){
  open (GTF, "zcat $input_gtf |") or die "Can't open $input_gtf: $!";
}
else{
  open (GTF, $input_gtf) or die "Can't open $input_gtf: $!";
}
while (<GTF>){
  chomp;
  if (/^#/){
    print OUTSP $_."\n";
    print OUTUN $_."\n";
  }
  else{
    my @cols = split("\t");
    foreach my $att (split(/;/, $cols[8])){
      if ($att =~ /transcript_id "(.+)"/){
        my $tid = $1;
        if (my $exon_count = $exon_count_by_tid{$tid}){
          if ($exon_count > 1){
            print OUTSP $_."\n";
          }
          else{
            print OUTUN $_."\n";
          }
        }
        else{
          die "Can't find transcript id $tid\n";
        }
        last;
      }
    }
  }
}
close (GTF);
close (OUTSP);
close (OUTUN);
print "DONE";
