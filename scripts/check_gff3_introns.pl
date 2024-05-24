#!/usr/bin/perl -w


#This script produces a breakdown of introns by splice site sequence in GFF3 annotation


use strict;
use Getopt::Long;
use Gencode::Default;


$|=1;

my $infile;
my $outfile;
my $species = "human";
my $assembly;

&GetOptions(
	    'in=s'      => \$infile,
	    'out=s'     => \$outfile,
            'species=s' => \$species,
            'asv=s'     => \$assembly,
           );



die "Please provide assembly version\n" unless $assembly;

#Database to provide sequence
my $db;
if ($species eq "human"){
  #$assembly = "GRCh38" unless $assembly;
  if ($assembly eq "GRCh37"){
    $db = Gencode::Default->dbconnect("ensdb-archive", 5304, "homo_sapiens_core_75_37", "ensro", undef);
  }
  elsif ($assembly eq "GRCh38"){
    $db = Gencode::Default->dbconnect("ensdb-archive", 5304, "homo_sapiens_core_80_38", "ensro", undef);
  }
}
elsif ($species eq "mouse"){
  $assembly = "GRCm38" unless $assembly;
  if ($assembly eq "GRCm38"){
    $db = Gencode::Default->dbconnect("ensdb-archive", 5304, "mus_musculus_core_80_38", "ensro", undef);
  }
}

print "\nUsing $species $assembly genome assembly version.\n\n";

my $sa = $db->get_SliceAdaptor;

open (OUT, ">$outfile") or die "Can't open $outfile\n";
my %trs;
open (GFF3, $infile) or die "Can't open $infile\n";
while (<GFF3>){
  next if /^#/;
  chomp;
  my @cols = split(/\s+/);
  if ($cols[2] eq "cDNA_match" or $cols[2] eq "exon"){
    my $chr = $cols[0];
    $chr =~ s/chr//;
    my $start = $cols[3];
    my $end = $cols[4];
    my $strand = $cols[6];
    $strand =~ s/\+/1/;
    $strand =~ s/\-/-1/;
    my $id = $cols[8];
    $trs{$id}{chr} = $chr;
    push(@{$trs{$id}{starts}}, $start);
    push(@{$trs{$id}{ends}}, $end);
    $trs{$id}{strand} = $strand;
  }
}
close (GFF3);

my %intron_sites;
foreach my $tid (keys %trs){
  my @starts = sort {$a <=> $b} @{$trs{$tid}{starts}};
  my @ends = sort {$a <=> $b} @{$trs{$tid}{ends}};
  for (my $i=0; $i<$#starts; $i++){
    #Get intron coordinates
    my $intron_start = $ends[$i]+1;
    my $intron_end = $starts[$i+1]-1;
    #Get intron sequence
    my $istart_slice = $sa->fetch_by_region("toplevel", $trs{$tid}{chr}, $intron_start, $intron_start+1, $trs{$tid}{strand});
    my $iend_slice = $sa->fetch_by_region("toplevel", $trs{$tid}{chr}, $intron_end-1, $intron_end, $trs{$tid}{strand});
    next unless $istart_slice and $iend_slice;
    my $isite_seq;
    if ($trs{$tid}{strand}==1){
      $isite_seq = $istart_slice->seq."..".$iend_slice->seq;
    }
    elsif ($trs{$tid}{strand}==-1){
      $isite_seq = $iend_slice->seq."..".$istart_slice->seq;
    }
    print OUT $tid."\t".$isite_seq."\n";
    $intron_sites{$isite_seq}++;
  }
}
close (OUT);

foreach my $site (sort {$intron_sites{$b} <=> $intron_sites{$a}} keys %intron_sites){
  print $site."\t".$intron_sites{$site}."\n";
}


