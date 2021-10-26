#!/usr/bin/perl -w

# This script gets the full-length supporting reads for each transcript from the CLS3 gff files.
# The 'longest_FL_supporters' attribute is expected, eg.
#chr1    tmerge  exon    1311599 1311924 0       -       .       gene_id "pacBio:Cshl:CapTrap:Corr0_HpreCap_0+_Brain01Rep1.NAM_TM_000000000148"; transcript_id "pacBio:Cshl:CapTrap:Corr0_HpreCap_0+_Brain01Rep1.NAM_TM_000000000148"; contains "m54333_190125_203311/35062307/ccs,m54333_190125_203311/73597686/ccs,m54333_190125_203311/33096297/ccs,m54333_190125_203311/37814484/ccs"; contains_count "4"; 3p_dists_to_3p "0,1,1,1"; 5p_dists_to_5p "0,0,0,2"; flrpm "2.011296"; longest "m54333_190125_203311/35062307/ccs"; longest_FL_supporters "m54333_190125_203311/35062307/ccs,m54333_190125_203311/33096297/ccs,m54333_190125_203311/37814484/ccs,m54333_190125_203311/73597686/ccs"; longest_FL_supporters_count "4"; mature_RNA_length "2115"; meta_3p_dists_to_5p "1,0.999527186761229,0.999527186761229,0.999527186761229"; meta_5p_dists_to_5p "0,0,0,0.000945626477541371"; rpm "2.011296"; spliced "1";

use strict;
use warnings;
use Getopt::Long;


my $infile;
my $dataset_name;
my $outfile;

&GetOptions(
            'infile=s'   => \$infile,
            'libname=s'  => \$dataset_name,
            'outfile=s'  => \$outfile,           
            );            


open (OUT, ">$outfile") or die "Can't open $outfile: $!";
my %seen_tids;
open (GTF, "gunzip -dc $infile |") or die "Can't open gtf file:$!";
while (<GTF>){
  chomp;
  if (/transcript_id \"(\S+)\"; contains \"(\S+)\";/){
    my $tid = $1;
    my $reads = $2;
    unless ($seen_tids{$tid}){
      $reads =~ s/\//_/g;
      print OUT $tid."\t".$dataset_name." : ".join(", ", split(/,/, $reads))." ;\n";
    }
    if (!$reads){
      warn "No FL reads for $tid\n";
    }
    $seen_tids{$tid} = 1;
  }
}
close (GTF);
close (OUT);


