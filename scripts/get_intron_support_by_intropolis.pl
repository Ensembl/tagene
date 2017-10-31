#!/usr/bin/perl -w

#This script provides a break down of transcripts in a database by level of intron support according to the Intropolis dataset.
#The intron support level of a transcript will be that of the less supported intron.

use strict;
use Getopt::Long;
use Gencode::Default;


$|=1;

my $infile;
my $outfile;
my $dbhost;
my $dbport;
my $dbuser;
my $dbpass;
my $dbname;
my $biotype;

&GetOptions(
            'in=s'      => \$infile, #Intropolis file with intron coordinates and number of experiments and reads
            'host=s'    => \$dbhost,
            'port=s'    => \$dbport,
            'user=s'    => \$dbuser,
            'pass=s'    => \$dbpass,
            'dbname=s'  => \$dbname,
            'biotype=s' => \$biotype,
            'out=s'     => \$outfile,
           );

#Read Intropolis file
my %intropolis;
if ($infile =~ /\.gz$/){
  open(IN, "gunzip -dc $infile |") or die "Cant open file $infile\n";
}
else{
  open(IN, $infile) or die "Can't open file $infile\n";
}
while (<IN>){
  #chr1    10113   10166   -       GT      AG      13028,17289,17307       1,5,1   chr1    10113   10166   -
  chomp;
  my @cols = split(/\t/);
  my $chr = $cols[8];
  my $start = $cols[9];
  my $end = $cols[10];
  my $strand = $cols[11];
  my $nreads;
  foreach (split(/,/, $cols[7])){
    $nreads += $_;
  }
  $intropolis{"$chr:$start-$end:$strand"} = $nreads;
}
close (IN);

open (OUT, ">$outfile") or die "Can't open $outfile\n";

my %intron_count_by_gtype;
my %intron_count_by_support;

#Find Intropolis support for each annotated intron
my $db = Gencode::Default->dbconnect($dbhost, $dbport, $dbname, $dbuser, $dbpass);
my $sa = $db->get_SliceAdaptor;
foreach my $slice (@{$sa->fetch_all('chromosome', undef, 0, 1)}){
  #next unless $slice->seq_region_name =~ /22/;
  print $slice->seq_region_name."...";
  my $chr = $slice->seq_region_name;
  $chr = "chr$chr" unless $chr =~ /^chr/;;
  $chr =~ s/(chr\w+)-38/$1/ if $chr =~ /chr\w+-38/;
  foreach my $gene (@{$slice->get_all_Genes}){
    foreach my $transcript (@{$gene->get_all_Transcripts}){
      if ($biotype){
        next unless $transcript->biotype eq $biotype;
      }
      my @introns = @{$transcript->get_all_Introns};
      next unless scalar @introns;
      #my $tsl = scalar(@{$transcript->get_all_Attributes('TSL')}) ? $transcript->get_all_Attributes('TSL')->[0]->value : "NA";
      print OUT $gene->stable_id."\t".$gene->biotype."\t".$transcript->stable_id."\t".$transcript->biotype."\t".$transcript->source."\t";
      my @supp_read_n;
      my $min_support;
      my @info;
      foreach my $intron (@introns){
        my $strand = Gencode::Default->strand_value($intron->seq_region_strand);
        my $coord = $chr.":".$intron->start."-".$intron->end.":".$strand;
        $intron_count_by_gtype{$gene->biotype}{$coord}++;
        if (my $val = $intropolis{$coord}){
          push(@supp_read_n, $val);
          push(@info, $coord."=".$val);
          $intron_count_by_support{get_category($val)}{$coord}++;
        }
        else{
          $min_support = "NA";
          push(@info, $coord."=NA");
          $intron_count_by_support{"none"}{$coord}++;
          #last;
        }
      }
      unless ($min_support){
        my @s = sort {$a <=> $b} @supp_read_n;
        $min_support = shift @s;
      }

      print OUT scalar(@introns)."\t".get_category($min_support)."\t".$min_support."\t".join("||", @info)."\n";
    }
  }
}
close (OUT);

print "\n\n";
foreach my $gtype (keys %intron_count_by_gtype){
  print $gtype."\t".scalar(keys %{$intron_count_by_gtype{$gtype}})."\n";
}
print "\n\n";
foreach my $categ (keys %intron_count_by_support){
  print $categ."\t".scalar(keys %{$intron_count_by_support{$categ}})."\n";
}
print "\n\n";


###

sub get_category {
  my $support_value = shift;
  my $categ;
  if ($support_value eq "NA"){
    $categ = "none";
  }
  elsif ($support_value >= 1 and $support_value < 100){
    $categ = "weak";
  }
  elsif ($support_value >= 100 and $support_value < 1000){
    $categ = "medium";
  }
  elsif ($support_value >= 1000){
    $categ = "strong";
  }
  return $categ;
}



