#!/usr/bin/perl -w

# Given a GTF file of long read models and 
# a list of Ensembl gene stable ids,
# this script prints out a GTF file with only transcripts
# that overlap the Ensembl genes in the list.
# The genes are looked up in the Havana database.

use strict;
use warnings;
use Getopt::Long;
use Gencode::Default;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


my $gene_list;
my $gtf_infile;
my $outfile;
my $dbhost;
my $dbport;
my $dbuser;
my $dbpass;
my $dbname;

&GetOptions(
            'list=s'     => \$gene_list,
            'infile=s'   => \$gtf_infile,
            'host=s'     => \$dbhost,
            'port=s'     => \$dbport,
            'user=s'     => \$dbuser,
            'pass=s'     => \$dbpass,
            'dbname=s'   => \$dbname,
            'outfile=s'  => \$outfile,
            );            

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                           -host    => $dbhost,
                                           -user    => $dbuser,
                                           -port    => $dbport,
                                           -dbname  => $dbname,
                                           ) or die "Can't connect to database";
my $ga = $db->get_GeneAdaptor();

#Store gene coordinates
my %gene_coords;
open (LIST, $gene_list) or die "Can't open $gene_list: $!";
while (<LIST>){
  chomp;
  my $gid = $_;
  my $gene = $ga->fetch_by_stable_id($gid);
  if ($gene){
    push(@{$gene_coords{$gene->seq_region_name}{Gencode::Default->strand_value($gene->strand)}}, [$gene->start, $gene->end]);
  }
  else{
    warn "Can't find gene $gid!\n";
  }     
}
close (LIST);


#Find long read models that overlap genes on the list
#Store transcript ids
#As the input GTF may contain exons only, do not print the output GTF now (some exons could be left out)
my %tids;
open (GTF, "gunzip -dc $gtf_infile |") or die "Can't open gtf file:$!";
while (<GTF>){
  next unless /\w+/;
  my ($fs, $attribs) = Gencode::Default->parse_gtf_line($_);
  if (my $tid = $attribs->{'transcript_id'}->[0]){
    next if $tids{$tid};
  }
  if ($fs->{'type'} eq "exon"){
    if ($fs->{'chr'} =~ /^chr/){
      $fs->{'chr'} =~ s/^chr//;
    }
    if (exists($gene_coords{$fs->{'chr'}}{$fs->{'strand'}})){
      foreach my $pair (@{$gene_coords{$fs->{'chr'}}{$fs->{'strand'}}}){
        my $gene_start = $pair->[0];
        my $gene_end = $pair->[1];
        if ($fs->{'start'} <= $gene_end and $fs->{'end'} >= $gene_start){
          if (my $tid = $attribs->{'transcript_id'}->[0]){
            $tids{$tid} = 1;
          }
        }
      }
    }
  }
}
close (GTF);
 

#Read the input file again and print the output GTF
open (OUT, ">$outfile") or die "Can't open $outfile: $!";
my %printed_lines;    
open (GTF, "gunzip -dc $gtf_infile |") or die "Can't open gtf file:$!";
while (<GTF>){
  next unless /\w+/;
  my ($fs, $attribs) = Gencode::Default->parse_gtf_line($_); 
  if ($fs->{'type'} eq "exon"){   
    if (my $tid = $attribs->{'transcript_id'}->[0]){
      if ($tids{$tid}){ 
        unless ($printed_lines{$_}){
          print OUT $_;
          $printed_lines{$_} = 1;           
        }
      }
    }
  }
}
close (GTF);
close (OUT);


