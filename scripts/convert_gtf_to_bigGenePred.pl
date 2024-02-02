#!/usr/bin/env perl

#This script converts annotation from a GTF file into a bigGenePred file by using UCSC tools

use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


my $gtf_file;
my $big_genepred_file;
my $host;
my $port;
my $user;
my $dbname;

&GetOptions(
             'gtf=s'     => \$gtf_file,
             'host=s'    => \$host,
             'port=i'    => \$port,
             'user=s'    => \$user,
             'dbname=s'  => \$dbname,
             'out=s'     => \$big_genepred_file,
          );

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
     -host   => $host,
     -port   => $port,
     -user   => $user,
     -pass   => undef,
     -dbname => $dbname,
     -driver => 'mysql',
);

my $sa = $db->get_SliceAdaptor();

my $meta_container = $db->get_MetaContainerAdaptor();
my $assembly_version = $meta_container->single_value_by_key('assembly.ucsc_alias');

system("gtfToGenePred -genePredExt $gtf_file file.gp");
system("genePredToBigGenePred file.gp stdout > file.bgpInput");
unless (-e "$assembly_version.chrom.sizes"){
  system("wget https://hgdownload.soe.ucsc.edu/goldenPath/$assembly_version/bigZips/$assembly_version.chrom.sizes");
}
check_chromosome_names("file.bgpInput", "$assembly_version.chrom.sizes");
system("LC_COLLATE=C sort -k1,1 -k2,2n file.bgpInput.2 > file.bgpInput");
print_bgp_definition_file("bigGenePred.as");
system("bedToBigBed -type=bed12+6 -tab -as=bigGenePred.as file.bgpInput $assembly_version.chrom.sizes $big_genepred_file");


###


sub check_chromosome_names {
  my ($bgp_file, $chr_file) = @_;

  my %chrs;
  open (C, $chr_file) or die "Can't open $chr_file: $!";
  while (<C>){
    my @cols = split("\t");
    my $chr = $cols[0];
    $chrs{$chr} = 1;
  }
  close (C);

  open (OUT, ">$bgp_file.2") or die "Can't open $bgp_file.2: $!";
  open (IN, $bgp_file) or die "Can't open $bgp_file: $!";
  while (<IN>){
    my @cols = split("\t");
    my $chr = $cols[0];
    if ($chrs{$chr}){
      print OUT $_;
    }
    else{
      my $slice = $sa->fetch_by_region('toplevel', $chr) or die "Can't get slice for $chr";
      my $ucsc_synonyms = $slice->get_all_synonyms('UCSC');
      if (scalar @$ucsc_synonyms){
        my $syn = $ucsc_synonyms->[0]->name;
        shift @cols;
        print OUT $syn."\t".join("\t", @cols);
      }
      else{
        print "No synonym for $chr\n";
      }
    }
  }
  close (OUT);

}


  
sub print_bgp_definition_file {
  my $file = shift;
  open ("T", ">$file") or die "Can't open $file: $!";
  my $text = <<TABLE;   #Extended from https://genome.ucsc.edu/goldenpath/help/bigGenePred.html
table bigGenePred
"genes"
   (
   string chrom;        "Reference sequence chromosome or scaffold"
   uint chromStart;     "Start position in chromosome"
   uint chromEnd;       "End position in chromosome"
   string name;         "Feature name"
   uint score;          "Score (0-1000)"
   char[1] strand;      "+ or - for strand"
   uint thickStart;     "Start of where display should be thick (start codon)"
   uint thickEnd;       "End of where display should be thick (stop codon)"
   uint reserved;       "RGB value"
   int blockCount;      "Number of blocks"
   int[blockCount] blockSizes; "Comma separated list of block sizes"
   int[blockCount] chromStarts; "Start positions relative to chromStart"
   string name2;        "Alternative/human readable name"
   string cdsStartStat; "Status of CDS start annotation (none, unknown, incomplete, or complete)"
   string cdsEndStat;   "Status of CDS end annotation (none, unknown, incomplete, or complete)"
   int[blockCount] exonFrames; "Exon frame {0,1,2}, or -1 if no frame for exon"
   string type;         "Transcript biotype"
   string geneName;     "Primary identifier for gene"
   string geneName2;    "Alternative/human readable gene name"
   string geneType;     "Gene biotype"
   )
TABLE
  print T $text."\n";
  close (T);
};



