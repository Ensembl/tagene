#!/usr/bin/env perl

#This script makes the basic structure of the trackhub
#It doesn't create the actual data files

use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $name;
my $annot_file;
my $input_file;
my $host = "ensembldb.ensembl.org";
my $port = "3306";
my $user = "anonymous";
my $dbname;

&GetOptions(
             'name=s'    => \$name,
             'annot=s'   => \$annot_file,
             'input=s'   => \$input_file,
             'host=s'    => \$host,
             'port=i'    => \$port,
             'user=s'    => \$user,
             'dbname=s'  => \$dbname,
          );


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host => $host,
    -port => $port,    
    -user => $user,
    -dbname => $dbname,
);
my $meta_container = $db->get_MetaContainerAdaptor();
my $assembly_version = $meta_container->single_value_by_key('assembly.ucsc_alias');

mkdir($name);
chdir($name);

open (H, ">hub.txt") or die "Can't open hub.txt: $!";
print H "hub TAGENE\n".
        "shortLabel TAGENE\n".
        "longLabel TAGENE annotation\n".
        "genomesFile genomes.txt\n".
        "email tagene\@ebi.ac.uk\n";
close (H);

open (G, ">genomes.txt") or die "Can't open genomes.txt: $!";
print G "genome $assembly_version\ntrackDb $assembly_version/trackDb.txt\n";
close (G);

mkdir($assembly_version);

my ($annot_filename) = $annot_file =~ /\/(\w+.bb)$/;
my ($input_filename) = $input_file =~ /\/(\w+.bb)$/;

open (T, ">$assembly_version/trackDb.txt") or die "Can't open trackDb.txt: $!";
print T "track TAGENE_annotation\n".
        "type bigBed 12 + 10\n".
        "bigDataUrl ../data/$annot_filename\n".
        "shortLabel havana_annotation\n".
        "longLabel Havana + TAGENE annotation\n".
        "noScoreFilter on\n".
        "itemRgb on\n".
        "visibility full\n\n";
        
print T "track input_annotation\n".
        "type bigBed 12 + 6\n".
        "bigDataUrl ../data/$input_filename\n".
        "shortLabel input_annotation\n".
        "longLabel Input gene annotation\n".
        "noScoreFilter on\n".
        "itemRgb on\n".
        "visibility full\n\n";
close (T);

mkdir("data");
system("cp -p $annot_file $input_file data/");



