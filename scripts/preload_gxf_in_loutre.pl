#!/usr/bin/perl -w

#This script splits a large input annotation file (in GTF or GFF format) into smaller chunks
#which are written in individual temporary files, 
#and then submits an individual job with the load_gxf_in_loutre.pl script for each of them.
#The script will group transcripts in clusters by simple genomic overlap by default,
#although it will only include all transcripts of the same cluster in the same file if
#the right option is given.
#It takes the same options as load_gxf_in_loutre.pl, which will be passed on to that script,
#in addition to some unique options:
#-the maximum number of transcripts for each file;
#-the maximum number of clusters per file;
#-the maximum number of files.
#Usage: perl preload_gxf_in_loutre.pl -file data.gtf -n 1000 -dataset human_test ...


use strict;
use warnings;
use Getopt::Long;
use List::Util qw(shuffle);
use File::Path qw(make_path);
use POSIX;
$| = 1;

my $file;
my $source_info;
my $write;
my $dataset_name;
my $author_name;
my $remark;
my $use_comp_pipe_biotype;
my $no_artifact_check;
my $analysis_name;
my $tsource;
my $assembly_version;
my $no_NFV;
my $do_not_add_cds;
my $no_intron_check;
my $host_biotype;
my $max_overlapped_loci;
my $filter_introns;
my $platinum;
my $only_chr;
my $max_trs_per_file;
my $max_clusters_per_file;
my $max_num_files;
my $tag;
my $dep_job;
my $tmpdir;
my $read_seq_dir;
my $cluster_by_gene_id;
my $no_shuffle;
my $job_limit = 20; #Max number of simultaneously running jobs in a job array

&GetOptions(
            'file=s'            => \$file,
            'dataset=s'         => \$dataset_name,
            'author=s'          => \$author_name,
            'source=s'          => \$source_info, #file with otter columns and read names supporting each transcript model
            'remark=s'          => \$remark, #remark to be added to every transcript model
            'readseqdir=s'      => \$read_seq_dir,
            'comp_pipe!'        => \$use_comp_pipe_biotype,
            'analysis=s'        => \$analysis_name,
            'tsource=s'         => \$tsource,
            'assembly=s'        => \$assembly_version,
            'no_check!'         => \$no_artifact_check,
            'no_NFV!'           => \$no_NFV,
            'no_CDS!'           => \$do_not_add_cds,
            'no_intron_check!'  => \$no_intron_check,
            'host_biotype=s'    => \$host_biotype,
            'max_ov_loc=i'      => \$max_overlapped_loci,
            'filter_introns!'   => \$filter_introns,
            'platinum!'         => \$platinum,
            'chr=s'             => \$only_chr,
            'write!'            => \$write,
            'max_trs=i'         => \$max_trs_per_file,
            'max_clus=i'        => \$max_clusters_per_file,
            'max_files=i'       => \$max_num_files,
            'tag=s'             => \$tag,
            'dep_job=s'         => \$dep_job,
            'tmpdir=s'          => \$tmpdir,
            'by_gene_id!'       => \$cluster_by_gene_id,
            'no_shuffle!'       => \$no_shuffle,
            'job_limit=i'       => \$job_limit,
            );


my $usage =  <<_USAGE_;

This script splits a gene annotation file (in gtf/gff3 format) in (clusters of) transcripts,
stores them in individual annotation files and submit an LSF job with the load_gxf_in_loutre.pl script for each file.

perl preload_gxf_in_loutre.pl -file ANNOTATION_FILE [OPTIONS]
 -file           annotation file in GTF or GFF3 format - only 'exon' lines will be read in
 -max_trs        maximum number of transcripts per file (transcripts will be randomly assigned to files irrespective of their clusters)
 -max_clus       maximum number of clusters per file
 -max_files      maximum number of files (used in conjunction with -max_clus)
 -tag            prefix for the individual file names and the LSF job names
 -dep_job        tag of another set of jobs that will be submitted simultaneously and that should end before this set of jobs can start
 -tmpdir         name of the directory that will contain the individual annotation files and LSF job output files
 -by_gene_id     use the 'gene_id' attribute instead of the 'transcript_id' attribute, to make sure that all transcripts under the same gene id are processed together
 -no_shuffle     do not shuffle the array of transcript clusters before creating the individual annotation files
 -job_limit      maximum number of LSF jobs to be run simultaneously

 OPTIONS TO BE PASSED ON TO THE load_gxf_in_loutre.pl SCRIPT:
 -dataset        dataset name (species) already known by the Otter system, ie. present in the species.dat file on the live server
 -source         tsv file with transcript id and transcript sources, ie. RNA-seq read names supporting the model or tissues where it was found
 -readseqdir     directory with supporting read sequences in Fasta format
 -author         author name in the corresponding loutre database
 -remark         annotation remark to be added to all transcripts, usually the ENA/GEO accession for the experiment that generated the models
 -comp_pipe      override the given biotypes and use "comp_pipe"
 -analysis       analysis logic name (default is "Otter")
 -tsource        transcript source (default is "havana")
 -assembly       assembly version of the annotation in the file, if not the default one
 -no_check       do no check for transcripts spanning a large number of genes
 -no_NFV         do not add the 'not for VEGA' attribute
 -no_CDS         do not try to add a CDS if the transcript falls in a coding gene
 -no_intron_check  allow transcripts with intron chains fully or partially identical to others in the database
 -host_biotype   restrict host genes by biotype (comma-separated list)
 -max_ov_loc     maximum number of existing loci that a novel transcript can overlap at the exon level (ignore the transcript if exceeded)
 -filter_introns assess introns and ignore transcript if at least an intron does not pass the filters
 -platinum       add a 'platinum' hidden remark to all transcripts
 -chr            restrict to annotation on this chromosome
 -write          store the gene annotation in the loutre database

_USAGE_

die $usage unless ($file and $dataset_name and $author_name);


if ($max_trs_per_file and $max_clusters_per_file or (!$max_trs_per_file and !$max_clusters_per_file)){
  die "Please provide one (and only one) of the -max_trs and -max_clus options\n";
}


#Double check in case of writing a production database
if ($write and ($dataset_name eq "human" or $dataset_name eq "mouse")){
  print STDERR "You are about to modify the $dataset_name production database.\n
                Do you want to proceed? yes | [no]\n";
  my $answer = <STDIN>;
  chomp $answer;
  unless ($answer eq "yes"){
    print STDERR "Stopping here...\n";
    exit(0);
  }
}


#Read input annotation file
my %transcript_data;
my %seen;
my $filetype;
if ($file =~ /\.gff3(.gz)?$/){
  $filetype = "gff3";
}
elsif ($file =~ /\.g(t|f)f(.gz)?$/){
  $filetype = "gtf";
}
else{
  die "Did not recognize the file format\n";
}

my %trs;
if ($file =~ /\.gz$/){
  open (IN, "gunzip -dc $file |") or die "Can't open file $file: $!\n";
}
else{
  open (IN, $file) or die "Can't open file $file: $!\n";
}
while (<IN>){
  next if /^#/;
  chomp;
  my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split(/\t/);
  next if (!$chr or $chr =~ /_/ or $chr eq "chrM");
  next unless $type eq "exon";
  my %attribs;
  foreach my $att (split(/;/, $attributes)){
    my ($name, $value);
    if ($filetype eq "gff3"){
      ($name, $value) = split(/=/, $att);
    }
    elsif ($filetype eq "gtf"){
      $att =~ s/^\s*//g;
      ($name, $value) = split(/\s+/, $att);
      $value =~ s/\"//g;
    }
    $attribs{$name} = $value;
  }
  my $transcript_id = $attribs{'transcript_id'};
  if (!$transcript_id){
    die "No transcript id found: $_";
  }
  if ($cluster_by_gene_id){
    if ($attribs{'gene_id'}){
      $transcript_id = $attribs{'gene_id'};
    }
    else{
      die "No gene id found: $_";
    }
  }
  push(@{$transcript_data{$transcript_id}}, $_);
  if ($trs{$chr}{$transcript_id}){
    if ($start < $trs{$chr}{$transcript_id}{'start'}){
      $trs{$chr}{$transcript_id}{'start'} = $start;
    }
    if ($end > $trs{$chr}{$transcript_id}{'end'}){
      $trs{$chr}{$transcript_id}{'end'} = $end;
    }
  }
  else{
    $trs{$chr}{$transcript_id}{'start'} = $start;
    $trs{$chr}{$transcript_id}{'end'} = $end;
  }
}
close (IN);

foreach my $chr (sort keys %trs){
  foreach my $tid (sort {$trs{$chr}{$a}{'start'} <=> $trs{$chr}{$b}{'start'}} keys %{$trs{$chr}}){
    #print join("\t", $tid, $chr, $trs{$chr}{$tid}{'start'}, $trs{$chr}{$tid}{'end'})."\n";
  }
}

#Clusters
my @clusters;
foreach my $chr (sort keys %trs){
  foreach my $tid (sort {$trs{$chr}{$a}{'start'} <=> $trs{$chr}{$b}{'start'}} keys %{$trs{$chr}}){
    my $clustered = 0;
    if (scalar @clusters){
      foreach my $cluster (@clusters){
        if ($cluster->{'chr'} eq $chr and $cluster->{'start'} <= $trs{$chr}{$tid}{'end'} and $cluster->{end} >= $trs{$chr}{$tid}{'start'}){
          if ($cluster->{'start'} > $trs{$chr}{$tid}{'start'}){
            $cluster->{'start'} = $trs{$chr}{$tid}{'start'};
          }
          if ($cluster->{'end'} < $trs{$chr}{$tid}{'end'}){
            $cluster->{'end'} = $trs{$chr}{$tid}{'end'};
          }
          push(@{$cluster->{'transcripts'}}, $tid);
          $clustered = 1;
          last;
        }
      }
    }
    if (scalar @clusters == 0 or !$clustered){
      my %cluster;
      $cluster{'chr'} = $chr;
      $cluster{'start'} = $trs{$chr}{$tid}{'start'};
      $cluster{'end'} = $trs{$chr}{$tid}{'end'};
      push(@{$cluster{'transcripts'}}, $tid);
      push(@clusters, \%cluster);
    }
  }
}

unless ($no_shuffle){
  @clusters = shuffle(@clusters);
}
my $c = 0;
my $n = 0;
foreach my $cluster (@clusters){
  print join("\t", "Cluster #".$c++, scalar @{$cluster->{'transcripts'}}, $cluster->{'chr'}, $cluster->{'start'}, $cluster->{'end'}, join(",", @{$cluster->{'transcripts'}}))."\n";
  $n = $n + scalar @{$cluster->{'transcripts'}};
}
print "Total transcripts = $n\n";



#Write temporary annotation files
if ($tmpdir){
  make_path("$tmpdir/gxf", "$tmpdir/out");
}
my $number_of_files = 0;
my $fn = 0;

##by transcript id
if ($max_trs_per_file){
  my $tn = 0;
  foreach my $transcript_id (keys %transcript_data){
    $tn++;
    if (($tn > 1 and $tn % $max_trs_per_file == 1) or $max_trs_per_file == 1){
      close (FILE);
    }
    if ($tn == 1 or $tn % $max_trs_per_file == 1 or $max_trs_per_file == 1){
      $fn++;
      open (FILE, ">>$tmpdir/gxf/$tag.$fn.$filetype") or die "Can't open file $tmpdir/gxf/$tag.$fn.$filetype: $!";
    }
    foreach my $line (@{$transcript_data{$transcript_id}}){
      print FILE $line."\n";
    }
  }
  close (FILE);
}

##by clusters
elsif ($max_clusters_per_file){
  if ($max_num_files){
    if (!$max_clusters_per_file or ($max_clusters_per_file and $max_clusters_per_file < int(scalar(@clusters)/$max_num_files))){
      $max_clusters_per_file = POSIX::ceil(scalar(@clusters)/$max_num_files);
    }
  }
  my $cn = 0;
  foreach my $cluster (@clusters){
    $cn++;
    if (($cn > 1 and $cn % $max_clusters_per_file == 1) or $max_clusters_per_file == 1){
      close (FILE);
    }
    if ($cn == 1 or $cn % $max_clusters_per_file == 1 or $max_clusters_per_file == 1){
      $fn++;
      open (FILE, ">>$tmpdir/gxf/$tag.$fn.$filetype") or die "Can't open file $tmpdir/gxf/$tag.$fn.$filetype: $!";
    }
    foreach my $transcript_id (@{$cluster->{'transcripts'}}){
      foreach my $line (@{$transcript_data{$transcript_id}}){
        print FILE $line."\n";
      }
    }
  }
  close (FILE);  
}

if ($fn > $number_of_files){
  $number_of_files = $fn;
}
print "\nNumber of files = ".$number_of_files."\n";



my $dependency = "";
if ($dep_job){
  $dependency = "-w 'ended(\"$dep_job"."[\%I]\")'";
}

  my $command = <<COM;
  bsub -q short -M1500 -R"select[mem>1500] rusage[mem=1500]" $dependency -J "$tag\[1-$number_of_files\]%$job_limit" -oo $tmpdir/out/f.$tag.\%I.out \\
    perl $ENV{LOUTRE_WRITE}/scripts/load_gxf_in_loutre.pl \\
      -file $tmpdir/gxf/$tag.\%I.$filetype \\
      -dataset $dataset_name \\
COM

  $command .= " -author $author_name \\\n"             if $author_name;
  $command .= " -readseqdir $read_seq_dir \\\n"        if $read_seq_dir;
  $command .= " -source $source_info \\\n"             if $source_info;
  $command .= " -remark \"$remark\" \\\n"              if $remark;
  $command .= " -comp_pipe \\\n"                       if $use_comp_pipe_biotype;
  $command .= " -analysis $analysis_name \\\n"         if $analysis_name;
  $command .= " -tsource $tsource \\\n"                if $tsource;
  $command .= " -assembly $assembly_version \\\n"      if $assembly_version;
  $command .= " -no_check \\\n"                        if $no_artifact_check;
  $command .= " -no_NFV \\\n"                          if $no_NFV;
  $command .= " -no_CDS \\\n"                          if $do_not_add_cds;
  $command .= " -no_intron_check \\\n"                 if $no_intron_check;
  $command .= " -host_biotype $host_biotype \\\n"      if $host_biotype;
  $command .= " -max_ov_loc $max_overlapped_loci \\\n" if $max_overlapped_loci;
  $command .= " -filter_introns \\\n"                  if $filter_introns;
  $command .= " -platinum \\\n"                        if $platinum;
  $command .= " -chr $only_chr \\\n"                   if $only_chr;
  $command .= " -write \\\n"                           if $write;

  print "Sending command:\n$command\n\n";
  `$command`; #execute command

