#!/usr/bin/perl -w

#This script splits a large input file into manageable chunks, 
#which are written in individual temporary files for the moment, 
#and submits the individual jobs calling the load_gxf_in_loutre script.
#It takes the same options as load_gxf_in_loutre, which will be passed on to this, 
#and the maximum number of transcripts per chromosome for each file.
#Usage: perl preload_gxf_in_loutre.pl -file data.gtf -n 1000 -dataset human_test ...


use strict;
use warnings;
use Getopt::Long;
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
my $max_n_tr;
my $tag;
my $dep_job;
my $tmpdir;
my $read_seq_dir;
my $cluster_by_gene_id;
my $job_limit = 20; #Max number of simultaneously running jobs in a job array

&GetOptions(
            'file=s'            => \$file,
            'dataset=s'         => \$dataset_name,
            'author=s'          => \$author_name,
            'source=s'          => \$source_info, #file with otter columns and read names supporting each transcript model
            'remark=s'          => \$remark, #remark to be added to every transcript model
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
            'max=i'             => \$max_n_tr,
            'tag=s'             => \$tag,
            'dep_job=s'         => \$dep_job,
            'tmpdir=s'          => \$tmpdir,
            'readseqdir=s'      => \$read_seq_dir,
            'by_gene_id!'       => \$cluster_by_gene_id,
            'job_limit=i'       => \$job_limit,
            );


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
  next if $type eq "gene";
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
}
close (IN);


#Write temporary annotation files
if ($tmpdir and !(-d $tmpdir)){
  mkdir($tmpdir);
}
my $number_of_files = 0;
my $tn = 0;
my $fn = 0;
foreach my $transcript_id (keys %transcript_data){
  $tn++;
  if (($tn > 1 and $tn % $max_n_tr == 1) or $max_n_tr == 1){
    close (FILE);
  }
  if ($tn == 1 or $tn % $max_n_tr == 1 or $max_n_tr == 1){
    $fn++;
    open (FILE, ">>$tmpdir/$tag.$fn.$filetype") or die "Can't open file $tmpdir/$tag.$fn.$filetype: $!";
  }
  foreach my $line (@{$transcript_data{$transcript_id}}){
    print FILE $line."\n";
  }
}
close (FILE);

if ($fn > $number_of_files){
  $number_of_files = $fn;
}


my $dependency = "";
if ($dep_job){
  $dependency = "-w 'ended(\"$dep_job"."[\%I]\")'";
}

  my $command = <<COM;
  bsub -q short -M1000 -R"select[mem>1000] rusage[mem=1000]" $dependency -J "$tag\[1-$number_of_files\]%$job_limit" -oo $tmpdir/f.$tag.\%I.out \\
    perl $ENV{LOUTRE_WRITE}/scripts/load_gxf_in_loutre.pl \\
      -file $tmpdir/$tag.\%I.$filetype \\
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


__END__
for (my $i=1; $i<=$number_of_files; $i++){
  my $dependency = "";
  if ($dep_job){
    $dependency = "-w 'ended(\"$dep_job"."[$i]\")'";
  }
#  if ($i == 1){
#    if ($dep_job){
#      $dependency = "-w 'ended(\"$dep_job"."[*]\")'";
#    }
#  }
#  if ($i > 1){
#    $dependency = "-w 'ended(\"$tag.".($i-1)."[*]\")'";
#  }
  my $command = <<COM;
  bsub -M3000 -R"select[mem>3000] rusage[mem=3000]" $dependency -J "$tag\[$i\]%$job_limit" -oo $tmpdir/f.$tag.\%I.out \\
    perl $ENV{LOUTRE_WRITE}/scripts/load_gxf_in_loutre.pl \\
      -file $tmpdir/$tag.$i.$filetype \\
      -dataset $dataset_name \\
      -readseqdir $read_seq_dir \\
COM


  $command .= " -author $author_name \\\n"             if $author_name;
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
}

