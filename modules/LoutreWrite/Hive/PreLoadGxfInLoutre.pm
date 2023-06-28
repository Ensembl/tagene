package LoutreWrite::Hive::PreLoadGxfInLoutre;

use strict;
use warnings;

use parent 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB';
use File::Path qw(remove_tree make_path);
use File::Spec::Functions qw(catfile catdir);
use List::Util qw(shuffle);

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
  }
}


sub fetch_input {
  my ($self) = @_;
  my $file = $self->param('input_file');

  unless ($self->param_is_defined('cluster_by_gene_id')){
    $self->param('cluster_by_gene_id', 0);
  }
  
  #Read input annotation file
  my %transcript_data;
  my %seen;
  my $filetype;
  if ($file =~ /\.gff3(.gz)?$/){
    $self->param('filetype' ,"gff3");
  }
  elsif ($file =~ /\.g(t|f)f(.gz)?$/){
    $self->param('filetype', "gtf");
  }
  else{
    $self->throw("Did not recognize the file format\n");
  }
  
  my %trs;
  if ($file =~ /\.gz$/){
    open (IN, "gunzip -dc $file |") or $self->throw("Can't open file $file");
  }
  else{
    open (IN, $file) or $self->throw("Can't open file $file");
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
      if ($self->param('filetype') eq "gff3"){
        ($name, $value) = split(/=/, $att);
      }
      elsif ($self->param('filetype') eq "gtf"){
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
    if ($self->param('cluster_by_gene_id')){
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
  close (IN) or $self->throw("Can't close file $file");


  $self->param('transcript_coords', \%trs);
  $self->param('transcript_data', \%transcript_data);  
}


sub run {
  my ($self) = @_;

  my %trs = %{$self->param('transcript_coords')};
  
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
  
  unless ($self->param('no_shuffle')){
    @clusters = shuffle(@clusters);
  }
  my $c = 0;
  my $n = 0;
  foreach my $cluster (@clusters){
    $self->say_with_header(join("\t", "Cluster #".$c++, scalar @{$cluster->{'transcripts'}}, $cluster->{'chr'}, $cluster->{'start'}, $cluster->{'end'}, join(",", @{$cluster->{'transcripts'}})));
    $n = $n + scalar @{$cluster->{'transcripts'}};
  }
  $self->say_with_header("Total transcripts = $n");
  
  
  $self->output(\@clusters);
}



sub write_output {
  my ($self) = @_;

  my %transcript_data = %{$self->param('transcript_data')};

  my $output_dir = $self->param('output_dir');
  my $log_dir = $self->param('log_dir');
  my $tag = $self->param('tag');
  my $filetype = $self->param('filetype');

  #Write temporary annotation files
  #if ($tmpdir){
  #  make_path("$tmpdir/gxf", "$tmpdir/out");
  #}
  
  #Delete files from a previous run
  #system("rm -rf $output_dir/*.$filetype");
  if (-d $output_dir){
    remove_tree($output_dir, {safe => 1, keep_root => 1});
  }
  if (-d $log_dir){
    remove_tree($log_dir, {safe => 1, keep_root => 1});
  }

  #Calculate expected number of files
  my $expected_number = 0;
  if (my $max_trs_per_file = $self->param('max_trs_per_file')){
    $expected_number = POSIX::ceil(scalar(keys %transcript_data)/$max_trs_per_file);
  }
  elsif (my $max_clusters_per_file = $self->param('max_clusters_per_file')){
    if (my $max_num_files = $self->param('max_num_files')){
      $expected_number = $max_num_files;
    }
    else{
      $expected_number = POSIX::ceil(scalar(@{$self->output})/$max_clusters_per_file);
    }
  }

  my $MAX_NUM_FILES_PER_DIR = 1000;
  my $number_of_files = 0;
  my $fn = 0;
  
  ##by transcript id (ignore clusters)
  if (my $max_trs_per_file = $self->param('max_trs_per_file')){
    my $tn = 0;
    foreach my $transcript_id (keys %transcript_data){
      $tn++;
      if (($tn > 1 and $tn % $max_trs_per_file == 1) or $max_trs_per_file == 1){
        close (FILE);
      }
      if ($tn == 1 or $tn % $max_trs_per_file == 1 or $max_trs_per_file == 1){
        $fn++;
        #open (FILE, ">>$output_dir/$tag.$fn.$filetype") or self->throw("Can't open file $output_dir/$tag.$fn.$filetype");
        #insert one level of subdirectories if the number of files is too high
        #assuming one level is enough if the limit is ~1000 files per directory
        if ($expected_number > $MAX_NUM_FILES_PER_DIR){ 
          my $subdir = POSIX::floor($fn/$MAX_NUM_FILES_PER_DIR);
          $output_dir = catdir($self->param('output_dir'), $subdir);
          $log_dir = catdir($self->param('log_dir'), $subdir);
        }
        my $output_file = catfile($output_dir, "$tag.$fn.$filetype");
        if (-e $output_file){
          unlink $output_file; #remove file in case it exists from a previous run
        }
        unless (-d $output_dir){
          make_path($output_dir);
        }
        unless (-d $log_dir){
          make_path($log_dir);
        }        
        open (FILE, ">$output_file") or $self->throw("Can't open file $output_file");
      }
      foreach my $line (@{$transcript_data{$transcript_id}}){
        print FILE $line."\n";
      }
    }
    close (FILE);
  }
  
  ##by cluster
  elsif (my $max_clusters_per_file = $self->param('max_clusters_per_file')){
    if (my $max_num_files = $self->param('max_num_files')){
      #if a 'max_num_files' value is given, increase the 'max_clusters_per_file' value as required
      if (!$max_clusters_per_file or ($max_clusters_per_file and $max_clusters_per_file < int(scalar(@{$self->output})/$max_num_files))){
        $max_clusters_per_file = POSIX::ceil(scalar(@{$self->output})/$max_num_files);
      }
    }
    my $cn = 0;
    foreach my $cluster (@{$self->output}){
      $cn++;
      if (($cn > 1 and $cn % $max_clusters_per_file == 1) or $max_clusters_per_file == 1){
        close (FILE);
      }
      if ($cn == 1 or $cn % $max_clusters_per_file == 1 or $max_clusters_per_file == 1){
        $fn++;
        #insert one level of subdirectories if the number of files is too high
        #assuming one level is enough if the limit is ~1000 files per directory
        if ($expected_number > $MAX_NUM_FILES_PER_DIR){ 
          my $subdir = POSIX::floor($fn/$MAX_NUM_FILES_PER_DIR);
          $output_dir = catdir($self->param('output_dir'), $subdir);
          $log_dir = catdir($self->param('log_dir'), $subdir);
        }
        my $output_file = catfile($output_dir, "$tag.$fn.$filetype");
        if (-e $output_file){
          unlink $output_file; #remove file in case it exists from a previous run
        }
        unless (-d $output_dir){
          make_path($output_dir);
        }
        unless (-d $log_dir){
          make_path($log_dir);
        }  
        open (FILE, ">$output_file") or $self->throw("Can't open file $output_file");
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
  $self->say_with_header("\nNumber of files = ".$number_of_files);



}


1;




  

