#!/usr/bin/env perl

#Filters out redundant reads from a file containing transcript ids and supporting reads by library, eg:
# transcript_1	library_1 : read_01, read_02, read_03, read_04 ; library_2 : read_5, read_6 ;

#Requires 
#a) input transcript source file,
#b) Ensembl core-like database with transcripts,
#c) read libraries (BAM files converted to BED format)

#Usage: perl filter_source_info.pl -in source_info.txt -host db_host -port db_port -user db_user -dbname db_name -analysis analysis_name -readfiles bed_file_1[, bed_file_2, ...] -out filtered_source_info.txt


use strict;
use warnings;
use Getopt::Long;
use List::Util qw(min max);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use File::Basename;



my $infile;
my $host = "mysql-ens-havana-prod-1";
my $port = 4581;
my $user = "ensro";
my $dbname;
my $analysis;
my $read_bed_path;
my $read_files;
my $read_list;
my $outfile;

&GetOptions(
            'in=s'          => \$infile,
            'host=s'        => \$host,
            'port=s'        => \$port,
            'user=s'        => \$user,
            'dbname=s'      => \$dbname,
            'analysis=s'    => \$analysis,
            'readbedpath=s' => \$read_bed_path,
            'readfiles=s'   => \$read_files,
            'readlist=s'    => \$read_list,
            'out=s'         => \$outfile,
           );


#Fetch database transcripts
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host => $host,
    -port => $port,    
    -user => $user,
    -dbname => $dbname,
);
my $sa = $db->get_SliceAdaptor();
my $ta = $db->get_TranscriptAdaptor();
my @transcripts = @{$ta->fetch_all_by_logic_name($analysis)};
print "Fetch all transcripts DONE\n";


#Read list of selected reads (if any)
my %selected_reads;
if ($read_list){
  open (RL, $read_list) or die "Cant open file $read_list\n";
  while (<RL>){
    chomp;
    my $read = $_;
    $read =~ s/\//_/g;
    $selected_reads{$read} = 1;
  }
  close (RL);
}
print "Select reads DONE\n";


#First pass read of source file
#Needed if the number of reads in the BED files is too large,
#so as not to store all BED lines in memory
my %all_source_info; my $x=0;
open (IN, $infile) or die "Cant open file $infile\n";
while (<IN>){
  chomp;
  my ($transcript_id, $source_info) = read_source_line($_);
  foreach my $library (keys %$source_info){
    foreach my $read (keys %{$source_info->{$library}}){
      #$read =~ s/\//_/g; #Otter issue with '/' in read names
      if ($read_list){
        next unless $selected_reads{$read};
      }
      $all_source_info{$library}{$read} = 1; $x++;
    }
  }
}
close (IN);
print "First pass source file DONE\n$x reads\n";


#Store read BED files in memory
unless (-d $read_bed_path){
  die "Can't find directory $read_bed_path";
}
my @read_files;
if ($read_files){
  @read_files = split(",", $read_files);
}
else{
  opendir (my($dh), $read_bed_path) or die "Couldn't open $read_bed_path: $!";
  @read_files = sort grep {$_ =~ /\.bed(.gz)?$/} readdir $dh;
  closedir $dh;
}
@read_files = map {$read_bed_path."/".$_} @read_files;
my $read_libraries = read_bed_files(\@read_files, \%all_source_info);
print "Read BED files DONE\n";

open (OUT, ">$outfile") or die "Can't open $outfile\n";

open (IN, $infile) or die "Cant open file $infile\n";
while (<IN>){
  chomp;
  my ($transcript_id, $source_info) = read_source_line($_);

  my @read_beds;
  my $n_reads;
  foreach my $library (keys %$source_info){
    foreach my $read (keys %{$source_info->{$library}}){
      if ($read_list){
        next unless $selected_reads{$read};
      }
      my $read_lines = find_bed_lines($read, $library);
      push(@read_beds, @{$read_lines}); 
      $n_reads++;  
    }  
  }
  print "\n$transcript_id ($n_reads reads)\n";
  
  my $transcript = get_transcript($transcript_id, \@transcripts);
  
  #For each read, 
  # -check there is only exon-exon and intron-intron overlap
  # -find % overlap with transcript
  my %overlaps;
  foreach my $read_bed (@read_beds){
    my $read_tr = make_transcript_from_bed_line($read_bed);
    if (are_compatible($read_tr, $transcript)){
      my $overlap_percent = get_overlap_percent($read_tr, $transcript);
      push(@{$overlaps{$overlap_percent}}, $read_tr);
    }  
  }
  
  #Sort reads by % overlap
  #Map reads onto transcript until every nucleotide is covered, 
  # discarding reads that don't provide extra coverage
  my $transcript_hash = trans_coords_to_hash($transcript);
  my %filtered_source;
  my $filtered_in = 0;
  my $c= 0;
  foreach my $ovpc (sort {$b <=> $a} keys %overlaps){
    foreach my $read_tr (@{$overlaps{$ovpc}}){
      #print $read_tr->stable_id."\t".$ovpc."\n";
      ($transcript_hash, $filtered_in) = map_read_to_transcript($read_tr, $transcript_hash); 
      if ($filtered_in){
        $filtered_source{$read_tr->stable_id} = 1;
        print ++$c.": ".$read_tr->stable_id." - ".sprintf("%.2f", $ovpc)."% overlap - ".coverage($transcript_hash)."% accumulated coverage\n";
      }
    }
  }
  
  #Go through original source info and delete non-essential reads
  foreach my $library (keys %{$source_info}){
    foreach my $read (keys %{$source_info->{$library}}){
      unless ($filtered_source{$read}){
        delete($source_info->{$library}{$read});
      }
    }
  }
  #Print filtered source info
  print OUT $transcript_id."\t";
  foreach my $library (keys %$source_info){
    if (scalar keys %{$source_info->{$library}}){
      print OUT $library." : ".join(", ", keys %{$source_info->{$library}})." ; ";
    }
  }
  print OUT "\n";
}
close (IN);
close (OUT);



######


sub read_bed_files {
  my ($files, $s_info) = @_;
#foreach my $a (keys %{$s_info}){
#  foreach my $b (keys %{$s_info->{$a}}){
#    print "$a  $b\n";
#  }
#}
  
  my %libraries;
  foreach my $file (@$files){
    my $lib_name;
    if ($file =~ /\.bed$/){
      $lib_name = basename($file, ".bed");
      open (BED, $file) or die "Can't open file $file\n";
    }
    else{
      $lib_name = basename($file, ".bed.gz"); print "LIB=".$lib_name."\n";
      open (BED, "gunzip -dc $file |") or die "Can't open file $file\n";
    }
    while (<BED>){
      chomp;
      my @fs = split(/\t/);
      my $read_name = $fs[3]; #print "Finding $read_name\n";
      if ($s_info->{$lib_name}{$read_name}){
        push(@{$libraries{$lib_name}{$read_name}}, $_); #print "Found $read_name\n";
      }
    }
    close (BED);
    if ($libraries{$lib_name}){
      print "Stored ".scalar(keys %{$libraries{$lib_name}})." reads from $lib_name\n";
    }
  }
  return \%libraries;
}


sub read_source_line {
  my $line = shift;
  my ($transcript_id, $sources) = split(/\t/, $line);
  my %source_info;
  foreach my $item (split(/ ; /, $sources)){
    my ($library, $reads) = split(/ : /, $item);
    foreach my $read (split(/,/, $reads)){
      $read =~ s/\//_/g; #Otter issue with '/' in read names
      $source_info{$library}{$read} = 1;  
    }  
  }
  return ($transcript_id, \%source_info);
}  
  
  
sub find_bed_lines {
  my ($read, $library) = @_;
  #There may be more than one read with the same name (multimappers?)
  my @read_lines;
  if ($read_libraries->{$library}{$read}){
    @read_lines = @{$read_libraries->{$library}{$read}};
  }
  return \@read_lines;
}


sub get_transcript {
  my ($tid, $trs) = @_;
  foreach my $tr (@$trs){
    if ($tr->stable_id eq $tid){
      return $tr;
    }
  }
}


sub make_transcript_from_bed_line {
  my $line = shift;
  my $transcript = new Bio::EnsEMBL::Transcript();
  my ($chr, $tstart, $tend, $name, $score, $strand, $cstart, $cend, $rgb, $nexons, $elengths, $estarts) = split(/\t/, $line);
  $chr =~ s/^chr//;
  $chr = "MT" if $chr eq "M"; #Some multimappers are aligned to chrM
  my $slice;
  $slice = $sa->fetch_by_region("toplevel", $chr);
  unless ($slice){
    $chr = "chr".$chr;
    foreach my $sl (@{$sa->fetch_all('toplevel')}){
      foreach my $synonym (@{$sl->get_all_synonyms('UCSC')}){
        if ($synonym->name eq $chr){
          $slice = $sl;
          last;
        }
      }
    }
  }
    
  my @exon_starts = split(/,/, $estarts);
  my @exon_lengths = split(/,/, $elengths);
  for (my $i=0; $i<$nexons; $i++){
    my $exon = new Bio::EnsEMBL::Exon(
      -START     => $tstart + $exon_starts[$i] + 1,
      -END       => $tstart + $exon_starts[$i] + $exon_lengths[$i],
      -STRAND    => ($strand eq "+" ? 1 : -1),
      -SLICE     => $slice,
    );
    $transcript->add_Exon($exon);  
  }
  $transcript->stable_id($name);
  return $transcript;
}


sub are_compatible {
  my ($read_tr, $db_tr) = @_;
  #The read must fit in the merged model (only exon-exon and intron-intron overlap)
  #Overhangs are allowed: some some reads were truncated because of small terminal exons, so the transcript models
  #can be shorter than the reads (this script uses the original reads)
  print $read_tr->slice->seq_region_name."\n";
  print $db_tr->slice->seq_region_name."\n";
  if ($read_tr->slice->seq_region_name ne $db_tr->slice->seq_region_name){
    return 0;
  }
  my $ex_ovlp = 0;
  foreach my $exon1 (@{$read_tr->get_all_Exons}){
    #Exon-exon overlap: required
    foreach my $exon2 (@{$db_tr->get_all_Exons}){
      if ($exon1->seq_region_start <= $exon2->seq_region_end and $exon1->seq_region_end >= $exon2->seq_region_start){
        $ex_ovlp = 1;
      }
    }
    #Exon-intron overlap: not allowed
    foreach my $intron2 (@{$db_tr->get_all_Introns}){
      if ($exon1->seq_region_start <= $intron2->seq_region_end and $exon1->seq_region_end >= $intron2->seq_region_start){
        return 0;
      }
    }
    #Exon overhang: not allowed
    #if ($exon1->seq_region_start < $db_tr->seq_region_start or $exon1->seq_region_end > $db_tr->seq_region_end){
     # return 0;
    #}
  }
  #Intron-exon overlap: not allowed
  foreach my $intron1 (@{$read_tr->get_all_Introns}){
    foreach my $exon2 (@{$db_tr->get_all_Exons}){
      if ($intron1->seq_region_start <= $exon2->seq_region_end and $intron1->seq_region_end >= $exon2->seq_region_start){
        return 0;
      }
    }
  }
  if ($ex_ovlp == 1){
    return 1;
  }
  return 0;
}


sub get_overlap_percent {
  my ($read_tr, $db_tr) = @_;
  my $overlap_length = 0;
  foreach my $exon2 (@{$db_tr->get_all_Exons}){
    foreach my $exon1 (@{$read_tr->get_all_Exons}){
      if ($exon1->seq_region_start <= $exon2->seq_region_end and $exon1->seq_region_end >= $exon2->seq_region_start){
        $overlap_length += min($exon1->seq_region_end, $exon2->seq_region_end) -  max($exon1->seq_region_start, $exon2->seq_region_start) + 1;
      }
    }
  }
  return $overlap_length/$db_tr->length*100;
}


sub trans_coords_to_hash {
  my $tr = shift;
  my %hash;
  foreach my $exon (@{$tr->get_all_Exons}){
    for (my $i=$exon->seq_region_start; $i<=$exon->seq_region_end; $i++){
      $hash{$i} = 0;
    }
  }
  return \%hash;
}


sub map_read_to_transcript {
  my ($read_tr, $tr_hash) = @_;
  my $valid = 0;
  my $has_new_pos;
  foreach my $exon (@{$read_tr->get_all_Exons}){
    for (my $i=$exon->seq_region_start; $i<=$exon->seq_region_end; $i++){
      if ($tr_hash->{$i} == 0){
        $tr_hash->{$i} = 1;
        $has_new_pos++;
      }
    }
  }
  if ($has_new_pos){
    $valid = 1;
  }
  return ($tr_hash, $valid);
}


sub coverage {
  my $hash = shift;
  my $t = 0;
  my $c = 0;
  foreach my $key (keys %{$hash}){
    $t++;
    if ($hash->{$key} == 1){
      $c++;
    } 
  }
  return sprintf("%.2f", $c/$t*100);
}


