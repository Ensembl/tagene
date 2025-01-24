#!/usr/bin/env perl

#For a given read supporting a transcript model, check if there are misalignments around any of the splice sites.
#Search for mismatches and indels in a window of X nucleotides before the donor site and after the acceptor site.
#Input:
# -GTF file
# -BAM files
# -List of reads supporting each transcript model
# -List of supporting reads and their number of introns (optional)
# -Parameters to be decided: window size around splice site and max number of allowed errors 
#Output:
# -Reads accepted/rejected for each transcript model, detailing failing introns in a log file

use strict;
use warnings;
use Getopt::Long;
use Bio::DB::HTS;
use Net::FTP;
use List::Util qw(min max uniq);

$|=1;

my $gtf_file;
my $only_chr;
my $bam_dir;
my $genome_fasta;
my $id_list;
my $read_support_list;
my $read_intron_number_list;
my $window_size = 10;
my $max_errors = 3;
my $outfile;
my $assembly_report_url;
my $keep_it_simple;
my $print_whole_alignment;

&GetOptions(
            'gtf=s'        => \$gtf_file,
            'bamdir=s'     => \$bam_dir,
            'fasta=s'      => \$genome_fasta,
            'readsup=s'    => \$read_support_list,
            'readintnum=s' => \$read_intron_number_list,
            'ids=s'        => \$id_list,
            'chr=s'        => \$only_chr,         
            'window=i'     => \$window_size,
            'maxerr=i'     => \$max_errors,
            'asreport=s'   => \$assembly_report_url,
            'outfile=s'    => \$outfile,
            'simple!'      => \$keep_it_simple,
            'print_ali!'   => \$print_whole_alignment,
           );



my $usage =  <<_USAGE_;

This script looks for RNA-seq-based transcripts containing splice sites supported by misaligned reads.
For each read supporting a spliced transcript, it searches a window of X nucleotides before the donor site
and after the acceptor site for mismatches and indels, assigning a score based on the distance from these
to the splice site. If the score is higher than a given threshold, the supporting read is considered invalid.
If a transcript is left without any valid supporting read, it is tagged for filtering.

perl check_splice_site_alignments.pl -gtf ANNOTATION_FILE -bamdir PATH_TO_BAM_DIR -fasta GENOME_FASTA_FILE -readsup READ_SUPPORT_FILE
 -gtf            annotation file in GTF format - only 'exon' lines and 'transcript_id' attributes are essential
 -bamdir         directory containing the BAM file(s)
 -fasta          genome sequence file in Fasta format
 -readsup        file with supporting reads for each transcript: tab-delimited file with columns 1) transcript name, 2) BAM file name (without .bam suffix) and 3) comma-separated list of read names
 -outfile        output file
 
  OPTIONAL PARAMETERS:
  -readintnum    list of number of introns of supporting reads: tab-delimited file with columns 1) BAM file name (without .bam suffix), 2) read name, 3) number of introns. It can help speed up the analysis
  -ids           only check the transcripts having the given IDs
  -chr           only check the annotation in the given chromosome
  -window        window size before donor splice site or after acceptor splice site (default = 10)
  -maxerr        maximum allowed number of misaligned nucleotides in the splice site window (default = 3)
  -asreport      genome assembly report from NCBI - only needed if chromosome or scaffold names differ between GTF and BAM files, eg they are synonymns
  -simple        stop checking a splice junction after finding a single valid supporting read
  -print_ali     print the read-genome alignment in the output file

_USAGE_

die $usage unless ($gtf_file and $bam_dir and $genome_fasta and $outfile);



#Connect to a remote or local directory and get a list of BAM files
my @bam_files;
if ($bam_dir =~ /^http/){
  my ($server, $dir) = $bam_dir =~ /^https:\/\/([a-z.]+)\/(.+)$/;
  my $ftp = Net::FTP->new($server, Debug => 0) or die "Cannot connect to $server: $@";
  $ftp->login("anonymous",'-anonymous@') or die "Cannot login ", $ftp->message;
  $ftp->cwd($dir) or die "Cannot change working directory ", $ftp->message;
  @bam_files = grep {$_ =~ /\.bam$/} $ftp->ls();
}
else{
  opendir (my($dh), $bam_dir) or die "Couldn't open $bam_dir: $!";
  @bam_files = sort grep {$_ =~ /\.bam$/} readdir $dh;
  closedir $dh;
}
my %bam_hts;
foreach my $bam (@bam_files){
  my $sam = Bio::DB::HTS->new( -bam => "$bam_dir/$bam", -fasta => $genome_fasta );
  $bam_hts{$bam} = $sam;
}


#Read list of supporting reads
my %supporting_reads;
open (RSL, $read_support_list) or die "Can't open file $read_support_list: $!";
while (<RSL>){
  chomp;
  #Expected format: tab-delimited file with transcript name, BAM file name (without .bam suffix) and comma-separated list of read names
  my ($tid, $file_name, $read_names) = split(/\t/);
  foreach my $read_name (split(/,/, $read_names)){
    push(@{$supporting_reads{$tid}}, $read_name."\t".$file_name);
  }
}
close (RSL);


#Read list of number of introns of supporting reads
my %read_intron_numbers;
if ($read_intron_number_list){
  open (RIN, $read_intron_number_list) or die "Can't open file $read_intron_number_list: $!";
  while (<RIN>){
    chomp;
    my ($sample_name, $read_name, $n_read_introns) = split(/\t/);
    $read_intron_numbers{$sample_name}{$read_name} = $n_read_introns;
  }
  close (RIN);
}


#Get UCSC scaffold synonyms
my %ucsc_synonyms;
if ($assembly_report_url){
  my $content = get($assembly_report_url);
  foreach my $line (split("\n", $content)){
    $line =~ s/\r//;
    next if $line =~ /^#/;
    my @cols = split("\t", $line);
    if ($cols[1] =~ /scaffold/){
      my $genbank = $cols[4];
      my $ucsc = $cols[9];
      $ucsc_synonyms{$genbank} = $ucsc;
    }
  }
}


#Read transcript id list
my @tids;
my %id_list;
if ($id_list){
  open (IN, $id_list) or die "Can't open $id_list: $!";
  while (<IN>){
    chomp;
    push(@tids, $_);
    $id_list{$_} = 1;
  }
  close (IN);
}


#Read GTF file
my %exon_coords_by_tid;
if ($gtf_file =~ /\.gz$/){
  open (GTF, "zcat $gtf_file |") or die "Can't open $gtf_file: $!";
}
else{
  open (GTF, "cat $gtf_file |") or die "Can't open $gtf_file: $!";
}
while (<GTF>){
  next if /^#/;
  chomp;
  #Read exon lines (do not assume that there are transcript lines)
  my @cols = split("\t");
  if ($cols[2] eq "exon"){
    if ($only_chr){
      unless ($cols[0] eq $only_chr or "chr".$cols[0] eq $only_chr or $cols[0] eq "chr".$only_chr){
        next;
      }
    }
    if ($cols[8] =~ /transcript_id \"(\S+)\";/){
      my $tid = $1;
      #Filter by ID list if any
      if ($id_list){
        next unless $id_list{$tid};
      }
      else{
        push(@tids, $tid);
      }
      #Store exon coordinates
      my ($chr, $start, $end, $strand) = @cols[0,3,4,6];
      $exon_coords_by_tid{$tid}{'chr'} = $chr;
      $exon_coords_by_tid{$tid}{'strand'} = $strand;
      push(@{$exon_coords_by_tid{$tid}{'starts'}}, $start);
      push(@{$exon_coords_by_tid{$tid}{'ends'}}, $end);  
    }
    else{
      die "No transcript_id field!";
    }
  }
}
close (GTF);


#Calculate intron coordinates for each transcript
my %intron_coords_by_tid;
foreach my $tid (keys %exon_coords_by_tid){
  if (scalar @{$exon_coords_by_tid{$tid}{'starts'}} > 1){
    $intron_coords_by_tid{$tid}{'chr'} = $exon_coords_by_tid{$tid}{'chr'};
    $intron_coords_by_tid{$tid}{'strand'} = $exon_coords_by_tid{$tid}{'strand'};
    if ($intron_coords_by_tid{$tid}{'strand'} eq "+"){
      my @starts = sort {$a<=>$b} @{$exon_coords_by_tid{$tid}{'starts'}};
      my @ends = sort {$a<=>$b} @{$exon_coords_by_tid{$tid}{'ends'}};
      @{$intron_coords_by_tid{$tid}{'donors'}} = map {$_ + 1} @ends[0..$#ends-1];
      @{$intron_coords_by_tid{$tid}{'acceptors'}} = map {$_ - 1} @starts[1..$#starts];
    }
    elsif ($intron_coords_by_tid{$tid}{'strand'} eq "-"){
      my @starts = sort {$b<=>$a} @{$exon_coords_by_tid{$tid}{'starts'}};
      my @ends = sort {$b<=>$a} @{$exon_coords_by_tid{$tid}{'ends'}};
      @{$intron_coords_by_tid{$tid}{'donors'}} = map {$_ - 1} @starts[0..$#starts-1];
      @{$intron_coords_by_tid{$tid}{'acceptors'}} = map {$_ + 1} @ends[1..$#ends];
    }
  }
}
#print Dumper(\%intron_coords_by_tid); 


open (OUT, ">$outfile") or die "Can't open $outfile: $!";
print OUT "#".join("\t", "transcript_ID",
                        "has_valid_support",
                        "number_of_introns",
                        "number_of_valid_reads",
                        "number_of_misaligned_reads",
                        "valid_read_names",
                        "misaligned_read_names",
                  )."\n";
open (LOG, ">$outfile.log") or die "Can't open $outfile.log: $!";


#Loop through transcripts
#foreach my $tid (keys %intron_coords_by_tid){
foreach my $tid (uniq @tids){
  print LOG "\n##".$tid."\n";
  my $timestamp = localtime(time);
  print LOG $timestamp."\n";
  my $transcript_has_valid_support = 0;
  my @valid_reads;
  my @misaligned_reads;
  if ($intron_coords_by_tid{$tid}{'donors'}){
    #Loop through supporting reads of the transcript
    READ:for (@{$supporting_reads{$tid}}){
      my ($read_name, $sample_name) = split(/\t/, $_);
      print LOG "\n#".$read_name." ($sample_name)\n";
      my $read_has_misaligned_splice_site = 0;

      #Remove read copy number that was appended to the read name in the GTF file
      #(when multiple reads/alignments under the same name are present in the BAM file)
      $read_name =~ s/_\d{1,2}$//;

      #If the read's number of introns is available from a list, compare it with the number of introns of the transcript
      #This saves the time that extracting the read from the BAM file would take
      if ($read_intron_numbers{$sample_name} and $read_intron_numbers{$sample_name}{$read_name}){
        my $n_read_introns = $read_intron_numbers{$sample_name}{$read_name};
        if ($n_read_introns != @{$intron_coords_by_tid{$tid}{'donors'}}){
          print LOG "  => Different number of introns: READ(L) = ".$n_read_introns." - TRANSCRIPT = ".scalar(@{$intron_coords_by_tid{$tid}{'donors'}})."\n\n";
          push(@misaligned_reads, $read_name);
          next;
        }
      }

      #Get the read alignment from the BAM file - search by both location and name (searching by name alone is very slow)
      my ($read_ali) = $bam_hts{$sample_name.".bam"}->get_features_by_location(
                                                      -seq_id => $ucsc_synonyms{$exon_coords_by_tid{$tid}{'chr'}} || $exon_coords_by_tid{$tid}{'chr'},
                                                      -start  => min(@{$exon_coords_by_tid{$tid}{'starts'}}),
                                                      -end    => max(@{$exon_coords_by_tid{$tid}{'ends'}}),
                                                      -name   => $read_name
                                                    );

      unless ($read_ali){
        print LOG "\nNo alignment for $read_name - $sample_name at ".($ucsc_synonyms{$exon_coords_by_tid{$tid}{'chr'}} || $exon_coords_by_tid{$tid}{'chr'}).":".min(@{$exon_coords_by_tid{$tid}{'starts'}})."-".max(@{$exon_coords_by_tid{$tid}{'ends'}})."\n\n";
        next;
      }
      #Skip this read if it's unspliced (i.e. there is no intron in the CIGAR string)
      unless ($read_ali->cigar_str =~ /\dN\d/){
        print LOG "\nUnspliced read: $read_name - $sample_name - CIGAR = ".$read_ali->cigar_str."\n\n";
        next;
      }

      #Skip this read if the number of introns doesn't match the number of transcript's introns
      my $n_read_introns = scalar(() = $read_ali->cigar_str =~ /\dN\d/g);
      unless ($n_read_introns == @{$intron_coords_by_tid{$tid}{'donors'}}){
        print LOG "  => Different number of introns: READ = ".$n_read_introns." - TRANSCRIPT = ".scalar(@{$intron_coords_by_tid{$tid}{'donors'}})."\n\n";
        push(@misaligned_reads, $read_name);
        next;
      }

      #Skip this read if the query sequence is not available
      unless ($read_ali->query->seq->seq =~ /[ATGC]/){
        print LOG "  => Alignment does not have query sequence\n\n";
        push(@misaligned_reads, $read_name);
        next;
      }        

      #Skip this read if it is not a primary alignment
      if ($read_ali->get_tag_values("FLAGS") =~ /NOT_PRIMARY/){
        print LOG "  => Alignment is not primary\n\n";
        push(@misaligned_reads, $read_name);
        next;
      }   

      #Using the padded alignment, find the genomic coordinates of mismatches, insertions and deletions
      #Store alignment output of each genomic position (exonic sequence only)
      my ($ref, $matches, $query) = $read_ali->padded_alignment;
      if ($print_whole_alignment){
        print LOG join("\n", "R: ".$ref, "   ".$matches, "Q: ".$query)."\n\n";
      }
      my $gpos = $read_ali->pos; print LOG "POS:".$gpos."\n";
      my @ref = split("", $ref);
      my @query = split("", $query);
      my %diff;
      for (my $i = 0; $i < $#ref; $i++){
        unless ($query[$i]){
          print LOG "MISSING SEQ AT $i\n";
        }
        if ($ref[$i] =~ /[ATGCN]/){
          $gpos++;
          if ($query[$i] eq "-"){
            $diff{$gpos} = "D"; #deletion
          }
          elsif ($ref[$i] ne $query[$i]){
            $diff{$gpos} = "M"; #match/mismatch
          }
        }
        elsif ($ref[$i] eq "-"){
          if ($ref[$i] ne $query[$i]){
            $diff{$gpos} .= "I"; #insertion
          }
          else{
            $gpos++; #intron
          }
        }
        #print LOG $gpos."->".$diff{$gpos}."\n";
      }

      #Using the CIGAR string, find intron coordinates
      #and get the number of mismatches/indels in the windows
      #before the donor splice site and after the acceptor splice site.
      #(Intron coordinates can also be found using the padded alignment)
      my $cigar = $read_ali->cigar_array; #arrayref
      print LOG $read_ali->cigar_str."\n\n";

      my %read_introns;
      
      my $pos = $read_ali->pos; #genomic coordinate of start of alignment
      my $rel = 0; #position relative to the start of the alignment
      my $intron_counter = 0;
      foreach my $bit (@$cigar){
        unless ($bit->[0] eq "H"){
          $rel += $bit->[1]; #STC593_SLO97_BN_IsoSeq_Pool_4_HQ_transcript/33424_7
        }
        if ($bit->[0] eq "M" or $bit->[0] eq "D"){ #match/mismatch or deletion
          $pos += $bit->[1];
        }
        elsif ($bit->[0] eq "N"){ #intron
          #Assess alignment quality around splice sites
          my $intron_ok = 1;
          my $istart = $pos + 1;
          $pos += $bit->[1];
          my $iend = $pos;

          #store read intron coordinates
          push(@{$read_introns{'starts'}}, $istart);
          push(@{$read_introns{'ends'}}, $iend);

          #check that intron is present in transcript model
          unless ($intron_coords_by_tid{$tid}{'strand'} eq "+" and scalar(grep {$_ == $istart} @{$intron_coords_by_tid{$tid}{'donors'}}) and scalar(grep {$_ == $iend} @{$intron_coords_by_tid{$tid}{'acceptors'}}) or
                  $intron_coords_by_tid{$tid}{'strand'} eq "-" and scalar(grep {$_ == $istart} @{$intron_coords_by_tid{$tid}{'acceptors'}}) and scalar(grep {$_ == $iend} @{$intron_coords_by_tid{$tid}{'donors'}})
          ){
            print LOG "Can't find intron $istart-$iend ($read_name) in transcript $tid\n";
          }

          #count number of sequence differences in windows around splice sites   
          my $n_diff_start = 0;
          my $closest_diff_start = $window_size + 1; #default value for closest mismatch to donor splice site
          my $score_start = 0;
          for (my $i = $istart-$window_size; $i < $istart; $i++){
            print LOG $i." ".($diff{$i}||" ")."\n";
            if ($diff{$i}){
              $n_diff_start++; #print LOG "AAA $i AAA\n";
              $closest_diff_start = $istart - $i;
              $score_start += $window_size - ($istart - $i) + 1; #misalignment score, directly proportional to the distance to the splice site
              if ($diff{$i} =~ /I/){
                print LOG $score_start."  ".$istart."  ".$i." ".$diff{$i}." ".scalar(() = $diff{$i} =~ /I/g)."\n";
                $score_start += scalar(() = $diff{$i} =~ /I/g) - 1; #increase misalignment score by adding insertion size
              }
            }
          }
          if ($closest_diff_start > $window_size){
            $closest_diff_start .= "+";
          }
          
          my $n_diff_end = 0;
          my $closest_diff_end = $window_size + 1; #default value for closest mismatch to acceptor splice site
          my $score_end = 0;
          for (my $i = $iend+$window_size; $i > $iend; $i--){
            print LOG $i." ".($diff{$i}||" ")."\n";
            if ($diff{$i}){
              $n_diff_end++; #print LOG "BBB $i BBB\n";
              $closest_diff_end = $i - $iend;
              $score_end += $window_size - ($i - $iend) + 1; #misalignment score, directly proportional to the distance to the splice site
              if ($diff{$i} =~/I/){
                print LOG $score_start."  ".$istart."  ".$i." ".$diff{$i}." ".scalar(() = $diff{$i} =~ /I/g)."\n";
                $score_end += scalar(() = $diff{$i} =~ /I/g) - 1; #increase misalignment score by adding insertion size
              }
            }
          }
          if ($closest_diff_end > $window_size){
            $closest_diff_end .= "+";
          }
          #unless ($n_diff_start <= $max_errors and $n_diff_end <= $max_errors){
          if ($score_start > $window_size or $score_end > $window_size){
            $intron_ok = 0;
            $read_has_misaligned_splice_site = 1;
            push(@misaligned_reads, $read_name);
          }
          
          #Print intron coordinates and score
          print LOG join("; ", "INTRON#".++$intron_counter."  ".$intron_coords_by_tid{$tid}{'chr'}.":$istart-$iend:".$intron_coords_by_tid{$tid}{'strand'},
                                "N_DIFFS:".($intron_coords_by_tid{$tid}{'strand'} eq "+" ? "donor=$n_diff_start,acceptor=$n_diff_end" : "donor=$n_diff_end,acceptor=$n_diff_start"),
                                "CLOSEST:".($intron_coords_by_tid{$tid}{'strand'} eq "+" ? "donor=$closest_diff_start,acceptor=$closest_diff_end" : "donor=$closest_diff_end,acceptor=$closest_diff_start"),
                                "SCORE:".($intron_coords_by_tid{$tid}{'strand'} eq "+" ? "donor=$score_start,acceptor=$score_end" : "donor=$score_end,acceptor=$score_start"),
                                $intron_ok ? "OK" : "FAILED",
                        )."\n";

          #Print sequence alignment around intron
           #DanaFarber_Liver_2_ppfoIcPb_12694765 (merged)
          my $a = $rel - $iend + $istart - 30;
          my $l = 29;
          if ($a < 0){
            $l = $l + $a;
            $a = 0;
          }
          #print LOG join("\n", "R: ".substr($ref, $rel-$iend+$istart-30, 29).lc(substr($ref, $rel-($iend-$istart+1), 2))."...".lc(substr($ref, $rel-2, 2)).substr($ref, $rel, 33),
                               #"   ".substr($matches, $rel-$iend+$istart-30, 31)."   ".substr($matches, $rel-2, 35),
                               #"Q: ".substr($query, $rel-$iend+$istart-30, 31)."...".substr($query, $rel-2, 35))."\n\n";          
          #print LOG "ISTART=$istart; IEND=$iend\n";
          print LOG join("\n", "R: ".substr($ref, $a, $l).lc(substr($ref, $a+$l, 2))."...".lc(substr($ref, $rel-2, 2)).substr($ref, $rel, 33),
                               "   ".substr($matches, $a, $l+2)."   ".substr($matches, $rel-2, 35),
                               "Q: ".substr($query, $a, $l+2)."...".substr($query, $rel-2, 35))."\n\n";
        }
      } #end cigar loop

      print LOG "  => ".($read_has_misaligned_splice_site ? "MISALIGNED READ" : "READ OK")."\n\n";

      #Check that the read covers all introns in the transcript model
      #Presence of each read intron in the transcript was checked above, so a total number check should suffice now
      my $intron_coord_match = 0;
      if (scalar @{$read_introns{'starts'}} == @{$intron_coords_by_tid{$tid}{'donors'}}){
        $intron_coord_match = 1;
      }
      else{
        print LOG "  => Different number of introns: READ = ".scalar(@{$read_introns{'starts'}})." - TRANSCRIPT = ".scalar(@{$intron_coords_by_tid{$tid}{'donors'}})."\n\n";
        push(@misaligned_reads, $read_name);
      }
      
      #Read is valid if it has the same splice sites as the transcript and none is misaligned
      if (!($read_has_misaligned_splice_site) and $intron_coord_match){
        $transcript_has_valid_support = 1;
        push(@valid_reads, $read_name);
        if ($keep_it_simple){
          last READ; #if 'simple' flag is set, stop after finding a single valid supporting read
        }
      }
    } #end read loop
  }
  else{
    print LOG "  ...skipping transcript: no introns\n";
  }
  print OUT join("\t", $tid,
                      ($transcript_has_valid_support ? "yes" : "no"),
                      ($intron_coords_by_tid{$tid}{'donors'} ? scalar(@{$intron_coords_by_tid{$tid}{'donors'}}) : "na"),
                      scalar(@valid_reads),
                      scalar(@misaligned_reads),
                      join(",", @valid_reads) || "na",
                      join(",", uniq @misaligned_reads) || "na",
                )."\n";
}
close (OUT);
close (LOG);



###################







