#!/usr/bin/perl -w


#1- Read BAM file
#2- Filter reads (mapping quality, exon size, ...)
#3- Write GFF3 file (this step also filters out short terminal exons)

#No splice site sequence filter implemented at the moment (see check_gff3_introns.pl)


#chr22   chr22.fa.gmap   cDNA_match      16124501        16124973        99      +       .       ID=HsBN_lib10:Read_468454-Barcode=BC256:length=1220.path1;Name=HsBN_lib10:Read_468454-Barcode=BC256:length=1220;Target=HsBN_lib10:Read_468454-Barcode=BC256:length=1220 1 473;Gap=M473
#chr22   chr22.fa.gmap   cDNA_match      16162397        16162487        98      +       .       ID=HsBN_lib10:Read_468454-Barcode=BC256:length=1220.path1;Name=HsBN_lib10:Read_468454-Barcode=BC256:length=1220;Target=HsBN_lib10:Read_468454-Barcode=BC256:length=1220 474 564;Gap=M91

use strict;
use Getopt::Long;
use Bio::DB::HTS;


my $quality_threshold = 10;
my $include_mono;
my $min_largest_exon_length = 30;
my $exclude_short_terminal_exons;
my $min_int_exon_length = 20;
my $min_ter_exon_length = 20;

my $bam_file;
my $fasta_file;
my $only_chr;
my $start;
my $end;
my $outfile;

&GetOptions(
	    'bam=s'      => \$bam_file,
	    'fasta=s'    => \$fasta_file,
            'chr=s'      => \$only_chr,
            'start=s'    => \$start,
            'end=s'      => \$end,
            'out=s'      => \$outfile,
            'qual=i'     => \$quality_threshold,
            'in_mono'    => \$include_mono,
            'min_lst_ex' => \$min_largest_exon_length,
            'clip_ends'  => \$exclude_short_terminal_exons,
            'min_int_ex' => \$min_int_exon_length,
            'min_ter_ex' => \$min_ter_exon_length,
           );


my $usage =  <<_USAGE_;

perl convert_bam_to_gff3_for_pasa.pl -bam file.bam -fasta seq.fa -chr chr -out file.gff3 [-in_mono] [-clip_ends]

 -bam         input BAM file (either path or URL) [required]
 -fasta       genome/chromosome sequence file in Fasta format (either path or URL; .gz accepted) [required]
 -chr         chromosome name [optional]
 -start       start coordinate (default: chromosome start) [optional]
 -end         end coordinate (default: chromosome end) [optional]
 -out         output GFF3 file [required]
 -qual        minimum quality score (default: 10) [optional]
 -min_lst_ex  minimum largest exon length (defaul: 30) [optional]
 -min_int_ex  minimum internal exon length (defaul: 20) [optional]
 -min_ter_ex  minimum terminal exon length (defaul: 20) [optional]
 -clip_ends   exclude terminal exons shorter than the minimum terminal exon length (default: no) [optional]
 -in_mono     include monoexonic features (default: no) [optional]

_USAGE_

die $usage unless ($bam_file and $fasta_file and $outfile);


open (LOG, ">$outfile.log") or die "can't open $outfile.log: $!";


#Read BAM file
#Input: BAM file or flat file with list of BAM files
#Output: array of Bio::DB::Sam alignment objects
my $alignments = read_bam_file($bam_file, $fasta_file);
print "N1=".scalar(@$alignments)." (Number of original reads) \n";


#Filter alignments
#Input: array of Bio::DB::Sam alignment objects
#Output: array of Bio::DB::Sam alignment objects
my $filtered_alignments = filter_alignments($alignments);
print "N2=".scalar(@$filtered_alignments)." (Number of reads after filtering) \n";

#Print alignment features
#Input: array of Bio::DB::Sam alignment objects, output file name
#Output: GFF3 file
print_gff3($filtered_alignments, $outfile);


close (LOG);



#####################

sub read_bam_file {
  my ($bam_file, $fasta_file) = @_;
  my $sam = Bio::DB::HTS->new (-bam           => $bam_file,
                               -fasta         => $fasta_file,
			       -autoindex     => 1,
			       -split_splices => 1,
                               -expand_flags  => 1,
                               -force_refseq  => 1,
                              );
 
  #my $header = $sam->header;
  #my @targets    = $sam->seq_ids;
  my @alignments = $sam->get_features_by_location(-seq_id => $only_chr,
                                                  -start  => $start,
                                                  -end    => $end
						 );
  #my $dna = $sam->seq($only_chr);
					
  return \@alignments;
}



sub filter_alignments {
  my $alignments = shift;
  my %filtered_alignments;
  my @filtered_alignments;

  while (my $a = shift @$alignments){
    my $read_name = $a->query->name;
    my @exons = $a->get_SeqFeatures; #NOTE: This method returns nothing unless there is an "N" in the cigar string (~ single-exon reads)
    if (scalar @exons == 0){
      push (@exons, $a);
    }

    #Filter out alignments on ERCC
    if ($a->seq_id =~ /^ERCC/){
      print LOG "UNWANTED_CHR\t$read_name\t".$a->cigar_str."\t".join("-", map {$_->length} @exons)."\n";
      next;
    }

    #Filter out low-quality alignments
    #NOTE: STAR alignments can only have 255, 3, 2, 1, 0 as mapping quality values, where 255 = uniquely mapped reads, 3 = read maps to 2 locations, and so on.
    if ($a->qual < $quality_threshold){
      print LOG "LOW_QUALITY\t$read_name\t".$a->cigar_str."\t".join("-", map {$_->length} @exons)."\n";
      next;
    } 

    #Exclude single-exon alignments
    if (scalar(@exons)==1 and !($include_mono)){
      print LOG "SINGLE_EXON\t$read_name\t".$a->cigar_str."\t".join("-", map {$_->length} @exons)."\n";
      next;
    }

    #Exon size
    #Discard alignments with all exons under 30bp in length
    if (scalar(@exons) == scalar(grep {$_->length < $min_largest_exon_length} @exons)){
      print LOG "MIN_LARGEST_EXON_LENGTH\t$read_name\t".$a->cigar_str."\t".join("-", map {$_->length} @exons)."\n";
      next;
    }

    my @int_exons = @exons[1..$#exons-1];
    #Discard alignments with any internal exon under 20bp in length
    if (scalar(grep {$_->length < $min_int_exon_length} @int_exons)){
      print LOG "MIN_INTERNAL_EXON_LENGTH\t$read_name\t".$a->cigar_str."\t".join("-", map {$_->length} @exons)."\n";
      next;
    }

    #Discard terminal exons under 20bp in length
    if ($exons[0]->length < $min_ter_exon_length or $exons[$#exons]->length < $min_ter_exon_length){ 
      #print LOG "MIN_TERMINAL_EXON_LENGTH: $read_name \t ".$a->cigar_str."\n";
      #I don't know how to remove features from alignment objects
    }


    #All but one alignments of a multimapper read should have a flag indicating secondary alignment (part of the bitwise flag). However, this doesn't happen for the SLR-seq alignments
      #NH:i: tag indicates the number of alignments for a given read in the BAM file
    #Alternative: filter by quality and then by aligned sequence length - Perhaps it should look at NM alignment tag instead? (choose smallest NM value - NOT ALWAYS!) 
    #unless ($filtered_alignments{$read_name}{ali} and ($filtered_alignments{$read_name}{ali}->qual > $a->qual or ($filtered_alignments{$read_name}{ali}->qual == $a->qual and length($filtered_alignments{$read_name}{ali}->dna) >= length($a->dna )))){

    if ($a->has_tag('SECONDARY') or #works for STAR
        #should work for GMAP
        ($filtered_alignments{$read_name}{ali} and ($filtered_alignments{$read_name}{ali}->qual > $a->qual or 
        ($filtered_alignments{$read_name}{ali}->qual == $a->qual and $filtered_alignments{$read_name}{ali}->get_tag_values('NM') < $a->get_tag_values('NM'))))
    ){
      print LOG "MULTIMAPPER\t$read_name\t".$a->cigar_str."\t".join("-", map {$_->length} @exons)."\n";
      next;
    }
    else{
      $filtered_alignments{$read_name}{ali} = $a;
    }
  }

  #Store filtered alignments
  foreach my $read_name (keys %filtered_alignments){
    push (@filtered_alignments, $filtered_alignments{$read_name}{ali});
  }
  
  return \@filtered_alignments;
}



sub print_gff3 {
  my ($filtered_alignments, $outfile) = @_;
  open (GFF3, ">$outfile") or die "can't open $outfile: $!";
  print GFF3 "##gff-version   3\n";


  while (my $fa = shift @$filtered_alignments){
    my $score = $fa->qual;
    my $chr = $fa->seq_id;
    my $strand; 
    if ($fa->strand==1){
      $strand = "+";
    }
    elsif ($fa->strand==-1){
      $strand = "-";
    }
    my $read_name = $fa->query->name;
    my $spliced_length = $fa->l_qseq; #$fa->query->length;


    #Deal with soft clipping (Bio::DB::Sam gets coordinates wrong)
    my $left_padding = 0;
    if ($fa->cigar_str =~ /^(\d+)(S|H)/){
      $left_padding = $1;
    }
    my $right_padding = 0;
    if ($fa->cigar_str =~ /(\d+)(S|H)$/){
      $right_padding = $1;
    }


    #Get overall percent id (it should be done at the exon level but this is quicker)
    #Number of matches
    my ($ref, $matches, $query) = $fa->padded_alignment;
    #print LOG "\n#$read_name\n$ref\n$matches\n$query\n\n";
    my @match_count = ($matches =~ /\|/g);
    #print GFF3 "\n\n$ref\n$matches\n$query\n".scalar(@match_count)."\n\n" if $read_name eq "HsBN_lib10:Read_432980-Barcode=BC008:length=1210";
    #my $percent_id = sprintf("%.0f", scalar(@match_count)/$spliced_length*100);
    #$percent_id = 100 if $percent_id > 100; #hack: $fa->query->length doesn't always give the right length (why?)

    #print $fa->qual."\t".$fa->l_qseq."\t".$fa->cigar2qlen."\n" if $read_name eq "HsBN_lib10:Read_432980-Barcode=BC008:length=1210";


####PACBIO####
    #Read length - remove clipped seq? Less conservative approach
    $spliced_length = $spliced_length - $left_padding - $right_padding; #Less stringent. Could also count deletions.
    if ($spliced_length == 0){
      print LOG "\nNULL_SPLICED_LENGTH\nname = $read_name\nspliced_length = ".$spliced_length."\nleft_padding = $left_padding\nright_padding = $right_padding\n\n";
      next;
    }
    #NOTE: number of matches could be calculated using the NM tag value (edit distance, which include ambiguous bases, so it is greater than the number of mismatches in the padded alignment)
    my $n_mismatches = $fa->get_tag_values('NM');
    my $percent_id = sprintf("%.0f", ($spliced_length - $n_mismatches)/$spliced_length*100);
####PACBIO####

####SLR-SEQ###
    #This works for the SLRseq GMAP alignments, but apparently underestimates id% for PacBio STAR alignments (l_qseq will include S and H bits)
#    my $percent_id = sprintf("%.0f", scalar(@match_count)/$spliced_length*100);
####SLR-SEQ###


    #if ($percent_id > 100){
      print LOG "\nname = $read_name\npercent_id = $percent_id\nquery_length = ".$fa->query->length."\nmatch_count = ".scalar(@match_count)."\nspliced_length = ".$spliced_length."\nleft_padding = $left_padding\nright_padding = $right_padding\n\n";
     #print LOG "\n#$read_name\n$ref\n$matches\n$query\n\n";
    #  $percent_id = 100; #hack
    #}

    my @features = $fa->get_SeqFeatures; #NOTE: This method doesn't return anything if there is no "N" in the cigar string (~ single-exon reads)
    if (scalar @features and $strand eq "-"){
      @features = sort {$b->start <=> $a->start} $fa->get_SeqFeatures;
    }
    if (scalar @features == 0){ #Single-exon reads
      push (@features, $fa);
    }

    my @lines_to_print;
    my ($gl_start, $gl_end);
    my $c=0;
    foreach my $feat (@features){
      my $feat_start = $feat->start;
      my $feat_end = $feat->end;
      if (scalar @features > 1){ #if there is one feature it will be the Alignment object, and this gives the correct start and stop despite hard/soft clipping
        $feat_start -= $left_padding unless (($c==0 and $strand eq "+") or ($c==$#features and $strand eq "-"));
        $feat_end -= $left_padding;
        $feat_end -= $right_padding if (($c==$#features and $strand eq "+") or ($c==0 and $strand eq "-"));
      } 
      my $feat_length = $feat_end - $feat_start + 1;

      #Exclude short terminal exons
      if ($exclude_short_terminal_exons and $feat_length < $min_ter_exon_length and ($c==0 or $c==$#features)){
        print LOG "MIN_TERMINAL_EXON_LENGTH\t$read_name\t".$fa->cigar_str."\n";
        next;
      }

      if ($gl_start and $gl_end){
        $gl_start = $gl_end + 1;
        $gl_end += $feat_length;
      }
      else{
        $gl_start = 1;
        $gl_end = $feat_length;
      }

      push (@lines_to_print, "$chr\tgff3_for_pasa\tcDNA_match\t".$feat_start."\t".$feat_end."\t$percent_id\t$strand\t.\tID=$read_name;Name=$read_name;Target=$read_name $gl_start $gl_end;Gap=".$fa->cigar_str);
#Gap value is not as in GMAP output, but apparently it is not used by PASA
      $c++;
    }
    #Print GFF3 lines - check number of features to avoid printing a single line for multi-exon reads where terminal exons have been removed leaving a single exon
    unless (scalar(@lines_to_print)==1 and !($include_mono)){
      print GFF3 join("\n", @lines_to_print)."\n";
    }
  }
  close (GFF3);
}


#chr22   chr22.fa.gmap   cDNA_match      16124501        16124973        99      +       .       ID=HsBN_lib10:Read_468454-Barcode=BC256:length=1220.path1;Name=HsBN_lib10:Read_468454-Barcode=BC256:length=1220;Target=HsBN_lib10:Read_468454-Barcode=BC256:length=1220 1 473;Gap=M473
#chr22   chr22.fa.gmap   cDNA_match      16162397        16162487        98      +       .       ID=HsBN_lib10:Read_468454-Barcode=BC256:length=1220.path1;Name=HsBN_lib10:Read_468454-Barcode=BC256:length=1220;Target=HsBN_lib10:Read_468454-Barcode=BC256:length=1220 474 564;Gap=M91






