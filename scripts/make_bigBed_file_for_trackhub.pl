#!/usr/bin/env perl

#This script dumps annotation from a database in GenePred format (extended BED 12+10).
#Then, it calls UCSC tools to convert it into a bigGenePred file.
#Transcripts will be coloured differently if there is a TAGENE_transcript remark.
#Optionally, a list of gene stable ids can be provided.
#If a path to an FTP directory is provided, it will upload the output file.


use strict;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use POSIX qw(strftime);
use Cwd;

my $host;
my $port;
my $user = 'ensro';
my $dbname;
my $gene_list;
my $outfile;
my $run_date; #YYYYmmdd
my $colour0 = "0,0,0";
my $colour1; 
my $colour2;
my $path_to_ftp;
my $ftp_user;

&GetOptions(
             'host=s'      => \$host,
             'port=i'      => \$port,
             'user=s'      => \$user,
             'dbname=s'    => \$dbname,
             'list=s'      => \$gene_list,
             'date=s'      => \$run_date,
             'colour1=s'   => \$colour1,
             'colour2=s'   => \$colour2,
             'out=s'       => \$outfile,
             'ftp_path=s'  => \$path_to_ftp,
             'ftp_user=s'  => \$ftp_user,
           );
$outfile ||= $dbname;

#Connect to database and get adaptors
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host => $host,
    -port => $port,    
    -user => $user,
    -dbname => $dbname
);
my $meta_container = $db->get_MetaContainer();
my $sa = $db->get_SliceAdaptor();
my $ta = $db->get_TranscriptAdaptor();

my @sources = ('ensembl', 'havana', 'ensembl_havana'); # fetch only these gene sources from the Havana database

#Read list of gene stable ids
my %selected_genes;
if ($gene_list){
  open (LIST, $gene_list) or die "can't open $gene_list";
  while (<LIST>){
    chomp;
    $selected_genes{$_} = 1;
  }
  close (LIST);
}


open BED, ">$outfile.bed" or die "could not open file $outfile.bed $!";

#Fetch havana db genes (include PAR regions, exclude non-reference regions)
#Sort chromosomes
foreach my $slice (sort {
                          if ($a->seq_region_name =~ /\d/ and 
                              $b->seq_region_name =~ /\d/){
                                $a->seq_region_name <=> $b->seq_region_name;
                          } 
                          else{
                            $a->seq_region_name cmp $b->seq_region_name;
                          }
                        }
              @{$sa->fetch_all('chromosome', undef, 0, 1)}){

  print $slice->seq_region_name."...";
  foreach my $gene (@{$slice->get_all_Genes}){
    if ($gene_list){
      next unless $selected_genes{$gene->stable_id};
    }
    next if $gene->biotype eq "artifact";
    #next if grep {$_->value eq "not for VEGA"} @{$gene->get_all_Attributes('remark')};
    next unless grep {$_ eq $gene->source} @sources;
    my $mod_num = 1; # number for appending to gene_name (for display in ensembl)
    
    foreach my $transcript (@{$gene->get_all_Transcripts}){
	    next if $transcript->biotype eq "artifact";
	    #next if grep {$_->value eq "not for VEGA"} @{$transcript->get_all_Attributes('remark')};
	    next unless grep {$_ eq $transcript->source} @sources;
		
	    #Assign colour based on TAGENE remark and created/modified date
	    my $colour;
      if (is_tagene($transcript)){
        if ($run_date){
          if (is_created_on_date($transcript, $run_date)){
            $colour = $colour1;
          }
          elsif (is_modified_on_date($transcript, $run_date)){
            $colour = $colour2;
          }
        }
        else{
          $colour = $colour1;
        }
      }
	    #Print BED line
	  	my $transcript_line = get_transcript_line($transcript->stable_id, $transcript->version, $mod_num, $colour);
      if ($transcript_line){
        print BED $transcript_line."\n";
        $mod_num++;
      }
	  }
  }	
}
close (BED);

print "\n";

#Sort BED file
system("sort -k1,1 -k2,2n $outfile.bed > $outfile.sorted.bed");

#Get chromosome size file
my $ucsc_name = $meta_container->single_value_by_key('assembly.ucsc_alias');
system("wget https://hgdownload.soe.ucsc.edu/goldenPath/$ucsc_name/bigZips/$ucsc_name.chrom.sizes");

#Create bigGenePred definition file
my $bgp_def_file = "genes.as";
print_bgp_definition_file($bgp_def_file);

#Make bigBed (bigGenePred) file
system("bedToBigBed -type=bed12+10 -as=$bgp_def_file -extraIndex=name $outfile.sorted.bed $ucsc_name.chrom.sizes $outfile.bb");

#Compress BED files
system("gzip $outfile.bed $outfile.sorted.bed");

#Upload file to FTP
my $cwd = getcwd;
if ($path_to_ftp){
  my $become_user = "";
  if ($ftp_user){
    $become_user = "become $ftp_user";
  }
  system("bsub -q datamover -I $become_user cp -p $cwd/$outfile.bb $path_to_ftp");
}



##########



sub exon_chain {
  my $transcript = shift;
  return join(":", map {$_->seq_region_start."-".$_->seq_region_end} @{$transcript->get_all_Exons});
}


sub cds_chain {
  my $transcript = shift;
  if ($transcript->translate){
    return join(":", map {$_->seq_region_start."-".$_->seq_region_end} @{$transcript->get_all_translateable_Exons});
  }
  else{
    return -1;
  }
}


sub is_tagene {
  my $transcript = shift;
  if (scalar grep {$_->value eq "TAGENE_transcript"} @{$transcript->get_all_Attributes('remark')}){
    return 1;
  }
  return 0;
}


sub is_created_on_date {
  my ($transcript, $date) = @_;
  my $created_date = strftime("%Y%m%d", localtime($transcript->created_date));
  if ($created_date eq $date){
    return 1;
  }
  return 0;
}


sub is_modified_on_date {
  my ($transcript, $date) = @_;
  my $modified_date = strftime("%Y%m%d", localtime($transcript->modified_date));
  if ($modified_date eq $date){
    return 1;
  }
  return 0;
}


sub get_transcript_line {
  my ($ltrans_id, $lversion, $mod_num, $colour) = @_;
  $mod_num += 1000;
  $mod_num =~ s/1/-/;
  $colour ||= $colour0;
  my $ltrans = $ta->fetch_by_stable_id_version($ltrans_id, $lversion);
  my $seq_region = $ltrans->seq_region_name;
  $seq_region = "chr".$seq_region;
  my ($cds_start_stat, $cds_end_stat) = ($ltrans->coding_region_start, $ltrans->coding_region_end);
  $cds_start_stat = $cds_start_stat ? 'cmpl' : 'none';
  $cds_end_stat   = $cds_end_stat   ? 'cmpl' : 'none';
  my $cds = $ltrans->get_all_CDS;
  # set default CDS start and end equal to feature start and end (they will be overwritten if there is a CDS)
  my ($cds_start, $cds_end) = ($ltrans->seq_region_start, $ltrans->seq_region_end);
  my $seq_strand = $ltrans->seq_region_strand;
  if (@{ $cds }) {
    ($cds_start, $cds_end) = $seq_strand eq '-1'
                         ? 
    ($cds->[-1]->seq_region_start, $cds->[0]->seq_region_end) 
                         :
    ($cds->[0]->seq_region_start, $cds->[-1]->seq_region_end);
  }
  my $trans_biotype = $ltrans->biotype;
  my $exons = $ltrans->get_all_ExonTranscripts; 
  my $num_exons = @{ $exons };
  my (@exon_lengths, @exon_starts, @exon_frames);
  foreach my $exon (@{ $exons }) {
    push @exon_lengths, $exon->length;
    push @exon_starts, ( $exon->start - $ltrans->seq_region_start);
    my $e_frame = $exon->phase;
    $e_frame = $e_frame eq q{.} ? "-1" : $e_frame;
    push @exon_frames, $e_frame;
  }
  my $exon_lengths = $seq_strand eq '-1' ? join",", reverse @exon_lengths : join',', @exon_lengths;
  my $exon_starts  = $seq_strand eq '-1' ? join",", reverse @exon_starts  : join',', @exon_starts;
  my $exon_frames  = $seq_strand eq '-1' ? join",", reverse @exon_frames  : join',', @exon_frames;
  unless ($exon_starts =~ /^0/){ #catch error that would prevent display in UCSC browser
    return;
  }
  my $tgene = $ltrans->get_Gene;
  my $gene_stid = $tgene->stable_id;
  my $gene_attribs = $tgene->get_all_Attributes;
  my ($gene_name, $mod_gene_name);
  my $trans_name = $ltrans->version ? $ltrans->stable_id.'.'.$ltrans->version : $ltrans->stable_id;
  foreach my $attrib (@{$gene_attribs}) {
    if ($attrib->name eq 'Name') {
      $gene_name = $attrib->value;
      $mod_gene_name = $gene_name . '-' . $trans_name;
    }
  }
  $seq_strand = $seq_strand == 1 ? '+' : '-';
  my $trans_attribs = get_transcript_attributes($ltrans);

  my $modif_date = strftime("%d-%b-%Y"."_%H:%M:%S", localtime($ltrans->modified_date));

  my $score = $tgene->biotype eq 'protein_coding' ? 100 : 1000;
  return join("\t", $seq_region, 
                    ($ltrans->seq_region_start - 1), 
                    $ltrans->seq_region_end, 
                    ($mod_gene_name || $trans_name),
                    $score, 
                    $seq_strand, 
                    ($cds_start -1), 
                    $cds_end, 
                    $colour, 
                    $num_exons, 
                    $exon_lengths,
                    $exon_starts, 
                    $trans_name,
                    $cds_start_stat, 
                    $cds_end_stat, 
                    $exon_frames, 
                    $ltrans->biotype, 
                    $gene_stid, 
                    ($gene_name || "no_name"), 
                    $tgene->biotype,
                    ($trans_attribs || "-"),
                    $modif_date
                    );
}


#Adapted from gencode/modules/Gencode/Ensembl2GTF_new.pm 
sub get_transcript_attributes {
  my $transcript = shift;
  my @attributes = ();

  my @att_values = ("alternative 3' UTR",
                    "alternative 5' UTR",
                    "bicistronic",
                    "CAGE supported TSS",
                    "dotter confirmed",
                    "downstream ATG",
                    "inferred exon combination",
                    "inferred transcript model",
                    "low sequence quality",
                    "MANE_plus",
                    "MANE_select",
                    "NAGNAG splice site",
                    "NMD exception",
                    "NMD likely if extended",
                    "[Nn]o.+org.+sup.+ted",
                    "non-ATG start",
                    "non canonical conserved",
                    "non canonical genome sequence error",
                    "non canonical other",
                    "non canonical polymorphism",
                    "non canonical TEC",
                    "non canonical U12",
                    "non-submitted evidence",
                    "not best-in-genome evidence",
                    "not for VEGA",
                    "not organism-supported",
                    "overlapping uORF",
                    "readthrough",
                    "retained intron CDS",
                    "retained intron final",
                    "retained intron first",
                    "RNA-Seq supported only",
                    "RNA-Seq supported partial",
                    "RP supported TIS",
                    "sequence error",
                    "TAGENE_transcript",
                    "upstream ATG",
                    "upstream uORF",
                    "3' nested supported extension",
                    "3' standard supported extension",
                    "5' nested supported extension",
                    "5' standard supported extension",
                    "454 RNA-Seq supported",
                    "nested 454 RNA-Seq supported"
                   );

  my $att_values_s = join("|", @att_values);

  foreach my $att (@{$transcript->get_all_Attributes}){
    my $val  = $att->value;
    my $code = $att->code;
    if ($code =~ /cds_end_NF|cds_start_NF|mRNA_end_NF|mRNA_start_NF/ and $val != 0){
      push(@attributes, $code);
    }
    elsif(($code eq "remark" or $code eq "hidden_remark" or $code eq "havana_cv") and ($val =~/^($att_values_s)$/)){
      $val =~ s/ /_/g;
      $val =~ s/'//g;
      $val =~ s/-/_/g;
      $val =~ s/^[Nn]o.+org.+sup.+ted$/not_organism_supported/;
      push(@attributes, $val);
    }
  }
  return join(",", @attributes);
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
   string tags;         "Transcript attributes"
   string modifDate;    "Modification date"
   )
TABLE
  print T $text."\n";
};

