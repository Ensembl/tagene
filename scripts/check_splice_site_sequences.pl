#!/usr/bin/env perl

#Check if all the introns of the transcripts in the list have allowed splice site sequences
#Provide allowed sequences in a comm-separated list:
#eg. -allowed GT-AG,GC-AG
#If reference annotation is provided, annotated splice sites are not checked

use strict;
use warnings;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::DB::HTS::Faidx;
use List::Util qw(uniq);

$|=1;

my $gtf_file;
my $fasta_file;
my $id_list;
my $only_chr;
my $ref_dbhost;
my $ref_dbport;
my $ref_dbuser;
my $ref_dbname;
my $allowed_seqs;
my $outfile;

&GetOptions(
            'gtf=s'       => \$gtf_file,
            'fasta=s'     => \$fasta_file,
            'idlist=s'    => \$id_list,
            'chr=s'       => \$only_chr,
            'refhost=s'   => \$ref_dbhost,
            'refport=i'   => \$ref_dbport,
            'refuser=s'   => \$ref_dbuser,
            'refdbname=s' => \$ref_dbname,
            'allowed=s'   => \$allowed_seqs,
            'out=s'       => \$outfile,
           );

my @allowed_splice_site_seqs;
if ($allowed_seqs){
  @allowed_splice_site_seqs = split(/,/, $allowed_seqs);
}

#Get intron coordinates from the reference annotation
my %ref_introns;
if ($ref_dbhost and $ref_dbport and $ref_dbuser and $ref_dbname){
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
            -host => $ref_dbhost,
            -port => $ref_dbport,    
            -user => $ref_dbuser,
            -dbname => $ref_dbname,
  );
  my $sa = $db->get_SliceAdaptor();
  foreach my $slice (@{$sa->fetch_all("toplevel")}){
    next unless $slice->coord_system->is_default;
    if ($only_chr){
      $only_chr =~ s/^chr//;
      next unless $slice->seq_region_name eq $only_chr;
    }
    foreach my $gene (@{$slice->get_all_Genes}){
      unless ($gene->biotype eq "artifact" or !($gene->stable_id =~ /^ENS/) or scalar(grep {$_->value eq "not for VEGA"} @{$gene->get_all_Attributes})){
        foreach my $transcript (@{$gene->get_all_Transcripts}){
          unless ($transcript->biotype eq "artifact" or !($transcript->stable_id =~ /^ENS/) or scalar(grep {$_->value eq "not for VEGA"} @{$transcript->get_all_Attributes})){
            foreach my $intron (@{$transcript->get_all_Introns}){
              my $coord = $intron->seq_region_name.":".$intron->seq_region_start."-".$intron->seq_region_end;
              $coord .= ":".($intron->seq_region_strand == 1 ? "+" : "-");
              $ref_introns{$coord} = 1;
            }
          }
        }
      }
    }
  }
}


#Read list of selected RNA-seq transcript ids, if any
my @selected_ids;
my %selected_id_list;
if ($id_list){
  open (IN, $id_list) or die "Can't open $id_list: $!";
  while (<IN>){
    chomp;
    push(@selected_ids, $_);
    $selected_id_list{$_} = 1;
  }
  close (IN);
}


#Read input GTF file
my %exon_coords_by_tid;
if ($gtf_file =~ /\.gz$/){
  open (GTF, "zcat $gtf_file |") or die "Can't open $gtf_file : $!";
}
else{
  open (GTF, "cat $gtf_file |") or die "Can't open $gtf_file : $!";
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
        next unless $selected_id_list{$tid};
      }
      else{
        push(@selected_ids, $tid);
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
    my @starts = sort {$a<=>$b} @{$exon_coords_by_tid{$tid}{'starts'}};
    my @ends = sort {$a<=>$b} @{$exon_coords_by_tid{$tid}{'ends'}};
    @{$intron_coords_by_tid{$tid}{'starts'}} = map {$_ + 1} @ends[0..$#ends-1];
    @{$intron_coords_by_tid{$tid}{'ends'}} = map {$_ - 1} @starts[1..$#starts];
  }
}


#Open index for the genome sequence file
my $index = Bio::DB::HTS::Faidx->new($fasta_file);


open (OUT, ">$outfile") or die "Can't open $outfile: $!";
print OUT "#".join("\t", "transcript_ID",
                        "has_valid_support",
                        "number_of_introns",
                        "splice_site_sequences",
                  )."\n";
 
#Loop through transcripts
foreach my $tid (uniq @selected_ids){
  my $transcript_has_valid_support = 0;
  my $passed = 1;
  my @splice_site_seqs;
  if ($intron_coords_by_tid{$tid}{'starts'}){
    for (my $i=0; $i < scalar(@{$intron_coords_by_tid{$tid}{'starts'}}); $i++){
      my $seq;
      my $coords = $intron_coords_by_tid{$tid}{'chr'}.":".$intron_coords_by_tid{$tid}{'starts'}[$i]."-".$intron_coords_by_tid{$tid}{'ends'}[$i].":".$intron_coords_by_tid{$tid}{'strand'};
      my $donor_coord = $intron_coords_by_tid{$tid}{'chr'}.":".$intron_coords_by_tid{$tid}{'starts'}[$i]."-".($intron_coords_by_tid{$tid}{'starts'}[$i]+1);
      my $acceptor_coord = $intron_coords_by_tid{$tid}{'chr'}.":".($intron_coords_by_tid{$tid}{'ends'}[$i]-1)."-".$intron_coords_by_tid{$tid}{'ends'}[$i];
      my $donor_seq = $index->get_sequence_no_length($donor_coord);
      my $acceptor_seq = $index->get_sequence_no_length($acceptor_coord);
      $seq = $donor_seq."-".$acceptor_seq;
      if ($intron_coords_by_tid{$tid}{'strand'} eq "-"){
        $seq =~ tr/ACGT/TGCA/;
        $seq = reverse $seq;
      }
      push(@splice_site_seqs, $seq);

      if (scalar @allowed_splice_site_seqs){
        if (grep {$_ eq $seq} @allowed_splice_site_seqs){
          $passed = 1;
        }
        else{
          if ($ref_introns{$coords}){ #ignore sequence if intron was in the reference annotation
            $passed = 1;
            $splice_site_seqs[-1] .= "(a)";
          }
          else{ 
            $passed = 0;
          }
        }
      }
      else{
        $passed = 1;
      }
    }
    print OUT join("\t", $tid, ($passed ? "yes" : "no"), scalar(@{$intron_coords_by_tid{$tid}{'starts'}}), join(",", @splice_site_seqs), $intron_coords_by_tid{$tid}{'strand'})."\n";
  }
}
close (OUT);





