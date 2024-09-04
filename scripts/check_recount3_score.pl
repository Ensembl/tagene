#!/usr/bin/env perl

#Check if all the introns of the transcripts in the list have a Recount3 score equal or greater than the given threshold.
use strict;
use warnings;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::DB::HTS::Tabix;
use List::Util qw(uniq);

$|=1;

my $gtf_file;
my $cls_id_list;
my $recount3_file;
my $cutoff = 50;
my $only_chr;
my $ref_dbhost;
my $ref_dbport;
my $ref_dbuser;
my $ref_dbname;
my $outfile;

&GetOptions(
            'gtf=s'      => \$gtf_file,
            'idlist=s'   => \$cls_id_list,
            'chr=s'      => \$only_chr,
            'recount3=s' => \$recount3_file,
            'cutoff=i'   => \$cutoff,
            'refhost=s'  => \$ref_dbhost,
            'refport=i'  => \$ref_dbport,
            'refuser=s'  => \$ref_dbuser,
            'refdbname=s' => \$ref_dbname,
            'out=s'      => \$outfile,
           );


#Connect to database with reference annotation
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
      next unless $slice->seq_region_name eq $only_chr;
    }
    foreach my $gene (@{$slice->get_all_Genes}){
      unless ($gene->biotype eq "artifact" or !($gene->stable_id =~ /^ENS/) or scalar(grep {$_->value eq "not for VEGA"} @{$gene->get_all_Attributes})){
        foreach my $transcript (@{$gene->get_all_Transcripts}){
          unless ($transcript->biotype eq "artifact" or !($transcript->stable_id =~ /^ENS/) or scalar(grep {$_->value eq "not for VEGA"} @{$transcript->get_all_Attributes})){
            foreach my $intron (@{$transcript->get_all_Introns}){
              $ref_introns{$intron->seq_region_name.":".$intron->seq_region_start."-".$intron->seq_region_end} = 1;
            }
          }
        }
      }
    }
  }
}



#Open Tabix for Recount3
my $recount3_tabix = Bio::DB::HTS::Tabix->new( filename => $recount3_file );


#Read list of CLS models
my @cls_ids;
my %id_list;
if ($cls_id_list){
  open (IN, $cls_id_list) or die "Can't open $cls_id_list: $!";
  while (<IN>){
    chomp;
    push(@cls_ids, $_);
    $id_list{$_} = 1;
  }
  close (IN);
}


#Read GTF file with CLS models
my %exon_coords_by_tid;
open (GTF, "zcat $gtf_file |") or die "Can't open $gtf_file: $!";
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
    if ($cols[8] =~ /transcript_id \"(\w+)\";/){
      my $tid = $1;
      #Filter by ID list if any
      if ($cls_id_list){
        next unless $id_list{$tid};
      }
      else{
        push(@cls_ids, $tid);
      }
      #Store exon coordinates
      my ($chr, $start, $end, $strand) = @cols[0,3,4,6];
      $exon_coords_by_tid{$tid}{'chr'} = "chr".$chr;
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

open (OUT, ">$outfile") or die "Can't open $outfile: $!";
print OUT "#".join("\t", "transcript_ID",
                        "has_valid_support",
                        "number_of_introns",
                        "intron_scores",
                  )."\n";
                  
#Loop through CLS transcripts
foreach my $tid (uniq @cls_ids){
  my $transcript_has_valid_support = 0;
  my $passed = 1;
  my @scores;
  if ($intron_coords_by_tid{$tid}{'starts'}){
    for (my $i=0; $i < scalar(@{$intron_coords_by_tid{$tid}{'starts'}}); $i++){
      my $score = 0;
      my $coords = $intron_coords_by_tid{$tid}{'chr'}.":".$intron_coords_by_tid{$tid}{'starts'}[$i]."-".$intron_coords_by_tid{$tid}{'ends'}[$i];
      my $iter = $recount3_tabix->query($coords);
      if ($iter){
        while (my $line = $iter->next){
          my @cols = split(/\t/, $line);
          if ($cols[1]+1 == $intron_coords_by_tid{$tid}{'starts'}[$i] and $cols[2] == $intron_coords_by_tid{$tid}{'ends'}[$i] and $cols[5] eq $intron_coords_by_tid{$tid}{'strand'}){
            $score = $cols[4];
            last;
          }
        }
      }
      push(@scores, $score);
      if ($score < $cutoff){
        if ($ref_introns{$coords}){ #ignore low score if intron was in the reference annotation
          $passed = 1;
          $scores[-1] .= "(a)";
        }
        else{ 
          $passed = 0;
        }
      }
    }
    print OUT join("\t", $tid, ($passed ? "yes" : "no"), scalar(@{$intron_coords_by_tid{$tid}{'starts'}}), join(",", @scores))."\n";
  }
  else{ #mouse liftover model with one exon
    print OUT join("\t", $tid, "na", 0, 0)."\n";
  }
}
$recount3_tabix->close;
close (OUT);

