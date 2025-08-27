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
my $id_list;
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
            'idlist=s'   => \$id_list,
            'chr=s'      => \$only_chr,
            'recount3=s' => \$recount3_file,
            'cutoff=i'   => \$cutoff,
            'refhost=s'  => \$ref_dbhost,
            'refport=i'  => \$ref_dbport,
            'refuser=s'  => \$ref_dbuser,
            'refdbname=s' => \$ref_dbname,
            'out=s'      => \$outfile,
           );

if ($only_chr){
  $only_chr =~ s/^chr//;
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



#Open Tabix for Recount3
my $recount3_tabix = Bio::DB::HTS::Tabix->new( filename => $recount3_file );


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
    my $chr = $cols[0];
    $chr =~ s/^chr//;
    if ($only_chr){
      next unless $chr eq $only_chr;
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
      my ($start, $end, $strand) = @cols[3,4,6];
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

open (OUT, ">$outfile") or die "Can't open $outfile: $!";
print OUT "#".join("\t", "transcript_ID",
                        "has_valid_support",
                        "number_of_introns",
                        "intron_scores",
                  )."\n";
                  
#Loop through transcripts
foreach my $tid (uniq @selected_ids){
  my $transcript_has_valid_support = 0;
  my $passed = 0;
  my @scores;
  my $failed_splice_sites = 0;
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
        $coords = $coords.":".$intron_coords_by_tid{$tid}{'strand'};
        if ($ref_introns{$coords}){ #ignore low score if intron was in the reference annotation
          $scores[-1] .= "(a)";
        }
        else{ 
          $failed_splice_sites++;
        }
      }
    }
    if ($failed_splice_sites == 0){
      $passed = 1;
    }
    print OUT join("\t", $tid, ($passed ? "yes" : "no"), scalar(@{$intron_coords_by_tid{$tid}{'starts'}}), join(",", @scores))."\n";
  }
  else{ #mouse liftover model with one exon
    print OUT join("\t", $tid, "na", 0, 0)."\n";
  }
}
$recount3_tabix->close;
close (OUT);

