#!/usr/bin/env perl

#Check for novel introns that overlap MANE_Select exons
#Introns will be flagged if both splice sites overlap the same exon

use strict;
use warnings;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use List::Util qw(uniq);

$|=1;

my $gtf_file;
my $id_list;
my $only_chr;
my $host;
my $port;
my $user;
my $dbname;
my $outfile;

&GetOptions(
            'gtf=s'      => \$gtf_file,
            'idlist=s'   => \$id_list,
            'chr=s'      => \$only_chr,
            'host=s'     => \$host,
            'port=i'     => \$port,
            'user=s'     => \$user,
            'dbname=s'   => \$dbname,
            'out=s'      => \$outfile,
           );

if ($only_chr){
  $only_chr =~ s/^chr//;
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

#Connect to reference annotation database and get MANE_Select exons
my %ref_exons;
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
            -host => $host,
            -port => $port,    
            -user => $user,
            -dbname => $dbname,
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
          my $is_reference = 0;
          if ($dbname =~ /homo_sapiens/ or $dbname =~ /human/){
            if (scalar(grep {$_->value eq "MANE_select"} @{$transcript->get_all_Attributes("remark")}) or
                scalar(@{$transcript->get_all_Attributes("MANE_Select")})
            ){
              $is_reference = 1;
            }
          }
          else{
            if ($transcript->is_canonical){
              $is_reference = 1;
            }
          }
          if ($is_reference){
            foreach my $exon (@{$transcript->get_all_Exons}){
              my $strand = $exon->seq_region_strand == 1 ? "+" : "-";
              $ref_exons{$exon->seq_region_name}{$strand}{$exon->seq_region_start."-".$exon->seq_region_end} = 1;
            }
          }
        }
      }
    }
  }
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
                        "intron_coord",
                        "MS_exon_coord",
                  )."\n";
                  
#Loop through transcripts
foreach my $tid (uniq @selected_ids){
  if ($intron_coords_by_tid{$tid}{'starts'}){
    INTRON:for (my $i=0; $i < scalar(@{$intron_coords_by_tid{$tid}{'starts'}}); $i++){
      my $score = 0;
      my $coords = $intron_coords_by_tid{$tid}{'chr'}.":".$intron_coords_by_tid{$tid}{'starts'}[$i]."-".$intron_coords_by_tid{$tid}{'ends'}[$i];
      my $chr = $intron_coords_by_tid{$tid}{'chr'};
      my $strand = $intron_coords_by_tid{$tid}{'strand'};
      my $intron_start = $intron_coords_by_tid{$tid}{'starts'}[$i];
      my $intron_end = $intron_coords_by_tid{$tid}{'ends'}[$i];
      #Compare novel intron coordinates with MANE Select exon coordinates
      foreach my $exon_coord (keys %{$ref_exons{$chr}{$strand}}){
        my ($exon_start, $exon_end) = split(/-/, $exon_coord);
        if ($intron_start > $exon_start and $intron_end < $exon_end){
          print OUT join("\t", $tid,
                              $chr.":".$intron_start."-".$intron_end.":".$strand,
                              $chr.":".$exon_start."-".$exon_end.":".$strand,
                        )."\n";
          last INTRON;
        }
      }
    }
  }
}
close (OUT);

