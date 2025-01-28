#!/usr/bin/env perl

#For each spliced transcript in the input GTF, find if it has a novel intron where both splice sites overlap a repeat element,
#Novel introns are those not annotated in the reference annotation database provided (optional).

use strict;
use warnings;
use Getopt::Long;
use Bio::DB::HTS::Tabix;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

$|=1;


my $repeat_file;
my $repeat_types;
my $input_gtf;
my $only_chr;
my $ref_dbhost;
my $ref_dbport;
my $ref_dbuser;
my $ref_dbname;
my $same_repeat_overlap;
my $outfile;

&GetOptions(
            'repeats=s'   => \$repeat_file,  #BED file with repeat class/type names in columns 7 and 8
            'types=s'     => \$repeat_types, #eg. -types "Tandem repeats,Type I Transposons/SINE"
            'gtf=s'       => \$input_gtf,
            'refhost=s'   => \$ref_dbhost,
            'refport=i'   => \$ref_dbport,
            'refuser=s'   => \$ref_dbuser,
            'refdbname=s' => \$ref_dbname,
            'same_repeat!' => \$same_repeat_overlap,
            'chr=s'       => \$only_chr,
            'out=s'       => \$outfile,
           );


#Open Tabix for the repeat file
my $repeat_tabix = Bio::DB::HTS::Tabix->new( filename => $repeat_file );


#Read GTF file with CLS models
print "Reading input GTF...\n";
my %exon_coords_by_tid;
my @tids;
if ($input_gtf =~ /\.gz$/){
  open (GTF, "zcat $input_gtf |") or die "Can't open $input_gtf: $!";
}
else{
  open (GTF, "cat $input_gtf |") or die "Can't open $input_gtf: $!";
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
      #Add to ID list
      unless ($exon_coords_by_tid{$tid}){
        push(@tids, $tid);
      }
      #Store exon coordinates
      my ($chr, $start, $end, $strand) = @cols[0,3,4,6];
      if ($chr =~ /^chr/){
        $chr =~ s/^chr//;
      }
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
print "Calculating intron coordinates...\n";
my %intron_coords_by_tid;
foreach my $tid (keys %exon_coords_by_tid){
  $intron_coords_by_tid{$tid}{'chr'} = $exon_coords_by_tid{$tid}{'chr'};
  $intron_coords_by_tid{$tid}{'strand'} = $exon_coords_by_tid{$tid}{'strand'};
  if (scalar @{$exon_coords_by_tid{$tid}{'starts'}} > 1){
    my @starts = sort {$a<=>$b} @{$exon_coords_by_tid{$tid}{'starts'}};
    my @ends = sort {$a<=>$b} @{$exon_coords_by_tid{$tid}{'ends'}};
    @{$intron_coords_by_tid{$tid}{'starts'}} = map {$_ + 1} @ends[0..$#ends-1];
    @{$intron_coords_by_tid{$tid}{'ends'}} = map {$_ - 1} @starts[1..$#starts];
  }
}


#Connect to database with reference annotation
print "Getting reference introns...\n";
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




print "Finding overlap with repeats...\n";
open (OUT, ">$outfile") or die "Can't open $outfile\n";

my @repeat_types;
if ($repeat_types){
  @repeat_types = split(/,/, $repeat_types);
}
my %seen_introns;
foreach my $tid (@tids){
  next unless $intron_coords_by_tid{$tid}{'starts'};
  my $repeat_overlap = 0;
  my $failed_intron_coord;
  my $repeat_feat_coord = "na";
  INT:for (my $i=0; $i<scalar(@{$intron_coords_by_tid{$tid}{'starts'}}); $i++){
    my $coord = $intron_coords_by_tid{$tid}{'chr'}.":".$intron_coords_by_tid{$tid}{'starts'}[$i]."-".$intron_coords_by_tid{$tid}{'ends'}[$i];
    my $coord2 = $intron_coords_by_tid{$tid}{'chr'}.":".$intron_coords_by_tid{$tid}{'starts'}[$i]."-".$intron_coords_by_tid{$tid}{'ends'}[$i].":".$intron_coords_by_tid{$tid}{'strand'};
    
    if ($ref_introns{$coord2}){ #ignore low score if intron was in the reference annotation
      #print $coord2."\n";
      next;
    }
      
    if (my $overlap = $seen_introns{$coord2}){
      if ($overlap eq "yes"){
        $repeat_overlap = 1;
        $failed_intron_coord = $coord2;
      }
    }
    else{
      my $iter = $repeat_tabix->query($coord);
      my $overlap = "no";
      if ($iter){
        my $overlapped_start;
        my $overlapped_end;
        while (my $line = $iter->next){
          my @cols = split(/\t/, $line);
          if (scalar @repeat_types){
            next unless scalar(grep {$_ eq $cols[6] or $_ eq $cols[7]} @repeat_types);
          }
          #next unless ($cols[6] eq "trf" or $cols[6] eq "Type I Transposons/SINE");
          if ($same_repeat_overlap){
            $overlapped_start = 0;
            $overlapped_end = 0;
            if ($cols[1]+1 <= $intron_coords_by_tid{$tid}{'starts'}[$i] and $cols[2] >= $intron_coords_by_tid{$tid}{'ends'}[$i]){
              $overlapped_start = 1;
              $overlapped_end = 1;
            }
          }
          else{
            if ($cols[1]+1 <= $intron_coords_by_tid{$tid}{'starts'}[$i] and $cols[2] >= $intron_coords_by_tid{$tid}{'starts'}[$i]){
              $overlapped_start = 1;
            }
            if ($cols[1]+1 <= $intron_coords_by_tid{$tid}{'ends'}[$i] and $cols[2] >= $intron_coords_by_tid{$tid}{'ends'}[$i]){
              $overlapped_end = 1;
            }
          }
          if ($overlapped_start and $overlapped_end){
            $overlap = "yes";
            $repeat_overlap = 1;
            $failed_intron_coord = $coord2;
            $repeat_feat_coord = $cols[3];
            last;
          }
        }
      }
      $seen_introns{$coord2} = $overlap;
    }

    if ($repeat_overlap){
      print OUT join("\t", $tid,
                          $failed_intron_coord,
                          $repeat_feat_coord,
                    )."\n";
      last INT;
    }
  }
}
close (OUT);

print "DONE\n";

