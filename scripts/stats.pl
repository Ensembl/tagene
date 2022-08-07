#!/usr/bin/perl -w


use strict;
use warnings;
use Getopt::Long;
use Bio::Otter::Server::Config;
use Bio::EnsEMBL::Utils::ConversionSupport;
$| = 1;

my $dataset_name;
my $outfile;
my $date;

&GetOptions(
            'dataset=s'  => \$dataset_name,
            'date=s'     => \$date,
            'out=s'      => \$outfile,
            );

open (OUT, ">$outfile") or die "Can't open $outfile:$!";
open (EXT, ">$outfile.tsv") or die "Can't open $outfile.tsv:$!";

print OUT join("\t", "gene_id",
                     "gene_symbol",
                     "gene_biotype",
                     "modified_date",
                     "",
                     "novel_transcripts",
                     "extended_transcripts",
                     "transcript_info",
                     "data_source",
              )."\n";

print EXT join("\t", "chromosome",
                     "start",
                     "end",
                     "gene_id",
                     "gene_symbol",
                     "gene_biotype",
                     "transcript_id",
                     "transcript_name",
                     "transcript_biotype",
                     "previous_transcript_biotype",
                     "novel/extended",
                     "data_source",
                     "CDS_end_not_found?",
                     "model_name",
                     "unique_CDS?",
                     "translation_length",
              )."\n";

#Connect to loutre database
#DataSet interacts directly with an otter database
my $dataset = Bio::Otter::Server::Config->SpeciesDat->dataset($dataset_name);
my $otter_dba = $dataset->otter_dba or die "can't get db adaptor\n";
$otter_dba->dbc->reconnect_when_lost(1);
my $sa = $otter_dba->get_SliceAdaptor();


my $new_gene_count = 0;
my $modified_gene_count = 0;
my $added_tr_count = 0;
my $extended_tr_count = 0;
my $multiple_tagene_ds_count = 0;
my $multiple_tagene_ms_count = 0;
foreach my $slice (@{$sa->fetch_all("chromosome")}){
  print $slice->seq_region_name."...";
  foreach my $gene (@{$slice->get_all_Genes}){
    next unless $gene->source =~ /havana/;
    my $gene_modif_date = Bio::EnsEMBL::Utils::ConversionSupport->date_format($gene->modified_date, "%y-%m-%d");
    if ($date){
      next unless $gene_modif_date ge $date;
    }
    #Get all gene's translations in order to report uniqueness of TAGENE CDS later
    my %tn_seqs;
    foreach my $transcript (@{$gene->get_all_Transcripts}){
      if ($transcript->translate){
        $tn_seqs{$transcript->stable_id} = $transcript->translation->seq;
      }
    }
    my $added_c= 0; #Added by TAGENE pipeline
    my $extended_c = 0; #Extended by TAGENE pipeline
    my $pre_c = 0; #Pre-existing transcript, not extended by TAGENE pipeline
    my $is_tagene_gene = 0; #Has a TAGENE transcript
    my %report;
    foreach my $transcript (@{$gene->get_all_Transcripts}){
      my $is_tagene = 0;
      my $tagene_datasets = 0;
      my $tagene_models = 0;
      my $extended_flag = 0;
      my %tsources;
      foreach my $att (@{$transcript->get_all_Attributes}){
        #Was the transcript created or extended by the TAGENE pipeline?
        if ($att->code eq "remark" and $att->value eq "TAGENE_transcript"
            and $transcript->transcript_author->name eq "tagene"){
          $is_tagene = 1;
          if ($date and 
              Bio::EnsEMBL::Utils::ConversionSupport->date_format($transcript->modified_date, "%y-%m-%d") lt $date){
            $is_tagene = 0;
          }
        }
      }
      if ($is_tagene){
        foreach my $att (@{$transcript->get_all_Attributes}){
          #Increase the count of different long read models that made up the transcript
          if ($att->code eq "hidden_remark" and $att->value =~ /^ID:.+(align|compmerge|NAM_TM|TM_)/){
            $tagene_models++;
          }
          #Increase the count of different long read datasets involved in making the transcript
          if ($att->code eq "remark" and 
              ($att->value =~ /^Assembled from PacBio CLS reads - .+\(GSE93848\)$/ or 
               $att->value eq "Assembled from PacBio CLS3 reads" or
               $att->value eq "Assembled from SLRseq reads (SRP049776)" or 
               $att->value eq "Assembled from RACEseq reads" or
               $att->value =~ /^Assembled from LRGASP.+reads.+/
              )){
            $tagene_datasets++;
            my ($dataset) = $att->value =~ /Assembled from (.+) reads/;
            $tsources{$dataset}++;
            $report{'sources'}{$dataset}++;
          }
        }

        #Was the transcript extended rather than created by TAGENE?  
        if ($date){
          if (Bio::EnsEMBL::Utils::ConversionSupport->date_format($transcript->created_date, "%y-%m-%d") lt $date
              and
              Bio::EnsEMBL::Utils::ConversionSupport->date_format($transcript->modified_date, "%y-%m-%d") ge $date
          ){
            $extended_flag = 1;
          }
        }
        else{
          #If attribute was not added by TAGENE, the transcript must have been created previously
          foreach my $att (@{$transcript->get_all_Attributes}){
            unless (($att->code eq "remark" and 
                    ($att->value =~ /^Assembled from PacBio CLS reads -.+\(GSE93848\)$/ or 
                     $att->value eq "Assembled from PacBio CLS3 reads" or
                     $att->value eq "Assembled from SLRseq reads (SRP049776)" or 
                     $att->value eq "Assembled from RACEseq reads" or
                     $att->value =~ /^Assembled from LRGASP.+reads.+/)) or
                    ($att->code eq "remark" and $att->value eq "not for VEGA") or
                    ($att->code eq "hidden_remark" and $att->value =~ /^ID: .+(align|compmerge|NAM_TM)/) or
                    ($att->code eq "hidden_remark" and $att->value =~ /^pacbio_capture_seq/) or
                    ($att->code eq "hidden_remark" and $att->value =~ /^SLR-seq/) or
                    ($att->code eq "hidden_remark" and $att->value =~ /^pacbio_raceseq/) or
                    ($att->code eq "hidden_remark" and $att->value =~ /^pacBio(SII)?:Cshl/) or 
                     $att->code eq "name" or 
                     $att->code eq "TAGENE_transcript" or                   
                    ($att->code eq "remark" and $att->value eq "TAGENE_transcript") or
                    ($att->code eq "remark" and $att->value eq "TAGENE_extended") or
                    ($att->code eq "hidden_remark" and $att->value eq "Platinum set")                   
                   ){
              $extended_flag = 1;
            }
          }
        }
      }
      if ($is_tagene){
        $is_tagene_gene = 1;
        #Multiple datasets
        if ($tagene_datasets > 1){
          $multiple_tagene_ds_count++;
        }
        if ($tagene_models > 1){
          $multiple_tagene_ms_count++; 
          print "MULTI: ".$transcript->stable_id."\n";
        }
        #Extended existing one? - look at evidence list
        if (scalar @{$transcript->evidence_list} or $extended_flag == 1){
          $extended_c++;
          $report{$transcript->stable_id} = "extended";
        }
        else{
          $added_c++;
          $report{$transcript->stable_id} = "novel";
        }
        $report{$transcript->stable_id} .= "-".$transcript->biotype;
      }
      else{
        $pre_c++;
      }
      
      if ($is_tagene){
        my $tn_length;
        my $is_unique_cds = "NA";
        if ($transcript->translation){
          $tn_length = $transcript->translation->length;
          my $tnseq = $transcript->translation->seq;
          $is_unique_cds = "yes";
          foreach my $tid (grep {$_ ne $transcript->stable_id} keys %tn_seqs){
            if ($tn_seqs{$tid} eq $tnseq){
              $is_unique_cds = "no";
              last;
            }
            elsif ($tn_seqs{$tid} =~ /$tnseq/){
              $is_unique_cds = "partial";
            }
          }
        }
        
        print EXT join("\t",
                        $transcript->seq_region_name,
                        $transcript->seq_region_start,
                        $transcript->seq_region_end,
                        $gene->stable_id, 
                        ($gene->get_all_Attributes('name')->[0]->value || "NA"),
                        $gene->biotype, 
                        $transcript->stable_id,
                        ($transcript->get_all_Attributes('name')->[0]->value || "NA"),
                        $transcript->biotype,
                        (previous_biotype($transcript) || "NA"),
                        ($extended_flag ? "extended" : "novel"),
                        join(", ", keys %tsources),
                        scalar(@{$transcript->get_all_Attributes('cds_end_NF')}) ? "cds_end_NF" : "NA",
                        join(", ", map {$_->value} grep {$_->value =~ /^ID:.+(align|compmerge|NAM_TM|TM_)/} @{$transcript->get_all_Attributes('hidden_remark')}),
                        $is_unique_cds,
                        ($tn_length || "NA")
                      )."\n";
      
      }
    }
    if ($is_tagene_gene){
      if ($added_c > 0 and $extended_c == 0 and $pre_c == 0){
        $new_gene_count++;
        $report{$gene->stable_id} = "new";
      }
      else{
        $modified_gene_count++;
        $report{$gene->stable_id} = "modified";
      }
    }
    $added_tr_count += $added_c;
    $extended_tr_count += $extended_c;
    if (scalar keys %report){
      print OUT join(", ", map {$_} grep {/ENS(MUS)?G|OTT(HUM|MUS)G/} keys %report)."\t".
                ($gene->get_all_Attributes('name')->[0]->value || " ")."\t".
                $gene->biotype."\t".
                $gene_modif_date."\t".
                join(", ", map {$report{$_}} grep {/ENS(MUS)?G/} keys %report)."\t".
                scalar(grep {/novel/} values(%report))."\t". 
                scalar(grep {/extended/} values(%report))."\t". 
                join(", ", map {$_.":".$report{$_}} grep {/ENS(MUS)?T/} keys %report)."\t". 
                join(", ", keys %{$report{'sources'}})."\n";
    }
  }
}

print OUT "\nNew genes: $new_gene_count\n".
          "Modified genes: $modified_gene_count\n".
          "Novel transcripts: $added_tr_count\n".
          "Extended transcripts: $extended_tr_count\n".
          "Transcripts with >1 TAGENE datasets: $multiple_tagene_ds_count\n".
          "Transcripts with >1 TAGENE models: $multiple_tagene_ms_count\n";

close (OUT);
close (EXT);

######

sub previous_biotype {
  my $transcript = shift;
  my $prev_biotype;
  my $sth;
  if ($date){
    $sth = $otter_dba->dbc->prepare("select biotype from transcript where stable_id=? and modified_date < ? order by modified_date desc limit 1");
    $sth->execute($transcript->stable_id, $date);
  }
  else{
    $sth = $otter_dba->dbc->prepare("select biotype from transcript where stable_id=? and version=?");
    $sth->execute($transcript->stable_id, $transcript->version - 1);
  }
  while (my @cols = $sth->fetchrow_array()) {
    $prev_biotype = $cols[0];
  }
  $sth->finish;
  return $prev_biotype;
}





