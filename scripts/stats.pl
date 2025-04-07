#!/usr/bin/perl -w


use strict;
use warnings;
use Getopt::Long;
use Bio::Otter::Server::Config;
use Bio::EnsEMBL::Utils::ConversionSupport;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
$| = 1;


my $dataset_name;
my $host;
my $port;
my $user;
my $dbname;
my $outfile;
my $date;
my $remark;
my $only_chr;

&GetOptions(
            'dataset=s'  => \$dataset_name,
            'host=s'     => \$host,
            'port=i'     => \$port,
            'user=s'     => \$user,
            'dbname=s'   => \$dbname,
            'date=s'     => \$date,
            'remark=s'   => \$remark,
            'chr=s'      => \$only_chr,
            'out=s'      => \$outfile,
            );

open (G, ">$outfile.gene.txt") or die "Can't open $outfile.gene.txt:$!";
open (T, ">$outfile.transcript.txt") or die "Can't open $outfile.transcript.txt:$!";

print G join("\t", "chromosome",
                     "start",
                     "end",
                     "strand",
                     "gene_id",
                     "gene_symbol",
                     "gene_biotype",
                     "modified_date",
                     "new/modified",
                     "single_exon",
                     "transcript_count",
                     "novel_transcripts",
                     "extended_transcripts",
                     "transcript_info",
                     "data_source",
                     "novel_exonic_seq",
                     "novel_cds_seq",
                     "novel_splice_sites",
                     "novel_junctions",     
              )."\n";

print T join("\t", "chromosome",
                     "start",
                     "end",
                     "strand",
                     "gene_id",
                     "gene_symbol",
                     "gene_biotype",
                     "transcript_id",
                     "transcript_name",
                     "transcript_biotype",
                     "previous_transcript_biotype",
                     "novel/extended",
                     "data_source",
                     "model_name",
                     "exon_count",
                     "novel_CDS",
                     "CDS_end_not_found",
                     "SQANTI3_cat",
                     "novel_exonic_seq",
                     "novel_cds_seq",
                     "novel_splice_sites",
                     "novel_junctions",
                     "translation_length",
                     #"translation_seq",
              )."\n";

#Connect to the havana database
my $otter_dba;
if ($dataset_name){
  #DataSet interacts directly with an otter database
  my $dataset = Bio::Otter::Server::Config->SpeciesDat->dataset($dataset_name);
  $otter_dba = $dataset->otter_dba or die "can't get db adaptor\n";
  $otter_dba->dbc->reconnect_when_lost(1);
}
else{
  $otter_dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host => $host,
    -port => $port,    
    -user => $user,
    -dbname => $dbname,
  );
}
my $sa = $otter_dba->get_SliceAdaptor();

my $new_gene_count = 0;
my $new_se_gene_count = 0;
my $modified_gene_count = 0;
my $added_tr_count = 0;
my $extended_tr_count = 0;
my $multiple_tagene_ds_count = 0;
my $multiple_tagene_ms_count = 0;
foreach my $slice (@{$sa->fetch_all("toplevel")}){
  next unless $slice->coord_system->is_default;
  print $slice->seq_region_name."...";
  if ($only_chr){
    next unless $slice->seq_region_name eq $only_chr;
  }
  foreach my $gene (@{$slice->get_all_Genes}){
    next unless $gene->source =~ /ensembl|havana/;
    my $gene_modif_date = Bio::EnsEMBL::Utils::ConversionSupport->date_format($gene->modified_date, "%y-%m-%d");
    if ($date){
      next unless $gene_modif_date ge $date;
    }
    #Total transcript count
    my $total_tr_count = 0;
    foreach my $transcript (@{$gene->get_all_Transcripts}){
      unless ($transcript->biotype eq "artifact" or grep {$_->value eq "not for VEGA"} @{$transcript->get_all_Attributes("remark")}){
        $total_tr_count++;
      }
    }
    


    my $added_c= 0; #Added by TAGENE pipeline
    my %added_c_b; #Added, by biotype
    my $extended_c = 0; #Extended by TAGENE pipeline
    my $pre_c = 0; #Pre-existing transcript, not extended by TAGENE pipeline
    my $is_tagene_gene = 0; #Has a TAGENE transcript
    my $is_single_exon_gene = 0;
    my %report;
    my @novel_transcripts;
    my @current_transcripts;
    foreach my $transcript (@{$gene->get_all_Transcripts}){
      my $is_tagene = 0;
      if (scalar(@{$gene->get_all_Transcripts}) == 1 and scalar(@{$transcript->get_all_Exons}) == 1){
        $is_single_exon_gene = 1;
      }
      #Was the transcript created or extended by the TAGENE pipeline?
      if (scalar grep {$_->value eq "TAGENE_transcript"} @{$transcript->get_all_Attributes('remark')}){
        $is_tagene = 1;
        #Check modification date if provided
        if ($date and Bio::EnsEMBL::Utils::ConversionSupport->date_format($transcript->modified_date, "%y-%m-%d") lt $date){
          $is_tagene = 0;
        }
        #Check TAGENE run remark if provided
        if ($remark and scalar(grep {$_->value eq $remark} @{$transcript->get_all_Attributes('remark')}) == 0){
          $is_tagene = 0;
        }
      }
      if ($is_tagene){
        push(@novel_transcripts, $transcript);
      }
      else{
        unless (scalar grep {$_->value eq "not for VEGA"} @{$transcript->get_all_Attributes("remark")} or $transcript->biotype eq "artifact"){
          push(@current_transcripts, $transcript);
        }
      }
    }

    #Get data from current annotation for comparison with novel transcripts
    #Exon coverage
    my %current_exonic_sequence;
    foreach my $current_tr (@current_transcripts){
      foreach my $exon (@{$current_tr->get_all_Exons}){
        for (my $pos = $exon->seq_region_start; $pos <= $exon->seq_region_end; $pos++){
          $current_exonic_sequence{$pos} = 1;
        }
      }
    }
    #CDS coverage
    my %current_cds_sequence;
    foreach my $current_tr (@current_transcripts){
      foreach my $cds (@{$current_tr->get_all_CDS}){
        for (my $pos = $cds->seq_region_start; $pos <= $cds->seq_region_end; $pos++){
          $current_cds_sequence{$pos} = 1;
        }
      }
    }
    #Intron chain
    my %current_intron_chains;
    foreach my $current_tr (@current_transcripts){
      my $intron_chain = join(":", map {$_->seq_region_start."-".$_->seq_region_end} @{$current_tr->get_all_Introns});
      $current_intron_chains{$intron_chain} = 1;
    }    
    #CDS exon chain
    my %current_cds_exon_chains;
    foreach my $current_tr (@current_transcripts){
      if ($current_tr->translate){
        my $cds_chain = join(":", map {$_->seq_region_start."-".$_->seq_region_end} @{$current_tr->get_all_CDS});
        $current_cds_exon_chains{$cds_chain} = 1;
      }
    }
    #Translations
    my %current_translation_seqs;
    foreach my $transcript (@current_transcripts){
      if ($transcript->translate){
        $current_translation_seqs{$transcript->translation->seq} = 1;
      }
    }  
    #Splice junctions
    my %known_junctions;
    foreach my $current_tr (@current_transcripts){
      foreach my $intron (@{$current_tr->get_all_Introns}){
        $known_junctions{$intron->seq_region_start."-".$intron->seq_region_end} = 1;
      }
    }
    #Splice sites
    my %known_donor_splice_sites;
    my %known_acceptor_splice_sites;
    foreach my $current_tr (@current_transcripts){
      foreach my $intron (@{$current_tr->get_all_Introns}){
        $known_donor_splice_sites{$intron->seq_region_start} = 1; #no need to check strand
        $known_acceptor_splice_sites{$intron->seq_region_end} = 1;
      }
    }  



    ####
    my $novel_exon_coverage = 0;
    my %novel_exon_coverage;
    my $novel_cds_coverage = 0;
    my %novel_cds_coverage;
    my %novel_splice_sites;
    my %novel_splice_junctions;


    foreach my $transcript (@novel_transcripts){
      my $tagene_datasets = 0;
      my $tagene_models = 0;
      my $extended_flag = 0;
      my %tsources;
      my $is_tagene = 1;
      if ($is_tagene){
        foreach my $att (@{$transcript->get_all_Attributes}){
          #Increase the count of different long read models that made up the transcript
          if ($att->code eq "hidden_remark" and $att->value =~ /^ID:.+(align|compmerge|NAM_TM|TM_|PB|anchIC)/){
            $tagene_models++;
          }
          #Increase the count of different long read datasets involved in making the transcript
          if ($att->code eq "remark" and 
              #($att->value =~ /^Assembled from PacBio CLS reads - .+\(GSE93848\)$/ or 
               #$att->value eq "Assembled from PacBio CLS3 reads" or
               #$att->value eq "Assembled from SLRseq reads (SRP049776)" or 
               #$att->value eq "Assembled from RACEseq reads" or
               #$att->value =~ /^Assembled from LRGASP.+reads.+/ or
               #$att->value =~ /^Assembled from PacBio reads/
               ($att->value =~ /Assembled from (.+) reads/ or
               $att->value eq "CLS3 project" or
               $att->value =~ /^(DanaFarber|ENCODE|LRGASP)$/)
              ){
            $tagene_datasets++;
            #my ($dataset) = $att->value =~ /Assembled from (.+) reads/;
            my $dataset = $1;
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
                     $att->value =~ /^Assembled from LRGASP.+reads.+/) or
                     $att->value =~ /^Assembled from PacBio reads/ or
                     $att->value eq "CLS3 project"
                     ) or
                    ($att->code eq "remark" and $att->value eq "not for VEGA") or
                    ($att->code eq "hidden_remark" and $att->value =~ /^ID:.+(align|compmerge|NAM_TM|anchIC)/) or
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
          #print "MULTI: ".$transcript->stable_id."\n";
        }
        #Extended existing one? - look at evidence list
        if ($extended_flag == 1){
          $extended_c++;
          $report{$transcript->stable_id} = "extended";
        }
        else{
          $added_c++;
          $report{$transcript->stable_id} = "novel";
          #count by biotype
          $added_c_b{$transcript->biotype}++;
        }
        $report{$transcript->stable_id} .= "-".$transcript->biotype;
      }
      else{
        $pre_c++;
      }
      
      if ($is_tagene){
        my $tn_length;
        my $is_unique_cds = "NA";
        my $tn_seq;
        if ($transcript->translation){
          $tn_length = $transcript->translation->length;
          $tn_seq = $transcript->translation->seq;
        }

        my $cds_str = join(":", map {$_->seq_region_start."-".$_->seq_region_end} @{$transcript->get_all_CDS});

        #Novel exon coverage
        my $t_novel_exon_coverage = 0;
        foreach my $exon (@{$transcript->get_all_Exons}){
          for (my $pos = $exon->seq_region_start; $pos <= $exon->seq_region_end; $pos++){
            unless ($current_exonic_sequence{$pos}){
              $t_novel_exon_coverage++;
              $novel_exon_coverage{$pos} = 1;
            }
          }
        }

        #Novel CDS coverage
        my $t_novel_cds_coverage = 0;
        foreach my $cds (@{$transcript->get_all_CDS}){
          for (my $pos = $cds->seq_region_start; $pos <= $cds->seq_region_end; $pos++){
            unless ($current_cds_sequence{$pos}){
              $t_novel_cds_coverage++;
              $novel_cds_coverage{$pos} = 1;
            }
          }
        }

        #Novel CDS (exon chain)
        my $is_novel_cds = 0;
        if ($transcript->translate){
          my $cds_chain = join(":", map {$_->seq_region_start."-".$_->seq_region_end} @{$transcript->get_all_CDS});
          unless ($current_cds_exon_chains{$cds_chain}){
            $is_novel_cds = 1;
          }
        }

        #Novel splice sites
        my $t_novel_splice_sites = 0;
        foreach my $intron (@{$transcript->get_all_Introns}){
          unless ($known_donor_splice_sites{$intron->seq_region_start}){
            $t_novel_splice_sites++;
            $novel_splice_sites{$intron->seq_region_start} = 1;
          }
          unless ($known_acceptor_splice_sites{$intron->seq_region_end}){
            $t_novel_splice_sites++;
            $novel_splice_sites{$intron->seq_region_end} = 1;
          }
        }
 
        #Novel splice junctions
        my $t_novel_splice_junctions = 0;
        foreach my $intron (@{$transcript->get_all_Introns}){
          unless ($known_junctions{$intron->seq_region_start."-".$intron->seq_region_end}){
            $t_novel_splice_junctions++;
            $novel_splice_junctions{$intron->seq_region_start."-".$intron->seq_region_end} = 1;
          }
        }
      

        
        print T join("\t",
                        $transcript->seq_region_name,
                        $transcript->seq_region_start,
                        $transcript->seq_region_end,
                        $transcript->seq_region_strand,
                        $gene->stable_id, 
                        ($gene->get_all_Attributes('name')->[0] ? $gene->get_all_Attributes('name')->[0]->value : "NA"),
                        $gene->biotype, 
                        $transcript->stable_id,
                        ($transcript->get_all_Attributes('name')->[0] ? $transcript->get_all_Attributes('name')->[0]->value : "NA"),
                        $transcript->biotype,
                        (previous_biotype($transcript) || "NA"),
                        ($extended_flag ? "extended" : "novel"),
                        join(", ", keys %tsources),
                        join(", ", map {$_->value} grep {$_->value =~ /^ID:.+(align|compmerge|NAM_TM|TM_|PB|anchIC|anchUC|.+)/} @{$transcript->get_all_Attributes('hidden_remark')}),
                        scalar(@{$transcript->get_all_Exons}),
                        $is_novel_cds ? "yes" : "no",
                        scalar(@{$transcript->get_all_Attributes('cds_end_NF')}) ? "cds_end_NF" : "NA",
                        sqanti3_category($transcript, \@current_transcripts),
                        $t_novel_exon_coverage,
                        $t_novel_cds_coverage,
                        $t_novel_splice_sites,
                        $t_novel_splice_junctions,
                        ($tn_length || "NA"),
                       # ($tn_seq || "NA"),
                      )."\n";
      
      }
    }
    if ($is_tagene_gene){
      if (scalar @current_transcripts == 0){
        $new_gene_count++;
        $report{$gene->stable_id} = "new";
        if ($is_single_exon_gene){
          $new_se_gene_count++;
        }
      }
      else{
        $modified_gene_count++;
        $report{$gene->stable_id} = "modified";
      }
    }
    $added_tr_count += $added_c;
    $extended_tr_count += $extended_c;
    if (scalar keys %report){
      print G join("\t",
                  $gene->seq_region_name,
                  $gene->seq_region_start,
                  $gene->seq_region_end,
                  $gene->seq_region_strand,
                  join(", ", map {$_} grep {/ENS([A-Z]{3})?G|OTT([A-Z]{3})G/} keys %report),
                  ($gene->get_all_Attributes('name')->[0]->value || " "),
                  $gene->biotype,
                  $gene_modif_date,
                  join(", ", map {$report{$_}} grep {/ENS([A-Z]{3})?G/} keys %report),
                  ($is_single_exon_gene ? "yes" : "no"),
                  $total_tr_count,
                  scalar(grep {/novel/} values(%report)),
                  scalar(grep {/extended/} values(%report)),
                  join(", ", map {$_.":".$report{$_}} grep {/ENS([A-Z]{3})?T/} keys %report),
                  join(", ", keys %{$report{'sources'}}),
                  scalar keys %novel_exon_coverage,
                  scalar keys %novel_cds_coverage,
                  scalar keys %novel_splice_sites,
                  scalar keys %novel_splice_junctions,
                  join("\t", $added_c_b{"protein_coding"}||0, $added_c_b{"nonsense_mediated_decay"}||0),
                )."\n";
    }
  }
}

print G "\nNew genes: $new_gene_count\n".
          "New single exon genes: $new_se_gene_count\n".
          "Modified genes: $modified_gene_count\n".
          "Novel transcripts: $added_tr_count\n".
          "Extended transcripts: $extended_tr_count\n".
          "Transcripts with >1 TAGENE datasets: $multiple_tagene_ds_count\n".
          "Transcripts with >1 TAGENE models: $multiple_tagene_ms_count\n";

close (G);
close (T);

######

sub previous_biotype {
  my $transcript = shift;
  my $prev_biotype;
  my $sth;
  if ($date){
    $sth = $otter_dba->dbc->prepare("select biotype from transcript join seq_region using (seq_region_id)
                                      where stable_id = ? and modified_date < ? and name = ? and seq_region_start < ? and seq_region_end > ?
                                      order by transcript_id desc limit 1");
    $sth->execute($transcript->stable_id, $date, $transcript->slice->seq_region_name, $transcript->seq_region_end, $transcript->seq_region_start);
  }
  else{
    $sth = $otter_dba->dbc->prepare("select biotype from transcript join seq_region using (seq_region_id)
                                      where stable_id=? and version=? and name = ? and seq_region_start < ? and seq_region_end > ?");
    $sth->execute($transcript->stable_id, $transcript->version - 1, $transcript->slice->seq_region_name, $transcript->seq_region_end, $transcript->seq_region_start);
  }
  while (my @cols = $sth->fetchrow_array()) {
    $prev_biotype = $cols[0];
  }
  $sth->finish;
  return $prev_biotype;
}


sub sqanti3_category {
  my ($transcript, $current_transcripts) = @_;
  my $category;
  my $subcategory;
  my %known_intron_chains;
  my %known_junctions;
  my %known_donor_splice_sites;
  my %known_acceptor_splice_sites;

  foreach my $current_tr (@$current_transcripts){
    my $c_introns = $current_tr->get_all_Introns;
    my $c_intron_chain = join(":", map {$_->seq_region_start."-".$_->seq_region_end} @{$c_introns});
    $known_intron_chains{$c_intron_chain} = 1;
    foreach my $intron (@$c_introns){
      $known_junctions{$intron->seq_region_start."-".$intron->seq_region_end} = 1;
      $known_donor_splice_sites{$intron->seq_region_start} = 1;
      $known_acceptor_splice_sites{$intron->seq_region_end} = 1;
    }
  }

  #Multi-exon transcripts
  if (scalar @{$transcript->get_all_Exons} > 1){
    my $introns = $transcript->get_all_Introns;
    my $intron_chain = join(":", map {$_->seq_region_start."-".$_->seq_region_end} @$introns);
    my %junctions;
    my %donor_splice_sites;
    my %acceptor_splice_sites;
    foreach my $intron (@$introns){
      $junctions{$intron->seq_region_start."-".$intron->seq_region_end} = 1;
      $donor_splice_sites{$intron->seq_region_start} = 1;
      $acceptor_splice_sites{$intron->seq_region_end} = 1;
    }

    if ($known_intron_chains{$intron_chain}){
      $category = "FSM";
    }
    else{
      foreach my $ic (keys %known_intron_chains){
        if ($known_intron_chains{$ic} =~ /$intron_chain/){
          $category = "ISM";
          if ($known_intron_chains{$ic} =~ /^$intron_chain/){
            $subcategory = "5' fragment";
          }
          elsif ($known_intron_chains{$ic} =~ /$intron_chain$/){
            $subcategory = "3' fragment";
          }
          else{
            $subcategory = "internal fragment";
          }
        }
        # ...ISM intron retention
      }
      unless ($category){
        foreach my $ss (keys %donor_splice_sites){
          unless ($known_donor_splice_sites{$ss}){
            $category = "NNC";
          }
        }
        foreach my $ss (keys %acceptor_splice_sites){
          unless ($known_acceptor_splice_sites{$ss}){
            $category = "NNC";
          }
        }
      }
      unless ($category){
        foreach my $sj (keys %junctions){
          unless ($known_junctions{$sj}){
            $category = "NIC";
            $subcategory = "combination of known splice sites"
          }
        }
      }
      unless ($category){
        $category = "NIC";
        $subcategory = "combination of known SJs"       
      }
      #...NIC intron retention
    }
  }
  #Single-exon transcripts
  else{
    if ($known_acceptor_splice_sites{$transcript->seq_region_start-1} and $known_donor_splice_sites{$transcript->seq_region_end+1}){
      foreach my $ic (keys %known_intron_chains){
        my ($start, $end) = split($ic, "-");
        if ($start > $transcript->seq_region_start and $end < $transcript->seq_region_end){
          $category = "NIC";
          $subcategory = "mono-exon by intron retention";
        }
      }
    }
    else{
      foreach my $current_tr (@$current_transcripts){
        foreach my $exon (@{$transcript->get_all_Exons}){
          if ($exon->seq_region_start == $current_tr->seq_region_end and $exon->seq_region_end == $current_tr->seq_region_start){
            $category = "ISM";
            $subcategory = "mono-exon";
          }
        }
      }
    }
    unless ($category){
      CT:foreach my $current_tr (@$current_transcripts){
        foreach my $exon (@{$transcript->get_all_Exons}){
          if ($exon->seq_region_start < $current_tr->seq_region_end and $exon->seq_region_end > $current_tr->seq_region_start){
            $category = "Genic genomic";
            last CT;
          }
        }
      }
    }
    unless ($category){
      CT:foreach my $current_tr (@$current_transcripts){
        foreach my $intron (@{$transcript->get_all_Introns}){
          if ($intron->seq_region_start <= $current_tr->seq_region_start and $intron->seq_region_end >= $current_tr->seq_region_end){
            $category = "Genic intron";
            last CT;
          }
        }
      }
    }
  }
  if ($subcategory){
    $category .= " - ".$subcategory;
  }
  return $category;
}

