#!/usr/bin/perl -w

#This script merges redundant Otter transcripts from a list of pairs of stable ids.
#It does not guess what transcripts need to be merged.
#In the input file, the first column corresponds to the transcript that will be deleted
#from the database as it is merged into the transcript indicated by the second column.
#Written originally to fix redundant annotation created by the TAGENE pipeline
#Usage: perl merge_redundant_transcripts .pl -dataset human_test -in list.txt


use strict;
use warnings;
use Getopt::Long;
use LoutreWrite::AnnotUpdate qw( exon_novelty intron_novelty can_be_merged merge_transcripts $OTTER_DBA $OTTER_SA $WRITE );
use Bio::Otter::Lace::Defaults;
use Bio::Otter::Server::Config;
use Bio::EnsEMBL::Attribute;
$| = 1;


my $dataset_name;
my $infile;
my $logfile;

&GetOptions(
            'dataset=s'  =>  \$dataset_name,
            'in=s'       =>  \$infile,
            'log=s'      =>  \$logfile,
            'write!'     =>  \$WRITE
            );


#Connect to loutre database
#DataSet interacts directly with an otter database
my $dataset = Bio::Otter::Server::Config->SpeciesDat->dataset($dataset_name);
$OTTER_DBA = $dataset->otter_dba or die "can't get db adaptor\n";
$OTTER_DBA->dbc->reconnect_when_lost(1);
$OTTER_SA = $OTTER_DBA->get_SliceAdaptor();
my $ga = $OTTER_DBA->get_GeneAdaptor();
my $ta = $OTTER_DBA->get_TranscriptAdaptor();

  
open (OUT, ">$logfile") or die "Can't open $logfile: $!";

open (IN, $infile) or die "Can't open $infile: $!";
while (<IN>){
  chomp;
  my ($tid1, $tid2) = split(/\t/);
  my ($tr1, $tr2);
  eval {
    $tr1 = $ta->fetch_by_stable_id($tid1);
    $tr2 = $ta->fetch_by_stable_id($tid2);
  };
  if ($@){
    if ($@ =~ /$tid1/){
      print OUT "Transcript $tid1 not found in database\n";
    }
    if ($@ =~ /$tid2/){
      print OUT "Transcript $tid2 not found in database\n";
    }
    next;
  }

  print OUT "\n\n";
  print OUT join("\t",
            $tr1->stable_id,
            scalar(@{$tr1->get_all_Exons}),
            (is_tagene($tr1) ? "TAGENE" : "non-TAGENE"),
            "e:".exon_novelty($tr1, $tr2),
            "i:".intron_novelty($tr1, $tr2),
            "m:".can_be_merged($tr1, $tr2),
            )."\n";
  print OUT join("\t",
            $tr2->stable_id,
            scalar(@{$tr2->get_all_Exons}),
            (is_tagene($tr2) ? "TAGENE" : "non-TAGENE"),
            "e:".exon_novelty($tr2, $tr1),
            "i:".intron_novelty($tr2, $tr1),
            "m:".can_be_merged($tr2, $tr1),
            )."\n";
  
  
  if (scalar(@{$tr1->get_all_Exons})==2 and scalar(@{$tr2->get_all_Exons})==2 and
      !($tr1->translate or $tr2->translate) and
      exon_novelty($tr1, $tr2) and !(intron_novelty($tr1, $tr2) or intron_novelty($tr2, $tr1)) and can_be_merged($tr1, $tr2)){
      

    #Get Otter region for the gene
    my $gene = $tr1->get_Gene;
    my $local_server = Bio::Otter::Server::Support::Local->new( otter_dba => $OTTER_DBA );
    $local_server->set_params(
        dataset => $dataset_name,
        chr     => $gene->seq_region_name,
        start   => $gene->start,   
        end     => $gene->end,
        cs      => "chromosome",
        csver   => "Otter"
    );

    my $region_action = Bio::Otter::ServerAction::Script::Region->new_with_slice($local_server);
    my $region = $region_action->get_region; #Bio::Vega::Region
    print "\nFOUND ".scalar($region->genes)." GENES IN REGION ".$gene->seq_region_name.":".$gene->start."\n";

    #Reset gene coordinates to the start of the region slice,
    #otherwise the gene coordinates in loutre will be wrong
    my $slice_offset = $region->slice->start - 1;
    print "OFFSET=$slice_offset\n";
    foreach my $tr (@{$gene->get_all_Transcripts}){
        foreach my $exon (@{$tr->get_all_Exons}){
            $exon->start($exon->start - $slice_offset);
            $exon->end(  $exon->end   - $slice_offset);
            $exon->slice($region->slice);
        }
    }
      
      
      my $merged_transcript = merge_transcripts($tr1, $tr2);

      foreach my $exon (@{$merged_transcript->get_all_Exons}){
        $exon->start($exon->start - $slice_offset); print $exon->start."\n";
        $exon->end(  $exon->end   - $slice_offset); print $exon->end."\n";
        $exon->slice($region->slice);
      }

      foreach my $g ($region->genes) {
        foreach my $ts (@{$g->get_all_Transcripts}) {
          if ($ts->stable_id eq $merged_transcript->stable_id){
            $ts->flush_Exons; 
            foreach my $exon (@{$merged_transcript->get_all_Exons}){
              $ts->add_Exon($exon);
            }
                
            #Add remarks (except 'not for VEGA') and hidden remarks from the novel transcript
            foreach my $att (@{$merged_transcript->get_all_Attributes}){
              next if ($att->code eq "remark" and $att->value eq "not for VEGA");
              next if ($att->code eq "name");
              #Append read names to existing remark
              if ((($att->code eq "hidden_remark" and $att->value =~ /^pacbio_capture_seq_\w+ : .+;/) or
                   ($att->code eq "hidden_remark" and $att->value =~ /^SLR-seq_\w+ : .+;/) or
                   ($att->code eq "hidden_remark" and $att->value =~ /^pacbio_raceseq_\w+ : .+;/)) and
                   (scalar(@{$ts->get_all_Attributes('TAGENE_transcript')}) > 0)){
                      my $appended = 0;
                      foreach my $att2 (@{$ts->get_all_Attributes('hidden_remark')}){
                        if ($att2->value =~ /^pacbio_capture_seq_\w+ : .+;/ or
                            $att2->value =~ /^SLR-seq_\w+ : .+;/ or
                            $att2->value =~ /^pacbio_raceseq_\w+ : .+;/){
                              $att2->value($att2->value." ".$att->value);
                              $appended = 1;
                              last;
                        }
                      }
                      #If no read name remark yet, create one
                      unless ($appended){
                        $ts->add_Attributes($att);
                      }
              }
              else{
                $ts->add_Attributes($att);
              }
            }

            #If coding gene, try to assign a CDS and change the biotype accordingly
            #assign_cds_to_transcripts($ts, $g, $slice_offset);
                                
            #Change authorship
            $ts->transcript_author($merged_transcript->transcript_author);
                                
            print OUT "TR: Will merge transcript ".$tr1->stable_id." (".$tr1->biotype.") with transcript ".$ts->stable_id." (".$ts->biotype.") in gene ".$g->stable_id."\n";

          }
          #Make the first TAGENE transcript non-current
          if ($ts->stable_id eq $tr1->stable_id){
            my $nfv_remark = Bio::EnsEMBL::Attribute->new(
                                -code => 'remark',
                                -value => 'not for VEGA'
                             );
            #$ts->add_Attributes($nfv_remark);
           # $ts->is_current(0);
            $g->remove_Transcript($ts);
            #$ts->biotype('obsolete');
          }
        }
      }
    #Write region
    my $g_msg = " ";
    if ($WRITE){
        $local_server->authorized_user($gene->gene_author->name); # preserve authorship
        my $n = 0;
        while (($g_msg eq " " or $g_msg =~ /write failed/) and $n<10){
            $g_msg = LoutreWrite::AnnotUpdate::write_gene_region($region_action, $region);
            sleep(int(rand(3)));
            $n++;
        }
        print OUT $g_msg."\n";
    }
  }
  else{
    print OUT "Not merging transcripts\n";
  }
}
close (IN);

close (OUT);


sub is_tagene {
  my $tr = shift;
  if (scalar grep {$_->value eq "TAGENE_transcript"} @{$tr->get_all_Attributes('remark')}){
    return 1;
  }
  return 0;
}


__END__

# foreach my $slice (@{$sa->fetch_all("chromosome", "Otter")}){
#   next unless $slice->seq_region_name =~ /-38$/;
#   if ($only_chr){
#     next unless $slice->seq_region_name eq "chr$only_chr-38";
#   }
#   print $slice->seq_region_name."\n";
#   foreach my $gene (@{$slice->get_all_Genes}){
#     if ($only_gene){
#       next unless $gene->stable_id eq $only_gene;
#     }
#     my @trs = @{$gene->get_all_Transcripts};
#     for (my $i=0; $i<scalar(@trs); $i++){
#       my $tr = $trs[$i];
#       if (scalar grep {$_->value eq "TAGENE_transcript"} @{$tr->get_all_Attributes('remark')}){
#         #Check for redundant transcripts (where one could be merged into the other)
#         for (my $j=$i+1; $j<scalar(@trs); $j++){
#           my $db_tr = $trs[$j];
#           next if $db_tr->biotype eq "artifact";
#           next if scalar(grep {$_->value eq "not for VEGA"} @{$db_tr->get_all_Attributes('remark')}) and $db_tr->biotype ne "comp_pipe";
#           next if scalar(grep {$_->value eq "comp_pipe_rejected"} @{$db_tr->get_all_Attributes('remark')});
#           next if scalar @{$db_tr->get_all_Introns} == 0;
#           
#           if ((exon_novelty($tr, $db_tr) or exon_novelty($db_tr, $tr)) and can_be_merged($tr, $db_tr)){
#             print OUT join("\t", "Redundant_transcript", 
#                                  $tr->stable_id, scalar(@{$tr->get_all_Exons}), $tr->biotype,
#                                  $db_tr->stable_id, scalar(@{$db_tr->get_all_Exons}), $db_tr->biotype,
#                                  $gene->stable_id, $gene->biotype)."\n";
#           }
#         }
#         
#         #Check for NMDs with an ORF < 35 aa
#         if ($tr->biotype eq "nonsense_mediated_decay" and $tr->translation->length < 35){
#           print OUT join("\t", "Short_NMD_ORF", $tr->stable_id, $tr->translation->length)."\n";
#         }
#       }
#     }
#   }
# }
# 
# close (OUT);
# 
# 
#     if ($mode eq "update"){
#         my $host_gene = $dba->get_GeneAdaptor->fetch_by_stable_id($gene->stable_id);
#         $padded_start = min($gene->start, $host_gene->start);
#         $padded_end = max($gene->end, $host_gene->end);
#     }
# 
#     #Get Otter region for the gene
#     my $local_server = Bio::Otter::Server::Support::Local->new( otter_dba => $dba );
#     $local_server->set_params(
#         dataset => $dataset_name,
#         chr     => $gene->seq_region_name,
#         start   => $padded_start || $gene->start,   
#         end     => $padded_end || $gene->end,
#         cs      => "chromosome",
#         csver   => "Otter"
#     );
# 
#     my $region_action = Bio::Otter::ServerAction::Script::Region->new_with_slice($local_server);
#     my $region = $region_action->get_region; #Bio::Vega::Region
#     print "\nFOUND ".scalar($region->genes)." GENES IN REGION ".$gene->seq_region_name.":".$gene->start."\n";
# 
#     #Reset gene coordinates to the start of the region slice,
#     #otherwise the gene coordinates in loutre will be wrong
#     my $slice_offset = $region->slice->start - 1;
#     print "OFFSET=$slice_offset\n";
#     foreach my $tr (@{$gene->get_all_Transcripts}){
#         foreach my $exon (@{$tr->get_all_Exons}){
#             $exon->start($exon->start - $slice_offset);
#             $exon->end(  $exon->end   - $slice_offset);
#             $exon->slice($region->slice);
#         }
#     }
# 
# 
# 
# 
#                    my $merged_transcript = merge_transcripts($tr, $sel_db_tr );
# #     $db_gene->flush_Transcripts;
# #                     #$db_gene->remove_Transcript($sel_db_tr);
# #     $db_gene->add_Transcript($merged_transcript);    
#                     #$slice_offset = $region->slice->start - 1;
# #     foreach my $tr2 (@{$db_gene->get_all_Transcripts}){
#                     foreach my $exon (@{$merged_transcript->get_all_Exons}){
#                         $exon->start($exon->start - $slice_offset);
#                         $exon->end(  $exon->end   - $slice_offset);
#                         $exon->slice($region->slice);
#                     }
# #     }
#                     foreach my $g ($region->genes) {
#                         foreach my $ts (@{$g->get_all_Transcripts}) {
#                             if ($ts->stable_id eq $merged_transcript->stable_id){
#                                 $ts->flush_Exons; #This empties the exon list. If the same exons are added again, 
#                                                   #those that are not associated with other transcripts will get a new stable id
#                                 foreach my $exon (@{$merged_transcript->get_all_Exons}){
#                                     $ts->add_Exon($exon);
#                                 }
#                 
#                                 #Add remarks (except 'not for VEGA') and hidden remarks from the novel transcript
#                                 foreach my $att (@{$merged_transcript->get_all_Attributes}){
#                                     #COMMENTING THIS OUT AS WE MAY WANT TO KEEP THE NOT FOR VEGA REMARK (NEEDS TO READ THE $no_NFV VARIABLE) - TO DO
#                                     #next if ($att->code eq "remark" and $att->value eq "not for VEGA");
#                                     #Append read names to existing remark
#                                     if ((($att->code eq "hidden_remark" and $att->value =~ /^pacbio_capture_seq_\w+ : .+;/) or
#                                         ($att->code eq "hidden_remark" and $att->value =~ /^SLR-seq_\w+ : .+;/) or
#                                         ($att->code eq "hidden_remark" and $att->value =~ /^pacbio_raceseq_\w+ : .+;/)) and
#                                         (scalar(@{$ts->get_all_Attributes('TAGENE_transcript')}) > 0)){
#                                         my $appended = 0;
#                                         foreach my $att2 (@{$ts->get_all_Attributes('hidden_remark')}){
#                                             if ($att2->value =~ /^pacbio_capture_seq_\w+ : .+;/ or
#                                                 $att2->value =~ /^SLR-seq_\w+ : .+;/ or
#                                                 $att2->value =~ /^pacbio_raceseq_\w+ : .+;/){
#                                                 $att2->value($att2->value." ".$att->value);
#                                                 $appended = 1;
#                                                 last;
#                                             }
#                                         }
#                                         #If no read name remark yet, create one
#                                         unless ($appended){
#                                             $ts->add_Attributes($att);
#                                         }
#                                     }
#                                     else{
#                                         $ts->add_Attributes($att);
#                                     }
#                                 }
#                                 #If 'comp_pipe' transcript, change biotype and remove 'not for VEGA' attribute
#                                 #unless 'comp_pipe_rejected' remark
#                                 if ($ts->biotype eq "comp_pipe" and scalar(grep {$_->value eq "comp_pipe_rejected"} @{$ts->get_all_Attributes('remark')}) == 0){
#                                     my $t_biotype = get_transcript_biotype($db_gene);
#                                     $ts->biotype($t_biotype);
#                                     if (scalar(grep {$_->value eq "not for VEGA"} @{$ts->get_all_Attributes('remark')})){
#                                        my $array = $ts->{'attributes'};
#                                        @$array = grep { $_->value ne "not for VEGA" } @$array;
#                                     }
#                                 }
# 
#                                 #If coding gene, try to assign a CDS and change the biotype accordingly
#                                 assign_cds_to_transcripts($ts, $g, $slice_offset);
#                                 
#                                 #Change authorship
#                                 $ts->transcript_author($merged_transcript->transcript_author);
#                                 
#                                 print "TR: $id: Will extend transcript ".$ts->stable_id." (".$ts->biotype.") in gene ".$g->stable_id."\n";
#                                 push (@log, "TR2: $id: Extended transcript ".$ts->stable_id." (".$ts->biotype.") in gene ".$g->stable_id." (".$g->biotype.")");
#                             }
#                         }
#                     }
# 
