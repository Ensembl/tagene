
package LoutreWrite::GeneFilter;

use strict;
use warnings;
#use Sys::Hostname;
#use Try::Tiny;
use List::Util qw[min max];
use List::MoreUtils qw(uniq);
#use LoutreWrite::Config;




=head2 check_artifact_transcripts

 Arg[1]    : list of Bio::Vega::Gene objects
 Function  : Check for transcripts unusually long and spanning several loci - probably alignment artifacts (NOT THOROUGHLY CHECKED YET)
 Returntype: Hashref - list of transcript names

=cut

sub check_artifact_transcripts {
  my ($self, $genes) = @_;
  my %list;
  my $min_num_genes = 4;

  foreach my $gene (@$genes){
    foreach my $transcript (@{$gene->get_all_Transcripts}){
      #Fetch database genes overlapping our transcript
#      my $tr_slice = $transcript->slice->adaptor->fetch_by_region("toplevel", $transcript->seq_region_name, $transcript->start, $transcript->end);
#      my @db_genes = grep {$_->seq_region_strand == $gene->seq_region_strand} @{$tr_slice->get_all_Genes};
      my @db_genes = @{get_valid_overlapping_genes($transcript)};
      if (scalar @db_genes > $min_num_genes){
        my $t_name = $transcript->stable_id || $transcript->get_all_Attributes('hidden_remark')->[0]->value;
        print "KILL: ".$t_name."  ".$transcript->start."-".$transcript->end."\n";
        $list{$t_name} = 1;
      }
    }
  }

  return \%list;
}



=head2 check_max_overlapped_loci

 Arg[1]    : list of Bio::Vega::Gene objects
 Arg[2]    : Integer
 Arg[3]    : 1/0 (stringent check, just simple exonic overlap, no clusters taken into account)
 Function  : Check for transcripts overlapping a number of existing genes at the exon level that exceeds the number given in the second argument.
             Only genes of certain biotypes (coding and lncRNA but not small RNAs or pseudogenes), and not tagged as readthrough, count towards that number.
 Returntype: Hashref - list of transcript names

=cut

sub check_max_overlapped_loci {
  my ($self, $genes, $max_ov_loci, $stringent_check) = @_;
  my %list;
  foreach my $gene (@$genes){
    foreach my $transcript (@{$gene->get_all_Transcripts}){
      my $message = "";
      #Fetch database genes overlapping our transcript
      my @db_genes = @{get_valid_overlapping_genes($transcript)};
      my @overlapped_db_genes; #genes to be considered for overlap count: if gene has MANE Select, exonic overlap with MANE Select; else, if non-readthrough gene, exonic overlap with any transcript
      my @overlapped_select_db_genes; #genes with MANE Select transcript that has exonic overlap with the input transcript
      my @overlapped_readthrough_db_genes; #readthrough genes with exonic overlap with the input transcript
      DBG:foreach my $db_gene (@db_genes){
        #Count only certain biotypes (coding and lncRNA broad biotypes)
        if ($db_gene->biotype =~ /antisense|bidirectional_promoter_lncrna|ig_gene|IG_V_gene|lincRNA|macro_lncrna|overlapping_ncrna|polymorphic_pseudogene|processed_transcript|protein_coding|sense_intronic|sense_overlapping|tec|tr_gene/i){
          my $mane_select;
          my $is_readthrough;
          foreach my $db_tr (@{$db_gene->get_all_Transcripts}){
            next if $db_tr->biotype eq "artifact";
            next if scalar grep {$_->value eq "not for VEGA"} @{$db_tr->get_all_Attributes('remark')};
            #Do not count readthrough genes as this may affect transcripts belonging to the overlapping non-readthrough genes
            #In any case, if the novel transcript is readthrough, it will be assigned to the readthrough gene later and subsequently rejected
            if (scalar grep {$_->value eq "readthrough"} @{$db_tr->get_all_Attributes("remark")}){
              $is_readthrough = 1;
              #next DBG;
            }
            #Store the MANE Select transcript, to be used later
            if (scalar grep {$_->value eq "MANE_select"} @{$db_tr->get_all_Attributes("remark")}){
              $mane_select = $db_tr;
            }
          }
          #Find if there is exonic overlap with the MANE Select transcript
          #If a gene has a MANE Select transcript, count it as an overlapped gene only if there is exonic overlap with that transcript
          #NOTE: this could be extended to all canonical transcripts if this info was available in the havana database (or obtained from the core database)
          if ($mane_select){
            E:foreach my $db_exon (@{$mane_select->get_all_Exons}){
              foreach my $exon (@{$transcript->get_all_Exons}){
                if ($db_exon->seq_region_start <= $exon->seq_region_end and $db_exon->seq_region_end >= $exon->seq_region_start){
                  push(@overlapped_select_db_genes, $db_gene);
                  push(@overlapped_db_genes, $db_gene);
                  $message .= "OVms: ".$db_gene->stable_id." ".$db_gene->biotype."  ";
                  last E;
                }
              }
            }
          }
#11-6-2025
#Changes to avoid adding readthrough transcripts:
# -either avoid exonic overlap with any transcript,
#  or at the very least avoid sharing splice site with non-MANE Select transcript
# -if pre-existing overlapping genes create false readthrough cases,
#  check if the transcript overlaps two non-overlapping genes at the exon level

#Possible algorithm:
#If any two exonic-overlapped genes don't overlap at the exon level themselves, then reject.
#If they already overlap but don't share splice sites, and the transcript shares splice sites with both, reject.
#In case of overlap with a readthrough gene and a non-readthrough gene, check that it doesn't overlap too much with the readthrough gene... (that's the 3kb limit below)
#Or check if the novel transcript shares a splice site with the readthrough gene that wasn't shared between non-readthrough and readthrough genes

#Stringent mode: strictly no exonic overlap with two genes 

#Difficult to make this work if the number of max ovlpd loci is not 1...

# -require CAGE support for 5' extensions?
# -for the 3' end, currently sharing MANE Select stop codon removes requirement for polyA support, which allows 3'end extensions of any length...



          #If there is no MANE Select, find if there is exonic overlap with any transcript
          else{
            E:foreach my $db_exon (@{$db_gene->get_all_Exons}){
              foreach my $exon (@{$transcript->get_all_Exons}){
                if ($db_exon->seq_region_start <= $exon->seq_region_end and $db_exon->seq_region_end >= $exon->seq_region_start){
                  if ($is_readthrough){
                    push(@overlapped_readthrough_db_genes, $db_gene);
                  }
                  else{
                    push(@overlapped_db_genes, $db_gene);
                  }
                  $message .= "OV: ".$db_gene->stable_id." ".$db_gene->biotype."  ";
                  last E;
                }
              }
            }
          }
        }
      }
   
      if ($stringent_check){
        my $t_name = $transcript->stable_id || $transcript->get_all_Attributes('hidden_remark')->[0]->value;
        if (scalar @overlapped_db_genes > $max_ov_loci){
          print "KILL_2: ".$t_name."  ".$transcript->start."-".$transcript->end." overlaps ".scalar(@overlapped_db_genes)." loci: $message\n";
          $list{$t_name} = 1;
        }
        #If exonic overlap with readthrough gene, even if there is exonic overlap with only other transcript,
        #calculate the ratio between the transcript span and both gene spans, or set a limit on how far upstream or downstream
        #of the host gene the novel transcript can go. 
        #This will be useful in case the readthrough gene overlaps only that gene (eg. STAG3L5P-PVRIG2P-PILRB and PILRB)
        #Note: implemented below for the most common use case: max. one overlapping locus allowed, but it wouldn't work for other cases (TO DO)
        else{
          if ($max_ov_loci == 1 and scalar @overlapped_db_genes == 1 and scalar @overlapped_readthrough_db_genes > 0){
            if ($transcript->seq_region_start < ($overlapped_db_genes[0]->seq_region_start - 3000) or
                $transcript->seq_region_end > ($overlapped_db_genes[0]->seq_region_end + 3000)
            ){
              print "KILL_2: ".$t_name."  ".$transcript->start."-".$transcript->end." overlaps ".scalar(@overlapped_db_genes)." loci and ".scalar(@overlapped_readthrough_db_genes)." readthroughs: $message\n";
              $list{$t_name} = 1;
            }
          }
        }
      }

  #NOTE: we could stop here, or if the limit is exceeded we could check if the overlapped genes
  #      are part of a single cluster in which the canonical transcripts overlap...
  #If 'stringent_check' is not set, the code below will determine which novel transcripts are killed.
  #If it is set, the code below will run anyway but it should be less stringent than the code above,
  #so it won't change the outcome.
          
      #Instead of counting overlapped genes, count clusters of overlapped genes.
      #This is fine for lncRNAs. For coding genes, the presence of readthrough transcripts that are not tagged as such
      #can lead to two coding genes being included in the same cluster and to the insertion of novel readthrough transcripts, 
      #which we don't want.
      #REFINE THE DEFINITION OF CLUSTERS: ONE GENE INCLUDED IN THE OTHER, OR SUBSTANTIAL OVERLAP, OR SHARED SPLICE SITES!
      my @clusters;
      foreach my $db_gene (sort {$a->seq_region_start <=> $b->seq_region_start} @overlapped_db_genes){
        my $clustered = 0;
        if (scalar @clusters){
          foreach my $cluster (@clusters){
            if ($cluster->{'start'} <= $db_gene->seq_region_end and $cluster->{'end'} >= $db_gene->seq_region_start){
              if ($cluster->{'start'} > $db_gene->seq_region_start){
                $cluster->{'start'} = $db_gene->seq_region_start;
              }
              if ($cluster->{'end'} < $db_gene->seq_region_end){
                $cluster->{'end'} = $db_gene->seq_region_end;
              }
              push(@{$cluster->{'genes'}}, $db_gene->stable_id);
              $clustered = 1;
              last;
            }
          }
        }
        if (scalar @clusters == 0 or !$clustered){
          my %cluster;
          $cluster{'start'} = $db_gene->seq_region_start;
          $cluster{'end'} = $db_gene->seq_region_end;
          push(@{$cluster{'genes'}}, $db_gene->stable_id);
          push(@clusters, \%cluster);
        }
      }
            
      if (scalar @clusters > $max_ov_loci){
        my $t_name = $transcript->stable_id || $transcript->get_all_Attributes('hidden_remark')->[0]->value;
        print "KILL_2a: ".$t_name."  ".$transcript->start."-".$transcript->end." overlaps ".scalar(@clusters)." loci: $message\n";
        $list{$t_name} = 1;
      }

###
      #In addition, check if the transcript shares splice sites with two genes
      #and these splice sites are beyond the other gene's boundaries
      #If one of them is a readthrough gene, it should be ignored (filtered previously)
      #but do not ignore genes where the overlap is with non-MANE_Select transcripts!
      my %shared_splice_sites;
      foreach my $intron (@{$transcript->get_all_Introns}){
#        foreach my $db_gene (@overlapped_db_genes){
foreach my $db_gene (@db_genes){          
print "OVLP=".$db_gene->stable_id."\n";
          foreach my $db_intron (@{$db_gene->get_all_Introns}){
            if ($intron->seq_region_start == $db_intron->seq_region_start){
              $shared_splice_sites{$intron->seq_region_start}{$db_gene->stable_id} = 1;
            }
            if ($intron->seq_region_end == $db_intron->seq_region_end){
              $shared_splice_sites{$intron->seq_region_end}{$db_gene->stable_id} = 1;
            }
          }
        }
      }
      my %genes_not_overlapping_splice_site;
      foreach my $ss_site (keys %shared_splice_sites){
print "SS=".$ss_site." ".join(",", keys %{$shared_splice_sites{$ss_site}})."\n";
        if (scalar keys %{$shared_splice_sites{$ss_site}} == 1){
#          foreach my $db_gene (grep {$_->stable_id ne $shared_splice_sites{$ss_site}} @overlapped_db_genes){
foreach my $db_gene (grep {$_->stable_id ne $shared_splice_sites{$ss_site}} @db_genes){            
            unless ($db_gene->seq_region_start < $ss_site and $db_gene->seq_region_end > $ss_site){
              $genes_not_overlapping_splice_site{$db_gene->stable_id} = 1;
            }
          }
        }
      }
      #if (scalar @overlapped_db_genes and scalar keys %genes_not_overlapping_splice_site == scalar @overlapped_db_genes){
      #  if (scalar @clusters > $max_ov_loci){
      if (scalar @overlapped_db_genes and scalar keys %genes_not_overlapping_splice_site > $max_ov_loci){ #TEST THIS
          my $t_name = $transcript->stable_id || $transcript->get_all_Attributes('hidden_remark')->[0]->value;
          print "KILL_2b: ".$t_name."  ".$transcript->start."-".$transcript->end." overlaps ".scalar(@overlapped_db_genes)." loci: $message; ".
                " it shares non-overlapping splice sites with ".join(", ", keys %genes_not_overlapping_splice_site)."\n";
          $list{$t_name} = 1;
      #  }
      }
###      
       
    }
  }

  return \%list;
}



=head2 check_multiple_gene_overlap

 Arg[1]    : list of Bio::Vega::Gene objects
 Arg[2]    : Boolean (true for stringent check, no exonic overlap with more than one gene)
 Function  : Check for transcripts overlapping multiple genes at the exon level to avoid creating readthrough transcripts.
             Only genes of certain biotypes (coding and lncRNA but not small RNAs or pseudogenes) are counted.
             A transcript can overlap two genes if these genes already overlap each other.
             A transcript can share splice sites with two genes if these genes share the same splice sites with each other.
 Returntype: Hashref - list of transcript names to be rejected

=cut

sub check_multiple_gene_overlap {
  my ($self, $genes, $stringent_check) = @_;
  my %list;
  foreach my $gene (@$genes){
    foreach my $transcript (@{$gene->get_all_Transcripts}){
      my $t_name = $transcript->stable_id || $transcript->get_all_Attributes('hidden_remark')->[0]->value;
      my $message = "";
      #Fetch database genes overlapping the input transcript
      my @db_genes = @{get_valid_overlapping_genes($transcript)};
      #Find which genes overlap with the input transcript at the exon level
      my @overlapped_db_genes;
      my %overlapped_readthrough_db_genes; #ids of readthrough genes overlapped at the exon level
      DBG:foreach my $db_gene (@db_genes){
        #Count only certain biotypes (coding and lncRNA broad biotypes)
        if ($db_gene->biotype =~ /antisense|bidirectional_promoter_lncrna|ig_gene|IG_V_gene|lincRNA|macro_lncrna|overlapping_ncrna|polymorphic_pseudogene|processed_transcript|protein_coding|sense_intronic|sense_overlapping|tec|tr_gene/i){
          foreach my $db_tr (@{$db_gene->get_all_Transcripts}){
            next if $db_tr->biotype eq "artifact";
            next if scalar grep {$_->value eq "not for VEGA"} @{$db_tr->get_all_Attributes('remark')};
            my $is_readthrough = 0;
            if (scalar grep {$_->value eq "readthrough"} @{$db_tr->get_all_Attributes('remark')}){
              $is_readthrough = 1;
            }
            E:foreach my $db_exon (@{$db_tr->get_all_Exons}){
              foreach my $exon (@{$transcript->get_all_Exons}){
                if ($db_exon->seq_region_start <= $exon->seq_region_end and $db_exon->seq_region_end >= $exon->seq_region_start){
                  push(@overlapped_db_genes, $db_gene);
                  $message .= "OV: ".$db_gene->stable_id." ".$db_gene->biotype."  "; print $message."\n";
                  if ($is_readthrough){
                    $overlapped_readthrough_db_genes{$db_gene->stable_id} = 1;
                  } 
                  next DBG;
                }
              }
            }
          }
        }
      }
      
      #If more than one overlapped gene, find if the overlapping genes also overlap with each other at the exon level
      if (scalar @overlapped_db_genes > 1){
        my %non_exonic_overlapping_genes;
        for (my $i = 0; $i < $#overlapped_db_genes; $i++){
          for (my $j = $i+1; $j <= $#overlapped_db_genes; $j++){
            my $found_overlap;
            TR:foreach my $db_tr1 (@{$overlapped_db_genes[$i]->get_all_Transcripts}){
              next if $db_tr1->biotype eq "artifact";
              next if scalar grep {$_->value eq "not for VEGA"} @{$db_tr1->get_all_Attributes('remark')};
              foreach my $db_tr2 (@{$overlapped_db_genes[$j]->get_all_Transcripts}){
                next if $db_tr2->biotype eq "artifact";
                next if scalar grep {$_->value eq "not for VEGA"} @{$db_tr2->get_all_Attributes('remark')};
                foreach my $db_exon1 (@{$db_tr1->get_all_Exons}){
                  foreach my $db_exon2 (@{$db_tr2->get_all_Exons}){
                    if ($db_exon1->seq_region_start <= $db_exon2->seq_region_end and $db_exon1->seq_region_end >= $db_exon2->seq_region_start){
                      $found_overlap = 1;
                      last TR;
                    }
                  }
                }
              }
            }
            unless ($found_overlap){
              $non_exonic_overlapping_genes{$overlapped_db_genes[$i]->stable_id."-".$overlapped_db_genes[$j]->stable_id} = 1;
            }
          }
        }
        #Skip transcript if it shares exonic sequence with two non-overlapping genes
        #This will also skip a transcript if the two genes had intronic overlap with each other
        if (scalar keys %non_exonic_overlapping_genes){
          print "KILL_2a: ".$t_name."  ".$transcript->start."-".$transcript->end." overlaps ".scalar(@overlapped_db_genes)." loci: $message; ".
                "it overlaps genes that do not overlap with each other at the exon level: ".join(", ", keys %non_exonic_overlapping_genes)."\n";
          $list{$t_name} = 1;
        }

        #Find if the transcript shares splice sites with the overlapping genes.
        #If it shares different splice sites with at least two genes, check if these genes already share splice sites with each other,
        #which should be alright if neither is a readthrough gene.
        #We can also check if the shared splice sites are located beyond the boundaries of other overlapping genes, although this can create false positives,
        #for instance when a new 5'-terminal exon is created and one the new splice sites is shared with an overlapping locus
        #Also, reject the transcript if it shares a splice site that is unique to a readthrough gene
        else{
          #Find shared splice sites
          my %shared_splice_sites; #splice sites of novel transcript shared with current genes
          my %splice_site_sharing_genes; #ids of current genes sharing a splice site with the novel transcript
          foreach my $intron (@{$transcript->get_all_Introns}){
            foreach my $db_gene (@overlapped_db_genes){
              foreach my $db_tr (@{$db_gene->get_all_Transcripts}){
                next if $db_tr->biotype eq "artifact";
                next if scalar grep {$_->value eq "not for VEGA"} @{$db_tr->get_all_Attributes('remark')};
                foreach my $db_intron (@{$db_tr->get_all_Introns}){
                  if ($intron->seq_region_start == $db_intron->seq_region_start){
                    $shared_splice_sites{$intron->seq_region_start}{$db_gene->stable_id} = 1;
                    $splice_site_sharing_genes{$db_gene->stable_id} = 1;
                  }
                  if ($intron->seq_region_end == $db_intron->seq_region_end){
                    $shared_splice_sites{$intron->seq_region_end}{$db_gene->stable_id} = 1;
                    $splice_site_sharing_genes{$db_gene->stable_id} = 1;
                  }
                }
              }
            }
          }

          if (scalar keys %splice_site_sharing_genes > 1){
            #Get gene objects of the genes sharing a splice site with the novel transcript
            my @splice_site_sharing_genes;
            foreach my $ov_gene (@overlapped_db_genes){
              if ($splice_site_sharing_genes{$ov_gene->stable_id}){
                push(@splice_site_sharing_genes, $ov_gene);
              }
            }
            #Find genes that overlap the novel transcript but do not share any splice site with each other
            my %non_splice_site_sharing_overlapping_genes; 
            for (my $i = 0; $i < $#splice_site_sharing_genes; $i++){
              for (my $j = $i+1; $j <= $#splice_site_sharing_genes; $j++){
                my $found_shared_splice_sites;
                TR:foreach my $db_tr1 (@{$splice_site_sharing_genes[$i]->get_all_Transcripts}){
                  next if $db_tr1->biotype eq "artifact";
                  next if scalar grep {$_->value eq "not for VEGA"} @{$db_tr1->get_all_Attributes('remark')};
                  foreach my $db_tr2 (@{$splice_site_sharing_genes[$j]->get_all_Transcripts}){
                    next if $db_tr2->biotype eq "artifact";
                    next if scalar grep {$_->value eq "not for VEGA"} @{$db_tr2->get_all_Attributes('remark')};
                    foreach my $db_intron1 (@{$db_tr1->get_all_Introns}){
                      foreach my $db_intron2 (@{$db_tr2->get_all_Introns}){
                        if ($db_intron1->seq_region_start == $db_intron2->seq_region_start or $db_intron1->seq_region_end == $db_intron2->seq_region_end){
                          $found_shared_splice_sites = 1;
                          last TR;
                        }
                      }
                    }
                  }
                }
                unless ($found_shared_splice_sites){
                  $non_splice_site_sharing_overlapping_genes{$splice_site_sharing_genes[$i]->stable_id."-".$splice_site_sharing_genes[$j]->stable_id} = 1;
                }
              }
            }
            #Reject the novel transcript if it shares splice sites with different genes that do not share splice sites with each other
            if (scalar keys %non_splice_site_sharing_overlapping_genes){
              print "KILL_2b: ".$t_name."  ".$transcript->start."-".$transcript->end." overlaps ".scalar(@overlapped_db_genes)." loci: $message; ".
                "it shares splice sites with genes that do not share splice sites with each other: ".join(", ", keys %non_splice_site_sharing_overlapping_genes)."\n";
              $list{$t_name} = 1;
            }
          }
        

          #Find splice sites shared with readthrough genes only
          foreach my $splice_site (keys %shared_splice_sites){
            if (scalar keys %{$shared_splice_sites{$splice_site}} == 1){
              foreach my $rtid (keys %overlapped_readthrough_db_genes){
                if ($shared_splice_sites{$splice_site}{$rtid}){
                  print "KILL_2d: ".$t_name."  ".$transcript->start."-".$transcript->end." overlaps ".scalar(@overlapped_db_genes)." loci: $message; ".
                        " it shares splice sites uniquely with a readthrough gene: $rtid\n";
                  $list{$t_name} = 1;
                }
              }
            }
          }

          #Find shared splice sites that are beyond the boundaries of other overlapped genes
          #Add a buffer area to avoid false positives caused by eg. novel 5' exons?
          #THIS CHECK IS ALMOST SUPERFLUOUS AFTER 2d - SOME CASES NOT PICKED UP BY 2d TOO LOOK WRONG - REMOVE IT?
          my %genes_not_overlapping_splice_site;
          foreach my $splice_site (keys %shared_splice_sites){
            print "SS=".$splice_site." ".join(",", keys %{$shared_splice_sites{$splice_site}})."\n";
            if (scalar keys %{$shared_splice_sites{$splice_site}} < scalar @overlapped_db_genes){
              foreach my $db_gene (@overlapped_db_genes){
                unless ($shared_splice_sites{$splice_site}{$db_gene->stable_id}){
                  unless ($db_gene->seq_region_start < $splice_site and $db_gene->seq_region_end > $splice_site){
                    $genes_not_overlapping_splice_site{$db_gene->stable_id} = 1;
                  }
                }
              }
            }
          }
          if (scalar keys %genes_not_overlapping_splice_site){
#            print "KILL_2c: ".$t_name."  ".$transcript->start."-".$transcript->end." overlaps ".scalar(@overlapped_db_genes)." loci: $message; ".
#                  " it shares splice sites that are beyond the boundaries of other overlapped genes: ".join(", ", keys %genes_not_overlapping_splice_site)."\n";
#            $list{$t_name} = 1;
          }


        }
      }
    }
  }

  return \%list;
}



=head2 check_overlapped_gene_biotypes

 Arg[1]    : list of Bio::Vega::Gene objects
 Arg[2]    : Hashref - list of unallowed biotypes
 Function  : Check for transcripts overlapping at the exon level a gene with an unallowed biotype that belongs to the list given in the second argument.
 Returntype: Hashref - list of transcript names

=cut

sub check_overlapped_gene_biotypes {
  my ($self, $genes, $biotypes) = @_;
  my %list;
  foreach my $gene (@$genes){
    foreach my $transcript (@{$gene->get_all_Transcripts}){
      my $message = "";
      #Fetch database genes overlapping our transcript
      my @db_genes = @{get_valid_overlapping_genes($transcript)};
      my @overlapped_db_genes;
      DBG:foreach my $db_gene (@db_genes){
        foreach my $db_exon (@{$db_gene->get_all_Exons}){
          foreach my $exon (@{$transcript->get_all_Exons}){
            if ($db_exon->seq_region_start <= $exon->seq_region_end and $db_exon->seq_region_end >= $exon->seq_region_start){
              if ($biotypes->{$db_gene->biotype}){
                my $t_name = $transcript->stable_id || $transcript->get_all_Attributes('hidden_remark')->[0]->value;
                print "KILL_2: ".$t_name."  ".$transcript->start."-".$transcript->end." overlaps gene ".$db_gene->stable_id." - ".$db_gene->biotype."\n";
                $list{$t_name} = 1;
                next DBG;
              }
            }
          }
        }
      }
    }
  }

  return \%list;
}



=head2 recluster_transcripts

 Arg[1]    : list of Bio::Vega::Gene objects
 Arg[2]    : Hashref - list of transcript names to be ignored
 Function  : Check if there are gaps between transcripts and split genes if appropriate
 Returntype: list of Bio::Vega::Gene objects

=cut

sub recluster_transcripts {
  my ($self, $genes, $kill_list) = @_;
  my @checked_genes;
    
  foreach my $gene (@$genes){
    my @clusters;    
    foreach my $transcript (sort {$a->start <=> $b->start} @{$gene->get_all_Transcripts}){
      my $t_name = $transcript->stable_id || $transcript->get_all_Attributes('hidden_remark')->[0]->value;
      next if $kill_list->{$t_name}; #Ignore transcript if in kill list
      #print "T:".$transcript->start."-".$transcript->end."  ".$transcript->stable_id."\n";
      my $clustered = 0;
      if (scalar @clusters){
        foreach my $cluster (@clusters){
          if ($cluster->{start} <= $transcript->end and $cluster->{end} >= $transcript->start){
            if ($cluster->{start} > $transcript->start){
              $cluster->{start} = $transcript->start;
            }
            if ($cluster->{end} < $transcript->end){
              $cluster->{end} = $transcript->end;
            }
            push(@{$cluster->{transcripts}}, $transcript);
            $clustered = 1;
            last;
          }
        }
      }
      if (scalar @clusters == 0 or !$clustered){
        my %cluster;
        $cluster{start} = $transcript->start;
        $cluster{end} = $transcript->end;
        push(@{$cluster{transcripts}}, $transcript);
        push(@clusters, \%cluster);
      }
    }
    print "Clusters: ".scalar(@clusters)."\n";
    #if (scalar @clusters == 1){
    #    push (@checked_genes, $gene);
    #}
    #else{
    foreach my $cluster (@clusters){
      my $new_gene = Bio::Vega::Gene->new(
                                    -stable_id  => $gene->stable_id,
                                    -biotype    => $gene->biotype,
                                    -analysis   => $gene->analysis,
                                    -source     => $gene->source,
                                    -status     => $gene->status
                                    );
      $new_gene->gene_author($gene->gene_author);
      $new_gene->add_Attributes(@{$gene->get_all_Attributes});
        
      foreach my $transcript (@{$cluster->{transcripts}}){ 
        #Check again that no transcript on the kill list goes through
        # as the first checkpoint seems to fail (?)
        my $t_name = $transcript->stable_id || $transcript->get_all_Attributes('hidden_remark')->[0]->value;
        unless ($kill_list->{$t_name}){
          $new_gene->add_Transcript($transcript);
        }
      }
      unless (scalar(@{$new_gene->get_all_Transcripts}) == 0){
        push (@checked_genes, $new_gene);
      }
    }
    #}
  }
    
  return \@checked_genes;
}



=head2 find_host_gene

 Arg[1]    : Bio::Vega::Gene object
 Function  : Find and return the most suitable existing gene (if any) to host the transcripts of the input gene
 Returntype: Bio::Vega::Gene object

=cut

sub find_host_gene {
  my ($self, $gene) = @_;
  my $host_gene;
    
  #Fetch database genes overlapping our gene
#   my $gene_slice = $gene->slice->adaptor->fetch_by_region("toplevel", $gene->seq_region_name, $gene->start, $gene->end);
#   my @db_genes = grep {$_->stable_id ne $gene->stable_id and $_->seq_region_strand == $gene->strand} @{$gene_slice->get_all_Genes};
  my @db_genes = @{get_valid_overlapping_genes($gene)};
  #Find host gene with the most matching introns with the input gene, otherwise the gene with the most exon overlap
  #Make sure that all transcripts of the input gene overlap the host gene or none (a readthrough transcript may bias 
  #the choice towards a host gene that is not overlapped by the remaining transcripts)
        
  #Check intron match - find gene with the most matching introns
  my @candidates;
  my @best_candidates;
  my $max_n_introns = 0;
  my %compatible_hosts;
  foreach my $db_gene (@db_genes){
    my $n_introns = 0; #Number of matching introns
    my %matching_introns;
    foreach my $db_intron (@{$db_gene->get_all_Introns}){
      foreach my $tr (@{$gene->get_all_Transcripts}){
        foreach my $intron (@{$tr->get_all_Introns}){
          if ($db_intron->seq_region_start == $intron->seq_region_start and $db_intron->seq_region_end == $intron->seq_region_end){
            $matching_introns{$db_intron->seq_region_start."-".$db_intron->seq_region_end} = 1;
            $compatible_hosts{$tr->get_all_Attributes('hidden_remark')->[0]->value}{$db_gene->stable_id} = 1; #Record if input transcripts are compatible with this host gene candidate (at least one shared intron). As this is mostly relevant for readthrough transcripts, it will be applied to the intron match analysis only
          }
        }
      }            
    }
    $n_introns = scalar keys %matching_introns;
    if ($n_introns > 0){ #Keep record of good candidates
      push(@candidates, $db_gene);
    }
    if ($n_introns > $max_n_introns){ #If better candidate, reset list
      @best_candidates = ();
      push(@best_candidates, $db_gene);
      $max_n_introns = $n_introns;
    }
    elsif ($n_introns == $max_n_introns){
      push(@best_candidates, $db_gene);        
    }        
  }
  #If only one...
  if (scalar(@best_candidates) == 1){
    $host_gene = $best_candidates[0];
  }
  #Else, check exon overlap - find host gene with the most exon overlap
  else{
    if (scalar(@best_candidates) == 0){
      @best_candidates = @db_genes;
    }
    my %gene_loc; #Store genomic positions of gene's exons
    foreach my $exon (@{$gene->get_all_Exons}){
      for (my $i = $exon->seq_region_start; $i <= $exon->seq_region_end; $i++){
        $gene_loc{$i} = 1;
      }
    }
    my $max = 0;
    foreach my $db_gene (@best_candidates){
      my %seen; #Store genomic positions where gene's exons overlaps host gene's exons
      foreach my $db_exon (@{$db_gene->get_all_Exons}){
        for (my $i = $db_exon->seq_region_start; $i <= $db_exon->seq_region_end; $i++){
          if ($gene_loc{$i}){
            $seen{$i} = 1;
          }
        }
      }
      if (scalar(keys %seen) > $max){
        $host_gene = $db_gene;
        $max = scalar(keys %seen);
      }
    }
  }
  if ($host_gene and scalar(@db_genes)>1){
    print "Selected ".$host_gene->stable_id." from ".join(", ", map {$_->stable_id} @db_genes)."\n";
    foreach my $tr (@{$gene->get_all_Transcripts}){
      my $tid = $tr->get_all_Attributes('hidden_remark')->[0]->value;
      #if ($compatible_hosts{$tid} and !($compatible_hosts{$tid}{$host_gene->stable_id})){
      print "\tHost gene not compatible with transcript ".$tid."\n";
      #}
    }
  }

  return $host_gene;
}



=head2 assign_host_gene

 Arg[1]    : Bio::Vega::Gene object
 Arg[2]    : preserve host gene biotype, eg. ignore non-coding host gene if transcript is coding
 Function  : Find the most suitable existing gene(s) to host the transcripts of the given gene, 
             split the input gene into as many genes as host genes were found and set the new gene stable id(s);
             return the transformed input gene
 Returntype: list of Bio::Vega::Gene objects

=cut

sub assign_host_gene {
  my ($self, $gene, $preserve_biotype) = @_;
#print "AAAAA:".$gene->seq_region_start."\n";
  #Fetch database genes overlapping our gene
#  my $gene_slice = $gene->slice->adaptor->fetch_by_region("toplevel", $gene->seq_region_name, $gene->start, $gene->end);
#  my @db_genes = grep {$_->stable_id ne $gene->stable_id and $_->seq_region_strand == $gene->strand} @{$gene_slice->get_all_Genes};
  my @db_genes = @{get_valid_overlapping_genes($gene, 1)};
  my %host_by_transcript;
  foreach my $transcript (@{$gene->get_all_Transcripts}){
#print "BBBBB:".$transcript->seq_region_start."\n";
    my $host_gene;
    #Check intron match - find gene with the most matching introns
    my @candidates;
    my $max_n_introns = 0;
    foreach my $db_gene (@db_genes){
#            #Exclude undesired host genes
#            if (!($db_gene->source =~ /(ensembl|havana)/) or
#                $db_gene->biotype eq "artifact" or
#                $db_gene->biotype eq "comp_pipe" or
#                scalar (grep {$_->value eq "not for VEGA"} @{$db_gene->get_all_Attributes('remark')}) == 1
#            ){
#              next;
#            }
      #EXPERIMENTAL: exclude pseudogenes so that spliced transcripts are assigned to new (lncRNA) genes
      if (scalar @{$transcript->get_all_Introns} > 1 and $db_gene->biotype =~ /pseudogene/){
        next;
      }
            
      if ($preserve_biotype){
        if ($transcript->translation and $db_gene->biotype ne "protein_coding"){
          next;
        }
      }
      my $n_introns = 0; #Number of matching introns
      foreach my $db_intron (sort {$a->start <=> $b->start} @{$db_gene->get_all_Introns}){
#print $db_gene->stable_id.".".$db_gene->version.": ".$db_intron->seq_region_start."-".$db_intron->seq_region_end."\n";            
        foreach my $intron (@{$transcript->get_all_Introns}){
			#print "CCCCC: ".$intron->start." vs ".$db_intron->seq_region_start."\n";
          if ($db_intron->seq_region_start == $intron->start and $db_intron->seq_region_end == $intron->end){
            $n_introns++; #print "NNN:".$n_introns."\n"; #print $db_gene->stable_id.": ".$db_intron->seq_region_start."-".$db_intron->seq_region_end."\n" if $transcript->start == 38024151;
          }
        }
      }
      if ($n_introns > $max_n_introns){
        @candidates = ();
        push(@candidates, $db_gene);
        $max_n_introns = $n_introns;  #print "0-AAA ".$transcript->start."-".$transcript->end."  ".$db_gene->stable_id."  ".$n_introns."\n";
      }
      elsif ($n_introns >= 1 and $n_introns == $max_n_introns){
        push(@candidates, $db_gene);  #print "0-BBB ".$transcript->start."-".$transcript->end."  ".$db_gene->stable_id."  ".$n_introns."\n";    
      }
    }
    #If only one...
    if (scalar(@candidates) == 1){
      $host_gene = $candidates[0]; #print "AAA ".$transcript->start."-".$transcript->end."\n";
    }
    #Else, check exon overlap - find host gene with the most exon overlap
    else{
      if (scalar(@candidates) == 0){
        @candidates = grep {!($_->biotype =~ /pseudogene/)} @db_genes; #print "BBB ".$transcript->start."-".$transcript->end."\n";
      }
      my %tr_loc; #Store genomic positions of transcript's exons
      foreach my $exon (@{$transcript->get_all_Exons}){
        for (my $i = $exon->seq_region_start; $i <= $exon->seq_region_end; $i++){
          $tr_loc{$i} = 1;
        }
      }
      my $max = 0;
      foreach my $db_gene (@candidates){
        my %seen; #Store genomic positions where transcript's exons overlap host gene's exons
        foreach my $db_exon (@{$db_gene->get_all_Exons}){
          for (my $i = $db_exon->seq_region_start; $i <= $db_exon->seq_region_end; $i++){
            if ($tr_loc{$i}){
              $seen{$i} = 1;
            }
          }
        }
        if (scalar(keys %seen) > $max){
          $host_gene = $db_gene;
          $max = scalar(keys %seen);
        }
      }                
    }
    #Store most suitable host gene for the transcript, if any
    if ($host_gene){
      my $tid = $transcript->get_all_Attributes('hidden_remark')->[0]->value;
      $host_by_transcript{$tid} = $host_gene->stable_id;    
      print "Assigning transcript $tid to host gene ".$host_gene->stable_id."\n";
    }
    my $t_name = $transcript->stable_id || $transcript->get_all_Attributes('hidden_remark')->[0]->value;
    print "XCXCXCXCX-2\n" if $t_name =~ /anchIC_000000168991/;      
  }
    
  #How many different host genes were found? Split input gene if necessary
  my @host_gene_ids = uniq(values %host_by_transcript);
  my @host_genes;
  foreach my $gid (@host_gene_ids){
    push(@host_genes, grep {$_->stable_id eq $gid} @db_genes);
  }
  #Multiple host genes: split the input gene into as many new genes as host genes there are
  if (scalar @host_gene_ids > 1){
    #Make new empty genes
    my @new_genes;
    foreach my $host_gene_id (@host_gene_ids){
      my $new_gene = Bio::Vega::Gene->new(
                                    -stable_id  => $host_gene_id,
                                    -biotype    => $gene->biotype,
                                    -analysis   => $gene->analysis,
                                    -source     => $gene->source,
                                    -status     => $gene->status
                                    );
      $new_gene->gene_author($gene->gene_author);
      $new_gene->add_Attributes(@{$gene->get_all_Attributes});            
      push (@new_genes, $new_gene);
    }
    #Add input transcripts to new genes
    foreach my $transcript (@{$gene->get_all_Transcripts}){
      my $tid = $transcript->get_all_Attributes('hidden_remark')->[0]->value;
      if ($host_by_transcript{$tid}){
        foreach my $new_gene (@new_genes){
          if ($new_gene->stable_id eq $host_by_transcript{$tid}){
            $new_gene->add_Transcript($transcript);
            last;
          }                
        }            
      }
      #If transcript had no host gene, take the closest or max overlap host gene
      else{
        my $host_gene;
        #First, find overlapping genes
        my @ovlp_host_genes = grep {!($_->biotype =~ /pseudogene/)} grep {$_->seq_region_start <= $transcript->seq_region_end and $_->seq_region_end >= $transcript->seq_region_start} @host_genes;
        if (scalar @ovlp_host_genes){
          my $max_overlap = 0;
          foreach my $db_gene (@ovlp_host_genes){
            my $overlap = min($transcript->seq_region_end, $db_gene->seq_region_end) - max($transcript->seq_region_start, $db_gene->seq_region_start) + 1;
            if ($overlap > $max_overlap){
              $host_gene = $db_gene;
              $max_overlap = $overlap;
            }
          }
        }
        else{
          my $min_distance = 10000000;
          foreach my $db_gene (@host_genes){
            my $distance = max($transcript->seq_region_start, $db_gene->seq_region_start) - min($transcript->seq_region_end, $db_gene->seq_region_end) + 1;
            if ($distance < $min_distance){
            #22-1-2024      $host_gene = $db_gene;  #Commented out because there must be overlap between transcript and host gene
              $min_distance = $distance;
            }
          }
        }
        if ($host_gene){
          foreach my $new_gene (@new_genes){
            if ($new_gene->stable_id eq $host_gene->stable_id){
              $new_gene->add_Transcript($transcript);
              last;
            }                
          }
        }
        else{
          print "No host gene for transcript $tid\n";
        }
      }
    }
    return \@new_genes;
  }
  #A single host gene: just set stable id and return
  elsif (scalar @host_gene_ids == 1){
    $gene->stable_id($host_gene_ids[0]);
    return [$gene];    
  }
  else{
    return [$gene];    
  }
    
  return undef;
}



=head2 assign_host_gene_multi

 Arg[1]    : list of Bio::Vega::Gene objects
 Arg[2]    : preserve host gene biotype, eg. ignore non-coding host gene if transcript is coding
 Arg[3]    : list of ignored host gene biotypes, i.e. transcripts will not be added to a gene with one of these
             biotypes but may form part of an overlapping novel gene
 Function  : Find the most suitable existing gene(s) to host the transcripts of the given genes, 
             split the input genes into as many genes as host genes were found and set the new gene stable id(s);
             return the transformed input genes
 Returntype: list of Bio::Vega::Gene objects

=cut

sub assign_host_gene_multi {
  my ($self, $genes, $preserve_biotype, $ignored_biotypes) = @_;
  #my %host_by_transcript;
  my @all_db_genes;
  my @new_genes;
  my %ignored_biotypes;
  if ($ignored_biotypes){
    %ignored_biotypes = map {$_=>1} split(/,/, $ignored_biotypes);
  }

  #Sort input genes by decreasing genomic span - hopefully it will improve the host gene assignment
  foreach my $gene (sort {$b->length <=> $a->length} @$genes){
    #Fetch database genes overlapping our gene
    #Exclude 'not for VEGA' genes, artifact genes, ...
    my @db_genes = @{get_valid_overlapping_genes($gene, 1)};
    #Add genes created for other transcripts in the same cluster
    push(@db_genes, @new_genes);
    print "DB_GENES = ".join(",", map {$_->stable_id} @db_genes)."\n";
    my @filt_db_genes;
    DBG:foreach my $db_gene (@db_genes){
      if ($db_gene->seq_region_strand == $gene->seq_region_strand){
        #Ignore readthrough genes as candidate host genes
        foreach my $db_tr (@{$db_gene->get_all_Transcripts}){
          if (scalar grep {$_->value eq "readthrough"} @{$db_tr->get_all_Attributes("remark")}){
            next DBG;
          }
        }
        #Ignore genes with certain biotypes as candidate host genes
        unless ($ignored_biotypes{$db_gene->biotype}){
          push(@filt_db_genes, $db_gene);
        }
      }
    }
    print "FILT_DB_GENES = ".join(",", map {$_->stable_id} @filt_db_genes)."\n";
    
    foreach my $transcript (@{$gene->get_all_Transcripts}){
      my $host_gene;
      my $tid = $transcript->get_all_Attributes('hidden_remark')->[0]->value;

      ###
      #Test to manage cases involving overlapping genes, i.e. where a readthrough transcript from an adjacent gene
      #may be chosen instead of the expected gene.
      #Give priority to candidate host genes that have MANE Select transcripts with exonic overlap with the novel transcript.
      #Would that be right?
      #This is potentially dangerous, especially for lncRNA genes which don't have MANE Select transcripts
      #and could be ignored if there is a small overlap with coding genes.
      #It should only be applied if the top two candidates have a MANE Select.
      #It should require at least a common splice site.
      
      #TO DO:
      #1-Find extent of exonic overlap, number of shared splice sites, overlap with MANE_Select, with all the overlapping genes
      #2-Filter overlapped genes: if all exon-overlapped genes have a MANE_Select, filter for those where the transcript exon-overlaps the MANE_Select
      #3-Find the best candidate
      ###

      #Alternative - more conservative approach: if all the overlapping genes with exonic overlap have a MANE Select,
      #filter for those genes whose MANE Select transcript is overlapped by the novel transcript at the exon level 

      my @exon_overlapping_genes;
      DBG:foreach my $db_gene (@filt_db_genes){
        foreach my $db_transcript (@{$db_gene->get_all_Transcripts}){
          foreach my $db_exon (@{$db_transcript->get_all_Exons}){
            foreach my $exon (@{$transcript->get_all_Exons}){
              if ($db_exon->seq_region_start <= $exon->seq_region_end and $db_exon->seq_region_end >= $exon->seq_region_start){
                push(@exon_overlapping_genes, $db_gene);                 
                next DBG;
              }
            }
          } 
        }
      }
      print "EXON_OVLP_GENES = ".join(",", map {$_->stable_id} @exon_overlapping_genes)."\n";
      if (scalar @exon_overlapping_genes){
        @filt_db_genes = @exon_overlapping_genes;
      }
      print "FILT_DB_GENES_2 = ".join(",", map {$_->stable_id} @filt_db_genes)."\n";
      
      my @db_genes_with_mane_select;
      DBG:foreach my $db_gene (@filt_db_genes){
        foreach my $db_transcript (@{$db_gene->get_all_Transcripts}){
          if (scalar grep {$_->value eq "MANE_select"} @{$db_transcript->get_all_Attributes("remark")}){
            push(@db_genes_with_mane_select, $db_gene);
            next DBG;
          }
        }
      }
      print "GENES_WITH_SELECT = ".join(",", map {$_->stable_id} @db_genes_with_mane_select)."\n";
      #if (scalar @filt_db_genes == scalar @db_genes_with_mane_select){
      #Compare the numbers of unique stable ids as some genes may be stored twice in the arrays
      #Check if all the genes overlapped at the exon level have a MANE Select transcript
      #If so, reset the host gene candidate set to include only genes whose MANE Select transcripts are overlapped at the exon level
      if (scalar uniq(map {$_->stable_id} @filt_db_genes) == scalar uniq(map {$_->stable_id} @db_genes_with_mane_select)){ 
        my @mane_select_overlapping_genes;
        DBG:foreach my $db_gene (@filt_db_genes){
          foreach my $db_transcript (@{$db_gene->get_all_Transcripts}){
            if (scalar grep {$_->value eq "MANE_select"} @{$db_transcript->get_all_Attributes("remark")}){
              foreach my $db_exon (@{$db_transcript->get_all_Exons}){
                foreach my $exon (@{$transcript->get_all_Exons}){
                  if ($db_exon->seq_region_start <= $exon->seq_region_end and $db_exon->seq_region_end >= $exon->seq_region_start){
                    push(@mane_select_overlapping_genes, $db_gene);                 
                    next DBG;
                  }
                }
              }
            } 
          }
        }
        if (scalar @mane_select_overlapping_genes){
          @filt_db_genes = @mane_select_overlapping_genes;
        }
      }


      #Check intron match - find gene with the most matching introns
      my @candidates;
      my $max_n_introns = 0;
      foreach my $db_gene (@filt_db_genes){
        #Exclude pseudogenes so that spliced transcripts are assigned to new (lncRNA) genes
        #(although pseudogene biotypes may already be included in the set of ignored biotypes)
        if (scalar @{$transcript->get_all_Introns} > 1 and $db_gene->biotype =~ /pseudogene/){
          next;
        }               
        #If input transcript has a translation, skip non-coding host genes
        if ($preserve_biotype){
          if ($transcript->translation and $db_gene->biotype ne "protein_coding"){
            next;
          }
        }
        my $n_introns = 0; #Number of matching introns
        my %seen_db_introns;
        foreach my $db_transcript (@{$db_gene->get_all_Transcripts}){
          next if scalar (grep {$_->value eq "not for VEGA"} @{$db_transcript->get_all_Attributes('remark')});
          next if $db_transcript->biotype eq "artifact";
          foreach my $db_intron (@{$db_transcript->get_all_Introns}){
            next if $seen_db_introns{$db_intron->seq_region_start."-".$db_intron->seq_region_end};
            $seen_db_introns{$db_intron->seq_region_start."-".$db_intron->seq_region_end} = 1;
            foreach my $intron (@{$transcript->get_all_Introns}){
              if ($db_intron->seq_region_start == $intron->seq_region_start and $db_intron->seq_region_end == $intron->seq_region_end){
                $n_introns++;
              }
            }
          }
        }
    
        if ($n_introns > $max_n_introns){
          @candidates = ();
          push(@candidates, $db_gene);
          $max_n_introns = $n_introns;
        }
        elsif ($n_introns >= 1 and $n_introns == $max_n_introns){
          push(@candidates, $db_gene);   
        }
      }

      #If only one host gene with the highest number of matching introns, choose that
      if (scalar(@candidates) == 1){
        $host_gene = $candidates[0];
      }
      #Else, select the host gene with the largest exonic overlap
      else{
        if (scalar(@candidates) == 0){
          @candidates = grep {!($_->biotype =~ /pseudogene/)} @filt_db_genes;
        }
        my @final_candidates;
        my %tr_loc; #Store genomic positions of transcript's exons
        foreach my $exon (@{$transcript->get_all_Exons}){
          for (my $i = $exon->seq_region_start; $i <= $exon->seq_region_end; $i++){
            $tr_loc{$i} = 1;
          }
        }
        my $max = 0;
        foreach my $db_gene (@candidates){
          my %seen_pos; #Store genomic positions where transcript's exons overlap host gene's exons
          foreach my $db_transcript (@{$db_gene->get_all_Transcripts}){
            next if scalar (grep {$_->value eq "not for VEGA"} @{$db_transcript->get_all_Attributes('remark')});
            next if $db_transcript->biotype eq "artifact";
            foreach my $db_exon (@{$db_transcript->get_all_Exons}){
              for (my $i = $db_exon->seq_region_start; $i <= $db_exon->seq_region_end; $i++){
                if ($tr_loc{$i}){
                  $seen_pos{$i} = 1;
                }
              }
            }
          }
          if (scalar(keys %seen_pos) > $max){
            @final_candidates = ($db_gene);
            $max = scalar(keys %seen_pos);
          }
          elsif (scalar(keys %seen_pos) == $max and $max > 0){
            push(@final_candidates, $db_gene);
          }
        }
        if (scalar @final_candidates == 1){
          $host_gene = $final_candidates[0];
        }
        else{
          #If two host genes are tied, prioritise the gene with the oldest ENS stable id
          #or the longest gene if none has ENS stable id
          if (scalar grep {$_->stable_id =~ /^ENS/} @final_candidates){
            ($host_gene) = sort {$a->stable_id cmp $b->stable_id} grep {$_->stable_id =~ /^ENS/} @final_candidates;
          }
          else{
            ($host_gene) = sort {$b->length <=> $a->length} @final_candidates;
          }
        }
        
      }
      #Store most suitable host gene for the transcript, if any
      my $found = 0;
      if ($host_gene){   
        print "Assigning transcript $tid to host gene ".$host_gene->stable_id."\n";
        
        #Search for new gene with this stable id
        #i.e. a new gene created for another input transcript which had the same host gene
        #my $found = 0;
        foreach my $new_gene (@new_genes){
          if ($new_gene->stable_id eq $host_gene->stable_id){
            $new_gene->add_Transcript($transcript);
            $found = 1;
            last;
          }
        }
        unless ($found){
          my $new_gene = Bio::Vega::Gene->new(
                                    -stable_id  => $host_gene->stable_id,
                                    -biotype    => $gene->biotype,
                                    -analysis   => $gene->analysis,
                                    -source     => $gene->source,
                                    -status     => $gene->status
                                  );
          $new_gene->gene_author($gene->gene_author);
          $new_gene->add_Attributes(@{$gene->get_all_Attributes});
          $new_gene->add_Transcript($transcript);   
          push (@new_genes, $new_gene);
        }
        #Remove transcript from old gene
        my $array = $gene->{_transcript_array};
        @$array = grep {$_->get_all_Attributes('hidden_remark')->[0]->value ne $tid} @$array;
      }
      else{
        #If no host gene, ...
        print "No host gene for transcript $tid\n";
        #Assign random identifier and add to db genes array so other transcripts can be assigned to this gene
        my $new_gene = Bio::Vega::Gene->new(
                                    -stable_id  => "tmp_".int(rand(1000)),
                                    -biotype    => $gene->biotype,
                                    -analysis   => $gene->analysis,
                                    -source     => $gene->source,
                                    -status     => $gene->status
                                  );
        $new_gene->gene_author($gene->gene_author);
        $new_gene->add_Attributes(@{$gene->get_all_Attributes});
        $new_gene->add_Transcript($transcript);   
        push (@new_genes, $new_gene);
        #Remove transcript from old gene
        my $array = $gene->{_transcript_array};
        @$array = grep {$_->get_all_Attributes('hidden_remark')->[0]->value ne $tid} @$array;
              
        push(@db_genes, $new_gene); 
              
      }

    }
    #Add the old gene to the set of new genes if it still has any unassigned transcript
    if (scalar @{$gene->get_all_Transcripts}){
      push(@new_genes, $gene);
    }
  }

  return \@new_genes;

}



=head2 get_valid_overlapping_genes

 Arg[1]    : Bio::Vega::Gene or Bio::Vega::Transcript object
 Arg[2]    : include database features with the same stable id as the input feature
 Function  : Returns a list of genes that overlap the input feature on the same strand and would normally be part of the Ensembl release, i.e. the genes that are not artifact or "not for VEGA" and have at least a transcript that is neither of these.
 Returntype: list of Bio::Vega::Gene objects

=cut

sub get_valid_overlapping_genes {
  my ($feat, $incl_self) = @_;
  my @genes = ();
  my $slice = $feat->slice->adaptor->fetch_by_region("toplevel", $feat->seq_region_name, $feat->start, $feat->end); #print "FEAT=".$feat->seq_region_name.":".$feat->start."-".$feat->end."\n";
  foreach my $gene (@{$slice->get_all_Genes}){
    #print "OVLP=".$gene->stable_id." ".$gene->seq_region_name." ".$gene->seq_region_start."\n";
    if ($gene->seq_region_strand == $feat->strand and
        $gene->source =~ /(ensembl|havana)/ and
        $gene->biotype ne "artifact" and
        $gene->biotype ne "comp_pipe" and        
        scalar (grep {$_->value eq "not for VEGA"} @{$gene->get_all_Attributes('remark')}) == 0){
      if ($gene->stable_id eq $feat->stable_id){
        next unless $incl_self;
      }
      TR:foreach my $transcript (@{$gene->get_all_Transcripts}){
        if ($transcript->source =~ /(ensembl|havana)/ and
            $transcript->biotype ne "artifact" and
            $transcript->biotype ne "comp_pipe" and        
            scalar (grep {$_->value eq "not for VEGA"} @{$transcript->get_all_Attributes('remark')}) == 0){	  
          push(@genes, $gene); #print "Adding gene ".$gene->stable_id."\n";
          last TR;
        }
      }
    }
  }
  return \@genes;
}



=head2 reassign_host_genes

 Arg[1]    : list of Bio::Vega::Gene objects
 Arg[2]    : include genes from the database
 Arg[3]    : hashref of allowed biotypes for the host genes
 Function  : Finds the most suitable gene(s) to host the transcripts of the given genes, taking into account the whole set of transcripts in all genes.
             This should improve the original gene assignment and minimise the number of overlapping novel genes that should be merged
 Returntype: list of Bio::Vega::Gene objects

=cut

sub reassign_host_genes {
  my ($self, $genes, $incl_dbgenes, $allowed_biotypes) = @_;

  #Find mergeable genes (must share a splice site)
  my %mergeable_genes;
  my @gene_clusters;
  foreach my $gene (@$genes){
    my @other_genes = grep {$_->stable_id ne $gene->stable_id and $_->seq_region_strand == $gene->seq_region_strand} @$genes;
    if ($incl_dbgenes){
      push(@other_genes, @{get_valid_overlapping_genes($gene)});
    }
    G2:foreach my $gene2 (@other_genes){
      if (scalar keys %{$allowed_biotypes}){
        next unless $allowed_biotypes->{$gene2->biotype};
      }
      foreach my $transcript (@{$gene->get_all_Transcripts}){
        foreach my $transcript2 (@{$gene2->get_all_Transcripts}){
          next if $transcript2->biotype eq "artifact" or $transcript->biotype eq "comp_pipe" or scalar (grep {$_->value eq "not for VEGA"} @{$transcript->get_all_Attributes('remark')}); 
          foreach my $intron (@{$transcript->get_all_Introns}){
            foreach my $intron2 (@{$transcript2->get_all_Introns}){
              if ($intron->seq_region_start == $intron2->seq_region_start or $intron->seq_region_end == $intron2->seq_region_end){
                my $seen = 0;
                foreach my $cluster (@gene_clusters){
                  if (grep {$_ eq $gene->stable_id or $_ eq $gene2->stable_id} @{$cluster}){
                    my @extended = uniq(@{$cluster}, $gene->stable_id, $gene2->stable_id);
                    $cluster = \@extended;
                    $seen = 1;
                  }
                }
                unless ($seen){
                  my @cluster = ($gene->stable_id, $gene2->stable_id);
                  push(@gene_clusters, \@cluster);
                }
                $mergeable_genes{$gene->stable_id}{$gene2->stable_id} = 1;
                print "MERGEABLE\t".$gene->stable_id."\t".$gene2->stable_id."\n";
                next G2;
              }
            }
          }
        }
      }
    }
  }

  foreach my $gene_id (keys %mergeable_genes){
    if (!($gene_id =~ /^ENS/)){ #do not merge an existing gene into another (yet)
      #if (scalar keys %{$mergeable_genes{$gene_id}} == 1){
        my ($gene) = grep {$_->stable_id eq $gene_id} @$genes;
        my ($target_gene_id) = sort {$b cmp $a} keys %{$mergeable_genes{$gene_id}}; #prioritise novel genes (tmp_*) vs existing ones (ENS*)
        my ($target_gene) = grep {$_->stable_id eq $target_gene_id} @$genes;
        if ($target_gene){
          next unless scalar @{$target_gene->get_all_Transcripts}; #do not return transcripts to a gene that was just emptied
        }
        print "Merging gene $gene_id into gene $target_gene_id\n";
        if ($target_gene){
          foreach my $transcript (@{$gene->get_all_Transcripts}){
            $target_gene->add_Transcript($transcript);
          }
          $gene->flush_Transcripts;
        }
        else{
          $gene->stable_id($target_gene_id);
        }
      #}
    }
  }
  #return non-empty genes only
  @$genes = grep {scalar @{$_->get_all_Transcripts}} @$genes;

  return $genes;
}



1;
