
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
#            my $tr_slice = $transcript->slice->adaptor->fetch_by_region("toplevel", $transcript->seq_region_name, $transcript->start, $transcript->end);
#            my @db_genes = grep {$_->seq_region_strand == $gene->seq_region_strand} @{$tr_slice->get_all_Genes};
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



=head2 check_overlapped_loci

 Arg[1]    : list of Bio::Vega::Gene objects
 Arg[2]    : Integer
 Function  : Check for transcripts overlapping a number of existing genes at the exon level that exceeds the number given in the second argument.
 Returntype: Hashref - list of transcript names

=cut

sub check_overlapped_loci {
    my ($self, $genes, $max_ov_loci) = @_;
    my %list;
    foreach my $gene (@$genes){
        foreach my $transcript (@{$gene->get_all_Transcripts}){
            my $message = "";
            #Fetch database genes overlapping our transcript
#            my $tr_slice = $transcript->slice->adaptor->fetch_by_region("toplevel", $transcript->seq_region_name, $transcript->start, $transcript->end);
#            my @db_genes = grep {$_->seq_region_strand == $gene->seq_region_strand} @{$tr_slice->get_all_Genes};
            #my $overlapped_genes_count = 0;
            my @db_genes = @{get_valid_overlapping_genes($transcript)};
            my @overlapped_db_genes;
            DBG:foreach my $db_gene (@db_genes){
#                next unless $db_gene->source =~ /(ensembl|havana)/;
#                next if $db_gene->biotype eq "artifact";
#                next if scalar(grep {$_->value eq "not for VEGA"} @{$db_gene->get_all_Attributes('remark')}) > 0;
               # next if $db_gene->biotype =~ /^(miRNA|misc_RNA|ribozyme|rRNA|rRNA_pseudogene|scaRNA|snoRNA|snRNA|sRNA|vault_rna)$/i;  #Ignore small RNA genes
               # next if $db_gene->biotype =~ /ig_pseudogene|processed_pseudogene|pseudoexon|pseudogene|transcribed_.+_pseudogene|tr_pseudogene|unitary_pseudogene/; #Ignore pseudogenes
                if ($db_gene->biotype =~ /antisense|bidirectional_promoter_lncrna|ig_gene|IG_V_gene|lincRNA|polymorphic_pseudogene|processed_transcript|protein_coding|sense_intronic|sense_overlapping|tec/i){
                    foreach my $db_exon (@{$db_gene->get_all_Exons}){
                        foreach my $exon (@{$transcript->get_all_Exons}){
                            if ($db_exon->seq_region_start <= $exon->seq_region_end and $db_exon->seq_region_end >= $exon->seq_region_start){
                                push(@overlapped_db_genes, $db_gene);
                                #$overlapped_genes_count++; 
                                $message .= "OV: ".$db_gene->stable_id." ".$db_gene->biotype."  ";
                                next DBG;
                            }
                        }
                    }
                }
            }
            #Instead of counting overlapped genes, count clusters of overlapped genes
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
                print "KILL_2: ".$t_name."  ".$transcript->start."-".$transcript->end." overlaps ".scalar(@clusters)." loci: $message\n";
                $list{$t_name} = 1;
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
#    my $gene_slice = $gene->slice->adaptor->fetch_by_region("toplevel", $gene->seq_region_name, $gene->start, $gene->end);
#    my @db_genes = grep {$_->stable_id ne $gene->stable_id and $_->seq_region_strand == $gene->strand} @{$gene_slice->get_all_Genes};
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
print "AAAAA:".$gene->seq_region_start."\n";
    #Fetch database genes overlapping our gene
#    my $gene_slice = $gene->slice->adaptor->fetch_by_region("toplevel", $gene->seq_region_name, $gene->start, $gene->end);
#    my @db_genes = grep {$_->stable_id ne $gene->stable_id and $_->seq_region_strand == $gene->strand} @{$gene_slice->get_all_Genes};
    my @db_genes = @{get_valid_overlapping_genes($gene)};
    my %host_by_transcript;
    foreach my $transcript (@{$gene->get_all_Transcripts}){
print "BBBBB:".$transcript->seq_region_start."\n";
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
                @candidates = @db_genes; #print "BBB ".$transcript->start."-".$transcript->end."\n";
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
                my @ovlp_host_genes = grep {$_->seq_region_start <= $transcript->seq_region_end and $_->seq_region_end >= $transcript->seq_region_start} @host_genes;
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
                            $host_gene = $db_gene;
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



=head2 get_valid_overlapping_genes

 Arg[1]    : Bio::Vega::Gene or Bio::Vega::Transcript object
 Function  : Returns a list of genes that overlap the input feature on the same feature and would normally be part of the Ensembl release, i.e. the genes are not artifact or "not for VEGA" and have at least a transcript that is neither of these.
 Returntype: list of Bio::Vega::Gene objects

=cut

sub get_valid_overlapping_genes {
  my ($feat) = @_;
  my @genes = ();
  my $slice = $feat->slice->adaptor->fetch_by_region("toplevel", $feat->seq_region_name, $feat->start, $feat->end);
  foreach my $gene (@{$slice->get_all_Genes}){
    if ($gene->stable_id ne $feat->stable_id and 
        $gene->seq_region_strand == $feat->strand and
        $gene->source =~ /(ensembl|havana)/ and
        $gene->biotype ne "artifact" and
        $gene->biotype ne "comp_pipe" and        
        scalar (grep {$_->value eq "not for VEGA"} @{$gene->get_all_Attributes('remark')}) == 0){
      TR:foreach my $transcript (@{$gene->get_all_Transcripts}){
        if ($transcript->source =~ /(ensembl|havana)/ and
            $transcript->biotype ne "artifact" and
            $transcript->biotype ne "comp_pipe" and        
            scalar (grep {$_->value eq "not for VEGA"} @{$transcript->get_all_Attributes('remark')}) == 0){	  
          push(@genes, $gene);
          last TR;
        }
      }
    }
  }
  return \@genes;
}


1;
