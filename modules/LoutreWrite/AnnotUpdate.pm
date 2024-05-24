
package LoutreWrite::AnnotUpdate;

use strict;
use warnings;
use Bio::Vega::Gene;
use Bio::Vega::Transcript;
use Bio::Vega::Exon;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Attribute;
use Bio::Otter::ServerAction::Script::Region;
use Bio::Otter::Server::Support::Local;
use Sys::Hostname;
use Try::Tiny;
use List::Util qw[min max];
use List::MoreUtils qw(uniq);
use LoutreWrite::CDSCreation qw(assign_cds_to_transcripts has_complete_cds);
use LoutreWrite::Config;
use LoutreWrite::GeneFilter qw(get_valid_overlapping_genes);
use base 'Exporter';

our @EXPORT = qw( $WRITE );
our @EXPORT_OK = qw( exon_novelty intron_novelty can_be_merged merge_transcripts has_polyA_site_support has_polyAseq_support );
our $WRITE = 0;
our $CP_BIOTYPE = "comp_pipe";





=head2 process_gene_2

 Arg[1]    : Vega gene object to be written
 Arg[2]    : mode - either 'add' or 'update'
 Arg[3]    : submode - either 'cons' or 'aggr'
 Arg[4]    : Otter dataset name ('human', 'mouse', 'human_test', etc)
 Arg[5]    : Bio::Vega::DBSQL::DBAdaptor
 Arg[6]    : Boolean - if true, do not check for completely or partially identical intron chains in the annotation
 Arg[7]    : Boolean - force 'comp_pipe' biotype
 Arg[8]    : Boolean - if true, the "not for VEGA" remark will not be added
 Arg[9]    : Boolean - if true, do not try to add a CDS
 Arg[10]   : Boolean - if true, add a "platinum" hidden remark
 Function  : get region with the gene coordinates; 'add' the gene to the region or 'update' pre-existing gene annotation; store changes in the database 
 Returntype: String (message with the outcome)

=cut

sub process_gene_2 {
    my ($self, $gene, $mode, $submode, $dataset_name, $dba, $no_intron_check, $force_cp_biotype, $no_NFV, $do_not_add_cds, $platinum) = @_;
    my @log;
print "SUBMODE: $submode\n";

    #Get default coordinate system
    my $csa = $dba->get_CoordSystemAdaptor;
    my $csver = $csa->fetch_top_level->version;

    #If in "update" mode, the region has to span the whole host gene to stop this from being "spliced out" as an incomplete gene. 
    #That would prevent the new annotation from being added to the host gene.
    my ($padded_start, $padded_end);
    if ($mode eq "update"){
        my $host_gene = $dba->get_GeneAdaptor->fetch_by_stable_id($gene->stable_id); unless ($host_gene){print "No gene for ".$gene->stable_id."\n";}
        $padded_start = min($gene->start, $host_gene->start);
        $padded_end = max($gene->end, $host_gene->end);
    }

    #Get Otter region for the gene
    my $local_server = Bio::Otter::Server::Support::Local->new( otter_dba => $dba );
    $local_server->set_params(
        dataset => $dataset_name,
        chr     => $gene->seq_region_name,
        start   => $padded_start || $gene->start,   
        end     => $padded_end || $gene->end,
        cs      => "toplevel",
        csver   => $csver
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

    if ($mode eq "add"){
        #Create a new gene name
        my $new_gene_name = $self->get_new_gene_name($region, $gene, $dataset_name, $dba) or return " ";
        my $name_att = Bio::EnsEMBL::Attribute->new(-code => 'name', -value => $new_gene_name);
        $gene->add_Attributes($name_att);
        foreach my $tr (@{$gene->get_all_Transcripts}){
            #Create a new transcript name
            my $new_tr_name = $self->get_new_transcript_name($gene, $dba);
            my $tr_name_att = Bio::EnsEMBL::Attribute->new(-code => 'name', -value => $new_tr_name);
            $tr->add_Attributes($tr_name_att);
            my $id = $tr->get_all_Attributes('hidden_remark')->[0]->value;
            if ($platinum){
              $tr->add_Attributes( Bio::EnsEMBL::Attribute->new(-code => 'hidden_remark', -value => 'Platinum set') );
            }
            print "TR: $id: Will add transcript $new_tr_name to novel gene $new_gene_name\n";
            push (@log, "TR2: $id: Added transcript $new_tr_name to novel gene $new_gene_name");
        }
        $gene = $self->assign_biotypes($gene);
        #Add new gene to region
        $region->add_genes($gene);
    }

    elsif ($mode eq "update"){
        #Update existing gene: host gene that provided its stable id to the gene being stored
        if (scalar(grep {$_->stable_id eq $gene->stable_id} $region->genes) == 1){
            my ($db_gene) = grep {$_->stable_id eq $gene->stable_id} $region->genes;
            
            #Compare new transcripts with existing ones
            TR:foreach my $tr (@{$gene->get_all_Transcripts}){
                my $id = $tr->get_all_Attributes('hidden_remark')->[0]->value;
                print "\n\n#Looking at transcript $id\n";
                my $add_transcript = 0;
                my @merge_candidates = ();
                my %tr_comp;
                print "\n\nComparing with existing transcripts:\n";
                print "GSID ".$db_gene->stable_id."\n";              
                #First round: look for micro-intron overlaps (false partial intron retention due to misalignment)
                foreach my $db_tr (@{$db_gene->get_all_Transcripts}){
                    #Exclude transcripts not yet stored in the database
  #                  next unless $db_tr->stable_id;
                    #Exclude artifacts
                    next if $db_tr->biotype eq "artifact";
                    #Exclude unwanted transcript biotypes which may be associated to unexpected gene biotypes
                    next if $db_tr->biotype =~ /^(known_ncrna|macro_lncrna|mirna|misc_RNA|ribozyme|rrna|scrna|snorna|snrna|vaultrna|tr_gene)$/;
                    #Ignore lingering OTT transcripts
                    next if $db_tr->stable_id =~ /^OTT/;
                    #Ignore "not for VEGA" transcripts unless they have a "comp_pipe" biotype or a "TAGENE_transcript" remark
                    if (scalar(grep {$_->value eq "not for VEGA"} @{$db_tr->get_all_Attributes('remark')})){
                      next unless ($db_tr->biotype eq "comp_pipe" or scalar(grep {$_->value eq "TAGENE_transcript"} @{$db_tr->get_all_Attributes('remark')}));
                    }
                    #Exclude transcripts having a "comp_pipe_rejected" remark
                    next if scalar(grep {$_->value eq "comp_pipe_rejected"} @{$db_tr->get_all_Attributes('remark')});
                    #Ignore single-exon transcripts (?)
                    next if scalar @{$db_tr->get_all_Introns} == 0;  #Do not merge with single-exon transcripts? (In case they are CLS Platinum...)
                    
                    #######---------########------#######-------#######  i=0, e=1, m=0 ; e=0
                    #######-----------------------#######-------#######  i=1, e=1, m=0
                                ############------#######                i=1, e=1, m=0
                    #######---------########------########
                
                    #Compare coordinates: transcript's start/end vs. db transcript's internal exon
                    #Search for a terminal exon matching an existing internal exon except for a small overhang at the end
                    #Trim transcript start/end if small overhang detected
                    #Skip this if Platinum? 
                   # unless ($platinum){
                      $tr = clip_ends($tr, $db_tr, 10, $slice_offset);
                   # }

                    #Loop over all transcripts in case a second db transcript supports the transcript's start/end (exclude known retained introns?)
                    
                }
                #IDEALLY, CREATE A METHOD FOR CHECKING NOVELTY SO THAT IT IS EASIER TO RE-RUN IT IF A TRANSCRIPT IS CLIPPED

                
                ###
                my $outcome;
                foreach my $db_tr (@{$db_gene->get_all_Transcripts}){
                    #Exclude transcripts not yet stored in the database
#                    next unless $db_tr->stable_id;
                    #Exclude artifacts
                    next if $db_tr->biotype eq "artifact";
                    #Exclude unwanted transcript biotypes which may be associated to unexpected gene biotypes
                    next if $db_tr->biotype =~ /^(known_ncrna|macro_lncrna|mirna|misc_RNA|ribozyme|rrna|scrna|snorna|snrna|vaultrna|tr_gene)$/;
                    #Ignore lingering OTT transcripts
                    next if $db_tr->stable_id =~ /^OTT/;
                    #Ignore "not for VEGA" transcripts unless they have a "comp_pipe" biotype or a "TAGENE_transcript" remark
                    if (scalar(grep {$_->value eq "not for VEGA"} @{$db_tr->get_all_Attributes('remark')})){
                      next unless ($db_tr->biotype eq "comp_pipe" or scalar(grep {$_->value eq "TAGENE_transcript"} @{$db_tr->get_all_Attributes('remark')}));
                    }
                    #Exclude transcripts having a "comp_pipe_rejected" remark
                    next if scalar(grep {$_->value eq "comp_pipe_rejected"} @{$db_tr->get_all_Attributes('remark')});
                    #Ignore single-exon transcripts unless TAGENE transcript is single-exon too?
                    #No: for the moment, just allow adding multi-exon transcripts to existing single-exon genes         
                    if (scalar @{$tr->get_all_Exons} > 1){
                  #    next if scalar @{$db_tr->get_all_Exons} == 1;  #Do not merge with single-exon transcripts? (In case they are CLS Platinum...)
                    }

                    my $i = intron_novelty($tr, $db_tr);
                    my $e = exon_novelty($tr, $db_tr);
                    my $m = can_be_merged($tr, $db_tr);
                    ###
                 #   my $c = can_be_clipped($tr, $db_tr, 15);
                    ###
                    my $db_tr_id = $db_tr->stable_id || $db_tr->get_all_Attributes('name')->[0]->value;
                    $tr_comp{$db_tr_id} = {'intron' => $i, 'exon' => $e, 'merge' => $m};
                    print "COMP: i=$i, e=$e, m=$m ".$db_tr_id." ".$db_tr->biotype."\n";
                }
                #Multi-exon transcripts
                if (scalar @{$tr->get_all_Exons} > 1){
                    #Conservative mode: a new transcript model is only stored if it has a novel intron structure or an existing intron structure with novel sequence in its terminal exons; existing transcripts are never updated.
                    if ($submode eq "cons"){
                        #Add transcript if intron chain shows novelty w.r.t. all existing transcripts
                        #Discard transcript if an existing transcript has identical intron chain
                        DBTR:foreach my $db_tr_id (keys %tr_comp){
                            if ($tr_comp{$db_tr_id}->{'intron'} == 1 or
                               ($tr_comp{$db_tr_id}->{'intron'} == -1 and $tr_comp{$db_tr_id}->{'exon'} == 1)){
                                $add_transcript = 1;
                            }
                            elsif (($tr_comp{$db_tr_id}->{'intron'} == 0) or 
                                   ($tr_comp{$db_tr_id}->{'intron'} == -1 and $tr_comp{$db_tr_id}->{'exon'} == 0)){
                                $add_transcript = 0;
                                last DBTR;
                            }
                        }
                    }
                    #Aggressive mode: existing transcripts can be modified by merging them with new models if they are compatible.
                    #  TO DO: Also, allows replacing existing transcript models with new ones if they have the same stable ids.
                    elsif ($submode eq "aggr"){
                        #Merge transcript if there is exon novelty w.r.t. an existing transcript and both can be merged
                        #Add transcript if intron novelty, or if exon novelty but transcripts can't be merged
                        #Action may differ for different existing transcripts: prioritise merging?  with the one with the most intronic/exonic overlap?
                        DBTR:foreach my $db_tr_id (keys %tr_comp){
                            if ($tr_comp{$db_tr_id}->{'exon'} == 1){
                                if ($tr_comp{$db_tr_id}->{'merge'} == 1){
                                    my ($db_tr) = grep {$_->stable_id and $_->stable_id eq $db_tr_id} @{$db_gene->get_all_Transcripts};
                                    #Commented out as this will fetch the transcript from the database, ignoring modifications done in this session and not stored yet
                                    #my $db_tr = $dba->get_TranscriptAdaptor->fetch_by_stable_id($db_tr_id);
                                    unless ($db_tr){ #if not found, probably because it doesn't have a stable id yet, fetch it using the transcript name
                                        $db_tr = fetch_from_gene_by_transcript_name($db_tr_id, $db_gene);
                                    }
                                    #
                                    #Do not merge with MANE transcripts
                                    if (scalar grep {$_->value eq "MANE_select" or $_->value eq "MANE_plus_clinical"} @{$db_tr->get_all_Attributes('remark')}){
                                        print "Skipping transcript as partially redundant with a MANE transcript\n";
                                        $add_transcript = 0;
                                        @merge_candidates = ();
                                        $outcome = "partially redundant with a MANE transcript";
                                        last DBTR;
                                    }
                                    #Do not merge with full-length CDS transcripts (TO DO)
                                    #unless ($db_tr->translate and has_complete_cds($db_tr)){ 
                                    #Do not merge with coding transcripts
                                    #Skip transcript if it could be merged with a coding transcript,
                                    # so as not to make a partially redundant transcript
                                    elsif ($db_tr->translate()){
                                        print "Skipping transcript as partially redundant with a coding transcript\n";
                                        $add_transcript = 0;
                                        @merge_candidates = ();
                                        $outcome = "partially redundant with a coding transcript";
                                        last DBTR;
                                    }
                                    #Do not merge with Platinum transcripts unless this transcript is also Platinum
                                    elsif (scalar grep {$_->value eq "Platinum set"} @{$db_tr->get_all_Attributes('hidden_remark')}){
                                        if ($platinum){
                                            push(@merge_candidates, $db_tr);
                                        }
                                        else{
                                            print "Skipping transcript as partially redundant with a Platinum transcript\n";
                                            $add_transcript = 0;
                                            @merge_candidates = ();
                                            $outcome = "partially redundant with a Platinum transcript";
                                            last DBTR;
                                        }
                                    }
                                    #Do not merge with pseudogene transcripts
                                    elsif ($db_tr->biotype =~ /pseudogene/){
                                        print "Not merging with a pseudogene transcript\n";
                                        next DBTR;
                                    }
                                    else{
                                        push(@merge_candidates, $db_tr);
                                    }
                                }
                                else{
                                    #check if start/end extends internal exon of another transcript
                                    #clip start/end if likely misalignment: 
                                    #ref: ########-----#######-----####
                                    #tr:    ######-----########
                                    $add_transcript = 1;
                                }
                            }
                            elsif ($tr_comp{$db_tr_id}->{'exon'} == 0 and $tr_comp{$db_tr_id}->{'intron'} == 1){ #eg. exon skipping
                                $add_transcript = 1;
                            }
                            elsif ($tr_comp{$db_tr_id}->{'exon'} == 0 and ($tr_comp{$db_tr_id}->{'intron'} == 0 or $tr_comp{$db_tr_id}->{'intron'} == -1)){
                                $add_transcript = 0;
                                @merge_candidates = ();
                                last DBTR;
                            }
                        }
                    }
                }
                #Single-exon transcripts: do not merge unless it includes and extends a single-exon non-coding transcript; else, add unless it can be merged
                else{
                    DBTR:foreach my $db_tr_id (keys %tr_comp){
                        if ($tr_comp{$db_tr_id}->{'exon'} == 1 and $tr_comp{$db_tr_id}->{'merge'} == 1){
                            my $db_tr = $dba->get_TranscriptAdaptor->fetch_by_stable_id($db_tr_id);
                            unless ($db_tr){ #if not found, probably because it doesn't have a stable id yet, fetch it using the transcript name
                                $db_tr = fetch_from_gene_by_transcript_name($db_tr_id, $db_gene);
                                print "AAA2=".$db_tr->seq_region_start."-".$db_tr->seq_region_end."\n";
                                print "BBB2=".$tr->seq_region_start."-".$tr->seq_region_end."\n";
                            }
                            if (scalar(@{$db_tr->get_all_Exons}) == 1 and !($db_tr->translate())){ 
                                #Do not merge with MANE transcripts
                                if (scalar grep {$_->value eq "MANE_select" or $_->value eq "MANE_plus_clinical"} @{$db_tr->get_all_Attributes('remark')}){
                                  print "Skipping transcript as partially redundant with a MANE transcript\n";
                                  $add_transcript = 0;
                                  @merge_candidates = ();
                                  $outcome = "partially redundant with a MANE transcript";
                                  last DBTR;
                                }                                                         
                                #Do not merge with Platinum transcripts unless this transcript is also Platinum
                                elsif (scalar grep {$_->value eq "Platinum set"} @{$db_tr->get_all_Attributes('hidden_remark')}){
                                    if ($platinum){
                                        push(@merge_candidates, $db_tr);
                                    }
                                    else{
                                        print "Skipping transcript as partially redundant with a Platinum transcript\n";
                                        $add_transcript = 0;
                                        @merge_candidates = ();
                                        $outcome = "partially redundant with a Platinum transcript";
                                        last DBTR;
                                    }
                                }
                                #Do not merge with pseudogene transcripts
                                elsif ($db_tr->biotype =~ /pseudogene/){
                                    print "Not merging with a pseudogene transcript\n";
                                    next DBTR;
                                }
                                #elsif ($tr->start <= $db_tr->start and $tr->end => $db_tr->end and
                                #    !($tr->start == $db_tr->start and $tr->end == $db_tr->end)){
                                else{
                                    push(@merge_candidates, $db_tr);
                                }
                            }
                        }
                        #Add as new transcript if there is exon novelty with a transcript but can't be merged (eg. retained intron)
                        elsif ($tr_comp{$db_tr_id}->{'exon'} == 1 and $tr_comp{$db_tr_id}->{'merge'} == 0){
                            $add_transcript = 1;
                        }
                        #Reject if no exon novelty with a transcript
                        elsif ($tr_comp{$db_tr_id}->{'exon'} == 0){ 
                            $add_transcript = 0;
                            @merge_candidates = ();
                            $outcome = "no exon novelty";
                            last DBTR;
                        }
                    }

#                     #Merge with current single-exon transcript/gene if it encompasses it and extends it
#                     #Skip if the current transcript is coding!
#                     #This rule needs to be refined: single-exon transcripts in multi-transcript genes are ignored now
#                     if (scalar keys %tr_comp == 1){
#                         my ($db_tr_id) = keys %tr_comp;
#                         my $db_tr = $dba->get_TranscriptAdaptor->fetch_by_stable_id($db_tr_id);
#                         if (scalar(@{$db_tr->get_all_Exons}) == 1 and !($db_tr->translate())){
#                             if ($tr->start <= $db_tr->start and $tr->end => $db_tr->end and
#                                 !($tr->start == $db_tr->start and $tr->end == $db_tr->end)){
#                                 push(@merge_candidates, $db_tr);
#                             }
#                         }
#                     }
#                     #If there is exon novelty and can't be merged with any existing transcript, add new transcript
#                     #We are ignoring transcripts that could extend 5' or 3' ends for the moment
#                     else{
#                         DBTR:foreach my $db_tr_id (keys %tr_comp){
#                             if ($tr_comp{$db_tr_id}->{'exon'} == 1 and $tr_comp{$db_tr_id}->{'merge'} == 0){
#                                 $add_transcript = 1;
#                             }
#                             else{
#                                 $add_transcript = 0;
#                                 last DBTR;
#                             }
#                         }
#                     }
                }
              
                #Take action
                #my $id = $tr->get_all_Attributes('hidden_remark')->[0]->value;
                if (scalar @merge_candidates > 0){
                    my $sel_db_tr;
                    if (scalar @merge_candidates > 1){
                        $sel_db_tr = pick_db_tr(\@merge_candidates, $tr);
                    }
                    else{
                        $sel_db_tr = $merge_candidates[0];
                    }
                    print "TR=".$tr->seq_region_start."-".$tr->seq_region_end.":".$tr->seq_region_strand."  ".
                          join(", ", map {$_->seq_region_start."-".$_->seq_region_end.":".$_->seq_region_strand} @{$tr->get_all_Exons})."  ".
                          join(", ", map {$_->start."-".$_->end.":".$_->strand} @{$tr->get_all_Exons}).
                          "\n";
                    print "DB_TR=".$sel_db_tr->seq_region_start."-".$sel_db_tr->seq_region_end.":".$sel_db_tr->seq_region_strand."  ".
                          join(", ", map {$_->seq_region_start."-".$_->seq_region_end.":".$_->seq_region_strand} @{$sel_db_tr->get_all_Exons})."  ".
                          join(", ", map {$_->start."-".$_->end.":".$_->strand} @{$sel_db_tr->get_all_Exons}).
                          "\n";
                    #If the db_tr was not stored yet (ie. no ENS stable id), its exons will still be on local region slices and the coordinates of the merged transcript will get messed up
                    #Transform the db_tr to toplevel before merging it with the new transcript
                #    unless ($sel_db_tr->stable_id =~ /^ENST/){
                        $sel_db_tr = $sel_db_tr->transform("toplevel");
                        print "SLICES=".$sel_db_tr->slice->name."  ".join(", ", map {$_->slice->name} @{$sel_db_tr->get_all_Exons})."\n";
                #    }
                    #Merge transcripts
                    my $merged_transcript = merge_transcripts($tr, $sel_db_tr );
                    print "CCC=".$merged_transcript->seq_region_start."-".$merged_transcript->seq_region_end.":".$merged_transcript->seq_region_strand."  ".
                          join(", ", map {$_->seq_region_start."-".$_->seq_region_end.":".$_->seq_region_strand} @{$merged_transcript->get_all_Exons})."  ".
                          join(", ", map {$_->start."-".$_->end.":".$_->strand} @{$merged_transcript->get_all_Exons}).
                          "\n";
#     $db_gene->flush_Transcripts;
#                     #$db_gene->remove_Transcript($sel_db_tr);
#     $db_gene->add_Transcript($merged_transcript);    
                    #$slice_offset = $region->slice->start - 1;
#     foreach my $tr2 (@{$db_gene->get_all_Transcripts}){
                    foreach my $exon (@{$merged_transcript->get_all_Exons}){
                        $exon->start($exon->start - $slice_offset);
                        $exon->end(  $exon->end   - $slice_offset);
                        $exon->slice($region->slice);
                    }
#     }
                    print "DDD=".$merged_transcript->seq_region_start."-".$merged_transcript->seq_region_end."  ".
                          join(", ", map {$_->seq_region_start."-".$_->seq_region_end.":".$_->seq_region_strand} @{$merged_transcript->get_all_Exons})."  ".
                          join(", ", map {$_->start."-".$_->end.":".$_->strand} @{$merged_transcript->get_all_Exons}).
                          "\n";
                    foreach my $g ($region->genes) {
                        foreach my $ts (@{$g->get_all_Transcripts}){
                            #print "ABCDE ".$ts->stable_id." ".$ts->get_all_Attributes('name')->[0]->value."\n";
                            if (($merged_transcript->stable_id =~ /^ENS/ and $ts->stable_id eq $merged_transcript->stable_id)
                                  or ($merged_transcript->get_all_Attributes('name')->[0] and $ts->get_all_Attributes('name')->[0] and $merged_transcript->get_all_Attributes('name')->[0]->value eq $ts->get_all_Attributes('name')->[0]->value)
                                  or $ts == $merged_transcript){
                                print "FGHIJ S1=".$ts->stable_id." N1=".$ts->get_all_Attributes('name')->[0]->value." S2=".$merged_transcript->stable_id." N2=".$merged_transcript->get_all_Attributes('name')->[0]->value."\n";
                                print "TS=".$ts->seq_region_start."-".$ts->seq_region_end.":".$ts->seq_region_strand."  ".join(", ", map {$_->seq_region_start."-".$_->seq_region_end.":".$_->seq_region_strand} @{$ts->get_all_Exons})."\n";
                                print "MT=".$merged_transcript->seq_region_start."-".$merged_transcript->seq_region_end.":".$merged_transcript->seq_region_strand."  ".join(", ", map {$_->seq_region_start."-".$_->seq_region_end.":".$_->seq_region_strand} @{$merged_transcript->get_all_Exons})."\n";
                                # $ts and $merged_transcript point to the same object, so flushing exons in one does the same to the other
                                unless ($ts == $merged_transcript){
                                    $ts->flush_Exons; #This empties the exon list. If the same exons are added again, 
                                                      #those that are not associated with other transcripts will get a new stable id
                                    foreach my $exon (@{$merged_transcript->get_all_Exons}){
                                        $ts->add_Exon($exon); #print "ADD_EXON\n";
                                    }
                                }
                                $ts->recalculate_coordinates;
                                print "TS2=".$ts->seq_region_start."-".$ts->seq_region_end.":".$ts->seq_region_strand."  ".join(", ", map {$_->seq_region_start."-".$_->seq_region_end.":".$_->seq_region_strand} @{$ts->get_all_Exons})."\n";                       
                                print "ATTRIBS=".scalar(@{$merged_transcript->get_all_Attributes})."  ".scalar(@{$ts->get_all_Attributes})."\n";
                                my $ac = 0;
                                unless ($ts == $merged_transcript){
                                    #Add remarks and hidden remarks from the novel transcript
                                    foreach my $att (@{$merged_transcript->get_all_Attributes}){
                                        print "AAA1 ".++$ac."  ".$att->code." ".$att->value."\n";
                                        print "ATTRIBS=".scalar(@{$merged_transcript->get_all_Attributes})."  ".scalar(@{$ts->get_all_Attributes})."\n";
                                        #Do not add the 'not for VEGA' remark if the no_NFV flag is on
                                        if ($no_NFV){
                                            next if ($att->code eq "remark" and $att->value eq "not for VEGA");
                                        }
                                        #Append read names to existing remark
                                        if ((($att->code eq "hidden_remark" and $att->value =~ /^pacbio_capture_seq_\w+ : .+;/) or
                                            ($att->code eq "hidden_remark" and $att->value =~ /^SLR-seq_\w+ : .+;/) or
                                            ($att->code eq "hidden_remark" and $att->value =~ /^pacbio_raceseq_\w+ : .+;/) or
                                            ($att->code eq "hidden_remark" and $att->value =~ /^pacBio(SII):Cshl:\S+ : .+;/)
                                            ) and
                                            (scalar(@{$ts->get_all_Attributes('TAGENE_transcript')}) > 0)){
                                            my $appended = 0;
                                            print "AAA2\n";
                                            foreach my $att2 (@{$ts->get_all_Attributes('hidden_remark')}){
                                                if ($att2->value =~ /^pacbio_capture_seq_\w+ : .+;/ or
                                                    $att2->value =~ /^SLR-seq_\w+ : .+;/ or
                                                    $att2->value =~ /^pacbio_raceseq_\w+ : .+;/ or
                                                    $att2->value =~ /^pacBio(SII):Cshl:\S+ : .+;/){
                                                    $att2->value($att2->value." ".$att->value);
                                                    $appended = 1;
                                                    last;
                                                }
                                            }
                                            #If no read name remark yet, create one
                                            unless ($appended){
                                                $ts->add_Attributes($att);
                                                print "AAA3\n";
                                            }
                                        }
                                        else{
                                            $ts->add_Attributes($att);
                                            print "AAA4\n";
                                        }
                                    }
                                }
                                print "KLMNO\n";                               
                                #Add remark that indicates that the transcript was extended
                               # unless (scalar(grep {$_->value eq "TAGENE_extended"} @{$ts->get_all_Attributes('hidden_remark')})){
                               #     $ts->add_Attributes( Bio::EnsEMBL::Attribute->new(-code => 'hidden_remark', -value => 'TAGENE_extended') );
                               # }

                                
                                #If flag is on, force comp_pipe biotype and "not for VEGA" remark (they should normally be associated)
                                if ($force_cp_biotype){
                                    $ts->biotype($CP_BIOTYPE);
                                    unless (scalar(grep {$_->value eq "not for VEGA"} @{$ts->get_all_Attributes('remark')})){
                                        $ts->add_Attributes( Bio::EnsEMBL::Attribute->new(-code => 'remark', -value => 'not for VEGA') );
                                    }
                                }
                                else{
                                    #If 'comp_pipe' transcript, change its biotype and remove the 'not for VEGA' attribute
                                    if ($ts->biotype eq $CP_BIOTYPE){
                                        if ($ts->translation){ #CDS from the gtf/gff3 file
                                            $ts->biotype("protein_coding"); #should look for NMD too!
                                        }
                                        else{
                                            my $t_biotype = get_transcript_biotype($db_gene);
                                            $ts->biotype($t_biotype);
                                        }
                                        if ($no_NFV){
                                            if (scalar(grep {$_->value eq "not for VEGA"} @{$ts->get_all_Attributes('remark')})){
                                                my $array = $ts->{'attributes'};
                                                @$array = grep { $_->value ne "not for VEGA" } @$array;
                                            }
                                        }
                                    }
                                }
                                
                                #Add 'Platinum set' hidden remark if indicated
                                if ($platinum){
                                    unless (scalar(grep {$_->value eq "Platinum set"} @{$ts->get_all_Attributes('hidden_remark')})){
                                        $ts->add_Attributes( Bio::EnsEMBL::Attribute->new(-code => 'hidden_remark', -value => 'Platinum set') );
                                    }
                                }
                                 

                                #Firstly, check for intron retention
                                #Else, if allowed, try to assign a CDS and change the biotype accordingly
                                if (LoutreWrite::CDSCreation::is_retained_intron($ts, $g, 5)){
                                    $ts->biotype("retained_intron"); print "STRAND=".$ts->seq_region_strand."\n";
                                }
                                else{
                                    unless ($do_not_add_cds){
                                        assign_cds_to_transcripts($ts, $g, $slice_offset);
                                    }
                                }
                                
                                #Change authorship
                                $ts->transcript_author($merged_transcript->transcript_author);
                                
                                print "TR: $id: Will extend transcript ".($ts->stable_id || $ts->get_all_Attributes("name")->[0]->value)." (".$ts->biotype.") in gene ".$g->stable_id."\n";
                                push (@log, "TR2: $id: Extended transcript ".($ts->stable_id || $ts->get_all_Attributes("name")->[0]->value)." (".$ts->biotype.") in gene ".$g->stable_id." (".$g->biotype.")");
                            }
                        }
                    }
    
                    #Update transcript
                    #$sel_db_tr->flush_Exons;   #This empties the exon list. If the same exons are added again, 
                                                #those that are not associated with other transcripts will get a new stable id
                    #foreach my $exon (@{$merged_transcript->get_all_Exons}){
                    #    $sel_db_tr->add_Exon($exon);
                    #}
                    #$db_gene->add_Transcript($sel_db_tr);
              #      print "TR: $id: Will extend transcript ".$merged_transcript->stable_id." in gene ".$db_gene->stable_id."\n";
              #      push (@log, "TR2: $id: Extended transcript ".$merged_transcript->stable_id." (".$tr->biotype.") in gene ".$db_gene->stable_id." (".$db_gene->biotype.")");
                }
                elsif ($add_transcript == 1){
                    #Create a name for the new transcript
                    #Before that, add a name to the gene if it doesn't have one
                    if (scalar @{$db_gene->get_all_Attributes('name')} == 0){
                      my $new_gene_name = $self->get_new_gene_name($region, $db_gene, $dataset_name, $dba) or return " ";
                      my $gene_name_att = Bio::EnsEMBL::Attribute->new(-code => 'name', -value => $new_gene_name);
                      $db_gene->add_Attributes($gene_name_att);
                    }
                    my $new_tr_name = $self->get_new_transcript_name($db_gene, $dba);
                    my $name_att = Bio::EnsEMBL::Attribute->new(-code => 'name', -value => $new_tr_name);
                    $tr->add_Attributes($name_att);
                    print "TR0_START=".$tr->seq_region_start."; TR0_END=".$tr->seq_region_end."; EXONS=".join(":", map {$_->seq_region_start."-".$_->seq_region_end} @{$tr->get_all_Exons})."\n";
                    $tr->slice($region->slice);
                    $tr->start($tr->start - $slice_offset);
                    $tr->end($tr->end - $slice_offset);
                    #foreach my $exon (@{$tr->get_all_Exons}){
                    #    $exon->start($exon->start - $slice_offset);
                    #    $exon->end(  $exon->end   - $slice_offset);
                    #    $exon->slice($region->slice);
                    #}
                    if ($platinum){
                        $tr->add_Attributes( Bio::EnsEMBL::Attribute->new(-code => 'hidden_remark', -value => 'Platinum set') );
                    }

                    #unless forced comp_pipe biotype, assign transcript biotype based on other transcripts
                    unless ($force_cp_biotype){
                        print "XXX".$tr->biotype."\n";
                        if ($tr->translation){ #CDS from the gtf/gff3 file
                       #     $tr->biotype("protein_coding"); #should look for NMD too!
                        }
                        else{
                            my $t_biotype = get_transcript_biotype($db_gene);
                            $tr->biotype($t_biotype);
                        }                  
                        if ($no_NFV){
                            if (scalar(grep {$_->value eq "not for VEGA"} @{$tr->get_all_Attributes('remark')})){
                                my $array = $tr->{'attributes'};
                                @$array = grep { $_->value ne "not for VEGA" } @$array;
                            }
                        }
                    }

                    print "TR1_START=".$tr->seq_region_start."; TR1_END=".$tr->seq_region_end."; EXONS=".join(":", map {$_->seq_region_start."-".$_->seq_region_end} @{$tr->get_all_Exons})."\n";
                    #Firstly, check for intron retention
                    #Else, if allowed, try to assign a CDS and change the biotype accordingly
                    if (LoutreWrite::CDSCreation::is_retained_intron($tr, $db_gene, 5)){
                        $tr->biotype("retained_intron");
                    }
                    else{
                        unless ($do_not_add_cds){
                            assign_cds_to_transcripts($tr, $db_gene, $slice_offset);
                        }
                    }
                    $db_gene->add_Transcript($tr);
                    print "TR: $id: Will add transcript $new_tr_name (".$tr->biotype.") to gene ".$db_gene->stable_id."\n";
                    push (@log, "TR2: $id: Added transcript $new_tr_name (".$tr->biotype.") to gene ".$db_gene->stable_id." (".$db_gene->biotype.")");
                }
                else{
                    if ($outcome){
                        print "TR: $id: Will reject transcript in host gene ".$db_gene->stable_id." as $outcome\n";
                        push (@log, "TR2: $id: Rejected transcript in host gene ".$db_gene->stable_id." as $outcome");
                    }
                    else{
                        print "TR: $id: Will reject transcript in host gene ".$db_gene->stable_id." as intron chain exists\n";
                        push (@log, "TR2: $id: Rejected transcript in host gene ".$db_gene->stable_id." as intron chain exists");
                    }
                }
            }
        }
        else{
            warn "More than one gene with stable id ".$gene->stable_id."!!!\n";
        }
    }

    #Gene biotype
    #Check Bio::Vega::Gene->set_biotype_status_from_transcripts


    #Write region
    my $g_msg = " ";
    if ($WRITE){
        $local_server->authorized_user($gene->gene_author->name); # preserve authorship
        my $n = 0;
        while (($g_msg eq " " or $g_msg =~ /write failed/ or $g_msg =~ /lock failed/) and $n<5){
            $g_msg = write_gene_region($region_action, $region);
            sleep(30);
            $n++;
        }
        print $g_msg."\n";
    }

    return ($g_msg, \@log);
}




=head2 get_new_gene_name

 Arg[1]    : Bio::Vega::Region
 Arg[2]    : Bio::Vega::Gene object
 Arg[3]    : Otter dataset name, eg. 'human', 'mouse', 'human_test', etc
 Arg[4]    : Bio::Vega::DBSQL::DBAdaptor
 Function  : Get clone-based name for a new gene
 Returntype: String

=cut

sub get_new_gene_name {
    my ($self, $region, $gene, $dataset_name, $dba) = @_;

    #For new databases that only have a 'primary_assembly' coordinate system name, return a genomic coordinate-based gene name
    if ($region->slice->coord_system_name eq "primary_assembly"){
      return join("_", $gene->seq_region_name, $gene->seq_region_start, $gene->seq_region_end);
    }

    #Older databases should have a clone coordinate system name. The clone-based gene name is obtained as described below.
    
    #In Otter/Zmap, the gene name is created by selecting 'New' (Ctrl+N) in the session window
    #The code can be found in ensembl-otter/modules/MenuCanvasWindow/SessionWindow.pm (subs_edit_new_subsequence and _region_name_and_next_locus_number)
    #Thus, the name has been created when the gene attributes are dumped in the XML and SQLite files

    #Find the clone that overlaps the gene's 3' end - according to ensembl-otter/modules/MenuCanvasWindow/SessionWindow.pm (sub _region_name_and_next_locus_number)
    #my $most_3prime = $gene->strand == 1 ? $gene->end : $gene->start;
    my $most_3prime = $gene->seq_region_strand == 1 ? $gene->seq_region_end : $gene->seq_region_start; #try genome coordinates as local coordinates may fail in some cases
    
#    my $clone_name;
#    my $tmp_slice = $dba->get_SliceAdaptor->fetch_by_region("toplevel", $gene->seq_region_name, $most_3prime, $most_3prime);
#    if ($tmp_slice->project("contig")->[0]){
#        $clone_name = $tmp_slice->project("contig")->[0]->to_Slice->seq_region_name;
#    }
    
    my $clone_name;
    my $actual_cs;
    foreach my $cs ($region->fetch_CloneSequences){ #Bio::Vega::Region method - gets clones trimmed down to the region boundaries
      #print "AAAAAAAAAAAAA  ".$cs->chr_start." ".$cs->chr_end."  ".$most_3prime."  ".$gene->start."-".$gene->end."  ".$gene->seq_region_start."-".$gene->seq_region_end."\n";
        if ($cs->chr_start <= $most_3prime and $cs->chr_end >= $most_3prime){
            $clone_name = $cs->clone_name;
            $actual_cs = $cs;
        }
    }
    if (!$clone_name){
        warn "No clone for gene on ".$gene->seq_region_name." ".$gene->start." ".$gene->end."\n";
        return undef;
    }
    #Trim sequence version from clone name
    $clone_name =~ s/\.\d+$//;

    #Find unused gene name - NOTE: this looks in the gene region only (subsequence of the clone?) - should it fetch all gene names in the clone?
    my $max = 0;
    
    #foreach my $gene ($region->genes){
    #    if (scalar(@{$gene->get_all_Attributes('name')})){
    #        my ($n) = $gene->get_all_Attributes('name')->[0]->value =~ /$clone_name\.(\d+)/;
    #        if ($n && $n > $max) {
    #            $max = $n;
    #        }
    #    }
    #}
    
    #Get Otter region for the whole clone
#    my $local_server = Bio::Otter::Server::Support::Local->new( otter_dba => $dba );
#    $local_server->set_params(
#        dataset => $dataset_name,
#        chr     => "chr".$actual_cs->chromosome."-38",
#        start   => $actual_cs->chr_start,
#        end     => $actual_cs->chr_end,
#        cs      => "chromosome",
#        csver   => "Otter"
#    );
print "\nZZZ ".$actual_cs->chromosome." ".$actual_cs->chr_start." ".$actual_cs->chr_end."  ".$actual_cs->clone_name." ".$actual_cs->contig_name." ".$actual_cs->accession."  ZZZ\n";

#    my $region_action = Bio::Otter::ServerAction::Script::Region->new_with_slice($local_server);
#    my $clone_region = $region_action->get_region; #Bio::Vega::Region
#    foreach my $gene ($clone_region->genes){
#        if (scalar(@{$gene->get_all_Attributes('name')})){
#            my ($n) = $gene->get_all_Attributes('name')->[0]->value =~ /$clone_name\.(\d+)/;
#            if ($n && $n > $max) {
#                $max = $n;
#            }
#        }
#    }    

    #my $contig_slice = $dba->get_SliceAdaptor->fetch_by_region("toplevel", $gene->seq_region_name, $actual_cs->chr_start, $actual_cs->chr_end);
    my $contig_slice = $dba->get_SliceAdaptor->fetch_by_region("toplevel", $actual_cs->contig_name);    #Error: MSG: Cannot deal with mis-lengthed mappings so far
    #A contig slice can be projected to more than one discontiguous chromosome segments, eg. "U62317.2.1.139934" in chr22-38
    foreach my $proj (@{$contig_slice->project('chromosome')}){ 
        my $chr_slice = $proj->to_Slice;
        #Genes can have a gene-symbol-based name whereas their transcripts have clone-based names
        #Then, if we look at the gene names only, we might come up with a gene name that has been used previously in other locus,
        #and the transcripts of the new gene might be given names already in use in the other locus (this will make otter crash)
        foreach my $gene (@{$chr_slice->get_all_Genes}){
            print "AAA ".$gene->stable_id."\n\n"; #."  ".$gene->get_all_Attributes('name')->[0]->value."\n\n";
            foreach my $transcript (@{$gene->get_all_Transcripts}){
                if (scalar @{$transcript->get_all_Attributes('name')}){
                    #$clone_name =~ s/\.\d+$//;
                    my ($n, $tmp) = $transcript->get_all_Attributes('name')->[0]->value =~ /$clone_name\.(\d+)-(\d{3,})$/;
                    if ($n and $n > $max) {
                        $max = $n;
                    }
                }
            }
        }
    }

    my $new_gene_name;
    my $name_in_db;
    do {
        $new_gene_name = $clone_name.".".(++$max);
        #Sanity check - avoid duplication with past or present transcripts
        my $dbcon = $dba->dbc;
        my $sth = $dbcon->prepare("SELECT t.stable_id FROM transcript t, transcript_attrib ta 
                                    WHERE t.transcript_id=ta.transcript_id AND ta.value LIKE \"$new_gene_name-\%\"");
        $sth->execute();
        my @tsi = $sth->fetchrow_array;
        $sth->finish();
        $name_in_db = scalar @tsi;
    } while ($name_in_db);
    
    print "New gene name: $new_gene_name\n";
    return $new_gene_name;
}



=head2 get_new_transcript_name

 Arg[1]    : Bio::Vega::Gene object
 Arg[2]    : Bio::Vega::DBSQL::DBAdaptor
 Function  : Get name for a new transcript - it should be an extension of a clone-based locus name
 Returntype: String

=cut

sub get_new_transcript_name {
    my ($self, $gene, $dba) = @_;
    
    #In Otter/Zmap, the transcript name is created by selecting 'Variant' (Ctrl+I) in the session window
    #The code can be found in ensembl-otter/modules/MenuCanvasWindow/SessionWindow.pm (sub _make_variant_subsequence)
    #Thus, the name has been created when the transcript attributes are dumped in the XML and SQLite files
#    while (scalar @{$gene->get_all_Attributes('name')} == 0){ print "WAITING FOR 'name'...\n"; sleep(60); }
    my $gene_name = $gene->get_all_Attributes('name')->[0]->value; 
    #The gene name attribute may be a gene symbol rather than a clone-based name
    #Use it only as a default name (eg. for new genes)
    #To imitate what otter does, get a list of existing (clone-based) transcript names 
    #in order to infer the clone-based gene name and the next transcript number
    my $max = 0;
    foreach my $tr (@{$gene->get_all_Transcripts}){
        if (scalar @{$tr->get_all_Attributes('name')}){
            my ($root, $n) = $tr->get_all_Attributes('name')->[0]->value =~ /^(.+)-(\d{3,})$/;
            if ($n and $n > $max){
                $max = $n;
            }
            if ($root){
                $gene_name = $root;
            }
        }
    }
    
    my $new_tr_name;
    my $name_in_db;
    do {
        $new_tr_name = sprintf("%s-%03d", $gene_name, ++$max);
        #Sanity check - avoid duplication with past or present transcripts
        my $dbcon = $dba->dbc;
        my $sth = $dbcon->prepare("SELECT t.stable_id FROM transcript t, transcript_attrib ta 
                                    WHERE t.transcript_id=ta.transcript_id AND ta.value=\"$new_tr_name\"");
        $sth->execute();
        my @tsi = $sth->fetchrow_array;
        $sth->finish();
        $name_in_db = scalar @tsi;
    } while ($name_in_db);    

    return $new_tr_name;
}




=head2 intron novelty

 Arg[1]    : Bio::Vega::Transcript object (novel transcript)
 Arg[2]    : Bio::Vega::Transcript object (database transcript)
 Function  : Look for novel splice sites in the first transcript. Returns '1' if there is intron novelty, '0' if intron chains are identical, and '-1' if the intron chain of the first transcript is included in the intron chain of the second one and the latter has additional introns
 Returntype: Integer 

=cut

sub intron_novelty {
    my ($tr, $db_tr) = @_;
    INT:foreach my $intron (@{$tr->get_all_Introns}){
        my $seen = 0;
        foreach my $db_intron (@{$db_tr->get_all_Introns}){
            if ($intron->seq_region_start == $db_intron->seq_region_start and $intron->seq_region_end == $db_intron->seq_region_end){
                $seen = 1;
                next INT;
            }
        }
        #Novel intron (may overlap existing intron but splice sites are different)
        if ($seen == 0){
            return 1;
        }
    }
    #Even if no intron novelty, report if novel transcript's intron chain is included in database transcript's intron chain
    if (scalar(@{$tr->get_all_Introns}) < scalar(@{$db_tr->get_all_Introns})){
      return -1;
    }
    return 0;
}



=head2 exon novelty

 Arg[1]    : Bio::Vega::Transcript object (novel transcript)
 Arg[2]    : Bio::Vega::Transcript object (database transcript)
 Function  : Look for exon sequence novelty in the first transcript.
 Returntype: Boolean

=cut

sub exon_novelty {
    my ($tr, $db_tr) = @_;
    EX:foreach my $exon (@{$tr->get_all_Exons}){
        my $seen = 0;
        my $ovlp = 0;
        foreach my $db_exon (@{$db_tr->get_all_Exons}){
            if ($exon->seq_region_start <= $db_exon->seq_region_end and $exon->seq_region_end >= $db_exon->seq_region_start){
                $ovlp = 1;
                #if exons overlap but the novel exon extends past 5' or 3' end of the database exon
                if ($exon->seq_region_start < $db_exon->seq_region_start or $exon->seq_region_end > $db_exon->seq_region_end){
                    return 1;
                }
                else{
                    next EX;
                }
            }
        }
        #if exon doesn't overlap any exon of the database transcript 
        if ($ovlp == 0){
            return 1;
        }
    }
    return 0;
}



=head2 can_be_merged

 Arg[1]    : Bio::Vega::Transcript object (novel transcript)
 Arg[2]    : Bio::Vega::Transcript object (database transcript)
 Function  : Returns true if transcripts can be merged, i.e. they have exon-exon overlap but no exon-intron overlap.
             If both transcripts are multi-exonic, they must share at least an intron.
 Returntype: Boolean

=cut

sub can_be_merged {
  my ($tr, $db_tr) = @_;
    
  #Exon-exon overlap: required
  my $exon_ovlp = 0;
  EXON:foreach my $exon (@{$tr->get_all_Exons}){
    foreach my $db_exon (@{$db_tr->get_all_Exons}){
      if ($exon->seq_region_start <= $db_exon->seq_region_end and $exon->seq_region_end >= $db_exon->seq_region_start){
        $exon_ovlp = 1;
        last EXON;
      }
    }
  }
  unless ($exon_ovlp){
    return 0;
  }
    
  #Exon-intron overlap: not allowed
  foreach my $intron (@{$tr->get_all_Introns}){
    foreach my $db_exon (@{$db_tr->get_all_Exons}){
      if ($intron->seq_region_start <= $db_exon->seq_region_end and $intron->seq_region_end >= $db_exon->seq_region_start){
        return 0;
      }
    }
  }
  foreach my $db_intron (@{$db_tr->get_all_Introns}){
    foreach my $exon (@{$tr->get_all_Exons}){
      if ($db_intron->seq_region_start <= $exon->seq_region_end and $db_intron->seq_region_end >= $exon->seq_region_start){
        return 0;
      }
    }
  }
  
  #Intron match required if both are multi-exonic transcripts
  my $intron_match = 0;
  INTRON:foreach my $intron (@{$tr->get_all_Introns}){
    foreach my $db_intron (@{$db_tr->get_all_Introns}){
      if ($intron->seq_region_start == $db_intron->seq_region_start and $intron->seq_region_end >= $db_intron->seq_region_end){
        $intron_match = 1;
        last INTRON;
      }
    }
  }

  if ($exon_ovlp == 1){
    if (scalar @{$tr->get_all_Introns} and scalar @{$db_tr->get_all_Introns}){
      if ($intron_match){
        return 1;
      }
      else{
        return 0;
      }
    }
    else{
      return 1;
    }
  }
}



=head2 pick_db_tr

 Arg[1]    : Array-ref of Bio::Vega::Transcript objects (database transcripts that can be merged with the novel transcript)
 Arg[2]    : Bio::Vega::Transcript object (novel transcript)
 Function  : Returns the most suitable transcript to be merged with the novel transcript, based on exon-exon overlap and number of shared introns
 Returntype: Bio::Vega::Transcript

=cut

sub pick_db_tr {
    my ($merge_c, $tr) = @_;
    my $selected_db_tr;
    my @candidates;
    my @best_candidates;
    my $max_n_introns = 0;
    #Prioritise number of shared introns
    foreach my $db_tr (@$merge_c){
        my $shared_intron_count = 0;
        I:foreach my $intron (@{$tr->get_all_Introns}){
            foreach my $db_intron (@{$db_tr->get_all_Introns}){
                if ($intron->seq_region_start == $db_intron->seq_region_start and $intron->seq_region_end == $db_intron->seq_region_end){
                    $shared_intron_count++;
                    next I;
                }
            }
        }
        if ($shared_intron_count > $max_n_introns){
            @best_candidates = ();
            push(@best_candidates, $db_tr);
            $max_n_introns = $shared_intron_count;
        }
        elsif ($shared_intron_count == $max_n_introns){
            push(@best_candidates, $db_tr);
        }
    }
    #If only one, pick it
    if (scalar(@best_candidates) == 1){
        $selected_db_tr = $best_candidates[0];
    }
    #Else, check exon overlap - find transcript with the most exonic overlap
    else{
        if (scalar(@best_candidates) == 0){
            @best_candidates = @$merge_c;
        }
        my %tr_loc; #Store genomic positions of transcript's exons
        foreach my $exon (@{$tr->get_all_Exons}){
          print "TR EXON: ".$exon->seq_region_start."-".$exon->seq_region_end."\n";
            for (my $i = $exon->seq_region_start; $i <= $exon->seq_region_end; $i++){
                $tr_loc{$i} = 1;
            }
        }
        my $max = 0;
        foreach my $db_tr (@best_candidates){
            my %seen; #Store genomic positions where novel transcript's exons overlap database transcript's exons
            foreach my $db_exon (@{$db_tr->get_all_Exons}){
              print "DB_TR EXON: ".$db_exon->seq_region_start."-".$db_exon->seq_region_end."\n";
                for (my $i = $db_exon->seq_region_start; $i <= $db_exon->seq_region_end; $i++){
                    if ($tr_loc{$i}){
                        $seen{$i} = 1;
                    }
                }
            }
            if (scalar(keys %seen) > $max){
                $selected_db_tr = $db_tr;
                $max = scalar(keys %seen);
            }
        }
    }
    return $selected_db_tr;
}



=head2 merge_transcripts

 Arg[1]    : Bio::Vega::Transcript object (novel transcript)
 Arg[2]    : Bio::Vega::Transcript object (database transcript)
 Function  : Returns a transcript resulting from merging the two input transcripts
 Returntype: Bio::Vega::Transcript

=cut

sub merge_transcripts {
    my ($tr, $db_tr) = @_;
    print "DB_TR=".$db_tr->seq_region_start."-".$db_tr->seq_region_end."   ".scalar(@{$db_tr->get_all_Attributes})."\n";
    print "TR=".$tr->seq_region_start."-".$tr->seq_region_end."   ".scalar(@{$tr->get_all_Attributes})."\n";
    print "BEFORE: ".join('+', map {$_->vega_hashkey} @{$db_tr->get_all_Exons})."\n";
    EX:foreach my $exon (@{$tr->get_all_Exons}){
        #Check if exon is completely new or overlaps an existing exon
        my $seen = 0;
        foreach my $db_exon (@{$db_tr->get_all_Exons}){
            if ($exon->seq_region_start == $db_exon->seq_region_start and $exon->seq_region_end == $db_exon->seq_region_end){
                $seen = 1;
                next EX;
            }
            elsif ($exon->seq_region_start <= $db_exon->seq_region_end and $exon->seq_region_end >= $db_exon->seq_region_start){
                $seen = 1;
                #If exons overlap, merge them, remove the old exon and add the merged exon to the database transcript.
                #NOTE: not checking if there is exon-intron overlap or if an exon overlaps two exons of the other transcript. 
                #It has been checked elsewhere.
                my $new_start = min($exon->seq_region_start, $db_exon->seq_region_start);
                my $new_end = max($exon->seq_region_end, $db_exon->seq_region_end);
                my $novel_exon = new Bio::Vega::Exon ( -start => $new_start, 
                                                       -end => $new_end,
                                                       -strand => $db_tr->seq_region_strand,
                                                       -slice => $db_tr->slice
                                                      );
                $novel_exon->phase(-1);
                $novel_exon->end_phase(-1);
                $db_tr->swap_exons($db_exon, $novel_exon);
                next EX;
            }
        }
        #If novel, make new Exon object and add it to the database transcript
        if ($seen == 0){
            my $novel_exon = new Bio::Vega::Exon ( -start => $exon->seq_region_start, 
                                                   -end => $exon->seq_region_end,
                                                   -strand => $db_tr->seq_region_strand,
                                                   -slice => $db_tr->slice
                                                  );
            $novel_exon->phase(-1);
            $novel_exon->end_phase(-1);
            $db_tr->add_Exon($novel_exon);
        }
    }
    print "AFTER: ".join('+', map {$_->vega_hashkey} @{$db_tr->get_all_Exons})."\n";

    #Add remarks (except 'not for VEGA') and hidden remarks from the novel transcript
    foreach my $att (@{$tr->get_all_Attributes}){
        unless (grep {$_->code eq $att->code and $_->value eq $att->value} @{$db_tr->get_all_Attributes}){
            $db_tr->add_Attributes($att);
        }
    }
    
    #Change authorship
    $db_tr->transcript_author($tr->transcript_author);
    print "AAA3=".$db_tr->seq_region_start."-".$db_tr->seq_region_end."  ".$db_tr->get_all_Exons->[0]->seq_region_start."-".$db_tr->get_all_Exons->[0]->seq_region_end."   ".scalar(@{$db_tr->get_all_Attributes})."\n";
    print join(", ", map {$_->code.": ".$_->value} @{$db_tr->get_all_Attributes})."\n";
    return $db_tr;
}



#Get transcript biotype from existing transcript biotypes
#Mainly for non-coding genes
sub get_transcript_biotype {
    my $db_gene = shift;
    my %biotypes;
    if ($db_gene->biotype =~ /^(processed_transcript|antisense|lincrna|lincRNA|sense_intronic|sense_overlapping|bidirectional_promoter_lncrna|tec|comp_pipe)$/){
        foreach my $db_tr (@{$db_gene->get_all_Transcripts}){
            if ($db_tr->biotype =~ /^(antisense|lincrna|lincRNA|sense_intronic|sense_overlapping|bidirectional_promoter_lncrna|3\'_overlapping_ncrna)$/){
                $biotypes{$db_tr->biotype}++;
            }
            #if ($biotypes{"antisense"} > 0 and !$biotypes{"lincrna"} and !$biotypes{"sense_intronic"} and !$biotypes{"sense_overlapping"} and !$biotypes{"bidirectional_promoter_lncrna"} and !$biotypes{"3'_overlapping_ncrna"}){
        }
        if (scalar keys %biotypes == 1){
            foreach (keys %biotypes){
                return $_;
            }
        }
    }

    return "processed_transcript";
}





=head2 is_retained_intron

 Arg[1]    : Bio::Vega::Transcript object
 Arg[2]    : Bio::Vega::Gene object
 Arg[3]    : integer (minimum overhang in partial intron retention)
 Function  : finds whether the novel transcript is a retained_intron model in the host gene
 Returntype: boolean

=cut

sub is_retained_intron {
  my ($transcript, $host_gene, $min_overhang) = @_;
  #i) Full intron retention: the model fully transcribes the sequence between the splice donor site and acceptor site of adjacent coding exons (i.e. a ‘CDS intron’) of a model in loutre. This could include scenarios where the TAGENE model retains multiple CDS introns of the same loutre model. 
  #ii) Partial intron retention: the model has a start or end coordinate within a CDS intron of a loutre model, but not (respectively) a splice donor site or splice acceptor site within that same intron. Do we need some kind of size cut-off here? I.e. for fuzzy edges that are misalignments rather than real retained introns. Should we clip these? 
  my $slice_offset = $host_gene->seq_region_start - 1; #print "OFFSET=$slice_offset\n";
  $min_overhang ||= 1;
  foreach my $tr (@{$host_gene->get_all_Transcripts}){
    if ($tr->translate){
      #print "TR:  ".$tr->stable_id." ".$tr->biotype."\n"; 
      my $cds_start = $tr->coding_region_start + $slice_offset;
      my $cds_end = $tr->coding_region_end + $slice_offset;
      if ($cds_start > $cds_end){
        ($cds_start, $cds_end) = ($cds_end, $cds_start);
      }
      #print "CDS: $cds_start-$cds_end\n";
      foreach my $intron (@{$tr->get_all_Introns}){
        #print "INTRON: ".$intron->seq_region_start."-".$intron->seq_region_end."\n";
        #Is CDS intron?
        if ($intron->seq_region_start > $cds_start and $intron->seq_region_end < $cds_end){
          #print "CDS INTRON: ".$intron->seq_region_start."-".$intron->seq_region_end."\n";
          #Full intron retention
          #  ref: #######-------#######------######
          #   tr:       #########
          foreach my $exon (@{$transcript->get_all_Exons}){
            #print "EXON: ".$exon->seq_region_start."-".$exon->seq_region_end."\n";
            if ($exon->seq_region_start < $intron->seq_region_start and $exon->seq_region_end > $intron->seq_region_end){
              #print "FOUND\n";
              return 1;
            }
          }
          #Partial intron retention
          #  ref: #######-------#######------######
          #   tr:          #######---####
          if ($transcript->seq_region_strand == 1 and 
              $transcript->start_Exon->seq_region_start >= $intron->seq_region_start and
              $transcript->start_Exon->seq_region_end > $intron->seq_region_end and
              ($transcript->start_Exon->seq_region_start + $min_overhang) <= $intron->seq_region_end){
                return 1;
          }
          if ($transcript->seq_region_strand == -1 and 
              $transcript->end_Exon->seq_region_start >= $intron->seq_region_start and
              $transcript->end_Exon->seq_region_end > $intron->seq_region_end and
              ($transcript->end_Exon->seq_region_start + $min_overhang) <= $intron->seq_region_end){
                #Check polyA site support
                unless (has_polyA_site_support($transcript, 500) or has_polyAseq_support($transcript, 500)){
                  print "No polyA site support - will make a retained_intron\n";
                  return 1;
                }
          }
          if ($transcript->seq_region_strand == 1 and 
              $transcript->end_Exon->seq_region_start < $intron->seq_region_start and
              $transcript->end_Exon->seq_region_end <= $intron->seq_region_end and
              ($transcript->end_Exon->seq_region_end - $min_overhang) >= $intron->seq_region_start){
                #Check polyA site support
                unless (has_polyA_site_support($transcript, 500) or has_polyAseq_support($transcript, 500)){
                  print "No polyA site support - will make a retained_intron\n";
                  return 1;
                }
          }
          if ($transcript->seq_region_strand == -1 and 
              $transcript->start_Exon->seq_region_start < $intron->seq_region_start and
              $transcript->start_Exon->seq_region_end <= $intron->seq_region_end and
              ($transcript->start_Exon->seq_region_end - $min_overhang) >= $intron->seq_region_start){
                return 1;
          }
        }
      }
    }
  }
  return 0;
}




=head2 has_polyA_site_support

 Arg[1]    : Bio::Vega::Transcript object
 Arg[2]    : integer (distance threshold)
 Function  : Returns true if there is a Havana-annotated polyA site within the distance to the transcript 3' end indicated by the second argument
 Returntype: none

=cut

sub has_polyA_site_support {
  my ($transcript, $threshold) = @_;
  my $transcript_end = $transcript->seq_region_strand == 1 ? $transcript->seq_region_end : $transcript->seq_region_start;
  my $sa = $DBA{'havana'}->get_SliceAdaptor();
  my $ext_slice = $sa->fetch_by_region("toplevel", $transcript->slice->seq_region_name, $transcript_end - $threshold, $transcript_end + $threshold);
  if (scalar grep {$_->seq_region_strand==$transcript->seq_region_strand} @{$ext_slice->get_all_SimpleFeatures('polya_site')}){
    return 1;
  }
  return 0;
}



=head2 has_polyAseq_support

 Arg[1]    : Bio::Vega::Transcript object
 Arg[2]    : integer (distance threshold)
 Function  : Returns true if there are polyA-seq features from at least four different tissues within the distance to the transcript 3' end indicated by the second argument
 Returntype: none

=cut

sub has_polyAseq_support {
  my ($transcript, $threshold) = @_;
  my $transcript_end = $transcript->seq_region_strand == 1 ? $transcript->seq_region_end : $transcript->seq_region_start;
  my $sa = $DBA{'core'}->get_SliceAdaptor();
  my $slice = $sa->fetch_by_region("toplevel", $transcript->seq_region_name, $transcript_end, $transcript_end);
  #Project slice to polyA-seq feature assembly version
  my $cs_version = $DBA{'polyAseq'}->get_CoordSystemAdaptor->get_default_version;
  my $projection = $slice->project("toplevel", $cs_version);
  if (scalar @$projection){
    my $segment = $projection->[0];
    my $psa = $DBA{'polyAseq'}->get_SliceAdaptor();    
    my $ext_slice = $psa->fetch_by_region("toplevel", $segment->to_Slice->seq_region_name, $segment->to_Slice->seq_region_end - $threshold, $segment->to_Slice->seq_region_end + $threshold); 
    my @polyAseq_features;
    if ($ext_slice){
      @polyAseq_features = grep {$_->seq_region_strand==$transcript->seq_region_strand} @{$ext_slice->get_all_SimpleFeatures()};
    }
    my %anames;
    foreach my $feat (@polyAseq_features){
      $anames{$feat->analysis->logic_name} = 1;
    }
    if (scalar keys %anames >= 4){
      return 1;
    }
  }
  return 0;
}



sub clip_ends {
  my ($transcript, $db_tr, $max_overhang, $slice_offset) = @_;
  $max_overhang ||= 1;
  my $clipped = 0;
  foreach my $db_exon (@{$db_tr->get_all_Exons}){
    #  ref:  #####-----#######-------#######------######
    #   tr:           ########-------########
          
    if ($db_exon != $db_tr->start_Exon and $db_exon != $db_tr->end_Exon){
      if ($transcript->seq_region_strand == 1 and 
          $transcript->start_Exon->seq_region_start < $db_exon->seq_region_start and
          ($transcript->start_Exon->seq_region_start + $max_overhang) >= $db_exon->seq_region_start and
          $transcript->start_Exon->seq_region_end == $db_exon->seq_region_end){
            $transcript->start_Exon->start($db_exon->seq_region_start - $slice_offset);
            $clipped = 1;
            print "Clipping transcript start to exon in ".$db_tr->stable_id."\n";
      }
      if ($transcript->seq_region_strand == -1 and 
          $transcript->end_Exon->seq_region_start < $db_exon->seq_region_start and
          ($transcript->end_Exon->seq_region_start + $max_overhang) >= $db_exon->seq_region_start and
          $transcript->end_Exon->seq_region_end == $db_exon->seq_region_end){
            $transcript->end_Exon->start($db_exon->seq_region_start - $slice_offset);
            $clipped = 1;
            print "Clipping transcript end to exon in ".$db_tr->stable_id."\n";
      }
      if ($transcript->seq_region_strand == 1 and 
          $transcript->end_Exon->seq_region_end > $db_exon->seq_region_end and
          ($transcript->end_Exon->seq_region_end - $max_overhang) <= $db_exon->seq_region_end and
          $transcript->end_Exon->seq_region_start == $db_exon->seq_region_start){
            $transcript->end_Exon->end($db_exon->seq_region_end - $slice_offset);
            $clipped = 1;
            print "Clipping transcript end to exon in ".$db_tr->stable_id."\n";
      }
      if ($transcript->seq_region_strand == -1 and 
          $transcript->start_Exon->seq_region_end > $db_exon->seq_region_end and
          ($transcript->start_Exon->seq_region_end - $max_overhang) <= $db_exon->seq_region_end and
          $transcript->start_Exon->seq_region_start == $db_exon->seq_region_start){
            $transcript->start_Exon->end($db_exon->seq_region_end - $slice_offset);
            $clipped = 1;
            print "Clipping transcript start to exon in ".$db_tr->stable_id."\n";
      }
    }
  }
  #print "PRECLIPPED=".$transcript->start."-".$transcript->end."  ".$transcript->seq_region_start."-".$transcript->seq_region_end."\n";
  #Recalculate transcript coordinates but, unlike with exons, keep start/end identical to seq_region_start/seq_region_end
  if ($clipped){
    $transcript->start(min(map{$_->seq_region_start} @{$transcript->get_all_Exons}));
    $transcript->end(max(map{$_->seq_region_end} @{$transcript->get_all_Exons}));
  }
  #print "CLIPPED=".$transcript->start."-".$transcript->end."  ".$transcript->seq_region_start."-".$transcript->seq_region_end."\n";
  #print join(":", map {$_->seq_region_start."-".$_->seq_region_end} @{$transcript->get_all_Exons})."\n";
  return $transcript;
}



=head2 assign_biotypes

 Arg[1]    : Bio::Vega::Gene object
 Function  : Assign biotypes to gene and transcripts
 Returntype: Bio::Vega::Gene object

=cut

sub assign_biotypes {
    my ($self, $gene) = @_;
    
    #Default gene biotype in the Ensembl API is protein_coding
    #Let this biotype be assigned if the gene has transcripts with translations
    $gene->biotype("processed_transcript");
    #Initial hypothesis: coding
    foreach my $transcript (@{$gene->get_all_Transcripts}){
        if ($transcript->translation){
            $transcript->biotype("protein_coding");
            $gene->biotype("protein_coding");
        }
        else{
            $transcript->biotype("processed_transcript");
        }            
    }
    
    #Alternative: long non-coding RNA
    if ($gene->biotype ne "protein_coding"){
        #Get gene slice and overlapping genes
        my @other_genes = @{LoutreWrite::GeneFilter::get_valid_overlapping_genes($gene)};
        #my $gene_slice = $gene->slice->adaptor->fetch_by_region("toplevel", $gene->seq_region_name, $gene->start, $gene->end);
        #my @other_genes = grep {$_->stable_id ne $gene->stable_id} @{$gene_slice->get_all_Genes};
        
        my $trans_count;
        my %readthrough;
        
        #1-Check same strand
        my $is_in_intron = 0;
        my %in_intron_genes;
        my $has_intronic_genes = 0;
        my %intronic_genes;
        
        #Look at pseudogenes overlapping gene
        foreach my $other_gene (grep {$_->seq_region_strand == $gene->seq_region_strand} @other_genes){
            if ($other_gene->biotype =~ /pseudo/){
                  #last; #?
            }
        }
        
        #Check for protein-coding genes in introns
        TRANSCRIPT1:
        foreach my $transcript (@{$gene->get_all_Transcripts}){
            $trans_count++;
            #Set default transcript biotype - can be changed later
            $transcript->biotype("processed_transcript");
            #Ignore readthrough transcripts
            if (scalar grep {$_->code eq "remark" and $_->value eq "readthrough"} @{$transcript->get_all_Attributes} ) {
                $readthrough{$transcript->stable_id} = 1;
                next TRANSCRIPT1;
            }
        
            foreach my $intron (@{$transcript->get_all_Introns}){
                my $intron_slice = $gene->slice->adaptor->fetch_by_region("toplevel", $intron->seq_region_name, $intron->seq_region_start, $intron->seq_region_end);
                my @intron_genes = grep {$_->stable_id ne $gene->stable_id and $_->seq_region_strand == $intron->seq_region_strand} @{$intron_slice->get_all_Genes_by_type("protein_coding")};
                foreach my $other_gene (@intron_genes){
                    #check overlap
                    if ($other_gene->seq_region_start > $intron->seq_region_start and $other_gene->seq_region_end < $intron->seq_region_end){
                        $has_intronic_genes = 1;
                        if (scalar @{$other_gene->get_all_Attributes('name')}){
                            $intronic_genes{$other_gene->get_all_Attributes('name')->[0]->value} = 1;
                        }
                        last TRANSCRIPT1;
                    }
                }
            }
        }

        if (scalar keys %readthrough == $trans_count){
            print "All Read-Throughs!\n";
            $has_intronic_genes = 0; #check this makes sense 
        }

        #Check if gene is in a protein-coding gene's intron
        GENE2:
        foreach my $other_gene (grep {$_->seq_region_strand == $gene->seq_region_strand and $_->biotype eq "protein_coding"} @other_genes){
            INTRON:
            foreach my $other_intron (@{$other_gene->get_all_Introns}){
                if ($other_intron->seq_region_start < $gene->seq_region_start and $other_intron->seq_region_end > $gene->seq_region_end){
                    $is_in_intron = 1;
                    if (scalar @{$other_gene->get_all_Attributes('name')}){
                        $in_intron_genes{$other_gene->get_all_Attributes('name')->[0]->value} = 1;
                    }
                    last INTRON;
                }
            }
            #make sure it's not overlapping any other exons
            foreach my $other_exon (@{$other_gene->get_all_Exons}){
                if ($other_exon->seq_region_end >= $gene->seq_region_start and $other_exon->seq_region_start <= $gene->seq_region_end){
                    $is_in_intron = 0;
                    last GENE2;
                }
            }
        }    
        
        
        
        #2-Check opposite strand
        #Antisense: transcripts overlapping one or more coding loci on the opposite strand, based on either exonic or intronic sequences.
        my $has_other_strand_genes = 0;
        my %other_strand_genes;            
         
        foreach my $other_gene (grep {$_->seq_region_strand != $gene->seq_region_strand and $_->biotype eq "protein_coding"} @other_genes){
            my $has_other_strand_transcripts = 0;
            TRANSCRIPT3:
            foreach my $other_trans (@{$other_gene->get_all_Transcripts}){
                #ignore non-coding transcripts
                if (!($other_trans->coding_region_start and $other_trans->coding_region_end)){
                    next TRANSCRIPT3;
                }
                #ignore read-thrus
                if (scalar grep {$_->code eq "remark" and $_->value eq "readthrough"} @{$other_trans->get_all_Attributes}){
                    next TRANSCRIPT3;
                }
                #check overlap
                if ($other_trans->seq_region_end > $gene->start and $other_trans->seq_region_start < $gene->end){
                    $has_other_strand_transcripts = 1;
                    last TRANSCRIPT3;
                }
            }

            if ($has_other_strand_transcripts){
                $has_other_strand_genes = 1;
                if (scalar @{$other_gene->get_all_Attributes('name')}){
        print $other_gene->stable_id."\n";
                    $other_strand_genes{$other_gene->get_all_Attributes('name')->[0]->value} = 1;
                }
            }
        }
        
    
        #Assign biotypes
        #Add 'overlapping locus' attribute if sense_intronic, sense_overlapping, 3'_overlapping_ncRNA
        #Antisense takes precedence over any other non-coding biotype with the exception of bidirectional_promoter_lncRNA
        my $biotype;
        if ($has_other_strand_genes){
            $biotype = "antisense";
        }
        elsif ($has_intronic_genes){
            $biotype = "sense_overlapping";
            $gene->add_Attributes(new Bio::EnsEMBL::Attribute(-code => 'remark', -value => 'overlapping locus'));
        }
        elsif ($is_in_intron){
            $biotype = "sense_intronic";
            $gene->add_Attributes(new Bio::EnsEMBL::Attribute(-code => 'remark', -value => 'overlapping locus'));            
        }
        else{
            $biotype = "lincRNA";
        }
        
        $gene->biotype($biotype);
        foreach my $transcript (@{$gene->get_all_Transcripts}){
            $transcript->biotype($biotype);
        }    
            
        #Add gene description, eg.
        #novel transcript, antisense to ARHGAP1 and DST
        #novel transcript, sense overlapping XBOX1
        #novel transcript, sense intronic to RPS3
        if (!$gene->description){
            if ($gene->biotype eq "lincRNA"){
                $gene->description("novel transcript");
            }
            elsif ($gene->biotype eq "antisense"){
                if (scalar keys %other_strand_genes){
                    $gene->description("novel transcript, antisense to ".join(" and ", keys %other_strand_genes)); 
                }
                else{
                    $gene->description("novel transcript");
                }	               
            }
            elsif ($gene->biotype eq "sense_overlapping"){
                if (scalar keys %intronic_genes){
                    $gene->description("novel transcript, sense overlapping ".join(" and ", keys %intronic_genes));   
                }
                else{
                    $gene->description("novel transcript");
                }             
            }    
            elsif ($gene->biotype eq "sense_intronic"){
                if (scalar keys %in_intron_genes){
                    $gene->description("novel transcript, sense intronic to ".join(" and ", keys %in_intron_genes));
                } 
                else{
                    $gene->description("novel transcript");
                }             
            }
        }
          
    }

    return $gene;
}



=head2 assign_status

 Arg[1]    : Bio::Vega::Gene object
 Function  : Assign status to gene and transcripts
 Returntype: Bio::Vega::Gene object

=cut

sub assign_status {
    my ($self, $gene) = @_;
    if ($gene->biotype =~ /^(lincRNA|sense_intronic|sense_overlapping|antisense)$/){
        $gene->status("UNKNOWN"); #This is the default status in the GTF file anyway
        foreach my $transcript (@{$gene->get_all_Transcripts}){
            $transcript->status("UNKNOWN");
        }
    }
    return $gene;
}



=head2 fetch_from_gene_by_transcript_name

 Arg[1]    : string (transcript name)
 Arg[2]    : Bio::Vega::Gene object
 Function  : Return a transcript having the query transcript name from the gene's transcript complement
 Returntype: Bio::Vega::Transcript object

=cut

sub fetch_from_gene_by_transcript_name {
  my ($tr_name, $gene) = @_;
  foreach my $transcript (@{$gene->get_all_Transcripts}){
    if ($transcript->get_all_Attributes('name')->[0] and $transcript->get_all_Attributes('name')->[0]->value eq $tr_name){
      return $transcript;
    }
  }
  return;
}



=head2 write_gene_region

 Arg[1]    : Bio::Otter::ServerAction::Script::Region
 Arg[2]    : Bio::Vega::Region
 Function  : Locks, edits and unlocks a region - Locking a region prevents it from being edited by other user - If a region is already locked by an annotator working on it, locking will fail
 Returntype: String (message with the outcome)

=cut

sub write_gene_region {
  my ($region_action, $region) = @_;

  my @msg;
  my $lock;
  try {
    $region_action->server->add_param( hostname => hostname );
    $lock = $region_action->lock_region;
    push @msg, 'lock ok';
  }
  catch {
    my ($err) = ($_ =~ m/^MSG: (Failed to lock.*)$/m);
    push @msg, "lock failed: '$err'";
    push @msg, "lock failed: '$_'";
  };

  if ($lock) {
    my $new_region;
    try {
      $region_action->server->set_params( data => $region, locknums => $lock->{locknums} );
      $new_region = $region_action->write_region; #NOTE: This resets slice to default (chromosome)
      push @msg, 'write ok';
    }
    catch {
      my $err = $_;
      chomp $err;
      push @msg, "write failed: '$err'";
    };
    $region_action->server->set_params( locknums => $lock->{locknums} );
    $region_action->unlock_region;
    push @msg, 'unlock ok';
    if ($new_region and scalar($new_region->genes)){
      push @msg, "gene ".join(",", map {$_->stable_id} $new_region->genes);
    }
  }

  return join(',', @msg);
}



1;
