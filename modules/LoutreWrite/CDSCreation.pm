
package LoutreWrite::CDSCreation;

use strict;
use warnings;
use base 'Exporter';
our @EXPORT = qw( sort_by_categ get_host_gene_cds_set assign_cds_to_transcripts get_appris_tag get_ccds_tag cds_exon_chain start_codon_fits has_complete_cds get_host_gene_start_codon_set is_retained_intron %HOST_CDS_SET %HOST_START_CODON_SET %HOST_STOP_CODON_SET );
use Bio::Vega::Translation;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(calculate_exon_phases);
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use LoutreWrite::AnnotUpdate qw( has_polyA_site_support has_polyAseq_support );
use LoutreWrite::Config;

our %HOST_CDS_SET;
our %HOST_START_CODON_SET;
our %HOST_STOP_CODON_SET;
our %CORE_TRANSCRIPTS;

#my $registry = 'Bio::EnsEMBL::Registry';
#$registry->load_registry_from_db(
#    -host => 'ensembldb.ensembl.org',
#    -user => 'anonymous'
#);

#my $core_slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );
#my $core_gene_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Gene' );


# my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
#     -host   => 'mysql-ens-havana-prod-1',
#     -port   => '4581',
#     -user   => 'ensro',
#     -dbname => 'homo_sapiens_core_98_38_status',
#     -driver => 'mysql'
#   );
# $db->dbc->disconnect_when_inactive(1);
# my $core_slice_adaptor = $db->get_SliceAdaptor();
# my $core_gene_adaptor = $db->get_GeneAdaptor();


#https://docs.google.com/document/d/17PIripFaSkGkHmZbwR1ZCRYQhyu5ov6dgBVNqVtDjt4/edit?ts=5b6b0748

#Take annotated transcripts in host gene
#Find all full-lenght CDSs and sort them by:
# -APPRIS level (principal 1 to 4)
# -CCDS (5'-most ATG first)
# -complete CDSs (5'-most ATG first)
# -complete NMD CDSs (5'-most ATG first) - only use complete NMD CDSs that are not also found as non-NMD CDSs in loutre
#Find all start codons and sort by 5'-most first

#Search for each of these CDSs in comp_pipe models until one is found (including start codon, and stop codon if available) - copy biotype ('known_CDS', 'novel_CDS', 'putative_CDS' or 'coding')

#If none is suitable for a comp_pipe model, try finding ORFs starting from known ATGs (CDS_end_NF allowed) sorted by:
# -APPRIS level (principal 1 to 3)
# -CCDS (5'-most ATG first)
# -other CDSs (5'-most ATG first)

#Search for ORFs starting from these start codons - if stop codon and not NMD, call it 'putative_CDS'
# if no stop codon, add 'end not found' tag

#Always check if NMD rules apply -
#  -stop codon is found 50bps or more upstream of a splice donor site
#  -NMD CDS must be greater than 35 aa

#Check for retained intron rules - no CDS will be assigned
# -full intron retention: complete CDS intron(s)
# -partial intron retention: model has a start or end coordinate within a CDS intron of a loutre model, but not (respectively) a splice donor site or splice acceptor site within that same intron

#Readthrough comp_pipe models do not get any CDS


#NOTE: needs to deal with merged transcripts where the CDS must be modified, 
# eg. if the pre-existing transcript was 'CDS end not found' and it has been extended.

#Ensure MANE transcripts are not modified.



=head2 assign_cds_to_transcripts

 Arg[1]    : Bio::Vega::Gene or Bio::Vega::Transcript object
 Arg[2]    : Bio::Vega::Gene object
 Function  : assign a CDS to every transcript in the first gene, based on the CDSs and start codons in the second gene
 Returntype: none

=cut

sub assign_cds_to_transcripts {
  my ($novel_obj, $host_gene, $slice_offset) = @_;
  my @novel_transcripts;
  if ($novel_obj->isa("Bio::Vega::Gene")){
    @novel_transcripts = @{$novel_obj->get_all_Transcripts};
  }
  elsif ($novel_obj->isa("Bio::Vega::Transcript")){
    @novel_transcripts = ($novel_obj);
  }
  else{
    die "First parameter must be a Bio::Vega::Gene or Bio::Vega::Transcript object!";
  }
  #my $slice_offset = $host_gene->seq_region_start - 1;
  my $cds_set;
  if ($HOST_CDS_SET{$host_gene->stable_id}){
    $cds_set = $HOST_CDS_SET{$host_gene->stable_id};
    #print "Using HOST_CDS_SET\n";
  }
  else{
    $cds_set = get_host_gene_cds_set($host_gene);
    $HOST_CDS_SET{$host_gene->stable_id} = $cds_set;
  }
  my $start_codon_set;
  if ($HOST_START_CODON_SET{$host_gene->stable_id}){
    $start_codon_set = $HOST_START_CODON_SET{$host_gene->stable_id};
    #print "Using HOST_START_CODON_SET\n";
  }
  else{
    $start_codon_set = get_host_gene_start_codon_set($host_gene);
    $HOST_START_CODON_SET{$host_gene->stable_id} = $start_codon_set;
  }
  my $stop_codon_set;
  if ($HOST_STOP_CODON_SET{$host_gene->stable_id}){
    $stop_codon_set = $HOST_STOP_CODON_SET{$host_gene->stable_id};
    #print "Using HOST_STOP_CODON_SET\n";
  }
  else{
    $stop_codon_set = get_host_gene_stop_codon_set($host_gene);
    $HOST_STOP_CODON_SET{$host_gene->stable_id} = $stop_codon_set;
  }  
  
        #print "HOST_CDS_SET=".scalar(keys %HOST_CDS_SET)."\n";
        #print "HOST_START_CODON_SET=".scalar(keys %HOST_START_CODON_SET)."\n";
  print "HOST GENE START = ".$host_gene->seq_region_start." HOST GENE END = ".$host_gene->seq_region_end."\n";
  TR:foreach my $transcript (@novel_transcripts){
    my $t_name = $transcript->stable_id || $transcript->get_all_Attributes('hidden_remark')->[0]->value;
    print "\nFinding a CDS for transcript $t_name\n";
    print "TR_START=".$transcript->seq_region_start."; TR_END=".$transcript->seq_region_end."\n";
    
    #Find out if the transcript has polyA support - otherwise it can't be made coding
#    unless (LoutreWrite::Default::has_polyA_site_support($transcript, 500) or LoutreWrite::Default::has_polyAseq_support($transcript, 500)){
#      print "No polyA site support found\n";
#      next TR;
#    }
    
    print "\nAvailable CDSs:\n";
    foreach my $unique_cds_tr (@$cds_set){
      print "CDS_START=".$unique_cds_tr->coding_region_start."; CDS_END=".$unique_cds_tr->coding_region_end."\n";
    }
    #Try to find a suitable CDS from the host gene
    foreach my $unique_cds_tr (@$cds_set){
      if (cds_fits($transcript, $unique_cds_tr, $slice_offset)){
        print "CDS FITS: CDS_START=".$unique_cds_tr->coding_region_start."; CDS_END=".$unique_cds_tr->coding_region_end."\n";
        my ($cds_start, $cds_end) = $transcript->seq_region_strand == 1 ? ($unique_cds_tr->coding_region_start, $unique_cds_tr->coding_region_end) : ($unique_cds_tr->coding_region_end, $unique_cds_tr->coding_region_start);
        $cds_start = $cds_start + $slice_offset;
        $cds_end = $cds_end + $slice_offset;
        create_cds($transcript, $cds_start, $cds_end);
        #Re-assign biotype and status
        if (predicted_nmd_transcript($transcript, $cds_end)){
          $transcript->biotype("nonsense_mediated_decay");
        }
        else{
          if ($unique_cds_tr->biotype eq "nonsense_mediated_decay"){
            $transcript->biotype("protein_coding");
            $transcript->status("PUTATIVE");
          }
          else{
            $transcript->biotype($unique_cds_tr->biotype);
            $transcript->status($unique_cds_tr->status);
          }
        }
        print $unique_cds_tr->stable_id." CDS matches ".$transcript->stable_id." (".$transcript->biotype."-".$transcript->status.")\n";
        next TR;
      }
    }
    #Else, try to create a CDS using a start codon from the host gene    
#    #But first, check for intron retention - SKIP THIS: DONE EARLIER
#    print "\nChecking for intron retention...\n";
#    if (is_retained_intron ($transcript, $host_gene, 5)){
#      $transcript->biotype("retained_intron");
#      print "Found retained_intron\n";
#      next TR;
#    }
    print "\nChecking start codons...\n";
    foreach my $unique_start_codon (@$start_codon_set){
      #Start codon coordinates are relative to the region slice: convert to genomic
      print "start_codon=".($unique_start_codon+$slice_offset)."; transcript=".$transcript->seq_region_start."-".$transcript->seq_region_end."\n";
      my ($cds_start, $cds_end, $fl_cds) = start_codon_fits($transcript, $unique_start_codon + $slice_offset); 
      if ($cds_start and $cds_end){
        #Found novel CDS from know start codon 
        #Now, check for intron retention
        #If stop codon is upstream of retained intron, do not make it retained_intron (it will probably be NMD)
        if (is_retained_intron($transcript, $host_gene, 5, $cds_end)){
          print "Found retained_intron before stop codon at $cds_end\n";
          $transcript->biotype("retained_intron");
          next TR;
        }
        #Full length CDS?
        if ($fl_cds){
          #PolyA support or know stop codon required if full-length CDS
          unless (LoutreWrite::AnnotUpdate::has_polyA_site_support($transcript, 500) or
                  LoutreWrite::AnnotUpdate::has_polyAseq_support($transcript, 500) or
                  LoutreWrite::AnnotUpdate::has_last_exon_polyA_site_support($transcript, $host_gene) or
                  is_known_stop_codon($cds_end, $slice_offset, $host_gene)){
            print "Stop codon at $cds_end has no polyA site support nor has been annotated before\n";
            next TR;
          }
        }
        #Annotate CDS
        if (create_cds($transcript, $cds_start, $cds_end)){
          #Re-assign biotype and status
          print "TTTTTTTT\n";
          if (predicted_nmd_transcript($transcript, $cds_end)){
            $transcript->biotype("nonsense_mediated_decay"); print "XXXXXXX\n";
          }
          else{
            $transcript->biotype("protein_coding");
            $transcript->status("PUTATIVE");
          }
          #Add end_NF attributes if 3'-incomplete CDS
          add_end_NF_attributes($transcript, $slice_offset);
          print "".($unique_start_codon + $slice_offset)." start codon matches ".$transcript->stable_id." (".$transcript->biotype."-".$transcript->status.")\n";
          next TR;
        }
      }
    }
    print "No match for transcript ".$transcript->stable_id."\t".$transcript->biotype."\n";
    #Non-coding transcript: check for intron retention 
    print "\nChecking for intron retention...\n";
    if (is_retained_intron($transcript, $host_gene, 5)){
      $transcript->biotype("retained_intron");
      print "Found retained_intron\n";
    }
  }
}



=head2 is_retained_intron

 Arg[1]    : Bio::Vega::Transcript object
 Arg[2]    : Bio::Vega::Gene object
 Arg[3]    : integer (minimum overhang in partial intron retention)
 Function  : finds whether the novel transcript is a retained_intron model in the host gene
 Returntype: boolean

=cut

sub is_retained_intron {
  my ($transcript, $host_gene, $min_overhang, $stop_codon_coordinate) = @_;
  #i) Full intron retention: the model fully transcribes the sequence between the splice donor site and
  #   acceptor site of adjacent coding exons (i.e. a ‘CDS intron’) of an annotated transcript. This could include scenarios
  #   where the TAGENE model retains multiple CDS introns of the same transcript. 
  #ii) Partial intron retention: the model has a start or end coordinate within a CDS intron of an annotated transcript,
  #   but not (respectively) a splice donor site or splice acceptor site within that same intron. Do we need some kind
  #   of size cut-off here? I.e. for fuzzy edges that are misalignments rather than real retained introns. Should we clip these? 
  my $slice_offset = $host_gene->seq_region_start - 1; #print "OFFSET=$slice_offset\n";
  $min_overhang ||= 1;
  my @comp_trs;
  my $mane_select;
  #print "\nLOOKING AT INTRON RETENTION...\n";
  foreach my $tr (@{$host_gene->get_all_Transcripts}){
    #Ignore lingering OTT transcripts
    next unless $tr->stable_id =~ /^ENS([A-Z]{3})?T000/;
    #Ignore "not for VEGA" transcripts unless they have a "comp_pipe" biotype or a "TAGENE_transcript" remark
    if (scalar(grep {$_->value eq "not for VEGA"} @{$tr->get_all_Attributes('remark')})){
      next unless ($tr->biotype eq "comp_pipe" or scalar(grep {$_->value eq "TAGENE_transcript"} @{$tr->get_all_Attributes('remark')}));
    }
    #Compare with MANE Select transcript first?...
    if (scalar(grep {$_->value eq "MANE_select"} @{$tr->get_all_Attributes('remark')})){
      $mane_select = $tr;
    }
    elsif ($tr->translate){
      push(@comp_trs, $tr);
    }
  }
  if ($mane_select){
    unshift(@comp_trs, $mane_select);
  }


  #Firstly, find which of the input transcript's exons exist in annotated coding transcripts
  #These must not be queried for intron retention
  my %annotated_exons;
  foreach my $tr (@comp_trs){
    if ($tr->translate){
      foreach my $exon1 (@{$tr->get_all_Exons}){
        my $rank = 0; #not using the 'rank' API method because it relies on exon stable ids, which may not have been assigned yet
        foreach my $exon2 (@{$transcript->get_all_Exons}){
          $rank++;
          if ($exon1->seq_region_start == $exon2->seq_region_start and $exon1->seq_region_end == $exon2->seq_region_end){
            $annotated_exons{$rank}++;
            #print "SEEN=".$rank."\n";
          }
        }
      }
    }
  }

  #Now loop through all the coding transcripts again and search for intron retention in the input transcript
  foreach my $tr (@comp_trs){
    if ($tr->translate){
      #print "\nTR:  ".$tr->stable_id." ".$tr->biotype."\n"; 
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
          my $rank = 0;
          foreach my $exon (@{$transcript->get_all_Exons}){
            #print "EXON: ".$exon->seq_region_start."-".$exon->seq_region_end."\n";
            $rank++;
            if ($annotated_exons{$rank}){
              #print "ANNOTATED\n";
              next;
            }
            if ($exon->seq_region_start < $intron->seq_region_start and $exon->seq_region_end > $intron->seq_region_end){
              #print "FOUND\n";
              if ($stop_codon_coordinate){ #if CDS end coordinate, do not make retained_intron if CDS end is upstream of intron
                if (($transcript->seq_region_strand == 1 and $stop_codon_coordinate > $intron->seq_region_start) or
                    ($transcript->seq_region_strand == -1 and $stop_codon_coordinate < $intron->seq_region_end)
                ){
                  return 1;
                }
              }
              else{
                return 1;
              }
            }
          }
          #Partial intron retention
          #  ref: #######-------#######------######
          #   tr:          #######---####
          #a) partial intron retention of the start exon (forward strand)
#print "\nEntering partial intron retention\n";
#print "START_EXON: ".$transcript->start_Exon->seq_region_start."-".$transcript->start_Exon->seq_region_end."\n";
#print "END_EXON: ".$transcript->end_Exon->seq_region_start."-".$transcript->end_Exon->seq_region_end."\n";
#print "INTRON: ".$intron->seq_region_start."-".$intron->seq_region_end."\n\n";
          if ($transcript->seq_region_strand == 1 and 
              $transcript->start_Exon->seq_region_start >= $intron->seq_region_start and
              $transcript->start_Exon->seq_region_end > $intron->seq_region_end and
              ($transcript->start_Exon->seq_region_start + $min_overhang) <= $intron->seq_region_end){
#print "-case 1\n";
            unless ($annotated_exons{1}){
              return 1;
            }
          }
          #b) partial intron retention of the end exon (reverse strand)
          if ($transcript->seq_region_strand == -1 and 
              $transcript->end_Exon->seq_region_start >= $intron->seq_region_start and
              $transcript->end_Exon->seq_region_end > $intron->seq_region_end and
              ($transcript->end_Exon->seq_region_start + $min_overhang) <= $intron->seq_region_end){
#print "-case 2\n";
                #Check polyA site support
                unless (LoutreWrite::AnnotUpdate::has_polyA_site_support($transcript, 500) or
                        LoutreWrite::AnnotUpdate::has_polyAseq_support($transcript, 500) or
                        LoutreWrite::AnnotUpdate::has_last_exon_polyA_site_support($transcript, $host_gene) 
                        ){
                  unless ($annotated_exons{scalar(@{$transcript->get_all_Exons})}){
                    print "No polyA site support - will make a retained_intron\n";
                    return 1;
                  }
                }
          }
          #c) partial intron retention of the end exon (forward strand)
          if ($transcript->seq_region_strand == 1 and 
              $transcript->end_Exon->seq_region_start < $intron->seq_region_start and
              $transcript->end_Exon->seq_region_end <= $intron->seq_region_end and
              ($transcript->end_Exon->seq_region_end - $min_overhang) >= $intron->seq_region_start){
#print "-case 3\n";
                #print "AAA - LAST INTRON\n";
                if (LoutreWrite::AnnotUpdate::has_polyAseq_support($transcript, 500)){
                  print "FOUND POLYASEQ\n";
                }
                #Check polyA site support
                unless (LoutreWrite::AnnotUpdate::has_polyA_site_support($transcript, 500) or
                        LoutreWrite::AnnotUpdate::has_polyAseq_support($transcript, 500) or
                        LoutreWrite::AnnotUpdate::has_last_exon_polyA_site_support($transcript, $host_gene)
                        ){
                  unless (scalar(@{$transcript->get_all_Exons})){
                    print "No polyA site support - will make a retained_intron\n";
                    return 1;
                  }
                }
          }
          #d) partial intron retention of the start exon (reverse strand)
          if ($transcript->seq_region_strand == -1 and 
              $transcript->start_Exon->seq_region_start < $intron->seq_region_start and
              $transcript->start_Exon->seq_region_end <= $intron->seq_region_end and
              ($transcript->start_Exon->seq_region_end - $min_overhang) >= $intron->seq_region_start){
#print "-case 4\n";
            unless ($annotated_exons{1}){
              return 1;
            }
          }
        }
      }
    }
  }
  return 0;
}



=head2 get_host_gene_cds_set

 Arg[1]    : Bio::Vega::Gene object
 Function  : returns transcripts representing the set of unique complete CDSs in the host gene, 
             sorted as defined in 'sort_by_categ'. For each unique CDS chain there is also a transcript.
 Returntype: arrayref of Bio::Vega::Transcript objects

=cut

sub get_host_gene_cds_set {
  my $gene = shift;
  my @cds_set;
  my %seen_chains;
  
  my @filtered_transcripts = grep {$_->translate and has_complete_cds($_)} @{$gene->get_all_Transcripts};
 
  foreach my $transcript (sort_by_categ(\@filtered_transcripts)){
    #Ignore lingering OTT transcripts
    next unless $transcript->stable_id =~ /^ENS([A-Z]{3})?T000/;
    #Ignore "not for VEGA" transcripts unless they have a "comp_pipe" biotype or a "TAGENE_transcript" remark
    if (scalar(grep {$_->value eq "not for VEGA"} @{$transcript->get_all_Attributes('remark')})){
      next unless ($transcript->biotype eq "comp_pipe" or scalar(grep {$_->value eq "TAGENE_transcript"} @{$transcript->get_all_Attributes('remark')}));
    }
    
    my $chain = cds_exon_chain($transcript);
    unless ($seen_chains{$chain}){
      #push(@cds_set, {chain => $transcript->get_all_translateable_Exons, id=> $transcript->stable_id, biotype => $transcript->biotype, status => $transcript->status});
      push(@cds_set, $transcript);
      $seen_chains{$chain} = 1;
    }
  }
  return \@cds_set;
}



=head2 sort_by_categ

 Arg[1]    : arrayref of Bio::Vega::Transcript objects
 Function  : sorts coding transcripts according to the categories below:
             1 - MANE Select,
             2 - APPRIS principal 1 to 4, 
             3 - CCDS (5'-most ATG first), 
             4 - Other complete CDS (5'-most ATG first),
             5 - Complete CDS NMD (5'-most ATG first)
 Returntype: arrayref of Bio::Vega::Transcript objects

=cut


sub sort_by_categ {
  my $transcripts = shift;

  my @sorted_transcripts = sort {
    my $a_mane_select = scalar grep {$_->value eq "MANE_select"} @{$a->get_all_Attributes("remark")};
    my $b_mane_select = scalar grep {$_->value eq "MANE_select"} @{$b->get_all_Attributes("remark")};
    my $a_appris = get_appris_tag($a);
    my $b_appris = get_appris_tag($b);
    my $a_ccds = get_ccds_tag($a);
    my $b_ccds = get_ccds_tag($b);
    my $a_cds_start = $a->coding_region_start;
    my $b_cds_start = $b->coding_region_start;

    if ($a_mane_select){
      return -1;
    }
    elsif ($b_mane_select){
      return 1;
    }    
    else{
      if ($a_appris and $b_appris and $a_appris=~/principal(1|2|3|4)/ and $b_appris=~/principal(1|2|3|4)/){
        return $a_appris cmp $b_appris;
      }
      elsif ($a_appris and $a_appris=~/principal(1|2|3|4)/){
        return -1;
      }
      elsif ($b_appris and $b_appris=~/principal(1|2|3|4)/){
        return 1;
      }
      else{
        if ($a_ccds and $b_ccds){
          if ($a->seq_region_strand == 1){
            return $a_cds_start <=> $b_cds_start;
          }
          else{
            return $b_cds_start <=> $a_cds_start;
          }
        }
        elsif ($a_ccds){
          return -1;
        }
        elsif ($b_ccds){
          return 1;
        }
        else{
          if ($a->biotype eq $b->biotype){
            if ($a->seq_region_strand == 1){
              return $a_cds_start <=> $b_cds_start;
            }
            else{
              return $b_cds_start <=> $a_cds_start;
            }
          }
          elsif ($a->biotype eq "protein_coding"){
             return -1;
          }
          elsif ($a->biotype eq "nonsense_mediated_decay"){
             return 1;
          }
        }
      }
    }
  } grep {$_->translate} @$transcripts;

  return @sorted_transcripts;
}



=head2 cds_fits

 Arg[1]    : Bio::Vega::Transcript object
 Arg[2]    : Bio::Vega::Transcript object
 Function  : returns true if the CDS of the second transcript is suitable for the first transcript
 Returntype: boolean

=cut

sub cds_fits {
  my ($acceptor_transcript, $donor_transcript, $slice_offset) = @_;
  #$slice_offset = $donor_transcript->seq_region_start - 1;
  #Check CDS introns, start_codon, stop codon
  my $cds_start = $donor_transcript->coding_region_start + $slice_offset;
  my $cds_end = $donor_transcript->coding_region_end + $slice_offset;
  if ($cds_start > $cds_end){
    ($cds_start, $cds_end) = ($cds_end, $cds_start);
  }
  print "AAA1:".$donor_transcript->coding_region_start."\n";
  print "AAA2:".$donor_transcript->coding_region_end."\n";
  print "BBB:".$donor_transcript->cdna_coding_start."\n";
  print "CCC:".$donor_transcript->stable_id."\n";
  print "CDS_START=$cds_start\n";
  print "CDS_END=$cds_end\n";


  #Start codon must fall in exon
  my $cds_start_present;
  foreach my $exon (@{$acceptor_transcript->get_all_Exons}){
    print "EXON_START=".$exon->seq_region_start."\n";
    if ($exon->seq_region_start <= $cds_start and $exon->seq_region_end >= $cds_start){
      $cds_start_present = 1; print "cds_start_present=1\n";
      last;
    }
  }
  #Stop codon must fall in exon
  my $cds_end_present;
  foreach my $exon (@{$acceptor_transcript->get_all_Exons}){
    if ($exon->seq_region_start <= $cds_end and $exon->seq_region_end >= $cds_end){
      $cds_end_present = 1; print "cds_end_present=1\n";
      last;
    }
  }
  #print "cds_start_present=$cds_start_present; cds_end_present=$cds_end_present\n";
  #CDS introns must coincide
  print "acceptor_intron_chain($cds_start, $cds_end)=".intron_chain($acceptor_transcript, $cds_start, $cds_end)."\n";
  print "donor_intron_chain($cds_start, $cds_end)=".intron_chain($donor_transcript, $cds_start, $cds_end)."\n";
  if ($cds_start_present and $cds_end_present and
      intron_chain($acceptor_transcript, $cds_start, $cds_end) eq intron_chain($donor_transcript, $cds_start, $cds_end)){
      print "cds_start, cds_end and intron chain match\n";
    #Avoid making an NMD transcript with a CDS smaller than 35 aa
    #The true CDS end is required now
    $cds_end = ($donor_transcript->seq_region_strand == 1) ? $cds_end : $cds_start;
    unless (predicted_nmd_transcript($acceptor_transcript, $cds_end) and $donor_transcript->translation->length < 35){
      return 1;
    }
  }
  return 0;
}



=head2 has_complete_cds

 Arg[1]    : Bio::Vega::Transcript object
 Function  : returns true if transcript has neither 'CDS start not found' nor 'CDS end not found' attributes
 Returntype: boolean

=cut

sub has_complete_cds {
  my $transcript = shift;
  if (!($transcript->translate)){
    return 0;
  }
  foreach my $attribute (@{$transcript->get_all_Attributes('cds_start_NF')}){
    if ($attribute->value == 1){
      return 0;
    }
  }
  foreach my $attribute (@{$transcript->get_all_Attributes('cds_end_NF')}){
    if ($attribute->value == 1){
      return 0;
    }
  }
  return 1;
}



=head2 has_cds_start

 Arg[1]    : Bio::Vega::Transcript object
 Function  : returns true if transcript has no 'cds start not found' attribute
 Returntype: boolean

=cut

sub has_cds_start {
  my $transcript = shift;
  if (!($transcript->translate)){
    return 0;
  }
  foreach my $attribute (@{$transcript->get_all_Attributes('cds_start_NF')}){
    if ($attribute->value == 1){
      return 0;
    }
  }
  return 1;
}


=head2 has_cds_end

 Arg[1]    : Bio::Vega::Transcript object
 Function  : returns true if transcript has no 'cds end not found' attribute
 Returntype: boolean

=cut

sub has_cds_end {
  my $transcript = shift;
  if (!($transcript->translate)){
    return 0;
  }
  foreach my $attribute (@{$transcript->get_all_Attributes('cds_end_NF')}){
    if ($attribute->value == 1){
      return 0;
    }
  }
  return 1;
}



=head2 get_appris_tag

 Arg[1]    : Bio::Vega::Transcript object
 Function  : gets the APPRIS attribute value (or undef if none) for a loutre transcript 
             by looking at their counterparts in the latest Ensembl core database.
             If the same transcript is not found, other transcript with the same CDS would suffice
 Returntype: string

=cut

sub get_appris_tag {
  my $transcript = shift;
  my $core_transcript = get_core_transcript($transcript);
  if ($core_transcript){
    #print $core_transcript->stable_id."\n";
    if (scalar @{$core_transcript->get_all_Attributes('appris')}){
      #print $core_transcript->get_all_Attributes('appris')->[0]->value."\n";
      return $core_transcript->get_all_Attributes('appris')->[0]->value;
    }
  }
  else{
    print "No core transcript for ".$transcript->stable_id."\n";
  }
  return undef;
}



=head2 get_ccds_tag

 Arg[1]    : Bio::Vega::Transcript object
 Function  : gets the CCDS attribute value (or undef if none) for a loutre transcript 
             by looking at their counterparts in the latest Ensembl core database
 Returntype: string

=cut

sub get_ccds_tag {
  my $transcript = shift;
  my $core_transcript = get_core_transcript($transcript);
  if ($core_transcript){
    if (scalar @{$core_transcript->get_all_Attributes('ccds_transcript')}){
      return $core_transcript->get_all_Attributes('ccds_transcript')->[0]->value;
    }
  }
  return undef;
}



=head2 cds_exon_chain

 Arg[1]    : Bio::Vega::Transcript object
 Function  : gets the exon coordinates of the coding portion of a transcript
 Returntype: string

=cut

sub cds_exon_chain {
  my $transcript = shift;
  if ($transcript->translate){
    my @cds_exons = @{$transcript->get_all_translateable_Exons};
    #return \@cds_exons;
    #return join("_", map{$_->seq_region_start."-".$_->seq_region_end} @cds_exons);
    return join("_", map{$_->seq_region_start."-".$_->seq_region_end} sort {$a->seq_region_start <=> $b->seq_region_start} @cds_exons);
  }
  return undef;
}



=head2 intron_chain

 Arg[1]    : Bio::Vega::Transcript object
 Arg[2]    : start coordinate (optional)
 Arg[3]    : end coordinate (optional)
 Function  : gets the intron coordinates of a transcript (or of a transcript's segment defined 
             by the genomic coordinates provided)
 Returntype: string

=cut

sub intron_chain {
  my ($transcript, $from, $to) = @_;
  if (($from and !$to) or ($to and !$from)){
    die "Need both start and end coordinates or none!";
  }
  unless ($from and $to){
    $from = $transcript->seq_region_start;
    $to = $transcript->seq_region_end;
  }
  my @introns = grep {$_->seq_region_start >= $from and $_->seq_region_end <= $to} @{$transcript->get_all_Introns};

  return join("_", map {$_->seq_region_start."-".$_->seq_region_end} sort {$a->seq_region_start <=> $b->seq_region_start} @introns);
}



=head2 get_core_transcript

 Arg[1]    : Bio::Vega::Transcript object
 Function  : gets an Ensembl core transcript having the same CDS as the query transcript
 Returntype: Bio::EnsEMBL:Transcript object

=cut

sub get_core_transcript {
  my $transcript = shift; 
  if ($CORE_TRANSCRIPTS{$transcript->stable_id}){
    return $CORE_TRANSCRIPTS{$transcript->stable_id};
  }
  else{
    print "Finding a core transcript for ".$transcript->stable_id."\n";
    my $sa = $DBA{'core'}->get_SliceAdaptor();
    my $core_slice = $sa->fetch_by_region("toplevel", $transcript->seq_region_name, $transcript->seq_region_start, $transcript->seq_region_end);
    if ($core_slice){
      foreach my $core_transcript (grep {$_->seq_region_strand == $transcript->seq_region_strand} @{$core_slice->get_all_Transcripts}){
        my $core_cds_exon_chain = cds_exon_chain($core_transcript);
        my $tr_cds_exon_chain = cds_exon_chain($transcript);
        if ($core_cds_exon_chain and $tr_cds_exon_chain and $core_cds_exon_chain eq $tr_cds_exon_chain){
          $CORE_TRANSCRIPTS{$transcript->stable_id} = $core_transcript;
          return $core_transcript;
        }
      }
    }
  }
  return undef;
}



=head2 get_host_gene_start_codon_set

 Arg[1]    : Bio::Vega::Gene object
 Function  : returns the set of unique CDS start positions (except for cds_start_NF transcripts) in the host gene, sorted as defined in 'sort_by_categ_2'.
 Returntype: arrayref of integers

=cut

sub get_host_gene_start_codon_set {
  my $gene = shift;
  #print "G_START=".$gene->seq_region_start."; G_END=".$gene->seq_region_end."\n";
  #print "G_START=".$gene->start."; G_END=".$gene->end."\n";
  my @start_codon_set;
  my %seen_sc;
  my @filtered_transcripts = grep {has_cds_start($_)} @{$gene->get_all_Transcripts};
  foreach my $transcript (sort_by_categ_2(\@filtered_transcripts)){
    #Ignore lingering OTT transcripts
    next unless $transcript->stable_id =~ /^ENS([A-Z]{3})?T000/;
    #Ignore "not for VEGA" transcripts unless they have a "comp_pipe" biotype or a "TAGENE_transcript" remark
    if (scalar(grep {$_->value eq "not for VEGA"} @{$transcript->get_all_Attributes('remark')})){
      next unless ($transcript->biotype eq "comp_pipe" or scalar(grep {$_->value eq "TAGENE_transcript"} @{$transcript->get_all_Attributes('remark')}));
    }
    
    #$transcript = $transcript->transform("chromosome");
    #print "T_START=".$transcript->seq_region_start."; T_END=".$transcript->seq_region_end."\n";
    my $cds_start = $transcript->seq_region_strand == 1 ? $transcript->coding_region_start : $transcript->coding_region_end;
    #print "CDS_START=$cds_start\n";
    unless ($seen_sc{$cds_start}){
      push(@start_codon_set, $cds_start);
      $seen_sc{$cds_start} = 1;
    }
  }
  return \@start_codon_set;
}



=head2 sort_by_categ_2

 Arg[1]    : arrayref of Bio::Vega::Transcript objects
 Function  : sorts coding transcripts according to the categories below:
             1 - MANE Select,
             2 - APPRIS principal 1 to 3, 
             3 - CCDS (5'-most ATG first), 
             4 - Other 5'-complete CDS (5'-most ATG first)
 Returntype: arrayref of Bio::Vega::Transcript objects

=cut

sub sort_by_categ_2 {
  my $transcripts = shift;

  my @sorted_transcripts = sort {
    my $a_mane_select = scalar grep {$_->value eq "MANE_select"} @{$a->get_all_Attributes("remark")};
    my $b_mane_select = scalar grep {$_->value eq "MANE_select"} @{$b->get_all_Attributes("remark")};
    my $a_appris = get_appris_tag($a);
    my $b_appris = get_appris_tag($b);
    my $a_ccds = get_ccds_tag($a);
    my $b_ccds = get_ccds_tag($b);
    my $a_cds_start = $a->coding_region_start;
    my $b_cds_start = $b->coding_region_start;

    if ($a_mane_select){
      return -1;
    }
    elsif ($b_mane_select){
      return 1;
    }
    else{
      if ($a_appris and $b_appris and $a_appris=~/principal(1|2|3)/ and $b_appris=~/principal(1|2|3)/){
        return $a_appris cmp $b_appris;
      }
      elsif ($a_appris and $a_appris=~/principal(1|2|3)/){
        return -1;
      }
      elsif ($b_appris and $b_appris=~/principal(1|2|3)/){
        return 1;
      }
      else{
        if ($a_ccds and $b_ccds){
          if ($a->seq_region_strand == 1){
            return $a_cds_start <=> $b_cds_start;
          }
          else{
            return $b_cds_start <=> $a_cds_start;
          }
        }
        elsif ($a_ccds){
          return -1;
        }
        elsif ($b_ccds){
          return 1;
        }
        else{
          if ($a->seq_region_strand == 1){
            return $a_cds_start <=> $b_cds_start;
          }
          else{
            return $b_cds_start <=> $a_cds_start;
          }
        }
      }
    }
  } grep {$_->translate} @$transcripts;

  return @sorted_transcripts;
}



=head2 start_codon_fits

 Arg[1]    : Bio::Vega::Transcript object
 Arg[2]    : integer (genomic position)
 Function  : returns true if a CDS beginning from the given genomic position can be built for the transcript
 Returntype: array of integers corresponding to the CDS start and end genomic coordinates

=cut

sub start_codon_fits {
  my ($transcript, $cds_start) = @_;

  #Start codon must fall in exon
  foreach my $exon (@{$transcript->get_all_Exons}){
    if ($exon->seq_region_start <= $cds_start and $exon->seq_region_end >= $cds_start){
      print "cds_start $cds_start in exon ".$exon->seq_region_start."-".$exon->seq_region_end."\n";
      #Find start position in transcript coordinates
      my @coords = $transcript->genomic2cdna($cds_start, $cds_start, $transcript->seq_region_strand);
      my $tr_cds_start = $coords[0]->start;
      #Search for a putative ORF in the transcript subsequence 
      #starting from the provided CDS start
      my $subseq = substr($transcript->seq->seq, $tr_cds_start-1);
      #print $transcript->stable_id."\t".$subseq."\n";
      #Does it contain a complete ORF?
      if ($subseq =~ /^(ATG([ATGC]{3})*?(TGA|TAA|TAG))/){
        print "*".$1."*".length($1)."\n";
        #Find CDS end in genomic coordinates
        my $cds_length = length($1);
        my @coords2;
        print "TR_CDS_START=$tr_cds_start; CDS_LENGTH=".$cds_length."nt,".(($cds_length-3)/3)."aa\n";
        #if ($transcript->seq_region_strand == 1){
          @coords2 = $transcript->cdna2genomic($tr_cds_start + $cds_length - 1, $tr_cds_start + $cds_length - 1);
        #}
        #else{
        #  @coords2 = $transcript->cdna2genomic($tr_cds_start - $cds_length + 1, $tr_cds_start - $cds_length + 1);
        #}
        my $cds_end = $coords2[0]->start;
        print "CDS_END=$cds_end\n";
        #Avoid making an NMD transcript with a CDS smaller than 35 aa
        unless (predicted_nmd_transcript($transcript, $cds_end) and ($cds_length-3)/3  <= 35){
          print $transcript->stable_id.": complete CDS of $cds_length bp $3\n";
          return ($cds_start, $cds_end, 1);
        }
      }
      #Or an open-ended ORF?
      elsif ($subseq =~ /^(ATG([ATGC]+))$/){
        my $cds_length = length($1);
        print $transcript->stable_id.": end-NF-CDS of $cds_length bp\n";
        my $cds_end = $transcript->seq_region_strand == 1 ? $transcript->seq_region_end : $transcript->seq_region_start;
        return ($cds_start, $cds_end, 0);
      }
    }
  }
  return undef;
}



=head2 predicted_nmd_transcript

 Arg[1]    : Bio::Vega::Transcript object
 Arg[2]    : integer (CDS end in genomic coordinates)
 Function  : returns true if the transcript would get a nonsense_mediated_decay biotype if a CDS was created with the given genome coordinates, that is, if there were 50 exonic bp or more between the stop codon and the last downstream splice site (according to the HAVANA rules)
 Returntype: boolean

=cut

sub predicted_nmd_transcript {
  my ($transcript, $cds_end) = @_;
  my $rank = 0; #use a separate variable instead of the exon_rank method in the Ensembl Core API as this seems to rely on exon stable ids, which will not be available until the annotation is stored in the database
  my $cds_end_exon_rank = 0;
  my $is_exonic = 0;
  my $distance = 0;
  print "Entering predicted_nmd_transcript - cds_end=$cds_end\n";
  foreach my $exon (@{$transcript->get_all_Exons}){
    $rank++; 
    print "rank".$rank.":".$exon->seq_region_start."-".$exon->seq_region_end."\n";
    #CDS end must lie in an exon
    if ($exon->seq_region_start <= $cds_end and $exon->seq_region_end >= $cds_end){
      $is_exonic = 1;
      $cds_end_exon_rank = $rank;
    }
    #Start counting distance downstream of stop codon
    if ($rank == $cds_end_exon_rank){
      $distance = $transcript->seq_region_strand == 1 ? $exon->seq_region_end - $cds_end :  $cds_end - $exon->seq_region_start;
    }
    #Keep counting to find distance between the CDS end position and the last splice site
    elsif ($is_exonic and $rank < scalar(@{$transcript->get_all_Exons})){
      $distance += $exon->length;
    }
  }
  if ($is_exonic and $cds_end_exon_rank != scalar(@{$transcript->get_all_Exons}) and $distance >= 50){
    print $transcript->stable_id.": nonsense_mediated_decay  - cds_end_exon=$cds_end_exon_rank DIST=$distance\n";
    return 1;
  }
  unless ($is_exonic){
    warn "Position $cds_end is not exonic in this transcript!\n";
  }
  return 0;
}



=head2 create_cds

 Arg[1]    : Bio::Vega::Transcript object
 Arg[2]    : integer (CDS start in genomic coordinates)
 Arg[3]    : integer (CDS end in genomic coordinates)
 Function  : Creates a Bio::Vega::Translation object associated to the transcript using the genomic coordinates provided
 Returntype: boolean

=cut

sub create_cds {
  my ($transcript, $cds_start, $cds_end) = @_;
  unless ($transcript and $cds_start and $cds_end){
    print "Not enough info for create_cds\n"; 
    return undef;
  }
print "create_cds => CDS_START: $cds_start ; CDS_END: $cds_end\n";

  my $start_exon;
  my $end_exon;
  my $seq_start;
  my $seq_end;
  #Find translation start and end exons
  foreach my $exon (@{$transcript->get_all_Exons}){
    if ($exon->seq_region_start <= $cds_start and $exon->seq_region_end >= $cds_start){
      $start_exon = $exon;
    }
    if ($exon->seq_region_start <= $cds_end and $exon->seq_region_end >= $cds_end){
      $end_exon = $exon;
    }
  }
  #Find positions within exons
  if ($transcript->seq_region_strand == 1){
    $seq_start = $cds_start - $start_exon->seq_region_start + 1;
    $seq_end = $cds_end - $end_exon->seq_region_start + 1;
  }
  else{
    ($cds_start, $cds_end) = ($cds_end, $cds_start) if $cds_start > $cds_end;
    $seq_end = $end_exon->seq_region_end - $cds_start + 1;
    $seq_start = $start_exon->seq_region_end - $cds_end + 1;
  }
  
  my $translation = Bio::Vega::Translation->new(-START_EXON => $start_exon,
                                                -END_EXON   => $end_exon,
                                                -SEQ_START  => $seq_start,
                                                -SEQ_END    => $seq_end
                                               );
  $transcript->translation($translation);


  #Assign phase and end_phase values to CDS exons
  #Assume start phase equals 0 in all cases
  calculate_exon_phases($transcript, 0);
  return 1;
}


=head2 get_host_gene_stop_codon_set

 Arg[1]    : Bio::Vega::Gene object
 Function  : returns the set of unique CDS end positions (except for cds_end_NF transcripts) in the host gene
 Returntype: arrayref of integers

=cut

sub get_host_gene_stop_codon_set {
  my $gene = shift;
  #print "G_START=".$gene->seq_region_start."; G_END=".$gene->seq_region_end."\n";
  #print "G_START=".$gene->start."; G_END=".$gene->end."\n";
  my @stop_codon_set;
  my %seen_sc;
  my @filtered_transcripts = grep {has_cds_end($_)} @{$gene->get_all_Transcripts};
  foreach my $transcript (@filtered_transcripts){
    #Ignore lingering OTT transcripts
    next unless $transcript->stable_id =~ /^ENS([A-Z]{3})?T000/;
    #Ignore "not for VEGA" transcripts unless they have a "comp_pipe" biotype or a "TAGENE_transcript" remark
    if (scalar(grep {$_->value eq "not for VEGA"} @{$transcript->get_all_Attributes('remark')})){
      next unless ($transcript->biotype eq "comp_pipe" or scalar(grep {$_->value eq "TAGENE_transcript"} @{$transcript->get_all_Attributes('remark')}));
    }
    
    #$transcript = $transcript->transform("chromosome");
    #print "T_START=".$transcript->seq_region_start."; T_END=".$transcript->seq_region_end."\n";
    my $cds_end = $transcript->seq_region_strand == 1 ? $transcript->coding_region_end : $transcript->coding_region_start;
    #print "CDS_START=$cds_start\n";
    unless ($seen_sc{$cds_end}){
      push(@stop_codon_set, $cds_end);
      $seen_sc{$cds_end} = 1;
    }
  }
  return \@stop_codon_set;
}


=head2 is_known_stop_codon

 Arg[1]    : integer (CDS end in genomic coordinates)
 Arg[2]    : integer (offset)
 Arg[3]    : Bio::Vega::Gene object (host gene)
 Function  : returns true if the stop codon ending at that position has been annotated before
 Returntype: boolean

=cut

sub is_known_stop_codon {
  my ($cds_end, $slice_offset, $gene) = @_;
  foreach my $stop_codon (@{$HOST_STOP_CODON_SET{$gene->stable_id}}){
    print "FFF=$stop_codon + $slice_offset = ".($stop_codon + $slice_offset)."; GGG=$cds_end\n";
    if ($cds_end == $stop_codon+$slice_offset){
      return 1;
    }
  }
  return 0;
}


=head2 add_end_NF_attributes

 Arg[1]    : Bio::Vega::Transcript object
 Function  : Add 'cds_end_NF' and 'mRNA_end_NF' attributes if appropriate
 Returntype: none

=cut

sub add_end_NF_attributes {
  my ($transcript, $slice_offset) = @_;
  #print "add_end_NF_attributes: ".join("==", $transcript->seq_region_start, $transcript->coding_region_start+$slice_offset, $transcript->coding_region_end+$slice_offset, $transcript->seq_region_end)."\n";
  if (($transcript->seq_region_strand == 1 and $transcript->coding_region_end+$slice_offset == $transcript->seq_region_end) or
    ($transcript->seq_region_strand == -1 and $transcript->coding_region_start+$slice_offset == $transcript->seq_region_start)
  ){
    unless ($transcript->seq->seq =~ /(TGA|TAG|TAA)$/){
      $transcript->add_Attributes(
        new Bio::EnsEMBL::Attribute(-code => 'cds_end_NF', -value => 1), 
        new Bio::EnsEMBL::Attribute(-code => 'mRNA_end_NF', -value => 1)
      );
      print "Added cds_end_NF/mRNA_end_NF attribute to transcript ".$transcript->stable_id."\n";
    }
  }
}




1;


