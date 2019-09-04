
package LoutreWrite::Default;

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
use LoutreWrite::CDSUtils qw(assign_cds_to_transcripts has_complete_cds);
use base 'Exporter';

our @EXPORT = qw( $WRITE );
our @EXPORT_OK = qw( exon_novelty can_be_merged );
our $WRITE = 0;
our $CP_BIOTYPE = "comp_pipe";

sub new {
  my $package = shift;
  return bless({}, $package);
}



=head2 parse_gxf_file

 Arg[1]    : GTF or GFF3 file
 Function  : parse gene annotation from GTF/GFF3 file - gene, transcript, exon lines
 Returntype: hash-ref representing data structure with multiple genes, their transcripts, their exons and their attributes

=cut

sub parse_gxf_file {
    my ($self, $file) = @_;
    my %data;
    if ($file =~ /\.gz$/){
        open (IN, "gunzip -dc $file |") or die "Can't open file $file\n";
    }
    else{
        open (IN, $file) or die "Can't open file $file\n";
    }
    while (<IN>){
        next if /^#/;
        chomp;
        my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split(/\t/);
        next if (!$chr or $chr =~ /_/ or $chr eq "chrM");
        my %attribs;
        foreach my $att (split(/;/, $attributes)){
            my ($name, $value);
            if ($file =~ /\.gff3(.gz)?$/){
                ($name, $value) = split(/=/, $att);
            }
            elsif ($file =~ /\.gtf(.gz)?$/){
                $att =~ s/^\s*//g;
                ($name, $value) = split(/\s+/, $att);
                $value =~ s/\"//g;
            }
            $attribs{$name} = $value;
        }
        my $gene_id = $attribs{'gene_id'};
        my $transcript_id = $attribs{'transcript_id'};
        if ($type eq "gene"){
            $data{$gene_id} = {
                                'gene_type'   => $attribs{'gene_type'},
                                'gene_name'   => $attribs{'gene_name'},
                                'gene_status' => $attribs{'gene_status'}
                                };
        }
        elsif ($type eq "transcript"){
            $data{$gene_id}{transcripts}{$transcript_id} = {
                                                            'transcript_type'   => $attribs{'transcript_type'},
                                                            'transcript_name'   => $attribs{'transcript_name'},
                                                            'transcript_status' => $attribs{'transcript_status'}
                                                            };
        }
        elsif ($type eq "exon"){
            my $exon_id = $attribs{'exon_id'} || $.; #gtf line number as default exon id
            $data{$gene_id}{transcripts}{$transcript_id}{exons}{$exon_id} = { 
                                                                                'chr'    => $chr,
                                                                                'start'  => $start,
                                                                                'end'    => $end,
                                                                                'strand' => $strand,
                                                                                   #'rank'   => $attribs{'exon_number'}
                                                                            };
            #If gtf has exon lines only, fetch gene and transcript attributes
            unless ($data{$gene_id}{'gene_name'}){
              if ($attribs{'gene_name'}){
                  $data{$gene_id}{'gene_name'} = $attribs{'gene_name'};
              }
              else{
                  $data{$gene_id}{'gene_name'} = $attribs{'gene_id'};
              }
            }
            unless ($data{$gene_id}{transcripts}{$transcript_id}{'transcript_name'}){
              if ($attribs{'transcript_name'}){
                  $data{$gene_id}{transcripts}{$transcript_id}{'transcript_name'} = $attribs{'transcript_name'};
              }
              else{
                  $data{$gene_id}{transcripts}{$transcript_id}{'transcript_name'} = $attribs{'transcript_id'};
              }
            }
                        
            #$data{$gene_id}{'gene_type'} ||= $attribs{'gene_type'};
            #$data{$gene_id}{'gene_name'} ||= $attribs{'gene_name'};
            #$data{$gene_id}{'gene_status'} ||= $attribs{'gene_status'};
            
            #$data{$gene_id}{transcripts}{$transcript_id}{'transcript_type'} ||= $attribs{'transcript_type'};
            #$data{$gene_id}{transcripts}{$transcript_id}{'transcript_name'} ||= $attribs{'transcript_name'};
            #$data{$gene_id}{transcripts}{$transcript_id}{'transcript_status'} ||= $attribs{'transcript_status'};
        }
        elsif ($type eq "CDS"){
            my $exon_id = $attribs{'exon_id'} || $.; #gtf line number as default exon id
            $data{$gene_id}{transcripts}{$transcript_id}{cds}{$exon_id} = { 
                                                                                'chr'    => $chr,
                                                                                'start'  => $start,
                                                                                'end'    => $end,
                                                                                'strand' => $strand,
                                                                                'phase'  => "p".$phase,
                                                                            };
        }
        elsif ($type eq "stop_codon"){ #Relevant for GENCODE gtf type files, where stop codon is not included in CDS
            my $exon_id = $attribs{'exon_id'} || $.; #gtf line number as default exon id
            $data{$gene_id}{transcripts}{$transcript_id}{stop_codon}{$exon_id} = { 
                                                                                'chr'    => $chr,
                                                                                'start'  => $start,
                                                                                'end'    => $end,
                                                                                'strand' => $strand,
                                                                                'phase'  => "p".$phase,
                                                                            };
        }
    }
    close (IN);

    return \%data;
}



=head2 make_vega_objects

 Arg[1]    : hash-ref of gene annotation (as provided by the parse_gxf_file method)
 Arg[2]    : Bio::Vega::DBSQL::DBAdaptor
 Arg[3]    : author name
 Arg[4]    : remark to be added to all transcripts
 Arg[5]    : force 'computationally_generated' biotype
 Arg[6]    : gene and transcript analysis logic name
 Arg[7]    : gene and transcript source
 Arg[8]    : Boolean - true if the default transcript attributes used in the comp_pipe models must be ignored
 Arg[9]    : Boolean - if true, the "not for VEGA" remark will not be added
 Arg[10]   : String - chromosome
 Function  : convert gene annotation into Vega gene, transcript and exon objects
 Returntype: list of Bio::Vega::Gene objects

=cut

sub make_vega_objects {
    my ($self, $genes, $dba, $author_name, $remark, $force_cp_biotype, $analysis_name, $source, $no_comp_pipe, $no_NFV, $only_chr) = @_;
    
    #Some common features
    $source ||= "havana";
    $analysis_name ||= "Otter";
    my $analysis = Bio::EnsEMBL::Analysis->new( -logic_name => $analysis_name );
    
    my $nfv_remark = Bio::EnsEMBL::Attribute->new(
                                -code => 'remark',
                                -value => 'not for VEGA'
                    );

    my $tagene_attrib = Bio::EnsEMBL::Attribute->new(
                                -code => 'TAGENE_transcript',
                                -value => '1'
                        );
    my $tagene_remark = Bio::EnsEMBL::Attribute->new(
                                -code => 'remark',
                                -value => 'TAGENE_transcript'
                        );

    my $source_remark_string;
    if ($remark){
        $source_remark_string = $remark;
    }
    else{
        $source_remark_string = 'Assembled from SLR-seq (SRP049776) and PacBio Capture-seq (GSE93848) reads';
    }
    my $source_remark = Bio::EnsEMBL::Attribute->new(
                                -code => 'remark',
                                -value => $source_remark_string
                        );
                
    my $author = Bio::Vega::Author->new(
                                -name   => $author_name,
                                -email  => $author_name."\@ebi.ac.uk"
                 );

    my $sa = $dba->get_SliceAdaptor();
    my %slices;
    my @objects;

    my %genes = %{$genes};
    GENE:foreach my $gid (keys %genes){
        my $gene = Bio::Vega::Gene->new(
                                        -stable_id  => ($gid =~ /^OTT/ ? $gid : ""),
                                        -biotype    => ($force_cp_biotype ? $CP_BIOTYPE : $genes{$gid}{'gene_type'}),
                                        -analysis   => $analysis,
                                        -source     => $source,
                                        -status     => $genes{$gid}{'gene_status'} || "UNKNOWN"
                                        );
        #$gene->add_Attributes(new Bio::EnsEMBL::Attribute(-code => 'name', -value => $genes{$gid}{'gene_name'})); #Let the script generate a name if needed
        $gene->gene_author($author);
        unless ($no_comp_pipe){
            $gene->add_Attributes($source_remark);
        }
        unless ($no_NFV){
            $gene->add_Attributes($nfv_remark);
        }
        #$gene->add_Attributes($tagene_attrib);
        $gene->add_Attributes(Bio::EnsEMBL::Attribute->new(-code => 'hidden_remark', -value => $genes{$gid}{'gene_name'}));
        
        #Add biotype-status combination sanity check!!!

        #Make transcript objects
        foreach my $tid (keys %{$genes{$gid}{transcripts}}){
            my $transcript = Bio::Vega::Transcript->new(
                                                        -stable_id  => ($tid =~ /^OTT/ ? $tid : ""),
                                                        -biotype    => ($force_cp_biotype ? $CP_BIOTYPE : $genes{$gid}{transcripts}{$tid}{'transcript_type'}),
                                                        -analysis   => $analysis,
                                                        -source     => $source,
                                                        -status     => $genes{$gid}{transcripts}{$tid}{'transcript_status'} || "UNKNOWN"
                                                        );
            #$transcript->add_Attributes(new Bio::EnsEMBL::Attribute(-code => 'name', -value => $genes{$gid}{transcripts}{$tid}{'transcript_name'})); #Let the script generate a name if needed
            $transcript->transcript_author($author);
            unless ($no_comp_pipe){
              $transcript->add_Attributes($source_remark);
            }
            unless ($no_NFV){
                $transcript->add_Attributes($nfv_remark);
            }
            $transcript->add_Attributes(Bio::EnsEMBL::Attribute->new(-code => 'hidden_remark', 
                                                                     -value => "ID: ".$genes{$gid}{transcripts}{$tid}{'transcript_name'}));

            if ($source_remark){
              $transcript->add_Attributes($source_remark);
            }
            #$transcript->add_Attributes($tagene_attrib);
            $transcript->add_Attributes($tagene_remark);

            #Make exon objects
            foreach my $exid (keys %{$genes{$gid}{transcripts}{$tid}{exons}}){
                my $exon = Bio::Vega::Exon->new(
                                                -start   => $genes{$gid}{transcripts}{$tid}{exons}{$exid}{'start'},
                                                -end     => $genes{$gid}{transcripts}{$tid}{exons}{$exid}{'end'},
                                                -strand  => ($genes{$gid}{transcripts}{$tid}{exons}{$exid}{'strand'} eq "+" ? 1 : -1),
                                                );
                my $chr = $genes{$gid}{transcripts}{$tid}{exons}{$exid}{'chr'};
                $chr =~ s/^chr//;
                next GENE if ($only_chr and $chr ne $only_chr);
                if (!($slices{$chr})){
                    my $slice = $sa->fetch_by_region("chromosome", "chr$chr-38");
                    $slices{$chr} = $slice;
                }
                $exon->slice($slices{$chr});
                $exon->phase(-1);
                $exon->end_phase(-1);
                
                #Assign exon phases from gxf file if exon starts within CDS
                #Compare also with stop_codon coordinates as stop codon is not included in CDS in GENCODE gtf files
                my $gtf_phase = $genes{$gid}{transcripts}{$tid}{cds}{$exid}{'phase'} || $genes{$gid}{transcripts}{$tid}{stop_codon}{$exid}{'phase'};
                print "\n$gtf_phase \n" if $gtf_phase;
                if ($gtf_phase and ($gtf_phase eq "p0" or $gtf_phase eq "p1" or $gtf_phase eq "p2")){
                #if ($gtf_phase){
                  print "\nA\n";
                  if ((($exon->start == $genes{$gid}{transcripts}{$tid}{cds}{$exid}{'start'} or $exon->start == $genes{$gid}{transcripts}{$tid}{stop_codon}{$exid}{'start'}) and $exon->strand == 1) or 
                      (($exon->end == $genes{$gid}{transcripts}{$tid}{cds}{$exid}{'end'} or $exon->end == $genes{$gid}{transcripts}{$tid}{stop_codon}{$exid}{'end'}) and $exon->strand == -1)
                     ){
                    $exon->phase(0) if $gtf_phase eq "p0";
                    $exon->phase(1) if $gtf_phase eq "p2";
                    $exon->phase(2) if $gtf_phase eq "p1"; print "\nB".$exon->phase."\n";
                    if ((($exon->end == $genes{$gid}{transcripts}{$tid}{cds}{$exid}{'end'} or $exon->end == $genes{$gid}{transcripts}{$tid}{stop_codon}{$exid}{'end'}) and $exon->strand == 1) or 
                        (($exon->start == $genes{$gid}{transcripts}{$tid}{cds}{$exid}{'start'} or $exon->start == $genes{$gid}{transcripts}{$tid}{stop_codon}{$exid}{'start'}) and $exon->strand == -1)
                       ){
                       $exon->end_phase(($exon->phase + $exon->length) % 3); print "\nC".$exon->end_phase."\n";
                    }
                  }
                  else{
                    #Half-coding start codon: no phase but end_phase
                    if ($exon->end == $genes{$gid}{transcripts}{$tid}{cds}{$exid}{'end'} and $exon->strand == 1){
                      $exon->end_phase(($exon->end - $genes{$gid}{transcripts}{$tid}{cds}{$exid}{'start'} + 1) % 3); print "\nD".$exon->end_phase."\n";
                    }
                    elsif ($exon->start == $genes{$gid}{transcripts}{$tid}{cds}{$exid}{'start'} and $exon->strand == -1){
                      $exon->end_phase(($genes{$gid}{transcripts}{$tid}{cds}{$exid}{'end'} - $exon->start + 1) % 3); print "\nE".$exon->end_phase."\n";
                    }
                  }
                }
                $transcript->add_Exon($exon);
            }
            
            #Make translation object
            my $cds_start;
            my $cds_end;
            foreach my $exid (keys %{$genes{$gid}{transcripts}{$tid}{cds}}){
                if (!($cds_start) or $cds_start > $genes{$gid}{transcripts}{$tid}{cds}{$exid}{'start'}){
                    $cds_start = $genes{$gid}{transcripts}{$tid}{cds}{$exid}{'start'};
                }
                if (!($cds_end) or $cds_end < $genes{$gid}{transcripts}{$tid}{cds}{$exid}{'end'}){
                    $cds_end = $genes{$gid}{transcripts}{$tid}{cds}{$exid}{'end'};
                }
            }
            foreach my $exid (keys %{$genes{$gid}{transcripts}{$tid}{stop_codon}}){
                if ($cds_end and $transcript->strand==1){
                  if ($cds_end < $genes{$gid}{transcripts}{$tid}{stop_codon}{$exid}{'end'}){
                    $cds_end = $genes{$gid}{transcripts}{$tid}{stop_codon}{$exid}{'end'};
                  }
                }
                elsif ($cds_start and $transcript->strand==-1){
                  if ($cds_start > $genes{$gid}{transcripts}{$tid}{stop_codon}{$exid}{'start'}){
                    $cds_start = $genes{$gid}{transcripts}{$tid}{stop_codon}{$exid}{'start'};
                  }
                }
            }
            if ($cds_start and $cds_end){
                my $start_exon;
                my $end_exon;
                my $seq_start;
                my $seq_end;
                foreach my $exon (@{$transcript->get_all_Exons}){
                  if ($exon->start <= $cds_start and $exon->end >= $cds_start){
                    if ($transcript->strand==1){
                      $start_exon = $exon;
                      $seq_start = $cds_start - $start_exon->start + 1;
                    }
                    else{
                      $end_exon = $exon;
                      $seq_end = $end_exon->end - $cds_start + 1;
                    }
                  }
                  if ($exon->start <= $cds_end and $exon->end >= $cds_end){
                    if ($transcript->strand==1){
                      $end_exon = $exon;
                      $seq_end = $cds_end - $end_exon->start + 1;
                    }
                    else{
                      $start_exon = $exon;
                      $seq_start = $start_exon->end - $cds_end + 1;
                    }
                  }
                }
               
                my $translation = Bio::Vega::Translation->new(
                                                              -start_exon => $start_exon,
                                                              -end_exon   => $end_exon,
                                                              -seq_start  => $seq_start,
                                                              -seq_end    => $seq_end
                                                              );
                $transcript->translation($translation);
                
                #Re-check exon phases
                my $check_message = $self->check_exon_phases($transcript);
                print $check_message."\n";
            }

            #Add transcript to gene
            $gene->add_Transcript($transcript);
        }
        
        

        #Update gene and transcript biotypes
        if (!($gene->biotype) or $gene->biotype eq "missing_biotype"){
            $gene = $self->assign_biotypes($gene); #This also updates description if appropriate 
        }
        
        #Update gene and transcript status - only lncRNA genes for now
        $gene = $self->assign_status($gene);
        
        push(@objects, $gene);
    }

    return \@objects;
}



=head2 process_gene

 Arg[1]    : Vega gene object to be written
 Arg[2]    : mode - either 'add' or 'update'
 Arg[3]    : Otter dataset name ('human', 'mouse', 'human_test', etc)
 Arg[4]    : Bio::Vega::DBSQL::DBAdaptor
 Arg[5]    : Boolean - if true, do not check for completely or partially identical intron chains in the annotation
 Function  : get region with the gene coordinates; 'add' the gene to the region or 'update' pre-existing gene annotation; store changes in the database 
 Returntype: String (message with the outcome)

=cut

sub process_gene {
    my ($self, $gene, $mode, $dataset_name, $dba, $no_intron_check) = @_;

    #If in "update" mode, the region has to span the whole host gene to stop this from being "spliced out" as an incomplete gene. 
    #This would prevent the new annotation from being added to the host gene.
    my ($padded_start, $padded_end);
    if ($mode eq "update"){
        my $host_gene = $dba->get_GeneAdaptor->fetch_by_stable_id($gene->stable_id);
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
        cs      => "chromosome",
        csver   => "Otter"
    );

    my $region_action = Bio::Otter::ServerAction::Script::Region->new_with_slice($local_server);
    my $region = $region_action->get_region; #Bio::Vega::Region
    print "\nFOUND ".scalar($region->genes)." GENES IN REGION ".$gene->seq_region_name.":".$gene->start."\n";

    #Reset gene coordinates to the start of the region slice,
    #otherwise the gene coordinates in loutre will be wrong
    my $slice_offset = $region->slice->start - 1;

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
            print "TR: $id: Added transcript $new_tr_name to novel gene $new_gene_name\n";
        }
        #Add new gene to region
        $region->add_genes($gene);
    }

    elsif ($mode eq "update"){
        #Update existing gene: host gene that provided its stable id to the gene being stored
        if (scalar(grep {$_->stable_id eq $gene->stable_id} $region->genes) == 1){
            my ($db_gene) = grep {$_->stable_id eq $gene->stable_id} $region->genes;
            #There are two options now:
            my $submode = "cons";
            #A) CONSERVATIVE: a new transcript model is only stored if it has a novel intron structure or an existing intron structure with novel sequence in its terminal exons; existing transcripts are never updated.
            if ($submode eq "cons"){
                #Fetch all intron structures in the gene
                my %db_gene_intr_str;
                foreach my $db_tr (@{$db_gene->get_all_Transcripts}){
                    my $db_intron_str = ":".join(':', map {$_->start."-".$_->end} @{$db_tr->get_all_Introns}).":";
                    $db_gene_intr_str{$db_intron_str} ++;
                }
                #Compare each new transcript with existing transcripts in database
                foreach my $tr (@{$gene->get_all_Transcripts}){
                    my $add_transcript = 0;
                    #Quick intron structure check
                    #Add transcript to database gene if intron structure is novel and not included in any database transcript
                    my $full_intron_chain_in_annot = 0;
                    my $partial_intron_chain_in_annot = 0;
                    my @introns = @{$tr->get_all_Introns};
                    my $intron_str = ":".join(':', map {$_->start."-".$_->end} @introns).":";
                    foreach my $db_intron_str (keys %db_gene_intr_str){
                        if ($db_intron_str eq $intron_str){
                            $full_intron_chain_in_annot = 1;
                            print "FULL_INTRON_CHAIN: $db_intron_str  vs  $intron_str\n";
                        }
                        elsif ($db_intron_str =~ /$intron_str/){
                            $partial_intron_chain_in_annot = 1;
                            print "PARTIAL_INTRON_CHAIN: $db_intron_str  vs  $intron_str\n";
                        }
                    }
                    
                    if ($full_intron_chain_in_annot){
                      $add_transcript = 0;
                    }
                    #If intron chain found as part of longer intron chain in annotation, 
                    #check for terminal exon extensions that may make this model unique
                    elsif ($partial_intron_chain_in_annot){
                        #Look at first and last introns only
                        my @t_introns = ($introns[0], $introns[-1]);
                        #Loop through ALL database transcripts 
                        foreach my $db_tr (@{$db_gene->get_all_Transcripts}){
                            my @db_introns = @{$db_tr->get_all_Introns};
                            my $db_intron_str = ":".join(':', map {$_->start."-".$_->end} @db_introns).":";
                            #Only look at cases where the intron chain is included in a longer one in the database
                            #If the chains are identical, the 5' and 3' ends must be compared (LATER, IN TO-DO LIST)
                            if ($db_intron_str =~ /$intron_str/ and $db_intron_str ne $intron_str){
                                my $exon_extension = 0;
                                #Find shared introns, compare flanking exons
                                foreach my $intron (@t_introns){
                                    my ($db_intron) = grep {$_->start == $intron->start and $_->end == $intron->end} @db_introns;
                                    my $prev_exon = $intron->prev_Exon;
                                    my $prev_db_exon = $db_intron->prev_Exon;
                                    my $next_exon = $intron->next_Exon;
                                    my $next_db_exon = $db_intron->next_Exon;
                                    my $mee = 3; #minimum exon extension into the db intron to be considered for novelty
                                    
                                    #Do not compare 5' and 3' transcript ends now (LATER, IN TO-DO LIST)
                                    if ($tr->seq_region_strand == 1){
                                        #db: 5'#####----####-----######3'
                                        #tr:        5'######-----####3'
                                        if ($prev_exon->start <= $prev_db_exon->start - $mee){
                                            unless ($prev_db_exon->start == $db_tr->start){
                                                $exon_extension = 1;
                                            }
                                        }
                                        #db: 5'#####----####-----######3'
                                        #tr:  5'####----######3'
                                        if ($next_exon->end >= $next_db_exon->end + $mee){
                                            unless ($next_db_exon->end == $db_tr->end){
                                                $exon_extension = 1;
                                            }
                                        }
                                    }
                                    elsif ($tr->seq_region_strand == -1){
                                        #db: 3'#####----####-----######5'
                                        #tr:  3'####----######5'
                                        if ($prev_exon->end >= $prev_db_exon->end + $mee){
                                            unless ($prev_db_exon->end == $db_tr->end){
                                                $exon_extension = 1;
                                            }
                                        }
                                        #db: 3'#####----####-----######5'
                                        #tr:        3'######-----####5'
                                        if ($next_exon->start <= $next_db_exon->start - $mee){
                                            unless ($next_db_exon->start == $db_tr->start){
                                                $exon_extension = 1;
                                            }
                                        }
                                    }
                                }
                                if ($exon_extension){
                                  $add_transcript = 1; #Provisionally, not final until all db transcripts checked
                                }
                                else{
                                  #Nothing unique in the model w.r.t. the database transcript: discard
                                  $add_transcript = 0;
                                  last;
                                }
                            }
                        }
                    }
                    #If intron chain not in annotation, annotate transcript
                    else{
                        #Set flag to indicate that the transcript must be annotated
                        $add_transcript = 1;
                    }
                    
                    #If 'no_intron_check' set, override intron check results
                    if ($no_intron_check){
                      $add_transcript = 1;
                    }
                    
                    my $id = $tr->get_all_Attributes('hidden_remark')->[0]->value;
                    if ($add_transcript == 1){
                        #Create a name for the new transcript
                        my $new_tr_name = $self->get_new_transcript_name($db_gene, $dba);
                        my $name_att = Bio::EnsEMBL::Attribute->new(-code => 'name', -value => $new_tr_name);
                        $tr->add_Attributes($name_att);
                        $tr->slice($region->slice);
                        $db_gene->add_Transcript($tr);
                        print "TR: $id: Added transcript $new_tr_name to gene ".$db_gene->stable_id."\n";
                    }
                    else{
                        if ($full_intron_chain_in_annot or $partial_intron_chain_in_annot){
                          print "TR: $id: Rejected transcript in host gene ".$db_gene->stable_id." as intron chain exists\n";
                        }
                    }
                }
            }

            #B) AGGRESSIVE: existing transcripts can be modified by merging them with new models if they are compatible.
            #   Also, allows replacing existing transcript models with new ones if they have the same stable ids.
                #NOTE: UNCLEAR WHEN THIS OPTION WOULD BE USED - 
                #      JUST WANT TO KEEP A RECORD OF THE CODE THAT UPDATES EXON COORDINATES WHILE PRESERVING ATTRIBUTES AND EVIDENCE
            elsif ($submode eq "db"){
                foreach my $tr (@{$gene->get_all_Transcripts}){
                    #If a transcript with the same stable id exists, update it if there are any changes to exon coordinates, CDS or biotype
                    my ($db_tr) = grep {$_->stable_id eq $tr->stable_id} @{$db_gene->get_all_Transcripts};
                    if ($db_tr){
                        #Compare exon coordinates between transcripts
                        if (join('+', map {$_->vega_hashkey} @{$tr->get_all_Exons}) ne join('+', map {$_->vega_hashkey} @{$db_tr->get_all_Exons})){
                            #Update transcript
                            $db_tr->flush_Exons;     #This empties the exon list. If the same exons are added again, 
                                                    #those that are not associated with other transcripts will get a new stable id
                            foreach my $exon (@{$tr->get_all_Exons}){
                                $db_tr->add_Exon($exon);
                            }
                        }
                        
                        #Compare CDS
                        if (($tr->translation and !($db_tr->translation)) or
                            ($tr->translation and $db_tr->translation and $tr->coding_region_start ne $db_tr->coding_region_start) or 
                            ($tr->translation and $db_tr->translation and $tr->coding_region_end ne $db_tr->coding_region_end)){
                            foreach my $db_exon (@{$db_tr->get_all_Exons}){
                                #set translation start exon
                                if ($db_exon->start == $tr->translation->start_Exon->start){
                                    $db_exon->translation->start_Exon($db_exon);
                                    $db_tr->translation->start($tr->translation->start); #position within the start exon
                                }
                                #set translation end exon
                                if ($db_exon->start == $tr->translation->end_Exon->start){
                                    $db_exon->translation->end_Exon($db_exon);
                                    $db_tr->translation->end($tr->translation->end); #position within the end exon
                                }
                            }
                        }
                        
                        #Update transcript biotype
                        if ($tr->biotype ne $db_tr->biotype){
                            $db_tr->biotype($tr->biotype);
                        }
                    }
                
                    #Else, add new transcript to gene 
                    #TO DO: CHECK THAT IT WON'T CREATE A TRANSCRIPT DUPLICATION
                    else{
                        #Create a name for the new transcript
                        my $new_tr_name = $self->get_new_transcript_name($db_gene, $dba);
                        my $name_att = Bio::EnsEMBL::Attribute->new(-code => 'name', -value => $new_tr_name);
                        $tr->add_Attributes($name_att);
                        $tr->slice($region->slice);
                        $db_gene->add_Transcript($tr);
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
    #if ($write){
        $local_server->authorized_user($gene->gene_author->name); # preserve authorship
        $g_msg = write_gene_region($region_action, $region);
        print $g_msg."\n";
    #}

    return $g_msg;
}




=head2 process_gene_2

 Arg[1]    : Vega gene object to be written
 Arg[2]    : mode - either 'add' or 'update'
 Arg[3]    : submode - either 'cons' or 'aggr'
 Arg[4]    : Otter dataset name ('human', 'mouse', 'human_test', etc)
 Arg[5]    : Bio::Vega::DBSQL::DBAdaptor
 Arg[6]    : Boolean - if true, do not check for completely or partially identical intron chains in the annotation
 Function  : get region with the gene coordinates; 'add' the gene to the region or 'update' pre-existing gene annotation; store changes in the database 
 Returntype: String (message with the outcome)

=cut

sub process_gene_2 {
    my ($self, $gene, $mode, $submode, $dataset_name, $dba, $no_intron_check) = @_;
    my @log;
print "SUBMODE: $submode\n";
    #If in "update" mode, the region has to span the whole host gene to stop this from being "spliced out" as an incomplete gene. 
    #That would prevent the new annotation from being added to the host gene.
    my ($padded_start, $padded_end);
    if ($mode eq "update"){
        my $host_gene = $dba->get_GeneAdaptor->fetch_by_stable_id($gene->stable_id);
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
                foreach my $db_tr (@{$db_gene->get_all_Transcripts}){
                    #Exclude transcripts not yet stored in the database
                    next unless $db_tr->stable_id;
                    #Exclude artifacts
                    next if $db_tr->biotype eq "artifact";
                    #Ignore "not for VEGA" transcripts unless they have a "comp_pipe" biotype or a "TAGENE_transcript" remark
                    if (scalar(grep {$_->value eq "not for VEGA"} @{$db_tr->get_all_Attributes('remark')})){
                      next unless ($db_tr->biotype eq "comp_pipe" or scalar(grep {$_->value eq "TAGENE_transcript"} @{$db_tr->get_all_Attributes('remark')}));
                    }
                    #Exclude transcripts having a "comp_pipe_rejected" remark
                    next if scalar(grep {$_->value eq "comp_pipe_rejected"} @{$db_tr->get_all_Attributes('remark')});
                    #Ignore single-exon transcripts (?)
                    next if scalar @{$db_tr->get_all_Introns} == 0;  #Do not merge with single-exon transcripts? (In case they are CLS Platinum...)
                    my $i = intron_novelty($tr, $db_tr);
                    my $e = exon_novelty($tr, $db_tr);
                    my $m = can_be_merged($tr, $db_tr);  
                    $tr_comp{$db_tr->stable_id} = {'intron' => $i, 'exon' => $e, 'merge' => $m};
                    print "COMP: i=$i, e=$e, m=$m ".$db_tr->stable_id." ".$db_tr->biotype."\n";
                }
                #Multi-exon transcripts
                if (scalar @{$tr->get_all_Exons} > 1){
                    #Conservative mode: a new transcript model is only stored if it has a novel intron structure or an existing intron structure with novel sequence in its terminal exons; existing transcripts are never updated.
                    if ($submode eq "cons"){
                        #Add transcript if intron chain shows novelty w.r.t. all existing transcripts
                        #Discard transcript if an existing transcript has identical intron chain
                        DBTR:foreach my $db_tr_id (keys %tr_comp){
                            if (($tr_comp{$db_tr_id}->{'intron'} == 1) or 
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
                                    my $db_tr = $dba->get_TranscriptAdaptor->fetch_by_stable_id($db_tr_id);
                                    #Do not merge with full-length CDS transcripts (TO DO)
                                    #unless ($db_tr->translate and has_complete_cds($db_tr)){ 
                                    #Do not merge with coding transcripts
                                    #Skip transcript if it could be merged with a coding transcript,
                                    # so as not to make a partially redundant transcript
                                    if ($db_tr->translate()){
                                      print "Skipping transcript as partially redundant with a coding transcript\n";
                                      $add_transcript = 0;
                                      @merge_candidates = ();
                                      last DBTR;
                                    }
                                    else{
                                      push(@merge_candidates, $db_tr);
                                    }
                                }
                                else{
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
                            if (scalar(@{$db_tr->get_all_Exons}) == 1 and !($db_tr->translate())){
                                if ($tr->start <= $db_tr->start and $tr->end => $db_tr->end and
                                    !($tr->start == $db_tr->start and $tr->end == $db_tr->end)){
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
                    my $merged_transcript = merge_transcripts($tr, $sel_db_tr );
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
                    foreach my $g ($region->genes) {
                        foreach my $ts (@{$g->get_all_Transcripts}) {
                            if ($ts->stable_id eq $merged_transcript->stable_id){
                                $ts->flush_Exons; #This empties the exon list. If the same exons are added again, 
                                                  #those that are not associated with other transcripts will get a new stable id
                                foreach my $exon (@{$merged_transcript->get_all_Exons}){
                                    $ts->add_Exon($exon);
                                }
                
                                #Add remarks (except 'not for VEGA') and hidden remarks from the novel transcript
                                foreach my $att (@{$merged_transcript->get_all_Attributes}){
                                    #COMMENTING THIS OUT AS WE MAY WANT TO KEEP THE NOT FOR VEGA REMARK (NEEDS TO READ THE $no_NFV VARIABLE) - TO DO
                                    #next if ($att->code eq "remark" and $att->value eq "not for VEGA");
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
                                #If 'comp_pipe' transcript, change biotype and remove 'not for VEGA' attribute
                                #unless 'comp_pipe_rejected' remark
                                if ($ts->biotype eq "comp_pipe" and scalar(grep {$_->value eq "comp_pipe_rejected"} @{$ts->get_all_Attributes('remark')}) == 0){
                                    my $t_biotype = get_transcript_biotype($db_gene);
                                    $ts->biotype($t_biotype);
                                    if (scalar(grep {$_->value eq "not for VEGA"} @{$ts->get_all_Attributes('remark')})){
                                       my $array = $ts->{'attributes'};
                                       @$array = grep { $_->value ne "not for VEGA" } @$array;
                                    }
                                }

                                #If coding gene, try to assign a CDS and change the biotype accordingly
                                assign_cds_to_transcripts($ts, $g, $slice_offset);
                                
                                #Change authorship
                                $ts->transcript_author($merged_transcript->transcript_author);
                                
                                print "TR: $id: Will extend transcript ".$ts->stable_id." (".$ts->biotype.") in gene ".$g->stable_id."\n";
                                push (@log, "TR2: $id: Extended transcript ".$ts->stable_id." (".$ts->biotype.") in gene ".$g->stable_id." (".$g->biotype.")");
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
                    my $new_tr_name = $self->get_new_transcript_name($db_gene, $dba);
                    my $name_att = Bio::EnsEMBL::Attribute->new(-code => 'name', -value => $new_tr_name);
                    $tr->add_Attributes($name_att);
                    print "TR0_START=".$tr->seq_region_start."; TR0_END=".$tr->seq_region_end."\n";
                    $tr->slice($region->slice);
                    $tr->start($tr->start - $slice_offset);
                    $tr->end($tr->end - $slice_offset);
                    #foreach my $exon (@{$tr->get_all_Exons}){
                    #    $exon->start($exon->start - $slice_offset);
                    #    $exon->end(  $exon->end   - $slice_offset);
                    #    $exon->slice($region->slice);
                    #}
                    #unless $use_comp_pipe_biotype, assign transcript biotype based on other transcripts
                    my $use_comp_pipe_biotype = 0; #PASS THIS AS A PARAMETER!!!
                    unless ($use_comp_pipe_biotype){
                        my $t_biotype = get_transcript_biotype($db_gene);
                        $tr->biotype($t_biotype);
                    }
                    #If coding gene, try to assign a CDS and change the biotype accordingly
                    print "TR1_START=".$tr->seq_region_start."; TR1_END=".$tr->seq_region_end."\n";
                    assign_cds_to_transcripts($tr, $db_gene, $slice_offset);
                    $db_gene->add_Transcript($tr);
                    print "TR: $id: Will add transcript $new_tr_name (".$tr->biotype.") to gene ".$db_gene->stable_id."\n";
                    push (@log, "TR2: $id: Added transcript $new_tr_name (".$tr->biotype.") to gene ".$db_gene->stable_id." (".$db_gene->biotype.")");
                }
                else{
                    print "TR: $id: Will reject transcript in host gene ".$db_gene->stable_id." as intron chain exists\n";
                    push (@log, "TR2: $id: Rejected transcript in host gene ".$db_gene->stable_id." as intron chain exists");
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
        while (($g_msg eq " " or $g_msg =~ /write failed/) and $n<10){
            $g_msg = write_gene_region($region_action, $region);
            sleep(int(rand(3)));
            $n++;
        }
        print $g_msg."\n";
    }

    return ($g_msg, \@log);
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
    
    #In Otter/Zmap, the gene name is created by selecting 'New' (Ctrl+N) in the session window
    #The code can be found in ensembl-otter/modules/MenuCanvasWindow/SessionWindow.pm (subs_edit_new_subsequence and _region_name_and_next_locus_number)
    #Thus, the name has been created when the gene attributes are dumped in the XML and SQLite files

    #Find the clone that overlaps the gene's 3' end - according to ensembl-otter/modules/MenuCanvasWindow/SessionWindow.pm (sub _region_name_and_next_locus_number)
    my $most_3prime = $gene->strand == 1 ? $gene->end : $gene->start;
    
#    my $clone_name;
#    my $tmp_slice = $dba->get_SliceAdaptor->fetch_by_region("toplevel", $gene->seq_region_name, $most_3prime, $most_3prime);
#    if ($tmp_slice->project("contig")->[0]){
#        $clone_name = $tmp_slice->project("contig")->[0]->to_Slice->seq_region_name;
#    }
    
    my $clone_name;
    my $actual_cs;
    foreach my $cs ($region->fetch_CloneSequences){ #Bio::Vega::Region method - gets clones trimmed down to the region boundaries
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
            print "AAA ".$gene->stable_id."  ".$gene->get_all_Attributes('name')->[0]->value."\n\n";
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

    my $gene_name = $gene->get_all_Attributes('name')->[0]->value; 
    #The gene name attribute may be a gene symbol rather than a clone-based name
    #Use it only as a defaut name (eg. for new genes)
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



=head2 assign_biotypes

 Arg[1]    : Bio::Vega::Gene object
 Function  : Assign biotypes to gene and transcripts
 Returntype: Bio::Vega::Gene object

=cut

sub assign_biotypes {
    my ($self, $gene) = @_;
    
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
        if ($gene->length >= 200){
            #Get gene slice and overlapping genes
            my $gene_slice = $gene->slice->adaptor->fetch_by_region("toplevel", $gene->seq_region_name, $gene->start, $gene->end);
            my @other_genes = grep {$_->stable_id ne $gene->stable_id} @{$gene_slice->get_all_Genes};
            
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
                            $intronic_genes{$other_gene->get_all_Attributes('name')->[0]->value} = 1;
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
                        $in_intron_genes{$other_gene->get_all_Attributes('name')->[0]->value} = 1;
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
                    $other_strand_genes{$other_gene->get_all_Attributes('name')->[0]->value} = 1;
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
                    $gene->description("novel transcript, antisense to ".join("and ", keys %other_strand_genes));                
                }
                elsif ($gene->biotype eq "sense_overlapping"){
                    $gene->description("novel transcript, sense overlapping ".join("and ", keys %intronic_genes));                
                }    
                elsif ($gene->biotype eq "sense_intronic"){
                    $gene->description("novel transcript, sense intronic to ".join("and ", keys %in_intron_genes));                
                }
            }
        }
        else{
            #What if gene length < 200bp ?
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



=head2 add_sources

 Arg[1]    : list of Bio::Vega::Gene objects
 Arg[2]    : source info file (output of extract_pasa_assembly_info.pl)
 Function  : Store Otter column name as transcript remark and read names as transcript hidden remark
 Returntype: list of Bio::Vega::Gene objects

=cut

sub add_sources {
    my ($self, $genes, $source_file) = @_;
    
    #Read source info
    my %sources;
    if ($source_file =~ /\.gz$/){
        open (IN, "gunzip -dc $source_file |") or die "Can't open file $source_file\n";
    }
    else{
        open (IN, $source_file) or die "Can't open file $source_file\n";
    }
    while (<IN>){
        chomp;
        my ($transc_id, $source_info) = split(/\t/);
        $sources{$transc_id} = $source_info;    
    }
    close (IN);
    
    foreach my $gene (@$genes){
        foreach my $transcript (@{$gene->get_all_Transcripts}){
            #Find original transcript id that matches that in the source file
            my $transcript_id;
            if (my ($t_name) = map {$_->value} grep {$_->value =~ /^ID:/} @{$transcript->get_all_Attributes('hidden_remark')}){
                $t_name =~ s/ID: //;
                if ($sources{$t_name}){
                    $transcript_id = $t_name;
                }
            }
            unless ($transcript_id){ #Alternatively, try current stable id
                if ($sources{$transcript->stable_id}){
                    $transcript_id = $transcript->stable_id;
                }
            }
            if ($transcript_id){    
                $transcript->add_Attributes(new Bio::EnsEMBL::Attribute(-code => 'hidden_remark', 
                                                                        -value => $sources{$transcript_id}));
            }
        }
    }
    
    return $genes;
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
    my $gene_slice = $gene->slice->adaptor->fetch_by_region("toplevel", $gene->seq_region_name, $gene->start, $gene->end);
    my @db_genes = grep {$_->stable_id ne $gene->stable_id and $_->seq_region_strand == $gene->strand} @{$gene_slice->get_all_Genes};

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

    #Fetch database genes overlapping our gene
    my $gene_slice = $gene->slice->adaptor->fetch_by_region("toplevel", $gene->seq_region_name, $gene->start, $gene->end);
    my @db_genes = grep {$_->stable_id ne $gene->stable_id and $_->seq_region_strand == $gene->strand} @{$gene_slice->get_all_Genes};

    my %host_by_transcript;
    foreach my $transcript (@{$gene->get_all_Transcripts}){
        my $host_gene;
        #Check intron match - find gene with the most matching introns
        my @candidates;
        my $max_n_introns = 0;
        foreach my $db_gene (@db_genes){
            #Exclude undesired host genes
            if ($db_gene->source ne "havana" or
                $db_gene->biotype eq "artifact" or
                scalar (grep {$_->value eq "not for VEGA"} @{$db_gene->get_all_Attributes('remark')}) == 1
            ){
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
                    if ($db_intron->seq_region_start == $intron->seq_region_start and $db_intron->seq_region_end == $intron->seq_region_end){
                        $n_introns++; #print $db_gene->stable_id.": ".$db_intron->seq_region_start."-".$db_intron->seq_region_end."\n" if $transcript->start == 38024151;
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
                        my $distance = max($transcript->seq_region_start >= $db_gene->seq_region_start) - min($transcript->seq_region_end >= $db_gene->seq_region_end) + 1;
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
            my $tr_slice = $transcript->slice->adaptor->fetch_by_region("toplevel", $transcript->seq_region_name, $transcript->start, $transcript->end);
            my @db_genes = grep {$_->seq_region_strand == $gene->seq_region_strand} @{$tr_slice->get_all_Genes};
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
            my $tr_slice = $transcript->slice->adaptor->fetch_by_region("toplevel", $transcript->seq_region_name, $transcript->start, $transcript->end);
            my @db_genes = grep {$_->seq_region_strand == $gene->seq_region_strand} @{$tr_slice->get_all_Genes};
            my $overlapped_genes_count = 0;
            DBG:foreach my $db_gene (@db_genes){
                next if $db_gene->source ne "havana";
                next if $db_gene->biotype eq "artifact";
                next if scalar(grep {$_->value eq "not for VEGA"} @{$db_gene->get_all_Attributes('remark')}) > 0;
                foreach my $db_exon (@{$db_gene->get_all_Exons}){
                    foreach my $exon (@{$transcript->get_all_Exons}){
                        if ($db_exon->seq_region_start <= $exon->seq_region_end and $db_exon->seq_region_end >= $exon->seq_region_start){
                            $overlapped_genes_count++; 
                            $message .= "OV: ".$db_gene->stable_id." ".$db_gene->biotype."  ";
                            next DBG;
                        }
                    }
                }
            }
            if ($overlapped_genes_count > $max_ov_loci){
                my $t_name = $transcript->stable_id || $transcript->get_all_Attributes('hidden_remark')->[0]->value;
                print "KILL_2: ".$t_name."  ".$transcript->start."-".$transcript->end." overlaps $overlapped_genes_count loci: $message\n";
                $list{$t_name} = 1;
            }
        }
    }

    return \%list;
}



=head2 check_exon_phases

 Arg[1]    : Bio::Vega::Transcript object
 Function  : 
 Returntype: String

=cut

sub check_exon_phases {
    my ($self, $transcript) = @_;
    my $message = " ";
    if ($transcript->translation and 
        scalar(grep {$_->value == 1} @{$transcript->get_all_Attributes('CDS_start_NF')}) == 0){
        foreach my $exon (@{$transcript->get_all_Exons}){
            my ($phase, $end_phase);
            #Calculate phase
            #5' UTR
            if (($exon->strand == 1 and $exon->start < $transcript->coding_region_start) or
                ($exon->strand == -1 and $exon->end > $transcript->coding_region_end)){
                $phase = -1;
            }
            #3' UTR
            elsif (($exon->strand == 1 and $exon->start > $transcript->coding_region_end) or
                   ($exon->strand == -1 and $exon->end < $transcript->coding_region_start)){
                $phase = -1;
            }
            #CDS
            else{
                $phase = ($exon->cdna_start($transcript) - $transcript->cdna_coding_start) % 3;
            }
            #Calculate end phase
            #5' UTR
            if (($exon->strand == 1 and $exon->end < $transcript->coding_region_start) or
                ($exon->strand == -1 and $exon->start > $transcript->coding_region_end)){
                $end_phase = -1;
            }
            #3' UTR
            elsif (($exon->strand == 1 and $exon->end > $transcript->coding_region_end) or
                   ($exon->strand == -1 and $exon->start < $transcript->coding_region_start)){
                $end_phase = -1;
            }
            #CDS
            else{
                if ($phase == -1){
                    $end_phase = ($exon->cdna_end($transcript) - $transcript->cdna_coding_start + 1) % 3;
                }
                else{
                    $end_phase = ($phase + $exon->length) % 3;
                }
            }
            if ($phase != $exon->phase){
                $message .= "Expected phase ($phase) does not match database phase (".$exon->phase.") for exon ".$exon->stable_id." ".$exon->start."-".$exon->end." on strand ".$exon->strand."\n";
            }
            if ($end_phase != $exon->end_phase){
                $message .= "Expected end phase ($end_phase) does not match database end phase (".$exon->end_phase.") for exon ".$exon->stable_id." ".$exon->start."-".$exon->end." on strand ".$exon->strand."\n";
            }
        }
    }
    return $message;
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
 Function  : Returns true if transcripts can be merged, i.e. they have exon-exon overlap but no exon-intron overlap
 Returntype: Boolean

=cut

sub can_be_merged {
    my ($tr, $db_tr) = @_;
    my $ex_ovlp = 0;
    EX:foreach my $exon (@{$tr->get_all_Exons}){
        #Exon-exon overlap: required
        foreach my $db_exon (@{$db_tr->get_all_Exons}){
            if ($exon->seq_region_start <= $db_exon->seq_region_end and $exon->seq_region_end >= $db_exon->seq_region_start){
                $ex_ovlp = 1;
            }
        }
        #Exon-intron overlap: not allowed
        foreach my $db_intron (@{$db_tr->get_all_Introns}){
            if ($exon->seq_region_start <= $db_intron->seq_region_end and $exon->seq_region_end >= $db_intron->seq_region_start){
              return 0;
            }
        }
    }
    #Intron-exon overlap: not allowed
    INT:foreach my $intron (@{$tr->get_all_Introns}){
        foreach my $db_exon (@{$db_tr->get_all_Exons}){
            if ($intron->seq_region_start <= $db_exon->seq_region_end and $intron->seq_region_end >= $db_exon->seq_region_start){
                return 0;
            }
        }
    }
    if ($ex_ovlp == 1){
      return 1;
    }
    return 0;
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
            for (my $i = $exon->seq_region_start; $i <= $exon->seq_region_end; $i++){
                $tr_loc{$i} = 1;
            }
        }
        my $max = 0;
        foreach my $db_tr (@best_candidates){
            my %seen; #Store genomic positions where novel transcript's exons overlaps database transcript's exons
            foreach my $db_exon (@{$db_tr->get_all_Exons}){
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
        $db_tr->add_Attributes($att);
    }
    
    #Change authorship
    $db_tr->transcript_author($tr->transcript_author);

    return $db_tr;
}



#Get transcript biotype from existing transcript biotypes
#Mainly for non-coding genes
sub get_transcript_biotype {
    my $db_gene = shift;
    my %biotypes;
    if ($db_gene->biotype =~ /^(processed_transcript|tec|comp_pipe)$/){
        foreach my $db_tr (@{$db_gene->get_all_Transcripts}){
            if ($db_tr->biotype =~ /^(antisense|lincrna|sense_intronic|sense_overlapping|bidirectional_promoter_lncrna|3\'_overlapping_ncrna)$/){
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



1;

