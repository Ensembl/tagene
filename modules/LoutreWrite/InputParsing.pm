
package LoutreWrite::InputParsing;

use strict;
use warnings;
use Bio::Vega::Gene;
use Bio::Vega::Transcript;
use Bio::Vega::Exon;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Attribute;
use LoutreWrite::AnnotUpdate;

our $CP_BIOTYPE = "comp_pipe";


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
            next unless $att =~ /\w/;
            my ($name, $value);
            if ($file =~ /\.gff3(.gz)?$/){
                ($name, $value) = split(/=/, $att);
            }
            elsif ($file =~ /\.g(t|f)f(.gz)?$/){
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
                                'gene_type'   => $attribs{'gene_type'} || "missing_biotype",
                                'gene_name'   => $attribs{'gene_name'},
                                'gene_status' => $attribs{'gene_status'}
                                };
        }
        elsif ($type eq "transcript"){
            $data{$gene_id}{transcripts}{$transcript_id} = {
                                                            'transcript_type'   => $attribs{'transcript_type'} || "missing_biotype",
                                                            'transcript_name'   => $attribs{'transcript_name'},
                                                            'transcript_status' => $attribs{'transcript_status'}
                                                            };
        }
        elsif ($type eq "exon"){
            my $exon_id = $attribs{'exon_id'} || $attribs{'exon_number'} || $.; #gtf line number as default exon id
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
            my $exon_id = $attribs{'exon_id'} || $attribs{'exon_number'} || $.; #gtf line number as default exon id
            $data{$gene_id}{transcripts}{$transcript_id}{cds}{$exon_id} = { 
                                                                                'chr'    => $chr,
                                                                                'start'  => $start,
                                                                                'end'    => $end,
                                                                                'strand' => $strand,
                                                                                'phase'  => "p".$phase,
                                                                           }; print "B: ".$exon_id." ".$start."\n";
        }
        elsif ($type eq "stop_codon"){ #Relevant for GENCODE gtf type files, where stop codon is not included in CDS
            my $exon_id = $attribs{'exon_id'} || $attribs{'exon_number'} || $.; #gtf line number as default exon id
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
 Arg[5]    : force 'comp_pipe' biotype
 Arg[6]    : gene and transcript analysis logic name
 Arg[7]    : gene and transcript source
 Arg[8]    : Boolean - if true, the "not for VEGA" remark will not be added
 Arg[9]    : String - chromosome (ignore genes not on this chromosome)
 Function  : convert gene annotation into Vega gene, transcript and exon objects
 Returntype: list of Bio::Vega::Gene objects

=cut

sub make_vega_objects {
    my ($self, $genes, $dba, $author_name, $remark, $force_cp_biotype, $analysis_name, $source, $no_NFV, $only_chr, $assembly_version) = @_;
    
    #Some common features
    $source ||= "havana";
    $analysis_name ||= "Otter";
    my $analysis = Bio::EnsEMBL::Analysis->new( -logic_name => $analysis_name );
    
    my $nfv_remark = Bio::EnsEMBL::Attribute->new(
                                -code => 'remark',
                                -value => 'not for VEGA'
                     );

    my $tagene_gene_attrib = Bio::EnsEMBL::Attribute->new(
                                -code => 'TAGENE_gene',
                                -value => '1'
                             );
    my $tagene_gene_remark = Bio::EnsEMBL::Attribute->new(
                                -code => 'remark',
                                -value => 'TAGENE_gene'
                        );                             
    my $tagene_attrib = Bio::EnsEMBL::Attribute->new(
                                -code => 'TAGENE_transcript',
                                -value => '1'
                        );
    my $tagene_remark = Bio::EnsEMBL::Attribute->new(
                                -code => 'remark',
                                -value => 'TAGENE_transcript'
                        );

    my $source_remark;
    if ($remark){
        $source_remark = Bio::EnsEMBL::Attribute->new(
                                -code => 'remark',
                                -value => $remark
                         );
    }            
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
                                        -biotype    => ($force_cp_biotype ? $CP_BIOTYPE : ($genes{$gid}{'gene_type'} || "missing_biotype")),
                                        -analysis   => $analysis,
                                        -source     => $source,
                                        -status     => $genes{$gid}{'gene_status'} || "UNKNOWN"
                                        );
        #$gene->add_Attributes(new Bio::EnsEMBL::Attribute(-code => 'name', -value => $genes{$gid}{'gene_name'})); #Let the script generate a name if needed
        $gene->gene_author($author);
        if ($source_remark){
            #$gene->add_Attributes($source_remark);
        }
        unless ($no_NFV){
            $gene->add_Attributes($nfv_remark);
        }
        $gene->add_Attributes($tagene_gene_remark);
        #$gene->add_Attributes(Bio::EnsEMBL::Attribute->new(-code => 'hidden_remark', -value => $genes{$gid}{'gene_name'}));
        
        #Add biotype-status combination sanity check!!!

        #Make transcript objects
        foreach my $tid (keys %{$genes{$gid}{transcripts}}){
            my $transcript = Bio::Vega::Transcript->new(
                                                        -stable_id  => ($tid =~ /^OTT/ ? $tid : ""),
                                                        -biotype    => ($force_cp_biotype ? $CP_BIOTYPE : ($genes{$gid}{transcripts}{$tid}{'transcript_type'} || "missing_biotype")),
                                                        -analysis   => $analysis,
                                                        -source     => $source,
                                                        -status     => $genes{$gid}{transcripts}{$tid}{'transcript_status'} || "UNKNOWN"
                                                        );
            #$transcript->add_Attributes(new Bio::EnsEMBL::Attribute(-code => 'name', -value => $genes{$gid}{transcripts}{$tid}{'transcript_name'})); #Let the script generate a name if needed
            $transcript->transcript_author($author);
            $transcript->add_Attributes($tagene_remark);
            if ($source_remark){
              $transcript->add_Attributes($source_remark);
            }
            unless ($no_NFV){
                $transcript->add_Attributes($nfv_remark);
            }
            $transcript->add_Attributes(Bio::EnsEMBL::Attribute->new(-code => 'hidden_remark', 
                                                                     -value => "ID: ".$genes{$gid}{transcripts}{$tid}{'transcript_name'}));
        
            #Make exon objects
            foreach my $exid (sort keys %{$genes{$gid}{transcripts}{$tid}{exons}}){
                my $exon = Bio::Vega::Exon->new(
                                                -start   => $genes{$gid}{transcripts}{$tid}{exons}{$exid}{'start'},
                                                -end     => $genes{$gid}{transcripts}{$tid}{exons}{$exid}{'end'},
                                                -strand  => ($genes{$gid}{transcripts}{$tid}{exons}{$exid}{'strand'} eq "+" ? 1 : -1),
                                                );
                my $chr = $genes{$gid}{transcripts}{$tid}{exons}{$exid}{'chr'};
                $chr =~ s/^chr//;
                next GENE if ($only_chr and $chr ne $only_chr);
                if (!($slices{$chr})){
                    my $slice = $sa->fetch_by_region("toplevel", $chr, undef, undef, undef, $assembly_version);
                    unless ($slice){
                      die "Can't find slice for $chr\n";
                    }
                    $slices{$chr} = $slice;
                }
                $exon->slice($slices{$chr});
                $exon->phase(-1);
                $exon->end_phase(-1);
#foreach my $exonid (sort keys %{$genes{$gid}{transcripts}{$tid}{cds}}){
#  print "B$exid: ".$exonid."  ".$genes{$gid}{transcripts}{$tid}{cds}{$exonid}{'start'}."\n";
#} 
foreach my $exonid (sort keys %{$genes{$gid}{transcripts}{$tid}{stop_codon}}){
    print "d$exid: ".$exonid."  ".$genes{$gid}{transcripts}{$tid}{stop_codon}{$exonid}{'start'}."\n";
}                
                #Assign exon phases from gxf file if exon starts within CDS
                #Compare also with stop_codon coordinates as stop codon is not included in CDS in GENCODE gtf files
                my $gtf_phase; 
                if ($genes{$gid}{transcripts}{$tid}{cds}{$exid}){
                    $gtf_phase = $genes{$gid}{transcripts}{$tid}{cds}{$exid}{'phase'};
                }
                elsif ($genes{$gid}{transcripts}{$tid}{stop_codon}{$exid}){
                    $gtf_phase = $genes{$gid}{transcripts}{$tid}{stop_codon}{$exid}{'phase'};
                }
                
# foreach my $exonid (sort keys %{$genes{$gid}{transcripts}{$tid}{cds}}){
#  print "C$exid: ".$exonid."  ".$genes{$gid}{transcripts}{$tid}{cds}{$exonid}{'start'}."\n";
#} 
                 
                print "\n$gtf_phase \n" if $gtf_phase;
                if ($gtf_phase and ($gtf_phase eq "p0" or $gtf_phase eq "p1" or $gtf_phase eq "p2")){
                #if ($gtf_phase){
                    print "\nA\n";
foreach my $exonid (sort keys %{$genes{$gid}{transcripts}{$tid}{stop_codon}}){
    print "e$exid: ".$exonid."  ".$genes{$gid}{transcripts}{$tid}{stop_codon}{$exonid}{'start'}."\n";
}  
                    if ((($exon->start == $genes{$gid}{transcripts}{$tid}{cds}{$exid}{'start'} or 
                          defined($genes{$gid}{transcripts}{$tid}{stop_codon}{$exid}) and 
                          $exon->start == $genes{$gid}{transcripts}{$tid}{stop_codon}{$exid}{'start'})                        
                          and $exon->strand == 1) or 
                        (($exon->end == $genes{$gid}{transcripts}{$tid}{cds}{$exid}{'end'} or 
                          defined($genes{$gid}{transcripts}{$tid}{stop_codon}{$exid}) and 
                          $exon->end == $genes{$gid}{transcripts}{$tid}{stop_codon}{$exid}{'end'})                     
                          and $exon->strand == -1)
                    ){
                        $exon->phase(0) if $gtf_phase eq "p0";
                        $exon->phase(1) if $gtf_phase eq "p2";
                        $exon->phase(2) if $gtf_phase eq "p1"; print "\nB".$exon->phase."\n";
foreach my $exonid (sort keys %{$genes{$gid}{transcripts}{$tid}{stop_codon}}){
    print "f$exid: ".$exonid."  ".$genes{$gid}{transcripts}{$tid}{stop_codon}{$exonid}{'start'}."\n";
} 
                  
                        if ((($exon->end == $genes{$gid}{transcripts}{$tid}{cds}{$exid}{'end'} or 
                              $exon->end == $genes{$gid}{transcripts}{$tid}{stop_codon}{$exid}{'end'}) 
                              and $exon->strand == 1) or 
                            (($exon->start == $genes{$gid}{transcripts}{$tid}{cds}{$exid}{'start'} or 
                              $exon->start == $genes{$gid}{transcripts}{$tid}{stop_codon}{$exid}{'start'}) 
                              and $exon->strand == -1)
                         ){
                           $exon->end_phase(($exon->phase + $exon->length) % 3); print "\nC".$exon->end_phase."\n";
foreach my $exonid (sort keys %{$genes{$gid}{transcripts}{$tid}{stop_codon}}){
    print "g$exid: ".$exonid."  ".$genes{$gid}{transcripts}{$tid}{stop_codon}{$exonid}{'start'}."\n";
}                        
                          }
                    }
                    else{
foreach my $exonid (sort keys %{$genes{$gid}{transcripts}{$tid}{stop_codon}}){
    print "h$exid: ".$exonid."  ".$genes{$gid}{transcripts}{$tid}{stop_codon}{$exonid}{'start'}."\n";
}                    
                        #Half-coding start codon: no phase but end_phase
                        if ($exon->end == $genes{$gid}{transcripts}{$tid}{cds}{$exid}{'end'} and $exon->strand == 1){
                            $exon->end_phase(($exon->end - $genes{$gid}{transcripts}{$tid}{cds}{$exid}{'start'} + 1) % 3); print "\nD".$exon->end_phase."\n";
                        }
                        elsif ($exon->start == $genes{$gid}{transcripts}{$tid}{cds}{$exid}{'start'} and $exon->strand == -1){
                            $exon->end_phase(($genes{$gid}{transcripts}{$tid}{cds}{$exid}{'end'} - $exon->start + 1) % 3); print "\nE".$exon->end_phase."\n";
                        }
foreach my $exonid (sort keys %{$genes{$gid}{transcripts}{$tid}{stop_codon}}){
    print "i$exid: ".$exonid."  ".$genes{$gid}{transcripts}{$tid}{stop_codon}{$exonid}{'start'}."\n";
}
                    } 
                }
foreach my $exonid (sort keys %{$genes{$gid}{transcripts}{$tid}{stop_codon}}){
    print "j$exid: ".$exonid."  ".$genes{$gid}{transcripts}{$tid}{stop_codon}{$exonid}{'start'}."\n";
} 
                $transcript->add_Exon($exon);
            }
           
            #Make translation object
            if (scalar(keys %{$genes{$gid}{transcripts}{$tid}{cds}})){
                my $cds_start = $transcript->end;
                my $cds_end = $transcript->start;
                foreach my $exid (sort keys %{$genes{$gid}{transcripts}{$tid}{cds}}){
 print "B: ".$exid."  ".$genes{$gid}{transcripts}{$tid}{cds}{$exid}{'start'}."\n";
                    if ($genes{$gid}{transcripts}{$tid}{cds}{$exid}{'start'} < $cds_start){
                        $cds_start = $genes{$gid}{transcripts}{$tid}{cds}{$exid}{'start'};
                    }
                    if ($genes{$gid}{transcripts}{$tid}{cds}{$exid}{'end'} > $cds_end){
                        $cds_end = $genes{$gid}{transcripts}{$tid}{cds}{$exid}{'end'};
                    }
                }
                print "CDS_START=$cds_start\nCDS_END=$cds_end\n";
#                foreach my $exonid (sort keys %{$genes{$gid}{transcripts}{$tid}{stop_codon}}){
#  print "D:".$exonid."  ".$genes{$gid}{transcripts}{$tid}{stop_codon}{$exonid}{'start'}."\n";
#}
                #Include stop codon coordinates in the CDS 
                if ($genes{$gid}{transcripts}{$tid}{stop_codon}){
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
                }
                if ($cds_start and $cds_end){
                print "CDS_START=$cds_start\nCDS_END=$cds_end\n";
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
               } 
               #Re-check exon phases
               my $check_message = $self->check_exon_phases($transcript);
               print $check_message."\n";
            }

            #Add transcript to gene
            $gene->add_Transcript($transcript);
        }
        
print "ZZZ".$gene->biotype."\n";        
print "ZZZ".join(", ", map {$_->biotype} @{$gene->get_all_Transcripts})."\n"; 
        #Update gene and transcript biotypes
        if (!($gene->biotype) or $gene->biotype eq "missing_biotype"){
            $gene = LoutreWrite::AnnotUpdate->assign_biotypes($gene); #This also updates description if appropriate 
        }
print "ZZZ".$gene->biotype."\n"; 
print "ZZZ".join(", ", map {$_->biotype} @{$gene->get_all_Transcripts})."\n";          
        #Update gene and transcript status - only lncRNA genes for now
        $gene = LoutreWrite::AnnotUpdate->assign_status($gene);
        
        push(@objects, $gene);
    }

    return \@objects;
}



=head2 check_exon_phases

 Arg[1]    : Bio::Vega::Transcript object
 Function  : 
 Returntype: String

=cut

sub check_exon_phases {
    my ($self, $transcript) = @_;
    my $message = " ";
    if ($transcript->translation and scalar(grep {$_->value == 1} @{$transcript->get_all_Attributes('CDS_start_NF')}) == 0){
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
                $message .= "Expected phase ($phase) does not match database phase (".$exon->phase.") for exon ".($exon->stable_id || "")." ".$exon->start."-".$exon->end." on strand ".$exon->strand."\n";
            }
            if ($end_phase != $exon->end_phase){
                $message .= "Expected end phase ($end_phase) does not match database end phase (".$exon->end_phase.") for exon ".($exon->stable_id || "")." ".$exon->start."-".$exon->end." on strand ".$exon->strand."\n";
            }
        }
    }
    return $message;
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





1;
