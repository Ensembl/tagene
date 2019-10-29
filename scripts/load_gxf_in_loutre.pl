#!/usr/bin/perl -w

#This script inserts new transcripts into loutre
#Usage: perl load_gxf_in_loutre.pl -file my_gene_data.gtf -dataset human_test [-write]


use strict;
use warnings;
use Getopt::Long;
use LoutreWrite::Config;
use LoutreWrite::Default;
use LoutreWrite::IntronFilter;
use Bio::Otter::Lace::Defaults;
use Bio::Otter::Server::Config;
$| = 1;

our %HOST_CDS_SET;
our %HOST_START_CODON_SET;
my $file;
my $source_info;
my $write;
my $dataset_name;
my $author_name;
my $remark;
my $use_comp_pipe_biotype;
my $no_artifact_check;
my $analysis_name;
my $tsource;
my $no_NFV;
my $no_comp_pipe;
my $no_intron_check;
my $host_biotype;
my $max_overlapped_loci;
my $filter_introns;
my $only_chr;

&GetOptions(
            'file=s'            => \$file,
            'dataset=s'         => \$dataset_name,
            'author=s'          => \$author_name,
            'source=s'          => \$source_info, #file with otter columns and read names supporting each transcript model
            'remark=s'          => \$remark, #remark to be added to every transcript model
            'comp_pipe!'        => \$use_comp_pipe_biotype,
            'analysis=s'        => \$analysis_name,
            'tsource=s'         => \$tsource,
            'no_check!'         => \$no_artifact_check,
            'no_NFV!'           => \$no_NFV,
            'no_comp_pipe!'     => \$no_comp_pipe,
            'no_intron_check!'  => \$no_intron_check,
            'host_biotype=s'    => \$host_biotype,
            'max_ov_loc=s'      => \$max_overlapped_loci,
            'filter_introns!'   => \$filter_introns,
            'chr=s'             => \$only_chr,
            'write!'            => \$WRITE,
            );

#die unless $dataset_name && $dataset_name =~ /test/;

#Find chromosome name if it has been passed as a job array index
$only_chr ||= $ENV{LSB_JOBINDEX};
if ($only_chr){
  $only_chr = "X" if $only_chr eq "23";
  $only_chr = "Y" if $only_chr eq "24";
}


my $usage =  <<_USAGE_;

This script stores the gene annotation from a gtf/gff3 file into an Otter (loutre) database.

perl load_gxf_in_loutre.pl -file ANNOTATION_FILE -source SOURCE_INFO_FILE -dataset DATASET_NAME -author AUTHOR_NAME -remark REMARK_STRING [-comp_pipe] [-write]
 -file           annotation file in GTF or GFF3 format - it can only contain exon lines - it will ignore anything but gene, transcript and exon lines
 -dataset        dataset name (species) already known by the Otter system, ie. present in the species.dat file on the live server
 -source         tsv file with transcript id and transcript sources, ie. RNA-seq read names supporting the model or tissues where it was found
 -author         author name in the corresponding loutre database
 -remark         annotation remark to be added to all transcripts, usually the ENA/GEO accession for the experiment that generated the models
 -comp_pipe      override the given biotypes and use "comp_pipe"
 -analysis       analysis logic name (default is "Otter")
 -tsource        transcript source (default is "havana")
 -no_check       do no check for transcripts spanning a large number of genes
 -no_NFV         do not add the 'not for VEGA' attributes
 -no_comp_pipe   do not add the comp_pipe attributes by default ('Assembled from ... reads', 'ID: ...' remarks)
 -no_intron_check  allow transcripts with intron chains fully or partially identical to others in the database
 -host_biotype   restrict host genes by biotype (comma-separated list)
 -max_ov_loc     maximum number of existing loci that a novel transcript can overlap at the exon level (ignore the transcript if exceeded)
 -filter_introns assess introns and ignore transcript if at least an intron does not pass the filters
 -chr            restrict to annotation on this chromosome
 -write          store the gene annotation in the loutre database


_USAGE_

die $usage unless ($file and $dataset_name and $author_name);

my $time0 = time();

#Connect to loutre database
#DataSet interacts directly with an otter database
my $dataset = Bio::Otter::Server::Config->SpeciesDat->dataset($dataset_name);
my $otter_dba = $dataset->otter_dba or die "can't get db adaptor\n";
$otter_dba->dbc->reconnect_when_lost(1);
my $sa = $otter_dba->get_SliceAdaptor();
my $ga = $otter_dba->get_GeneAdaptor();
my $ta = $otter_dba->get_TranscriptAdaptor();
$DBA{'otter'} = $otter_dba;


#Read data file
my $genes = LoutreWrite::Default->parse_gxf_file($file);

#Make gene objects
my $gene_objects = LoutreWrite::Default->make_vega_objects($genes, $otter_dba, $author_name, $remark, $use_comp_pipe_biotype, $analysis_name, $tsource, $no_comp_pipe, $no_NFV, $only_chr);


#Long artifact transcripts, spanning multiple real loci, make long artificial genes
unless ($no_artifact_check){
  #Generate a kill list of "artifact" transcripts (those with long introns spanning multiple genes, likely caused by misalignments)
  my $kill_list = LoutreWrite::Default->check_artifact_transcripts($gene_objects);

  #Remove artifacts and split genes left with genomic gaps
  $gene_objects = LoutreWrite::Default->recluster_transcripts($gene_objects, $kill_list);
}

#Remove transcripts overlapping more than the allowed number of existing loci at the exon level
my $gene_objects_2;
if ($max_overlapped_loci){
  my $kill_list_2 = LoutreWrite::Default->check_overlapped_loci($gene_objects, $max_overlapped_loci);
  $gene_objects_2 = LoutreWrite::Default->recluster_transcripts($gene_objects, $kill_list_2);
  foreach my $tid (keys %{$kill_list_2}){
    print "TR2: $tid: max allowed number of overlapped loci exceeded\n";
  }
}

#Add source info
#my $gene_objects_3;
if ($source_info){
    $gene_objects_2 = LoutreWrite::Default->add_sources($gene_objects_2, $source_info);
}

my $count = 1;
my $total = scalar(@$gene_objects_2);
#Process gene objects
foreach my $gene_obj (@$gene_objects_2){
    print "##########\n\nLOOKING AT GENE IN ".$gene_obj->seq_region_name.":".$gene_obj->start."-".$gene_obj->end."  ($count/$total)\n\n";
    $count++;

    #Check whether add as novel gene or merge with existing gene
    #This method returns the most suitable gene to be merged with, or undef if there is none (ie. novel gene)
    #my $host_gene = LoutreWrite::Default->find_host_gene($gene_obj);
    #if ($host_gene){
    #    print "HOST: ".$host_gene->stable_id."\n\n";
    #}

   #IT HAS TO BE DONE BY SPLITTING GENES    
    #Readthrough genes -> should be split and assign host genes independently by transcript
    #Why not replace "find_host_gene" with "assign_host_gene" which would return an array of gene objects, each one with their new stable id if they have a host gene?
    my $new_gene_objects = LoutreWrite::Default->assign_host_gene($gene_obj, 1);
    my %allowed_biotypes;
    if ($host_biotype){
      %allowed_biotypes = map {$_=>1} split(/,/, $host_biotype);
    }
    GENE:foreach my $new_gene_obj (@$new_gene_objects){
        print "HOST: ".($new_gene_obj->stable_id || "NONE")."\n";

        #Find biotype of host gene
        if ($new_gene_obj->stable_id =~ /^OTT/){
            my $wrong_host;
            my $host_gene = $ga->fetch_by_stable_id($new_gene_obj->stable_id);
            if ($host_biotype and !($allowed_biotypes{$host_gene->biotype})){
                print "Gene ".$host_gene->stable_id." will be ignored as its biotype (".$host_gene->biotype.") is not allowed\n";
                $wrong_host = 1;
            }
            #Ignore if gene has an ASB_protein_coding remark
            elsif ($host_biotype and !($host_biotype =~ /protein_coding/) and scalar(grep {$_->value eq "ASB_protein_coding"} @{$host_gene->get_all_Attributes('remark')})){
                print "Gene ".$host_gene->stable_id." will be ignored as it has an ASB_protein_coding remark\n";
                $wrong_host = 1;
            }
            #Double check that the gene does not have a coding transcript
            #There are a few cases in loutre_human
            else{
                foreach my $tr (@{$host_gene->get_all_Transcripts}){
                    if ($tr->translate and $host_biotype and !($host_biotype =~ /protein_coding/)){
                        print "Gene ".$host_gene->stable_id." will be ignored as it has a coding transcript\n";
                        $wrong_host = 1;
                        last;
                    }
                }
            }
            if ($wrong_host){
                foreach my $transcript (@{$new_gene_obj->get_all_Transcripts}){
                    my $t_name = $transcript->stable_id || $transcript->get_all_Attributes('hidden_remark')->[0]->value;
                    print "TR2: $t_name: host gene biotype not allowed\n";
                }
                next GENE;
            }
            
            #$new_gene_obj = LoutreWrite::Default->assign_cds_to_transcripts($new_gene_obj, $host_gene);
        }
        
        #Predict intron outcome
        #Only if novel intron?
      if ($filter_introns){
        my %PASSED_INTRONS;
        TR:foreach my $transcript (@{$new_gene_obj->get_all_Transcripts}){
            my $t_name = $transcript->stable_id || $transcript->get_all_Attributes('hidden_remark')->[0]->value;
            foreach my $intron (@{$transcript->get_all_Introns}){
                my $pass;
                if ($PASSED_INTRONS{$intron->seq_region_start."-".$intron->seq_region_end}){
                  $pass = 1;
print "Using PASSED_INTRONS for intron at ".$intron->seq_region_start."-".$intron->seq_region_end."\n";
                }
                else {
                  $pass = LoutreWrite::IntronFilter->predict_outcome($intron, $transcript, 1);
                  $PASSED_INTRONS{$intron->seq_region_start."-".$intron->seq_region_end} = $pass;
print "Predicting outcome for intron at ".$intron->seq_region_start."-".$intron->seq_region_end."\n";
                }
                #Give another chance to HiSeq-supported introns
                if (($file =~ /_HiSeq/ or $remark =~ /_HiSeq/) and $pass == 0){
                  $pass = LoutreWrite::IntronFilter->exonerate_support($intron, $transcript);
print "Testing exonerate support for intron at ".$intron->seq_region_start."-".$intron->seq_region_end."\n";
                }
                print "PASS: ".($pass ? "YES" : "NO")."\n";
                if (!$pass){
                    print "TR2: $t_name: intron did not pass the filter\n";
                    if (scalar(@{$new_gene_obj->get_all_Transcripts}) > 1){
                        #$new_gene_obj->remove_Transcript($transcript);
                        #Ensembl API remove_Transcript method requires a dbID, but transcript objects do not have one at this stage 
                        #Reimplement method using transcript name as identifier
                        my $array = $new_gene_obj->{_transcript_array};
                        @$array = grep { $_->get_all_Attributes('hidden_remark')->[0]->value ne $t_name } @$array;
                        $ga->update_coords($new_gene_obj);
                        next TR;
                    }
                    else{
                        next GENE;
                    }
                }
            }
        }
        if (scalar(@{$new_gene_obj->get_all_Transcripts}) == 0){
            next GENE;
        }
      }
    
        #Store gene
        #if ($write){
            my $mode;
            #if ($host_gene){
            #    $gene_obj->stable_id($host_gene->stable_id);
            #    $mode = "update";
            #}
            #else{
            #    $mode = "add";
            #}
            if ($new_gene_obj->stable_id =~ /^OTT/){
                $mode = "update";
            }
            else{
                $mode = "add";
            }        
            my $novel_gene = $new_gene_obj->stable_id ? 0 : 1;
            print "\nMODE: $mode\n";
            my ($msg, $log) = LoutreWrite::Default->process_gene_2($new_gene_obj, $mode, "aggr", $dataset_name, $otter_dba, $no_intron_check);
        #print "HOST_CDS_SET_B=".scalar(keys %HOST_CDS_SET)."\n";
        #print "HOST_START_CODON_SET_B=".scalar(keys %HOST_START_CODON_SET)."\n";
            if ($msg =~ /lock ok,write ok,unlock ok(,(.+))?/){
                if ($2){
                    if ($novel_gene){
                        print "\nRESULT: novel ".$2." was created\n";
                        print join("\n", @$log)."\n";
                    }
                    else{
                        print "\nRESULT: ".$2." was modified\n";
                        print join("\n", @$log)."\n";
                    }
                }
                else{
                    print "\nRESULT: gene ".$new_gene_obj->stable_id." checked but not modified\n";
                    print join("\n", @$log)."\n";
                }
            }
            else{
                if ($WRITE){
                    print "\nRESULT: gene ".$new_gene_obj->stable_id." could not be unlocked\n";
                }
                else{
                    print "\nRESULT: gene ".$new_gene_obj->stable_id." not modified (WRITE = 0)\n";
                }
            }
        #}
    }
}


print "\nDONE - Elapsed time: ".format_time(time()-$time0)."\n\n";



sub format_time {
  my $time = shift; #Time in seconds
  my $string = "";
  $string .= int($time/3600)."h ";
  $time = $time % 3600;
  $string .= int($time/60)."m ".($time % 60)."s";
  return $string;
}



