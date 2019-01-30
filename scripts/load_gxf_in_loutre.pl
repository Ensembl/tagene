#!/usr/bin/perl -w

#This script inserts new transcripts into loutre
#Usage: perl load_gxf_in_loutre.pl -file my_gene_data.gtf -dataset human_test [-write]


use strict;
use warnings;
use Getopt::Long;
use LoutreWrite::Default;
use Bio::Otter::Lace::Defaults;
use Bio::Otter::Server::Config;
$| = 1;

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
my $no_comp_pipe;
my $no_intron_check;

&GetOptions(
            'file=s'        => \$file,
            'dataset=s'     => \$dataset_name,
            'author=s'      => \$author_name,
            'source=s'      => \$source_info, #file with otter columns and read names supporting each transcript model
            'remark=s'      => \$remark, #remark to be added to every transcript model
            'comp_pipe!'    => \$use_comp_pipe_biotype,
            'analysis=s'    => \$analysis_name,
            'tsource=s'     => \$tsource,
            'no_check!'     => \$no_artifact_check,
            'no_comp_pipe!' => \$no_comp_pipe,
            'no_intron_check!' => \$no_intron_check,
            'write!'        => \$write,
            );

#die unless $dataset_name && $dataset_name =~ /test/;

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
 -no_comp_pipe   do not add the comp_pipe attributes by default ('not_for_VEGA', 'Assembled from ... reads', 'ID: ...' remarks)
 -no_intron_check  allow transcripts with intron chains fully or partially identical to others in the database
 -write          store the gene annotation in the loutre database


_USAGE_

die $usage unless ($file and $dataset_name and $author_name);



#Connect to loutre database
#DataSet interacts directly with an otter database
my $dataset = Bio::Otter::Server::Config->SpeciesDat->dataset($dataset_name);
my $otter_dba = $dataset->otter_dba or die "can't get db adaptor\n";
my $sa = $otter_dba->get_SliceAdaptor();
my $ga = $otter_dba->get_GeneAdaptor();
my $ta = $otter_dba->get_TranscriptAdaptor();


#Read data file
my $genes = LoutreWrite::Default->parse_gxf_file($file);

#Make gene objects
my $gene_objects = LoutreWrite::Default->make_vega_objects($genes, $otter_dba, $author_name, $remark, $use_comp_pipe_biotype, $analysis_name, $tsource, $no_comp_pipe);

#Long artifact transcripts, spanning multiple real loci, make long artificial genes
unless ($no_artifact_check){
  #Generate a kill list of "artifact" transcripts (those with long introns spanning multiple genes, likely caused by misalignments)
  my $kill_list = LoutreWrite::Default->check_artifact_transcripts($gene_objects);

  #Remove artifacts and split genes left with genomic gaps
  $gene_objects = LoutreWrite::Default->recluster_transcripts($gene_objects, $kill_list);
}

#Add source info
if ($source_info){
    $gene_objects = LoutreWrite::Default->add_sources($gene_objects, $source_info);
}

#Process gene objects
foreach my $gene_obj (@$gene_objects){
    print "##########\n\nLOOKING AT GENE IN ".$gene_obj->seq_region_name.":".$gene_obj->start."-".$gene_obj->end."\n\n";

    #Check whether add as novel gene or merge with existing gene
    #This method returns the most suitable gene to be merged with, or undef if there is none (ie. novel gene)
    #my $host_gene = LoutreWrite::Default->find_host_gene($gene_obj);
    #if ($host_gene){
    #    print "HOST: ".$host_gene->stable_id."\n\n";
    #}

   #IT HAS TO BE DONE BY SPLITTING GENES    
    #Readthrough genes -> should be split and assign host genes independently by transcript
    #Why not replace "find_host_gene" with "assign_host_gene" which would return an array of gene objects, each one with their new stable id if they have a host gene?
    my $gene_objects = LoutreWrite::Default->assign_host_gene($gene_obj, 1);
    foreach my $new_gene_obj (@$gene_objects){
        print "HOST: ".$new_gene_obj->stable_id."\n";
        
        #Find biotype of host gene
        if ($new_gene_obj->stable_id =~ /^OTT/){
          my $host_gene = $ga->fetch_by_stable_id($new_gene_obj->stable_id);
          if ($host_gene->biotype eq "protein_coding"){
            #$new_gene_obj = LoutreWrite::Default->assign_cds_to_transcripts($new_gene_obj, $host_gene);
          }
        }

        #Store gene
        if ($write){
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
        
            print "\nMODE: $mode\n";
            my $msg = LoutreWrite::Default->process_gene($new_gene_obj, $mode, $dataset_name, $otter_dba, $no_intron_check);
            if ($msg =~ /lock ok,write ok,unlock ok(,(.+))?/){
                if ($2){
                    print "\nRESULT: ".$2." was modified\n";
                }
                else{
                    print "\nRESULT: gene ".$new_gene_obj->stable_id." checked but not modified\n";
                }
            }
            else{
                print "\nRESULT: gene ".$new_gene_obj->stable_id." could not be unlocked\n";
            }
        }
    }
}




