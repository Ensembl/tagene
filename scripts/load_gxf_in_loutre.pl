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

&GetOptions(
            'file=s'      => \$file,
            'dataset=s'   => \$dataset_name,
            'author=s'    => \$author_name,
            'source=s'    => \$source_info, #file with otter columns and read names supporting each transcript model
            'remark=s'    => \$remark, #remark to be added to every transcript model
            'comp_pipe!'  => \$use_comp_pipe_biotype,
            'write!'      => \$write,            
            );

die unless $dataset_name && $dataset_name =~ /test/;

my $usage =  <<_USAGE_;

This script stores the gene annotation from a gtf/gff3 file into an Otter (loutre) database.

perl load_gxf_in_loutre.pl -file ANNOTATION_FILE -source SOURCE_INFO_FILE -dataset DATASET_NAME -author AUTHOR_NAME -remark REMARK_STRING [-comp_pipe] [-write]
 -file      annotation file in GTF or GFF3 format - it can only contain exon lines - it will ignore anything but gene, transcript and exon lines
 -dataset   dataset name (species) already known by the Otter system, ie. present in the species.dat file on the live server
 -source    tsv file with transcript id and transcript sources, ie. RNA-seq read names supporting the model or tissues where it was found
 -author    author name in the corresponding loutre database
 -remark    annotation remark to be added to all transcripts, usually the ENA/GEO accession for the experiment that generated the models
 -comp_pipe override the given biotypes and use "comp_pipe"
 -write     store the gene annotation in the loutre database


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
my $gene_objects = LoutreWrite::Default->make_vega_objects($genes, $otter_dba, $author_name, $remark, $use_comp_pipe_biotype);

#Long artifact transcripts, spanning multiple real loci, make long artificial genes
#Generate a kill list of "artifact" transcripts (those with long introns spanning multiple genes, likely caused by misalignments)
my $kill_list = LoutreWrite::Default->check_artifact_transcripts($gene_objects);

#Remove artifacts and split genes left with genomic gaps
$gene_objects = LoutreWrite::Default->recluster_transcripts($gene_objects, $kill_list);

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
    my $gene_objects = LoutreWrite::Default->assign_host_gene($gene_obj);
    foreach my $new_gene_obj (@$gene_objects){
        print "HOST: ".$new_gene_obj->stable_id."\n";

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
            my $msg = LoutreWrite::Default->process_gene($new_gene_obj, $mode, $dataset_name, $otter_dba);
            if ($msg =~ /lock ok,write ok,unlock ok/){
                print "\nRESULT: Gene ".$new_gene_obj->stable_id." written\n";
            }
            else{
                print "\nRESULT: Gene ".$new_gene_obj->stable_id." failed\n";
            }
        }
    }
}




