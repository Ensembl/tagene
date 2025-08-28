# README

This software allows the user to import gene annotation from a gtf or gff3 file into a MySQL database with the Ensembl/Otter schema.
It is part of the gene annotation workflow followed by the Ensembl Havana team.

The TAGENE pipeline:

* reads gene annotation from an input file in gtf or gff3 format, normally derived from long RNA-seq read data;
* filters the input annotation to remove transcripts containing potentially spurious splice sites;
* compares the input annotation with the reference annotation from a database with the Ensembl/Otter schema;
* stores in the database the input annotation that provides some novelty, either by adding novel transcripts or by extending existing transcripts using novel exonic sequence.


## Software
The pipeline code is mostly written in Perl. It relies heavily on the Ensembl and Otter Perl APIs.
The target database is similar to an Ensembl core database, with some extensions to track annotation authorship and avoid concurrent annotation of the same region. 
Database access is granted via an Ensembl private repository (server-config).
The pipeline manager is eHive (https://ensembl-hive.readthedocs.io).


### Dependencies

* ensembl (https://github.com/Ensembl/ensembl)
* ensembl-otter (https://github.com/Ensembl/ensembl-otter)
* server-config (https://github.com/Ensembl/server-config) [private]
* ensembl-analysis (https://github.com/Ensembl/ensembl-analysis)
* ensembl-hive (https://github.com/Ensembl/ensembl-hive)



### Installation

Install software dependencies:
```sh
cd $HOME/software
git clone https://github.com/Ensembl/tagene.git
git clone https://github.com/Ensembl/ensembl.git
git clone https://github.com/Ensembl/ensembl-otter.git
git clone https://github.com/Ensembl/ensembl-analysis.git
git clone https://github.com/Ensembl/ensembl-hive.git
```
 
Install cpanm to have other Perl dependencies installed:
```sh
curl -L http://cpanmin.us | perl - App::cpanminus
#Add "$HOME/perl5/bin" to the $PATH environment variable
export PATH="$HOME/perl5/bin:$PATH"
```

Install other known Perl dependencies (this list may not be comprehensive):
```sh
cpanm Try::Tiny
cpanm Readonly
cpanm Lingua::EN::Inflect
cpanm Test::MockObject
cpanm Data::Rmap
cpanm Hash::Merge::Simple
cpanm Crypt::JWT
cpanm Config::IniFiles
cpanm Text::sprintfn
cpanm DBD::mysql
cpanm DBI
```

Set your PERL5LIB variable so that it includes the paths to the installed dependencies:
```sh
export PERL5LIB="$PERL5LIB:$HOME/perl5/lib/perl5:$HOME/software/tagene/modules:$HOME/software/ensembl/modules:$HOME/software/ensembl-otter/modules:$HOME/software/ensembl-analysis/modules:$HOME/software/ensembl-hive/modules"
```

Your PATH variable must include the path to ensembl-hive/scripts:
```sh
export PATH="$HOME/software/ensembl-hive/scripts:$PATH"
```

Set up a basic Otter client configuration so that the pipeline can pick up dataset names from the server-config/species.dat file for database access:
```sh
#Create an .otter directory and clone the server-config repository there
mkdir $HOME/.otter && cd $HOME/.otter
git clone https://github.com/Ensembl/server-config

#Create environment variable
ANACODE_SERVER_CONFIG=${HOME}/.otter/server-config
#Restrict read permission to group
chmod 0750 $HOME/.otter/server-config
```

### Test
[TO DO]


## Pipeline steps
### 1. Input file pre-filter (splice-site-based filter of input annotation)

The input annotation file can be pre-filtered at the start of the pipeline to remove transcripts with unreliable splice sites. This step consists of a set of filters that check for:

- introns having a Recount3 score below certain threshold (eg. 50) unless already present in the reference annotation
- introns with any obvious misalignments around the splice sites (requires access to BAM file with supporting reads)
- splice site overlap with repeat features
- overlap with processed pseudogenes
- overlap with coding genes in the opposite strand

The final output is a list of input transcripts that contain at least an intron that failed any of the filters. Transcripts in the list will be ignored in subsequent steps.

This step is optional and can be skipped if:

- the filter list from a previous pipeline run using the same input annotation file is provided
- the filters are carried out during the annotation step
- no intron filters are required

### 2. Reset the target database

For pipeline tests, the target database can be synchronised with another reference database before the annotation step. Basically, this step consists of replacing the target database with a source database, which will normally be a database with the current reference annotation. This step is optional but advisable when the pipeline has been previously run against the same target database. However, this step must never be done when the target database is one of the production databases.

### 3. Split the input annotation file

The input gtf/gff3 file is normally split in smaller input files to be processed independently. Each smaller input file contains one or more genomic clusters of transcripts (a cluster may contain transcripts on both strands but no transcriptless gap). The original input file doesn't have to be sorted by genomic location: the pipeline will sort and group the input transcripts in clusters.

### 4. Annotation

The annotation in each individual input file is parsed and converted into Bio::Vega objects. The pipeline then searches the reference annotation database for possible host genes for each input transcript and determines which gene is the most suitable one. If there is no host gene for a transcript, it will create a new gene for that transcript. At the same time, several filters are applied: input transcripts can be rejected based on the biotype or protected status of the chosen host gene, or on the number of genes with exonic overlap with the transcript.

Input transcripts are subsequently compared with each existing transcript in the host gene and the pipeline determines if the input transcript provides any exon or intron novelty and whether both transcripts could be merged. Based on all these comparisons, the pipeline chooses whether to merge the input transcript with an existing transcript or to add the input transcript as a novel transcript.

Then, the pipeline attempts to create a CDS if the host gene is coding, and chooses an adequate biotype for the transcript.

(...)

Eventually, the pipeline will try to store the new annotation (one gene at a time) in the target database. As this process uses the Otter manual annotation tool logic, it will try to lock the clones corresponding to the genome region being updated. If the clones can't be locked, the annotation won't succeed. In this case, depending on the chosen pipeline options, the job can either finish "successfully" (i.e. fail silently) or fail with a proper error message.

### 5. Stats

When all the annotation jobs are done, the pipeline will extract some statistics of the updated annotation from the log files and the target database.

### 6. Track hub

The pipeline makes bigBed files for the input annotation and the updated annotation and creates a track hub directory containing those files. This track hub directory can then be uploaded to our ftp site (manually).


## Running the pipeline

Run the pipeline in a tmux or screen session.

### Software environment
The PERL5LIB and PATH variables must include the paths to the software dependencies as described above.

### Initialisation

Initialise the pipeline using the LoutreWrite::Hive::LoutreWrite_conf module.
```sh
init_pipeline.pl LoutreWrite::Hive::LoutreWrite_conf \
  [OPTIONS] \
  --hive_force_init 1
```


### Pipeline options
| Name | Description | Default |
| ------ | ------ | ------ |
| pipe_db_host | eHive pipeline database host	| |
| pipe_db_port	| eHive pipeline database port	| |
| pipe_db_user	| eHive pipeline database user	| |
| pipeline_name	| eHive pipeline database name	| |
| species	| species name	| |
| dataset	| species name assigned to the target database in server-config/species.dat	| |
| reset_target_db	| set to '1' if the target database should be replaced with a copy of a 'source' reference database (usually done for test databases, NEVER for production databases)	| 0 |
| source_db_name	| if reset_target_db was set to '1', name of the database that will be copied	| |
| target_db_name	| reference annotation database name	| |
| output_dir 	| directory where all the output files (split input files, log files, stats files, etc) will be stored	| |
| input_file	| input annotation file in GTF or GFF3 format (zipped files accepted)	| |
| max_clus	| maximum number of transcript clusters per input file	| 1 |
| chr	| only annotate this chromosome	| |
| tag	| prefix for output files...	| |
| host_biotype	| accepted host gene biotypes, i.e. genes where input transcripts can be added (comma-separated list)	| |
| no_overlap_biotype	| input transcripts that overlap these gene biotypes are rejected (comma-separated list)	| |
| protected_loci	| path to a file containing a list of gene stable ids that cannot be edited by the pipeline, so input transcripts that would be added to these genes are rejected	| |
| no_artifact_check	| | |
| no_intron_check	| do not apply filters in the annotation step	| |
| max_ov_loc	| maximum number of genes that an input transcript is allowed to overlap at the exon level without being rejected	| |
| no_NFV	| do not add a 'not for VEGA' transcript remark to new transcripts	| 0 |
| platinum	| add a 'Platinum' transcript remark	| 0 |
| remark	| transcript remark to be added to all transcripts created or extended in the current pipeline run	| 'Assembled using long reads' |
| no_CDS	| do not try to add a CDS to any transcript	| 0 |
| registry	| Ensembl registry file containing connection details of the auxiliary databases (recount3, polyA sites, etc)	| |
| source	| file containing a list of read names supporting each input transcript	| |
| die_locked	| let a pipeline job die if the region to be edited is locked in the target database (this requires manual intervention to unlock the region but the pipeline would not be aware of incomplete jobs otherwise)	| |
| scripts_dir	| path to auxiliary scripts	| |
| job_limit	| maximum number of jobs allowed to run simultaneously	| 20 |
| write	| edit the annotation in the target database	| |
| prefilter_input_file	| set to '1' if filters must be run at the start of the pipeline	| 0 |
| intron_score_file	| path to file with Recount3 intron scores (in case an equivalent database is not provided)	| |
| intron_score_cutoff	| minimum intron score required to pass the filter	| 50 |
| repeat_feature_file	| path to file with repeat elements (in case an equivalent database is not provided)	| |
| bam_directory	| path to directory with read alignments, required for the splice site misalignment filter	| |
| genome_fasta_file	| path to genome sequence file, required for the splice site misalignment filter	| |
| read_support_file	| path to supporting read file, required for the splice site misalignment filter (not the same file as in the -source parameter)	| |
| filter_file	| file containing a list of input transcripts that failed the filters and should be ignored (useful if the filters were run previously)	| |

NOTE: all these parameters require a value (default: '1') or eHive will parse them incorrectly.

The init_pipeline.pl command will create the pipeline database. Note that the final pipeline database name will be the value of the -pipeline_name option, prefixed with the user name and suffixed with 'pipe' (see example below).



### Start the pipeline

Start the pipeline providing the URL of the pipeline database.
```sh
PIPELINE_DATABASE_URL="mysql://"${DBUSER}":"${DBPASS}"@"${DBHOST}":"${DBPORT}"/"${DBNAME}
beekeeper.pl \
  -url $PIPELINE_DATABASE_URL \
  -loop
```


### GuiHive

Check the pipeline status in a web browser using a URL with the following pattern:
```sh
http://guihive.ebi.ac.uk:8080/versions/97/?driver=mysql&username=${DBUSER}&host=${DBHOST}&port=${DBPORT}&dbname=${DBNAME}&passwd=xxxxx
```





