# README

This software allows the user to import gene annotation from a gtf or gff3 file into a database with Ensembl/Otter schema.

The pipeline:

* takes gene annotation from an input file (in gtf or gff3 format) and a reference annotation database (with Ensembl/Otter schema);
* filters the input annotation to remove transcripts containing probably spurious splice sites;
* compares the input annotation with the reference annotation;
* stores in the database the input annotation that provides some novelty, either by adding novel transcripts or by extending existing transcripts.


## Dependencies

* ensembl (https://github.com/Ensembl/ensembl)
* ensembl-otter (https://github.com/Annosoft/ensembl-otter)
* server-config (https://github.com/Ensembl/server-config)
* ensembl-hive (https://github.com/Ensembl/ensembl-hive)
* try-tiny (https://github.com/nothingmuch/try-tiny)


## Installation

Your PERL5LIB variable must include the paths to the following:
* loutre_write/modules
* ensembl/modules
* ensembl-otter/modules
* ensembl-hive/modules

Your PATH variable must include the path to:
* ensembl-hive/scripts

**TO DO**: describe how to set up a basic Otter client configuration so that the pipeline can pick up dataset names from server-config/species.dat.
