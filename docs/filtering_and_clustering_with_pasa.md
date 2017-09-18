
## Summary
This pipeline takes one or more bam files with long reads (from SLRseq or PacBio), filters out spurious reads and makes use of the PASA pipeline to assemble the remaining reads into transcript models and cluster these into loci, generating a gff3 file that can be fed into the LoutreWrite software.

## Dependencies
- PASA (https://pasapipeline.github.io)

  NOTE: this summary doesn't cover the configuration of PASA pipeline. See https://pasapipeline.github.io/ for instructions.


## Pipeline

### 1. Filter reads and print gff3 files individually
  The script *convert_bam_to_gff3_for_pasa.pl* reads the bam file, filters the reads (by mapping quality, exon size, ...) and writes the gff3 output file. The format of this file is tweaked so that it mimicks the output of the alignment step of the PASA pipeline. 
 
    bsub -M12000 -R"select[mem>12000] rusage[mem=12000]" -oo f.file.out \
    perl convert_bam_to_gff3_for_pasa.pl \
      -in_mono \
      -clip_ends \
      -qual 0 \
      -bam ftp://path/to/file.bam \
      -fasta /path/to/genome.fa \
      -out file.gff3
 
  NOTE: the script was coded to deal with the slightly different outputs of the GMAP and STAR aligners, but in its current status there is no option to switch between them. It worked by commenting certain lines in and out. THIS HAS TO BE FIXED.


### 2. Merge gff3 files
  If several bam files are processed in step 1, all the output gff3 files must be merged into a single one.

    for file in *.gff3; do cat $file >> merged.gff3; done


### 3. Write non-redundant read sequence fasta file
  It is probably easier to start from the original bam files, as it doesn't matter if the sequences of the discarded reads are present in the fasta file.

    for file in [bam_file_list]; do
      samtools view ftp://path/to/$file.s.bam | cut -f1,10 | \
      perl -ne '($name, $seq) = split; print ">$name\n$seq\n" unless $seen{$name}; $seen{$name} = 1' >> merged_reads_seq.fasta
    done


### 4. Make or edit the PASA config file (alignAssembly.config)
  Indicate the PASA database name (MYSQLDB variable) and set the PASA alignment filtering paramaters.

  For example:
  
    #script validate_alignments_in_db.dbi
    validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=75
    validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=95
    validate_alignments_in_db.dbi:--NUM_BP_PERFECT_SPLICE_BOUNDARY=0

    #script subcluster_builder.dbi
    subcluster_builder.dbi:-m=50


### 5. Run the PASA pipeline
  This step only runs the PASA alignment assembly, as nothing else is needed.
  
    bsub -M10000 -R"select[mem>10000] rusage[mem=10000]" -q basement -oo f.pasa.out \
    $PASAHOME/scripts/Launch_PASA_pipeline.pl \
      -c alignAssembly.config \
      -C \
      -R \
      -g /path/to/genome.fa \
      -t merged_reads_seq.fasta \
      --IMPORT_CUSTOM_ALIGNMENTS_GFF3 merged.gff3

  NOTE: this job may take 2-3 days


### 6. Intron QC and stats
  The PASA pipeline produces several sets of files in bed, gtf and gff3 format.
   -failed_custom_alignments.bed
   -valid_custom_alignments.bed
   -pasa_assemblies.bed

  Check how many reads were discarded by PASA by counting the bed file lines.

    wc -l *.bed
  
  The *check_gff3_introns.pl* script generates a breakdown of introns by splice site sequence in the gff3 alignment assembly files.

    perl check_gff3_introns.pl \
      -i failed_custom_alignments.gff3 \
     -out intron_qc.failed \
      -asv GRCh37

    cut -f2 intron_qc.failed | sort | uniq -c | sort -nr -k1,1 | head


## For more information...
* ~/long_read_pipeline/README
* ~/long_read_pipeline/README_SLRseq
* /nfs/th_group/jmg/long_read_pipeline/using_pasa/SLRseq/README\*
* /nfs/th_group/jmg/long_read_pipeline/using_pasa/Capture_seq/README\*
* /nfs/th_group/jmg/long_read_pipeline/using_pasa/mouse_Capture_seq/README\*








