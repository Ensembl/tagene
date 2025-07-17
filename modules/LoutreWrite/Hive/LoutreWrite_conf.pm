package LoutreWrite::Hive::LoutreWrite_conf;

use strict;
use warnings;


use File::Spec::Functions qw(catfile catdir);
use DateTime;

use parent ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');


=head2 default_options

 Arg [1]    : None
 Description: Create default hash for this configuration file
 Returntype : Hashref
 Exceptions : None

=cut

sub default_options {
  my ($self) = @_;

  my $date = DateTime->now()->strftime("%Y_%m_%d");
  my $date2 = DateTime->now()->strftime("%Y-%m-%d");
  my $date3 = DateTime->now()->strftime("%Y%m%d");

  return {
    %{$self->SUPER::default_options()},

    user_r => 'ensro',
    user => 'ensadmin',
    password => $ENV{EHIVE_PASS},
    pipe_db_host => 'mysql-ens-havana-prod-2',
    pipe_db_port => 4682,    
    pipeline_name => join('_', 'loutre_write', $self->o('dataset'), $date, $self->o('db_prefix')),
    email => $ENV{'USER'}.'@ebi.ac.uk',
    hive_capacity => 0, # 0 means unlimited

    species => 'homo_sapiens',
    dataset => 'human_tagene_test_6',
    db_prefix => 1,
    input_file => $self->o('input_file'),

    reset_target_db => 0,
    source_db_name => '',
    source_db_host => 'mysql-ens-havana-prod-2',
    source_db_port => 4682,
    source_db_user => $self->o('user_r'),
    source_db_driver => $self->o('hive_driver'),
    source_db => {
      -dbname => $self->o('source_db_name'),
      -host   => $self->o('source_db_host'),
      -port   => $self->o('source_db_port'),
      -user   => $self->o('source_db_user'),
      -driver => $self->o('source_db_driver'),
    },
    target_db_name => '',
    target_db_host => 'mysql-ens-havana-prod-2',
    target_db_port => 4682,
    target_db_user => $self->o('user'),
    target_db_password => $self->o('password'),
    target_db_driver => $self->o('hive_driver'),
    target_db => {
      -dbname => $self->o('target_db_name'),
      -host   => $self->o('target_db_host'),
      -port   => $self->o('target_db_port'),
      -user   => $self->o('target_db_user'),
      -pass   => $self->o('target_db_password'),
      -driver => $self->o('target_db_driver'),
    },
    
    author => 'tagene',
    source => '',
    remark => 'Assembled using long reads',
    read_seq_dir => '',
    comp_pipe => 0,
    analysis => '',
    tsource => '',
    assembly_version => '',
    no_NFV => 0,
    no_CDS => 0,
    no_intron_check => 0,
    host_biotype => '',
    no_overlap_biotype => '',
    tr_biotypes => '',
    complete_cds => 0,
    mane_select_start => 0,
    mane_select_stop => 0,
    no_novel_genes => 0,
    protected_loci => '',
    protected_regions => '',
    max_ov_loc => 0,
    no_extensions => 0,
    filter_introns => 0,
    platinum => 0,
    die_locked => 0,
    chr => '',
    write => 0,
    max_trs => 0,
    max_clus => 1,
    max_files => 0,
    tag => 'split',
    dep_job => '',
    output_dir => '',
    filter_file => 'all_failed_non_redundant.txt',
    by_gene_id => 0,
    no_shuffle => 0,
    job_limit => 20,

    prefilter_input_file => 0,
    intron_score_file => '',
    intron_score_cutoff => 50,
    repeat_feature_file => '',
    bam_directory => '',
    genome_fasta_file => '',
    read_support_file => '',
    
    scripts_dir => '',
    recount3_score_script => catfile($self->o('scripts_dir'), 'check_recount3_score.pl'),
    repeat_overlap_script => catfile($self->o('scripts_dir'), 'check_repeat_overlap.pl'),
    pseudogene_overlap_script => catfile($self->o('scripts_dir'), 'check_pseudogene_overlap.pl'),
    opp_strand_mismap_script => catfile($self->o('scripts_dir'), 'check_opposite_strand_mismapping.pl'),
    splice_site_misali_script => catfile($self->o('scripts_dir'), 'check_splice_site_alignments.pl'),
    load_gxf_script => catfile($self->o('scripts_dir'), 'load_gxf_in_loutre.pl'),
    log_stats_script => catfile($self->o('scripts_dir'), 'report'),
    db_stats_script => catfile($self->o('scripts_dir'), 'stats.pl'),
    trackhub_script => catfile($self->o('scripts_dir'), 'make_trackhub_directory.pl'),
    annot_to_bgp_script => catfile($self->o('scripts_dir'), 'make_bigBed_file_for_trackhub.pl'),
    gtf_to_bgp_script => catfile($self->o('scripts_dir'), 'convert_gtf_to_bigGenePred.pl'),
    date => $date2,
  };
}

sub pipeline_wide_parameters {
  my ($self) = @_;
  return {
    %{$self->SUPER::pipeline_wide_parameters},
  }
}

sub pipeline_analyses {
  my ($self) = @_;

  #Check that 'reset_target_db' is not being used with a production database
  if (($self->o('dataset') =~ /^(human|mouse)$/ or $self->o('target_db') =~ /^havana_(human|mouse)$/) and $self->o('reset_target_db') == 1){
    die "Attempting to reset the ".$self->o('dataset')." production database!";
  }
  
  my @analyses = (
    {
      -logic_name => 'create_working_directory',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'for DIR in #file_dir# #log_dir# #th_dir# #filter_dir#; do mkdir -p $DIR; done',
        use_bash_errexit => 1,
        file_dir => catdir($self->o('output_dir'), 'gxf'),
        log_dir => catdir($self->o('output_dir'), 'log'),
        th_dir => catdir($self->o('output_dir'), 'trackhub'),
        filter_dir => catdir($self->o('output_dir'), 'filters'),
      },
      -input_ids => [{
        species => $self->o('species'),
      }],
      -flow_into => {
        '1->A' => ['do_reset_database','do_prefilter_input_file'],
        'A->1' => ['split_annotation_file'],
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'do_reset_database',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'if [ #reset_target_db# -eq 0 ]; then exit 42; fi',
        return_codes_2_branches => {
          42 => 2,
        },
        reset_target_db => $self->o('reset_target_db'),
      },
      -rc_name => 'default',
      -flow_into => {
        1 => ['reset_database'],
      },
    },

    {
      -logic_name => 'reset_database',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('source_db'),
        target_db => $self->o('target_db'),
        force_drop => 1,
        ignore_dna => 1,
        create_type => 'copy',
      },
      -rc_name => 'default',
      -flow_into => {
        1 => ['remove_slice_locks'],
      },
    },
    
    {
      -logic_name => 'remove_slice_locks',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('target_db'),
        sql => [
          'TRUNCATE slice_lock',
        ],
      },            
      -rc_name => 'default',
    },
    
    {
      -logic_name => 'do_prefilter_input_file',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'if [ #prefilter# -eq 0 ]; then exit 42; fi',
        return_codes_2_branches => {
          42 => 2,
        },
        prefilter => $self->o('prefilter_input_file'),
      },
      -rc_name => 'default',
      -flow_into => {
        1 => ['get_chr_list'],
      },
    },

    {
      -logic_name => 'get_chr_list',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        inputcmd => 'zcat #input_file# | grep -v ^# | cut -f1 | sort -u',
        input_file => $self->o('input_file'),
        column_names => ['chr'],
      },
      -flow_into => {
        '2->A' => ['recount3_score_filter'],
        'A->1' => ['make_filter_list'],
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'recount3_score_filter',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl #recount3_score_script#'
              .' -gtf #input_file#'
              .' -chr #chr#'
              .' -recount3 #recount3_file#'
              .' -cutoff #recount3_cutoff#'
              .' -out #output_file#',
        recount3_score_script => $self->o('recount3_score_script'),
        input_file => $self->o('input_file'),
        recount3_file => $self->o('intron_score_file'),
        recount3_cutoff => $self->o('intron_score_cutoff'),
        filter_dir => catdir($self->o('output_dir'), 'filters'),
        output_file => catfile('#filter_dir#', 'recount3_score_#chr#.txt'),
      },
      -rc_name => '5GB',
      -flow_into => {
        1 => ['repeat_overlap_filter'],
      },
    },

    {
      -logic_name => 'repeat_overlap_filter',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl #repeat_overlap_script#'
              .' -gtf #input_file#'
              .' -chr #chr#'
              .' -repeats #repeat_file#'
              .' -out #output_file#',
        repeat_overlap_script => $self->o('repeat_overlap_script'),
        input_file => $self->o('input_file'),
        repeat_file => $self->o('repeat_feature_file'),
        filter_dir => catdir($self->o('output_dir'), 'filters'),
        output_file => catfile('#filter_dir#', 'repeat_overlap_#chr#.txt'),
      },
      -rc_name => '5GB',
      -flow_into => {
        1 => ['pseudogene_overlap_filter'],
      },
    },

    {
      -logic_name => 'pseudogene_overlap_filter',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl #pseudogene_overlap_script#'
              .' -gtf #input_file#'
              .' -chr #chr#'
              .' -dataset '.$self->o('dataset')
              .' -utr5 180'
              .' -utr3 561'
              .' -out #output_file#',
        pseudogene_overlap_script => $self->o('pseudogene_overlap_script'),
        input_file => $self->o('input_file'),
        filter_dir => catdir($self->o('output_dir'), 'filters'),
        output_file => catfile('#filter_dir#', 'pseudogene_overlap_#chr#.txt'),
      },   
      -rc_name => '5GB',
      -flow_into => {
        1 => ['opposite_strand_mismapping_filter'],
      },
    },

    {
      -logic_name => 'opposite_strand_mismapping_filter',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl #opp_strand_mismap_script#'
              .' -gtf #input_file#'
              .' -chr #chr#'
              .' -dataset '.$self->o('dataset')
              .' -out #output_file#',
        opp_strand_mismap_script => $self->o('opp_strand_mismap_script'),
        input_file => $self->o('input_file'),
        filter_dir => catdir($self->o('output_dir'), 'filters'),
        output_file => catfile('#filter_dir#', 'opp_strand_mismap_#chr#.txt'),
      },
      -rc_name => '5GB',
      -flow_into => {
        1 => ['splice_site_misalignment_filter'],
      },
    },

    {
      -logic_name => 'splice_site_misalignment_filter',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl #splice_site_misali_script#'
              .' -gtf #input_file#'
              .' -chr #chr#'
              .' -bamdir #bam_directory#'
              .' -fasta #genome_fasta_file#'
              .' -readsup #read_support_file#'
              .' -window 10'
              .' -simple'
              .' -out #output_file#',
        splice_site_misali_script => $self->o('splice_site_misali_script'),
        input_file => $self->o('input_file'),
        bam_directory => $self->o('bam_directory'),
        genome_fasta_file => $self->o('genome_fasta_file'),
        read_support_file => $self->o('read_support_file'),
        filter_dir => catdir($self->o('output_dir'), 'filters'),
        output_file => catfile('#filter_dir#', 'splice_site_misali_#chr#.txt'),
      },
      -rc_name => '5GB',
    },

    {
      -logic_name => 'make_filter_list',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'cd #filter_dir#;
                cat recount3_score_*.txt | grep -P "\tno\t" | cut -f1 > failed_recount3_score.txt;
                cat repeat_overlap_*.txt | cut -f1 > failed_repeat_overlap.txt;
                cat pseudogene_overlap_*.txt | cut -f1 > failed_pseudogene_overlap.txt;
                cat opp_strand_mismap_*.txt | cut -f1 > failed_opposite_strand_mismapping.txt;
                cat splice_site_misali_*.txt | grep -P "\tno\t" | cut -f1 > failed_splice_site_misalignment.txt;
                cat failed_*.txt | sort -u > #output_file#',
        filter_dir => catdir($self->o('output_dir'), 'filters'),
        output_file => $self->o('filter_file'),
      },
      -rc_name => 'default',
    },


    {
      -logic_name => 'split_annotation_file',
      -module     => 'LoutreWrite::Hive::PreLoadGxfInLoutre',
      -parameters => {
        input_file => $self->o('input_file'),
        filter_output_file => catfile(catdir($self->o('output_dir'), 'filters'), $self->o('filter_file')),
        output_dir => catdir($self->o('output_dir'), 'gxf'),
        log_dir => catdir($self->o('output_dir'), 'log'),
        max_trs_per_file => $self->o('max_trs'),
        max_clusters_per_file => $self->o('max_clus'),
        max_num_files => $self->o('max_files'),
        tag => $self->o('tag'),
        cluster_by_gene_id => $self->o('by_gene_id'),
        no_shuffle => $self->o('no_shuffle'),
      },
      -flow_into => {
        1 => ['create_individual_jobs'],
      },
      -rc_name => '40GB',
    },

    {
      -logic_name => 'create_individual_jobs',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        inputcmd => 'find #gxf_dir# -name "*.gtf" -or -name "*.gff3" | replace #gxf_dir# / | sed "s/^\/\///"',
        column_names => ['filename'],
        gxf_dir => catdir($self->o('output_dir'), 'gxf'),
      },
      -flow_into => {
        '2->A' => ['load_gxf_in_loutre'],
        'A->1' => ['get_stats_from_logs'],
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'load_gxf_in_loutre',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl #load_gxf_script#'
          .' -file #infile#'
          .' -dataset '.$self->o('dataset')
          .' -author '.$self->o('author')
          .' #source#'
          .' -remark "'.$self->o('remark').'"'
          .' #readseqdir#'
          .' #comp_pipe#'
          .' #analysis#'
          .' #tsource#'
          .' #assembly#'
          .' #registry#'
          .' #no_nfv#'
          .' #no_cds#'
          .' #no_intron_check#'
          .' -host_biotype '.$self->o('host_biotype')
          .' -no_overlap_biotype '.$self->o('no_overlap_biotype')
          .' #tr_biotypes#'
          .' #mane_select_start#'
          .' #mane_select_stop#'
          .' #complete_cds#'
          .' #no_novel_genes#'
          .' -max_ov_loc '.$self->o('max_ov_loc')
          .' #protected_loci#'
          .' #protected_regions#'
          .' #no_extensions#'
          .' #filter_introns#'
          .' #platinum#'
          .' #die_locked#'
          .' #chr#'
          .' #write#'
          .' > #output_file#'
          .' 2> #output_file#.err',
        load_gxf_script => $self->o('load_gxf_script'),
        gxf_dir => catdir($self->o('output_dir'), 'gxf'),
        infile => catfile('#gxf_dir#', '#filename#'),
        log_dir => catdir($self->o('output_dir'), 'log'),
        output_file => catfile('#log_dir#', '#filename#.log'),
        source => ($self->o('source') ? '-source '.$self->o('source') : ''),
        readseqdir => ($self->o('read_seq_dir') ? '-readseqdir '.$self->o('read_seq_dir') : ''),
        comp_pipe => ($self->o('comp_pipe') ? '-comp_pipe' : ''),
        analysis => ($self->o('analysis') ? '-analysis '.$self->o('analysis') : ''),        
        tsource => ($self->o('tsource') ? '-tsource '.$self->o('tsource') : ''), 
        assembly => ($self->o('assembly_version') ? '-assembly_version '.$self->o('assembly_version') : ''),
        registry => ($self->o('registry') ? '-registry '.$self->o('registry') : ''),
        chr => ($self->o('chr') ? ' -chr '.$self->o('chr') : ''),
        no_intron_check => ($self->o('no_intron_check') ? '-no_intron_check ' : ''),
        no_nfv => ($self->o('no_NFV') ? '-no_NFV' : ''),
        no_cds => ($self->o('no_CDS') ? '-no_CDS' : ''),
        tr_biotypes => ($self->o('tr_biotypes') ? '-tr_biotypes '.$self->o('tr_biotypes') : ''),
        mane_select_start => ($self->o('mane_select_start') ? ' -mane_select_start' : ''),
        mane_select_stop => ($self->o('mane_select_stop') ? ' -mane_select_stop' : ''),
        complete_cds => ($self->o('complete_cds') ? '-complete_cds' : ''),
        no_novel_genes => ($self->o('no_novel_genes') ? '-no_novel_genes' : ''),
        protected_loci => ($self->o('protected_loci') ? '-protected_loci '.$self->o('protected_loci') : ''),
        protected_regions => ($self->o('protected_regions') ? '-protected_regions '.$self->o('protected_regions') : ''),
        no_extensions => ($self->o('no_extensions') ? '-no_extensions' : ''),
        filter_introns => ($self->o('filter_introns') ? '-filter_introns' : ''),
        platinum => ($self->o('platinum') ? '-platinum' : ''),
        die_locked => ($self->o('die_locked') ? '-die_locked' : ''),
        write => ($self->o('write') ? '-write' : ''),
      },
      -max_retry_count => 0,
      -rc_name => '5GB',
      -analysis_capacity => $self->o('job_limit'),
      -flow_into => {
        'MEMLIMIT' => ['load_gxf_in_loutre_high_mem'],
        'RUNLIMIT' => ['load_gxf_in_loutre_high_mem'],
      },
    },

    {
      -logic_name => 'load_gxf_in_loutre_high_mem',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl #load_gxf_script#'
          .' -file #infile#'
          .' -dataset '.$self->o('dataset')
          .' -author '.$self->o('author')
          .' #source#'
          .' -remark "'.$self->o('remark').'"'
          .' #readseqdir#'
          .' #comp_pipe#'
          .' #analysis#'
          .' #tsource#'
          .' #assembly#'
          .' #registry#'
          .' #no_nfv#'
          .' #no_cds#'
          .' #no_intron_check#'
          .' -host_biotype '.$self->o('host_biotype')
          .' -no_overlap_biotype '.$self->o('no_overlap_biotype')
          .' #tr_biotypes#'
          .' #mane_select_start#'
          .' #mane_select_stop#'
          .' #complete_cds#'
          .' #no_novel_genes#'
          .' -max_ov_loc '.$self->o('max_ov_loc')
          .' #protected_loci#'
          .' #protected_regions#'
          .' #no_extensions#'
          .' #filter_introns#'
          .' #platinum#'
          .' #die_locked#'
          .' #chr#'
          .' #write#'
          .' > #output_file#'
          .' 2> #output_file#.err',
        load_gxf_script => $self->o('load_gxf_script'),
        gxf_dir => catdir($self->o('output_dir'), 'gxf'),
        infile => catfile('#gxf_dir#', '#filename#'),
        log_dir => catdir($self->o('output_dir'), 'log'),
        output_file => catfile('#log_dir#', '#filename#.log'),
        source => ($self->o('source') ? '-source '.$self->o('source') : ''),
        readseqdir => ($self->o('read_seq_dir') ? '-readseqdir '.$self->o('read_seq_dir') : ''),
        comp_pipe => ($self->o('comp_pipe') ? '-comp_pipe' : ''),
        analysis => ($self->o('analysis') ? '-analysis '.$self->o('analysis') : ''),        
        tsource => ($self->o('tsource') ? '-tsource '.$self->o('tsource') : ''), 
        assembly => ($self->o('assembly_version') ? '-assembly_version '.$self->o('assembly_version') : ''),
        registry => ($self->o('registry') ? '-registry '.$self->o('registry') : ''),
        chr => ($self->o('chr') ? ' -chr '.$self->o('chr') : ''),
        no_intron_check => ($self->o('no_intron_check') ? '-no_intron_check ' : ''),
        no_nfv => ($self->o('no_NFV') ? '-no_NFV' : ''),
        no_cds => ($self->o('no_CDS') ? '-no_CDS' : ''),
        tr_biotypes => ($self->o('tr_biotypes') ? '-tr_biotypes '.$self->o('tr_biotypes') : ''),
        mane_select_start => ($self->o('mane_select_start') ? ' -mane_select_start' : ''),
        mane_select_stop => ($self->o('mane_select_stop') ? ' -mane_select_stop' : ''),
        complete_cds => ($self->o('complete_cds') ? '-complete_cds' : ''),
        no_novel_genes => ($self->o('no_novel_genes') ? '-no_novel_genes' : ''),
        protected_loci => ($self->o('protected_loci') ? '-protected_loci '.$self->o('protected_loci') : ''),
        protected_regions => ($self->o('protected_regions') ? '-protected_regions '.$self->o('protected_regions') : ''),
        no_extensions => ($self->o('no_extensions') ? '-no_extensions' : ''),
        filter_introns => ($self->o('filter_introns') ? '-filter_introns' : ''),
        platinum => ($self->o('platinum') ? '-platinum' : ''),
        die_locked => ($self->o('die_locked') ? '-die_locked' : ''),
        write => ($self->o('write') ? '-write' : ''),
      },
      -max_retry_count => 0,
      -rc_name => '20GB',
      -analysis_capacity => $self->o('job_limit'),
    },

    {
      -logic_name => 'get_stats_from_logs',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'bash #stats_script#'     
          .' #log_dir#'
          .' > #output_file#',
        stats_script => $self->o('log_stats_script'),
        log_dir => catdir($self->o('output_dir'), 'log'),
        output_file => catfile($self->o('output_dir'), join('_', 'log', 'stats', $self->o('dataset'), $self->o('date'))),
      },
      -max_retry_count => 0,
      -flow_into => {
        1 => ['get_stats_from_db'],
      },
      -rc_name => '5GB',
    },

    {
      -logic_name => 'get_stats_from_db',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl #stats_script#'     
          .' -dataset '.$self->o('dataset')
          .' -date '.$self->o('date')
          .' -remark "'.$self->o('remark').'"'
          .' -out #output_file#',
        stats_script => $self->o('db_stats_script'),
        output_file => catfile($self->o('output_dir'), join('_', 'db', 'stats', $self->o('dataset'), $self->o('date'))),
      },
      -max_retry_count => 0,
      -flow_into => {
        1 => ['make_trackhub_annotation_file'],
      },
      -rc_name => '5GB',
    },  

    {
      -logic_name => 'make_trackhub_annotation_file',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'cd #outdir# && perl #annot_to_bgp_script#'     
          .' -host '.$self->o('target_db_host')
          .' -port '.$self->o('target_db_port')
          .' -user '.$self->o('user_r')
          .' -dbname '.$self->o('target_db_name')
          .' -date '.$self->o('date')
          .' -colour1 22,74,216'
          .' -colour2 220,21,21'
          .' -out #output_file#',
        annot_to_bgp_script => $self->o('annot_to_bgp_script'),
        outdir => catdir($self->o('output_dir'), "trackhub"),
        output_file => catfile('#outdir#', $self->o('dataset')),
      },
      -max_retry_count => 0,
      -flow_into => {
        1 => ['make_trackhub_input_file'],
      },
      -rc_name => '5GB',
    },

    {
      -logic_name => 'make_trackhub_input_file',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'cd #outdir# && perl #gtf_to_bgp_script#'
          .' -gtf '.$self->o('input_file')
          .' -host '.$self->o('target_db_host')
          .' -port '.$self->o('target_db_port')
          .' -user '.$self->o('user_r')
          .' -dbname '.$self->o('target_db_name')
          .' -out #output_file#',
        gtf_to_bgp_script => $self->o('gtf_to_bgp_script'),
        outdir => catdir($self->o('output_dir'), "trackhub"),
        output_file => catfile('#outdir#', "input_annot.bb"),
      },
      -max_retry_count => 0,
      -flow_into => {
        1 => ['make_trackhub_directory'],
      },
      -rc_name => '5GB',
    },

    {
      -logic_name => 'make_trackhub_directory',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'cd #outdir# && perl #trackhub_script#'
          .' -name '.$self->o('dataset').'_'.$self->o('date')
          .' -host '.$self->o('target_db_host')
          .' -port '.$self->o('target_db_port')
          .' -user '.$self->o('user_r')
          .' -dbname '.$self->o('target_db_name')
          .' -input #input_file#'
          .' -annot #annot_file#',
        trackhub_script => $self->o('trackhub_script'),
        outdir => catdir($self->o('output_dir'), "trackhub"),
        input_file => catfile('#outdir#', "input_annot.bb"),
        annot_file => catfile('#outdir#', $self->o('dataset').".bb"),
      },
      -max_retry_count => 0,
      -rc_name => 'default',
    },

  );

  foreach my $analysis (@analyses) {
    $analysis->{'-max_retry_count'} = 0 unless (exists $analysis->{'-max_retry_count'});
    $analysis->{'-hive_capacity'} = $self->o('hive_capacity') if ($self->o('hive_capacity'));
    if ($analysis->{'-logic_name'} eq 'pseudogene_overlap_filter'){
      if ($self->o('prefilter_input_file') and $self->o('reset_target_db')){
        $analysis->{'-wait_for'} = 'reset_database';
      }
    }
  }

  return \@analyses;
}


sub resource_classes {
  my ($self) = @_;

  return {
    'default' => {SLURM => '--time=3:00:00 --mem=1000', LSF => '-q short -M 1000 -R"select[mem>1000] rusage[mem=1000]"'},   
    '5GB' => {SLURM => '--time=6:00:00 --mem=5000', LSF => '-q short -M 5000 -R"select[mem>5000] rusage[mem=5000]"'},
    '10GB' => {SLURM => '--time=8:00:00 --mem=10000', LSF => '-q standard -M 10000 -R"select[mem>10000] rusage[mem=10000]"'},
    '20GB' => {SLURM => '--time=10:00:00 --mem=20000', LSF => '-q standard -M 20000 -R"select[mem>20000] rusage[mem=20000]"'},
    '40GB' => {SLURM => '--time=12:00:00 --mem=40000', LSF => '-q short -M 40000 -R"select[mem>40000] rusage[mem=40000]"'},
  };
}

1;

