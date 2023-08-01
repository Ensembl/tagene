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
    no_check => 0,
    no_NFV => 0,
    no_CDS => 0,
    no_artifact_check => 0,
    no_intron_check => 0,
    host_biotype => '',
    no_overlap_biotype => '',
    max_ov_loci => 0,
    filter_introns => 0,
    platinum => 0,
    chr => '',
    write => 0,
    max_trs => 0,
    max_clus => 1,
    max_files => 0,
    tag => 'split',
    dep_job => '',
    output_dir => '',
    by_gene_id => 0,
    no_shuffle => 0,
    job_limit => 20,

    scripts_dir => '',
    load_gxf_script => catfile($self->o('scripts_dir'), 'load_gxf_in_loutre.pl'),
    log_stats_script => catfile($self->o('scripts_dir'), 'report'),
    db_stats_script => catfile($self->o('scripts_dir'), 'stats.pl'),
    trackhub_script => catfile($self->o('scripts_dir'), 'make_trackhub_files.pl'),
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

  my @analyses = (
    {
      -logic_name => 'create_working_directory',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'for DIR in #file_dir# #log_dir# #th_dir# ; do mkdir -p $DIR; done',
        use_bash_errexit => 1,
        file_dir => catdir($self->o('output_dir'), 'gxf'),
        log_dir => catdir($self->o('output_dir'), 'log'),
        th_dir => catdir($self->o('output_dir'), 'trackhub'),
      },
      -input_ids => [{
        species => $self->o('species'),
      }],
      -flow_into => {
        #1 => ['fan_reset_database'],
        '1->A' => ['do_reset_database'],
        'A->1' => ['split_annotation_file'],
      },
      -rc_name => 'default',
    },

    #{
      #-logic_name => 'fan_reset_database',
      #-module => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      #-parameters => {
      #},
      #-rc_name => 'default',
      #-flow_into => {
        #'2->A' => ['do_reset_database'],
        #'A->1' => ['split_annotation_file'],
      #},
    #},

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
      -logic_name => 'split_annotation_file',
      -module     => 'LoutreWrite::Hive::PreLoadGxfInLoutre',
      -parameters => {
        input_file => $self->o('input_file'),
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
      -rc_name => '2GB',
    },

    {
      -logic_name => 'create_individual_jobs',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        inputcmd => 'find #gxf_dir# -name "*.gtf" -or -name "*.gff3" | replace #gxf_dir# / | sed "s/^\/\///"',
        #inputcmd => 'ls -1 #gxf_dir#',
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
       #   .' -source '.$self->o('source')
          .' -remark "'.$self->o('remark').'"'
       #   .' -readseqdir '.$self->o('read_seq_dir')
          .' #readseqdir#'
       #   .' -comp_pipe '.$self->o('comp_pipe')
       #   .' -analysis '.$self->o('analysis')
       #   .' -tsource '.$self->o('tsource')
       #   .' -assembly '.$self->o('assembly_version')
          .' #assembly#'
       #   .' -no_check '.$self->o('no_check')
          .' #no_artifact_check#'
       #   .' -no_NFV '.$self->o('no_NFV')
          .' #no_nfv#'
       #   .' -no_CDS '.$self->o('no_CDS')
          .' #no_cds#'
       #   .' -no_intron_check '.$self->o('no_intron_check')
          .' -host_biotype '.$self->o('host_biotype')
          .' -no_overlap_biotype '.$self->o('no_overlap_biotype')
          .' -max_ov_loc '.$self->o('max_ov_loci')
      #    .' -filter_introns '.$self->o('filter_introns')
          .' #filter_introns#'
          #.' -platinum '.$self->o('platinum')
          .' #platinum#'
      #    .' -chr '.$self->o('chr')
          .' #chr#'
          #.' -write '.$self->o('write')
          .' #write#'
          .' > #output_file#'
          .' 2> #output_file#.err',
        load_gxf_script => $self->o('load_gxf_script'),
        gxf_dir => catdir($self->o('output_dir'), 'gxf'),
        infile => catfile('#gxf_dir#', '#filename#'),
        log_dir => catdir($self->o('output_dir'), 'log'),
        output_file => catfile('#log_dir#', '#filename#.log'),
        #log_dir_2 => join("", grep {$_ =~ /\//} splitpath('#output_file#')),
        readseqdir => ($self->o('read_seq_dir') ? '-readseqdir '.$self->o('read_seq_dir') : ''),
        assembly => ($self->o('assembly_version') ? '-assembly_version '.$self->o('assembly_version') : ''),
        chr => ($self->o('chr') ? ' -chr '.$self->o('chr') : ''),
        no_artifact_check => ($self->o('no_artifact_check') ? ' -no_check ' : ''),
        no_nfv => ($self->o('no_NFV') ? '-no_NFV' : ''),
        no_cds => ($self->o('no_CDS') ? '-no_CDS' : ''),
        filter_introns => ($self->o('filter_introns') ? '-filter_introns' : ''),
        platinum => ($self->o('platinum') ? '-platinum' : ''),
        write => ($self->o('write') ? '-write' : ''),
      },
      -max_retry_count => 0,
      -rc_name => 'default',
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
       #   .' -source '.$self->o('source')
          .' -remark "'.$self->o('remark').'"'
       #   .' -readseqdir '.$self->o('read_seq_dir')
          .' #readseqdir#'
       #   .' -comp_pipe '.$self->o('comp_pipe')
       #   .' -analysis '.$self->o('analysis')
       #   .' -tsource '.$self->o('tsource')
       #   .' -assembly '.$self->o('assembly_version')
          .' #assembly#'
       #   .' -no_check '.$self->o('no_check')
          .' #no_artifact_check#'
       #   .' -no_NFV '.$self->o('no_NFV')
          .' #no_nfv#'
       #   .' -no_CDS '.$self->o('no_CDS')
          .' #no_cds#'
       #   .' -no_intron_check '.$self->o('no_intron_check')
          .' -host_biotype '.$self->o('host_biotype')
          .' -no_overlap_biotype '.$self->o('no_overlap_biotype')
          .' -max_ov_loc '.$self->o('max_ov_loci')
      #    .' -filter_introns '.$self->o('filter_introns')
          .' #filter_introns#'
          #.' -platinum '.$self->o('platinum')
          .' #platinum#'
      #    .' -chr '.$self->o('chr')
          .' #chr#'
          #.' -write '.$self->o('write')
          .' #write#'
          .' > #output_file#'
          .' 2> #output_file#.err',
        load_gxf_script => $self->o('load_gxf_script'),
        gxf_dir => catdir($self->o('output_dir'), 'gxf'),
        infile => catfile('#gxf_dir#', '#filename#'),
        log_dir => catdir($self->o('output_dir'), 'log'),
        output_file => catfile('#log_dir#', '#filename#.log'),
        #log_dir_2 => join("", grep {$_ =~ /\//} splitpath('#output_file#')),
        readseqdir => ($self->o('read_seq_dir') ? '-readseqdir '.$self->o('read_seq_dir') : ''),
        assembly => ($self->o('assembly_version') ? '-assembly_version '.$self->o('assembly_version') : ''),
        chr => ($self->o('chr') ? ' -chr '.$self->o('chr') : ''),
        no_artifact_check => ($self->o('no_artifact_check') ? ' -no_check ' : ''),
        no_nfv => ($self->o('no_NFV') ? '-no_NFV' : ''),
        no_cds => ($self->o('no_CDS') ? '-no_CDS' : ''),
        filter_introns => ($self->o('filter_introns') ? '-filter_introns' : ''),
        platinum => ($self->o('platinum') ? '-platinum' : ''),
        write => ($self->o('write') ? '-write' : ''),
      },
      -max_retry_count => 0,
      -rc_name => '6GB',
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
      -rc_name => 'default',
    },

    {
      -logic_name => 'get_stats_from_db',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl #stats_script#'     
          .' -dataset '.$self->o('dataset')
          .' -date '.$self->o('date')
          .' -out #output_file#',
        stats_script => $self->o('db_stats_script'),
        output_file => catfile($self->o('output_dir'), join('_', 'db', 'stats', $self->o('dataset'), $self->o('date'))),
      },
      -max_retry_count => 0,
      -flow_into => {
        1 => ['make_trackhub_files'],
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'make_trackhub_files',
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
        trackhub_script => $self->o('annot_to_bgp_script'),
        outdir => catdir($self->o('output_dir'), "trackhub"),
        input_file => catfile('#outdir#', "input_annot.bb"),
        annot_file => catfile('#outdir#', $self->o('dataset')),
      },
      -max_retry_count => 0,
      -flow_into => {
        1 => ['make_trackhub_annotation_file'],
      },
      -rc_name => 'default',
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
      -rc_name => 'default',
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
          .' -species '.$self->o('species')
          .' -out #output_file#',
        gtf_to_bgp_script => $self->o('gtf_to_bgp_script'),
        outdir => catdir($self->o('output_dir'), "trackhub"),
        output_file => catfile('#outdir#', "input_annot.bb"),
      },
      -max_retry_count => 0,
      -rc_name => 'default',
    },

  );

  foreach my $analysis (@analyses) {
    $analysis->{'-max_retry_count'} = 0 unless (exists $analysis->{'-max_retry_count'});
    $analysis->{'-hive_capacity'} = $self->o('hive_capacity') if ($self->o('hive_capacity'));
  }

  return \@analyses;
}


sub resource_classes {
  my ($self) = @_;

  return {
    'default' => { LSF => '-q short -M 1000 -R"select[mem>1000] rusage[mem=1000]"'},
    '2GB' => { LSF => '-q short -M 2000 -R"select[mem>2000] rusage[mem=2000]"'},    
    '4GB' => { LSF => '-q short -M 4000 -R"select[mem>4000] rusage[mem=4000]"'},
    '6GB' => { LSF => '-q standard -M 6000 -R"select[mem>6000] rusage[mem=6000]"'},
  };
}

1;

