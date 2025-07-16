package LoutreWrite::Hive::Prefilter3_conf;

use strict;
use warnings;


use File::Spec::Functions qw(catfile catdir splitpath);
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

  return {
    %{$self->SUPER::default_options()},

    user_r => 'ensro',
    user => 'ensadmin',
    password => $ENV{EHIVE_PASS},
    pipe_db_host => 'mysql-ens-havana-prod-2',
    pipe_db_port => 4682,    
    pipeline_name => join('_', 'prefilter', $self->o('dataset'), $date, $self->o('db_prefix')),
    email => $ENV{'USER'}.'@ebi.ac.uk',
    hive_capacity => 0, # 0 means unlimited

    species => 'homo_sapiens',
    dataset => 'human_tagene_test_6',
    db_prefix => 1,
    input_file => $self->o('input_file'),

    ref_db_name => '',
    ref_db_host => 'mysql-ens-havana-prod-2',
    ref_db_port => 4682,
    ref_db_user => $self->o('user'),
    ref_db_driver => $self->o('hive_driver'),
    ref_db => {
      -dbname => $self->o('ref_db_name'),
      -host   => $self->o('ref_db_host'),
      -port   => $self->o('ref_db_port'),
      -user   => $self->o('ref_db_user'),
      -driver => $self->o('ref_db_driver'),
    },
    
    chr => '',
    
    output_dir => '',
    recount3_score_dir => 'recount3_score',
    repeat_overlap_dir => 'repeat_overlap',
    canonical_splice_site_dir => 'canonical_splice_site',
    pseudogene_overlap_dir => 'pseudogene_overlap',
    opp_strand_mismap_dir => 'opp_strand_mismap',
    mane_exon_overlap_dir => 'mane_exon_overlap',
    terminal_exon_size_dir => 'terminal_exon_size',
    splice_site_alignment_dir => 'splice_site_alignment',
    
    filter_file => '',

    intron_score_file => '',
    intron_score_cutoff => 50,
    repeat_feature_file => '',
    bam_directory => '',
    genome_fasta_file => '',
    read_support_file => '',
    allowed_splice_site_seq => '',
    min_terminal_exon_size => 15,
    
    scripts_dir => '',
    recount3_score_script => catfile($self->o('scripts_dir'), 'check_recount3_score.pl'),
    repeat_overlap_script => catfile($self->o('scripts_dir'), 'check_repeat_overlap.pl'),
    canonical_splice_site_script => catfile($self->o('scripts_dir'), 'check_splice_site_sequences.pl'),
    pseudogene_overlap_script => catfile($self->o('scripts_dir'), 'check_pseudogene_overlap.pl'),
    opp_strand_mismap_script => catfile($self->o('scripts_dir'), 'check_opposite_strand_mismapping.pl'),
    mane_exon_overlap_script => catfile($self->o('scripts_dir'), 'check_MANE_exon_overlap.pl'),
    terminal_exon_size_script => catfile($self->o('scripts_dir'), 'check_terminal_exons.pl'),
    splice_site_alignment_script => catfile($self->o('scripts_dir'), 'check_splice_site_alignments.pl'),
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
        cmd => 'for DIR in #filter_dir#; do mkdir -p $DIR; done',
        use_bash_errexit => 1,
        filter_dir => catdir($self->o('output_dir'), 'filters', $self->o('recount3_score_dir')).' '.
                      catdir($self->o('output_dir'), 'filters', $self->o('repeat_overlap_dir')).' '.
                      catdir($self->o('output_dir'), 'filters', $self->o('canonical_splice_site_dir')).' '.
                      catdir($self->o('output_dir'), 'filters', $self->o('pseudogene_overlap_dir')).' '.
                      catdir($self->o('output_dir'), 'filters', $self->o('opp_strand_mismap_dir')).' '.
                      catdir($self->o('output_dir'), 'filters', $self->o('mane_exon_overlap_dir')).' '.
                      catdir($self->o('output_dir'), 'filters', $self->o('terminal_exon_size_dir')).' '.           
                      catdir($self->o('output_dir'), 'filters', $self->o('splice_site_alignment_dir'))      
      },
      -input_ids => [{
        species => $self->o('species'),
      }],
      -flow_into => {
        1 => ['get_chr_list', 'prepare_input_ids'],
      },
      -rc_name => 'default',
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
        '2->A' => ['recount3_score_filter',
                    'repeat_overlap_filter',
                    'canonical_splice_site_filter',
                    'pseudogene_overlap_filter',
                    'opposite_strand_mismapping_filter',
                    'mane_exon_overlap_filter',
                    'terminal_exon_size_filter',
                   # 'splice_site_alignment_filter'
                  ], 
        'A->1' => ['make_filter_list'],
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'recount3_score_filter',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('recount3_score_script')
              .' -gtf '.$self->o('input_file')
              .' -chr #chr#'
              .' -recount3 '.$self->o('intron_score_file')
              .' -cutoff '.$self->o('intron_score_cutoff')
              .' -refhost '.$self->o('ref_db_host')
              .' -refport '.$self->o('ref_db_port')
              .' -refuser '.$self->o('ref_db_user')
              .' -refdbname '.$self->o('ref_db_name')
              .' -out #output_file#',
        filter_dir => catdir($self->o('output_dir'), 'filters', $self->o('recount3_score_dir')),
        output_file => catfile('#filter_dir#', 'recount3_score_#chr#.txt'),
      },
      -rc_name => '5GB',
      -analysis_capacity => 20,
    },

    {
      -logic_name => 'repeat_overlap_filter',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('repeat_overlap_script')
              .' -gtf '.$self->o('input_file')
              .' -chr #chr#'
              .' -repeats '.$self->o('repeat_feature_file')
              .' -refhost '.$self->o('ref_db_host')
              .' -refport '.$self->o('ref_db_port')
              .' -refuser '.$self->o('ref_db_user')
              .' -refdbname '.$self->o('ref_db_name')
              .' -out #output_file#',
        filter_dir => catdir($self->o('output_dir'), 'filters', $self->o('repeat_overlap_dir')),
        output_file => catfile('#filter_dir#', 'repeat_overlap_#chr#.txt'),
      },
      -rc_name => '5GB',
      -analysis_capacity => 20,
    },

    {
      -logic_name => 'canonical_splice_site_filter',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('canonical_splice_site_script')
              .' -gtf '.$self->o('input_file')
              .' -chr #chr#'
              .' -fasta '.$self->o('genome_fasta_file')
              .' -allowed '.$self->o('allowed_splice_site_seq')
              .' -refhost '.$self->o('ref_db_host')
              .' -refport '.$self->o('ref_db_port')
              .' -refuser '.$self->o('ref_db_user')
              .' -refdbname '.$self->o('ref_db_name')
              .' -out #output_file#',
        filter_dir => catdir($self->o('output_dir'), 'filters', $self->o('canonical_splice_site_dir')),
        output_file => catfile('#filter_dir#', 'canonical_splice_sites_#chr#.txt'),
      },
      -rc_name => '5GB',
      -analysis_capacity => 20,
    },

    {
      -logic_name => 'pseudogene_overlap_filter',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('pseudogene_overlap_script')
              .' -gtf '.$self->o('input_file')
              .' -chr #chr#'
              .' -host '.$self->o('ref_db_host')
              .' -port '.$self->o('ref_db_port')
              .' -user '.$self->o('ref_db_user')
              .' -dbname '.$self->o('ref_db_name')
              .' -utr5 180'
              .' -utr3 561'
              .' -out #output_file#',
        filter_dir => catdir($self->o('output_dir'), 'filters', $self->o('pseudogene_overlap_dir')),
        output_file => catfile('#filter_dir#', 'pseudogene_overlap_#chr#.txt'),
      },   
      -rc_name => '5GB',
      -analysis_capacity => 20,
    },

    {
      -logic_name => 'opposite_strand_mismapping_filter',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('opp_strand_mismap_script')
              .' -gtf '.$self->o('input_file')
              .' -chr #chr#'
              .' -host '.$self->o('ref_db_host')
              .' -port '.$self->o('ref_db_port')
              .' -user '.$self->o('ref_db_user')
              .' -dbname '.$self->o('ref_db_name')
              .' -out #output_file#',
        filter_dir => catdir($self->o('output_dir'), 'filters', $self->o('opp_strand_mismap_dir')),
        output_file => catfile('#filter_dir#', 'opp_strand_mismap_#chr#.txt'),
      },
      -rc_name => '5GB',
      -analysis_capacity => 20,
    },

    {
      -logic_name => 'mane_exon_overlap_filter',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('mane_exon_overlap_script')
              .' -gtf '.$self->o('input_file')
              .' -chr #chr#'
              .' -host '.$self->o('ref_db_host')
              .' -port '.$self->o('ref_db_port')
              .' -user '.$self->o('ref_db_user')
              .' -dbname '.$self->o('ref_db_name')
              .' -out #output_file#',
        filter_dir => catdir($self->o('output_dir'), 'filters', $self->o('mane_exon_overlap_dir')),
        output_file => catfile('#filter_dir#', 'mane_exon_overlap_#chr#.txt'),
      },
      -rc_name => '5GB',
      -analysis_capacity => 20,
    },

    {
      -logic_name => 'terminal_exon_size_filter',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('terminal_exon_size_script')
              .' -gtf '.$self->o('input_file')
              .' -chr #chr#'
              .' -minsize '.$self->o('min_terminal_exon_size')
              .' -out #output_file#',
        filter_dir => catdir($self->o('output_dir'), 'filters', $self->o('terminal_exon_size_dir')),
        output_file => catfile('#filter_dir#', 'terminal_exon_size_#chr#.txt'),
      },
      -rc_name => '5GB',
      -analysis_capacity => 20,
    },    

    {
      -logic_name => 'prepare_input_ids',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'cd #filter_dir# && mkdir idfiles && zcat #input_file# | awk \'{print $12}\' | sort -u | perl -pe\'s/\"|;//g\' | shuf - | split -l 1000 -a 4 -d - idfiles/split_',
        input_file => $self->o('input_file'),
        filter_dir => catdir($self->o('output_dir'), 'filters', $self->o('splice_site_alignment_dir')),
      },
      -flow_into => {
        1 => ['create_jobs'],
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'create_jobs',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        inputcmd => 'cd #id_file_dir# && find . -type f -name "split_*" | cut -d"/" -f2',
        column_names => ['iid'],
        id_file_dir => catdir($self->o('output_dir'), 'filters', $self->o('splice_site_alignment_dir'), 'idfiles'),
      },
      -flow_into => {
        '2->A' => ['splice_site_alignment_filter'],
        'A->1' => ['make_filter_list'],
      },
      -rc_name => 'default',
    },    
 
    {
      -logic_name => 'splice_site_alignment_filter',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('splice_site_alignment_script')
              .' -gtf '.$self->o('input_file')
            #  .' -chr #chr#'
              .' -ids #id_file#'
              .' -bamdir '.$self->o('bam_directory')
              .' -fasta '.$self->o('genome_fasta_file')
              .' -readsup '.$self->o('read_support_file')
              .' -window 10'
              .' -same_introns'
              .' -simple'
              .' -out #output_file#',
        filter_dir => catdir($self->o('output_dir'), 'filters', $self->o('splice_site_alignment_dir')),
        id_file_dir => catdir($self->o('output_dir'), 'filters', $self->o('splice_site_alignment_dir'), 'idfiles'),
        id_file => catfile('#id_file_dir#', '#iid#'),
        output_file => catfile('#filter_dir#', 'splice_site_align_'.'#iid#'.'.txt'),
      },
      -rc_name => '10GB',
      -analysis_capacity => 25,
    },

    {
      -logic_name => 'make_filter_list',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'cd #filter_dir#;
                cd '.$self->o('recount3_score_dir').';
                cat recount3_score_*.txt | grep -P "\tno\t" | cut -f1 > failed_recount3_score.txt;
                cd ../'.$self->o('repeat_overlap_dir').';
                cat repeat_overlap_*.txt | cut -f1 > failed_repeat_overlap.txt;
                cd ../'.$self->o('canonical_splice_site_dir').';
                cat canonical_splice_sites_*.txt | grep -P "\tno\t" | cut -f1 > failed_canonical_splice_sites.txt;
                cd ../'.$self->o('pseudogene_overlap_dir').';
                cat pseudogene_overlap_*.txt | cut -f1 > failed_pseudogene_overlap.txt;
                cd ../'.$self->o('opp_strand_mismap_dir').';
                cat opp_strand_mismap_*.txt | cut -f1 > failed_opposite_strand_mismapping.txt;
                cd ../'.$self->o('mane_exon_overlap_dir').';
                cat mane_exon_overlap_*.txt | grep -v ^# | cut -f1 > failed_mane_exon_overlap.txt;
                cd ../'.$self->o('terminal_exon_size_dir').';
                cat terminal_exon_size_*.txt | grep -P "\tno\t" | cut -f1 > failed_terminal_exon_size.txt;                
                cd ../'.$self->o('splice_site_alignment_dir').';
                cat splice_site_align_*.txt | awk \'$2=="no" && $3>0\' | cut -f1 > failed_splice_site_alignment.txt;
                cd ..;
                find . -name "failed_*.txt" | xargs cat | sort -u > #output_file#',
        filter_dir => catdir($self->o('output_dir'), 'filters'),
        output_file => $self->o('filter_file'),
      },
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
    'default' => {SLURM => '--time=1:00:00 --mem=1000', LSF => '-q short -M 1000 -R"select[mem>1000] rusage[mem=1000]"'},   
    '5GB' => {SLURM => '--time=2:00:00 --mem=5000', LSF => '-q short -M 5000 -R"select[mem>5000] rusage[mem=5000]"'},
    '10GB' => {SLURM => '--time=4:00:00 --mem=10000', LSF => '-q short -M 10000 -R"select[mem>10000] rusage[mem=10000]"'},    
    '10GB_long' => {SLURM => '--time=3-00:00 --mem=10000', LSF => '-q long -M 10000 -R"select[mem>10000] rusage[mem=10000]"'},
  };
}

1;

