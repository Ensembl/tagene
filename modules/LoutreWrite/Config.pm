
package LoutreWrite::Config;

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

use base 'Exporter';
our @EXPORT = qw(%DBA);

our %DBA;

my $loutre_db = get_loutre_db_adaptor();
$DBA{'loutre'} = $loutre_db;
my $pipe_db = get_pipe_db_adaptor();
$DBA{'pipe'} = $pipe_db;
my $core_db = get_core_db_adaptor();
$DBA{'core'} = $core_db;
my $intron_db = get_intropolis_db_adaptor();
$DBA{'intron'} = $intron_db;
my $polyAseq_db = get_polyAseq_db_adaptor();
$DBA{'polyAseq'} = $polyAseq_db;



#Connect to the loutre database (the live/production one, not necessarily the otter database that is being edited)
sub get_loutre_db_adaptor {
  my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
     -host   => 'mysql-ens-havana-prod-1',
     -port   => 4581,
     -user   => 'ensro',
     -pass   => undef,
     -dbname => 'loutre_human',
     -driver => 'mysql',
  );
  $db->dbc->reconnect_when_lost(1);
  return $db;
}


#Connect to the Otter pipeline database
sub get_pipe_db_adaptor { 
  my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host   => 'mysql-ens-havana-prod-1',
    -port   => 4581,
    -user   => 'ensro',
    -pass   => undef,
    -dbname => 'pipe_human',
    -driver => 'mysql',
  );
  $db->dbc->reconnect_when_lost(1);
  return $db;
}


#Connect to the Ensembl core database
sub get_core_db_adaptor {
  my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
     -host   => 'mysql-ens-havana-prod-1',
     -port   => 4581,
     -user   => 'ensro',
     -pass   => undef,
     -dbname => 'homo_sapiens_core_101_38',
     -driver => 'mysql',
  );
  $db->dbc->reconnect_when_lost(1);
  #DNA db
  my $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
     -host   => 'mysql-ens-havana-prod-1',
     -port   => 4581,
     -user   => 'ensro',
     -pass   => undef,
     -dbname => 'homo_sapiens_core_101_38',
  );
  $db->dnadb($dnadb);
  return $db;
}


#Connect to the Intropolis database
sub get_intropolis_db_adaptor { 
  my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host   => 'mysql-ens-havana-prod-1',
    -port   => 4581,
    -user   => 'ensro',
    -pass   => undef,
    -dbname => 'gencode_sf5_human_introns',
    -driver => 'mysql',
  );
  $db->dbc->reconnect_when_lost(1);
  return $db;
}


#Connect to the (Merck) polyA-seq database
sub get_polyAseq_db_adaptor { 
  my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host   => 'mysql-ens-havana-prod-1',
    -port   => 4581,
    -user   => 'ensro',
    -pass   => undef,
    -dbname => 'gencode_polyAseq',
    -driver => 'mysql',
  );
  $db->dbc->reconnect_when_lost(1);
  return $db;
}





1;


