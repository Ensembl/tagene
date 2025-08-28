
package LoutreWrite::Config;

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;

use base 'Exporter';
our @EXPORT = qw( %DBA $SPECIES );

our $SPECIES;
our %DBA;

if (defined($SPECIES)){
  get_db_adaptors();
}


sub get_db_adaptors_from_registry {
  my $registry_file = shift;

  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_all($registry_file);
  $registry->set_reconnect_when_lost();

  $DBA{'havana'} = $registry->get_DBAdaptor('havana_'.$SPECIES, 'core');
  unless ($DBA{'havana'}){
    $DBA{'havana'} = $DBA{'otter'};
  }
  $DBA{'core'} = $registry->get_DBAdaptor('core_'.$SPECIES, 'core');
  $DBA{'intron'} = $registry->get_DBAdaptor('intron_'.$SPECIES, 'core');
  $DBA{'polyAseq'} = $registry->get_DBAdaptor('polyAseq_'.$SPECIES, 'core');
  $DBA{'polyAdb'} = $registry->get_DBAdaptor('polyAdb_'.$SPECIES, 'core');  
}
  

sub get_db_adaptors{
  $DBA{'havana'} = get_havana_db_adaptor();
  unless ($DBA{'havana'}){
    $DBA{'havana'} = $DBA{'otter'};
  }
  $DBA{'core'} = get_core_db_adaptor();
  $DBA{'intron'} = get_intron_db_adaptor();
  $DBA{'polyAseq'} = get_polyAseq_db_adaptor();
}


#Connect to the havana database (the live/production one, not necessarily the otter database that is being edited)
sub get_havana_db_adaptor {
  my $host = '';
  my $port = '';
  my $user = '';
  my $dbname = "havana_".$SPECIES;
  if ($SPECIES eq "human"){
    $host = '';
    $port = '';
    $dbname = "havana_".$SPECIES."_backup";
    #$dbname = "havana_".$SPECIES."_tagene_test_1"; #temporarily to avoid failures when the backup database is being updated
  }
  elsif ($SPECIES eq "rat"){
    $dbname = "havana_rattus_norvegicus";
  }
  my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
     -host   => $host,
     -port   => $port,
     -user   => $user,
     -pass   => undef,
     -dbname => $dbname,
     -driver => 'mysql',
  );
  $db->dbc->reconnect_when_lost(1);
  return $db;
}


#Connect to the Ensembl core database
sub get_core_db_adaptor {
  my $host = '';
  my $port = '';
  my $user = '';
  my $dbname;
  if ($SPECIES eq "human"){
    $dbname = "homo_sapiens_core_112_38";
  }
  elsif ($SPECIES eq "mouse"){
    $dbname = "mus_musculus_core_112_39";
  }
  elsif ($SPECIES eq "rat"){
    $host = '';
    $port = '';
    $user = '';
    $dbname = "rattus_norvegicus_core_111_72";
  } 
  my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
     -host   => $host,
     -port   => $port,
     -user   => $user,
     -pass   => undef,
     -dbname => $dbname,
     -driver => 'mysql',
  );
  $db->dbc->reconnect_when_lost(1);
  #DNA db
  #my $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
     #-host   => $host,
     #-port   => $port,
     #-user   => $user,
     #-pass   => undef,
     #-dbname => $dbname,
  #);
  #$db->dnadb($dnadb);
  return $db;
}


#Connect to the intron database
sub get_intron_db_adaptor {
  my $host = '';
  my $port = '';
  my $user = '';
  my $dbname;
  if ($SPECIES eq "human"){
    $dbname = "gencode_snaptron";
  }
  elsif ($SPECIES eq "mouse"){
    $dbname = "gencode_mouse_recount3";
  }
  elsif ($SPECIES eq "rat"){
    $host = '';
    $port = '';
    $user = '';
    $dbname = "rattus_norvegicus_core_111_72";
  } 
  my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
     -host   => $host,
     -port   => $port,
     -user   => $user,
     -pass   => undef,
     -dbname => $dbname,
     -driver => 'mysql',
  );
  $db->dbc->reconnect_when_lost(1);
  return $db;
}


#Connect to the (Merck) polyA-seq database
sub get_polyAseq_db_adaptor {
  my $host = '';
  my $port = '';
  my $user = '';
  my $dbname;
  if ($SPECIES eq "human"){
    $dbname = "gencode_polyAseq";
  }
  elsif ($SPECIES eq "mouse"){
    $dbname = "gencode_mouse_polyAseq_mm9_assembly_schema_73";
  }
  elsif ($SPECIES eq "rat"){
    $host = '';
    $port = '';
    $user = '';
    $dbname = "rattus_norvegicus_core_111_72";
  } 
  my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
     -host   => $host,
     -port   => $port,
     -user   => $user,
     -pass   => undef,
     -dbname => $dbname,
     -driver => 'mysql',
  );
  $db->dbc->reconnect_when_lost(1);
  return $db;
}


1;


