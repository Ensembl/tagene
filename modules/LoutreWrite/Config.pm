
package LoutreWrite::Config;

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

use base 'Exporter';
our @EXPORT = qw( %DBA $SPECIES );

our $SPECIES;
our %DBA;

if (defined($SPECIES)){
  get_db_adaptadors();
}


sub get_db_adaptadors{
  $DBA{'havana'} = get_havana_db_adaptor();
  $DBA{'pipe'} =  get_pipe_db_adaptor();
  $DBA{'core'} = get_core_db_adaptor();
  $DBA{'intron'} = get_intron_db_adaptor();
  $DBA{'polyAseq'} = get_polyAseq_db_adaptor();
}


#Connect to the havana database (the live/production one, not necessarily the otter database that is being edited)
sub get_havana_db_adaptor {
  my $host = 'mysql-ens-havana-prod-1';
  my $port = 4581;
  my $user = 'ensro';
  my $dbname = "havana_".$SPECIES;
  if ($SPECIES eq "human"){
    $host = 'mysql-ens-havana-prod-2';
    $port = 4682;
    #$dbname = "havana_".$SPECIES."_backup";
    $dbname = "havana_".$SPECIES."_tagene_test_1"; #temporarily to avoid failures when the backup database is being updated
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


#Connect to the Otter pipeline database
sub get_pipe_db_adaptor {
  my $host = 'mysql-ens-havana-prod-1';
  my $port = 4581;
  my $user = 'ensro';
  my $dbname = "pipe_".$SPECIES;
  if ($SPECIES eq "rat"){
    $dbname = "pipe_rattus_norvegicus";
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
  my $host = 'mysql-ens-havana-prod-2';
  my $port = 4682;
  my $user = 'ensro';
  my $dbname;
  if ($SPECIES eq "human"){
    $dbname = "homo_sapiens_core_112_38";
  }
  elsif ($SPECIES eq "mouse"){
    $dbname = "mus_musculus_core_112_39";
  }
  elsif ($SPECIES eq "rat"){
    $host = 'mysql-ens-mirror-1';
    $port = 4240;
    $user = 'ensro';
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
  my $host = 'mysql-ens-havana-prod-2';
  my $port = 4682;
  my $user = 'ensro';
  my $dbname;
  if ($SPECIES eq "human"){
    $dbname = "gencode_snaptron";
  }
  elsif ($SPECIES eq "mouse"){
    $dbname = "gencode_mouse_recount3";
  }
  elsif ($SPECIES eq "rat"){
    $host = 'mysql-ens-mirror-1';
    $port = 4240;
    $user = 'ensro';
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
  my $host = 'mysql-ens-havana-prod-2';
  my $port = 4682;
  my $user = 'ensro';
  my $dbname;
  if ($SPECIES eq "human"){
    $dbname = "gencode_polyAseq";
  }
  elsif ($SPECIES eq "mouse"){
    $dbname = "gencode_mouse_polyAseq_mm9_assembly_schema_73";
  }
  elsif ($SPECIES eq "rat"){
    $host = 'mysql-ens-mirror-1';
    $port = 4240;
    $user = 'ensro';
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


