#!/usr/bin/perl -w

#This script extracts the read source information for each PASA assembled transcript in a GTF file
#Usage: perl extract_pasa_assembly_info.pl \
#                              -gtf my_gene_data.gtf \
#                              -pasa gencode_pasa_Captureseq_2.pasa_assemblies_described.txt.gz \
#                              -lib hsAll_Cap1_Undeter_primAlign.gff3.gz,hsAll_Cap1_Testes_primAlign.gff3.gz,hsAll_Cap1_Brain_primAlign.gff3.gz

use strict;
use warnings;
use Getopt::Long;
$| = 1;

my $gtf_file;
my $pasa_file;
my @lib_files;
my $outfile;

&GetOptions(
			'gtf=s'      => \$gtf_file,
			'pasa=s'     => \$pasa_file,
			'lib=s'      => \@lib_files,
			'outfile=s'  => \$outfile,		
			);


#Find original library/column for each read
@lib_files = split(/,\s*/, join(',', @lib_files));
my %read_libs;
foreach my $lib_file (@lib_files){
	print "Reading $lib_file ...\n";
	my $lib_name;
	if ($lib_file =~ /hsAll_Cap1_(\w+)_primAlign/){
		$lib_name = "pacbio_capture_seq_hs".$1;
	}
	elsif ($lib_file =~ /([A-Za-z0-9_-]+)_comb_withLaneID/){ #SLRseq library names are included in the read names, anyway ...
		$lib_name = "SLR-seq_".$1;
		$lib_name =~ s/-DNA_A01-LRAAD-01//;
	}
	else{
	  ($lib_name) = $lib_file =~ /\/(\w+)\.gff3(\.gz)?/;
	}
	
	if ($lib_file =~ /\.gz/){
		open (LIB, "zcat $lib_file | ") or die "Can't open file $lib_file: $!";
	}
	else{
		open (LIB, $lib_file) or die "Can't open file $lib_file: $!";
	}
	while (<LIB>){
		#chr17   gff3_for_pasa   cDNA_match      16381902        16382740        98      -       .       ID=m150425_202400_42137_c100813572550000001823175210081502_s1_p0_96135_ccs;Name=m150425_202400_42137_c100813572550000001823175210081502_s1_p0_96135_ccs;Target=m150425_202400_42137_c100813572550000001823175210081502_s1_p0_96135_ccs 1 839;Gap=97S37M2D49M717N60M1D53M2D93M1D78M1D52M1D147M1D26M1D165M1D56M1D17M1D81M123S
		my $read_name;
		my @fs = split(/;/);
		foreach my $f (@fs){
			if ($f =~ /Name=(.+)/){ #Eg. align_id:1326054|asmbl_45957
				$read_name = $1;
				last;
			}
		}
		if ($read_name){
			if ($read_libs{$read_name} and $read_libs{$read_name} ne $lib_name){
				warn "$read_name found in ".$read_libs{$read_name}." and $lib_name\n";
			}
			$read_libs{$read_name} = $lib_name;
		}	
	}
	close (LIB);
}



#Find read names for each transcript (PASA assembly)
my %assemblies;
if ($pasa_file =~ /\.gz/){
	open (PASA, "zcat $pasa_file | ") or die "Can't open file $pasa_file: $!";
}
else{
	open (PASA, $pasa_file) or die "Can't open file $pasa_file: $!";
}
while (<PASA>){
	#chr1    3       asmbl_22        m150415_052801_42137_c100793412550000001823172010081557_s1_p0_9754_ccs,m150415_115530_42137_c100793272550000001823172010081550_s1_p0_116644_ccs,m150415_203101_42137_c100793272550000001823172010081552_s1_p0_99804_ccs,m150418_110041_42137_c100813522550000001823175210081553_s1_p0_27937_ccs,m150420_082649_42137_c100793212550000001823172010081515_s1_p0_123713_ccs,m150420_170847_42137_c100793212550000001823172010081517_s1_p0_77763_ccs,m150424_115814_42137_c100813372550000001823175210081563_s1_p0_145508_ccs,m150426_134437_42137_c100813572550000001823175210081506_s1_p0_90639_ccs,m150426_180445_42137_c100813572550000001823175210081507_s1_p0_73799_ccs       orient(a-/s-) align: 266119(2441)-268204(356)>CT....AC<268667(355)-268816(206)>CT....AC<289266(205)-289370(101)>CT....AC<710635(100)-710734(1)
	my @fs = split(/\s+/);
	my ($assembly_id, $read_names) = ($fs[2], $fs[3]);
	$assemblies{$assembly_id} = $read_names;
}
close (PASA);


open (OUT, ">$outfile") or die "Can't open file $outfile: $!";

#Read GTF file - extract transcript ids
my %seen;
if ($gtf_file =~ /\.gz/){
	open (GTF, "zcat $gtf_file | ") or die "Can't open file $gtf_file: $!";
}
else{
	open (GTF, $gtf_file) or die "Can't open file $gtf_file: $!";
}
while (<GTF>){
	#19      pacbio  transcript      3216669 3218031 .       -       .       gene_id "PASA_cluster_10799"; transcript_id "align_id:1326054|asmbl_45957"; gene_type "missing_biotype"; gene_status "UNKNOWN"; gene_name "PASA_cluster_10799"; transcript_type "missing_biotype"; transcript_status "UNKNOWN"; transcript_name "align_id:1326054|asmbl_45957"; level 2;
	next if /^#/;
	my @fs = split(/;/);
	my $tid;
	foreach my $f (@fs){
		if ($f =~ /transcript_id \"(.+)\"/){ #Eg. align_id:1326054|asmbl_45957
			$tid = $1;
			last;
		}
	}
	next unless $tid;
	if (!$seen{$tid}){
		if ($tid =~ /(asmbl_\d+)/){
			my $assembly_id = $1;
			if ($assemblies{$assembly_id}){
				my %read_columns; #Otter columns
				my @reads = split(/,/, $assemblies{$assembly_id});
				foreach my $read (@reads){
					if (my $colname = $read_libs{$read}){
						$read_columns{$colname}{$read} = 1;
					}
				}
				print OUT $tid."\t";
				foreach my $colname (keys %read_columns){
					print OUT $colname." : ".join(", ", keys %{$read_columns{$colname}})." ; ";
				}
				print OUT "\n";
			}
		}
		$seen{$tid} = 1;
	}
}
close (GTF);

close (OUT);




