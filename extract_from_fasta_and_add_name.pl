#!/usr/bin/perl -w

use strict;

use Bio::SeqIO;
use Data::Dumper;

# Iker Irisarri, Univeristy of Konstanz, January 2015
# modification of extract_from_fasta_by_name.pl to add gene names to fasta headers

# example of usage:
# files=*.fa
# for f in $files; do perl /computation/iker/software/scripts/phylogm/extract_from_fasta_and_add_name.pl $f lati >> all_lati_queries.fa; done

my $usage = "extractSeqFromFasta fasta_file query_file\n";
my $fasta = $ARGV[0] or die $usage;
my $query = $ARGV[1] or die $usage;

# read fasta file with SeqIO
my $seqio_obj = Bio::SeqIO->new('-file' => "<$fasta",
                		        '-format' => "fasta");

# get gene name
$fasta =~ /(\w+)\.fa\w*/;
my $name = $1;
                		        
# declare hash
my %queries;

#open query file, chomp line and save it into the array with push
#the hash will contain gene names (variable $line) as keys and 1 as a value
open (IN , $query) or die "Can't open $query, $!\n";

while (my $line = <IN>){
	chomp $line;
	$queries{$line} = 1;	
	}

#check structure of hash	
#print Dumper \%queries, "\n";

while (my $seq_obj = $seqio_obj->next_seq){

    my $seqname = $seq_obj->primary_id;
    my $description = $seq_obj->description;

	if (exists($queries{$seqname})) {

                print ">",  $seq_obj->primary_id, "_$name\n";
		print $seq_obj->seq, "\n";
	}
}



__END__

Structure of %queries

$VAR1 = {
          'comp161255_c2_seq1' => 1,
          'comp177440_c0_seq1' => 1,
          'comp178046_c1_seq9' => 1,
          'comp170200_c0_seq5' => 1,
          'comp166596_c2_seq13' => 1,
          'comp149867_c0_seq1' => 1,
          'comp174145_c0_seq7' => 1,
          'comp165635_c0_seq1' => 1,
          'comp167515_c0_seq5' => 1,
				[...]
        };
$VAR2 = '
';
