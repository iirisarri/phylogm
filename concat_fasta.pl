#!/usr/bin/perl
 
use strict;
use warnings;
use Bio::SeqIO;
use Data::Dumper;

# concat_fasta.pl
# Iker Irisarri. Jul 2016.
# This script will concatenate the files in the order given in the command line
# requirements: all taxa in individual fasta files should have matching headers

my $usage = "concat_fasta.pl infile1.fa infile2. etc. > outfile.fa";
my @infiles = @ARGV or die $usage; 

# declare hashes
my %new;
my %concat;

# loop through the list of infiles
foreach my $infile (@infiles) {


	# read file in with seqio
	my $seqio_obj = Bio::SeqIO->new('-file' => "<$infile", 
									'-format' => "fasta", 
					);

	# store sequences of new file into hash %new
	while (my $inseq = $seqio_obj->next_seq) {

		$new{$inseq->primary_id}=$inseq->seq;
	}

	# add sequences to hash %concat
	foreach my $key ( keys %new ) {

		# create new taxa if not existing in %concat
		if ( !exists $concat{$key} ) {
	
			$concat{$key} = $new{$key};
		}
		else {
	
			# concatenate
			my $seq2 = $concat{$key} . $new{$key};
			# reassign
			$concat{$key} = $seq2;
		}
	}
}

# print out concatenated sequence
foreach my $taxa ( sort keys %concat ) {

	print ">$taxa\n";
	print $concat{$taxa}, "\n";
}

print STDERR "\ndone!\n\n";