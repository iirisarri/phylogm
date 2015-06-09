#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;


# Iker Irisarri. University of Konstanz, June 2015
# random_sample_aln_columns.pl
# Subsamples alignment (fasta format) by randomly picking N columns (without replacement)


my $usage = "random_sample_aln_columns.pl in.fasta num_final_cols(int) > output.fa\n";
my $in_fasta = $ARGV[0] or die $usage;
my $num_out_pos = $ARGV[1] or die $usage;


## Store alignment in hash of arrays
open(IN, "<", $in_fasta);

my %fasta_hash;
my $header;
my @seq;

while ( my $line =<IN> ) {
	chomp $line;
	if ($line =~ /(>.+)/ ) {
		$header = $1;
	}
	else {
		@seq = split (//, $line);
	} 
	$fasta_hash{$header} = [@seq];
}

#print Dumper \%fasta_hash;


## Get N random numbers

# get total number of positions in aln & check if seqs are aligned
my $in_aln_total = 0;

foreach my $sp ( keys %fasta_hash ) {
	#print "$sp: ", scalar @{ $fasta_hash{$sp} }, "\n";
	if ( $in_aln_total != 0 && $in_aln_total != scalar @{ $fasta_hash{$sp} } ) {
		print STDERR "WARN: sequence lengths are unequal.
		 They might not be aligned: check your input!\n";
	}
	else {
		$in_aln_total = scalar @{ $fasta_hash{$sp} };
	}
}

## Generate N random numbers (they must be different!)

#my $num_out_pos = 1;
my %rand_numbers;
my $pos = 0;

while ( $pos < $num_out_pos ) {

	my $num = int(rand($in_aln_total));
	
	# avoid numbers already in %rand_numbers (aln cols are sampled without replacement)
	if ( !exists $rand_numbers{$num} ) {
		$pos++;
		$rand_numbers{$pos} = $num;
	}
}

#print Dumper \%rand_numbers;
	
## Write out aln consisting of N number of columns

# numbers from %rand_numbers are used to extract a given aln column
# by using the index of the @seq array in which sequences are stored

foreach my $key ( sort keys %fasta_hash ) {
	print "$key\n";
	foreach my $n ( keys %rand_numbers ) { 
		print @{ $fasta_hash{$key} }[$rand_numbers{$n}];
	}
	print "\n";
}

