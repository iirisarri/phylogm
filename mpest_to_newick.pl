#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;


## Iker Irisarri, University of Konstanz. Jun 2015 ##

# mpest_to_newick.pl
# simple script to parse MP-EST results and print final tree in newick format
# MP-EST produces a nexus file where the last tree is the MP-EST tree
# Taxa should also be translated using nexus translate block

my $usage = "mpest_to_newick.pl mpest_result.nex.tre > STDOUT\n";
my $intree = $ARGV[0] = shift or die "$usage";

my %hash;
my $mpest_tree;

open (IN, "<", $intree);

while ( my $line =<IN> ) {

	chomp $line;
	# parse taxa from translate block (line start with \t)
	if ( $line =~ /^\t/ ) {
		
		$line =~ /^\t(\d+)\s(\w+)/;
		my $num = $1;
		my $taxa =$2;
		
		$hash{$num} = $taxa;	
	}
	# parse final mpest tree
	if ( $line =~ /^\s\stree mpest/ ) {
		
		my @lines = split (" ", $line);
		
		$mpest_tree = $lines[4];
		last;
	}
}

# substitute numbers by species using info in %hash

foreach my $key ( keys %hash ) {

	# to make sure only the complete numbers are replaced,
	# anchor them with "(" or "," and ":"
	$mpest_tree =~ s/\($key:/\($hash{$key}:/;
	$mpest_tree =~ s/,$key:/\,$hash{$key}:/;
	#print "$mpest_tree\n";
}

print "$mpest_tree\n";

print STDERR "done!\n";