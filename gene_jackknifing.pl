#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

# Iker Irisarri, January 2015, University of Konstanz
# This script takes a file with gene file names and subsamples it without replacement (gene jackknifing)
# the number (-n) and size (-s) of the replicates should be given in the command line (integers)


my $usage = "gene_jackknifing.pl -f list_of_gene_files -s size_of_replicate -n number of replicates > out\n";
my $genes_in; 
my $replicate_size;
my $num_replicates;

GetOptions (
    "file=s" => \$genes_in,  # string
    "size=i" => \$replicate_size,  # integer
    "number=i" => \$num_replicates
    ) or die ("$usage\n");

my @genes_all;
my @subsample;

# read-in gene file names from file (-f)
open (IN, "<", $genes_in);

while ( my $line =<IN> ) {
    chomp $line;
    push (@genes_all, $line);
}



for (my $i=0; $i<$num_replicates; $i++) {

    # start with complete list of gene files every replicate
    # and empty the list of subsampled genes
    my @genes = @genes_all;
    my @replicate = ();

    for (my $j=0; $j<$replicate_size; $j++) {

	#get one element at random from array
	my $random_elem = $genes[rand @genes];

	#save it to new array
	push (@replicate, $random_elem);

	#remove that element from origina array (to avoid sampling this element again)
	#jackknife is sampling without replacement
	@genes = grep { $_ ne $random_elem } @genes;

	#print "$random_elem\n";
	#print join ('--', @genes), "\n";

    }

    print "replicate number $i:\n";
    print join ("\n", @replicate), "\n\n";

}
