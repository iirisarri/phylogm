#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

# Iker Irisarri, January 2015, University of Konstanz
# This script takes a file with gene file names and subsamples it without replacement (gene jackknifing)
# the number (-n) and size (-s) of the replicates should be given in the command line (integers)
# bash commands to sort files into new folders, ready for concatenation

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

    #add folder name to each array element and sort
    foreach my $elem ( @replicate ) {
	$elem = $elem . " replicate.\$i";
    }

    my @replicate_sorted = sort @replicate;

    print "i=$i\n\n";
    print "mkdir replicate.\$i\n\ncp ../";
    print join ("\ncp ../", @replicate_sorted), "\n\n";

}


__END__

print statements will print something like:

"i=0

 mkdir replicate.$i
 cp ../G15821_all.aln.trim.fas replicate.$i
 cp ../G15823_all.aln.trim.fas replicate.$i
 cp ../G15837_all.aln.trim2.fas replicate.$i
 cp ../G15845_all.aln.trim2.fas replicate.$i
 cp ../G15851_all.aln.trim.fas replicate.$i
 cp ../G15858_all.aln.trim.fas replicate.$i
 cp ../G15867_all.aln.trim.fas replicate.$i
 cp ../G15869_all.aln.trim.fas replicate.$i


i=1

 mkdir replicate.$i
 cp ../G15807_all.aln.trim.fas replicate.$i
 cp ../G15820_all.aln.trim.fas replicate.$i
 cp ../G15823_all.aln.trim.fas replicate.$i
 cp ../G15832_all.aln.trim.fas replicate.$i
 cp ../G15837_all.aln.trim2.fas replicate.$i
 cp ../G15844_all.aln.trim.fas replicate.$i
 cp ../G15845_all.aln.trim2.fas replicate.$i

i=2

etc"
