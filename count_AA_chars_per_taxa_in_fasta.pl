#!/usr/bin/perl                                                                                                                                                     

use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;

##########################################################################################                                                                          
#                                                                                                                                                                   
# #################### Iker Irisarri. May 2018. Uppsala University ##################### #                                                                          
#                                                                                                                                                                   
# It counts the number AA chars in a fasta file                                                                                                               
#                                                                                                                                                                   
##########################################################################################                                                                          


my $usage = "perl count_AA_chars_per_taxa_in_fasta.pl *.fasta\n";
my @infiles = @ARGV or die $usage;


my %results;

foreach my $infile ( @infiles ) {

    my $seqio_obj = Bio::SeqIO->new('-file' => "<$infile",
                                    '-format' => "fasta",
        );

    while (my $inseq = $seqio_obj->next_seq) {

		my $aa = "0";                                                                                                                                           
		my $other = "0";                                                                                                                                          
		my $total;
		my %others = (); # will save characters that are counted as not DNA

        my $header = $inseq->primary_id;
        my $sequence = $inseq->seq;
        my @sequences = split ("", $sequence);

        foreach my $elem ( @sequences ) {

            if ( $elem eq "A" || $elem eq "R" || $elem eq "N" || $elem eq "D" ||
            	 $elem eq "C" || $elem eq "Q" || $elem eq "E" || $elem eq "G" ||
            	 $elem eq "H" || $elem eq "I" || $elem eq "L" || $elem eq "K" ||
            	 $elem eq "M" || $elem eq "F" || $elem eq "P" || $elem eq "S" ||
            	 $elem eq "T" || $elem eq "W" || $elem eq "Y" || $elem eq "V" ) {

                $aa++;
                $total++;
            }
            else {

                $other++;
                
                if ( !exists $others{$elem} ) {
                	$others{$elem} = 1;
	                $total++;
                }
                else {
	                $others{$elem}++;
	                $total++;
                }
            }
        }
        my $others_ref = \%others;
        $results{$header} = [$aa, $other, $total, $others_ref];
    }
}

# print out results                                                                                                                                                 

foreach my $key ( sort keys %results ) {

	print "$key\tDNA: $results{$key}[0] OTHER: $results{$key}[1] TOTAL: $results{$key}[2]\n";

=pod
	foreach my $key ( sort keys %others ) {

		print "\t$key: $others{$key}\n";
	}
=cut
}

print "\ndone!\n\n";

