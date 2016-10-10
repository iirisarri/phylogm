#!/usr/bin/perl
 
use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Data::Dumper;

#######################################################
#
# gene_jackknife_from_fasta_collection.pl
# merging of concat_fasta.pl and gene_jackknifing.pl
#
# Iker Irisarri. Jul 2016.
#
#######################################################

# This script will generate gene jackknife replicates by concatenating genes 
# sampled without replacement until the specified replicate length is reached.
#
# Infiles are taken from current directory, reading only those matching the provided file extension.
# WATCH OUT: running the script multiple times on the same folder can produce 
# incorporation of already generated replicates as single genes if same extension is used.
#
# Requirements: all taxa in individual fasta files should have matching headers.
# Besides replicates, a summary file is printed out including basic information.

my $usage = "gene_jackknife_from_fasta_collection.pl input_file_extension num_replicates num_final_positions outfile_pattern\ne.g. gene_jackknife_from_fasta_collection.pl fasta 10 10000 myreplicate\n";
my $input_file_extension = $ARGV[0] or die $usage;
my $num_replicates = $ARGV[1] or die $usage;
my $replicate_size = $ARGV[2] or die $usage;
my $outfile_pattern = $ARGV[3] or die $usage;

# open summary file and print out stuff
my $summary = "_" . $outfile_pattern . ".summary.txt";
open (SUM, ">", $summary);
print SUM "Num. replicates = $num_replicates\n";
print SUM "Min. replicate size = $replicate_size\n";
print SUM "Input files = \.$input_file_extension\n";
print SUM "Outfile name = $outfile_pattern\n\n";

# read input files from current directory
# only files with the provided extension will be read
opendir(DIR, ".");
my @genes_all = grep(/\.$input_file_extension$/,readdir(DIR));
closedir(DIR);


# loop through jackknife replicates
for (my $i=0; $i<$num_replicates; $i++) {

	# print out to summary file
	print SUM "replicate: $i\n";

    # start with complete list of gene files every replicate
    # and empty the list of subsampled genes
    my @genes = @genes_all;
    my @replicate = ();

	# declare hashes (new every jackknife replicate)
	my %new = ();
	my %concat = ();
	# measure total length of concatenation
	my $total_length = 0;

	# loop through genes until $replicate_size is reached
	while ( $total_length < $replicate_size ) {

		# print out error when gene alignments are finished
		if ( scalar @genes == 0 ) {
	
			print STDERR "ERR: not enough genes to reach specified replicate length $replicate_size!\n";
		}
		# choose random gene from @genes
		my $gene = $genes[rand @genes];

		my %new = ();
		my $new_seq_length = 0;

		# read file in with seqio
		my $seqio_obj = Bio::SeqIO->new(
			'-file' => "<$gene", 
			'-format' => "fasta", 
		);
					
		# store sequences of new file into hash %new
		while (my $inseq = $seqio_obj->next_seq) {

			# scape empty sequences in present
			next if ( $inseq->seq =~ /^[-XNn?]*$/ );

			$new{$inseq->primary_id}=$inseq->seq;
			$new_seq_length = $inseq->length;
		}

		# add sequences to hash %concat
		foreach my $k ( keys %new ) {

			# create new taxa if not existing in %concat
			if ( !exists $concat{$k} ) {

			# prior to concatenation, fill previous positions with "?s" for $total_length
			# for gene #1 this $total_concat_length == 0
			my $seq1 = '?' x $total_length . $new{$k};

			$concat{$k} = $seq1;

			}
			else {
	
				# concatenate
				my $seq2 = $concat{$k} . $new{$k};
				# reassign
				$concat{$k} = $seq2;
			}
		}
		
		# fill with ?s any taxa present in %concat but not in %new
		foreach my $j ( keys %concat ) {
	
			if ( !exists $new{$j} ) {
		
				my $seq3 = $concat{$j} . '?' x $new_seq_length;
				$concat{$j} = $seq3;
			}
		}

		# add sequence length of new alignment
		$total_length += $new_seq_length;
	
		# remove gene from @genes_all
		@genes = grep { $_ ne $gene } @genes;

		# print out to summary file
		print SUM "\t$gene $new_seq_length\n";

	}
	# check that replicates reached the min number of positions
	if ( $total_length < $replicate_size && scalar @genes == 0 ) {
	
		print STDERR "ERR: concatenation of all genes does not reach the 
						specified replicate size of $replicate_size!\n";
	}
	# print out concatenated sequence to file
	my $outfile = $outfile_pattern . "_rep${i}.fa";

	open (OUT, ">", $outfile);
	
	foreach my $taxa ( sort keys %concat ) {

		print OUT ">$taxa\n";
		print OUT $concat{$taxa}, "\n";
	}
	close(OUT);
	print SUM "\ttotal length: $total_length\n\n";

	# got to next replicate
	next;
}
close(SUM);

print STDERR "\nNum. replicates generated = $num_replicates\n";
print STDERR "Min. number of positions per replicate = $replicate_size\n";
print STDERR "\ndone!\n\n";

