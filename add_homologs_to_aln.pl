#!/usr/bin/perl

##########################################################################################
#
# #################### Iker Irisarri. Feb 2018. Uppsala University ##################### #
#
# Pipeline that identifies homologs from database and adds them to alignment
#
# 1. Read in input alignment (protein)
# 2. Blast all sequences in alignment against protein dataset (eg. TransDecoder-predicted 
#		ORFs from transcriptome)
# 3. parse BLASTP report
# 4. Extract N best BLAST hits to file
#
# A few settings can be modified (see below usage)
#
# Fasta parsing in subroutine
#	-Alignment: it will de-align it to *.nogaps for BLAST
#	-Database (transcriptome): it will make headers unique (eg. TransDecoder can predict 
#	multiple ORFs per transcript). It will append seq length and "***" to headers
#
# In addition to main output file, the script will create intermediate files that will be
#	retained for later use in "for" loops (*.nogaps, *.flags,blast DB to be removed manually)
#
# Typical use in for loops:
#	Single transcriptome to multiple alignments:
#		for f in *alns; do perl add_homologs_to_aln.pl $f transcriptome.fa; done	
#	Multiple transcriptomes to multiple alignments:
#	for f in *alns; do for d in transcriptome.*fa; do perl add_homologs_to_aln.pl $f $d; done	
#
##########################################################################################

use strict;
use warnings;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Index::Fasta; 
use Data::Dumper;

my $usage = "perl add_homologs_to_aln.pl alignment.fa database.fa > out to: alignment.hits.fa\n";
my $alignment = $ARGV[0] or die $usage;
my $database = $ARGV[1] or die $usage;
my $num_hits = "2"; # max number of final hits to retain per database
my $evalue = "1e-6"; # e-value threshold for BLASTP
my $threads = "4"; # number of threads for BLASTP


# 1. Input files

# ALIGNMENT
my $alignment_unaligned = $alignment . ".nogaps";

# check if unaligned file exists already (-z checks if file has 0 size)
if ( !-z "$alignment_unaligned") {

	system("sed '/>/! s/-//g' $alignment > $alignment_unaligned");
}

# DATABASE (TRANSCRIPTOME)
# check if unaligned file exists already
my $database_flags = $database . ".flags";
my %database;

if ( !-z "$database_flags") {

	# reformat alignment by removing gaps and adding flags to header (also made unique)
	print STDERR "\nReformating input alignment: $alignment\nrm gaps and add flags to unique headers...\n";

	%database = read_fasta($database);

	open(DB, ">", $database_flags) or die "Can't create $database_flags: $!";

	foreach my $d ( sort keys %database ) {
	
		print DB ">$d\n";
		print DB "$database{$d}\n";
	}
	close(DB);
}

# 2. blast each sequence in alignment against database

# check if blast db already exists, otherwise create
if ( !-z "$database_flags.phr") {

	print STDERR "\nCreating BLASTP database from: $database_flags...\n";
	system("makeblastdb -in $database_flags -dbtype prot -parse_seqids");
}
else {
	
	print STDERR "\nBLASTP database already exists for: $database_flags...\n";
}

# blastp
print STDERR "\nRunning BLASTP...\n";

my $blast_out = $alignment_unaligned . "_vs_" . $database_flags . ".blastp";

system("blastp -query $alignment_unaligned -db $database_flags -evalue $evalue -num_descriptions 1 -num_alignments 1 -num_threads $threads -out $blast_out");

# 3. parse blast report (save x best hits)
print STDERR "\nParsing BLASTP, retained $num_hits hits...\n";

my $report = new Bio::SearchIO(	-format => "blast", 
				-file   => "<$blast_out",
	);

my %hits;

while ( my $result = $report->next_result ) {

    while ( my $hit = $result->next_hit ) {

		my $hit_name = $hit->name;
		my $evalue = $hit->significance;

		# fill hash until $num_hits
		next if ( exists $hits{$hit_name} );
		
		if ( scalar keys %hits < $num_hits ) {
			
			$hits{$hit_name} = $evalue;
		}
		# replace hit if new hit has lower evalue
		else {
			foreach my $key ( keys %hits ) {
				
				if ( $evalue < $hits{$key} ) {
					
					delete $hits{$key};
					$hits{$hit_name} = $evalue;
				}
			}	
		}
	}
}

# 4. extract hits and create new output file

# extract hits and generate output
my $out_hits = $alignment . ".hits.fa";

# for checks: alternative print hits to individual files per taxa
# my $out_hits = $database . ".hits.fa";

print STDERR "\nExtracting best hits from: $database_flags\nOutput appended to: $out_hits\n";

open(OUT, ">>", $out_hits) or die "Can't create $out_hits: $!\n'";

foreach my $h ( keys %hits ) {
	
	if ( exists $database{$h} ) {
		
		print OUT ">$h\n";
		print OUT "$database{$h}\n";
	}
}
close(OUT);


# these three files can be reused if not removed every time
#system("rm $database_flags");
#system("rm $alignment_unaligned");
#system("rm $database_flags.p*"); # removes blast DB

# when running pipeline for single species, we can add new sequences into the alignment and
# re-align with mafft. When doing it for multiple species, it is easier to create one file
# per alignment where we add all hits in fasta, later to be added to the alignment
=pod
# add to alignment
my $new_alignment = $alignment . ".hits";
my $new_mafft = $new_alignment . ".mafft";

system("cat $alignment_unaligned $out_hits > $new_alignment");


# 5. realign

system("mafft $new_alignment > $new_mafft");

# clean up the temp files
system("rm $alignment_unaligned");
system("rm $new_alignment");
system("rm $database_flags*"); # removes fasta and blast DB
=cut

print STDERR "\n\ndone!\n\n";


##########################
sub read_fasta {

	my $fasta = $database;
	my %seq_hash;
	my $repeated = 0;

	my $seqio_obj = Bio::SeqIO->new('-file' => "<$fasta",
									'-format' => "fasta");

	while (my $seq_obj = $seqio_obj->next_seq) {

		my $seqname = $seq_obj->primary_id;
		my $sequence = $seq_obj->seq;
		my $length = $seq_obj->length;
		
		# add sequence length and "***" to new sequences
		my $seqname_with_flag = $seqname . "-len=$length-***";
		
		# remove gaps to sequence, if present
		if ( $sequence =~ /\w*[-]+\w*/ ) {

			$sequence =~ s/-//g;
		}
		# rename if duplicated headers are left (only if two sequences have same length)
    	if ( !exists $seq_hash{$seqname_with_flag} ) {

			$seq_hash{$seqname_with_flag} = $sequence;
	    }
		else {
	
        $repeated++;

		my $new_seqname = $seqname_with_flag . "_" . $repeated;

		$seq_hash{$new_seqname}	= $sequence;
		}
	}
	return %seq_hash;
}
