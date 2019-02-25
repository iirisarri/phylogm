#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;

# Iker Irisarri, November 2016. University of Konstanz.
# 
# Obtains basic stats from a fasta alignment: variable, invariable and parsimony-informative sites
# 
# NOTES:
# By default, gaps considered as missing data (useful for concatenated matrices)
# If we want gaps to count as variability uncomment line num 126.
# Parsimony informative sites: min 2 chars, each present in at least 2 sequences
# 	this requires at least 4 sequences. If not present, these site category are ignored
# 	when calculating parsimony-informative sites, we assume each sequence is a different taxa 
# Alignment column numbers are printed as 0-based

my $usage = "aln_stats_var_pi_sites.pl in.fasta DNA/PROT> STDOUT\n";
my $in_fasta = $ARGV[0] or die $usage;
my $data_type = $ARGV[1] or die $usage;



# define ambiguity
my $ambiguous = "n"; # default DNA

if ( $data_type eq "PROT" ) {

	$ambiguous = "?";
}
elsif ( $data_type ne "DNA" ) {

	die $usage;
}

## Store alignment in hash of arrays
open(IN, "<", $in_fasta);

my %fasta_hash;
my $header;
my $num_seqs = 0;
my @seq;
my $length = 0;

# STORE SEQUENCES
=pod
# read fasta file with SeqIO
my $seqio_obj = Bio::SeqIO->new('-file' => "<$in_fasta",
								'-format' => "fasta");

while ( my $seq_obj = $seqio_obj->next_seq ) {

    $header = $seq_obj->primary_id;
    my $sequence = $seq_obj->seq;
    
    @seq = split (//, $sequence);
    
	$fasta_hash{$header} = [@seq];

}
=cut

# fasta parser without bioperl
while ( my $line =<IN> ) {
	
	chomp $line;
	
	if ($line =~ />(.+)/ ) {
	
		$header = $1;
		$num_seqs++;
	}
	else {
		# convert everything to lowercase
		$line = lc $line;
		
		@seq = split (//, $line);

		# check that sequences are aligned
		if ( $length == 0 ) {

			$length = scalar @seq;
		} 
		else {

			die "Sequences are not aligned, check $in_fasta!" if ( $length != scalar @seq );
		}
	} 
	$fasta_hash{$header} = [@seq];

}

#print Dumper \%fasta_hash;

# DECLARE PATTERN TYPES
my @indet_only;
my @indet_invar;
my @indet_var;
my @indet_pi;
my @indet_total;
my @det_invar;
my @det_var;
my @det_pi;
my @det_total;

# SAVE PATTERNS

# loop through each aln position
for ( my $i=0; $i<$length + 1; $i++ ) {

	my %pattern = ();
	my $first_char = ();

	# process pattern of column $i by counting the frequency of each char
	foreach my $head ( keys %fasta_hash ) {

		my $nt = $fasta_hash{$head}[$i];
		
		# skip undefined $nt (for some reason at the end of the aln I had an undefined $nt in all sequences)
		next if ( !defined $nt );
		
		# COUNT GAPS AS MISSING DATA (DEFAULT)
		# OTHERWISE, UNCOMMENT NEXT LINE
		$nt = $ambiguous if ( $nt eq "-");
		
		if ( !exists $pattern{$nt} ) {
		
			$pattern{$nt} = 1;
		}
		else {
			$pattern{$nt}++;
		}
	}	
	#print Dumper \%pattern;
	
	# CLASSIFY PATTERN & STORE IT
	
	# sort keys by its value (higher to lower)
	foreach my $l ( sort { $pattern{$b} <=> $pattern{$a} } keys %pattern ) {

		if ( exists $pattern{$ambiguous} ) {

			# INDETERM (with Ns)
			# only indet
			if ( $pattern{$ambiguous} == $num_seqs ) {
			
				push (@indet_only, $i);
				push (@indet_total, $i);
				last;
			}
			# indet invariable
			if ( scalar keys %pattern == 2 ) {
			
				push (@indet_invar, $i);
				push (@indet_total, $i);
				last;
			}
			# indet parsimony-informative (min 2 chars, each present in at least 2 sequences)
			# requires the presence of at least 4 sequences
			# skip ambiguous (n or ?) character for the identification of further sites
			# ( gaps "-" are considered as an extra state)
			next if ( $l eq $ambiguous);
			if ( !defined $first_char && $pattern{$l} > 1 && $num_seqs > 3 ) {
				
				$first_char = $pattern{$l};
				next;
			} 
			if ( defined $first_char && $pattern{$l} > 1 && $num_seqs > 3 ) {
			
				push (@indet_pi, $i);
				push (@indet_var, $i);
				push (@indet_total, $i);
				last;
			}
			# indet variable
			else {
			
				push (@indet_var, $i);
				push (@indet_total, $i);
				last;
			}
		}
		else {
		
			# DETETERM	(without Ns)
			# det invariable
			if ( scalar keys %pattern == 1 ) {
			
				push (@det_invar, $i);
				push (@det_total, $i);
				last;
			}
			# det parsimony-informative
			if ( !defined $first_char && $pattern{$l} > 1 && $num_seqs > 3 ) {
				
				$first_char = $pattern{$l};
				next;
			}
			if ( defined $first_char && $pattern{$l} > 1 && $num_seqs > 3 ) {
			
				#print "DET_PI:";
				#print Dumper \%pattern;
				push (@det_pi, $i);
				push (@det_var, $i);
				push (@det_total, $i);
				last;
			}
			# det variable
			else {
			
				push (@det_var, $i);
				push (@det_total, $i);
				last;
			}
		}
	}
}

# add up some counts
my $total_pos = scalar @det_total + scalar @indet_total;
my $det_total = scalar @det_total;
my $det_invar = scalar @det_invar;
my $det_var = scalar @det_var;
my $det_pi = scalar @det_pi;
my $indet_total = scalar @indet_total;
my $indet_invar = scalar @indet_invar;
my $indet_var = scalar @indet_var;
my $indet_pi = scalar @indet_pi;
my $indet_only = scalar @indet_only;
my $tot_invar = $det_invar + $indet_invar;
my $tot_var = $det_var + $indet_var;
my $tot_pi = $det_pi + $indet_pi;

my $perc_det_total = 0;
my $perc_det_invar = 0;
my $perc_det_var = 0;
my $perc_det_pi = 0;
my $perc_indet_total = 0;
my $perc_indet_invar = 0;
my $perc_indet_var = 0;
my $perc_indet_pi = 0;
my $perc_indet_only = 0;

if ( $det_total != 0 ) {
	$perc_det_total = ($det_total / $total_pos) * 100;
	$perc_det_invar = ($det_invar / $det_total) * 100;
	$perc_det_var = ($det_var / $det_total) * 100;
	$perc_det_pi = ($det_pi / $det_total) * 100;
}
if ( $indet_total != 0 ) {
	$perc_indet_total = ($indet_total / $total_pos) * 100;
	$perc_indet_invar = ($indet_invar / $indet_total) * 100;
	$perc_indet_var = ($indet_var / $indet_total) * 100;
	$perc_indet_pi = ($indet_pi / $indet_total) * 100;
	$perc_indet_only = ($indet_only / $indet_total) * 100;
}
my $perc_tot_invar = ($tot_invar / $total_pos) * 100;
my $perc_tot_var = ($tot_var / $total_pos) * 100;
my $perc_tot_pi = ($tot_pi / $total_pos) * 100;
my $perc_tot_indet_only = ($indet_only / $total_pos) * 100;

# print out statsistics
print "\nTotal number of positions:\t"; print "$total_pos\n";
print "Total number of sequences:\t$num_seqs\n\n";

#=pod
# print out counts
print "Complete positions:\t\t\t"; print scalar @det_total, " ($perc_det_total% of total)\n";
print "\tinvariable:\t\t\t\t"; print scalar @det_invar, " ($perc_det_invar%  of complete)\n";
print "\tvariable:\t\t\t\t"; print scalar @det_var, " ($perc_det_var%  of complete)\n";
print "\tparsimony-informative:\t"; print scalar @det_pi, " ($perc_det_pi%  of complete)\n";

print "Incomplete positions:\t\t"; print scalar @indet_total, " ($perc_indet_total% of total)\n";
print "\tinvariable:\t\t\t\t"; print scalar @indet_invar, " ($perc_indet_invar% of incomplete)\n";
print "\tvariable:\t\t\t\t"; print scalar @indet_var, " ($perc_indet_var% of incomplete)\n";
print "\tparsimony-informative:\t"; print scalar @indet_pi, " ($perc_indet_pi% of incomplete)\n";
print "\tindetermined only:\t\t"; print scalar @indet_only, " ($perc_indet_only% of incomplete)\n";

print "\nAll positions:\t\t\t\t"; print "$total_pos\n";
print "\tinvariable:\t\t\t\t"; print $tot_invar, " ($perc_tot_invar% of total)\n";
print "\tvariable:\t\t\t\t"; print $tot_var, " ($perc_tot_var% of total)\n";
print "\tparsimony-informative:\t"; print $tot_pi, " ($perc_tot_pi% of total)\n";
print "\tindetermined only:\t\t"; print scalar @indet_only, " ($perc_tot_indet_only% of total)\n";
#=cut

=pod
# print out aln columns (counts)
print "Complete positions:\t\t"; print join (",", @det_total), "\n";
print "\tinvariable:\t\t"; print join (",", @det_invar), "\n";
print "\tvariable:\t\t"; print join (",", @det_var), "\n";
print "\tparsimony-informative:\t"; print join (",", @det_pi), "\n";

print "Incomplete positions:\t\t";print join (",", @indet_total), "\n";
print "\tinvariable:\t\t"; print join (",", @indet_invar), "\n";
print "\tvariable:\t\t"; print join (",", @indet_var), "\n";
print "\tparsimony-informative:\t"; print join (",", @indet_pi), "\n";
print "\tindetermined only:\t"; print join (",", @indet_only), "\n";
=cut

print STDERR "\ndone!\n\n";

