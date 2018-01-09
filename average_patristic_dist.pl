#!/usr/bin/perl
 
use strict;
use warnings;
use Data::Dumper;

##########################################################################################
#
# #################### Iker Irisarri. Jan 2018. Uppsala University ##################### #
#
# Average pairwise patristic distances using .mldist files form IQTREE as input
#	Lower-triangle values are selected for computation (including diagonal "0" or
#	duplicated values in upper-triangle changes average value
#
##########################################################################################

my $usage = "average_patristic_distances.pl infile.mldist > STDOUT\n";
my $infile = shift or die $usage; 

my @data;
my $data_lines;

open (IN, "<", $infile);

while ( my $line =<IN> ) {

	chomp $line;
	
	# skip first line
	next if $line =~ /^\d+$/;
	
	my @lines = split (" ", $line);

	# remove everything in the upper triangle
	# I9169_1_80 0 10 20 30 -> 1st line ($data_lines=1), splice will not extract any number -> skip
	# I9129_1_40 10 0 10 20	-> 2nd line: we retain from element 1 to 1
	# I9172_1_83 20 10 0 10	-> 3rd: we retain 1-2
	# I9097_1_8_ 30 20 10 0 -> 4th: we retain 1-4

	$data_lines++;
	my $last_elem_in_array = $data_lines - 1;
	
	# skip first line
	next if $last_elem_in_array == 0;

	# extract lower-triangle values
	@lines = splice @lines, 1, $last_elem_in_array;

	# checking
	#print "$line\n";
	#print "LINE: $data_lines\tELEM: 1 - $last_elem_in_array:\t";
	#print join "--", @lines, "\n";
	
	# save rest of elements in array
	foreach my $elem ( @lines ) {
	
		push ( @data, $elem );
	}
}

# calculate average patristic distances and print out result
my $num_points = scalar @data;

my $total;

foreach my $d ( @data ) {

	$total += $d;
}

my $average = $total / $num_points;

print "$infile\t$average\n";


