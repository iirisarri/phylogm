#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::SearchIO;
use Data::Dumper;
use Bio::DB::Fasta;


# Iker Irisarri. Konstanz, Jul 2015  
# modified version of parse_blast_RBHH_to_prot.pl
# it checks for reciprocal blast best hits
# option to use a filter for an evalue threshold
# IMPORTANT: It will output hits from $blast_report_1 (order in which reports are given is important!)
# originally for getting comps of reciprocal best hits to extract from fasta files (transcriptomes)
# It does not check for frame shifts, so it can be used for all flavours of BLAST
# Allows repeated queries (e.g. several ORF estiamted by transdecoder on the same transcfript) by appending a number

my $usage = "parse_blast_RBBH.pl blast_report_1 blast_report_2 evalue_threshold > output.fa\n";
my $blast_report_1 = $ARGV[0] or die $usage;
my $blast_report_2 = $ARGV[1] or die $usage;
my $evalue_threshold = $ARGV[2] or die $usage;

# process best hits from blast reports by sending it the subroutine (= blast_to_best_hit_prot.pl)
print STDERR "\nprocessing blast report #1...\n";

my $blast_1_ref = blast_to_hash($blast_report_1);

print STDERR "processing blast report #2...\n";
my $blast_2_ref = blast_to_hash($blast_report_2);

# de-reference hashes
my %blast_1 = %{ $blast_1_ref };
my %blast_2 = %{ $blast_2_ref };

#print Dumper \%blast_1;
#print Dumper \%blast_2;

# count original number of queries in blast report #1
my $original_query_num = scalar keys %blast_1;

#print Dumper \$blast_1_ref;
#print Dumper \$blast_2_ref;

my %RBBH_results;

print STDERR "searching for RBBH...\n";

foreach my $query_1 ( keys %blast_1 ) {

    my $best_hit_1 = $blast_1{$query_1}{'blast_info'}[0];
    
    # find BH1 as query in second blast report
    if ( exists $blast_2{$best_hit_1} ) {

	# if BH2 is not equal to query in first blast report, remove that key value pair

	if ( $query_1 ne $blast_2{$best_hit_1}{'blast_info'}[0] ) {  
	
	    delete $blast_1{$query_1};

	}
    }
}


#print Dumper \%blast_1;

# print out table of RBBH

foreach my $key ( keys %blast_1 ) {

    #print $key, "\t", $blast_1{$key}{'blast_info'}[0], "\n";
    print $key, "\t", $blast_1{$key}{'blast_info'}[0], "\t", $blast_1{$key}{'blast_info'}[2], "\n";

}

# count number of queries in blast report #1 after deleting non-RBBH 
my $query_num_after_delete = scalar keys %blast_1;
my $percent = ($query_num_after_delete / $original_query_num ) * 100;

print STDERR "original queries: $original_query_num\n";
print STDERR "number of RBBHs: $query_num_after_delete ($percent %)\n";


### SUBROUTINES ###

sub blast_to_hash {

    my $blast = shift;
    # appends number if queries are repeated
    my $query_count = 0;

    # read blast report
    my $report = new Bio::SearchIO(-format => "blast", 
				   -file   => "<$blast"
	);

    my %hash;
    #my @frames;

	# initialize evalues as 1, frames as 0
	my $significance1 = 1;
	my $evalue1 = 1;
	#my $frame_null = 0;
	#my $frame;

	while( my $result = $report->next_result ) {

		# get query name
		my $query = $result->query_name;
		
		# reasign 1 to evalue for each new result 
		$significance1 = 1;

		while( my $hit = $result->next_hit ) {

			# get evalue for that hit
			my $significance = $hit->significance;
			
			if ( $significance < $evalue_threshold ) {

			    # compare evalue of this hit with that from previous hits
			    # for first hit, compares against evalue of 1
			    # for next hits, compares with evalue of previous hit
			    if ( $significance < $significance1 ) {

				# assing new evalue to $value1 (for comparison)
				$significance1 = $significance;

				my $hit_name = $hit->name;
				# to remove "lcl|" that appears in hit name
				#$hit->name =~ /lcl\|(.+)/;
				#my $hit_name = $1;
				my $bits1 = $hit->bits;

				# empty $frame1 & @frames for each new hit
				#$frame_null = 0;
				#my $frame = 0;
				#@frames = ();

				$evalue1 = 1;

				# get reading frame
		#		while ( my $hsp = $hit->next_hsp ) {
		
		#		my $evalue = $hsp->evalue;
		#		if ( $evalue < $evalue1 ) {
		#		    $evalue1 = $evalue;

				# get reading frames for all hsps

				# process first hsp for each hit ( only time $frame_null will be 0 )
				#	if ( $frame_null == 0 ) {
					# get reading frame (1,2,3) multiplied by strand
				#		$frame = ( $hsp->hit->frame + 1 ) * $hsp->strand('hit');
				#		$frame_null = $frame;
						# store frame into array @frames
				#		push ( @frames, $frame );
				#	}
					# get further $frames in the loop
				#	else {

				#	$frame = ( $hsp->hit->frame + 1 ) * $hsp->strand('hit');
					
				#	}

				# make sure all frames are the same for the different hsp within the same hit
				#if ( $frame != $frame_null ) {
				# will print as many times as reading frame shifts are present for each hit
				#	print STDERR "Warning! Frameshift for",
				#	" query $query and hit $hit_name\n";

					# add additional frames to array @frames
				#	push ( @frames, $frame );
			
				#}
			
				# store info into hash of hashes
				# APPEND NUM TO KEY IF QUERIES ARE REPEATED
				if ( !exists $hash{$query} ) {
				    
				    $hash{$query}{'blast_info'} = [$hit_name,$bits1,$significance];
				}
				else {
				
				    $query_count++;
				    my $query_repeated = $query . "_" . $query_count;
				    $hash{$query_repeated}{'blast_info'} = [$hit_name,$bits1,$significance];
				 }

				# hash 'frames' would allow translating the same query in multiple frames
				# to try to correct from frameshifts
				#$hash{$query}{'frames'} = [@frames];
			    }
			}
		}
	}

	# return parsed info of best hits as a hash ref
    return \%hash;

}

__END__
