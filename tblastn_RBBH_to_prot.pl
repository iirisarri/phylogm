#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::SearchIO;
use Data::Dumper;
use Bio::DB::Fasta;


# Iker Irisarri. Konstanz, January 2015  
# modified version of parse_blast_best_hit_prot.pl
# it identifies RBBH and translates nt sequence into aa
# REQUIRES: two blast reports (e.g. prot vs nucl –tblastn– and nucl vs prots – blastx–), source fasta and taxon name to output in fasta files
# OPTIONAL use a filter for an evalue threshold. If undefined as an option, evalue=1e-3
# OUTPUT: RBBH as new fasta files, one for each orthogroup

# IMPORTANT: It will output hits from $blast_report_1 (order in which reports are given is important!)
# originally for getting comps of reciprocal best hits to extract from fasta files (transcriptomes)
# and translate them using the correct reading frame
# IMPORTANT: define fasta header to be printed into output files

my $usage = "tblastn_RBBH_to_prot.pl tblastn_rep_1 blastx_rep_2 source_fasta species_name (evalue threshold) > output.fa\n";
my $blast_report_1 = $ARGV[0] or die $usage;
my $blast_report_2 = $ARGV[1] or die $usage;
my $in_fasta = $ARGV[2] or die $usage;
my $header_out = $ARGV[3] or die $usage;
my $evalue_threshold =  $ARGV[4];


if ( !defined($evalue_threshold) ) {

    $evalue_threshold = 1e-3;

}


# process best hits from blast reports by sending it the subroutine (= blast_to_best_hit_prot.pl)
print STDERR "\nprocessing blast report #1 $blast_report_1...\n";

my $blast_1_ref = tblastn_to_hash($blast_report_1);

print STDERR "\nprocessing blast report #2 $blast_report_2...\n";
my $blast_2_ref = blastx_to_hash($blast_report_2);

# de-reference hashes
my %blast_1 = %{ $blast_1_ref };
my %blast_2 = %{ $blast_2_ref };

my $original_query_num = `grep -c "Query=" $blast_report_1`;

# count original number of queries that have significant hits in blast report #1                            
my $signif_query_num = scalar keys %blast_1;

#print Dumper \$blast_1_ref;
#print Dumper \$blast_2_ref;

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

# count number of queries in blast report #1 after deleting non-RBBH 
my $query_num_after_delete = scalar keys %blast_1;
my $percent = ($query_num_after_delete / $signif_query_num ) * 100;

print STDERR "\noriginal number of queries: $original_query_num\n";
print STDERR "of which they have significant hits: $signif_query_num\n";
print STDERR "number of significant RBBHs: $query_num_after_delete ($percent %)\n\n";
print STDERR "extracting RBBH (blast report #1) from source fasta...\n";

###########################################
# extract_from_fasta_by_name && translate #
###########################################

# read fasta file with SeqIO
my $seqio_obj = Bio::SeqIO->new('-file' => "<$in_fasta",
         			'-format' => 'fasta');

my $prot_obj = Bio::SeqIO->new('-format' => "fasta");

my $perl_frame;

#my $header_out;

# extract best hits from source fasta & translate
while (my $seq_obj = $seqio_obj->next_seq){

    my $seqname = $seq_obj->primary_id;
    #my $description = $seq_obj->description;

    foreach my $key ( keys %blast_1 ) {
	# get best hits from source fasta
	if ( $seqname eq $blast_1{$key}{'blast_info'}[0] ) {
	    # get all elements from @frames
	    foreach my $real_frame ( @{ $blast_1{$key}{'frames'} } ) {

			# get "real" reading frame and convert into perl reading frames by doing -1 (0,1,2)
			$real_frame =~ /-*(\d)/;
			$perl_frame = $1 - 1;
	#		$header_out = "$key";

			# revcomp sequence if negative frame
			if ( $real_frame =~ /-\d/ ) {
		    
			    my $revcomp_obj = $seq_obj->revcom;
		    
			    # codontable_id 1 is for standard genetic code
			    $prot_obj = $revcomp_obj->translate(-codontable_id => 1,
								-frame => $perl_frame);
			}
			else {

		    	$prot_obj = $seq_obj->translate(-codontable_id => 1,
							-frame => $perl_frame);
			}

			# get ortholog group
	   	 	$key =~ /(EOG\w+)/;
			my $orthogroup = $1;
			my $outfile = $orthogroup . "_" . $header_out . ".fa";

			open (OUT, ">", "$outfile");

			print OUT ">$header_out\n";
			print OUT $prot_obj->seq, "\n";
			# scape foreach loop of %blast_1 once the sequence is found
			next;
			close(OUT);
	    }
	}
    }
}

print STDERR "\ndone!\n\n";

### SUBROUTINES ###

sub tblastn_to_hash {

    my $blast = shift;

    # read blast report
    my $report = new Bio::SearchIO(-format => "blast", 
				   -file   => "<$blast");

    my %tblastn;
    my @frames;

    #initialize evalues as 1, frames as 0
    my $significance1 = $evalue_threshold;
    my $evalue1 = $evalue_threshold;
    my $frame_null = 0;
    my $frame;

    while( my $result = $report->next_result ) {

	# get query name
	my $query = $result->query_name;

	# reasign 1 to evalue for each new result 
	$significance1 = $evalue_threshold;

	while( my $hit = $result->next_hit ) {

		# get evalue for that hit
		my $significance = $hit->significance;

		# compare evalue of this hit with that from previous hits
		# for first hit, compares against evalue of 1-e3
		# for next hits, compares with evalue of previous hit
		if ( $significance < $significance1 ) {

			# assing new evalue to $value1 (for comparison)
			$significance1 = $significance;
			my $bits1 = $hit->bits;

			# to remove "lcl|" that appears in hit name
			my $hit_name;
			if ( $hit->name =~ /lcl\|.+/ ) {
			    $hit->name =~ /lcl\|(.+)/;
			    $hit_name = $1;
			}
                
			else {

			    $hit_name = $hit->name;

			}


			# empty $frame1 & @frames for each new hit
			$frame_null = 0;
			my $frame = 0;
			@frames = ();

			$evalue1 = $evalue_threshold;

			# get reading frame
			while ( my $hsp = $hit->next_hsp ) {
		
			# get reading frames for all hsps

			# process first hsp for each hit ( only time $frame_null will be 0 )
			    if ( $frame_null == 0 ) {
				# get reading frame (1,2,3) multiplied by strand
				$frame = ( $hsp->hit->frame + 1 ) * $hsp->strand('hit');
				$frame_null = $frame;
				# store frame into array @frames
				push ( @frames, $frame );
			    }
			    # get further $frames in the loop
			    else {

				$frame = ( $hsp->hit->frame + 1 ) * $hsp->strand('hit');
					
			    }

			    # make sure all frames are the same for the different hsp within the same hit
			    if ( $frame != $frame_null ) {
			    
				#will print as many times as reading frame shifts are present for each hit
				print STDERR "Warning! Frameshift for",
				" query $query and hit $hit_name\n";

				# add additional frames to array @frames
				push ( @frames, $frame );
			
			    }
			
			    # store info into hash of hashes
			    $tblastn{$query}{'blast_info'} = [$hit_name,$bits1,$significance,$frame];
			    # hash 'frames' would allow translating the same query in multiple frames
			    # to try to correct from frameshifts
			    $tblastn{$query}{'frames'} = [@frames];
		       	}
		}
	}
    }

    # return parsed info of best hits as a hash ref
    return \%tblastn;

}


sub blastx_to_hash {

    my $blast = shift;

    # read blast report
    my $report = new Bio::SearchIO(-format => "blast", 
				   -file   => "<$blast"
	);

    my %blastx;

	# initialize evalues as threshold
	my $significance1 = $evalue_threshold;
	my $evalue1 = $evalue_threshold;

	while( my $result = $report->next_result ) {

	    my $query = $result->query_name;

	    # reasign evalue to the threshold for each new result 
	    $significance1 = $evalue_threshold;

	    while( my $hit = $result->next_hit ) {

		# get evalue for that hit
		my $significance = $hit->significance;

		# compare evalue of this hit with that from previous hits
		# for first hit, compares against evalue of 1
		# for next hits, compares with evalue of previous hit
		if ( $significance < $significance1 ) {

		    # assing new evalue to $value1 (for comparison)
		    $significance1 = $significance;

		    # to remove "lcl|" that appears sometimes in hit name
		    # blastp output of callorhinchus vs queries will have lcl| before the hits 
				
		    my $hit_name;
		    if ( $hit->name =~ /lcl\|.+/ ) {
			$hit->name =~ /lcl\|(.+)/;
			$hit_name = $1;
		    }
		    else {
			
			$hit_name = $hit->name;
			
		    }
			
		    my $bits1 = $hit->bits;

		    $evalue1 = $evalue_threshold;

		    # store info into hash of hashes
		    $blastx{$query}{'blast_info'} = [$hit_name,$bits1,$significance];

		}
	    }
	}

    # return parsed info of best hits as a hash ref
    return \%blastx;

}

__END__
