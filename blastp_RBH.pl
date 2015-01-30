#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::SearchIO;
use Data::Dumper;
use Bio::DB::Fasta;


# Iker Irisarri. Konstanz, January 2015  
# modified version of parse_blast_best_hit_prot.pl to parse blastp reports (it does not translate hits)
# it checks for reciprocal blast best hits
# option to use a filter for an evalue threshold. If undefined as an option, evalue=1e-3
# OUTPUT: new output file created for each orthogroup containing the ortholog of the blasted species

# IMPORTANT: It will output hits from $blast_report_1 (order in which reports are given is important!)
# IMPORTANT: check how query and hit names are handled (perl -d) cause for different fasta files primary_id and hit->name will be different
# originally for calculating RBBH of vert queries vs. callorhinchus proteome and extracting seq from proteome


my $usage = "blastp_RBH.pl blast_report_1 blast_report_2 source_fasta (evalue threshold)> output.fa\n";
my $blast_report_1 = $ARGV[0] or die $usage;
my $blast_report_2 = $ARGV[1] or die $usage;
my $in_fasta = $ARGV[2] or die $usage;
my $evalue_threshold =  $ARGV[3];

if ( !defined($evalue_threshold) ) {
    
    $evalue_threshold = 1e-3;

}

# process best hits from blast reports by sending it the subroutine (= blast_to_best_hit_prot.pl)
print STDERR "\nprocessing blast report #1 $blast_report_1...\n";

my $blast_1_ref = blast_to_hash($blast_report_1);

print STDERR "processing blast report #2 $blast_report_2...\n";
my $blast_2_ref = blast_to_hash($blast_report_2);

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


# extract best hits from source fasta
while (my $seq_obj = $seqio_obj->next_seq){
    
    # primary id would be 'gi|632933636|ref|XP_007888002.1|'
    $seq_obj->primary_id =~ /gi\|\d+\|(ref\|\w+\.*\w*\|)/;
    my $seqname = $1;

#    foreach my $key ( sort { $a cmp $b } keys %blast_1 ) {
    foreach my $key ( keys %blast_1 ) {
    # get best hits from source fasta
	if ( $seqname eq $blast_1{$key}{'blast_info'}[0] ) {

	    # get ortholog group
	    $key =~ /(EOG\w+)/;
	    my $orthogroup = $1;
	    my $outfile = $orthogroup . "_callorhinchus.fa";

	    open (OUT, ">", "$outfile");

            print OUT ">Callorhinchus_milii\n";
	    print OUT $seq_obj->seq, "\n";
	    # scape foreach loop of %blast_1 once the sequence is found
	    next;
	    close(OUT);
	}
    }
}

print STDERR "\ndone!\n\n";

### SUBROUTINES ###

sub blast_to_hash {

    my $blast = shift;

    # read blast report
    my $report = new Bio::SearchIO(-format => "blast", 
				   -file   => "<$blast"
	);

    my %hash;

	# initialize evalues as threshold
	my $significance1 = $evalue_threshold;
	my $evalue1 = $evalue_threshold;

	while( my $result = $report->next_result ) {

	    # get query name and extract ref|| name if present
            # blastp output of queries vs callorhinchus will have gi|\d+|ref|\w+|
	    # but hits in callorhinchus vs queries will have only |ref|\w+|
	    my $query;
            if ( $result->query_name =~ /gi\|\d+\|ref\|\w+\.*\w*\|/ ) {

		$result->query_name =~ /gi\|\d+\|(ref\|\w+\.*\w*\|)/;
#		$result->query_name =~ /(ref.+)/;
		$query = $1;

	    }
	    else {

		$query = $result->query_name;

	    }

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
			    $hash{$query}{'blast_info'} = [$hit_name,$bits1,$significance];

			}
		}
	}

    # return parsed info of best hits as a hash ref
    return \%hash;

}

__END__
