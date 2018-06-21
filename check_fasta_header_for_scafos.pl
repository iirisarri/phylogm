#!/usr/bin/perl                                                                                                                                                     

use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;

##########################################################################################                                                                          
#                                                                                                                                                                   
# ##### Iker Irisarri. June 2018. Museo Nacional de Ciencias Naturales University ########
#                                                                                                                                                                   
# It checks that fasta headers do not contain weird characters that are not compatible
#	with Scafos
#
# Not accepted: \ | ( ) [ ] { } ^ $ * + ? ~ ! % & = ' ; : , < > #                                                                                                       
#                                                                                                                                                                   
##########################################################################################                                                                          


my $usage = "perl check_fasta_header_for_scafos.pl *.fasta\n";
my @infiles = @ARGV or die $usage;

open (OUTERR, ">>check_fasta_header_for_scafos.err") or die "Can't create check_fasta_header_for_scafos.err\n";

foreach my $infile ( @infiles ) {

    my $seqio_obj = Bio::SeqIO->new('-file' => "<$infile",
                                    '-format' => "fasta",
        );

	print OUTERR "FILE: $infile\n\n";
	
	my %fasta = ();

	my $duplicated = "0"; # counter to make headers unique in case they are repeated
	
    while (my $inseq = $seqio_obj->next_seq) {

		my $seqname = $inseq->primary_id;
		my $sequence = $inseq->seq;
        
        if ( $seqname =~ /.*?\\.*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\\/_/;
	        print OUTERR "\t$seqname\n\n";     
        }
        if ( $seqname =~ /.*?\/.*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\\/_/;
	        print OUTERR "\t$seqname\n\n";     
        }
        if ( $seqname =~ /.*?\|.*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\|/_/;
	        print OUTERR "\t\t$seqname\n\n";     
        }
        if ( $seqname =~ /.*?\(.*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\(/_/;
	        print OUTERR "\t$seqname\n\n";     
        }
        if ( $seqname =~ /.*?\).*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\)/_/;
	        print OUTERR "\t$seqname\n\n";     
        }
        if ( $seqname =~ /.*?\[.*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\[/_/;
	        print OUTERR "\t$seqname\n\n";     
        }
        if ( $seqname =~ /.*?\].*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\]/_/;
	        print OUTERR "\t$seqname\n\n";     
        }
        if ( $seqname =~ /.*?\{.*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\{/_/;
	        print OUTERR "\t$seqname\n\n";     
        }
        if ( $seqname =~ /.*?\}.*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\}/_/;
	        print OUTERR "\t$seqname\n\n";     
        }
        if ( $seqname =~ /.*?\^.*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\^/_/;
	        print OUTERR "\t$seqname\n\n";     
        }
        if ( $seqname =~ /.*?\$.*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\$/_/;
	        print OUTERR "\t$seqname\n";     
        }
        if ( $seqname =~ /.*?\*.*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\*/_/;
	        print OUTERR "\t$seqname\n\n";     
        }
        if ( $seqname =~ /.*?\+.*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\+/_/;
	        print OUTERR "\t$seqname\n\n";     
        }
        if ( $seqname =~ /.*?\?.*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\?/_/;
	        print OUTERR "\t$seqname\n\n";     
        }
        if ( $seqname =~ /.*?\~.*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\~/_/;
	        print OUTERR "\t$seqname\n\n";     
        }
        if ( $seqname =~ /.*?\!.*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\!/_/;
	        print OUTERR "\t$seqname\n\n";     
        }
        if ( $seqname =~ /.*?\%.*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\%/_/;
	        print OUTERR "\t$seqname\n\n";     
        }
        if ( $seqname =~ /.*?\&.*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\&/_/;
	        print OUTERR "\t$seqname\n\n";     
        }
        if ( $seqname =~ /.*?\=.*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\=/_/;
	        print OUTERR "\t$seqname\n\n";     
        }
        if ( $seqname =~ /.*?\'.*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\'/_/;
	        print OUTERR "\t$seqname\n\n";     
        }
        if ( $seqname =~ /.*?\;.*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\;/_/;
	        print OUTERR "\t$seqname\n\n";     
        }
        if ( $seqname =~ /.*?\:.*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\:/_/;
	        print OUTERR "\t$seqname\n\n";     
        }
        if ( $seqname =~ /.*?\,.*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\,/_/;
	        print OUTERR "\t$seqname\n\n";     
        }
        if ( $seqname =~ /.*?\<.*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\</_/;
	        print OUTERR "\t$seqname\n\n";     
        }
        if ( $seqname =~ /.*?\>.*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\>/_/;
	        print OUTERR "\t$seqname\n\n";     
        }
        if ( $seqname =~ /.*?\#.*?/g ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname =~ tr/\#/_/;
	        print OUTERR "\t$seqname\n\n"; 
        }
        # remove if header finishes is "_"
        if ( $seqname =~ /(.+?)_+$/ ) {
        
        	print OUTERR "\t$seqname\n\n\t\treplaced by:\n\n";
        	$seqname = $1;
	        print OUTERR "\t$seqname\n\n"; 
        }
        # save sequence with new header
        if ( exists $fasta{$seqname} ) {
        
    		print OUTERR "\tWARN: duplicated header $seqname\n";
			$duplicated++;
			$seqname .= "_$duplicated";
			print OUTERR "\t\tsaved as: $seqname\n\n";

    	}    
		$fasta{$seqname} = $sequence;
	}
	# print out file
	my $outfile = $infile . ".out";
	
	open(OUT, ">", $outfile) or die "Can't create $outfile\n";
	
	foreach my $key ( sort keys %fasta ) {
	
		print OUT ">$key\n";
		print OUT "$fasta{$key}\n";
	}
}

