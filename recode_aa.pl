#!/usr/bin/perl
 
use strict;
use warnings;
use Bio::SeqIO;
use Data::Dumper;

############################################################################
#
# Iker Irisarri, Sep 2018. Museo Nacional de Ciencias Naturales MNCN-CSIC
# 
# It recodes amino acids to different schemes: dayhoff, 
# Susko & Roger 2007's 6-category and Kosiol et al. 2004's 6-category
# 
#
############################################################################

my $usage = "recode_aa.pl infile.fa\n";
my @infiles = @ARGV or die $usage; 

#dayhoff (AGPST) (DENQ) (HKR) (MIVL) (WFY) (C)
#SR-6 (APST) (DENG) (QKR) (MIVL) (WC) (FYH)
#KGB-6 (AGPS) (DENQHKRT) (MIL) (W) (FY) (CV)

# loop through the list of infiles
foreach my $infile (@infiles) {

	my $new_seq_length = 0;
	my %DAYH6 = ();
	my %SR6 = ();
	my %KGB6 = ();
	
	# read file in with seqio
	my $seqio_obj = Bio::SeqIO->new(
		'-file' => "<$infile", 
		'-format' => "fasta", 
		);

	while (my $inseq = $seqio_obj->next_seq) {

		my $sequence = $inseq->seq;
		
		my $dayh6 = uc $sequence;
		my $sr6 = uc $sequence;
		my $kgb6 = uc $sequence;

		# recode
		$dayh6 =~ tr/AGPSTDENQHKRMIVLWFYC/11111222233344445556/; # Dayhoff
		$sr6   =~ tr/APSTDENGQKRMIVLWCFYH/11112222333444455666/; # SR-6
		$kgb6  =~ tr/AGPSDENQHKRTMILWFYCV/11112222222233345566/; # KGB-6

		$DAYH6{$inseq->primary_id}=$dayh6;
		$SR6{$inseq->primary_id}=$sr6;
		$KGB6{$inseq->primary_id}=$kgb6;
	}

	# print out output files

	my $outfile1 = $infile . ".dayhoff6.fa";
	my $outfile2 = $infile . ".SR6.fa";
	my $outfile3 = $infile . ".KGB6.fa";
	
	open(OUT1, ">", $outfile1);
	open(OUT2, ">", $outfile2);
	open(OUT3, ">", $outfile3);

	foreach my $k1 ( sort keys %DAYH6 ) {

		print OUT1 ">$k1\n";
		print OUT1 $DAYH6{$k1}, "\n";
	}
	foreach my $k2 ( sort keys %SR6 ) {

		print OUT2 ">$k2\n";
		print OUT2 $SR6{$k2}, "\n";
	}
	foreach my $k3 ( sort keys %KGB6 ) {

		print OUT3 ">$k3\n";
		print OUT3 $KGB6{$k3}, "\n";
	}
}

print STDERR "\ndone!\n\n";