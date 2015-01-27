#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Bio::DB::Fasta;
use Bio::SeqIO;

# Iker Irisarri, January 2015, University of Konstanz
# Modification of filter_orthodb_single_copy.pl
# Added: extract prot sequences from ensembl proteome files and create fasta files for each orthogroup
# OrthoDB is based on Ensembl v75, so use those proteomes!!
# If the protein is not found in the fasta file it will print the following to STDERR:
# Use of uninitialized value $seq in print at /computation/iker/software/scripts/phylogm/filter_orthodb_2fasta.pl line 155.


# From filter_orthodb_single_copy.pl :
# Script for parsing OrthoDB (v.8) information containing ortholog groups, species and genes
# info available @ ftp://cegg.unige.ch/OrthoDB8/Eukaryotes/Genes_to_OGs/
# It prints out (i) ortholog group, (ii) species name and (iii) ensembl ID for single-copy genes that are present in the reference set of species
# It requires the orthoDB file with the relationships of orthologs and a text file with the names of species for which we require genes to be present as single copy

## OPTIONS ##
# There is the possibility of requiring single copy genes to be present in (i) ALL ref species
# or (ii) only in a subset (e.g. >70%, >80%, <90%). This option is hard-coded. 
# To select option 1 or 2, go the specific line and comment out the if statement that that does not apply

my $usage = "filter_orthodb_single_copy.pl orthodb_file ref_species fasta_source_file > stdout\n";

my $orthodb_file = $ARGV[0];
my $reference_set = $ARGV[1];
my $fasta = $ARGV[2];

# create fasta DB && read fasta file with SeqIO
my $db = Bio::DB::Fasta->new("$fasta");

my $seqio_obj = Bio::SeqIO->new('-file' => "<$fasta",
                             '-format' => "fasta");


# read in reference species & save it in %hash
open (DB, "<", $reference_set) or die $usage;

my %ref_species;

while ( my $line0 = <DB> ) {
    chomp $line0;

    # change spaces by underescores and remove parentheses in species names
    $line0 =~ s/\s/_/g;
    $line0 =~ s/[\(\)]//g;

    $ref_species{$line0} = 1;
}

# get number of referenc species (number of keys in hash)
my $ref_species_length = scalar keys %ref_species;


# read-in orthoDB file & save it in %hash
open (IN, "<", $orthodb_file) or die $usage;

my %orthoDB;
my $orthogroup;
my $prot_id;
my $species;

while ( my $line1 = <IN> ) {
    chomp $line1;

    next if ($line1 =~ /odb8_level/ || $line1 =~ /odb8_og_id/ || $line1 =~ /protein_id/ || $line1 =~ /odb8_prot_id/ || $line1 =~ /organism/ || $line1 =~ /uniprot_description/ );

    my @lines = split ("\t", $line1);

    $orthogroup = $lines[1];
    $prot_id = $lines[2];
    $species = $lines[4];
    # change spaces by underescores and remove parentheses in species names
    $species =~ s/\s/_/g;
    $species =~ s/[\(\)]//g;

    if ( !exists $orthoDB{$orthogroup}{$species} ) {

	$orthoDB{$orthogroup}{$species} = [$prot_id];

    }
    else {
	
	# if the species is already present, push the new gene_id to the annonymous array
	push ( @{ ${ $orthoDB{$orthogroup} }{$species} }, $prot_id );

    }
}


#print Dumper \%orthoDB;


# check single-copy genes for each orthogroup
my @single_copy;
my @single_copy_in_ref;
my $orthogroup_count = 0;

foreach my $group ( keys %orthoDB ) {

# WE ARE NOW INSIDE EACH ORTHOGROUP

    # empty list of single copy genes for each orthogroup
    @single_copy = ();
    @single_copy_in_ref = ();

    foreach my $sp ( keys %{ $orthoDB{$group} } ) {

	# filter species have a single copy
	if ( scalar @{ ${ $orthoDB{$group} }{$sp} } == 1 ) {

	    #save species with single gene into an array
	    push ( @single_copy, $sp );

	}

    }

    # push into new array the species that have a single copy && are present in the reference set
    foreach my $s ( @single_copy ) {
	if ( exists $ref_species{$s} ) { 
	    
	    push ( @single_copy_in_ref, $s );

	}
    }

    # get number of genes that are single copy for current orthogroup
    my $single_copy_in_ref_length = scalar @single_copy_in_ref;

    # if current orthogroup has copies for all (or e.g. 80%) of the species in the reference set
    # (i.e. # elems in @single_copy_in_ref is >= than a given % of species present in the reference set)
    # comment out one option or the other to require all species or only a % of them to be present

    # SINGLE COPY FOR ALL SPECIES IN THE REF SET
    #if (  $single_copy_in_ref_length == $ref_species_length ) {

    # SINGLE COPY FOR AT LEAST >70% OF THE SPECIES IN THE REF SET
    if (  $single_copy_in_ref_length > ($ref_species_length * 0.7)  ) {

	# count number of groups
	$orthogroup_count++;

	open (OUT, ">", "$group.fa");

	# print out ortholog group, species and ensembl ID
	foreach my $spec ( @single_copy_in_ref ) {

            foreach my $elem ( @{ ${ $orthoDB{$group} }{$spec} } ) {

		my $seq=$db->seq($elem);

		print OUT ">$spec\n";
		print OUT $seq,"\n";
	
	    }
	}
    }

    close(OUT);

}

print STDERR "\n$orthogroup_count ortholog groups found!\n";
print STDERR "\ndone!\n\n";
