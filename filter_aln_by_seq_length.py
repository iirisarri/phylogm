#!/usr/bin/env python3

import sys 
from Bio import SeqIO
import re

'''

Iker Irisarri. LIB, Nov 2023

Filters alignments by % of missing data. Default filter requires unambiguous characters 
	for more than 30% of the total alignment length.

Defaults are for amino acids and 30% of the sequence length, but might be modified

Sequences are filtered individually by the amount of undetermined characters. Thus, for
	filtering alignments, sequences must be all of the same length!

usage: filter_aln_by_seq_length.py input.fasta > filtered.fasta

'''

input_file = sys.argv[1]
threshold = float("0.3") # 30%

for record in SeqIO.parse(input_file, "fasta"):

	length = len(record)
	total_AA = record.seq.count("A") + record.seq.count("R") + record.seq.count("N") + record.seq.count("D") + record.seq.count("C") + record.seq.count("E") + record.seq.count("Q") + record.seq.count("G") + record.seq.count("H") + record.seq.count("I") + record.seq.count("L") + record.seq.count("K") + record.seq.count("M") + record.seq.count("F") + record.seq.count("P") + record.seq.count("S") + record.seq.count("T") + record.seq.count("W") + record.seq.count("Y") + record.seq.count("V")
	# for nucleotides
	#total_AA = record.seq.count("A") + record.seq.count("C") + record.seq.count("G") + record.seq.count("T")
	#total_gaps = record.seq.count("-")
	#total_undet = record.seq.count("X ") + record.seq.count("?")
	#print(record.id, length, total_AA)

	# filter by % of unambiguous characters
	if ( total_AA > length * threshold ):
	
		#print(record.id, length, total_AA)
		print(">", record.id, sep="")
		print(record.seq)
		
	
