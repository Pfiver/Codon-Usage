#!/usr/bin/python

from Bio.Data import CodonTable 
from Bio.GenBank import FeatureParser

def main():

	input_files = [ "EcoliK12.gb", "genome_sso.gb" ]

	standard_forward_table = CodonTable.standard_dna_table.forward_table
	counting_forward_table = LookupCountingForwardTable(standard_forward_table)
	
	codon_table = CodonTable.CodonTable(forward_table=counting_forward_table)

	for fil in input_files:

		record = FeatureParser().parse(open(fil))

		complement = record.seq.complement()

		counting_forward_table.reset_counts()

		for feat in record.features:

			if feat.type == "gene":

				bSeq = feat.strand == 1 and record.seq or complement

				start = feat.location.start.position
				end = feat.location.end.position

				pSeq = bSeq[start:end].translate(codon_table)

		print counting_forward_table.counts

class LookupCountingForwardTable(dict):

	counts = {}
	real_table = {}

	def __init__(self, table):
		self.real_table = table
		self.reset_counts()

	def __getitem__(self, key):
		self.counts[key] += 1
		return self.real_table[key]

	def reset_counts(self):
		self.counts = dict((key, 0) for key in self.real_table)


main()
