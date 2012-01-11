#!/usr/bin/python

from Bio.Data import CodonTable 
from Bio.GenBank import FeatureParser


input_files = [ "EcoliK12.gb", "genome_sso.gb" ]


def main():

	standard_forward_table = CodonTable.standard_dna_table.forward_table

	codons = standard_forward_table.keys()
	amino_acids = set(standard_forward_table.values())

	print "Translation table"
	print standard_forward_table
	print "Number of codons / amino_acids: %d / %d" % (len(codons), len(amino_acids))
	print

	for fil in input_files:

		counts = count_usage(fil)

		hz = {}

		# calculate
		for aa in amino_acids:
			hz[aa] = homozygosity(counts, aa)

		# display
		print "homozygosities:"
		for fam in (1, 2, 3, 4, 6):
			print "fam:", fam
			for aa in aminoacids(fam):
				print " %s: %f" % (aa, hz[aa])

		print
		print "absolute codon usage N_c:", absolute_codon_usage(hz)
		print


def family(aa):
	
	table = CodonTable.standard_dna_table.forward_table

	return len(filter(lambda it: it[1] == aa, table.iteritems()))


def codons(aa):

	table = CodonTable.standard_dna_table.forward_table

	return list(it[0] for it in table.iteritems() if it[1] == aa)


def aminoacids(fam):

	table = CodonTable.standard_dna_table.forward_table
	amino_acids = set(table.values())

	return filter(lambda aa: family(aa) == fam, amino_acids)


def homozygosity(counts, aa):

	fam = family(aa)

	if fam == 1:
		return 1

	sum_ = 0
	cdns = codons(aa)
	n = sum(counts[codon] for codon in cdns)

	for i in range(fam):
		p = float(counts[cdns[i]]) / n
		sum_ += p * p

	return (n * sum_ - 1) / (n - 1)


def absolute_codon_usage(hz):

	avg = {}
	for fam in (2, 3, 4, 6):
		hzs = list(hz[aa] for aa in aminoacids(fam))
		avg[fam] = sum(hzs) / len(hzs)

	return 2 + (9/avg[2]) + (1/avg[3]) + (5/avg[4]) + (3/avg[6])


def count_usage(input_file):

	standard_forward_table = CodonTable.standard_dna_table.forward_table
	counting_forward_table = LookupCountingForwardTable(standard_forward_table)
	codon_table = CodonTable.CodonTable(forward_table=counting_forward_table)

	record = FeatureParser().parse(open(input_file))

	complement = record.seq.complement()

	for feat in record.features:

		if feat.type == "gene":

			bSeq = feat.strand == 1 and record.seq or complement

			start = feat.location.start.position
			end = feat.location.end.position

			pSeq = bSeq[start:end].translate(codon_table)

#	print "Usage counts for ", input_file
#	print counting_forward_table.counts

	return counting_forward_table.counts


class LookupCountingForwardTable(dict):

	counts = {}
	real_table = {}

	def __init__(self, table):
		self.real_table = table
		self._init_counts()

	def __getitem__(self, key):
		self.counts[key] += 1
		return self.real_table[key]

	def _init_counts(self):
		self.counts = dict((key, 0) for key in self.real_table)


main()
