# em1
A 1-locus EM implementation to estimate haplotype frequencies.

usage: em1.groovy <input file> <frequencies output> <individual output>
  e.g., em1.groovy sout_interps.txt output/em-haps.txt output/em-ids.txt

The input file is a tab-separted 3-column file. Each line represents an  individual's haplotype list:
 column 1: an individual ID
 column 2: not currently used; reserved for allelic information
 column 3: a haplotype list (e.g., 1+1|2+3); all ambiguity is in the '|'

The frequencies output is a tab-separated 2-column file. Each line contains the haplotype name and its estimated frequency.

The individual output is a tab-separated 3-column file. Each line contains the individual ID, its estimated haplotype-pair list (which may be ambiguous), and a single haplotype pair, which is the same as column 2 if the haplotype list is unambiguous or a randomly chosen haplotype pair from column 2 if the haplotype list is ambiguous.
