Statistical Analysis of the Codon-Usage in _Sulfolobus solfataricus_ and _E. Coli_
==================================================================================

by Michel Heiniger and Patrick Pfeifer

Standard Methods
----------------

Various research papers published during the past three decades are focusing on this particular problem.
One paper by Frank Wright that seems to be widely cited and was published on the 1st of March 1990 in
Volume 87, Issue 1, Pages 23–29 of Elsevier's "Gene" Journal is titled "The ‘effective number of codons’
used in a gene". It is available at <http://dx.doi.org/10.1016/0378-1119(90)90491-9> (and, although illegally,
at <http://crocodoc.com/2N8JTIc> -- choose "Download original" there).

### Student's $t$-test

### $\chi^2$-test

### The "effective number of codons", $N_c$, as defined by Wright

"The codon usage table of a gene can be subdivided
according to the number of synonymous codons belonging
to each aa. Thus for a gene using the 'universal' code, there
are 2 aa with only one codon choice, 9 with two, 1 with
three, 5 with four, and 3 with six. These represent five SF
types, designated SF types 1, 2, 3, 4, and 6 according to
their respective number of synonymous codons."
(Wright, 1990) (aa = amino-acid)

Then, for each aa, calculate:

\begin{align}\hat{F} &= \frac{n \displaystyle\sum_{i = 1}^{k}{(p_i^2 - 1)} }{n - 1}\end{align}

And for a complete gene, calculate:

\begin{align}\hat{N}_c &= 2 + 9 / \overline{\hat{F}_2} + 1 / \hat{F}_3 + 5 / \overline{\hat{F}_4} + 3 / \overline{\hat{F}_6} \end{align}

"where $\overline{\hat{F_i}}$, is the average homozygosity estimate for SF type
i, and $F_i$ for each aa are calculated using Eqn. (1)." (Wright, 1990)

Our Methods
-----------

We will (a) carry out $t$- and $\chi^2$-tests for the usage of each codon in respect to the possible codons for Wright's SF2, SF3, SF4 and SF6 families of aa and (b) calculate $\hat{N}_c$ for the whole genome.

### Student's $t$-test

#### Null-Hypothesis ($H_0$):

We calculate the hypotetical probability of the codon usage for the amino acids of
all genes by assuming the following null hypothesis:

Assuming there are $n$ codons coding for a given amino acid $AS$,
then the probability for any of those codons $C_x$ to be actually used, would be:

$P(C_x, AS) = \frac{1}{n}$

The opposite Probability, that the codon $C_x$ would *not* be used, would be:

$P(\overline{C_x}, AS) = 1 - \frac{1}{n}$

#### Test

We examine the actual codon usage in both organisms. For a given amino-acid,
e.g. Glycin: CGA, GGC, GGG, GGT, we examine if codon-usage satisfies the
null-hypothesis.

e.g. $H_0$: $P(C_x, Glycin) = \frac{1}{4}$; $P(\overline{C_x}, Glycin) = \frac{3}{4}$

Assuming the null-hypothesis, we can now calculate the probability that out
of $m$ identical amino-acids $n$ _or more_ are coded
