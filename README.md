# Clivage-site-prediction

Each entry in these files consists in exactly three lines. This is an example of a protein description:

55 2SS8 HELAN 25 ALBUMIN 8 PRECURSOR (METHIONINE-RICH 2S PROTEIN) (SFA8).
MARFSIVFAAAGVLLLVAMAPVSEASTTTIITTIIEENPYGRGRTESGCYQQMEE
SSSSSSSSSSSSSSSSSSSSSSSSSCMMMMMMMMMMMMMMMMMMMMMMMMMMMMM

- The first line gives a description of the protein and its content may be ignored.
- The second line gives the beginning of the primary structure of the protein, with letters encoding the amino acids. Note that the rst amino acid is almost always M
(methionine), as it corresponds to the start codon.
- The third line is an annotation corresponding to the localization of the cleavage site. For each amino acid on the previous line, a letter indicates whether this amino acid belongs to the signal peptide (annotated with a letter S) or to the mature protein (annotated with a letter M). The first amino acid of the mature protein is annotated with a letter C.

Over a set of N sequences, with a known cleavage site for each, we first count $c(a,i)$ the number of occurrences of each amino acid $a \in A$, at every position $i \in [-p; ...; q-1] $, relative to the corresponding cleavage site. Then, for each a and i, let define $f(a, i) = c(a, i)/N$, the observed frequency of amino acid a at the relative position $i$. 

In a same way, by counting over the whole length of given sequences, one can compute the observed general background frequency $g(a$) of amino acid a in the given set, regardless of the position $i$. However, it must be noticed that the very first amino acid at the beginning of a sequence is almost always an $M$, because it corresponds to the transcription of the start codon. Also, one will not count letters on this first position to avoid a bias. 

These frequencies will be used as estimated probabilities to compute the probability of a given word to be located at a cleavage site, under an independent model. We rather use the logarithm of probabilities to go on additive calculations. 

Then, $\forall a \in A, \forall i \in [-p; ...; q-1]$ we define $s(a, i) = log(f(a, i)) - log(g(a))$. Also, as zero counts may occur, we use additive smoothing. Finally, for any word $w = a_{0}a_{1}...a_{p+q-1}$ the $q-1$ score defined as $\sum\limits_{i=-p}^{q-1} s(a_{p+i},i)$ may tell whether $w$ is the neighborhood of a cleavage site or not.

Then we tune a threshold to define a binary classifier.
