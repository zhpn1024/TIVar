# TIVar
Translation Initiation Variation

Predict translation initiation (TI) efficiency for potential start codons, based on the context sequence near the start codon. Given SNP/Indel variation, this tools can predict changes of TI efficiencies between ref and alt alleles.

# INSTALL

Python version >= 3.4.

**Requirements**

[NumPy](https://numpy.org/)

[PyTorch](https://pytorch.org/)

**Install from source**

`git clone https://github.com/zhpn1024/TIVar`

`python setup.py install`

or

`python setup.py install --user`

# Usage

**predict**

This module can and calculate TI efficiency scores from given sequences.

Fasta sequence file as input:

`tivar predict -S test1.fa -o out1.txt`

Provide sequence in the parameter:

`tivar predict -s aaaaaacaaaaaaaTGTACAATGGATGCATTGAAATTATATGTAATTGTATAAATGGTGCAACA -o out1.txt`

Provide transcript annotation and genome sequence:

`tivar predict -g hg38_gc31.gtf.gz -f hg38.fa -o out1.txt`

The output is like:

|SeqID|Pos|StartSeq|EffScore|
|-----|-----|-----|-----|
|Seq|13|aacaaaaaa-aTG-TACA|0.08535|
|Seq|20|aaaTGTACA-ATG-GATG|0.34153|


**diff**

This module predict TI changes caused by sequence variation.

`tivar diff -i test.vcf -g hg38_gc31.gtf.gz -f hg38.fa -o out2.txt`

The output is like:

|Gid|Tid|Var|GenoPos|Strand|Pos|RefSeq|AltSeq|EffeRef|EffeAlt|Diff|FC|Type|
|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
|ENSG00000134262.13|ENST00000369569.6|chr1:113895309:A>AC|113895310|-|2056|ACCCTCCAG-ATG-GCTC|ACCCTCCAG-AGT-GGCT|0.31412|0.0|-0.3141|0.0|TI_decreased|
|ENSG00000134262.13|ENST00000369569.6|chr1:113895309:A>AC|113895310|-|2056|ACCCTCCAG-ATG-GCTC|CCCTCCAGA-GTG-GCTC|0.31412|0.04335|-0.2708|0.138|TI_decreased|

