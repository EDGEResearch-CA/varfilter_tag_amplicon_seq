# varfilter_tag_amplicon_seq

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15446855.svg)](https://doi.org/10.5281/zenodo.15446855)

Data processing and filtering steps specific to mutation calling in tagged amplicon sequencing data. 

Evaluating support metrics like amplicon-read counts and variant allele frequency (VAF) to determine if a variant is a real mutation.

# Script Summary
Custom mutation confidence scoring system. 

Application of thresholds to determine the confidence of each variant being a real mutation.

It considers formalin-induced artifacts, which can falsely look like mutations.

It assumes specific columns names - edit if needed

It contains two example data frames in the input folder

# Summary of custom rules
Minimum of 2 supporting amplicons

VAF > 2%

For 2 supporting amplicons: Read support for amplicon > 20, except when the alteration is C>T/G>A then it must be >50

For 3 supporting amplicons: Read support for amplicon > 3, except when the alteration is C>T/G>A then it must be >20
