# Allele match

This is the Pysam implementation of the variant allele matching algorithm described in this [preprint](https://doi.org/10.1101/2021.01.06.425550).

A variant can be represented in multiple formats, making annotating not straightforward.
This algorithm constructs local haplotype around a variant and its nearby cohort variants and then perform exact matching.
Thereby, variants can be annotated accurately regardless of representation.

## Pipeline
```
ls cohort-chr*.release_missing2ref.no_calls.vcf.gz | sort -V > cohort.list
bcftools concat -f cohort.list -O z -o cohort.release_missing2ref.no_calls.vcf.gz; tabix cohort.release_missing2ref.no_calls.vcf.gz
```
