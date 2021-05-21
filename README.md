# Matchall

`matchall` (MATCH ALLele) annotates a VCF with the information in another VCF. 
This is a Pysam implementation of the variant allele matching algorithm described in this [preprint](https://doi.org/10.1101/2021.01.06.425550).

Example use case: 

- filling the population allele frequency (`AF`) field in a variant set using information in a reference panel.
- comparing two VCFs to find shared and private variants

A variant can be represented in multiple formats. An example in the table below shows a variant in two forms. The ambiguity in variant representation can confound annotating and result in errors.

| POS | REF | ALT |
|-----|-----|-----|
| 1   | GAC | GA  |
| 2   | AC  | A   |

The `matchall` algorithm solves this issue by comparing a variant and a set of queried variants from another VCF using re-constructed local haplotypes.
This algorithm annotates variants accurately regardless of representation.


## Installation
### Dependencies - required
- Python (3.6+)
- pysam (0.15.3)

### Dependencies. - optional (used in our data processing pipelines)
- bcftools (1.12)
- tabix (1.12)

### Download matchall
```
https://github.com/milkschen/matchall.git
```

## Usage
### Annotate
```
python annotate.py -v target.vcf.gz -q query.vcf.gz -r ref.fa -o out.vcf.gz
```
- [Detailed tutorial for matchall-annotate](tutorials/annotate.md)

### Compare
```
python src/compare.py -v A.vcf.gz -q B.vcf.gz -op A_0-B_1 -m annotate,private,isec -o A_0-B_1.vcf.gz -r ref.fa
```
- [Detailed tutorial for matchall-compare](tutorials/compare.md)

## Test
```
python src/test_matchall.py
python src/test_end_to_end.py
```
Or
```
sh test_all.sh
```
