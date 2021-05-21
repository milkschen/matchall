# Matchall - compare

`matchall-compare` uses the haplotype-based allele matching algorithm to support accurate allele/genotype matching.

### Arguments
- `-v`: VCF_0
- `-q`: VCF_1
- `-r`: reference FASTA
- `-op`: output prefix for `isec` and `private` results
- `-o`: annotated output VCF
- `-gt`: set to evaluate in genotype resolution

*Usage:*

Allele-resolution:
```
python src/compare.py -v A.vcf.gz -q B.vcf.gz -op out-A_0-B_1 -m annotate,private,isec -o out-A_0-B_1.vcf.gz -r ref.fa
```

Genotype-resolution:
```
python src/compare.py -v A.vcf.gz -q B.vcf.gz -op out-A_0-B_1 -m annotate,private,isec -o out-A_0-B_1.vcf.gz -r ref.fa -gt
```

*Difference between allele- and genotype-resolution:*
- allele-resolution: for a HET variant where one allele matches and the other doesn't, split the genotype. Report the matched allele in `isec` and the unmatched allele in `private`. 
- genotype-resolution: report a match genotype when all alleles are matched.

### Examples

The examples in this section are not handled by `bcftools isec`.

#### Alternative VCF representation
An insertion is represented in different forms:
```
A.vcf
chr20	590301	.	AAC	AACAC	.	PASS	GT	1/0
B.vcf
chr20	590301	.	A	AAC	.	PASS		GT	1|0
```

#### Multi-allelic results
VCF B uses two records to represent one genotype and VCF A uses one:
```
A.vcf
chr20	201932	.	CA	C,CAAAAA	.	PASS	.	GT	1/2
B.vcf
chr20	201932	.	C	CAAAA	.	PASS	GT	1|0
chr20	201932	.	CA	C	.	PASS	GT	0|1
```

```
A.vcf
chr20	201932	.	CA	C,CAAAAA	.	PASS	MATCH=1,1	GT	1/2

B.vcf
chr20	201932	.	C	CAAAA	.	PASS	GT	1|0
chr20	201932	.	CA	C	.	PASS	GT	0|1
```


### Allelic-resolution matching
```
A.vcf
(1)	chr20	156877	.	GTA	G	.	PASS	.	GT	1/1:20:26:0,25:0.961538:38,20,0
B.vcf
(2)	chr20	156877	.	GTA	G,*	.	PASS	GT	2|1
(3)	chr20	156877	.	GTATA	G	.	PASS	GT	1|.
```

In this case both VCFs agree with the 156877_GTA_G deletion, but VCF B additionally calls a deletion (156877_GTATA_G).
`bcftools isec` considers the genotypes and reports this as private to both VCFs.
`matchall` further splits the alleles in B and reports 156877_GTA_G as matched and 156877_GTATA_G as unique to B.
