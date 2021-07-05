# Matchall - compare

`matchall-compare` uses the haplotype-based allele matching algorithm to support accurate allele/genotype matching.

#### Usage:

Allele-resolution:
```
python src/compare.py -v A.vcf.gz -q B.vcf.gz -op out-A_0-B_1 -m annotate,private,isec -o out-A_0-B_1.vcf.gz -r ref.fa
```

**Example using a test case:**
```
python src/compare.py -r test_data/chr20_1_580000.fa -v test_data/chr20_156877-A.vcf.gz -q test_data/chr20_156877-B.vcf.gz -m isec,private -op test_data/test-compare-chr20_156877
```

Genotype-resolution:
```
python src/compare.py -v A.vcf.gz -q B.vcf.gz -op out-A_0-B_1 -m annotate,private,isec -o out-A_0-B_1.vcf.gz -r ref.fa -gt
```

**Example using a test case:**
```
python src/compare.py -r test_data/chr20_1_580000.fa -v test_data/chr20_156877-A.vcf.gz -q test_data/chr20_156877-B.vcf.gz -m isec,private -gt -op test_data/test-compare-chr20_156877-gt
```

#### Difference between allele- and genotype-resolution:
- allele-resolution: for a HET variant where one allele matches and the other doesn't, split the genotype. Report the matched allele in `isec` and the unmatched allele in `private`. 
- genotype-resolution: report a match genotype when all alleles are matched.


### Arguments
```
usage: compare.py [-h] -v VCF -q QUERY_VCF -r REF [-o OUT] [-op OUT_PREFIX]
                  [--padding PADDING] [-m MODE] [-gt] [--happy] [--debug]

optional arguments:
  -h, --help            show this help message and exit
  -v VCF, --vcf VCF     Path to target VCF. [required]
  -q QUERY_VCF, --query-vcf QUERY_VCF
                        Path to query VCF (TBI or CSI indexes are required).
                        [required]
  -r REF, --ref REF     Path to the reference FASTA (FAI index is required).
                        [required]
  -o OUT, --out OUT     Path to output VCF. Set to "-" to print to stdout.
                        ["-"]
  -op OUT_PREFIX, --out-prefix OUT_PREFIX
                        Prefix to output VCF for "isec" and "private" modes.
                        [None]
  --padding PADDING     Size of paddings when performing allele matching. [10]
  -m MODE, --mode MODE  Mode: ["annotate", "isec", "private"]. Multiple modes
                        can be toggled at once, e.g. "annotate,private,isec".
                        ["annotate"]
  -gt, --gt-resolution  Set to use genotype-level resolution (less precise).
                        [False]
  --happy               Set for hap.py VCFs. Will only consider variants in
                        confident regions. [False]
  --debug               Set to print debug messages. [False]
```


### Examples

We show examples that can be correctly matched by `matchall-compare`, but often missed using typical variant matching tools.

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


### Allelic-resolution matching
```
A.vcf
(1)	chr20	156877	.	GTA	G	.	PASS	.	GT	1/1:20:26:0,25:0.961538:38,20,0
B.vcf
(2)	chr20	156877	.	GTA	G,*	.	PASS	GT	2|1
(3)	chr20	156877	.	GTATA	G	.	PASS	GT	1|.
```

In the `allele-resolution` mode (default), `matchall` splits the alleles in B and reports 156877_GTA_G as matched and 156877_GTATA_G as unique to B.
In this case both VCFs agree with the 156877_GTA_G deletion, but VCF B additionally calls a deletion (156877_GTATA_G).
In the `gt-resolution` mode, `matchall` reports this as private to both VCFs.

