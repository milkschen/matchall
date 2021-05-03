# Allele match
 Ex
`allele_match` annotates a VCF with the information in another VCF. 
Example use case: filling the population allele frequency (`AF`) field in a called variant set using information in a reference panel, such as the 1000 Genomes Project.
This is a Pysam implementation of the variant allele matching algorithm described in this [preprint](https://doi.org/10.1101/2021.01.06.425550).

A variant can be represented in multiple formats. An example in the table below shows a variant in two forms. The ambiguity in variant representation can confound annotating and result in errors.

| POS | REF | ALT |
|-----|-----|-----|
| 1   | GAC | GA  |
| 2   | AC  | A   |

The `allele_match` algorithm solves this issue by comparing a variant and a set of queried variants from another VCF using re-constructed local haplotypes.
This algorithm annotates variants accurately regardless of representation.

Currently, we only support annotating the `AF` (population allele frequency) field, but it's supposedly capable of supporting other tags. Please file an issue or pull request if there's a need.

## Dependencies
### Required
- Python (3.6+)
- pysam (0.15.3)

### Optional (used in our data processing pipelines)
- bcftools (1.12)
- tabix (1.12)




## Cohort variants
A DeepVariant-GLnexus-based call set for the 1000 Genomes Project (2504 samples) is avaialable [here](https://console.cloud.google.com/storage/browser/brain-genomics-public/research/cohort/1KGP/cohort_dv_glnexus_opt/v3_missing2ref;tab=objects?prefix=&forceOnObjectsSortingFiltering=false).
You can download it through the website, or use the following command:

```
for i in $(seq 1 22); do wget https://storage.googleapis.com/brain-genomics-public/research/cohort/1KGP/cohort_dv_glnexus_opt/v3_missing2ref/cohort-chr${i}.release_missing2ref.no_calls.vcf.gz; done
for i in $(seq 1 22); do wget https://storage.googleapis.com/brain-genomics-public/research/cohort/1KGP/cohort_dv_glnexus_opt/v3_missing2ref/cohort-chr${i}.release_missing2ref.no_calls.vcf.gz.tbi; done
```

You may also use other reference panels, such as the GRCh38-based 1000 Genomes calls, gnomAD, etc. 
We require a VCF format where allele frequency information is provided as an `AF` tag in the `INFO` field.

## Usage
We first concatenate cohort VCFs for each contig as a unified one. If you're cohort call set is already unified, you can skip this step:
```
ls cohort-chr*.release_missing2ref.no_calls.vcf.gz | sort -V > cohort.list
bcftools concat -f cohort.list -O z -o cohort.release_missing2ref.no_calls.vcf.gz; tabix cohort.release_missing2ref.no_calls.vcf.gz
```

Annotate the population allele frequency tags (`AF`) for a VCF file using the reference panel:
```
REF=<grch38.fa> # reference FASTA for VCF
VCF=<vcf> # VCF to be annotated
VCF_AF=<annotated.vcf> # path to the output annotated VCF

python src/allele_match.py -p cohort.release_missing2ref.no_calls.vcf.gz -v ${VCF} -r ${REF} -o - | bgzip > ${VCF_AF}; tabix ${VCF_AF}
```

After annotating, we can split the VCF file based on an allele frequency cutoff:
```
AF_CUTOFF=0.05 # allele frequency cutoff; we'll generate a VCF with AF > `AF_CUTOFF` and a VCF with AF <= `AF_CUTOFF`
VCF_COMMON=<annotated.common.vcf> # path to the output annotated common VCF
VCF_RARE=<annotated.rare.vcf> # path to the output annotated rare VCF

bcftools view -O z -i "AF>${AF_CUTOFF}" -o ${VCF_COMMON} ${VCF_AF}; tabix ${VCF_COMMON}
bcftools view -O z -i "AF<=${AF_CUTOFF}" -o ${VCF_RARE} ${VCF_AF}; tabix ${VCF_RARE}
```

## Test
```
python src/test-allele_match.py
```
