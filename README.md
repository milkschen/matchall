# Allele match

This is a Pysam implementation of the variant allele matching algorithm described in this [preprint](https://doi.org/10.1101/2021.01.06.425550).

A variant can be represented in multiple formats, making annotating not straightforward.
This algorithm constructs local haplotype around a variant and its nearby cohort variants and then perform exact matching.
Thereby, variants can be annotated accurately regardless of representation.

Currently, we only support annotating the `AF` (population allele frequency) field, but it's not difficult to support other tags as long as they are provided in the cohort call set.

## Dependencies
- pysam (0.15.3)
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

## Pipeline
We first concatenate cohort VCFs for each contig as a unified one. If you're cohort call set is already unified, you can skip this step:
```
ls cohort-chr*.release_missing2ref.no_calls.vcf.gz | sort -V > cohort.list
bcftools concat -f cohort.list -O z -o cohort.release_missing2ref.no_calls.vcf.gz; tabix cohort.release_missing2ref.no_calls.vcf.gz
```

```
VCF=<vcf> # VCF to be annotated
REF=<grch38.fa> # reference FASTA for VCF
AF_CUTOFF="0.05" # allele frequency cutoff; we'll only a VCF with AF > `AF_CUTOFF` and a VCF with AF < `AF_CUTOFF`
VCF_AF=<annotated.vcf> # path to the output annotated VCF
VCF_COMMON=<annotated.common.vcf> # path to the output annotated common VCF
VCF_RARE=<annotated.rare.vcf> # path to the output annotated rare VCF

python src/allele_match.py -p cohort.release_missing2ref.no_calls.vcf.gz -v ${VCF} -r ${REF} -o - | bgzip > ${VCF_AF}; tabix ${VCF_AF}
bcftools view -O z -i 'AF>"${AF_CUTOFF}"' -o ${VCF_COMMON} ${VCF_AF}; tabix ${VCF_COMMON}
bcftools view -O z -i 'AF<="${AF_CUTOFF}"' -o ${VCF_RARE} ${VCF_AF}; tabix ${VCF_RARE}
```

## Test
```
python src/test-allele_match.py
```
