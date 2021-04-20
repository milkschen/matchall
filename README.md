# Allele match

This is the Pysam implementation of the variant allele matching algorithm described in this [preprint](https://doi.org/10.1101/2021.01.06.425550).

A variant can be represented in multiple formats, making annotating not straightforward.
This algorithm constructs local haplotype around a variant and its nearby cohort variants and then perform exact matching.
Thereby, variants can be annotated accurately regardless of representation.

Currently, we only support annotating the `AF` (population allele frequency) field, but it's not difficult to support other tags as long as they are provided in the cohort call set.

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
```
ls cohort-chr*.release_missing2ref.no_calls.vcf.gz | sort -V > cohort.list
bcftools concat -f cohort.list -O z -o cohort.release_missing2ref.no_calls.vcf.gz; tabix cohort.release_missing2ref.no_calls.vcf.gz
python src/allele_match.py -p cohort.release_missing2ref.no_calls.vcf.gz -v <vcf> -r <grch38.fa> -o - | bgzip > <af.vcf>; tabix <af.vcf>
```
