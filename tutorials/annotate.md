# Matchall - annotate

### Usage
```
python src/annotate.py -q cohort.release_missing2ref.no_calls.vcf.gz -v ${VCF} -r ${REF} -o ${VCF_AF} --info {INFO}
```

**Example using a test case:**
```
python src/annotate.py -r test_data/chr20_1_580000.fa -v test_data/HG00733-hifi_deepvariant-chr20_568936_571052.vcf.gz -q test_data/chr20_560000_580000.cohort.vcf.gz --info AF
```

### Cohort variants
A DeepVariant-GLnexus-based call set for the 1000 Genomes Project (2504 samples) is avaialable [here](https://console.cloud.google.com/storage/browser/brain-genomics-public/research/cohort/1KGP/cohort_dv_glnexus_opt/v3_missing2ref;tab=objects?prefix=&forceOnObjectsSortingFiltering=false).
You can download it through the website, or use the following command:

```
for i in $(seq 1 22); do wget https://storage.googleapis.com/brain-genomics-public/research/cohort/1KGP/cohort_dv_glnexus_opt/v3_missing2ref/cohort-chr${i}.release_missing2ref.no_calls.vcf.gz; done
for i in $(seq 1 22); do wget https://storage.googleapis.com/brain-genomics-public/research/cohort/1KGP/cohort_dv_glnexus_opt/v3_missing2ref/cohort-chr${i}.release_missing2ref.no_calls.vcf.gz.tbi; done
```

You may also use other reference panels, such as the [GRCh38-based 1000 Genomes calls](https://www.internationalgenome.org/announcements/Variant-calls-from-1000-Genomes-Project-data-on-the-GRCh38-reference-assemlby/), [gnomAD](https://gnomad.broadinstitute.org/downloads), etc. 
We require a VCF format where allele frequency information is provided as an `AF` tag in the `INFO` field.


We first concatenate cohort VCFs for each contig as a unified one. If you're cohort call set is already unified, you can skip this step:
```
ls cohort-chr*.release_missing2ref.no_calls.vcf.gz | sort -V > cohort.list
bcftools concat -f cohort.list -O z -o cohort.release_missing2ref.no_calls.vcf.gz; tabix cohort.release_missing2ref.no_calls.vcf.gz
```

Annotate the population allele frequency tags (`AF`) for a VCF file using the reference panel:
```
REF=<grch38.fa> # reference FASTA for VCF
VCF=<vcf.gz> # VCF to be annotated
VCF_AF=<annotated.vcf.gz> # path to the output annotated VCF

python src/annotate.py -q cohort.release_missing2ref.no_calls.vcf.gz -v ${VCF} -r ${REF} -o ${VCF_AF}
tabix ${VCF_AF}
```

After annotating, we can split the VCF file based on an allele frequency cutoff:
```
AF_CUTOFF=0.05 # allele frequency cutoff; we'll generate a VCF with AF > `AF_CUTOFF` and a VCF with AF <= `AF_CUTOFF`
VCF_COMMON=<annotated.common.vcf> # path to the output annotated common VCF
VCF_RARE=<annotated.rare.vcf> # path to the output annotated rare VCF

bcftools view -O z -i "AF>${AF_CUTOFF}" -o ${VCF_COMMON} ${VCF_AF}
bcftools view -O z -i "AF<=${AF_CUTOFF}" -o ${VCF_RARE} ${VCF_AF}
# index (optional)
tabix ${VCF_COMMON}
tabix ${VCF_RARE}
```