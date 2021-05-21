# Matchall - annotate

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