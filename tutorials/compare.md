# Matchall - compare

Matchall supports accurate 


### Complex matching results
```
A.vcf
(1)	chr20	156877	.	GTA	G	38	PASS	.	GT:GQ:DP:AD:VAF:PL	1/1:20:26:0,25:0.961538:38,20,0
B.vcf
(2)	chr20	156877	.	GTA	G,*	427.84	PASS	AC=1,1;AN=2;DP=26;MQ=40;MQ0=0;NS=1	GT:GQ:DP:MQ:PS:PQ:FT	2|1:210:26:40:156813:99:PASS
(3)	chr20	156877	.	GTATA	G	190.31	PASS	AC=1;AN=1;DP=26;MQ=40;MQ0=0;NS=1	GT:GQ:DP:MQ:PS:PQ:FT	1|.:190:26:40:156813:99:PASS
```

In this case both VCFs agree with the 156877_GTA_G deletion, but VCF B additionally calls a deletion (156877_GTATA_G).
`bcftools isec` considers the genotypes and reports this as private to both VCFs.
`matchall` further splits the alleles in B and reports 156877_GTA_G as matched and 156877_GTATA_G as unique to B.
