
BEDtools intersect analysis
```
### Pol III tRNA
bedtools intersect -s -F 0.25 -wao -a ARPE19_UNINF_24h-4_dorado.0.6.0_bwa_W13k6T20.HG38.merged.primary.sorted.bed -b ../GTF/hg38-tRNAs.bed | awk '$7 ~ /chr/' > intersect/ARPE19_UNINF_24h-4_dorado.0.6.0_bwa_W13k6T20.HG38.merged.primary.sorted.hg38-tRNAs.intersect.bed

### Pol III primary only
bedtools intersect -s -F 0.25 -wao -a ARPE19_UNINF_24h-4_dorado.0.6.0_bwa_W13k6T20.HG38.merged.primary.sorted.bed -b ../GTF/gencode.v45.PolIIItx.primaryOnly.bed | awk '$7 ~ /chr/' > intersect/ARPE19_UNINF_24h-4_dorado.0.6.0_bwa_W13k6T20.HG38.merged.primary.sorted.PolIIItx.primaryOnly.intersect.bed

### Pol III psuedogene only
bedtools intersect -s -F 0.25 -wao -a ARPE19_UNINF_24h-4_dorado.0.6.0_bwa_W13k6T20.HG38.merged.primary.sorted.bed -b ../GTF/gencode.v45.PolIIItx.pseudogeneOnly.bed | awk '$7 ~ /chr/' > intersect/ARPE19_UNINF_24h-4_dorado.0.6.0_bwa_W13k6T20.HG38.merged.primary.sorted.PolIIItx.pseudoOnly.intersect.bed

### Pol III SNAR only
bedtools intersect -s -F 0.25 -wao -a ARPE19_UNINF_24h-4_dorado.0.6.0_bwa_W13k6T20.HG38.merged.primary.sorted.bed -b ../GTF/SNARs.bed | awk '$7 ~ /chr/' > intersect/ARPE19_UNINF_24h-4_dorado.0.6.0_bwa_W13k6T20.HG38.merged.primary.sorted.SNARs.intersect.bed

### no Pol III 
bedtools intersect -s -F 0.25 -wao -a ARPE19_UNINF_24h-4_dorado.0.6.0_bwa_W13k6T20.HG38.merged.primary.sorted.bed -b ../GTF/gencode.v45.NoPolIIItx.bed | awk '$7 ~ /chr/' > intersect/ARPE19_UNINF_24h-4_dorado.0.6.0_bwa_W13k6T20.HG38.merged.primary.sorted.NoPolIIItx.intersect.bed

### Full Anno
bedtools intersect -s -F 0.25 -wao -a ARPE19_UNINF_24h-4_dorado.0.6.0_bwa_W13k6T20.HG38.merged.primary.sorted.bed -b ../GTF/gencode.v45.fullAnnotation_inc_PolIIItx_tRNAs_SNARs.bed | awk '$7 ~ /chr/' > intersect/ARPE19_UNINF_24h-4_dorado.0.6.0_bwa_W13k6T20.HG38.merged.primary.sorted.fullAnnotation_inc_PolIIItx_tRNAs_SNARs.intersect.bed

### no overlap at all
bedtools intersect -s -v -a ARPE19_UNINF_24h-4_dorado.0.6.0_bwa_W13k6T20.HG38.merged.primary.sorted.bed -b ../GTF/gencode.v45.fullAnnotation_inc_PolIIItx_tRNAs_SNARs.bed > intersect/ARPE19_UNINF_24h-4_dorado.0.6.0_bwa_W13k6T20.HG38.merged.primary.sorted.NoOverlaps.intersect.bed
```

Processing of intersect bed outputs into count files
```
cut -f10 intersect/ARPE19_UNINF_24h-4_dorado.0.6.0_bwa_W13k6T20.HG38.merged.primary.sorted.hg38-tRNAs.intersect.bed | sort | uniq -c | sed 's/^ *//g' | sed 's/ /\t/g' > count/ARPE19-UNINF-24h-4_bwa_W13k6T20.SR.primary-merged.HG38.sorted.hg38-tRNAs.count.txt

cut -f10 intersect/ARPE19_UNINF_24h-4_dorado.0.6.0_bwa_W13k6T20.HG38.merged.primary.sorted.PolIIItx.primaryOnly.intersect.bed | sort | uniq -c | sed 's/^ *//g' | sed 's/ /\t/g' > count/ARPE19-UNINF-24h-4_bwa_W13k6T20.SR.primary-merged.HG38.sorted.PolIIItx.primaryOnly.count.txt

cut -f10 intersect/ARPE19_UNINF_24h-4_dorado.0.6.0_bwa_W13k6T20.HG38.merged.primary.sorted.PolIIItx.pseudoOnly.intersect.bed | sort | uniq -c | sed 's/^ *//g' | sed 's/ /\t/g' > count/ARPE19_UNINF_24h-4_dorado.0.6.0_bwa_W13k6T20.HG38.merged.primary.sorted.PolIIItx.pseudoOnly.count.txt

cut -f10 intersect/ARPE19_UNINF_24h-4_dorado.0.6.0_bwa_W13k6T20.HG38.merged.primary.sorted.SNARs.intersect.bed | sort | uniq -c | sed 's/^ *//g' | sed 's/ /\t/g' > count/ARPE19_UNINF_24h-4_dorado.0.6.0_bwa_W13k6T20.HG38.merged.primary.sorted.SNARs.count.txt

cut -f10 ARPE19_UNINF_24h-4-SR_bwa_W13k6T20.merged.HG38.primary.sorted.noPolIIItx.intersect.bed | sort | uniq -c | sed 's/^ *//g' | sed 's/ /\t/g' > count/ARPE19_UNINF_24h-4-SR_bwa_W13k6T20.merged.HG38.primary.sorted.noPolIIItx.counts.txt






```
