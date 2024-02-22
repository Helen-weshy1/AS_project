# **Regulation of endocrine cell alternative splicing revealed by single-cell RNA sequencing in type 2 diabetes pathogenesis**


---
The prevalent RNA alternative splicing (AS) contributes to molecular diversity, which has been demonstrated in cellular function regulation and disease pathogenesis. However, the contribution of AS in pancreatic islets during diabetes progression remains unclear. Here, we reanalyzed the full-length single-cell RNA sequencing data from the deposited database to investigate AS regulation across human pancreatic endocrine cell types in non-diabetic (ND) and type 2 diabetic (T2D) individuals.

## **Table of contents**
---
1. Initialization
2. Data availability
3. Data preprocessing
    1. AS quantification
    2. Gene expression quantification
4. Downstream analysis
    1. AS profile
    2. Detect differential splicing events
    3. De novo motif analysis

## **Initialization**
---
Prerequisite programming languages and related software:

* Perl
* Python
* R

## **Data availability**
---
Two single-cell sequencing FASTQ data produced by The Jackson Laboratory and Regeneron Pharmaceuticals were downloaded from NCBI Sequence Read Archive (SRA) under accession numbers SRP075970 and SRP075377, respectively. Another scRNA-seq FASTQ data by Department of Genetics and Genome Sciences was downloaded from NCBI SRA under accession number GSE101207.

## **Data preprocessing**
---
### **AS quantification**
#### **Quantas**
Islet cell RNA-seq data were firstly mapped by OLego (v1.1.2) to the reference genome (hg19).

We used the Quantas pipeline (http://zhanglab.c2b2.columbia.edu/index.php/Quantas) to quantify AS based on the number of exon junction reads.The level of inclusion of alternative exons was represented by percent spliced in (PSI). We required exons to have junction read coverage ≥20 when estimated PSI.

``` 
perl /summarize_splicing_wrapper.pl -c ./cache -v -big -conf /Quantas/index_annot/hg19/annotation/hg19.conf -dbkey hg19 -cass test.uniq.bed /count_AS/test
``` 
Then created the PSI matrix with R, which code saved in ./QUANTAS/PSI quantification.R. Seurat was used to analyze PSI matrix in ./QUANTAS/AS seuratobject.R.

#### **MARVEL**
For MARVEL pipline, fastq data was mapped by STAR (v2.7.10b) to the reference genome (hg19) and used rMATS to identify SE, MXE, RI, A5SS and A3SS splicing events. 

``` 
#STAR mapping
STAR --runThreadN 10 --runMode alignReads --readFilesCommand zcat --quantMode TranscriptomeSAM GeneCounts --outSAMattrIHstart 0 --outSAMstrandField intronMotif --outSAMattributes NH HI AS nM XS --twopassMode Basic --outSAMtype BAM Unsorted --outSAMunmapped None --genomeDir /STAR_index/hg19/ --readFilesIn $i --outFileNamePrefix /bam/SRR 
    
#rMATs annotation
rmats.py --b1 ./rMATS/ASanno/SRR/BAM_fls.txt --gtf /gencode.v38lift37.annotation.gtf --od ./rMATS/ASanno/SRR --tmp ./rMATS/ASanno/SRR/tmp -t single --readLength 75 --variable-read-length --nthread 10 --statoff

#samtools sort bam
samtools sort -@8 /bam/SRR.Aligned.out.bam -o /bam/SRR.sort.bam
    
#create intron_count
Rscript ./MARVEL/rscript_bedtools_input.R SRR
bedtools coverage -g /marvel/intron/SRR/hg19.chrom.sizes.txt -split -sorted -a /SRR/RI_Coordinates_sorted.bed -b /bam/SRR.sort.bam > /SRR/intron_count.txt -d
``` 

We first prepared input file in linux, which code was saved in ./MARVEL/pre_input.R. Then create the MARVEL object in ./MARVEL/marvel.R.

### **Gene expression quantification**
Gene expression quantification was quantified using the Quantas pipline.

``` 
perl /summarize_expression_wrapper.pl -big -c ../cache -exon /Quantas/index_annot/hg19/annotation/hg19.exon.uniq.core.bed -e2g /Quantas/index_annot/hg19/annotation/hg19.exon.uniq.core.id2gene2symbol -v test.uniq.bed /gene_exp/test
```

Then created the seurat object with ./QUANTAS/Gene seuratobject.R


## **Downstream analysis**
### **AS profile**
---
We analyzed splicing profile through Seurat pipline. To filter unquantifiable exons, we maintained exons (4,616 exons in Lawlor data and 2974 in Xin data passed this filtering) with junction read coverage ≥20 in ≥10% of 972 cells. The missing values in y matrix were replaced by an extract value far from all y. The details were described in  ./QUANTAS/Gene seuratobject.R.

### **Detect differential splicing events**
---
Differential splicing exons were performed using the Quantas pipeline.

``` 
perl /test_splicing_diff.pl -type cass -v --min-cov 20 --id2gene2symbol /Quantas/index_annot/hg19/annotation/Hs.seq.all.AS.chrom.can.id2gene2symbol test.conf test.diff.txt
```
we performed one vs. other comparisons between each splicing cluster and defined the exons with the following criteria as significant: junction read coverage ≥ 20, False discovery rate (FDR, BH correction) ≤ 0.05, exon quantifiable in ≥ 10% of cells in different groups and |PSI| ≥ 0.1. This part was mainly described in ./Figure 3/Fig3_code.R.


### **De novo motif analysis**
---
We used rMAPS2 to analyze de novo motif and infer RBP activity in differential splicing exons. The input files were prepared in ./Figure 5/Figure 5 prepare motif analysis.R and get the result, which saved as ./Figure 5/A Enriched_motif_predicted_RBP_&exp.csv. 








