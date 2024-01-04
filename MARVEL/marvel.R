library(devtools)
install_github("wenweixiong/MARVEL")

# Load MARVEL package
library(MARVEL)

# Load adjunct packages for selected MARVEL features
# General data processing, plotting
library(ggnewscale)
library(ggrepel)
library(parallel)
library(reshape2)
library(stringr)
library(textclean)

# Dimension reduction analysis
library(factoextra)
library(FactoMineR)

# Modality analysis
library(fitdistrplus)

# Differential splicing analysis
library(kSamples)
library(twosamples)

# Gene ontology analysis
library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)

# Nonsense-mediated decay (NMD) analysis
library(Biostrings)
library(BSgenome)
#BiocManager::install('BSgenome.Hsapiens.NCBI.GRCh38')
library(BSgenome.Hsapiens.NCBI.GRCh38)

# Load adjunct packages for this tutorial
library(data.table)
library(ggplot2)
library(gridExtra)


######################### Create MARVEL object #################################
sample_meta=read.table('/AS/lawlor/marvel/input_file/sample_meta_res.txt',
                       sep="\t",header = T)
sj_res=read.table('/AS/lawlor/marvel/input_file/sj_res.txt',
                  sep="\t",header = T)
intron_res=read.table('/AS/lawlor/marvel/input_file/intron_res.txt',
                      sep="\t",header = T)
gene_FPKM_res=read.table('/AS/lawlor/marvel/input_file/gene_FPKM_res.txt',
                         sep="\t",header = T)
gene_FPKM_res[1:4,1:4]
colnames(gene_FPKM_res)[1]='gene_id'
gene_FPKM_res=gene_FPKM_res[!duplicated(gene_FPKM_res$gene_id),]
gene_meta=read.table('/AS/lawlor/marvel/input_file/metadata_res.txt',
                     sep="\t",header = T)
gene_meta=gene_meta[!duplicated(gene_meta$gene_id),]
gtf <- as.data.frame(data.table::fread('/index/gtf/gencode.v38lift37.annotation.gtf', 
                                       sep="\t", header=FALSE, stringsAsFactors=FALSE, quote=""))

df.feature.se=read.table('/AS/lawlor/marvel/input_file/SE_res.txt',
                         sep="\t",header = T)
df.feature.se=df.feature.se[!duplicated(df.feature.se$tran_id),]
df.feature.mxe=read.table('/AS/lawlor/marvel/input_file/MXE_res.txt',
                          sep="\t",header = T)
df.feature.mxe=df.feature.mxe[!duplicated(df.feature.mxe$tran_id),]
df.feature.ri=read.table('/AS/lawlor/marvel/input_file/RI_res.txt',
                         sep="\t",header = T)
df.feature.ri=df.feature.ri[!duplicated(df.feature.ri$tran_id),]
df.feature.a5ss=read.table('/AS/lawlor/marvel/input_file/A5SS_res.txt',
                           sep="\t",header = T)
df.feature.a5ss=df.feature.a5ss[!duplicated(df.feature.a5ss$tran_id),]
df.feature.a3ss=read.table('/AS/lawlor/marvel/input_file/A3SS_res.txt',
                           sep="\t",header = T)
df.feature.a3ss=df.feature.a3ss[!duplicated(df.feature.a3ss$tran_id),]
df.feature.list <- list(df.feature.se,
                        df.feature.mxe,
                        df.feature.ri,
                        df.feature.a5ss,
                        df.feature.a3ss
)
names(df.feature.list) <- c("SE", "MXE", "RI", "A5SS", "A3SS")


marvel <- CreateMarvelObject(SpliceJunction=sj_res,
                             SplicePheno=sample_meta,
                             SpliceFeature=df.feature.list,
                             IntronCounts=intron_res,
                             GeneFeature=gene_meta,
                             Exp=gene_FPKM_res,
                             GTF=gtf
)
class(marvel)


########################## Detect additional events ############################
# Detect AFE
marvel <- DetectEvents(MarvelObject=marvel,
                       min.cells=50,
                       min.expr=1,
                       track.progress=FALSE,
                       EventType="AFE"
)

# Detect ALE
marvel <- DetectEvents(MarvelObject=marvel,
                       min.cells=50,
                       min.expr=1,
                       track.progress=FALSE,
                       EventType="ALE"
)
saveRDS(marvel,file='/AS/lawlor/marvel/input_file/marvel_20230831.rds')


############################# Compute PSI ######################################
# Check splicing junction data
marvel <- CheckAlignment(MarvelObject=marvel, level="SJ")

marvel$SpliceFeature$SE=df.feature.se
marvel$SpliceFeature$MXE=df.feature.mxe
marvel$SpliceFeature$RI=df.feature.ri
marvel$SpliceFeature$A5SS=df.feature.a5ss
marvel$SpliceFeature$A3SS=df.feature.a3ss
# Validate, filter, compute SE splicing events
marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     UnevenCoverageMultiplier=10,
                     EventType="SE"
)

# Validate, filter, compute MXE splicing events    
marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     UnevenCoverageMultiplier=10,
                     EventType="MXE"
)

# Validate, filter, compute RI splicing events      
marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     EventType="RI",
                     thread=4
)

# Validate, filter, compute A5SS splicing events  
marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     EventType="A5SS"
)

# Validate, filter, compute A3SS splicing events  
marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     EventType="A3SS"
)

# Validate, filter, compute AFE splicing events     
marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     EventType="AFE"
)

# Validate, filter, compute ALE splicing events      
marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     EventType="ALE"
)


############################## Transform expression values #####################
marvel <- TransformExpValues(MarvelObject=marvel,
                             offset=1,
                             transformation="log2",
                             threshold.lower=1
)

########################### Check matrices and metadata ########################
# Check splicing data
marvel <- CheckAlignment(MarvelObject=marvel, level="splicing")

# Check gene data
marvel <- CheckAlignment(MarvelObject=marvel, level="gene")

rownames(gene_meta)=gene_meta$gene_id
gene_meta=gene_meta[gene_FPKM_res$gene_id,]
match(gene_FPKM_res$gene_id,gene_meta$gene_id)
head(gene_meta)
rownames(gene_meta)=NULL
marvel$GeneFeature=gene_meta

marvel <- CheckAlignment(MarvelObject=marvel, level="gene")

# Cross-check splicing and gene data
marvel <- CheckAlignment(MarvelObject=marvel, level="splicing and gene")


############################ Overview of splicing events #######################
library(ggplot2)
# Retrieve sample metadata
df.pheno <- marvel$SplicePheno
head(df.pheno)
# Define sample ids
sample.ids=df.pheno$sample.id
sample.ids <- df.pheno[which(df.pheno$Status=="ND"), "sample.id"]

# Tabulate expressed events
marvel <- CountEvents(MarvelObject=marvel,
                      sample.ids=sample.ids,
                      min.cells=10
)

# Output (1): Plot
marvel$N.Events$Plot
ggsave(file='/AS/lawlor/marvel/figure/1.ND_overview.pdf',
       width=3,height=3)
# Output (2): Table
marvel$N.Events$Table
write.table(marvel$N.Events$Table,file='/AS/lawlor/marvel/figure/1.ND_overview_table.txt',
            sep='\t',quote=F,row.names = F)


############################### Modality analysis ##############################
sample.ids=df.pheno$sample.id
sample.ids <- df.pheno[which(df.pheno$Status=="T2D"), "sample.id"]
# Assign modality
marvel <- AssignModality(MarvelObject=marvel,
                         sample.ids=sample.ids,
                         min.cells=10,
                         seed=10
)

marvel$Modality$Results[1:5, c("tran_id", "event_type", "gene_id", "gene_short_name", "modality.bimodal.adj")]
write.table(marvel$Modality$Results[, c("tran_id", "event_type", "gene_id", "gene_short_name", "modality.bimodal.adj")],
            file='/AS/lawlor/marvel/figure/2.T2D_Modality_table.txt',
            sep='\t',quote=F,row.names = F)
# Tabulate modality proportion (overall)
marvel <- PropModality(MarvelObject=marvel,
                       modality.column="modality.bimodal.adj",
                       modality.type="extended",
                       event.type=c("SE", "MXE", "RI", "A5SS", "A3SS", "AFE", "ALE"),
                       across.event.type=FALSE
)

marvel$Modality$Prop$DoughnutChart$Plot
ggsave(file='/AS/lawlor/marvel/figure/2.T2D_Modality.pdf',
       width=3,height=3)

write.table(marvel$Modality$Prop$DoughnutChart$Table,
            file='/AS/lawlor/marvel/figure/2.T2D_Modality_Prop_table.txt',
            sep='\t',quote=F,row.names = F)

# Tabulate modality proportion (by event type)
marvel <- PropModality(MarvelObject=marvel,
                       modality.column="modality.bimodal.adj",
                       modality.type="extended",
                       event.type=c("SE", "MXE", "RI", "A5SS", "A3SS", "AFE", "ALE"),
                       across.event.type=TRUE,
                       prop.test="chisq",
                       prop.adj="fdr",
                       xlabels.size=8
)

marvel$Modality$Prop$BarChart$Plot
ggsave(file='/AS/lawlor/marvel/figure/2.T2D_Modality_Bar.pdf',
       width=6,height=3)

write.table(marvel$Modality$Prop$BarChart$Table,
            file='/AS/lawlor/marvel/figure/2.T2D_Modality_Prop_Bar_table.txt',
            sep='\t',quote=F,row.names = F)

############################# Differential analysis ############################
##### Differential gene expression analysis #####
# Define cell groups
# Retrieve sample metadata
df.pheno <- marvel$SplicePheno
head(df.pheno)

# Cell group 1 (reference)
cell.group.g1 <- df.pheno[which(df.pheno$celltype !="Beta"), "sample.id"]

# Cell group 2
cell.group.g2 <- df.pheno[which(df.pheno$celltype =="Beta"), "sample.id"]

# DE analysis
marvel <- CompareValues(MarvelObject=marvel,
                        cell.group.g1=cell.group.g1,
                        cell.group.g2=cell.group.g2,
                        min.cells=3,
                        method="wilcox",
                        method.adjust="fdr",
                        level="gene",
                        show.progress=FALSE
)

marvel$DE$Exp$Table[1:5, ]

# Differential splicing analysis
marvel <- CompareValues(MarvelObject=marvel,
                        cell.group.g1=cell.group.g1,
                        cell.group.g2=cell.group.g2,
                        min.cells=25,
                        method=c("ad", "dts"),
                        method.adjust="fdr",
                        level="splicing",
                        event.type=c("SE", "MXE", "RI", "A5SS", "A3SS", "ALE", "AFE"),
                        show.progress=FALSE
)

head(marvel$DE$PSI$Table[["ad"]])


############################# Principal component analysis #####################
# # Retrieve non-DE gene_ids
# results.de.exp <- marvel$DE$Exp$Table
# index <- which(results.de.exp$p.val.adj > 0.10 )
# gene_ids <- results.de.exp[, "gene_id"]

# Retrieve DE tran_ids
method <- c("ad", "dts")

tran_ids.list <- list()

for(i in 1:length(method)) {
  
  results.de.psi <- marvel$DE$PSI$Table[[method[i]]]
  index <- which(results.de.psi$p.val.adj < 0.10 & results.de.psi$outlier==FALSE)
  tran_ids <- results.de.psi[index, "tran_id"]
  tran_ids.list[[i]] <- tran_ids
  
}

tran_ids <- unique(unlist(tran_ids.list))


head(marvel$DE$PSI$Table$ad)
tran_ids <- unique(marvel$DE$PSI$Table$ad$tran_id)

table(df.pheno$celltype)
marvel <- RunPCA.PSI(MarvelObject=marvel,
                     cell.group.column="celltype",
                     cell.group.order=unique(df.pheno$celltype),
                     cell.group.colors=NULL,
                     min.cells=25,
                     features=tran_ids,
                     point.size=2.5,
                     method.impute="random",
                     seed=1
)

marvel$PCA$PSI$Plot
ggsave(file='/AS/lawlor/marvel/figure/3.PCA_all_PSI.pdf',
       width=7,height=5)


# Retrieve DE genes
# Retrieve DE result table
results.de.exp <- marvel$DE$Exp$Table    
gene_ids <- results.de.exp[, "gene_id"]


# Reduce dimension
marvel <- RunPCA.Exp(MarvelObject=marvel,
                     cell.group.column="celltype",
                     cell.group.order=unique(df.pheno$celltype),
                     cell.group.colors=NULL,
                     min.cells=25,
                     features=gene_ids,
                     point.size=2.5
)

marvel$PCA$Exp$Plot
ggsave(file='/AS/lawlor/marvel/figure/3.PCA_all_Exp.pdf',
       width=7,height=5)


write.table(marvel$PSI$SE,
            file='/AS/lawlor/marvel/input_file/SE_PSI.txt',
            sep='\t',quote=F)


################### 20231007 ####################################################
### Differential analysis ###
# Define cell groups
# Retrieve sample metadata
df.pheno <- marvel$SplicePheno
head(df.pheno)

# Cell group 1 (reference)
cell.group.g1 <- df.pheno[which(df.pheno$celltype =="Beta" &
                                  df.pheno$Status=='ND'), "sample.id"]

# Cell group 2
cell.group.g2 <- df.pheno[which(df.pheno$celltype =="Beta" &
                                  df.pheno$Status=='T2D'), "sample.id"]

# DE analysis
marvel <- CompareValues(MarvelObject=marvel,
                        cell.group.g1=cell.group.g1,
                        cell.group.g2=cell.group.g2,
                        min.cells=25,
                        method=c("ad", "dts"),
                        method.adjust="fdr",
                        level="splicing",
                        event.type=c("SE"
                                     #, "MXE", "RI", "A5SS", "A3SS", "ALE", "AFE"
                        )
)

head(marvel$DE$PSI$Table[["ad"]])
write.table(marvel$DE$PSI$Table[["ad"]],
            file='/AS/lawlor/marvel/input_file/marvel_DE_PSI_ad.txt',
            sep='\t',quote=F)
write.table(marvel$DE$PSI$Table[["dts"]],
            file='/AS/lawlor/marvel/input_file/marvel_DE_PSI_dts.txt',
            sep='\t',quote=F)
# Differential gene expression analysis
marvel <- CompareValues(MarvelObject=marvel,
                        cell.group.g1=cell.group.g1,
                        cell.group.g2=cell.group.g2,
                        min.cells=3,
                        method="wilcox",
                        method.adjust="fdr",
                        level="gene",
                        show.progress=FALSE
)

marvel$DE$Exp$Table[1:5, ]
write.table(marvel$DE$Exp$Table,
            file='/AS/lawlor/marvel/input_file/marvel_DE_Exp.txt',
            sep='\t',quote=F)

# Differential (spliced) gene analysis
## Next, we will perform differential gene expression analysis only on the 
## differentially spliced genes. This will enable us to investigate the 
## gene-splicing relationship between iPSCs and endoderm cells downstream.
marvel <- CompareValues(MarvelObject=marvel,
                        cell.group.g1=cell.group.g1,
                        cell.group.g2=cell.group.g2,
                        psi.method=c("ad", "dts"),
                        psi.pval=c(0.10, 0.10),
                        psi.delta=0,
                        method.de.gene="wilcox",
                        method.adjust.de.gene="fdr",
                        downsample=FALSE,
                        level="gene.spliced"
)

head(marvel$DE$Exp.Spliced$Table)
write.table(marvel$DE$Exp.Spliced$Table,
            file='/AS/lawlor/marvel/input_file/marvel_DE_Exp.Spliced_SE.txt',
            sep='\t',quote=F)



marvel=readRDS('/AS/lawlor/marvel/input_file/marvel_20230831.rds')
# Volcano plot: Spliced genes
# Plot: Annotate top genes
results <- marvel$DE$Exp.Spliced$Table

index <- which((results$log2fc > 3 | results$log2fc < -3) & -log10(results$p.val.adj) > 15)
gene_short_names <- results[index, "gene_short_name"]

marvel <- PlotDEValues(MarvelObject=marvel,
                       method=c("ad", "dts"),
                       psi.pval=c(0.10, 0.10),
                       psi.delta=0,
                       gene.pval=0.10,
                       gene.log2fc=0.25,
                       point.size=0.1,
                       xlabel.size=8,
                       level="gene.spliced",
                       # anno=TRUE,
                       # anno.gene_short_name=gene_short_names
)

marvel$DE$Exp.Spliced$Summary
marvel$DE$Exp.Spliced$Plot

# Plot DE results
marvel <- PlotDEValues(MarvelObject=marvel,
                       pval=0.10,
                       log2fc=0.25,
                       point.size=0.1,
                       xlabel.size=8,
                       level="gene.global",
                       anno=FALSE
)

marvel$DE$Exp.Global$Plot
marvel$DE$Exp.Global$Summary
head(marvel$DE$Exp.Global$Table[,c("gene_id", "gene_short_name", "sig")])
# Plot DE results with annotation of selected genes
# Retrieve DE output table
results <- marvel$DE$Exp$Table

# Retrieve top genes
index <- which(results$log2fc > 2 | results$log2fc < -2)
gene_short_names <- results[index, "gene_short_name"]

# Plot
marvel <- PlotDEValues(MarvelObject=marvel,
                       pval=0.10,
                       log2fc=0.25,
                       point.size=0.1,
                       xlabel.size=10,
                       level="gene.global",
                       anno=TRUE,
                       anno.gene_short_name=gene_short_names
)

marvel$DE$Exp.Global$Plot
ggsave(file='/AS/lawlor/marvel/figure/4.DE_gene_vol.pdf',
       width=4,height=4)


######################## plot MYL6 #############################################
# Example 1
# Define sample groups
# Retrieve sample metadata
df.pheno <- marvel$SplicePheno

# Cell group 1 (reference)
cell.group.g1 <- df.pheno[which(df.pheno$celltype =="Beta" &
                                  df.pheno$Status=='ND'), "sample.id"]

# Cell group 2
cell.group.g2 <- df.pheno[which(df.pheno$celltype =="Beta" &
                                  df.pheno$Status=='T2D'), "sample.id"]

# Merge
cell.group.list <- list("ND"=cell.group.g1,
                        "T2D"=cell.group.g2
)
# Gene
df.feature <- marvel$GeneFeature
gene_id <- df.feature[which(df.feature$gene_short_name=="MYL6"), "gene_id"]

marvel <- PlotValues(MarvelObject=marvel,
                     cell.group.list=cell.group.list,
                     feature=gene_id,
                     maintitle="gene_short_name",
                     xlabels.size=7,
                     level="gene"
)

plot.1_gene <- marvel$adhocPlot$Exp

# Splicing
a=marvel$DE$PSI$Table[["ad"]]
tran_id <- a[a$gene_short_name=='MYL6'&
               a$p.val<0.05,]$tran_id
tran_id

marvel <- PlotValues(MarvelObject=marvel,
                     cell.group.list=cell.group.list,
                     feature=tran_id,
                     xlabels.size=7,
                     level="splicing",
                     min.cells=25
)

plot.1_splicing <- marvel$adhocPlot$PSI

plot.1_gene+plot.1_splicing
ggsave(file='/AS/lawlor/marvel/figure/5.MYL6_beta.pdf',
       width=4.5,height=2.5)


####################### GO #####################################################
marvel <- BioPathways(MarvelObject=marvel,
                      method=c("ad", "dts"),
                      pval=0.10,
                      species="human"
)

head(marvel$DE$BioPathways$Table)
# Plot top pathways
df <- marvel$DE$BioPathways$Table
go.terms <- df$Description[c(1:10)]

marvel <- BioPathways.Plot(MarvelObject=marvel,
                           go.terms=go.terms,
                           y.label.size=10,
                           x.axis = "pval"
)

marvel$DE$BioPathways$Plot
ggsave(file='/AS/lawlor/marvel/figure/6.GO_beta_T2DvsND.pdf',
       width=5,height=5)
write.csv(marvel$DE$BioPathways$Table,
          file='/AS/lawlor/marvel/input_file/6.GO_beta_T2DvsND_result.csv',
          quote=F)


saveRDS(marvel,file='/AS/lawlor/marvel/input_file/marvel_20230831.rds')
