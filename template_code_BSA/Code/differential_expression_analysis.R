####################################
# Differential Expression Analysis #
####################################


# Da terminale creazione cartella per salvataggio tabelle e grafici
cd ..
mkdir -p results_de_analysis


# Ci spostiamo ora sulla console
# Settaggio working directory
setwd("/workspace/class-rnaseq/analysis_tutoring01")

# Caricamento delle librerie
library(DESeq2)
library(tximport)
library(tidyverse)
library(pheatmap)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)



############################
# Preparazione del dataset #
############################


# Creazione dei metadati
dataset <- tibble(
  sample = c("sample_01",
             "sample_02",
             "sample_03",
             "sample_04",
             "sample_05",
             "sample_06"),
  condition = c(rep("control", 3),
                rep("case", 3))
)



#######################################
# Lettura dei file di quantificazione #
#######################################


# Creazione di un vettore con i path dei file di quantificazione
# file.path() --> funzione per costruire i path dei file
# paste0() --> funzione per combinare il nome di ogni campione con il suffisso .quant
# quant.sf --> nome del file di quantificazione presente in ogni sottocartella
files <- file.path("/workspace/class-rnaseq/analysis_tutoring01/reads", paste0(dataset$sample,".quant"), "quant.sf")
files
names(files)


# Associare il file con la quantificazione al nome del campione
names(files) <- dataset$sample
files
names(files)


# Lettura del file di associazione tra ID trascritti e ID geni
tx2gene <- read_tsv("/workspace/class-rnaseq/datasets_reference_only/trascriptome/gencode.v29.transcripts_no-vers_chr21_tx2gene.txt")
head(tx2gene)


# Importazione dati di quantificazione genica ottenuti con Salmon
# tximport() --> funzione che permette di aggregare i dati da livello di trascritto a livello di gene
# files --> vettore con i path dei file di quantificazione e nomi dei campioni assegnati precedentemente
# type = "salmon" --> specifica che la quantificazione proviene da Salmon
# tx2gene --> file di associazione tra ID trascritti e ID geni
# txi --> lista contenente varie componenti (counts, abundance and length)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)


# Verifica della funzionalità di tximport prendendo come esempio il gene ENSG00000277117
counts <- as.data.frame(txi$counts) %>% 
  rownames_to_column(var = "gene_name")


# Un punto importante per evitare errori nei passaggi successivi è allineare i nome delle righe dei metadati con i nomi delle colonne di txi$counts
colnames(txi$counts) # Verifica nomi delle colonne di txi$counts
rownames(dataset)    # Verifica nomi delle righe dei metadati


# Allineare nomi delle righe dei metadati con nomi delle colonne di txi$counts
# Collegamento tra i metadati (dataset) e i dati di espressione genica (txi)
rownames(dataset) <- colnames(txi$counts)
colnames(txi$counts)
rownames(dataset)


# La funzione DESeqDataSetFromTximport viene utilizzata per creare un oggetto DESeqDataSet (dds), a partire dai dati di espressione genica e dai metadati associati
# Input principali:
# - txi: oggetto prodotto da tximport
# - dataset: dataframe di metadati che descrive le condizioni sperimentali dei campioni
# - ~ condition: formula che specifica la variabile sperimentale (nel nostro caso "condition"). Questo permette a DESeq2 di sapere quali gruppi comparare
dds <- DESeqDataSetFromTximport(txi, dataset, ~ condition)


# Ispezione dell'oggetto dds

colData(dds) # metadati associati ai campioni

design(dds) # design del modello

head(counts(dds))


# Salvataggio immagine in caso di chiusura di gitpod
save.image("tutoring01.RData")
# Per ricaricare l'mmagine comando load("tutoring01.RData")



#########################################
# Filtraggio in base al numero di reads #
#########################################


# La funzione counts() estrae la matrice di conteggi dal DESeqDataSet
# rowSums() calcola la somma dei conteggi per ciascun gene su tutti i campioni
# Vengono mantenuti solo i geni con una somma di conteggi >= 10
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


# Impostazione del livello di riferimento per la condizione sperimentale
# relevel() consente di definire il livello di base (controllo) per la variabile condition
# Importante per l'interpretazione dei risultati, poiché i confronti saranno fatti rispetto al controllo
dds$condition
dds$condition <- relevel(dds$condition, ref = "control")
dds$condition



###########################
# Differential expression #
###########################

# Utilizzo del pacchetto DESeq2 per eseguire un'analisi di espressione differenziale
dds <- DESeq(dds)



##########################################
# Analisi esplorativa di Quality Control #
##########################################


# Normalizzazione dei conteggi, molto utile per visualizzazioni downstream (heatmap e PCA, per esempio)
# In DESeq2 esistono poi anche trasformazioni più complesse per stabilizzare la varianza (vst e rlog)
ntd <- normTransform(dds)

# Selezione dei 20 geni con la media di espressione più alta
# (counts(dds, normalized = TRUE)) --> estrae la matrice dei conteggi normalizzati dal dds
# rowMeans() --> calcola la media di ciascuna riga (cioè per ogni gene) nei vari campioni
# order() --> restituisce un vettore di indici ordinati in modo decresecente (decreasing = TRUE) in base alla media
# [1:20] --> seleziona i primi 20 geni dal vettore ordinato (quindi i 20 con la media più alta)
select <- order(rowMeans(counts(dds, normalized = TRUE)),
                decreasing = TRUE)[1:20]

# Creazione di un dataframe con le informazioni sulla condizione
# colData(dds) --> estrae i metadati dal dds
# [,c("condition")] --> seleziona solo la colonna "condition"
df <- as.data.frame(colData(dds)[,c("condition")])

# Creazione di una heatmap dei 20 geni più espressi
# assay(ntd) --> estrae i conteggi trasformati
# [select,] --> seleziona solo le righe corrispondenti ai 20 geni selezionati
# cluster_cols = FALSE --> disabilita il clustering delle colonne
# annotation_col = df$condition --> colora le colonne in base alla condizione
pdf("results_de_analysis/heatmap.pdf")
pheatmap(assay(ntd)[select,],
         cluster_cols = FALSE, annotation_col = df$condition)
dev.off()

# Creazione del grafico di PCA per visualizzare la variabilità tra i campioni utilizzando la variabile condition
pdf("results_de_analysis/plotPCA.pdf")
plotPCA(ntd, intgroup = c("condition"))
dev.off()

# Visualizzazione le stime di dispersione per i dati del modello DESeq2
pdf("results_de_analysis/plotDispEsts.pdf")
plotDispEsts(dds)
dev.off()

# Salvataggio immagine in caso di chiusura di gitpod
save.image("tutoring01.RData")



############################
# Estrazione dei risultati #
############################


# Estrazione dei risultati dell'analisi di espressione differenziale
res <- results(dds)
res

# Ordina in modo cresecente i risultati in base al valore di p-value
# I geni con p-value più basso (più significativi) saranno in cima alla lista
resOrdered <- res[order(res$pvalue),]
resOrdered


# Plot esplorativi dei risultati


# Crea un plot di tipo MA per visualizzare il log2foldchange rispetto alla media normalizzata
pdf("results_de_analysis/plotMA.pdf")
plotMA(res, ylim = c(-3, 3))
dev.off()


# Crea un grafico dei conteggi normalizzati per il gene con il valore di padj più basso
# gene = which.min(res$padj) individua e restituisce l'indice del gene con il valore di padj più basso (più significativo)
# `intgroup` specifica la variabile del design che deve essere utilizzata nel grafico
# gene ENSG00000160191
pdf("results_de_analysis/plotCounts.pdf")
plotCounts(dds, gene = which.min(res$padj), intgroup = "condition")
dev.off()

ENSG00000160191 <- counts %>% 
  filter(gene_name == "ENSG00000160191")



#############################
# Salvataggio dei risultati #
#############################


# Converte i risultati ordinati in un tibble
resdata <- as_tibble(resOrdered)

# Aggiunge una colonna `gene` contenente i nomi dei geni, presi dai nomi delle righe dell'oggetto `resOrdered`
resdata$gene <- rownames(resOrdered)

# Mettere gene come prima colonna
resdata <- resdata %>%
  relocate(gene, .before = everything())

# Salvataggio risultati come file TSV (tab-separated values)
write_tsv(resdata, "results_de_analysis/results.tsv")

# Estrazione dei geni significativi (padj < 0.05)
significant_genes <- as_tibble(resdata %>%
                                 filter(padj < 0.05))

# Salvataggio risultati significativi
write_tsv(significant_genes, "results_de_analysis/significant_genes.tsv")

# Salvataggio immagine in caso di chiusura di gitpod
save.image("tutoring01.RData")




################################################
# Analisi di enrichment con Gene Ontology (GO) #
################################################


# Estrazione geni significativi (come fatto sopra) ma questa volta sotto forma di vettore
# [which(resdata$padj<0.05)] --> restituisce gli indici delle righe di resdata in cui la condizione padj < 0.05 è soddisfatta
# resdata$gene --> restituisce i nomi dei geni corrispondenti a tali indici
sig_genes <- resdata$gene[which(resdata$padj < 0.05)]

# Analisi di arricchimento che consente di identificare i termini GO (Gene Ontology) rappresentati tra i geni significativi
# in particolare ci concentriamo sui termini GO relativi ai processi biologici (BP), alle funzioni molecolari (MF) e ai componenti cellulari (CC)
ego_BP <- enrichGO(gene = sig_genes,
                   universe = unique(tx2gene$GENEID),
                   OrgDb = org.Hs.eg.db,
                   keyType = 'ENSEMBL',
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05)
ego_BP

ego_MF <- enrichGO(gene = sig_genes,
                   universe = unique(tx2gene$GENEID),
                   OrgDb = org.Hs.eg.db,
                   keyType = 'ENSEMBL',
                   ont = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05)
ego_MF

ego_CC <- enrichGO(gene = sig_genes,
                   universe = unique(tx2gene$GENEID),
                   OrgDb = org.Hs.eg.db,
                   keyType = 'ENSEMBL',
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05)
ego_CC

# Creazione di un dotplot che consente la visualizzazione dei primi 30 termini GO più significativi in base al valore di padj
pdf("results_de_analysis/dotplot_ego.pdf")
dotplot(ego_MF, showCategory = 10)
dev.off()

# Creazione di un gene-concept network per visualizzare la connessione tra i geni e i termini GO arricchiti
pdf("results_de_analysis/cnetplot_ego.pdf")
cnetplot(ego_MF, foldChange = resdata$log2FoldChange[which(resdata$padj<0.5)], showCategory = 10)
dev.off()



#####################
# DisGeNET analysis #
####################


# DisGeNET è un database di associazioni tra malattie e geni
# Utile per vedere se alcuni dei nostri geni significativi sono già stati descritti in alcune patologie

# Lettura del file di associazione tra malattie e geni
gda <- read_tsv(gzfile("/workspace/class-rnaseq/datasets_reference_only/trascriptome/all_gene_disease_associations.tsv.gz"))

# Creazione dei dataframe per l'arricchimento
# disease2gene --> associazione tra malattie e geni
# disease2name --> associazione tra malattie e nomi delle malattie
disease2gene = gda[, c("diseaseId", "geneId")]
disease2name = gda[, c("diseaseId", "diseaseName")]

# Creazione di un universo di geni che contanga per ogni gene annotato nel genoma umano (org.Hs.eg.db):
# ENTREZID --> ID assegnato da NCBI Entrez Gene
# SYMBOL --> simbolo del gene
# ENSEMBL --> ID del gene nel database Ensembl
# ENSEMBLTRANS --> ID del trascritto del gene nel database Ensembl
universe <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = keys(org.Hs.eg.db),
                                  columns = c('ENTREZID','SYMBOL','ENSEMBL','ENSEMBLTRANS'),
                                  keytype = 'ENTREZID')
                        
# Filtraggio geni significativi dall'universo basandosi sugli ID ENSEMBL
# unique() --> rimuove i duplicati
# which() --> restituisce gli indici degli elementi che soddisfano la condizione
# universe$ENSEMBL %in% sig_genes --> verifica se gli ID ENSEMBL sono presenti tra i geni significativi
# $ENTREZID --> restituisce gli ID Entrez corrispondenti ai geni significativi
entrez_genes_sig <- unique(universe[which(universe$ENSEMBL %in% sig_genes),]$ENTREZID)

# Arricchimento funzionale con enrcher(), funzione del pacchetto clusterProfiler
# entrez_genes_sig --> vettore di ID dei geni significativi
disgnet = enricher(entrez_genes_sig, 
                   TERM2GENE = disease2gene, 
                   TERM2NAME = disease2name)


# Creazione di un gene-concept network per visualizzare le relazioni tra le malattie e i geni significativi associati
# foldChange = resdata$log2FoldChange[which(resdata$padj<0.5)] --> specifica il log2FoldChange per ciascun gene significativo (padj < 0.5)
pdf("results_de_analysis/cnetplot_disgnet.pdf")
cnetplot(disgnet, foldChange = resdata$log2FoldChange[which(resdata$padj<0.5)], showCategory = 10)
dev.off()



##########
# GitHub #
##########


# Prima di pushare su GitHub rimuoviamo prima le cartelle che contengono file pesanti come le reads
# Assicuriamoci di essere in /workspace/class-rnaseq sul terminale (non sulla console)
rm -rf */.git
rm -r dataset_tutoring_rnaseq01
rm -r analysis_tutoring01/reads
rm analysis_tutoring01/tutoring01.RData 

git add *
git commit -m "analisi RNA-Seq primo tutorato"
git push
