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
files <- file.path("/workspace/class-rnaseq/analysis_tutoring01/reads/", paste0(dataset$sample,".quant"), "quant.sf")


# Assegnazione dei nomi dei campioni come nomi dei campioni in files
names(files) <- dataset$sample


# Lettura del file di associazione tra ID trascritti e ID geni
tx2gene <- read_tsv("/workspace/class-rnaseq/datasets_reference_only/trascriptome/gencode.v29.transcripts_no-vers_chr21_tx2gene.txt")


# Importazione dati di quantificazione genica ottenuti con Salmon
# tximport() --> funzione che permette di aggregare i dati da livello di trascritto a livello di gene
# files --> vettore con i path dei file di quantificazione. Nomi dei campioni assegnati con il precedente comando
# type = "salmon" --> specifica che la quantificazione proviene da Salmon
# tx2gene --> file di associazione tra ID trascritti e ID geni
# txi --> lista contenente varie componenti (counts, abundance and length)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)


# Verifica della funzionalità di tximport
inspecition <- as.data.frame(txi$counts) %>% 
  rownames_to_column(var = "gene_name")


colnames(txi$counts) # Verifica nomi delle colonne di txi$counts
rownames(dataset)    # Verifica nomi delle righe di dataset


# Allineare nomi delle righe di dataset con nomi delle colonne di txi$counts
# Collegamento tra i metadati (dataset) e i dati di espressione genica (txi)
rownames(dataset) <- colnames(txi$counts)
rownames(dataset)


# La funzione DESeqDataSetFromTximport viene utilizzata per creare un oggetto DESeqDataSet (dds), a partire dai dati di espressione genica e dai metadati associati
# Input principali:
# - txi: oggetto prodotto da tximport
# - dataset: dataframe di metadati che descrive le condizioni sperimentali dei campioni
# - ~ condition: formula che specifica la variabile sperimentale (nel nostro caso "condition"). Questo permette a DESeq2 di sapere quali gruppi comparare
dds <- DESeqDataSetFromTximport(txi, dataset, ~ condition)



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
dds$condition <- relevel(dds$condition, ref = "control")



###########################
# Differential expression #
###########################

# Utilizzo del pacchetto DESeq2 per eseguire un'analisi di espressione differenziale
# La funzione `DESeq` applica il modello per identificare i geni differenzialmente espressi
dds <- DESeq(dds)



############################
# Estrazione dei risultati #
############################


# Estrazione dei risultati dell'analisi di espressione differenziale
res <- results(dds)

# Ordina in modo cresecente i risultati in base al valore di p-value
# I geni con p-value più basso (più significativi) saranno in cima alla lista
resOrdered <- res[order(res$pvalue),]

# Crea un plot di tipo MA per visualizzare il log2foldchange rispetto alla media normalizzata
plotMA(res, ylim = c(-3, 3))

# Visualizzazione le stime di dispersione per i dati del modello DESeq2
plotDispEsts(dds)

# Crea un grafico dei conteggi normalizzati per il gene con il valore di padj più basso
# gene = which.min(res$padj) individua e restituisce l'indice del gene con il valore di padj più basso (più significativo)
# `intgroup` specifica la variabile del design che deve essere utilizzata nel grafico
plotCounts(dds, gene = which.min(res$padj), intgroup = "condition")



#########################
# Salvataggio risultati #
#########################


# Converte i risultati ordinati in un tibble
resdata <- as_tibble(resOrdered)

# Aggiunge una colonna `gene` contenente i nomi dei geni, presi dai nomi delle righe dell'oggetto `resOrdered`
resdata$gene <- rownames(resOrdered)

# Mettere gene come prima colonna
resdata <- resdata %>%
  relocate(gene, .before = everything())

# Salvataggio come file TSV (tab-separated values)
write_tsv(resdata, "analysis_results.tsv")



##############
# CLUSTERING #
##############


# Normalizzazione dei conteggi, molto utile per visualizzazioni downstream (heatmap e PCA, per esempio)
# In DESeq2 esistono poi anche trasformazioni più complesse per stabilizzare la varianza (vst e rlog)
ntd <- normTransform(dds)

# Selezione dei 20 geni con la media di espressione più alta
# rowMeans() --> calcola la media per ciascun gene nei vari campioni
# (counts(dds, normalized = TRUE)) --> estrae i conteggi normalizzati dal dds
# order() --> ordina in modo decresecente i geni in base alla media
# [1:20] --> seleziona i primi 20 geni (quindi i 20 con la media più alta)
select <- order(rowMeans(counts(dds, normalized = TRUE)),
                decreasing = TRUE)[1:20]

# Creazione di un dataframe con le informazioni sulla condizione
# colData(dds) --> estrae i metadati dal dds
# [,c("condition")] --> seleziona solo la colonna "condition"
df <- as.data.frame(colData(dds)[,c("condition")])

# Creazione di una heatmap dei 20 geni più espressi
# assay(ntd) --> estrae i conteggi normalizzati
# [select,] --> seleziona solo le righe corrispondenti ai 20 geni selezionati
# cluster_cols = FALSE --> disabilita il clustering delle colonne
# annotation_col = df$condition --> colora le colonne in base alla condizione
pheatmap(assay(ntd)[select,],
         cluster_cols = FALSE, annotation_col = df$condition)

# Creazione del grafico di PCA per visualizzare la variabilità tra i campioni utilizzando la variabile condition
plotPCA(ntd, intgroup = c("condition"))



####################################################
# Estrazione dei geni statisticamente significanti #
####################################################

# Creazione dell' "universo", un dataframe che associa identificatori di geni (ENTREZID, SYMBOL, ENSEMBL, ENSEMBLTRANS) utilizzando il database org.Hs.eg.db
# Ogni chiave di tipo 'ENTREZID' è mappata con i suoi corrispondenti valori
universe <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = keys(org.Hs.eg.db),
                                  columns = c('ENTREZID','SYMBOL','ENSEMBL','ENSEMBLTRANS'),
                                  keytype = 'ENTREZID')

# Estrazione dei geni con padj < a 0.05.
sig_genes <- resdata$gene[which(resdata$padj<0.05)]

# Conversione degli identificatori ENSEMBL dei geni significativi in identificatori ENTREZID
# Si ottiene un vettore unico di ENTREZID corrispondente ai geni significativi
entrez_genes_sig <- unique(universe[which(universe$ENSEMBL %in% sig_genes),]$ENTREZID)

# Creazione di un vettore di padj per i geni ENSEMBL significativi
# I nomi del vettore sono gli identificatori ENSEMBL dei geni significativi
pvalue_ens_genes <- resdata$padj[which(resdata$padj<0.05)]
names(pvalue_ens_genes) <- sig_genes

# Creazione di un vettore di padj per i geni ENTREZID significativi
# I nomi del vettore sono gli identificatori ENTREZID dei geni significativi
pvalue_entrez_genes <- resdata$padj[which(resdata$padj<0.05)]
names(pvalue_entrez_genes) <- entrez_genes_sig




################################################
# Analisi di enrichment con Gene Ontology (GO) #
################################################


ego <- enrichGO( gene = sig_genes,
                 universe = unique(tx2gene$GENEID),
                 OrgDb = org.Hs.eg.db,
                 keyType = 'ENSEMBL',
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05)

pdf("plots_ego.pdf")
dotplot(ego, showCategory = 30)
dev.off()

pdf("plots_network-ego.pdf")
cnetplot(ego, foldChange = resdata$log2FoldChange[which(resdata$padj<0.5)])
dev.off()



####################
# DISGNET ANALYSIS #
####################


# DisGeNet è un database di associazioni tra malattie e geni

# Lettura del file di associazione tra malattie e geni
gda <- read_tsv(gzfile("/workspace/class-rnaseq/datasets_reference_only/trascriptome/all_gene_disease_associations.tsv.gz"))

# Creazione dei dataframe per l'arricchimento
# disease2gene --> associazione tra malattie e geni
# disease2name --> associazione tra malattie e nomi delle malattie
disease2gene = gda[, c("diseaseId", "geneId")]
disease2name = gda[, c("diseaseId", "diseaseName")]

# Arricchimento funzionale con enrcher(), funzione del pacchetto clusterProfiler
# entrez_genes_sig --> vettore di ID dei geni significativi
disgnet = enricher(entrez_genes_sig, 
                   TERM2GENE = disease2gene, 
                   TERM2NAME = disease2name)

# Creazione di un concept network plot per visualizzare le relazioni tra le malattie e i geni significativi associati
# foldChange = resdata$log2FoldChange[which(resdata$padj<0.5)] --> specifica il log2FoldChange per ciascun gene significativo (padj < 0.5)
cnetplot(disgnet, foldChange = resdata$log2FoldChange[which(resdata$padj<0.5)])