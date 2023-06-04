library("xlsx")
library(data.table)
library(readxl)

Discoba_xlsx <- read_excel("C:/Users/micka/OneDrive/Bureaublad/ITAP/Discoba.xlsx.xls")
Discoba.dt <- data.table(Discoba_xlsx)
#Find gene IDs that don't match the manually curated labels
Metazoa.dt <- data.table(Metazoa_xlsx)
setnames(Metazoa.dt, "GENE_ID", "MAN_ID")
Metazoa.dt[, GENE_ID := gsub(".*\\|","",`# BLASTP 2.11.0+`)]
GENE_ID.translate.dt <- Metazoa.dt[GENE_ID != MAN_ID,.(GENE_ID,MAN_ID)]
Metazoa.dt[GENE_ID != MAN_ID,GENE_ID := MAN_ID]
setnames(Metazoa.dt, "GENE_ID", "Gene_ID")
Metazoa.dt$MAN_ID <- NULL
Metazoa.dt$`# BLASTP 2.11.0+` <- NULL
#
Discoba.dt[, GENE_ID := gsub(".*\\|","",`# BLASTP 2.11.0+`)]
Discoba.trans.dt <- merge(Discoba.dt, GENE_ID.translate.dt, by = "GENE_ID", all.x = T)
Discoba.trans.dt[!is.na(MAN_ID) & GENE_ID != MAN_ID, GENE_ID := MAN_ID]
setnames(Discoba.trans.dt, "GENE_ID", "Gene_ID")
Discoba.trans.dt$MAN_ID <- NULL
Discoba.trans.dt$`# BLASTP 2.11.0+` <- NULL
#


nITAP_quality <- read.xlsx("C:/Users/micka/OneDrive/Bureaublad/ITAP/ITAP_quality.xlsx",sheetIndex = 1)
nITAP_quality$Metazoa <- NULL
nITAP_quality$Discoba <- NULL
nITAP_quality.dt <- data.table(nITAP_quality)[!is.na(Gene_ID)]
nITAP.metaz.dt <- merge(nITAP_quality.dt, Metazoa.dt, by = "Gene_ID", all = T)
nITAP.metaz.disco.dt <- merge(nITAP.metaz.dt, Discoba.trans.dt, by = "Gene_ID", all = T)


GeneID_list<-read.xlsx("C:/Users/micka/OneDrive/Bureaublad/ITAP/GeneID list.xlsx", sheetIndex = 1)
GeneID_list.dt <- data.table(GeneID_list)[!is.na(Gene_ID)]
GeneID_list.dt[,Gene:=gsub(x= Gene_ID, pattern = "_.*", replacement = "" )]
#$Gene <- gsub(x= GeneID_list$Gene_ID, pattern = "_.*", replacement = "" )
functionITAP.dt <- merge(nITAP.metaz.disco.dt, GeneID_list.dt, all.y = T)
functionITAP.clean.dt <- functionITAP.dt[ Gene != "VAPA" ]
#nITAP <- data.table(functionITAP[order(functionITAP$Function),])[!is.na(Gene_ID)]
colnames(functionITAP.clean.dt) <- as.character(colnames(functionITAP.clean.dt))


library(reshape2)
ITAPQual_GGPLOT <- melt.data.table(functionITAP.clean.dt, id.vars = c("Gene_ID","Gene", "Function"), variable.name = "Taxon", value.name = "Evalue")
class(ITAPQual_GGPLOT$Evalue)
#evalue is read as a character object, meaning it r believes it to be text, so we convert it to "numeric"
ITAPQual_GGPLOT$Evalue <- as.numeric(ITAPQual_GGPLOT$Evalue)
#let's check if that worked
class(ITAPQual_GGPLOT$Evalue)

#Let us get rid of NA's, these will just not be drawn in the plot. 
#Also gets rid of rows that are all NAs
ITAPQual_GGPLOT.noNA <- ITAPQual_GGPLOT[!is.na(ITAPQual_GGPLOT$Evalue),]

ecutoff <- 0.00001
ITAPQual_GGPLOT.noNA.merged.dt <- ITAPQual_GGPLOT.noNA
ITAPQual_GGPLOT.noNA.merged.dt[,pres := 0]
ITAPQual_GGPLOT.noNA.merged.dt[Evalue < ecutoff,pres := 1]
ITAPQual_GGPLOT.noNA.merged.dt[,anypres:= max(pres), by =.(Gene)]
ITAPQual_GGPLOT.perGene.dt <- unique(ITAPQual_GGPLOT.noNA.merged.dt[,.(Gene,Taxon,Function,Present = anypres)])
ITAPQual_GGPLOT.perGene.dt$Gene <- factor(ITAPQual_GGPLOT.perGene.dt$Gene, levels =  unique(ITAPQual_GGPLOT.perGene.dt[order(Function,decreasing = T)]$Gene))
ITAPQual_GGPLOT.perGene.dt[Taxon == "Chloroplastida", Taxon := "Chlorophyta"]
ITAPQual_GGPLOT.perGene.dt[Taxon == "higher_plants", Taxon := "Streptophyta"]
Grouporder.v <- c(
  "Rhizaria",
  "Alveolata",
  "Stramenopiles",
  "Telonemia",
  "Katablepharida",
  "Palpitomonas",
  "Cryptophyta",
  "Haptophyta",
  "Centrohelida",
  "Glaucophyta",
  "Streptophyta",
  "Chlorophyta",
  "Rhodophyta",
  "Fungi",
  "Amoebozoa",
  "Metazoa",
  "Discoba",
  "Malawimonadida",
  "Ancyromonadidae",
  "Bacteria",
  "Archaea"
  )

ITAPQual_GGPLOT.perGene.dt$Taxon <- factor(
  ITAPQual_GGPLOT.perGene.dt$Taxon,
  levels =  Grouporder.v
  )
Telonemia.dt <- data.table("Gene" = "APC", "Taxon" = "Telonemia", "Function" = "Organization", "Present" = 0)
ITAPQual_GGPLOT.perGene.plusT.dt <- rbind(ITAPQual_GGPLOT.perGene.dt, Telonemia.dt)
library(ggplot2)
ggplot(ITAPQual_GGPLOT.perGene.plusT.dt, aes(x = `Gene_ID`, y = Taxon, fill = Evalue)) +
  geom_tile() + 
  scale_fill_gradient(low = "blue", high="white") +
  theme_light() + #just my preference, feel free to get rid of it
  theme(axis.text.x = element_text(angle = 90))

ggplot(ITAPQual_GGPLOT.perGene.plusT.dt, aes(x = Gene, y = Taxon, fill = Present)) +
  geom_tile(color="black") + 
  scale_fill_gradient(low = "white", high="blue") +
  theme_light() + #just my preference, feel free to get rid of it
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_discrete(position = "top")+
  xlab("") +
  ylab("")

