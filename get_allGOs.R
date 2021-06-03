if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("annotate")

library(annotate)
library(GO.db)

roots = c("GO:0003673","GO:0003674","GO:0005575","GO:0008150","all")

tabla_GOmi=read.table("~/Downloads/Micelioptera/para_GOterms/Anot_unida_Micel.txt", header=TRUE, sep="\t")
head(tabla_GOmi)

# unique GO terms
gos = unique(tabla_GOmi$GO)

head(gos)
# get Ontology, remove obsolete
onts = Ontology(gos)
obs = is.na(onts)
gos = gos[!obs]
onts = onts[!(obs)]

# Obtain all ancestors for our gos

ances = c(
   as.list(GOBPANCESTOR[keys(GOBPANCESTOR) %in% gos]),
   as.list(GOMFANCESTOR[keys(GOMFANCESTOR) %in% gos]),
   as.list(GOCCANCESTOR[keys(GOCCANCESTOR) %in% gos])
)
# remove roots
ances = lapply(ances, function(x) x[!x %in% roots])
# add actual GO into the list
allGos = lapply(names(ances), function(x) c(x, ances[[x]]))
names(allGos) = names(ances)

# Process annotation table
# remove obsolete
tabla_GOmi = tabla_GOmi[tabla_GOmi$GO %in% gos,]
head(tabla_GOmi)
# split by gene, expand GO terms
gosByID = split(tabla_GOmi$GO, tabla_GOmi$ID)
gosByIDexp = lapply(gosByID, function(x) unique(unlist(allGos[x])))

# New table, with the expanded list of GOs
allTab = data.frame("ID"=rep(names(gosByIDexp), sapply(gosByIDexp, length)),
               "GO"=unlist(gosByIDexp), stringsAsFactors=FALSE)
allTab$Ontology = Ontology(allTab$GO)
allTab$Term = Term(allTab$GO)

# Save, split and save
bpTab = allTab[allTab$Ontology == "BP",]
mfTab = allTab[allTab$Ontology == "MF",]
ccTab = allTab[allTab$Ontology == "CC",]

write.table(allTab, "GOALL.tab", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(bpTab, "GOBP.tab", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(mfTab, "GOMF.tab", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(ccTab, "GOCC.tab", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
