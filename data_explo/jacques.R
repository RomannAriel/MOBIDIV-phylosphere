# Charger un fichier RDS nommé "data.rds" dans une variable nommée "data"
#file<-"C:/Users/david/Documents/RTRA et Campus/PPR 2019/Mobidiv/stage Romann Charbonnier/sequences/tax_final.rds"


if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

# BiocManager::install("msa")

# Charger le package seqinr
library(seqinr)
library(ape)
library(msa)
library(pegas)
library(stringr)


# Etape 1 : données brutes
## lecture du fichier issus de l'assignation
file<-"C:/Users/david/Documents/RTRA et Campus/PPR 2019/Mobidiv/stage Romann Charbonnier/sequences/tax_euka2207final.rds"

data <- readRDS(file)
data<-as.data.frame(data)
dim(data)

# compter les reads de chacun des clusters
file<-"C:/Users/david/Documents/RTRA et Campus/PPR 2019/Mobidiv/stage Romann Charbonnier/sequences/seqtab_final.rds"

reads <- readRDS(file)
reads<-as.data.frame(reads)
dim(reads)

TOTreads<-colSums(reads[1:ncol(reads)])
length(TOTreads)

data$Totreads<-TOTreads

# Analyse des reads par origine
levels(as.factor(data$tax.Kingdom))

data<-data[grep("k__Fungi",data[,1] ),]
dim(data)
names(data)
fam<-levels(as.factor(data$tax.Family))

B<-data[grep("thora",data$tax.Family),c(4:7,15)]
rownames(B)<-1:nrow(B)
B

C<-aggregate(data$Totreads, list(phylum=data$tax.Genus, species=data$tax.Species), sum )
print(C[order(C$x, decreasing=T),])
hist(C$x)

data<-data[grep("g__Blumeria",data[,6] ),]
dim(data)

data<-data[grep("s__graminis",data[,7] ),]
dim(data)

hist(data$boot.Genus)
hist(data$boot.Species)

plot(data$boot.Genus, data$boot.Species)
dim(data)


## construire le fichier texte contenant les noms et les séquences
names <- data[, 6]
seqs <- rownames(data)

lines<-c()
for (i in 1:length(names) )
{ lines[2*(i-1)+1]<-paste0(">",names[i],i)
lines[2*i]<-seqs[i]
}

## Ecrire le fichier fasta avec les noms et les séquences
write(lines, file= "sequences.fasta")

# si on utilise le package biocmanager 
# https://bioconductor.org/packages/release/bioc/vignettes/msa/inst/doc/msa.pdf

## lecture du fichier fasta brut de départ
mySequences <- readDNAStringSet("sequences.fasta", format="fasta")
mySequences
names(mySequences)

## alignement CLUSTALW par defaut
myFirstAlignment <- msa(mySequences)

## ca peut aussi etre msa(mySequences, "Muscle") si on veut utliser d autres algos

## conversion
align_obj <- msaConvert(myFirstAlignment, type = "seqinr::alignment")

## Puis convertir l'objet alignment en un objet DNAbin avec la fonction as.DNAbin du package ape
dnabin_obj <- as.DNAbin(align_obj)

# Detecter les ratons laveurs  
## Calculer la matrice de distances
dist <- dist.dna(dnabin_obj, model = "JC69")

##  Estimer un arbre avec la fonction nj
tree <- nj(dist)
plot(tree)

hist(dist)
D<-as.matrix(dist)
dim(D)
hist(colSums(D))
liste<-names(which(colSums(D)<100))


# Eliminer ces séquences du fichier de départ
##  lecture du fichier

mySequences <- readDNAStringSet("sequences.fasta", format="fasta")
##  Éliminer des  séquences de la liste  avec subset

if(length(liste)!=0) {aa <- mySequences[which(names(mySequences) %in%  liste)]
}else {aa <- mySequences}

length(aa)
width(aa)

## re - alignement CLUSTALW par defaut
myFirstAlign <- msa(aa)

## conversion
bb <- msaConvert(myFirstAlign, type = "seqinr::alignment")
# Puis convertir l'objet alignment en un objet DNAbin avec la fonction as.DNAbin du package ape
bb <- as.DNAbin(bb)

## Exporter l'objet dnabin au format fasta
write.dna(bb, file = "bb.fasta", format = "fasta")

# Etape 2 : Détection des Indels et modification des séquences 

## lecture du fichier sans ratons laveurs
mySeque_bb <- readDNAStringSet("bb.fasta",format="fasta")
mySeque_bb
names(mySeque_bb)
length(mySeque_bb)
width(mySeque_bb)

seq <- as.character(mySeque_bb[names(mySeque_bb)])

## Compter le nombre de "-" par position

mat<-
  sapply(seq, function(x) {
    v <- strsplit(x, "")[[1]]
    return(v)})

Count <-apply(mat,1, function(x) {sum(str_count(x, "-")) } )
plot(Count)

# on ne garde que les bases qui n'ont pas d'indels
liste <-which(Count==0)

c<-  sapply(seq, function(x) { y<-str_sub(x, start=liste, end=liste)
y<-paste(y, collapse="")
return(y)} )

# on remet dans l'alignement 
mySeque_bb[names(mySeque_bb)]<-c

# Exporter l'objet dnabin au format fasta
write.dna(mySeque_bb[names(mySeque_bb)], file = "Racoon_and_indel_free.fasta", format = "fasta")

mySeque_bb <- readDNAStringSet("Racoon_and_indel_free.fasta", format="fasta")

car<-names(mySeque_bb)
car<-str_replace(car,"g__Blumeria","")
names(mySeque_bb)<-car

# conversion
rm(align)
align<- msa(mySeque_bb)

align_obj <- msaConvert(align, type = "seqinr::alignment")


# Puis convertir l'objet alignment en un objet DNAbin avec la fonction as.DNAbin du package ape
rm(dnabin_obj)
dnabin_obj <- as.DNAbin(align_obj)

# Calculer la matrice de distances

D <- dist.dna(dnabin_obj)

D<-as.matrix(D)
dim(D)

D[2,1]

# recherche des paires identiques
hist(D)
LIGNE<-1+round(which(D==0)/454)
#Colonne
COLONNE<-which(D==0) %% 454

COORD <- cbind(LIGNE, COLONNE)
COORD<-as.data.frame(COORD)

COORD <- COORD[which(COORD$LIGNE!=COORD$COLONNE),]
COORD[1,]

colnames(D)[COORD[1,1]]
rownames(D)[COORD[1,2]]

heatmap(D)

# Estimer un arbre avec la fonction nj
tree <- nj(D)

# Visualiser l'arbre avec la fonction plot
plot(tree,cex=0.1)

write.tree(tree, file = "mon_arbre.newick")


# Construction du réseau d'haplotypes

# Créer un objet haplotype à partir de l'objet DNAbin filtré
h <- haplotype(dnabin_obj)

# Inférer et représenter un réseau d'haplotypes avec la méthode TCS
net <- haploNet(h)
plot(net, cex=0.1)


# https://bioconductor.org/packages/release/bioc/manuals/fastreeR/man/fastreeR.pdf

# Clustering 
D <- dist.dna(dnabin_obj)

D<-as.matrix(D)

# Étape 3 : Clustering hiérarchique
hc <- hclust(as.dist(D), method = "complete")

# Couper le dendrogramme pour obtenir le nombre souhaité de clusters
num_clusters <- 3
cluster_ids <- cutree(hc, k = num_clusters)


table(cluster_ids)
cluster_ids<-as.data.frame(cluster_ids)

cluster_ids$count<-data$Totreads[as.numeric(rownames(cluster_ids))]

A<-aggregate(cluster_ids$count, list(Cluster=cluster_ids$cluster_ids), sum)

A$x[1]/sum(A$x)
A$x[2]/sum(A$x)
A$x[3]/sum(A$x)

datanames(mySeque_bb)

merge()

# Étape 4 : Évaluation des clusters
# Utilisez des mesures d'évaluation appropriées pour évaluer la qualité de vos clusters
# Convertir les distances en similarités
max_dist <- max(D)
similarity_matrix <- 1 - (D / max_dist)

library(cluster) # Assurez-vous d'avoir installé le package cluster



## en reseau
library(igraph)

# Calculer la distance génétique entre les séquences
dist_matrix <- D

# Convertir les distances en similarités (ou en différences) pour construire le réseau
# Par exemple, vous pouvez utiliser une transformation de similarité comme la distance maximale moins les distances
max_dist <- max(dist_matrix)
similarity_matrix <- 1 - (dist_matrix / max_dist)

# Créer le réseau à partir de la matrice de similarité avec un seuil pour les arêtes
max_dist <- max(dist_matrix)
similarity_matrix <- 1 - (dist_matrix / max_dist)

# Seuillage manuel pour déterminer les arêtes du réseau
threshold <- 0.7# Ajustez ce seuil selon vos besoins
network_adj <- similarity_matrix > threshold

# Créer le réseau à partir de la matrice d'adjacence
network <- graph_from_adjacency_matrix(network_adj, mode = "undirected", diag = FALSE)

network <-graph_from_adjacency_matrix(similarity_matrix)

# Attribuer les noms de séquence comme étiquettes de noeud
V(network)$label <- rownames(similarity_matrix)

# Visualiser le réseau
plot(network, layout = layout_with_fr(network))  # Utilisation de l'algorithme de Fruchterman-Reingold pour disposer les noeuds

plot(network, vertex.label.cex = 0.8, edge.arrow.size = 0.5)


# Créer le dendrogramme basé sur la distance génétique
dendrogram <- hclust(as.dist(dist_matrix), method = "average")  # Méthode d'agrégation utilisée pour le regroupement hiérarchique (ici, average linkage)
library("dendextend")

# Modifier les étiquettes des branches terminales du dendrogramme
labels(dendrogram) <- custom_labels

# Afficher le dendrogramme avec des branches proportionnelles aux distances
plot(dendrogram, hang = -1, main = "Dendrogramme des séquences", xlab = "Séquences", ylab = "Distance")

# bootstrap

fun <- function(x) as.phylo(hclust(dist.dna(x), "average")) # upgma() in phangorn
tree <- fun(dnabin_obj)

## get 100 bootstrap trees:
bstrees <- boot.phylo(tree, dnabin_obj, fun, trees = TRUE)$trees

## get proportions of each clade:
clad <- prop.clades(tree, bstrees, rooted = TRUE)

## get proportions of each bipartition:
boot <- prop.clades(tree, bstrees)
layout(1)
par(mar = rep(2, 4))
plot(tree, main = "Bipartition vs. Clade Support Values")
drawSupportOnEdges(boot)
nodelabels(clad)
legend("bottomleft", legend = c("Bipartitions", "Clades"), pch = 22,
       pt.bg = c("green", "lightblue"), pt.cex = 2.5)
