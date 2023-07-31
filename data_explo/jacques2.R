
# Charger le package seqinr
library(seqinr)
library(ape)
library(msa)
library(pegas)
library(stringr)


# Etape 1 : données brutes
## lecture du fichier issus de l'assignation
file<-"E:/Cours/SupAgro/Stage M1/Analyses/Romann/Data/tax_euka2207final3.rds"
data <- readRDS(file)
data<-as.data.frame(data)
dim(data)

# compter les reads de chacun des clusters
file<-"E:/Cours/SupAgro/Stage M1/Analyses/Romann/Data/seqtab_final3.rds"
reads <- readRDS(file)
reads<-as.data.frame(reads)
dim(reads)

# on fait la somme des colonnes
TOTreads<-colSums(reads[1:ncol(reads)])
length(TOTreads)

# on fait confiance à JF que les deux fichiers ont les séquences dans le meme ordre
data$Totreads<-TOTreads

# calcul du nombre de fois ou on voit un variant avec au moins 1 read
# on met 1 si il y a un read au moins , 0 sinon
R<-(reads>!0)*1 + (reads==0)*0
R<-as.data.frame(R)
R<-apply(R,2,as.numeric)
dim(R)

# comme la fonctions colSums a refusé de marcher, on en fait une
TOTR<-apply(R,2,sum)

length(TOTR)
hist(TOTR)

# on fait encore confiance à JF 
data$TOTACC<-TOTR

hist(data$TOTACC)

# on ne garde que blumeria
data<-data[grep("g__Blumeria",data[,6] ),]
dim(data)

data<-data[grep("s__graminis",data[,7] ),]
dim(data)

# on ne garde que les variants qui ont plus de 1 librairies
data<-data[which(data$TOTACC>1),]
dim(data)


## construire le fichier texte contenant les noms et les séquences
## important car on va devoir retrouver cette info dans le reseau d'haplotype
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

## Exporter l'objet dnabin au format fasta
write.dna(dnabin_obj, file = "bb.fasta", format = "fasta")

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

# on relit le fichier sans indel
mySeque_bb <- readDNAStringSet("Racoon_and_indel_free.fasta", format="fasta")

# on recupere les numeros de lignes pour pouvoir recuperer les comptages dans data
# car les séquences ne sont pas dans le même ordre a cause du fait qu'on les a alignees
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

# par curiosite Calculer la matrice de distances

D <- dist.dna(dnabin_obj)
D<-as.matrix(D)
dim(D)

heatmap(D)


 # Construction du réseau d'haplotypes
# voir cette doc http://ape-package.ird.fr/pegas/PlotHaploNet.pdf

# Créer un objet haplotype à partir de l'objet DNAbin filtré
# attention les noms ne sont pas dans le meme ordre que dans data...
noms <-   attributes(dnabin_obj)$dimnames[[1]]
noms

h <- haplotype(dnabin_obj, labels=1:23 )

## the indices of the individuals belonging to the 1st haplotype:
attr(h, "index")

net <- haploNet(h)

# fq donne l'effectif des individus dans chaque haplotype ici 1
fq <- attr(net, "freq")

# on met le nombre de reads dans fq en donnant les bonnes coordonnées dans data
fq<-log10(data$Totreads[as.numeric(noms)])

# on fait le reseau
plot(net, size = fq, cex=0.6, font=3, scale.ratio = 1,
     col='blue', bg="yellow", asp=1, lwd=1, lty=1,
     show.mutation =2,      fast=F)

# et on le retrace manuellement 
o <- replot()

## les dfferences entre haplotypes
Mutation_matrix<-matrix(NA, nrow = 23, ncol = 23)
for (i in 1:22) {
  for (j in (i+1):23) {
    Mutation_matrix[i,j]<-nrow(diffHaplo(h, a = i, b = j, strict = TRUE)) }}

# l'arbre en fonction des accessions touchées
# on met le nombre de reads dans fq
fq<-data$TOTACC[as.numeric(noms)]

plot(net, size = fq, cex=0.6, font=3, scale.ratio = 1,
     col='blue', bg="yellow", asp=1, lwd=1, lty=1,
     show.mutation =2,      fast=F)
o <- replot()

# la relation entre le nombre de reads et le nombre d'accessions touchées
plot(log10(data$Totreads), data$TOTACC)

