library(readr)
library(HGNChelper) #Used version 0.8.14, May 2024 updated gene names
library(dplyr)
library(caret)
library(preprocessCore)
library(rrcov)
library(Hmisc)
library(factoextra)
library(M3C)
library(GEOquery)
library(biomaRt)
library(tidyverse)
library(sva)


########Preprocess the Immport dataset##########

gene_pbmc_normalized <- read_csv("gene_pbmc_normalized.csv") #Downloaded from 10K immunomes
save(gene_pbmc_normalized, file = "gene_pbmc_normalized.RData") #Saved as R object

madat <- data.frame(gene_pbmc_normalized[,9:ncol(gene_pbmc_normalized)])
(row.names(madat) <- gene_pbmc_normalized$subject_accession)

checkedgenenames <- checkGeneSymbols(colnames(madat), unmapped.as.na = TRUE, map = NULL, species = "human")
checkedgenenames$Suggested.Symbol.2 <- sapply(strsplit(checkedgenenames$Suggested.Symbol, "/\\s*"), tail, 1)

madat.2.df <- as.data.frame(t(madat))
madat.2.df$ID.orig <- row.names(madat.2.df)

madat.2.df$ID.update <- checkedgenenames$Suggested.Symbol.2[match(madat.2.df$ID.orig, checkedgenenames$x)]
madat.3  <- subset(madat.2.df, !is.na(madat.2.df$ID.update))

sum(is.na(madat.2.df$ID.update)) #1527 

dim(madat.2.df) #[1] 19254   167
dim(madat.3) #[1] 17727      167

sum(duplicated(madat.3$ID.update)) #74

madat.4 <- madat.3 %>% group_by(ID.update) %>%
  summarise(across(where(is.numeric), median))

madat.4.df <- as.data.frame(madat.4)
row.names(madat.4.df) <- madat.4.df$ID.update

madat.5 <- subset(madat.4.df, select = -ID.update) 

transf.params.madat.cent <- preProcess(t(madat.5), method="center")
madat.cent <- t(predict(transf.params.madat.cent, t(madat.5)))
madat.cent.QN <- normalize.quantiles(madat.cent, keep.names = T)

set.seed(1)
pcaHub <- PcaHubert(t(madat.cent.QN), mcd=F, scale = T)
(outliers.pca <- names(which(pcaHub$flag=='FALSE')))
madat.df <- madat.cent.QN[,colnames(madat.cent.QN) %nin% outliers.pca]
dim(madat.cent.QN.out)

fviz_pca_ind(prcomp(t(madat.df)), habillage = subj.age.sex$study_accession, addEllipses=T)

set.seed(1)
tsne(madat.df, labels=subj.age.sex$study_accession)

boxplot(madat.df)

plot(density(t(madat.df)),
     frame = FALSE, col = "blue",main = "Density plot")


subj.age.sex <- subset(gene_pbmc_normalized, subject_accession %in% colnames(madat.df))

subj.age.sex$age.sex <- paste(floor(subj.age.sex$age), substring(subj.age.sex$gender,1,1), 
                              subj.age.sex$subject_accession, sep="_")

subj.age.sex.ord <- subj.age.sex[order(subj.age.sex$age.sex),]

colnames(madat.df) <- subj.age.sex.ord$age.sex[match(colnames(madat.df), subj.age.sex.ord$subject_accession)]
madat.df.sorted <- madat.df[,order(colnames(madat.df))]
dim(madat.df.sorted)
colnames(madat.df.sorted)

dim(madat.df.sorted)# 17653   129

mean(subj.age.sex.ord$age) #45.81736
sd(subj.age.sex.ord$age) # 20.85259

table(subj.age.sex.ord$gender)
table(subj.age.sex.ord$gender)/ncol(madat.df.sorted)*100

table(subj.age.sex.ord$race)
table(subj.age.sex.ord$race)/ncol(madat.df.sorted)*100


########Preprocess the GSE3365 External Validation set##########

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

tables <- listAttributes(mart)
tables[grep('affy', tables[,1]),]

affyHGU133p2 <- getBM(attributes = c('affy_hg_u133_plus_2',"hgnc_symbol"), mart = mart, useCache = FALSE)

gse.GSE3365 <- getGEO("GSE3365", GSEMatrix = TRUE, getGPL = T)
show(gse.GSE3365)
#GPL96	[HG-U133A] Affymetrix Human Genome U133A Array

pData.GSE3365 <- pData(gse.GSE3365$GSE3365_series_matrix.txt.gz)

pData.GSE3365.normalonly <- pData.GSE3365[grep("^Normal", pData.GSE3365$title),]

pData.GSE3365.normalonly$age <- parse_number(pData.GSE3365.normalonly$description.2)
pData.GSE3365.normalonly$sex <- sub(".*: ", "", pData.GSE3365.normalonly$description.4)
pData.GSE3365.normalonly$race <- sub(".*: ", "", pData.GSE3365.normalonly$description.3)

table(pData.GSE3365.normalonly$age)
#25 27 28 31 33 35 36 39 40 42 44 45 46 47 48 49 50 51 52 53 54 56 57 60 
#1  1  1  3  1  1  1  1  3  6  1  2  1  4  1  1  1  2  1  2  4  1  1  1 
mean(pData.GSE3365.normalonly$age)
sd(pData.GSE3365.normalonly$age)

table(pData.GSE3365.normalonly$sex)
#female   male 
#18       24 

table(pData.GSE3365.normalonly$race)

eset.GSE3365 <- exprs(gse.GSE3365$GSE3365_series_matrix.txt.gz)
eset.GSE3365.df <- data.frame(eset.GSE3365)[, pData.GSE3365.normalonly$geo_accession]
eset.GSE3365.df_log <- log2(eset.GSE3365.df)
dim(eset.GSE3365.df_log) #22283    42

length(intersect(affyHGU133p2$affy_hg_u133_plus_2, row.names(eset.GSE3365.df_log))) #20258

rownamesof_GSE3365 <- checkGeneSymbols(affyHGU133p2$hgnc_symbol, unmapped.as.na = TRUE, map = NULL, species = "human")
rownamesof_GSE3365$Suggested.Symbol.2 <- sapply(strsplit(rownamesof_GSE3365$Suggested.Symbol, "/\\s*"), tail, 1)
rownamesof_GSE3365 <- rownamesof_GSE3365[!duplicated(rownamesof_GSE3365), ]

eset.GSE3365.df_log$affy_hg_u133_plus_2 <- row.names(eset.GSE3365.df_log)

eset.GSE3365.df_2 <- merge(eset.GSE3365.df_log, affyHGU133p2, by="affy_hg_u133_plus_2", all.x = T)

dim(eset.GSE3365.df) #22283
dim(eset.GSE3365.df_2) #25649
dim(rownamesof_GSE3365) #41063

eset.GSE3365.df_3 <- subset(eset.GSE3365.df_2, !is.na(hgnc_symbol))
dim(eset.GSE3365.df_3) #23624    44

sum(is.na(eset.GSE3365.df_2$hgnc_symbol)) #2025

eset.GSE3365.df_3$ID.update <- rownamesof_GSE3365$Suggested.Symbol.2[match(eset.GSE3365.df_3$hgnc_symbol, rownamesof_GSE3365$x)]
eset.GSE3365.df_4  <- subset(eset.GSE3365.df_3, !is.na(eset.GSE3365.df_3$ID.update))
dim(eset.GSE3365.df_4) #22325    45

sum(duplicated(eset.GSE3365.df_3$ID.update)) #9465

eset.GSE3365.df_4$ID.update[duplicated(eset.GSE3365.df_4$ID.update)==TRUE]

eset.GSE3365.df_5 <- eset.GSE3365.df_4 %>% group_by(ID.update) %>%
  summarise(across(where(is.numeric), median))

dim(eset.GSE3365.df_4) #22325
dim(eset.GSE3365.df_5) #14158

eset.GSE3365.df_6 <- as.data.frame(eset.GSE3365.df_5)

row.names(eset.GSE3365.df_6) <- eset.GSE3365.df_6$ID.update

eset.GSE3365.df_7 <- subset(eset.GSE3365.df_6, select = -ID.update) 

dim(eset.GSE3365.df_7) #14158    42

plot(density(t(eset.GSE3365.df_7)),
     frame = FALSE, col = "blue",main = "Density plot")

boxplot(eset.GSE3365.df_7)

(pData.GSE3365.normalonly$replacecolnames <- paste0(pData.GSE3365.normalonly$age, "_", 
                                                    substr(pData.GSE3365.normalonly$sex, 1, 1), "_",
                                                    pData.GSE3365.normalonly$geo_accession))

names(eset.GSE3365.df_7) <- pData.GSE3365.normalonly$replacecolnames[match(names(eset.GSE3365.df_7), pData.GSE3365.normalonly$geo_accession)]
names(eset.GSE3365.df_7)

GSE3365.df.ord <- eset.GSE3365.df_7[,order(names(eset.GSE3365.df_7))]
colnames(GSE3365.df.ord)


########## Subset to use common genes only ##########

common_genes <- intersect(row.names(madat.df.sorted), row.names(GSE3365.df.ord))
length(common_genes) ##12417

madat.df.sorted.filt <- madat.df.sorted[common_genes,]
dim(madat.df.sorted.filt)

GSE3365.filt <- GSE3365.df.ord[row.names(GSE3365.df.ord) %in% common_genes,]


#######randomly split data into training & test set #########

madat_ages <- parse_number(colnames(madat.df.sorted.filt))

set.seed(12345)
trainIndex <- createDataPartition(y=madat_ages, p = .8, times=1)

madat.train <- madat.df.sorted.filt[,trainIndex$Resample1]
madat.test <- madat.df.sorted.filt[,-trainIndex$Resample1]

#######ensure the training & test set demographics are comparable#########

mean(parse_number(colnames(madat.train))) #45.34286
sd(parse_number(colnames(madat.train))) #20.55239

(Age.vec.train <- parse_number(colnames(madat.train)))
(Age.vec.test <- parse_number(colnames(madat.test)))

range(Age.vec.train) #21 90
range(Age.vec.test) #22 87

mean(Age.vec.train) #45.34286
mean(Age.vec.test) #47.58333
sd(Age.vec.train) #20.55239
sd(Age.vec.test) #22.42071

t.test(Age.vec.train, Age.vec.test, var.equal =T) #p=0.6365


tab.sex <- matrix(c(sum(str_sub(colnames(madat.train), 4, 4)=="F"),
                     sum(str_sub(colnames(madat.train), 4, 4)=="M"),
                     sum(str_sub(colnames(madat.test), 4, 4)=="F"),
                     sum(str_sub(colnames(madat.test), 4, 4)=="M")), 
                   ncol=2, byrow=T)

colnames(tab.sex) <- c("F", "M")
row.names(tab.sex) <- c("Train", "Test")

tab.sex 
#F  M
#Train 61 44
#Test  16  8

chisq.test(tab.sex) #p=0.588

table(subj.age.sex.ord$race)
table(subj.age.sex.ord$race)/ncol(madat.df.sorted)*100

train.race <- subset(subj.age.sex.ord, subject_accession %in% str_sub(colnames(madat.train), 6, 14))
test.race <- subset(subj.age.sex.ord, subject_accession %in% str_sub(colnames(madat.test), 6, 14))

train.race$race.2 <- ifelse(train.race$race=="Unknown"|train.race$race=="Other", "Other/unk",
                            train.race$race)

test.race$race.2 <- ifelse(test.race$race=="Unknown"|test.race$race=="Other", "Other/unk",
                            test.race$race)
table(train.race$race.2)
table(test.race$race.2)

table(train.race$race.2)/ncol(madat.train)*100
table(test.race$race.2)/ncol(madat.test)*100

race.tab <- as.matrix(rbind(as.vector(table(train.race$race.2)), as.vector(table(test.race$race.2))))
colnames(race.tab) <- c("Asian", "Black/AfricanAmerican", "Other/unk", "White")
row.names(race.tab) <- c("train", "test")

fisher.test(race.tab) #0.7215


#### External validation set - batch removal & quantile normalization for comparability with training set####

GSE3365.m <- cbind(madat.train, GSE3365.filt)
GSE3365.batch <- c(rep("trainset", ncol(madat.train)), rep("GSE3365", ncol(GSE3365.filt)))
dim(GSE3365.m) 

GSE3365.m.combat <- ComBat(dat=GSE3365.m, batch=GSE3365.batch, ref.batch = "trainset", 
                           mod=NULL, par.prior=F, prior.plots=FALSE)
sum(is.na(GSE3365.m.combat)) #0

targetQN <- normalize.quantiles.determine.target(as.matrix(madat.train), target.length=NULL,subset=NULL)
GSE3365.QN <-normalize.quantiles.use.target(x=GSE3365.m.combat, target=targetQN, copy=TRUE,subset=NULL)
colnames(GSE3365.QN) <- colnames(GSE3365.m)
row.names(GSE3365.QN) <- row.names(GSE3365.m)

GSE3365.nomalized <- GSE3365.QN[, GSE3365.batch=="GSE3365"]


#### standardize training, test, external validation sets according to the training set's distribution ####

transf.params.sc <- preProcess(t(madat.train), method=c("center", "scale"))
trainset <- t(predict(transf.params.sc, t(madat.train)))
testset <- t(predict(transf.params.sc, t(madat.test)))
extvalset <- t(predict(transf.params.sc, t(GSE3365.nomalized)))

summary(as.vector(data.matrix(trainset)))
sd(as.vector(data.matrix(trainset)))

summary(as.vector(data.matrix(testset)))
sd(as.vector(data.matrix(testset)))

summary(as.vector(data.matrix(GSE3365.ext)))
sd(as.vector(data.matrix(GSE3365.ext)))

######
write.csv(trainset, file="trainset.csv")
write.csv(testset, file="testset.csv")
write.csv(extvalset, file="extvalset.csv")
