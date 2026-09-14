#BiocManager::install(c("affy","GEOquery","Biobase","hgu133plus2.db","affyPLM","hgu133acdf","hgu133a.db","hgu133plus2cdf","genefilter")) 

library(affy)
library(GEOquery)
library(Biobase)
library(AnnotationDbi)
#devtools::install_github("https://github.com/bioc/simpleaffy")
library(simpleaffy)
library(affyPLM)
library(hgu133plus2.db)
#library(hgu133acdf)
#library(hgu133a.db)
library(hgu133plus2cdf)
library(genefilter)

## Call the data in R
gse<- list.celfiles("GSE45827_RAW", full.names=T)
class(gse)
gse

affy.data = ReadAffy(filenames=gse)
class(affy.data)
affy.data

#get pheno data
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 4)
gset <- getGEO(GEO="GSE45827",GSEMatrix =TRUE)
class(gset)
names(gset)
gset
data.gse <- exprs(gset[[1]])
class(data.gse)
dim(data.gse)
head(data.gse)

#BOXPLOT before pre-processing
#dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE45827", '/', annotation(gset[[1]]), " Dataset", sep ='')
boxplot(exprs(gset[[1]]), boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)

pheno <- pData(phenoData(gset[[1]]))
varLabels(phenoData(gset[[1]]))
dim(pheno)
class(pheno)

#menyimpan data pheno, sesuaikan dengan laptop masing2 cara menyimpan nya
#ganti tulisan 'N/A' dengan 'None (normal)'
#pada variabel/kolom 'diagnosis:ch1'

write.csv(pheno,"phe.csv")

#memanggil data pheno setelah revisi, sesuaikan dengan laptop masing2
phenoa <- read.csv("phe.csv")
class(phenoa)
phenoc<-data.frame(phenoa)
phenob<-phenoc[,-1]
read.csv
class(phenob)
dim(phenob)
rownames(phenob)<-rownames(pheno)
colnames(phenob)
#Take alook the Pheno
table(phenob$`diagnosis.ch1`)
table(phenob$`tumor.subtype.ch1`)

#which(phenob$`tumor.subtype.ch1`=="Bcc")
#pheno$`diagnosis:ch1`[pheno$`diagnosis:ch1` == 'NA'] <- 'Breast cancer'


#remove subgroup SHH OUTLIER and U
#phenob<-pheno
#phenob
#dim(phenob)

#table(phenob$`geo_accession`)


#pie chart diagnosis
des <-table(phenob$`diagnosis.ch1`)
des
persen <- round(des/sum(des)*100)
des <- as.data.frame(des)
lbls <- paste(des$Var1,'-',persen, '%', sep='')
#quartz()
pie(des$Freq, label= lbls, col=c('blue',"pink"))

#pie chart tumor subtype
des <-table(phenob$`tumor.subtype.ch1`)
des
persen <- round(des/sum(des)*100)
des <- as.data.frame(des)
lbls <- paste(des$Var1,'-',persen, '%', sep='')
#quartz()
pie(des$Freq, label= lbls, col=c('blue',"pink","green","black","maroon","yellow"))


#quartz()
#pie(des$Freq, label= lbls, col=c(1:7))

# New data based on phenob
gseb<- gse
gseb
affy.datab = ReadAffy(filenames=gseb)

affy.datab
#pre processing data
eset.dChip=threestep(affy.datab,background.method = "RMA.2", normalize.method="quantile",summary.method="median.polish")
class(eset.dChip)
dim(eset.dChip)
eset.dChip

Ekspres <- exprs(eset.dChip)
class(Ekspres)
dim(Ekspres)
head(Ekspres)

#BOXPLOT after pre-processing
#dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(phenob$geo_accession))/2),4,2,1))
title <- paste ("GSE45827", '/',annotation(gset[[1]]), " Dataset", sep ='')
boxplot(Ekspres, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)


#filtering
filterdataa <- nsFilter(eset.dChip, require.entrez =T,var.func = IQR, 
                        remove.dupEntrez = T,var.cutoff = 0.5, feature.exclude = "^AFFX")
log <-filterdataa$filter.log
eset <- filterdataa$eset
featureNames(eset) <- make.names(featureNames(eset))

#Filtering results (eset)
dim(eset) 
#head(eset)
class(eset)
eset

# create the matrix/data frame
databaru <- exprs(eset)
dim(databaru)
class(databaru)
head(databaru)

#View(databaru)

#There are 4 Subgroup classes: G4, G3, SHH and WNT
phenob$`diagnosis.ch1`
datacl <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
datacl

#Filtering with multtest
# devtools::install_github("https://github.com/bioc/multtest")
library(multtest, verbose = FALSE)
datattest <- mt.teststat(databaru,datacl,test="t") 
class(datattest)
length(datattest)
qqnorm(datattest)
qqline(datattest)

#Adjusted p-value  
rawp = 2 * (1 - pnorm(abs(datattest))) 
prosedur = c("Bonferroni", "Holm", "Hochberg", "BH", "BY") 
adjusted = mt.rawp2adjp(rawp, prosedur) 
data <- adjusted$adjp[,] 
data1 <- data[order(adjusted$index), ] 
head(data1)
dim(data1)

#Take the Bonferroni adjusted p-value data column
ffs <- data1[,2] 
class(ffs)
length(ffs)
#ffs[1 : 10]
#View(ffs)

#Adjusted rawp 
datarawp <- data.frame(databaru, ffs) 
row.names(datarawp) <- row.names(databaru)
class(datarawp)
#head(datarawp)
dim(datarawp) 

library(dplyr) 
#datadatarawpfilter <- filter(datarawp, ffs < 0.0005) 
#class(datarawpfilter) dimensi = 96 39
datarawpfilterfinal <- subset(datarawp, ffs < 0.0000001)
#rownames(datarawpfilterfinal)
class(datarawpfilterfinal)
dim(datarawpfilterfinal) # 21 x 179
head(datarawpfilterfinal)

## Define a new data after filtering
datadef <- datarawpfilterfinal[,1:178]
head(datadef)
dim(datadef)
summary(datadef)
colnames(datadef)

#Preparing data for classification
#install.packages(c("e1071","pROC"))
library(e1071)
library(pROC)

data   = as.data.frame (t((datadef)))
dim(data)
#head(data)
#dataY = as.factor(datacl) 
dataY = datacl
dataY

#if it is possible add the demographical variables from phenob data: `age:ch1`, `Sex:ch1`, `m stage:ch1`, `subgroup:ch1`, `ethnic:ch1`
datause = as.data.frame(cbind(data,dataY))
#datause
dim(datause)
#head(datause)
dim(data)
write.table(datause,"C:/Users/Lenovo/Downloads/datause.txt")
write.csv(datause,"C:/Users/Lenovo/Downloads/datause2.csv")

#Classification with SVM kernel Linier
set.seed(123)
rasio = 0.8
train = sample(length(dataY), size = floor(rasio*length(dataY)))
datatrain<- datause[train,]
dim(datatrain)
#datatrain
datatest <- datause[-train,]
dim(datatest)

#tuning (best parameter)
tunaslin <- tune(svm, dataY~. ,data = datatrain, kernel="linear",types = "C-clasification",ranges= list( cost = c(0.1,0.01, 0.001,1 , 10 , 100)))
summary(tunaslin)


# ========================================================================================================
# ========================================================================================================

  #SVM

# ========================================================================================================
# ========================================================================================================
#Build the model with the best tuning parameter

model <- svm(dataY~. , datatrain, kernel = "linear", cost=0.1,scale=F, types = "C-clasification",decision.value=T)
predictionste <- predict(model, datatest)
predictionste <- ifelse(predictionste < 1, 0, 1) 
table(predictionste,datatest$dataY)
mean(predictionste == datatest$dataY)
predictionstr <- predict(model, datatrain) 
predictionstr <- ifelse(predictionstr < 1, 0, 1) 
table(predictionstr,datatrain$dataY) 
mean(predictionstr == datatrain$dataY)

#plot(model,datatrain, datatrain$gender~datatrain$ethnic, slice = list(ethnic = 4, gender = 4))
#weighted
w <- t(model$coefs) %*% model$SV        
w <- apply(w, 2, function(v){sqrt(sum(v^2))})   
w <- sort(w, decreasing = T)
x<- as.data.frame(w)
#View(x)
all(colnames(w)==colnames(datatrain))

#See the genes: gene yang paling berpengaruh
id <- substring(as.character(head(rownames(x), n=25)),2) 
id
#View(id)

#GENE Profiling Gene Ontology
AnnotationDbi::select(hgu133plus2.db, id, c("SYMBOL","GENENAME","ENTREZID","ONTOLOGY"), "PROBEID")

#gene-gene interaction
datagambar<-data
dim(datagambar)
library(qgraph)
#pdf("paper_jurnal.pdf")
qgraph(cor(datagambar[,1:21]), layout="spring", posCol="darkgreen", negCol="darkmagenta")

#Heatmap
idb<-as.character(tail(rownames(x),n=21 ))
idb
dim(datarawpfilterfinal)
#View(id)
datakuh<-datarawpfilterfinal[,-156]
dim(datakuh)
colnames(datakuh)
head(datakuh)
aaa<-rownames(datakuh)
aaa
length(aaa)
idi<-which(is.element(aaa,idb)==TRUE)
idi
#1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21
rna<-aaa[idi]
rna
dataheatmap<-datakuh[idi,]
dim(dataheatmap)
#dataheatmap
row.names(dataheatmap)<-rna
colnames(dataheatmap)<-datacl
head(dataheatmap)
colnames(dataheatmap)
class(dataheatmap)
summary(dataheatmap)
datacl
#heatmap(as.matrix(dataheatmap),scale="none")
#legend(x = "topright", legend = c("1: Vehicle","2: ONECUT2 overexpression construct","3: shRNA scramble","4: shONECUT2"), cex = 1)


#heatmap(as.matrix(dataheatmap), col=topo.colors(100))
library(gplots)
library(RColorBrewer)
color.map <- function(datacl) { if (datacl=="1") "#FF0000" else  "#0000FF" }
patientcolors <- unlist(lapply(datacl, color.map))
patientcolors
colMain <- colorRampPalette(brewer.pal(8, "Greens"))(25)
heatmap.2(as.matrix(dataheatmap), col=topo.colors(75), scale="none", ColSideColors=patientcolors,
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)

heatmap.2(as.matrix(dataheatmap), col=colMain , scale="none", ColSideColors=patientcolors,
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)

heatmap.2(as.matrix(dataheatmap), col=redgreen(75) , scale="none", ColSideColors=patientcolors,
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
#scale could be row, column or none
#One subtle point in the previous examples is that the heatmap function has automatically 
#scaled the colours for each row (i.e. each gene has been individually normalised across patients). 
#This can be disabled using scale="none", which you might want to do 
#if you have already done your own normalisation (or this may not be appropriate for your data):
#Jadi, kl sdh dinormalisasi lebih baik menggunakan "none"


# ========================================================================================================
# ========================================================================================================

  #RANDOM FOREST
# ========================================================================================================
# ========================================================================================================

#libraries/packages yang diperlukan 
library(ROSE)
library(randomForest)
library(caret)

#memastikan bahwa column/variabel dalam data set train - test sama
#colnames(datatrain)<-colnames(datatest)

#tuning parameters
set.seed(123)
randomfrst<-tuneRF(datatrain[,-22], datatrain[,22],stepFactor=2, improve=0.05,trace=TRUE,plot=TRUE,doBest=FALSE)
randomfrst

#gunakan mtry yang sesuai yaitu yang OOB nya paling kecil
#Dicoba dgn beberapa mtry
Dataranfrst1<-randomForest(datatrain[,-22], datatrain[,22],importance=T,ntree=50,mtry=7)
Dataranfrst2<-randomForest(datatrain[,-22], datatrain[,22],importance=T,ntree=100,mtry=7)
Dataranfrst3<-randomForest(datatrain[,-22], datatrain[,22],importance=T,ntree=200,mtry=7)
Dataranfrst4<-randomForest(datatrain[,-22], datatrain[,22],importance=T,ntree=300,mtry=7)
Dataranfrst5<-randomForest(datatrain[,-22], datatrain[,22],importance=T,ntree=400,mtry=7)
Dataranfrst6<-randomForest(datatrain[,-22], datatrain[,22],importance=T,ntree=500,mtry=7)
Dataranfrst1
Dataranfrst2
Dataranfrst3
Dataranfrst4
Dataranfrst5
Dataranfrst6

plot(Dataranfrst1)
plot(Dataranfrst2)
plot(Dataranfrst3)
plot(Dataranfrst4)
plot(Dataranfrst5)
plot(Dataranfrst6)

#Dipilih model dengan akurasi terbaik, kalau disini misal menggunakan:
#Dataranfrst6

Predtrainrf<-predict(Dataranfrst3,datatrain[,-22])
Predtrainrf <- ifelse(Predtrainrf < 1, 0, 1) 
confusionMatrix(as.factor(Predtrainrf),as.factor(datatrain[,22]))
Predtestrf<-predict(Dataranfrst3,datatest[,-22])
Predtestrf <- ifelse(Predtestrf < 1, 0, 1)
confusionMatrix(as.factor(Predtestrf),as.factor(datatest[,22]))

#menghitung auc nya
library(pROC)
foraucRF<-as.numeric(predict(Dataranfrst3,datatest[,-22]))
aucRF<-multiclass.roc(datatest$dataY,foraucRF)
aucRF$auc

#plotting importance variabel terbaik 
#View(importance(Dataranfrst2),"GENENAME")
varImpPlot(Dataranfrst3,pch=19,main="importance variable optimum")

