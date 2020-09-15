
################################################################################
#    &&&....&&&    % Project: Mapping KO ids to gene symbols  #
#  &&&&&&..&&&&&&  % Authors: Xiner Nie, Bo Li, Jingxin Tao, Hao He            #
#  &&&&&&&&&&&&&&  % Date: Aug. 18th, 2020                                     #
#   &&&&&&&&&&&&   %                                                           #
#     &&&&&&&&     % Environment: R version 3.6.0;                             #
#       &&&&       % Platform: x86_64-w64-mingw32/x64 (64-bit)                 #
#        &         %                                                           #
################################################################################

### ****************************************************************************
### code chunk number 10: Evaluation Criteria for data clustering.
### ****************************************************************************

download.file("http://rest.kegg.jp/list/ko", "ko.txt")

anno <- readLines("ko.txt")

strsplit(anno[1], "\t")[[1]][1]
strsplit(anno[1], "\t")[[1]][2]

strsplit(strsplit(anno[1], "\t")[[1]][2], "; ")[[1]]

strsplit(strsplit(strsplit(anno[1], "\t")[[1]][2], "; ")[[1]][2], " \\[")[[1]][1]

gsub("EC", "\\[EC", strsplit(strsplit(strsplit(anno[1], "\t")[[1]][2], "; ")[[1]][2], " \\[")[[1]][2])





# ============================
#     评估聚类               #
# ============================

# 引入fpc包(cluster.stats)
library(fpc)

# 引入包库(clara、fanny)
library(cluster)


#=====调用聚类算法=======================================================

# 确定簇心个数
cluster_num <- 3

# 读取数据
data <- read.csv("data.csv",header = T)

# 调用kmeans算法
km <- kmeans(data,cluster_num)

# 调用(clara)算法
cl <- clara(data,cluster_num)

# 调用模糊C-Means聚类算法
fan <- fanny(data,cluster_num) 

#=====调用聚类算法=======================================================


# 聚类评价统计量
km_stats <- cluster.stats(dist(data), km$cluster)
cl_stats <- cluster.stats(dist(data), cl$cluster)
fcm_stats <- cluster.stats(dist(data), fan$clustering)

# 信息数据框表化
info <- data.frame(
  Algorithm = c("KMeans", "Clara", "FCM"), 
  Silwidth = c(km_stats$avg.silwidth,cl_stats$avg.silwidth,fcm_stats$avg.silwidth),
  AverageWithin = c(km_stats$average.within, cl_stats$average.within, fcm_stats$average.within),
  averageBetween = c(km_stats$average.between, cl_stats$average.between, fcm_stats$average.between),
  ch = c(km_stats$ch, cl_stats$ch, fcm_stats$ch)
)

# 重命名字段
names(info)[2:5] <- c("轮廓系数","簇内平均距离","簇间平均聚类","Calinski和Harabasz指数")



library("biomaRt")                                                                                                                   
listMarts()                                                                                                                        
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")                                                                         
filters = listFilters(ensembl)                                                                                                        
entrezgene = ("105371919")             

genes <- getBM(filters="entrezgene_id", attributes=c("ensembl_gene_id","entrezgene_id"), values=entrezgene, mart=ensembl)                                                                                                                 
print(genes)


library(caret)
library(FSinR)

data(mtcars)


evaluator_1 <- filterEvaluator('determinationCoefficient')
evaluator_2 <- filterEvaluator('ReliefFeatureSetMeasure')





####### GA ##################

library(FSinR)

filter_evaluator <- filterEvaluator('binaryConsistency')
#Generates a filter function to be used as an evaluator in the feature selection proccess. 

search_method <- searchAlgorithm('geneticAlgorithm',
                                 list(popSize = 200, 
                                      pcrossover = 0.8,
                                      pmutation = 0.1, 
                                      maxiter = 10,
                                      run = 10))

load("benchmark data_GSE89.RData")

#list()……modify parameters

res.GA <- featureSelection(input, 'sam.lab', search_method, filter_evaluator)

#data_x,data_y

GA.features<-res.GA$bestFeatures

GADEG <- t(GA.features)

GADEG <- apply(GADEG, 1, sum)

GADEG <- names(GADEG)[rev(order(GADEG))]

deg.GA <- GADEG[1:10]

deg.GA <- na.omit(deg.GA)


if (length(deg.GA) < 10) {
  
  deg.cer <- deg.GA
  
} else {
  
  deg.cer <- deg.GA[1:10]
  
}



sam.iris <- input[, c("sam.lab", deg.cer)]

X <- sam.iris[, names(sam.iris) != "sam.lab"]

y <- as.character(sam.iris$sam.lab)

folds <- 5

test.fold <- split(sample(1:length(y)), 1:folds) # ignore warning


all.pred.tables <- lapply(1:folds, function(i) {
  
  test.id <- test.fold[[i]]
  
  X.train <- X[-test.id, ]
  
  y.train <- as.factor(y[-test.id])
  
  model <- svm(X.train, y.train, kernel = "radial", prob = TRUE, cross = 5) # some tuning may be needed
  
  predict.test <- predict(model, X[test.id, ], prob = TRUE)
  
  prob.benign <- attr(predict.test, "probabilities")[, 2]
  
  data.frame(y.test = y[test.id], y.pred = prob.benign) # returning this
  
})

full.pred.table <- do.call(rbind, all.pred.tables)

res.roc <- roc(full.pred.table$y.test, 
               full.pred.table$y.pred, 
               plot = TRUE, 
               legacy.axes = TRUE)
GA.auc.value <- auc(res.roc)









ind<-sample(2,nrow(input),replace = TRUE, prob=c(0.7,0.3))

train <-input[ind==1,]
test <-input[ind==2,]

set.seed(100)

RF <- randomForest(sam.lab ~ ., data = train, mtry = 2, ntree=100, importance=TRUE)

deg.RF<-as.data.frame(RF$importance)

deg.RF<-deg.RF[rev(order(deg.RF$MeanDecreaseAccuracy)),]

deg.RF<-rownames(deg.RF)[1:10]

set.seed(5)

deg.cer<-ifelse(length(deg.RF)<10, deg.RF, deg.RF[1:10])

sam.iris <- input[, c("sam.lab", deg.cer)]

X <- sam.iris[, names(sam.iris) != "sam.lab"]

y <- as.character(sam.iris$sam.lab)

test.fold <- split(sample(1:length(y)), 1:folds) # ignore warning

all.pred.tables <- lapply(1:folds, function(i) {
  
  test.id <- test.fold[[i]]
  
  X.train <- X[-test.id, ]
  
  y.train <- as.factor(y[-test.id])
  
  model <- svm(X.train, y.train, kernel = "radial", prob = TRUE, cross = 5) # some tuning may be needed
  
  predict.test <- predict(model, X[test.id, ], prob = TRUE)
  
  prob.benign <- attr(predict.test, "probabilities")[, 2]
  
  data.frame(y.test = y[test.id], y.pred = prob.benign) # returning this
  
})

full.pred.table <- do.call(rbind, all.pred.tables)

res.roc <- roc(full.pred.table$y.test, 
               full.pred.table$y.pred, 
               plot = TRUE, 
               legacy.axes = TRUE)
RF.auc.value <- auc(res.roc)

####### end ##################






