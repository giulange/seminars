# https://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm

## WRITE HERE THE R VERSION TO USE: R???

## *** INSTALL PACKAGES ***
install.packages("tree",dependencies=T)
install.packages("randomForest",dependencies=T)
install.packages("rpart",dependencies=T)
install.packages("partykit",dependencies=T)
install.packages("e1071",dependencies=T)
install.packages("party",dependencies=T)
install.packages("caret",dependencies=T)
# manual install from local zip:
install.packages("/home/giuliano/R/strucchange_1.5-1.tar.gz",repos=NULL)
install.packages("/home/giuliano/R/mvtnorm_1.0-6.tar.gz",repos=NULL)
install.packages("/home/giuliano/R/coin_1.1-3.tar.gz",repos=NULL)
install.packages("/home/giuliano/R/modeltools_0.2-21.tar.gz",repos=NULL)
install.packages("/home/giuliano/R/party_1.2-2.tar.gz",repos=NULL)
install.packages("/home/giuliano/R/partykit_1.1-1.tar.gz",repos=NULL)
install.packages("/home/giuliano/R/ModelMetrics_1.1.0.tar.gz",repos=NULL)
install.packages("/home/giuliano/R/pbkrtest_0.4-7.tar.gz",repos=NULL)
install.packages("/home/giuliano/R/minqa_1.2.3.tar.gz",repos=NULL)
install.packages("/home/giuliano/R/nloptr_1.0.4.tar.gz",repos=NULL)
install.packages("/home/giuliano/R/RcppEigen_0.3.2.9.0.tar.gz",repos=NULL)
install.packages("/home/giuliano/R/pbkrtest_0.4-6.tar.gz",repos=NULL)
install.packages("/home/giuliano/R/lme4_1.1-11.tar.gz",repos=NULL)
install.packages("/home/giuliano/R/car_2.1-3.tar.gz",repos=NULL)
install.packages("/home/giuliano/R/randomForest_4.6-10.tar.gz",repos=NULL)
install.packages("/home/giuliano/R/e1071_1.6-7.tar.gz",repos=NULL)
install.packages("/home/giuliano/R/caret_6.0-73.tar.gz",repos=NULL)


## *** PRE ***
# Load libraries
require(tree)
require(randomForest)
require(rpart)
require(partykit)
require(e1071)
require(party)
library(caret)
# Help on ramdonForest package and function
library(help=randomForest)
help(randomForest)


## *** INPUT ***
# Read data
D<-read.table(file="~/R/randomForests/D_cl-o-r_corr.txt",header = T, sep=",")
names(D)
## *** DATA RESAMPLING ***
# Both development and validation samples have 
# similar target variable distribution:
sample.ind <- sample(2, 
                     nrow(D),
                     replace = T,
                     prob = c(0.7,0.3))
tr <- D[sample.ind==1,]
te <- D[sample.ind==2,]
# check the type of response variable:
class(tr$Al05Fe)

## *** Make Formula ***
varNames <- names(tr)
# Exclude ID or Response variable
varNames <- varNames[!varNames %in% c("Al05Fe")]
# add + sign between exploratory variables
varNames1 <- paste(varNames, collapse = "+")
# Add response variable and convert to a formula object
formula.all <- as.formula(paste("Al05Fe", varNames1, sep = " ~ "))

## *** BUILD TREE ***
mod.tree <- tree(formula.all, data=tr)
plot(mod.tree)
text(mod.tree, cex=.75)
deciles <- quantile(tr$Al05Fe, 0:10/10)
cut    <- cut(tr$Al05Fe, deciles, include.lowest=TRUE)
plot(tr$Easting_32632, tr$Northing_32632, col=grey(10:2/11)[cut], pch=20, xlab="Easting",ylab="Northing")
partition.tree(mod.tree, ordvars=c("Easting_32632","Northing_32632"), add=TRUE)
summary(mod.tree)
te$predicted.tree <- predict(mod.tree,te)
plot(te$predicted.tree,te$Al05Fe)
cor(te$predicted.tree,te$Al05Fe,use="pairwise.complete.obs",method="pearson")
# cross-validation & pruning
cv.model <- cv.tree(mod.tree)
plot(cv.model)
cv.model$dev # gives the deviance for each K (small is better)
best.size <- cv.model$size[which(cv.model$dev==min(cv.model$dev))] # which size is better?
best.size
# let's refit the tree model (the number of leafs will be no more than best.size)
cv.model.pruned <- prune(mod.tree, best=best.size)
summary(cv.model.pruned)
te$predicted.cv.model.pruned <- predict(cv.model.pruned,te)
plot(te$predicted.cv.model.pruned,te$Al05Fe)
cor(te$predicted.cv.model.pruned,te$Al05Fe,use="pairwise.complete.obs",method="pearson")

mod.tree2 <- tree(formula.all, data=tr, mindev=0.001)
plot(mod.tree2)
text(mod.tree2, cex=.75)
summary(mod.tree2)
te$predicted.tree2 <- predict(mod.tree2,te)
plot(te$predicted.tree2,te$Al05Fe)
cor(te$predicted.tree2,te$Al05Fe,use="pairwise.complete.obs",method="pearson")

mod.tree.pruned <- prune.tree(mod.tree, best=4)
plot(mod.tree.pruned)
text(mod.tree.pruned)
te$predicted.tree.pruned <- predict(mod.tree.pruned,te)
plot(te$predicted.tree.pruned,te$Al05Fe)
cor(te$predicted.tree.pruned,te$Al05Fe,use="pairwise.complete.obs",method="pearson")

# rpart:
rpart.tree <- rpart(formula.all, data=tr)
plot(rpart.tree, uniform=TRUE, branch=0.6, margin=0.05)
text(rpart.tree, all=TRUE, use.n=TRUE)
title("Training Set's Regression Tree")
predictions <- predict(rpart.tree, te)
table(te$Al05Fe, predictions)
prune.rpart.tree <- prune(rpart.tree, cp=0.02) # pruning the tree
plot(prune.rpart.tree, uniform=TRUE, branch=0.6)
text(prune.rpart.tree, all=TRUE, use.n=TRUE)
# lmat | I should understand "lmat" and better define it!!
lmat <- matrix(c(0,1,2,
                 1,0,100,
                 2,100,0), ncol = 3)
lmat

rpart.tree <- rpart(formula.all, data=tr, parms = list(loss = lmat))
te$predicted.rparttree <- predict(rpart.tree, te)
table(te$predicted.rparttree, predictions)
plot(rpart.tree)
text(rpart.tree)
plot(te$predicted.rparttree,te$Al05Fe)
cor(te$predicted.rparttree,te$Al05Fe,use="pairwise.complete.obs",method="pearson")

## Define a plotting function with decent defaults
plot.rpart.obj <- function(rpart.obj, font.size = 0.8) {
    ## plot decision tree
    plot(rpart.obj,
         uniform   = T,    # if 'TRUE', uniform vertical spacing of the nodes is used
         branch    = 1,    # controls the shape of the branches from parent to child node
         compress  = F,    # if 'FALSE', the leaf nodes will be at the horizontal plot
         nspace    = 0.1,
         margin    = 0.1, # an extra fraction of white space to leave around the borders
         minbranch = 0.3)  # set the minimum length for a branch

    ## Add text
    text(x      = rpart.obj,   #
         splits = T,           # If tree are labeled with the criterion for the split
         all    = T,           # If 'TRUE', all nodes are labeled, otherwise just terminal nodes
         use.n  = T,           # Use numbers to annotate
         cex    = font.size)   # Font size
}

plot.rpart.obj(rpart.tree, 0.7)

# partykit
rparty.tree <- as.party(rpart.tree)
rparty.tree
fit.rparttree <- rpart(formula.all, method="anova", data=tr)
printcp(fit.rparttree) # display the results
plotcp(fit.rparttree) # visualize cross-validation results
summary(fit.rparttree) # detailed summary of splits
# create additional plots
par(mfrow=c(1,2)) # two plots on one page
rsq.rpart(fit.rparttree) # visualize cross-validation results  
par(mfrow=c(1,1)) 
# plot tree
plot(fit.rparttree, uniform=TRUE, main="Regression Tree for Al+0.5Fe ")
text(fit.rparttree, use.n=TRUE, all=TRUE, cex=.8)
# create attractive postcript plot of tree
post(fit, file = "rtree.ps", title = "Regression Tree for Al+0.5Fe ")

## *** BUILD RANDOM FOREST ***
mod.rf <- randomForest(formula.all,
                       tr,
                       ntree=500,
                       importance=T)
plot(mod.rf)

## *** IMPORTANCE ***
# Variable Importance Plot
varImpPlot(mod.rf,
           sort = T,
           main="Variable Importance",
           n.var=5)
varImpPlot(mod.rf)
# Variable Importance Table
var.imp <- data.frame(importance(mod.rf,type=2))
# make row names as columns
var.imp$Variables <- row.names(var.imp)
var.imp[order(var.imp$MeanDecreaseGini,decreasing = T),]

## *** Predict Response Variable Value using Random Forest ***
# Predicting response variable
tr$predicted.rf <- predict(mod.rf,tr)
plot(tr$predicted.rf,tr$Al05Fe)
cor(tr$predicted.rf,tr$Al05Fe,use="pairwise.complete.obs",method="pearson")
# Predictions at Te
te$predicted.rf <- predict(mod.rf,te)
plot(te$predicted.rf,te$Al05Fe)
cor(te$predicted.rf,te$Al05Fe,use="pairwise.complete.obs",method="pearson")

## *** Confusion Matrix ***
# Load Library or packages
library(e1071)
library(caret)
# Create Confusion Matrix
confusionMatrix(data=tr$predicted.rf,
                reference=tr$Al05Fe,
                positive='yes')

# Predicting response variable
cross.sell.val$predicted.response <- predict(cross.sell.rf ,cross.sell.val)
# Create Confusion Matrix
confusionMatrix(data=cross.sell.val$predicted.response,
                reference=cross.sell.val$y,
                positive='yes')

#_____________________________________________________
#
# SEE :: http://trevorstephens.com/kaggle-titanic-tutorial/r-part-5-random-forests/
#_____________________________________________________
#
library(party)
set.seed(415)
varImpPlot(mod.rf)
fit <- cforest(formula.all,
               tr,
               controls=cforest_unbiased(ntree=500, mtry=3))
te$predicted.fitcfor <- predict(fit, te, OOB=TRUE, type = "response")
cor(te$predicted.fitcfor,te$Al05Fe,use="pairwise.complete.obs",method="pearson")
plot(te$predicted.fitcfor,te$Al05Fe)


