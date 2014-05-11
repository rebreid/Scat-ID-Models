#Load the data
Scat<-read.csv(file="DNAid.csv",sep=",",header=TRUE) 

#Subset data to exclude mixed, missing and spotted skunk
Scat2<-subset(Scat, Species!="spotted skunk")
Scat2<-subset(Scat2, Species!="mixed")
Scat2<-subset(Scat2, Species!="failed")
Scat2<-subset(Scat2,select=-c(18:19)) #drop out the isotope data for now
Scat2<-subset(Scat2,select=-c(13)) #drop out the C:N ratio for now
Scat2<-subset(Scat2,select=-c(1,3:4)) #drop out the sample id, month and site columns for now

#Reset the Species variable such that it doesn't retain the removed subset labels
Scat2$Species<-factor(Scat2$Species)

#Before we can do any modeling, it's important that any categorical variables are seen as factors, not as numeric variables. check the encoding with:
str(Scat2)  #if categorical variables are not factors, use as.factor() to change.
Scat2$Age<-factor(Scat2$Age)

#Now we'll try building a conditional inference tree using the package "Party".
library(party)
speciesID<-ctree(Scat2$Species ~ ., data = Scat2)
speciesID #outputs a description of the model
plot(speciesID) #shows the tree

#to test the fit of the model, plot a table of the prediction errors
table(predict(speciesID), Scat2$Species)

#Estimated class probabilities
tr.pred<-predict(speciesID, type="prob")
tr.pred

## what if we change the tree control type to "MonteCarlo"?
Scat3<-na.omit(Scat2)
speciesID<-ctree(Scat3$Species ~ ., data = Scat3, controls = ctree_control(testtype = c("MonteCarlo"), nresample=9999))
plot(speciesID)
table(predict(speciesID), Scat3$Species)

##Now change the control type to "Bonferroni"
speciesID<-ctree(Scat2$Species ~ ., data = Scat2, controls = ctree_control(testtype = c("Univariate"), mincriterion = 0.95, minbucket = 11, savesplitstats = TRUE))
plot(speciesID)
table(predict(speciesID), Scat2$Species)

#Now let's try creating a Random Forest Model
#start by checking out the data structure, be sure that categorical variables are listed as factors and binary variables are listed as integers.
str(Scat2)

#Then split the data into a training and a test set
ind<-sample(2,nrow(Scat3), replace=TRUE, prob=c(0.7,0.3))
trainData<-Scat3[ind==1,]
testData<-Scat3[ind==2,]

#Generate random forest learning trees
library(randomForest)
set.seed(1385)
Scat_rf<-randomForest(Species~.,data=trainData,ntree=1000,proximity=TRUE)
table(predict(Scat_rf),trainData$Species)
print(Scat_rf)
attributes(Scat_rf)
plot(Scat_rf)
importance(Scat_rf)
varImpPlot(Scat_rf)

#extract the predicted classes from the random forest model
rfTrainPred <- predict(Scat_rf, trainData, type = "prob")
head(rfTrainPred)
trainData$RFprobCoy <- rfTrainPred[,"coyote"]
trainData$RFprobBob <- rfTrainPred[,"bobcat"]
trainData$RFprobFox <- rfTrainPred[,"gray fox"]
trainData$RFclass <- predict(Scat_rf, trainData)

#extract the confusion matrix
confusionMatrix(data = trainData$RFclass, reference = trainData$Species)

#Test the built random forest for the test data
ScatPred<-predict(Scat_rf,newdata=testData)
table(ScatPred, testData$Species)
plot(margin(Scat_rf,testData$Species))

outlier(Scat_rf)
plot(outlier(Scat_rf), type="h",
     col=c("red", "green", "blue")[as.numeric(Scat$Species)])

#extract the predicted classes from the random forest model
rfTestPred <- predict(Scat_rf, testData, type = "prob")
head(rfTestPred)
testData$RFprobCoy <- rfTestPred[,"coyote"]
testData$RFprobBob <- rfTestPred[,"bobcat"]
testData$RFprobFox <- rfTestPred[,"gray fox"]
testData$RFclass <- predict(Scat_rf, testData)

#extract the confusion matrix
confusionMatrix(data = testData$RFclass, reference = testData$Species)

#Plot probability distributions of scat ids by species
library(ggplot2)
require(gridExtra)
bobcat<-testData$RFprobBob
coyote<-testData$RFprobCoy
grayfox<-testData$RFprobFox
predictions<-data.frame(testData$Species,bobcat,coyote,grayfox)
colnames(predictions)<-c("Species","Bobcat","Coyote","Fox")
plot1<-ggplot(predictions, aes(x=Bobcat, fill=Species)) + geom_density(alpha=0.2) + xlab("Probability - Bobcat")
plot2<-ggplot(predictions, aes(x=Coyote, fill=Species)) + geom_density(alpha=0.2) + xlab("Probability - Coyote")
plot3<-ggplot(predictions, aes(x=Fox, fill=Species)) + geom_density(alpha=0.2) + xlab("Probability - Gray Fox")
grid.arrange(plot1,plot2,plot3, ncol=3)

###Try the same thing, but using a different R package:
### This is the same type of model, just using a different package to execute it.
#begin by setting the controls for the random forest
data.controls<-cforest_unbiased(ntree=2000, mtry=3.5) #the ntree argument controls the overall number of trees in the forest, the mtry argument controls the number of randomly preselected predictor variables for each split. The square root of the number of variables is the suggested default value for mtry in the literature, we have 9 variables, so 3 is an appropriate choice.

#set a random seed for the forest run, any random number will do:
set.seed(894)
mycforest<-cforest(Species~., data=trainData, controls=data.controls, savesplitstats=TRUE)
mycforest

#Look at some trees in the forest (code from Strobl et al. 2009 supplementary material)
xgr <- 2
grid.newpage()
cgrid <- viewport(layout = grid.layout(xgr, xgr), name = "cgrid")
pushViewport(cgrid)
for (j in 1:xgr) {
  pushViewport(viewport(layout.pos.col = i, layout.pos.row = j)) 
  tr <- party:::prettytree(mycforest@ensemble[[i + j * xgr]], names(mycforest@data@get("input")))
  plot(new("BinaryTree", tree = tr, data = mycforest@data, responses = mycforest@responses), newpage = FALSE, pop = FALSE, type="simple") 
  upViewport()
}

#compute and plot the permutation importance of each predictor variable
set.seed(560)
myvarimp<-varimp(mycforest)
myvarimp
dotchart(sort(myvarimp))

#conditional = true
myvarimp<-varimp(mycforest, conditional=TRUE)
dotchart(sort(myvarimp))
###



