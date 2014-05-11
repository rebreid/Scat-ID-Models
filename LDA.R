#Load the data
Scat<-read.csv(file="DNAid.csv",sep=",",header=TRUE) 

#Subset data to exclude mixed, missing and spotted skunk
Scat2<-subset(Scat, Species!="spotted skunk")
Scat2<-subset(Scat2, Species!="mixed")
Scat2<-subset(Scat2, Species!="failed")
Scat2<-subset(Scat2,select=-c(18:19)) #drop out the isotope data for now
Scat2<-subset(Scat2,select=-c(13)) #drop out the C:N ratio for now
              
#Reset the Species variable such that it doesn't retain the removed subset labels
Scat2$Species<-factor(Scat2$Species)
Scat2<-na.omit(Scat2)
str(Scat2) #check the data structure 

#assess assumptions of LDA before starting
#Are the variance/covariance matrices for each group equivalent (they should be)?
#Separate the data by groups
bobcats<-Scat2[Scat2$Species=="bobcat",]
coyotes<-Scat2[Scat2$Species=="coyote",]
foxes<-Scat2[Scat2$Species=="gray fox",]

#standardize the continuous variables in each group
StandBob<-as.data.frame(scale(bobcats[7:12]))
StandCoy<-as.data.frame(scale(coyotes[7:12]))
StandFox<-as.data.frame(scale(foxes[7:12]))

#Look at the correlation matrices
BobCor<-cor(StandBob)
CoyCor<-cor(StandCoy)
FoxCor<-cor(StandFox)
BobCor
CoyCor
FoxCor
#These are not exactly equivalent - may need to proceed with Quadratic Discriminant Analysis, 
#which allows for these to differ. Start with LDA, though, and compare! Also note that taper and TI are highly correlated
#probably best to include just one or the other in the model, not both.

### Linear Discriminant Analysis ###
#Begin by log transforming the variable that needs it (determined during data inspection).
logMass<-log(Scat2$Mass)
#put the transformed data into the same data frame as the other numeric variables
NumScat<-Scat2[7:12]
NumScat<-data.frame(NumScat,logMass)

#scale the data
StanScat <- as.data.frame(scale(NumScat))
#add the sample id column back
StanScat$ID<-Scat2$ID
#now put the standardised numeric data back together with the other data
Scat2<-subset(Scat2, select=-c(7:12)) #drops the unstandardized versions of the numerical columns from the dataset Scat2
Scat3<-merge(Scat2,StanScat, by="ID")
str(Scat3)

#Compute Wilk's lambda to evaluate whether or not the species are significantly different
X<-as.matrix(Scat3[,11:17]) #isolate the standardized continuous variables
scat.manova<-manova(X~Scat3$Species)
scat.wilks<-summary(scat.manova,test="Wilks")
scat.wilks

#Split the data into a training and a test set
ind<-sample(2,nrow(Scat3), replace=TRUE, prob=c(0.7,0.3))
trainData<-Scat3[ind==1,]
testData<-Scat3[ind==2,]

#Perform linear discriminant analysis
library(MASS)
#Note - set prior probability based on the structure of the dataset. This dataset is about 1/2 bobcat scats and 1/4 each coyote and fox scats.
Scatlda<-lda(Species ~ Number + Length + Diameter + Taper + logMass, prior = c(1/2, 1/4, 1/4), data=trainData, CV=F) 
Scatlda
plot(Scatlda)

#plot the correlation coefficients between the first two axes 
#and each of the original variables
correlation<-cor(x = trainData[, c(11:15,17)], y = predict(Scatlda)$x)
plot(c(-0.8, 0.8), c(-0.8, 0.8), type = "n", 
     xlab = "Discriminant function 1", ylab = "Discriminant function 2")
lines(c(-1, 1), c(0, 0), lty = 2)
lines(c(0, 0), c(-1, 1), lty = 2)
library(ade4)
s.arrow(correlation[, 1:2], grid = F, add.plot = T, addaxes = T, clabel = 1.5)
correlation

Scatlda$scaling[,1] #these are the loadings for the 1st discriminant function
Scatlda.values<-predict(Scatlda)

#View the separation between the groups on each of the LDA axes
ldahist(data = Scatlda.values$x[,1], g=trainData$Species) #LDA1
ldahist(data = Scatlda.values$x[,2], g=trainData$Species) #LDA2
par(op)

#Plot (same as the plot produced earlier, but with numbers for the groups instead of species names...)
LD1<-Scatlda.values$x[,1]
LD2<-Scatlda.values$x[,2]
plot(LD1,LD2,xlab="first linear discriminant",ylab="second linear discriminant",type="n")
text(cbind(LD1,LD2),labels=unclass(trainData$Species))

#extract the predicted classes from the qda model for the training data
ldaTrainPred <- predict(Scatlda, trainData)
names(ldaTrainPred)
head(ldaTrainPred$class)
head(ldaTrainPred$posterior)
trainData$LDAprobCoy <- ldaTrainPred$posterior[,"coyote"] #probability that a scat in the trainData is identified by the model as coyote
trainData$LDAprobBob <- ldaTrainPred$posterior[,"bobcat"] #probability that a scat in the trainData is identified by the model as bobcat
trainData$LDAprobFox <- ldaTrainPred$posterior[,"gray fox"] #probability that a scat in the trainData is identified by the model as gray fox

#Calculate the confusion matrix for the train data
confusionMatrix(data = ldaTrainPred$class, reference = trainData$Species)

#Plot with the centroids drawn in (95% confidence level)
library(ggplot2)
datPred<-data.frame(Species=trainData$Species,Scatlda.values$x)
library(ellipse) 
dat_ell <- data.frame() 
for(g in levels(datPred$Species)){ 
  dat_ell <- rbind(dat_ell, cbind(as.data.frame(with(datPred[datPred$Species==g,], ellipse(cor(LD1, LD2), 
                                                      scale=c(sd(LD1),sd(LD2)), 
                                                      centre=c(mean(LD1),mean(LD2))))),Species=g)) 
} 

ggplot(datPred, aes(x=LD1, y=LD2, col=Species) ) + geom_point( size = 4, aes(color = Species))+theme_bw()+geom_path(data=dat_ell,aes(x=x,y=y,color=Species),size=1,linetype=2) 

#Calculate the centroids by group...
bobcatx<-sum(LD1*(trainData$Species=="bobcat"))/sum(trainData$Species=="bobcat")
bobcaty<-sum(LD2*(trainData$Species=="bobcat"))/sum(trainData$Species=="bobcat")
coyotex<-sum(LD1*(trainData$Species=="coyote"))/sum(trainData$Species=="coyote")
coyotey<-sum(LD2*(trainData$Species=="coyote"))/sum(trainData$Species=="coyote")
foxx<-sum(LD1*(trainData$Species=="gray fox"))/sum(trainData$Species=="gray fox")
foxy<-sum(LD2*(trainData$Species=="gray fox"))/sum(trainData$Species=="gray fox")
bobcent<-c(bobcatx,bobcaty)
coycent<-c(coyotex,coyotey)
foxcent<-c(foxx,foxy)

#Plot probability distributions of scat ids by species
library(ggplot2)
require(gridExtra)
Scatlda.values<-predict(Scatlda)
bobcat<-Scatlda.values$posterior[,1]
coyote<-Scatlda.values$posterior[,2]
grayfox<-Scatlda.values$posterior[,3]
predictions<-data.frame(trainData$Species,bobcat,coyote,grayfox)
colnames(predictions)<-c("Species","Bobcat","Coyote","Fox")
plot1<-ggplot(predictions, aes(x=Bobcat, fill=Species)) + geom_density(alpha=0.2) + xlab("Probability - Bobcat")
plot2<-ggplot(predictions, aes(x=Coyote, fill=Species)) + geom_density(alpha=0.2) + xlab("Probability - Coyote")
plot3<-ggplot(predictions, aes(x=Fox, fill=Species)) + geom_density(alpha=0.2) + xlab("Probability - Gray Fox")
grid.arrange(plot1,plot2,plot3, ncol=3)

#Group effect
summary(manova(as.matrix(trainData[, c(11:15,17)]) ~ trainData$Species), "Wilks")
summary(manova(as.matrix(trainData[, c(11:15,17)]) ~ trainData$Species), "Pillai")
summary(manova(as.matrix(trainData[, c(11:15,17)]) ~ trainData$Species), "Hotelling-Lawley")

#Now assess the accuracy using the test dataset
test<-predict(Scatlda, newdata=testData) 
ct2<-table(testData$Species,test$class)
diag(prop.table(ct2,1))
sum(diag(prop.table(ct2)))
class.table<-table(testData$Species,test$class)
class.table

#or another way...
ldaTestPred <- predict(Scatlda, testData)
names(ldaTestPred)
head(ldaTestPred$class)
head(ldaTestPred$posterior)
testData$LDAprobCoy <- ldaTestPred$posterior[,"coyote"] #probability that a scat in the testData is identified by the model as coyote
testData$LDAprobBob <- ldaTestPred$posterior[,"bobcat"] #probability that a scat in the testData is identified by the model as bobcat
testData$LDAprobFox <- ldaTestPred$posterior[,"gray fox"] #probability that a scat in the testData is identified by the model as gray fox

#Calculate the confusion matrix for the test data
confusionMatrix(data = ldaTestPred$class, reference = testData$Species)

#Plot probability distributions of scat ids by species
library(ggplot2)
require(gridExtra)
bobcat<-test$posterior[,1]
coyote<-test$posterior[,2]
grayfox<-test$posterior[,3]
predictions<-data.frame(testData$Species,bobcat,coyote,grayfox)
colnames(predictions)<-c("Species","Bobcat","Coyote","Fox")
plot1<-ggplot(predictions, aes(x=Bobcat, fill=Species)) + geom_density(alpha=0.2) + xlab("Probability - Bobcat")
plot2<-ggplot(predictions, aes(x=Coyote, fill=Species)) + geom_density(alpha=0.2) + xlab("Probability - Coyote")
plot3<-ggplot(predictions, aes(x=Fox, fill=Species)) + geom_density(alpha=0.2) + xlab("Probability - Gray Fox")
grid.arrange(plot1,plot2,plot3, ncol=3)

#Extract means and standard deviations of LDA1 and LDA2 by group
#first make a data frame containing the species assignments and the LDA values
LDs<-data.frame(Scatlda.values$class, Scatlda.values$x[,1], Scatlda.values$x[,2]) 
colnames(LDs)<-c("Species","LD1","LD2")

#Are the centroids different?
aov1<-aov(LD1~Species, data=LDs)
aov2<-aov(LD2~Species, data=LDs)
summary(aov1)
summary(aov2)

par(mfrow=c(1,2))
library(gplots)
boxplot2(LD1~Species, data=LDs)
boxplot2(LD2~Species, data=LDs)

TukeyHSD(aov1)
TukeyHSD(aov2)

#attach posterior predictions to the scat dataset
Scat4<-predict(Scatlda, newdata=Scat3) #make predictions for the entire dataset
Post<-data.frame(Scat3$ID,Scat4$posterior[,1],Scat4$posterior[,2],Scat4$posterior[,3])
colnames(Post)<-c("ID","Bobcat","Coyote","Fox")
Scat5<-merge(Scat3,Post, by="ID")
str(Scat5)

#Seprate by groups
bobcats<-Scat5[Scat5$Species=="bobcat",]
coyotes<-Scat5[Scat5$Species=="coyote",]
foxes<-Scat5[Scat5$Species=="gray fox",]

#Pull out the highly identifiable scats
posbobs<-subset(bobcats, Bobcat>=0.8)
poscoys<-subset(coyotes, Coyote>=0.8)
posfox<-subset(foxes, Fox>=0.8)
negfox<-subset(foxes, Fox<=0.8)
negbob<-subset(bobcats, Bobcat<=0.8)
negcoy<-subset(coyotes, Coyote<=0.8)

#take a look at some of them
posbobs
poscoys

##Create the lift curve using the "lift" function in the "caret" package
library(caret)

#Do this for coyotes first - the classification variable should have only two levels, so need a new column with "coyote" or "not coyote"
coyote<-testData$Species
coyote[coyote=="bobcat"]<-"gray fox"
testData<-data.frame(testData,coyote)
liftCurve1 <- lift(coyote ~ testData$LDAprobCoy, data = testData)
liftCurve1
plot(liftCurve1)

#Then for the bobcats; start by creating a column "bobcat" with records labeled as bobcat or gray fox
bobcat<-testData$Species
bobcat[bobcat=="coyote"]<-"gray fox"
testData<-data.frame(testData,bobcat)
liftCurve2 <- lift(bobcat ~ testData$LDAprobBob, data = testData)
liftCurve2
plot(liftCurve2)

#And finally for gray foxes; start by creating a column "grayfox" with records labeled as gray fox and coyote
grayfox<-testData$Species
grayfox[grayfox=="bobcat"]<-"coyote"
testData<-data.frame(testData,grayfox)
liftCurve3 <- lift(grayfox ~ testData$LDAprobFox, data = testData)
liftCurve3 #####Still not quite right, this is event "coyote" rather than event "gray fox"
plot(liftCurve3) 

##Compare the LDA model with a QDA model
Scatqda<-qda(Species ~ Number + Length + Diameter + Taper + logMass, prior = c(1/2, 1/4, 1/4), data=trainData, CV=F) 
Scatqda

#extract the predicted classes from the qda model for the training data
qdaTrainPred <- predict(Scatqda, trainData)
names(qdaTrainPred)
head(qdaTrainPred$class)
head(qdaTrainPred$posterior)
trainData$QDAprobCoy <- qdaTrainPred$posterior[,"coyote"] #probability that a scat in the trainData is identified by the model as coyote
trainData$QDAprobBob <- qdaTrainPred$posterior[,"bobcat"] #probability that a scat in the trainData is identified by the model as bobcat
trainData$QDAprobFox <- qdaTrainPred$posterior[,"gray fox"] #probability that a scat in the trainData is identified by the model as gray fox

#Calculate the confusion matrix for the train data
confusionMatrix(data = qdaTrainPred$class, reference = trainData$Species)

#Do the same for the test data
qdaTestPred <- predict(Scatqda, testData)
names(qdaTestPred)
head(qdaTestPred$class)
head(qdaTestPred$posterior)
testData$QDAprobCoy <- qdaTestPred$posterior[,"coyote"] #probability that a scat in the trainData is identified by the model as coyote
testData$QDAprobBob <- qdaTestPred$posterior[,"bobcat"] #probability that a scat in the trainData is identified by the model as bobcat
testData$QDAprobFox <- qdaTestPred$posterior[,"gray fox"] #probability that a scat in the trainData is identified by the model as gray fox

#Calculate the confusion matrix for the test data
confusionMatrix(data = qdaTestPred$class, reference = testData$Species)

##Create the lift curve using the "lift" function in the "caret" package
library(caret)

#Do this for coyotes first - the classification variable should have only two levels, so need a new column with "coyote" or "not coyote"
liftCurve1 <- lift(coyote ~ testData$QDAprobCoy, data = testData)
liftCurve1
plot(liftCurve1)

#Then for the bobcats; start by creating a column "bobcat" with records labeled as bobcat or gray fox
liftCurve2 <- lift(bobcat ~ testData$QDAprobBob, data = testData)
liftCurve2
plot(liftCurve2)

#compare the new class probabilities
library(caret)
calCurve <- calibration(testData$coyote ~ testData$QDAprobCoy + testData$LDAprobCoy, data = testData)
xyplot(calCurve, auto.key=TRUE)
