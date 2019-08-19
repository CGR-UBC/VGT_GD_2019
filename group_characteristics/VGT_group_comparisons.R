# Eve Limbrick-Oldfield 14/08/2019
# Statistics for Table 1 [pubmed link]

library(pastecs)

rfromwilcox<-function(wilcoxModel,N){
  z<-qnorm(wilcoxModel$p.value/2)
  r<-z/sqrt(N)
  cat(wilcoxModel$data.name,"Effect size, r = " , r)
}

# edit
setwd("***")
data=read.csv("***.csv")

data$Group <- ifelse(data$participant > 200, 1,0)
data<-within (data, {Group<-factor(Group,levels=0:1,labels = c("HC","GD"))})
data$Druguser<-ifelse(data$DAST > 0, 1,0)
data$Druguser[is.na(data$Druguser)] <- 0 # Replace NA with zero
data_smokers<-subset(data,Smoker==1)
data_drug_users<-subset(data,DAST>0)


by(data$PGSI,data$Group,stat.desc,basic=TRUE,norm=TRUE)

by(data$Win_coins,data$Group,stat.desc,basic=TRUE,norm=TRUE)
by(data$Loss_coins,data$Group,stat.desc,basic=TRUE,norm=TRUE)


################ NON-PARAMETRIC #######################

by(data$Win_coins,data$Group,stat.desc,basic=TRUE,norm=TRUE)
Win_coins_model<-wilcox.test(data$Win_coins~data$Group,paired=FALSE,exact=FALSE)
Win_coins_model
rfromwilcox(Win_coins_model,nrow(data))

by(data$Loss_coins,data$Group,stat.desc,basic=TRUE,norm=TRUE)
Loss_coins_model<-wilcox.test(data$Loss_coins~data$Group,paired=FALSE,exact=FALSE)
Loss_coins_model
rfromwilcox(Loss_coins_model,nrow(data))

by(data$Age,data$Group,stat.desc,basic=TRUE,norm=TRUE)
Agemodel<-wilcox.test(data$Age~data$Group,paired=FALSE,exact=FALSE)
Agemodel
rfromwilcox(Agemodel,nrow(data))

by(data$Age,data$Group,stat.desc,basic=TRUE,norm=TRUE)
Agemodel<-wilcox.test(data$Age~data$Group,paired=FALSE,exact=FALSE)
Agemodel
rfromwilcox(Agemodel,nrow(data))

by(data$DASS,data$Group,stat.desc,basic=TRUE,norm=TRUE)
DASS_TOTALmodel<-wilcox.test(data$DASS_TOTAL~data$Group,paired=FALSE,exact=FALSE)
DASS_TOTALmodel
rfromwilcox(DASS_TOTALmodel,nrow(data))

by(data$Total_GRCS,data$Group,stat.desc,basic=TRUE,norm=TRUE)
Total_GRCSmodel<-wilcox.test(data$Total_GRCS~data$Group,paired=FALSE,exact=FALSE)
Total_GRCSmodel
rfromwilcox(Total_GRCSmodel,nrow(data))

by(data$AUDIT,data$Group,stat.desc,basic=TRUE,norm=TRUE)
AUDITmodel<-wilcox.test(data$AUDIT~data$Group,paired=FALSE,exact=FALSE)
AUDITmodel
rfromwilcox(AUDITmodel,nrow(data))

by(data_smokers$Fagerstrom,data_smokers$Group,stat.desc,basic=TRUE,norm=TRUE)
Fagerstrommodel<-wilcox.test(data_smokers$Fagerstrom~data_smokers$Group,paired=FALSE,exact=FALSE)
Fagerstrommodel
rfromwilcox(Fagerstrommodel,nrow(data_smokers))

by(data_drug_users$DAST,data_drug_users$Group,stat.desc,basic=TRUE,norm=TRUE)
DASTmodel<-wilcox.test(data_drug_users$DAST~data_drug_users$Group,paired=FALSE,exact=FALSE)
DASTmodel
rfromwilcox(DASTmodel,nrow(data_drug_users))

################ PARAMETRIC #######################

by(data$Verbal.IQ,data$Group,stat.desc,basic=TRUE,norm=TRUE)

Verbal.IQmodel<-t.test(Verbal.IQ~Group,data=data,paired=FALSE)
Verbal.IQmodel
t<-Verbal.IQmodel$statistic[[1]]
df<-Verbal.IQmodel$parameter[[1]]
r<-sqrt(t^2/(t^2+df))
r

# chi-squre tests
temp<-subset(data,Gender<3)
Gender_table <-table(temp$Group,temp$Gender)
chisq.test(Gender_table)

Smoker_table <-table(data$Group,data$Smoker)
chisq.test(Smoker_table)

Druguser_table <-table(data$Group,data$Druguser)
chisq.test(Druguser_table)