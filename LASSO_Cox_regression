set.seed(1)

genelist = "A+B+C+D+E" # Gene list -> Input your interested genes
genes = unlist(strsplit(genelist, split="+",fixed=T))
matchgene=intersect(genes,colnames(tmp))
genes=matchgene
print(length(genes))
library(glmnet)
S.OS = Surv(time = as.numeric(data[,time]),event= as.numeric(data[,event]))

y<-cbind(time=data[,time],status=data[,event])
#y=Surv(type='right',time=data[,time],event=data[,event])
coxlasso<-glmnet(x=data[,genes],y=y,family='cox') # 0 때문에 오류 -> 0보다 큰 값으로 진행할 것 (filtering 으로 제거)
dev.new()
plot(coxlasso,label=T)
# 1. Standardization 문제 -> False 사용해도 될지?
# 2. 재현성 문제 -> LOOCV 로 해결해야할지? -> LOOCV 사용 but 이에 따른 문제점도 발생할 수 있음 (overfit)
# y 값으로 data[,event], data[,time] 이용하기
cv.lasso_OS<-cv.glmnet(x=as.matrix(data[,genes]),y,family='cox',alpha=1,nfolds=149,parallel=T) # 수정 필요 -> 재현성 이슈가 있어서 LOOCV
# Ref:https://stats.stackexchange.com/questions/572211/questions-related-to-survival-analysis-and-lasso-cox-regression
cv.lasso_OS
coef(cv.lasso_OS,s="lambda.min") # lambda.min, lambda.1se (1se is more strict than lambda.min)

plot(coxlasso,'lambda')
#plot(coxlasso,'norm',label=TRUE) # abs(Coefficient) 를 모두 더한 값이 L1 norm
plot(cv.lasso_OS)
#gene300=fread('Gene_Signature_LASSO_300genes.txt',stringsAsFactors = F)
#genelst<-str_c(gene300$Gene,collapse='+')
#genelst_2<-str_c(gene300$Gene)

### type measure = 'C' -> glmnet 한번 찾아보기
cv.lasso_OS
coef(cv.lasso_OS,s='lambda.min')
predict.glmnet(cv.lasso_OS,as.matrix(data[,genes])) # why error
predict(cv.lasso_OS,as.matrix(data[,genes]))

#predict.coxph(cv.lasso_OS_fit,data[,genes])

fit<-glmnet(x=data[,genes],y=y,family='cox')

#cv.lasso_OS_fit<-cv.glmnet(x=as.matrix(data[,genes]),y,family='cox',standardize=F,alpha=1,nfolds=5,parallel=T,type.measure='C') # 수정 필요
#cv.lasso_OS_fit #-> LASSO -> k-fold CV 대신 LOOCV 로 사용 ( 재현성 문제 )
# Error 발생하면 survival:::predict.coxph 로 사용하기 
#coxlasso<-glmnet(x=data[,genes],y=y,family='cox',type.measure='C',standardize = F)
#plot(coxlasso,'lambda')
#coef(cv.lasso_OS_fit)
#plot(cv.lasso_OS_fit)
#tmp_coeffs<-coef(cv.lasso_OS,s='lambda.1se')
tmp_coeffs<-coef(cv.lasso_OS,s='lambda.min') # 더 많이나옴 feature
res<-data.frame(name=tmp_coeffs@Dimnames[[1]][tmp_coeffs@i+1],coefficient=tmp_coeffs@x)
res


risk_lst<-c()

risk_table<-data[,res$name]*res$coefficient

for (i in seq(1,149)){
  risk_lst<-c(risk_lst,sum(risk_table[i,]))
}
data$risk_score<-risk_lst
#write.csv(data,'20230718_Risk_score_calculated.csv')
#write.csv(data[,c('sampleID','risk_score',genes)],'20230718_Parsed_risk_score_calculated.csv')
low_group<-data[data$risk_score<median(data$risk_score),c('sampleID','risk_score',genes)]
high_group<-data[data$risk_score>=median(data$risk_score),c('sampleID','risk_score',genes)]
#write.csv(low_group,'20230718_lowrisk_group.csv')
#write.csv(high_group,'20230718_highrisk_group.csv')

label_risk<-c()
for (i in seq(1,149)){
  if (data[i,'risk_score']<median(data$risk_score)){
    label_risk<-c(label_risk,0)
  } else{
    label_risk<-c(label_risk,1)
  }
}
data$high_risk<-label_risk
