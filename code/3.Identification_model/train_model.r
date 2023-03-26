ls()
rm(list=ls(all=TRUE))

library(splines)
library(survival)
library(survminer)

#==data1.Clinic Information
PRJEB_Ann_train <- read.delim("/work/dy/project/2020CircRNA/temp_github/data/SampleClinicInformation.txt")
all_preID <- as.character(PRJEB_Ann_train[PRJEB_Ann_train$PrePost=="PRE",]$SRR)
PD1_preID <- as.character(PRJEB_Ann_train[PRJEB_Ann_train$Treatment=="PD1" & PRJEB_Ann_train$PrePost=="PRE",]$SRR)
ipiPD1_preID <- as.character(PRJEB_Ann_train[PRJEB_Ann_train$Treatment=="ipiPD1" & PRJEB_Ann_train$PrePost=="PRE",]$SRR)

sur_survival=PRJEB_Ann_train[match(all_preID,PRJEB_Ann_train$SRR,nomatch=0),]

#=data2.
all_CircDataM88sample_89204Circ_exp_sub5350 <- readr::read_rds("/work/dy/project/2020CircRNA/temp_github/data/circRNA_exp.rds")

#=data3.
FPS_MK_SIG=read.table("/work/dy/project/2020CircRNA/temp_github/data/pre_PFS25_CoxMK_SIG.tab")
FPS_MK_SIG=as.character(FPS_MK_SIG[,1])

pre_exp_PFS=all_CircDataM88sample_89204Circ_exp_sub5350[FPS_MK_SIG,all_preID]
pre_exp_PFSsur=cbind(sur_survival[,c("PFST", "PFS")],t(pre_exp_PFS))

#=================== lasso
#=================== lasso
library(glmnet)
x <- t(as.matrix(pre_exp_PFS))
y <- sur_survival[,c("PFST", "PFS")]
names(y) <- c("time","status")
y$time <- as.double(y$time)
y$status <- as.double(y$status)
y <- as.matrix(survival::Surv(y$time,y$status))

lasso <- glmnet(x, y, family = "cox", alpha = 1)
p2 <-plot(lasso, xvar = "lambda", label = TRUE)
set.seed(123)
fitCV <- cv.glmnet(x, y, family = "cox",
                   type.measure = "deviance",
                   nfolds = 5)
plot(fitCV)

coefficient <- coef(fitCV, s= fitCV$lambda.min)
Active.Index <- which(as.numeric(coefficient) !=0)
active.coefficients <- as.numeric(coefficient)[Active.Index]
sig_gene_multi_cox <- rownames(coefficient)[Active.Index]
sig_gene_multi_cox

###= cox 
Mulcox <- coxph(Surv(PFST,PFS) ~ `chr1_236178413_236183906_+_GPR137B` + 
									`chr2_202755537_202759353_+_FAM117B` + 
									`chr12_88148287_88176319_+_TMTC3` + 									
									`chr11_6941570_6955782_+_ZNF215` , data =  pre_exp_PFSsur)
Mulcox_coef <- coef(Mulcox)
Mulcox_Coxp <- signif(as.numeric(summary(Mulcox)$coefficients[,c( "Pr(>|z|)")]),digit=4)
Mulcox_coef_used=Mulcox_coef[Mulcox_Coxp<0.05]

options(scipen = 1)
ggforest(
  model = Mulcox, data = pre_exp_PFSsur,
  main = "Hazard ratio",
  cpositions = c(0.10, 0.22, 0.4),
  fontsize = 1.0,
  refLabel = "1", noDigits = 4
)

temp_PFS_coxresult=summary(Mulcox)$coefficients[summary(Mulcox)$coefficients[,5]<0.05,]
rownames(temp_PFS_coxresult)=substr(rownames(temp_PFS_coxresult),2,(nchar(rownames(temp_PFS_coxresult))-1))

#=====================score 
tempexp=all_CircDataM88sample_89204Circ_exp_sub5350[,all_preID]
temp_sur=PRJEB_Ann_train[match(all_preID,PRJEB_Ann_train$SRR,nomatch=0),]		

PFS_cox_model_lasso_score=apply(tempexp[rownames(temp_PFS_coxresult),],2,function(x){sum(as.numeric(x)*as.numeric(temp_PFS_coxresult[rownames(temp_PFS_coxresult),1]),na.rm=T)})
temp_exp=PFS_cox_model_lasso_score
gene <- temp_exp
gene <- unlist(gene)
group <- ifelse(gene > median(gene,na.rm = T),'high', 'low')
survival_dat <- data.frame(group = group, 
							 groupHR = relevel(factor(group), ref ="low"), 
                             status = temp_sur$PFS,
                             time = temp_sur$PFST,
							 gene = gene,						 
                             stringsAsFactors = F)
		
# harzad ratio
model1 <- survival::coxph(Surv(time, status) ~ gene, data = survival_dat,  na.action=na.exclude)
coef = signif(as.numeric(summary(model1)$coefficients[1,c( "coef")]),digit=4)		
HR = signif(as.numeric(summary(model1)$coefficients[1,c( "exp(coef)")]),digit=4)
Coxp <- signif(as.numeric(summary(model1)$coefficients[1,c( "Pr(>|z|)")]),digit=4)
HR_detail = summary(model1)
CILow =  HR_detail$conf.int[,"lower .95"]
CIHigh =  HR_detail$conf.int[,"upper .95"]
CI =  paste(signif(HR_detail$conf.int[,"lower .95"],2),'-',signif(HR_detail$conf.int[,"upper .95"],2),sep=" ") 			
fit <- survfit(Surv(time, status) ~ group, data = survival_dat)	
   
maxstat_test_result=maxstat.test(Surv(time, status) ~ gene, data=survival_dat, smethod="LogRank", pmethod="HL")
survival_dat$max_group <- ifelse(survival_dat$gene > maxstat_test_result$estimate,"High","Low")
# maxstat MK 
logRankTestmaxstat=survdiff(Surv(time, status) ~ max_group, data = survival_dat, na.action=na.exclude) ##log-rank test using survdiff(formula)
KM_pvaluemaxstat <- signif(1-pchisq(logRankTestmaxstat[[5]],1))		
		
fit <- survfit(Surv(time, status) ~ max_group, data = survival_dat)	
p_sur_MK_max=ggsurvplot(fit, data = survival_dat,palette = c("red", "blue"),
                    pval = F, title=paste('log-rank'," p = ",signif(as.numeric(KM_pvaluemaxstat),digits = 2),
                                          "\nHR ",HR, " [95% CI:",CI,"]",sep=""),
                    fun = function(y) y*100,legend = c(0.2, 0.2), legend.title=F,font.legend=10,font.title=10,
                    legend.labs = c(paste("High (n=",table(survival_dat$max_group)[[1]],")",sep=""),
                                    paste("Low (n=",table(survival_dat$max_group)[[2]],")",sep="")))

p_sur_MK_max$plot <- p_sur_MK_max$plot + xlab("Months")+ theme(legend.title = element_blank(),legend.background = element_blank())

#=====================score -- response
survival_dat$Benefit=temp_sur$Benefit
ggboxplot(survival_dat, x = "Benefit", y = "gene",color='Benefit',add = "jitter",ylab='score',x.text.angle= 0)+theme(legend.position='none')+
          stat_compare_means(method = "wilcox") + scale_color_npg() #scale_fill_manual(values=c('red','blue')) + labs(ylab='score') #kruskal.test

#=data4.
tempExp=readRDS('/work/dy/project/2020CircRNA/temp_github/data/exp_CD274_PDCD1.rds')
#========cox 
survival_dat=data.frame(survival_dat,t(tempExp[,rownames(survival_dat)]))
survival_dat$gender=temp_sur$Sex	
survival_dat$Age=as.numeric(temp_sur$Age)
survival_dat$Treatment=temp_sur$Treatment	
final_score_coxmodel=coxph(Surv(time,status ) ~ gene + Age + gender + Treatment + CD274 + PDCD1, data =  survival_dat)

#forest
options(scipen = 1)
ggforest(
			model = final_score_coxmodel, data = survival_dat,
			main = "Hazard ratio",
			cpositions = c(0.10, 0.22, 0.4),
			fontsize = 1.0,
			refLabel = "1", noDigits = 4
			)	

#==================timeROC 
#==================timeROC 
library(timeROC)
ROC<- timeROC(T=survival_dat$time,
                      delta=survival_dat$status,
                      marker=survival_dat$gene,
                      cause=1,
                      weighting="marginal",
                      times=c(12,24),
                      ROC = TRUE,
                      iid = TRUE)
        
        auc_1 = ROC$AUC[[1]]
        auc_2 = ROC$AUC[[2]]
        
        dat = data.frame(tpr_1 = ROC$TP[,1],
                         fpr_1 = ROC$FP[,1],
                         tpr_2 = ROC$TP[,2],
                         fpr_2 = ROC$FP[,2]
                         )
        
        library(ggplot2)        
        ggplot() + 
          geom_smooth(data = dat,aes(x = fpr_1, y = tpr_1),color = "blue",se = F) + 
          geom_smooth(data = dat,aes(x = fpr_2, y = tpr_2),color = "red",se = F)+
          geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
          theme(panel.background=element_rect(colour=NA,fill="white"),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                axis.line.y = element_line(colour = "black"),
                axis.line.x = element_line(colour = "black"))+
          annotate("text",x = .75, y = .25,
                   label = paste("AUC of 12 mouths = ",round(auc_1,2)),color = "blue")+
          annotate("text",x = .75, y = .15,label = paste("AUC of 24 mouths = ",round(auc_2,2)),color = "red")+
          scale_x_continuous(name  = "1-Specificity")+
          scale_y_continuous(name = "Specificity")
        
		